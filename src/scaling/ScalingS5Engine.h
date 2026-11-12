#ifndef SCALING_S5_ENGINE_H
#define SCALING_S5_ENGINE_H

// =============================================================================
// ScalingS5Engine -- the paper-faithful §5 implementation of Edmonds' search for
// max-weight PERFECT matching (Duan-Pettie-Su arXiv:1411.1919, §5). Built as a
// SEPARATE, self-contained engine so the proven ScalingMatching::compute()
// (40/40 vs LEMON) stays the correctness ORACLE we gate against.
//
// The three §5 data structures (docs/Scaling_S5_DataStructures.md):
//   1. UNION-FIND over OUTER blossoms: find(v) = the outer node (super-vertex)
//      containing v. Replaces the mBParent pointer-walk of the old engine. union
//      on blossom formation. This makes "same outer blossom?" / "same tree?" O(a)
//      and ALWAYS valid under merges -- the thing the bolted-on BucketPQ lacked,
//      which is why deferred events corrupted there.
//   2. SPLIT-FINDMIN over INNER blossoms: maintains, per inner vertex, the min
//      slack of an incident edge from an outer vertex (for dissolution min-edges).
//      (Stage 3.3; correctness baseline scans instead.)
//   3. PRIORITY QUEUE (BucketPQ) of grow / blossom-formation / dissolution events.
//      (Stage 3.2; correctness baseline scans for the min event instead.)
//
// BUILD ORDER (each gated vs the oracle):
//   3.1 union-find super-nodes + shrink/augment rewritten on them; duals
//       materialized; NEXT-EVENT found by an O(m) scan (isolates UF correctness
//       from the PQ). <-- THIS FILE STARTS HERE.
//   3.2 replace the scan with the BucketPQ event scheduler.
//   3.3 replace the inner min-edge scan with split-findmin.
//
// Integer weights/duals (long long). Maximization direction (Property 1/Crit 1).
// C++98 / -ansi. Vertices 0..n-1; blossom super-nodes get ids >= n.
// =============================================================================

#include <vector>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <utility>
#include <algorithm>
#include "Types.h"
#include "BucketPQ.h"
#ifdef PS5_SFM
#include "SplitFindmin.h"   // §5.1#2 split-findmin (ConcatenableQueue-backed) for no-relabel dissolution
#endif
#ifdef PS5_CQ_EVENTS
#include "ConcatenableQueue.h"   // K_Scope Stage 1: per-blossom O(log n) event queues (shadow, then driver)
#endif
// ---------------------------------------------------------------------------------------------
// THE DEFAULT STACK (2026-07-12 DEFAULT FLIP -- the N_Scope arc's endpoint). Flagless builds get
// the full validated 10-flag configuration: gates 348/594 green, LEMON-int64-exact at every
// refereed scale incl. the certified 1M (wt=4434326679609, 210.0s), fastest-ever on all three
// regimes (W=4 wall dead: 244.6s -> 6.4s @51200). Opt out wholesale: -DPS5_NO_BEST_STACK
// (restores the pre-flip base engine); individual opt-outs: the usual PS5_NO_* / omitting flags
// via PS5_NO_BEST_STACK + explicit -D picks.
#ifndef PS5_NO_BEST_STACK
#ifndef PS5_VTX_BEST
#define PS5_VTX_BEST 1           // DPS Fig-9: keep the least-slack edge per target, discard others
#endif
#if !defined(PS5_EVT_HORIZON) && !defined(PS5_NO_EVT_HORIZON)
#define PS5_NO_EVT_HORIZON 1     // horizon OFF in the stack (marker rescans backfire on giant blossoms)
#endif
#ifndef PS5_NODE_Y
#define PS5_NODE_Y 1             // faithful §5.2 member-y through the outer blossom's role
#endif
#ifndef PS5_FLIP_DIRTY
#define PS5_FLIP_DIRTY 1         // event-driven role sync (O(#flips), the persistent-clock enabler)
#endif
#if !defined(PS5_LEAF_NOVEC) && !defined(PS5_LEAF_SHADOW) && !defined(PS5_LEAF_SPLICE)
#define PS5_LEAF_NOVEC 1         // Appx-C leaf fringe as intrusive lists; vectors deleted
#endif
#ifndef PS5_PERSIST_SHRINK
#define PS5_PERSIST_SHRINK 1     // persistent boundaries extend to shrink phases (teardowns gone)
#endif
#ifndef PS5_NESTED_INTACT
#define PS5_NESTED_INTACT 1      // in-place relink through nested + intact sub-blossoms
#endif
#ifndef PS5_DEEP_INTACT
#define PS5_DEEP_INTACT 1        // uncapped-depth struct-eligibility with per-batch memo
#endif
#ifndef PS5_INTACT_POLE
#define PS5_INTACT_POLE 1        // on-arc intact rotation to the MATCHED pole (the exactness fix)
#endif
#ifndef PS5_OUTERCACHE
#define PS5_OUTERCACHE 1         // epoch-memoized outerNode + cross-blossom slackZ fast path
#endif
#endif // PS5_NO_BEST_STACK
// ---------------------------------------------------------------------------------------------
#if defined(PS5_LEAF_NOVEC) && !defined(PS5_LEAF_SPLICE)
#define PS5_LEAF_SPLICE   // NOVEC (drop vector maintenance) requires the list readers
#endif
#if defined(PS5_LEAF_NOVEC) && defined(PS5_LEAF_SHADOW)
#error "PS5_LEAF_SHADOW (vector referees) cannot combine with PS5_LEAF_NOVEC (no vectors)"
#endif
#if (defined(PS5_LEAF_SPLICE) || defined(PS5_LEAF_SHADOW)) && !defined(PS5_LEAF_LIST)
#define PS5_LEAF_LIST     // the list-maintenance layer
#endif
#ifdef PS5_NESTED_INTACT
#define PS5_NESTED_RELINK   // INTACT-nested builds on the z=0 recursive relink + dissolve uf-fix
#endif

// EVENT-HORIZON DEFER: DEFAULT-ON (2026-07-05). QSTAT-measured on the scaling loop: 92% of event
// inserts were never popped (phases end within <=3 ticks; inserts spread to 100K+ ticks, far ones
// through the allocating multiset overflow). Defer events firing more than H ticks past the window
// base; keep ONE per-node marker (kind 3) at the deferred MIN fire time; popping a marker rescans
// the node's CURRENT top-level blossom with the window re-based at the marker time (the outer-walk
// is required: no-relabel shrink relies on absorbed members' events surviving -- n=4 rep=0 trace).
// Marker time <= every deferred fire time, so the clock can never overrun a deferred edge.
// Interleaved A/B (ring+3n W=1000): loop -29%/-35% at n=51200/102400, single-shot -6%, gates
// 594/594 + 348/348 exact. Opt out: -DPS5_NO_EVT_HORIZON; tune: -DPS5_EVT_HORIZON=<ticks>.
#if !defined(PS5_EVT_HORIZON) && !defined(PS5_NO_EVT_HORIZON)
#define PS5_EVT_HORIZON 8
#endif

// PERSISTENT CLOCK + FOREST: DEFAULT-ON (2026-07-06). Stage 0-2 of docs/H_Scope_PersistentClock.md,
// all gates exact (348+594, flagged and flagless). Parity at n<=102400 (2.35/6.03s vs 2.34/6.21)
// and 2.7x at n=204800 (30.2s vs 81.8s) -- the per-phase forest regrowth it skips is what explodes
// on giant-blossom boundaries at scale (the 204800 cliff: 1.95B slackZW calls, ~95-edge boundary
// scans x 2480 phases). Opt out: -DPS5_NO_PERSIST (flips both together; PS5_PERSIST_SHRINK stays
// opt-in -- correct but structurally slow, see 320a6c9).
#ifndef PS5_NO_PERSIST
#ifndef PS5_PERSIST_CLOCK
#define PS5_PERSIST_CLOCK 1
#endif
#ifndef PS5_PERSIST_FOREST
#define PS5_PERSIST_FOREST 1
#endif
#endif

namespace Matchbox {

#ifdef PS5_ZSHADOW
static long long gZShadowB = -1;   // watch a single blossom (or -1 = first divergence)
#endif
#ifdef MATEGUARD
// writer-tag (relink-orientation-bug method): abort at the FIRST vertex-slot write of a node id.
#define MG_CHK(site, a, b) do { if ((Sz)(a) >= mN || (Sz)(b) >= mN) { \
  std::fprintf(stderr,"[MATEGUARD %s a=%ld b=%ld mN=%lu]\n",(site),(long)(a),(long)(b),(unsigned long)mN); \
  std::fflush(stderr); std::abort(); } } while(0)
#else
#define MG_CHK(site, a, b)
#endif
#ifdef PS5_UF_CHECK
long long gUfCheckFired = 0;
#endif
#ifdef S5DBG
static long long gS5Shrinks = 0, gS5Dissolves = 0, gS5Pops = 0;
static long long gS5Searches = 0, gS5EdgeScans = 0, gS5DissolveRelabels = 0, gS5PhaseRelabels = 0;
static long long gS5VbSup = 0, gS5VbRepairs = 0, gS5VbRepairScans = 0;   // v5 attribution
static long long gS5PersistPhases = 0;   // phase boundaries that took the PERSISTENT path (no relabel/reseed)
static long long gS5Augs = 0;   // augments in the no-relabel multi-augment path (augments/phase probe)
static long long gS5DisScanWork = 0, gS5SyncWork = 0, gS5Iters = 0;
// CLIFF profiling (2026-07-02): the per-line cost inside shrink + relinkOddCycleS5 (callgrind said
// these two are ~64% of the run at n=12800 uniform). gS5RelinkLeafScan = leaf-array scan work in
// relink's entry/attach lookups (loops over mNodeLeaves of every cycle member). gS5RelinkChildScan =
// the "for c in all mNumNodes" child-search sweep. gS5ShrinkLeafCopy = shrink's mNodeLeaves[B] rebuild
// copy work (append every leaf of every cycle member). gS5RelinkCalls/gS5ShrinkCalls = invocation counts.
static long long gS5RelinkLeafScan = 0, gS5RelinkChildScan = 0, gS5ShrinkLeafCopy = 0;
static long long gS5RelinkCalls = 0, gS5ShrinkCalls = 0;
// distinguish "inherent giant blossom" from "redundant re-copy under nesting": gS5MaxBlossom = largest
// leaf-set of any blossom formed; gS5SubBlossomShrinks = shrinks whose cycle contains a sub-blossom (nesting).
static long long gS5MaxBlossom = 0, gS5SubBlossomShrinks = 0;
static long long gS5MaxDistinct = 0;   // distinct-vertex count of the biggest blossom at its peak (dup check)
static long long gS5InternalScans = 0;   // edges scanned that were INTERNAL (same outer blossom) = pure waste
static long long gS5BigNodeCalls = 0, gS5BigNodeLeafSum = 0, gS5BigNodeThresh = 1000;  // scheduleNodeEdges on >1000-leaf nodes
static long long gS5BigViaMarker = 0, gS5BigViaReseed = 0;  // big-node scans via horizon-marker pop vs full reseed
static long long gS5BigRedundant = 0, gS5BigNecessary = 0;  // big-node scans that were provably redundant (memoizable) vs necessary
static long long gS5BigCacheHit = 0;  // big-node scans served from the event cache (PS5_BIGNODE_CACHE)
static long long gS5IncrSkip = 0;     // STAGE A: big-node scans that SKIPPED re-emission (events still live)
static long long gS5ReplayTotal = 0, gS5ReplayWchanged = 0;  // prediction test: replayed events vs those whose w-label changed
// distinct-weight cliff probe: gS5SyncCalls = # lazySyncRoles invocations (each O(n)); gS5LazyBeginCalls = # phases.
static long long gS5SyncCalls = 0, gS5LazyBeginCalls = 0;
// dirty-set feasibility: gS5SyncChanged = # vertices that ACTUALLY changed role across all lazySyncRoles
// calls (the rest of the O(n) scan per call is wasted no-ops). If changed << syncWork, a dirty-set wins.
static long long gS5SyncChanged = 0;
static long long gS5SyncFull = 0, gS5SyncIncr = 0;   // dirty-set: # syncs that took the full-scan vs incremental path
// M==E scan-restart: how many times solve() fell back to the O(nm) pure-scan re-solve (the DISTINCT-weight cliff).
static long long gS5ScanRestarts = 0;
static long long gS5OuterSteps = 0, gS5OuterCalls = 0, gS5SlackZSteps = 0, gS5SlackZCalls = 0, gS5UfSteps = 0, gS5UfCalls = 0;
// REDUNDANCY profiling (answer: are the 2.49M schedule() calls mostly re-inserting the SAME
// pending event?). gS5SchedCalls = scheduleNodeEdges invocations; gS5SchedInserts = mPQ.schedule
// calls; gS5SchedRedundant = inserts where this directed edge (u,w) was ALREADY scheduled at the
// SAME tNow since the last structural relabel (provably no new information). gS5NodeReSched =
// scheduleNodeEdges calls on a node already scheduled at this tNow/epoch.
static long long gS5SchedCalls = 0, gS5SchedInserts = 0, gS5SchedRedundant = 0, gS5NodeReSched = 0;
#endif
#ifdef RELINKDBG
static bool gS5NrdSeen = false;   // set true after the first in-place (no-relabel) dissolution
#endif

class ScalingS5Engine {
public:
#ifdef PS5_IDX32
  // INDEX DIET (2026-07-12, the HPC enabler): 32-bit node/vertex indices halve every index
  // vector and improve cache density in the memory-bound walk/scan loops. Valid for
  // n < 2^31 (mNumNodes <= 2n+1 and NIL()=(Sz)-1=0xFFFFFFFF stay representable); init()
  // enforces the bound. Weights, times, tokens, and epochs remain 64-bit (LL) throughout --
  // tokens/epochs exceed 2^32 within a single run (measured 1.3B+ lcaOuter calls).
  typedef std::uint32_t Sz;
#else
  typedef std::size_t Sz;
#endif
  typedef long long LL;
  static Sz NIL() { return (Sz)-1; }
  enum { CRIT1 = 1, CRIT2 = 2, CRIT3 = 3 };   // eligibility criteria (public: driver sets via setCrit)

  ScalingS5Engine() : mLastBoundedCalls(0), mPhaseBudget(0), mLastPhaseCount(0), mTickBudget(0), mTicksUsed(0), mN(0), mNReal(0), mNumNodes(0), mMaxW(0), mUseEvent(false), mLazyActive(false), mLazyVerify(false), mLazyEvent(false), mNoRelabelShrink(false), mNoRelabelAug(false), mNoRelabelDissolve(false), mCrit(CRIT1), mIgnoreDummies(false), mAllowOvertight(false), mLcaTok(0)
    , mDbgEpoch(0)
#ifdef PS5_SFM
    , mSFMon(false), mSfmChecks(0), mSfmMismatch(0)
#endif
#ifdef PS5_CQ_EVENTS
    , mCqChecks(0), mCqMismatch(0)
#endif
#ifdef PS5_FLIP_DIRTY
    , mFlipAll(true)
#endif
#ifdef PS5_INVKEY
    , mIkTok(0)
#endif
#ifdef PS5_ANCESTOR_HINT
    , mAhTok(0)
#endif
#ifdef PS5_VTX_BEST
    , mVbEnabled(true)
#endif
#if defined(PS5_VTX_BEST) && defined(PS5_VB_DEFER_REPAIR)
    , mVbRepTok(0)
#endif
#ifdef PS5_DEEP_INTACT
    , mSEMemoTok(0)
#endif
#ifdef PS5_OUTERCACHE
    , mOuterEpoch(0)
#endif
  { mEventStuck = false; mSyncDirtyAll = true;
#ifdef PS5_VTX_BEST
    mVBestEpoch = 0;
#endif
#ifdef PS5_ZSHADOW
    mZShadowOk = false;
#endif
  }
#ifdef PS5_SFM
  Sz sfmChecks() const { return mSfmChecks; }
  Sz sfmMismatch() const { return mSfmMismatch; }
#endif
  void init(Sz n) {
#ifdef PS5_IDX32
    // the diet's validity bound: mNumNodes can reach 2n+1 and NIL()=0xFFFFFFFF must stay free
    if ((unsigned long long)n >= 0x7FFFFFFFULL) { std::fprintf(stderr,"PS5_IDX32: n=%llu exceeds 2^31-1\n",(unsigned long long)n); std::abort(); }
#endif
    mN = n;
    mNReal = n;                       // default: all vertices real (bridge overrides via setReals)
    mNumNodes = n;                    // was UNINITIALIZED until the first reset/solve wrote it; any
                                      // earlier reader (the liquidation fold was the first) saw
                                      // garbage -- the n=130 "hang" was a loop over ~1.5e19 nodes.
    mAdjV.assign(n, std::vector<Sz>());
    mAdjW.assign(n, std::vector<LL>());
    mMaxW = 0;
#ifdef PS5_QSTAT
    qstatReset();
#endif
#ifdef PS5_STAGE0
    s0Reset();
#endif
    // NOTE: do NOT reset mUseEvent here -- callers may set it before init().
  }
  void addEdge(Sz u, Sz v, LL w) {
    if (u == v) return;
    mAdjV[u].push_back(v); mAdjW[u].push_back(w);
    mAdjV[v].push_back(u); mAdjW[v].push_back(w);
    if (w > mMaxW) mMaxW = w;
  }

  // Returns max-weight perfect matching weight; fills mate (0-based). -1 if none.
  LL solve(std::vector<Sz>* mateOut);
  // foreign-index overload (driver passes vector<Size>; under PS5_IDX32 that differs from Sz)
  template <class I>
  LL solve(std::vector<I>* mateOut) {
    std::vector<Sz> mm; LL g = solve(&mm);
    if (mateOut) { mateOut->assign(mm.size(), (I)-1);
      for (Sz v = 0; v < mm.size(); ++v) (*mateOut)[v] = (mm[v] == NIL()) ? (I)-1 : (I)mm[v]; }
    return g;
  }
  Sz mateOf(Sz v) const { return mMate[v]; }

  // ---- PERSISTENT-S5 DRIVER INTERFACE (D path Option 1) -----------------------------------------
  // The scaling loop (ScalingMatching::computePersistS5) drives ONE S5 across all scales, doing the
  // per-scale Scaling/Liquidate/Perfection ops DIRECTLY on S5's state between searchEvent calls
  // (S5's lazy duals are local to each search; between searches mY/mZ are materialized). These
  // accessors let the driver read/write that state. (S5 + the driver are coupled by design -- the
  // §5 doc says "build them together".)
  Sz nVerts() const { return mN; }
  Sz numNodes() const { return mNumNodes; }
  std::vector<LL>& y() { return mY; }              // per-vertex dual (size >= mN)
  std::vector<LL>& z() { return mZ; }              // per-node dual (blossoms >= mN)
  std::vector<Sz>& mate() { return mMate; }
  Sz lastBoundedCalls() const { return mLastBoundedCalls; }   // diagnostic for searchBounded
  // Phase budget: bound a searchEvent to at most `k` completed MV multi-augment phases (the tau
  // bound that makes per-scale Phase-A O(tau*m) and LEAVES free reals for park-on-dummies, instead
  // of one searchEvent draining the whole matching). 0 = unbounded (default; all other callers).
  void setPhaseBudget(Sz k) { mPhaseBudget = k; }
  Sz lastPhaseCount() const { return mLastPhaseCount; }       // diagnostic: phases the last searchEvent ran
  // TICK budget: bound the TOTAL dual adjustment (clock advance, summed across phases) a search may
  // perform before stopping -- the paper's tau semantics (Fig 5 step 7 = tau SearchOne calls = tau
  // unit dual adjustments per scale). 0 = unbounded. Unlike the phase budget (whose unit -- an MV
  // phase -- can contain MANY adjustments), this bounds exactly what the paper bounds.
  void setTickBudget(LL b) { mTickBudget = b; mTicksUsed = 0; }
  LL ticksUsed() const { return mTicksUsed; }                 // diagnostic
  std::vector<std::vector<LL> >& edgeW() { return mAdjW; }   // mAdjW[u][k] mirrors mAdjV[u][k]
  std::vector<std::vector<Sz> >& adjV() { return mAdjV; }
  // Initialize the search state for a fresh scale's FVR: empty matching, fresh forest (no blossoms).
  // Duals (mY/mZ) are KEPT (carried). Sets mNumNodes=mN. (The driver calls this when the loop empties
  // M per scale -- spec step 1.) ensureCap so later searches can grow blossom ids.
  void resetForestKeepDuals() {
    Sz cap = 2 * mN + 1; ensureCap(cap);
#ifndef PS5_NO_LIQ_FOLD
    // DEFAULT-ON (2026-07-05, 594-gated; opt out with -DPS5_NO_LIQ_FOLD, matching the driver's
    // z-doubling flag -- the two must flip together).
    // PAPER-FAITHFUL LIQUIDATION (Fig 5 steps 3+5): before dropping a blossom, fold its dual into
    // its members -- Liquidate(B): y(u) += z(B)/2 for all u in B. Dropping z WITHOUT the fold (the
    // old behavior below) strictly lowers yz on all edges inside former blossoms and violates
    // Property 2 (measured: the violation appears at the first blossom-forming scale and doubles
    // every scale thereafter -- SCALECHK, n=10 rep=4). The driver doubles BOTH y and z at each
    // scale start (spec step 2), so the fold here is uniformly z/2: equal to the spec's Liquidate
    // on the doubled z at scale boundaries, and to liquidate-now semantics at the Finalization
    // reset. Nested blossoms fold once per containing node via each node's own leaf list.
    { Sz hi = mNumNodes;                       // defensive clamp: never trust a count past the arrays
      if (hi > mZ.size()) hi = mZ.size();
      if (hi > mNodeLeaves.size()) hi = mNodeLeaves.size();
      for (Sz b = mN; b < hi; ++b) {
        if (mZ[b] != 0) {
#ifdef PS5_LEAF_NOVEC
          // NOVEC: vectors are empty for shrink-created blossoms -- fold via the SEGMENT walk
          // (works for interior nodes too: their (head,tail) delimits their members).
          if (b < mLfHead.size() && mLfHead[b] != NIL()) {
            Sz e = mLfTail[b];
            for (Sz u = mLfHead[b], nx = NIL(); u != NIL(); u = nx) {
              nx = (u == e) ? NIL() : mLfNext[u];
              if (u < mN) mY[u] += mZ[b] / 2;
            }
          }
#else
          for (Sz j = 0; j < mNodeLeaves[b].size(); ++j) {
            Sz u = mNodeLeaves[b][j];
            if (u < mN) mY[u] += mZ[b] / 2;
          }
#endif
        }
      }
    }
#endif
    mNumNodes = mN;
#ifdef PS5_ANCESTOR_HINT
    ++mAhTok; ahEnsure();   // node ids get reused from mN up -> all hints invalid
#endif
    mMate.assign(mN, NIL());
    mBParent.assign(cap, NIL()); mBBase.assign(cap, NIL());
    mBCycle.assign(cap, std::vector<Sz>()); mBCycEU.assign(cap, std::vector<Sz>()); mBCycEV.assign(cap, std::vector<Sz>());
    mNodeLeaves.assign(cap, std::vector<Sz>());
    mLabel.assign(cap, LBL_FREE); mParent.assign(cap, NIL()); mParEU.assign(cap, NIL()); mParEV.assign(cap, NIL()); mRoot.assign(cap, NIL());
    mUf.assign(mN, 0); mFreeBlossoms.clear();
    for (Sz v = 0; v < mN; ++v) { mBBase[v] = v; mNodeLeaves[v].assign(1, v); mUf[v] = v; }
#ifdef PS5_LEAF_LIST
    lfRebuildFromVectors();
#endif
    // z for the (now-gone) blossoms is stale; keep mY (vertex duals) which the scale carries.
    for (Sz b = mN; b < cap; ++b) mZ[b] = 0;
  }
  // Run the FVR search to a perfect matching on the CURRENT state (used as a scale's search). Reuses
  // searchToPerfect (already added). Returns false if no PM at this scale's duals/weights.

  // ---- D PATH (a) STATE BRIDGE ----------------------------------------------------------------
  // Load a per-scale state from the loop engine (ScalingMatching) so S5 can run the FVR search on
  // it with its sqrt(n) per-level batching, then read the result back. The loop and S5 share the
  // SAME blossom-forest model (mBParent/mBBase/mBCycle/mBCycEU(loop:mBCycEdgeU)/mNodeLeaves), so this
  // is a direct array copy + a per-scale weight refresh (S5 mAdjW = loop scale-weight / 2, since the
  // loop's mW is the doubled-even 2*w and S5 slackZW uses 2*mAdjW).
  // Caller passes the loop's current arrays. numNodes = loop mNumNodes (verts + live blossoms).
  // Templated on the caller's index element type (the driver uses Matchbox::Size; under
  // PS5_IDX32 the engine's Sz is narrower -- these bridges convert element-wise; transfers
  // happen once per scale, so the copy cost is irrelevant). NIL sentinels map (I)-1 <-> NIL().
  template <class I>
  void loadState(Sz numNodes,
                 const std::vector<I>& mate, const std::vector<LL>& y, const std::vector<LL>& z,
                 const std::vector<I>& bparent, const std::vector<I>& bbase,
                 const std::vector<std::vector<I> >& bcycle,
                 const std::vector<std::vector<I> >& bcycU, const std::vector<std::vector<I> >& bcycV,
                 const std::vector<std::vector<I> >& nodeLeaves,
                 const std::vector<I>& uf) {
    Sz cap = 2 * mN + 1; if (numNodes + 1 > cap) cap = numNodes + 1;
    ensureCap(cap);
    mNumNodes = numNodes;
    mMate.assign(mN, NIL());
    for (Sz v = 0; v < mN && v < mate.size(); ++v) mMate[v] = (mate[v] == (I)-1) ? NIL() : (Sz)mate[v];
    mY.assign(mN, 0); for (Sz v = 0; v < mN && v < y.size(); ++v) mY[v] = y[v];
    mZ.assign(cap, 0); for (Sz b = 0; b < numNodes && b < z.size(); ++b) mZ[b] = z[b];
    mBParent.assign(cap, NIL()); for (Sz b = 0; b < numNodes && b < bparent.size(); ++b) mBParent[b] = (bparent[b] == (I)-1) ? NIL() : (Sz)bparent[b];
    mBBase.assign(cap, NIL());   for (Sz b = 0; b < numNodes && b < bbase.size(); ++b) mBBase[b] = (bbase[b] == (I)-1) ? NIL() : (Sz)bbase[b];
    mBCycle.assign(cap, std::vector<Sz>());      for (Sz b = 0; b < numNodes && b < bcycle.size(); ++b) mBCycle[b].assign(bcycle[b].begin(), bcycle[b].end());
    mBCycEU.assign(cap, std::vector<Sz>());       for (Sz b = 0; b < numNodes && b < bcycU.size(); ++b) mBCycEU[b].assign(bcycU[b].begin(), bcycU[b].end());
    mBCycEV.assign(cap, std::vector<Sz>());       for (Sz b = 0; b < numNodes && b < bcycV.size(); ++b) mBCycEV[b].assign(bcycV[b].begin(), bcycV[b].end());
    mNodeLeaves.assign(cap, std::vector<Sz>());   for (Sz b = 0; b < numNodes && b < nodeLeaves.size(); ++b) mNodeLeaves[b].assign(nodeLeaves[b].begin(), nodeLeaves[b].end());
    mUf.assign(mN, 0); for (Sz v = 0; v < mN; ++v) mUf[v] = (v < uf.size()) ? (Sz)uf[v] : v;
    // forest labels are search-internal scratch -- relabel() (first searchEvent iter) rebuilds them.
    mLabel.assign(cap, LBL_FREE); mParent.assign(cap, NIL()); mParEU.assign(cap, NIL());
    mParEV.assign(cap, NIL()); mRoot.assign(cap, NIL());
    mFreeBlossoms.clear();
#ifdef PS5_LEAF_LIST
    lfRebuildFromVectors();   // AFTER mBParent/mNumNodes restored
#endif
  }
  // Refresh S5's per-scale edge weights from the loop's doubled-even scale weight (loop passes
  // 2*w as scaleW2[u][k] mirroring its mAdj[u][k]; S5 stores w = scaleW2/2). mAdjV order must match
  // the loop's mAdj order -- the bridge builds S5 adjacency in that exact order (see ScalingMatching).
  void setScaleWeights(const std::vector<std::vector<LL> >& scaleW2) {
    mMaxW = 0;
    for (Sz u = 0; u < mN && u < scaleW2.size(); ++u)
      for (Sz k = 0; k < mAdjW[u].size() && k < scaleW2[u].size(); ++k) {
        mAdjW[u][k] = scaleW2[u][k] / 2;
        if (mAdjW[u][k] > mMaxW) mMaxW = mAdjW[u][k];
      }
  }
  template <class I>
  void readMate(std::vector<I>* out) const {
    if (!out) return; out->assign(mMate.size(), (I)-1);
    for (Sz v = 0; v < mMate.size(); ++v) (*out)[v] = (mMate[v] == NIL()) ? (I)-1 : (I)mMate[v];
  }
  void readDuals(std::vector<LL>* y, std::vector<LL>* z) const { if (y) *y = mY; if (z) { z->assign(mNumNodes, 0); for (Sz b = 0; b < mNumNodes; ++b) (*z)[b] = mZ[b]; } }
  // Read back the BLOSSOM FOREST that S5 built during the search (blossoms formed in the FVR must
  // be carried into the next scale -- their z>0 is part of the carried dual state). Without this the
  // loop loses the blossoms S5 created (the sumZ-divergence bug).
  Sz readNumNodes() const { return mNumNodes; }
  // Dissolve all z==0 ROOT blossoms (any label) before readback -- the loop treats a node with z!=0
  // as live and z==0 as dead, but S5 may leave z==0 blossoms (formed in the search, never grew, or
  // EVEN-roots that ended at 0) with intact sub-structure. Carrying them bloats mNumNodes and leaves
  // vertices "inside" dead blossoms. Flatten them so their vertices become outer (clean carry).
  void dissolveZeroRoots() {
    bool any = true; Sz guard = 0;
    while (any && ++guard <= mNumNodes + 4) {
      any = false;
      for (Sz b = mN; b < mNumNodes; ++b)
        if (mBParent[b] == NIL() && mZ[b] == 0 && !mBCycle[b].empty()) { dissolve(b); any = true; }
    }
  }
#ifdef PS5_LEAF_SPLICE
  // Stage 3 move 1: the leaf-list state transfer (O(n) arrays; NO synthesized vectors).
  // Call readLeafLists after readForest; call loadLeafLists after loadState (it overwrites the
  // vector-rebuilt lists with the transferred truth -- under SHADOW the sweep then validates
  // the transfer itself against the vectors).
  template <class I>
  void readLeafLists(std::vector<I>* nx, std::vector<I>* pv,
                     std::vector<I>* hd, std::vector<I>* tl, std::vector<I>* cn) const {
    if (nx) { nx->assign(mLfNext.size(), (I)-1); for (Sz v = 0; v < mLfNext.size(); ++v) (*nx)[v] = (mLfNext[v]==NIL())?(I)-1:(I)mLfNext[v]; }
    if (pv) { pv->assign(mLfPrev.size(), (I)-1); for (Sz v = 0; v < mLfPrev.size(); ++v) (*pv)[v] = (mLfPrev[v]==NIL())?(I)-1:(I)mLfPrev[v]; }
    if (hd) { hd->assign(mNumNodes, (I)-1); for (Sz b = 0; b < mNumNodes && b < mLfHead.size(); ++b) (*hd)[b] = (mLfHead[b]==NIL())?(I)-1:(I)mLfHead[b]; }
    if (tl) { tl->assign(mNumNodes, (I)-1); for (Sz b = 0; b < mNumNodes && b < mLfTail.size(); ++b) (*tl)[b] = (mLfTail[b]==NIL())?(I)-1:(I)mLfTail[b]; }
    if (cn) { cn->assign(mNumNodes, 0);     for (Sz b = 0; b < mNumNodes && b < mLfCnt.size();  ++b) (*cn)[b] = (I)mLfCnt[b]; }
  }
  template <class I>
  void loadLeafLists(const std::vector<I>& nx, const std::vector<I>& pv,
                     const std::vector<I>& hd, const std::vector<I>& tl, const std::vector<I>& cn) {
    if (nx.empty() && hd.empty()) return;   // nothing transferred (first scale): keep the rebuild
    mLfNext.assign(mN, NIL()); for (Sz v = 0; v < mN && v < nx.size(); ++v) mLfNext[v] = (nx[v]==(I)-1)?NIL():(Sz)nx[v];
    mLfPrev.assign(mN, NIL()); for (Sz v = 0; v < mN && v < pv.size(); ++v) mLfPrev[v] = (pv[v]==(I)-1)?NIL():(Sz)pv[v];
    Sz cap = mNodeLeaves.size();
    mLfHead.assign(cap, NIL()); for (Sz b = 0; b < mNumNodes && b < hd.size(); ++b) mLfHead[b] = (hd[b]==(I)-1)?NIL():(Sz)hd[b];
    mLfTail.assign(cap, NIL()); for (Sz b = 0; b < mNumNodes && b < tl.size(); ++b) mLfTail[b] = (tl[b]==(I)-1)?NIL():(Sz)tl[b];
    mLfCnt.assign(cap, 0);      for (Sz b = 0; b < mNumNodes && b < cn.size(); ++b) mLfCnt[b] = (Sz)cn[b];
#ifdef PS5_LEAF_SHADOW
    lfSweepAll("loadLeafLists");
#endif
  }
#endif
  template <class I>
  void readForest(std::vector<I>* bparent, std::vector<I>* bbase,
                  std::vector<std::vector<I> >* bcycle,
                  std::vector<std::vector<I> >* bcU, std::vector<std::vector<I> >* bcV,
                  std::vector<std::vector<I> >* nodeLeaves, std::vector<I>* uf) const {
    if (bparent) { bparent->assign(mNumNodes, (I)-1); for (Sz b = 0; b < mNumNodes; ++b) (*bparent)[b] = (mBParent[b]==NIL())?(I)-1:(I)mBParent[b]; }
    if (bbase)   { bbase->assign(mNumNodes, (I)-1);   for (Sz b = 0; b < mNumNodes; ++b) (*bbase)[b] = (mBBase[b]==NIL())?(I)-1:(I)mBBase[b]; }
    if (bcycle)  { bcycle->assign(mNumNodes, std::vector<I>());     for (Sz b = 0; b < mNumNodes; ++b) (*bcycle)[b].assign(mBCycle[b].begin(), mBCycle[b].end()); }
    if (bcU)     { bcU->assign(mNumNodes, std::vector<I>());        for (Sz b = 0; b < mNumNodes; ++b) (*bcU)[b].assign(mBCycEU[b].begin(), mBCycEU[b].end()); }
    if (bcV)     { bcV->assign(mNumNodes, std::vector<I>());        for (Sz b = 0; b < mNumNodes; ++b) (*bcV)[b].assign(mBCycEV[b].begin(), mBCycEV[b].end()); }
    if (nodeLeaves) { nodeLeaves->assign(mNumNodes, std::vector<I>()); for (Sz b = 0; b < mNumNodes; ++b) (*nodeLeaves)[b].assign(mNodeLeaves[b].begin(), mNodeLeaves[b].end()); }
    if (uf)      { uf->assign(mN, 0); for (Sz v = 0; v < mN; ++v) (*uf)[v] = (I)mUf[v]; }
  }
  // Run searchEvent to a perfect matching on the ALREADY-LOADED state (does NOT re-init like solve).
  // For the bridge: after loadState + setScaleWeights, drive the FVR search. Returns false if no PM.
  bool searchToPerfect() {
    for (Sz guard = 0; guard <= 2 * mN + 2; ++guard) {
      bool perfect = true;
      for (Sz v = 0; v < mNReal; ++v) if (mMate[v] == NIL()) { perfect = false; break; }
      if (perfect) return true;
      if (!searchEvent()) return false;
    }
    for (Sz v = 0; v < mNReal; ++v) if (mMate[v] == NIL()) return false;
    return true;
  }
  // BOUNDED: run at most `maxCalls` searchEvent calls (the tau-bounded FVR -- do NOT drive to a real
  // perfect matching, leave reals parked for Finalization). Stops early if perfect.
  void searchBounded(Sz maxCalls) {
    mLastBoundedCalls = 0;
    for (Sz c = 0; c < maxCalls; ++c) {
      bool perfect = true;
      for (Sz v = 0; v < mNReal; ++v) if (mMate[v] == NIL()) { perfect = false; break; }
      if (perfect) return;
      ++mLastBoundedCalls;
      if (!searchEvent()) return;
    }
    // (mLazyActive is reset to false inside searchEvent on every augment-materialize, so mY/mZ are
    //  current after the loop. A call returning via !anyRoot also leaves the last materialize intact.)
  }

private:
  Sz mLastBoundedCalls;             // diagnostic: searchEvent calls made by the most recent searchBounded
  Sz mPhaseBudget;                  // 0 = unbounded; >0 = stop searchEvent after this many completed MV-phases (the tau bound)
  Sz mLastPhaseCount;               // diagnostic: MV-phases the most recent searchEvent completed
  LL mTickBudget;                   // 0 = unbounded; >0 = stop when total clock advance exceeds this (paper tau)
  LL mTicksUsed;                    // accumulated clock advance since setTickBudget
  Sz mN;
  std::vector<std::vector<Sz> > mAdjV;
  std::vector<std::vector<LL> > mAdjW;
  LL mMaxW;

  // ---- node model: 0..mN-1 vertices, mN..mNumNodes-1 blossom super-nodes ----
  Sz mNumNodes;
  std::vector<Sz> mMate;            // mMate[v] = matched vertex (real verts only) or NIL

  // duals (materialized; Stage 3.1). y per vertex, z per blossom node (even >=0).
  std::vector<LL> mY;               // size mN
  std::vector<LL> mZ;               // size mNumNodes (0 for vertices)

  // ---- UNION-FIND over OUTER blossoms (structure #1) ----
  // ufFind(v) returns the representative VERTEX of v's current outer blossom.
  // We keep the union-find over VERTICES; a blossom's representative vertex is its
  // base. ufParent/ufRank are the standard DSU. ufReset re-points all to self.
  std::vector<Sz> mUf;              // union-find parent (over vertices 0..mN-1)
#ifdef DISS_TRACE
  std::vector<Sz> gDissKids;        // last dissolved blossom's cycle (for tracing)
#endif
#ifdef SHADOW_DUAL
  std::vector<LL> mShadowY;         // parallel materialized y (by outer label) to catch lazy drift
#endif
#ifdef PS5_ZSHADOW
  // PERSISTENT materialized shadow z, updated by the EXACT per-step rule (+-2*delta for ROOT blossoms) and
  // NEVER re-seeded across lazyBegin/dissolution. Seeded once (mZShadowOk). At each dual advance we bump the
  // root-blossom z and assert zLazy(b)==mZShadow[b] -> the FIRST divergence pinpoints the drift. This is the
  // paper-faithful z; if it stays == zLazy everywhere, zLazy is right and the bug is elsewhere.
  std::vector<LL> mZShadow; bool mZShadowOk;
  std::vector<LL> mYShadow;   // persistent materialized y (never re-seeded across lazyBegin) -- catches y drift
  void zShadowEnsure(Sz b) { if (b >= mZShadow.size()) mZShadow.resize(b + 16, 0); }
#endif
#ifdef PS5_OUTERCACHE
  // outerNode(v) cache: valid iff mOuterEpochOf[v]==mOuterEpoch. mOuterEpoch bumps at every shrink +
  // dissolve (the only mBParent mutators), invalidating all stale entries in O(1). mutable: outerNode
  // is const but memoizes.
  mutable std::vector<Sz> mOuterCache;
  mutable std::vector<LL> mOuterEpochOf;
  mutable LL mOuterEpoch;
  void outerCacheBump() { ++mOuterEpoch; }
  void outerCacheEnsure() {
    if (mOuterCache.size() < mN) { mOuterCache.assign(mN, 0); mOuterEpochOf.assign(mN, -1); }
  }
#endif
#ifdef PS5_DEEP_INTACT
  // memo for relinkStructEligible: per-blossom result valid within one dissolution batch (mSEMemoTok).
  // 0 = unknown, 1 = eligible, 2 = ineligible. Avoids the O(depth^2) re-walk (same subtree re-checked
  // at every ancestor). Keyed by blossom id; mSEMemoStamp[B]==mSEMemoTok means mSEMemo[B] is current.
  std::vector<char> mSEMemo; std::vector<LL> mSEMemoStamp; LL mSEMemoTok;
#endif
#ifdef PS5_DISHEAP
  // DISSOLUTION MIN-HEAP (replaces the O(iters x blossoms) per-iteration scan for the next dissolution).
  // A blossom's dissolution time tNow+z/2 is INVARIANT while it stays ODD-root (as tNow advances by d,
  // zLazy drops by 2d). So we key by that stable time. Entries go STALE on role change -> lazy-validate
  // on pop (discard if the blossom is no longer an ODD-root at the stored time). Push when a blossom
  // becomes ODD-root; the rare relabel-fallback rebuilds the heap by a one-time scan.
  std::vector<std::pair<LL,Sz> > mDisHeap;   // (dissolutionTime, blossomId) min-heap on .first
  void disHeapPush(Sz b, LL t) {
    mDisHeap.push_back(std::make_pair(t, b));
    Sz i = mDisHeap.size() - 1;
    while (i > 0) { Sz p = (i - 1) / 2; if (mDisHeap[p].first <= mDisHeap[i].first) break;
      std::swap(mDisHeap[p], mDisHeap[i]); i = p; }
  }
  void disHeapPopMin() {
    Sz n = mDisHeap.size(); if (n == 0) return;
    mDisHeap[0] = mDisHeap[n - 1]; mDisHeap.pop_back(); --n;
    Sz i = 0; for (;;) { Sz l = 2*i+1, r = 2*i+2, s = i;
      if (l < n && mDisHeap[l].first < mDisHeap[s].first) s = l;
      if (r < n && mDisHeap[r].first < mDisHeap[s].first) s = r;
      if (s == i) break; std::swap(mDisHeap[s], mDisHeap[i]); i = s; }
  }
  // schedule b's dissolution (call when b becomes an ODD root). tNow supplied by caller.
  void disHeapSchedule(Sz b, LL tNow) {
    if (b < mN) return;
    LL zb = mLazyActive ? zLazy(b) : mZ[b];
    disHeapPush(b, tNow + zb / 2);
  }
  // rebuild from scratch (rare: the relabel fallback). O(#odd-root blossoms).
  void disHeapRebuild(LL tNow) {
    mDisHeap.clear();
    for (Sz b = mN; b < mNumNodes; ++b)
      if (mBParent[b] == NIL() && mLabel[b] == LBL_ODD) disHeapSchedule(b, tNow);
  }
  // pop the earliest VALID dissolution: returns blossom id + sets *outT, or NIL if none.
  // A blossom that is still an ODD root is a valid dissolution candidate; if its dissolution time DRIFTED
  // (its z changed via a different role interval), the stored key is stale -> re-push at the CORRECT time
  // and re-heapify (NOT discard -- discarding drops real dissolutions -> relabel-fallback storm). Only
  // truly-dead entries (b no longer an ODD root) are discarded.
  Sz disHeapNext(LL tNow, LL* outT) {
    while (!mDisHeap.empty()) {
      LL t = mDisHeap[0].first; Sz b = mDisHeap[0].second;
      if (b >= mNumNodes || mBParent[b] != NIL() || mLabel[b] != LBL_ODD) { disHeapPopMin(); continue; }
      LL cur = tNow + (mLazyActive ? zLazy(b) : mZ[b]) / 2;
      if (cur == t) { *outT = t; return b; }         // key current -> this is the min valid dissolution
      disHeapPopMin(); disHeapPush(b, cur);          // key drifted -> re-key and re-sift
    }
    return NIL();
  }
#endif
#ifdef PS5_SFM
  // ConcatenableQueue-backed blossom MEMBERSHIP (the cliff fix): a blossom = a contiguous sublist of
  // member vertices. merge() on shrink, split() on dissolve -> structurally-correct O(log n) nested
  // dissolution (the cycle-walk relinkOddCycleS5 can't do nested; this can). key(u)=min incident slack
  // from outside (split-findmin findmin) gives the next-eligible edge for freed-region scheduling.
  SplitFindmin mSFM;
  bool mSFMon;                      // SFM consistency/maintenance active for the current search
#endif
#ifdef PS5_CQ_EVENTS
  // K_Scope Stage 1: per-EVEN-blossom event queue keyed by fire-time. Stage 1a = REBUILT each reseed
  // + refereed against sfmBruteNextFire (validates the per-blossom decomposition + CQ query). Stage 1b
  // will maintain it incrementally via merge (shrink) / split (dissolve). NOT driving yet.
  ConcatenableQueue mCQ;
  Sz mCqChecks; Sz mCqMismatch;
  // build every EVEN top-level blossom's boundary-edge fire-times into one CQ queue per blossom, then
  // assert the global min matches the brute next-fire. tNow-relative; O(m log n) (shadow cost only).
  void cqCheckAfterReseed(LL tNow);
  Sz cqChecks() const { return mCqChecks; }
  Sz cqMismatch() const { return mCqMismatch; }
  // STAGE 1b: LIVE per-vertex membership shadow via SplitFindmin (mSFM). One sublist per top-level node;
  // split at dissolve, re-init at shrink boundaries. Referee: sameList(u,v) iff outerNode(u)==outerNode(v).
  // This isolates the merge/split STRUCTURAL correctness (separate from keying) that Stage 2 depends on.
  std::vector<Sz> mCqOrder;         // last init's vertex order (for split pivots)
  std::vector<std::vector<Sz> > mCqPivots;   // mCqPivots[B] = the split pivot leaves recorded when B merged
                                             // (last leaf of each sub-run, in MERGE order) -- replayed at
                                             // dissolve, immune to rotateToBase reordering mBCycle.
  void cqInitMembership();          // init mSFM: vertices grouped by current top-level node, contiguous
  void cqCheckMembership(const char* where);   // assert sameList == same-outer-node for sampled pairs
#endif
#ifdef PS5_UF_CHECK
  // assert every vertex's uf-rep == its outer node's base (uf and mBParent agree). Fires at FIRST violation
  // with a tag naming the operation that just created it. gUfCheckFired stops after the first.
  void ufCheck(const char* where) {
    extern long long gUfCheckFired; if (gUfCheckFired) return;
    for (Sz v = 0; v < mN; ++v) {
      Sz on = outerNode(v); Sz wantRep = (on < mN) ? on : mBBase[on];
      if (ufFind(v) != ufFind(wantRep)) {
        std::fprintf(stderr,"[UF-DESYNC @%s v=%lu ufFind(v)=%lu outerNode(v)=%lu base=%lu ufFind(base)=%lu]\n",
          where,(unsigned long)v,(unsigned long)ufFind(v),(unsigned long)on,(unsigned long)wantRep,(unsigned long)ufFind(wantRep));
        gUfCheckFired = 1; return;
      }
    }
  }
#endif
#ifdef PS5_UF_RESYNC
  // rebuild union-find from the AUTHORITATIVE mBParent tree: every leaf of an outer node points to that
  // node's base. Fixes the uf<->mBParent desync (a nested dissolution leaves absorbed leaves with stale
  // separate reps -> the search sees an INTERNAL blossom edge as a live cross-blossom formation -> overrun
  // + suboptimal, e.g. n=100 W50 s3 edge (36,28): outerNode says same blossom, ufFind says different).
  void ufResyncFromTree() {
    for (Sz v = 0; v < mN; ++v) mUf[v] = v;
    for (Sz b = mN; b < mNumNodes; ++b) if (mBParent[b] == NIL()) {
      Sz rep = mBBase[b];
      std::vector<Sz>& lv = mNodeLeaves[b];
      for (Sz j = 0; j < lv.size(); ++j) mUf[lv[j]] = rep;
      mUf[rep] = rep;
    }
  }
#endif
  Sz ufFind(Sz v) {
#ifdef S5DBG
    ++gS5UfCalls;
#endif
    while (mUf[v] != v) {
#ifdef S5DBG
      ++gS5UfSteps;
#endif
      mUf[v] = mUf[mUf[v]]; v = mUf[v]; } return v; }
  void ufUnion(Sz a, Sz b, Sz newRep) {
    // point both roots at newRep (newRep must be one of the two roots' vertices).
    Sz ra = ufFind(a), rb = ufFind(b);
    mUf[ra] = newRep; mUf[rb] = newRep; mUf[newRep] = newRep;
  }

  // ---- blossom forest (for expansion/augment bookkeeping) ----
  // Each outer blossom super-node B (>= mN) records its odd cycle of CHILD nodes
  // and the edges linking them, exactly like the reference engine, so we can
  // expand (dissolve) and rotate the base for augmentation.
  std::vector<Sz> mBParent;         // immediate enclosing blossom node, or NIL
  std::vector<Sz> mBBase;           // base vertex of node (vertex id for a vertex)
  std::vector<std::vector<Sz> > mBCycle;   // child nodes in cycle order (base first)
  std::vector<std::vector<Sz> > mBCycEU;   // linking edge endpoint (in child i)
  std::vector<std::vector<Sz> > mBCycEV;   // linking edge endpoint (in child i+1)
  std::vector<std::vector<Sz> > mNodeLeaves; // original vertices under a node

  // ---- alternating forest over OUTER nodes ----
  enum { LBL_FREE = 0, LBL_EVEN = 1, LBL_ODD = 2 };
  std::vector<int> mLabel;          // per super-node (indexed by representative-as-node)
  // STRUCTURE EPOCH: bumps on every relabel / dissolution / shrink-relink (every event-affecting
  // structural change). Used by S5DBG redundancy stats AND by the PS5_BIGNODE_CACHE validity key.
  // Unconditional (cheap counter) so the cache invalidates correctly even without S5DBG.
  long long mDbgEpoch;
#ifdef S5DBG
  // redundancy detector: per vertex, the (epoch, tNow) at which its edges were last scheduled.
  std::vector<long long> mDbgSchedEpoch, mDbgSchedT;
  std::vector<long long> mBigLastEpoch; std::vector<Sz> mBigLastLeaves;  // memoizability probe (big-node redundancy)
  // N_Scope Stage 0: per-node reseed-epoch stamp -- mRsEp[i]==e means node i was top-level EVEN
  // at reseed #e. A node stamped e-1 at reseed #e is a REPEAT (its events likely survived the
  // boundary; label-only proxy, see docs/N_Scope_InvariantKeys.md).
  std::vector<long long> mRsEp;
#endif
#ifdef PS5_BIGNODE_CACHE
  struct BcEvent { long long u, w; int kind; LL off; int wlbl; };   // cached event (fire=tNow+off); wlbl = w's outer label when cached
  std::vector<long long> mBcEpoch;  std::vector<Sz> mBcLeaves;  // per-node cache validity key (epoch, leafcount)
  std::vector<std::vector<BcEvent> > mBcEvents;               // per-node cached event list
#ifdef PS5_INCR_EVENTS
  // STAGE A (docs/I_Scope_IncrementalBigBlossomEvents.md): inverse index w -> the (bigNode E, event
  // index) pairs referencing external node w. Built O(1) as events are cached; on a marker-pop HIT
  // with the queue generation UNCHANGED (events still live), re-emit ONLY the events whose w was
  // dirtied since the cache build, skipping the 99.96% unchanged. mWref[w] size <= deg(w) (bounded).
  struct WrefEntry { Sz node; Sz idx; };                     // an event: mBcEvents[node][idx]
  std::vector<std::vector<WrefEntry> > mWref;                // per external node -> its events
  std::vector<long long> mBcGen;                             // queue generation when node E's cache built
#endif
#endif
#ifdef PS5_MARKER_ADAPT
  std::vector<unsigned> mMarkerPops;   // per-node marker-pop count (giant-blossom rescan-thrash detector)
  std::vector<char> mNoDeferNode;      // per-node: once marker-thrashing, stop deferring this node
#endif
  std::vector<Sz>  mParent;         // parent OUTER node in the tree
  std::vector<Sz>  mParEU, mParEV;  // the (childVtx, parentVtx) edge that grew it
  std::vector<Sz>  mRoot;           // tree root super-node

  // map a representative VERTEX to its outer NODE id and back. For a bare vertex
  // the node id == the vertex; for a blossom the node id is the super-node whose
  // base's uf-rep is this vertex. We store nodeOfRep so labels live on node ids.
  // Simplest faithful approach: labels/parent/root are indexed by NODE id; to get
  // the outer NODE of vertex v we walk mBParent to the top (still O(depth) but
  // only used where structure is stable). ufFind gives same-blossom test in O(a).
  Sz outerNode(Sz v) const {
#ifdef PS5_OUTERCACHE
    // cache outerNode(v), invalidated per structural epoch (bumped at every shrink/dissolve, the only
    // mBParent mutators). Turns the O(depth) mBParent-walk into O(1) amortized for repeat queries in the
    // same epoch -- the residual-cliff fix (slackZW/yLazy/zLazy/outerNode all share this deep walk).
    if (mOuterEpochOf.size() > v && mOuterEpochOf[v] == mOuterEpoch) {
#ifdef S5DBG
      ++gS5OuterCalls;   // count the query (0 steps -- cache hit)
#endif
      return mOuterCache[v];
    }
#endif
    Sz x = v;
#ifdef PS5_ANCESTOR_HINT
    // O_Scope: jump to the last-known top-level ancestor and continue upward. Any live hinted
    // ancestor is a valid waypoint (shrink only adds parents ABOVE old tops; dissolution is the
    // only membership remover and it marks the id dead; ids are only reused across forest
    // resets, which bump mAhTok). Write-back below = path-compression semantics kept OUTSIDE
    // mBParent (which stays exact for dissolve/rotate).
    if (v < mAhStamp.size() && mAhStamp[v] == mAhTok) {
      Sz h = mAhHint[v];
      if (h != NIL() && h < mAhDead.size() && !mAhDead[h]) {
#ifdef PS5_AH_SHADOW
        // referee: the hint must be v's ancestor (walk from v must pass through h)
        { Sz c = v; bool seen = (c == h);
          while (!seen && mBParent[c] != NIL()) { c = mBParent[c]; if (c == h) seen = true; }
          if (!seen) { std::fprintf(stderr,"[AH-NOT-ANCESTOR v=%lu h=%lu]\n",(unsigned long)v,(unsigned long)h); std::fflush(stderr); std::abort(); } }
#endif
        x = h;
      }
    }
#endif
#ifdef S5DBG
    Sz g = 0; ++gS5OuterCalls;
#endif
    while (mBParent[x] != NIL()) { x = mBParent[x];
#ifdef S5DBG
      ++gS5OuterSteps;
      if (++g > mNumNodes + 5) { std::fprintf(stderr,"  [outerNode CYCLE from %lu at %lu]\n",(unsigned long)v,(unsigned long)x); std::fflush(stderr); return x; }
#endif
    }
#ifdef PS5_OUTERCACHE
    if (mOuterEpochOf.size() > v) { mOuterCache[v] = x; mOuterEpochOf[v] = mOuterEpoch; }
#endif
#ifdef PS5_ANCESTOR_HINT
    if (x != v && v < mAhStamp.size()) { mAhHint[v] = x; mAhStamp[v] = mAhTok; }
#endif
    return x;
  }

  // is vertex v a leaf of blossom node B? (ungated -- used by the intact-rotate external-match fix)
  bool leafOfNode(Sz v, Sz B) const {
    if (B < mN) return v == B;
#ifdef PS5_LEAF_SPLICE
    // Paper semantics (Gabow #2.3/Appx C): membership = ancestry in the laminar tree, not a list
    // scan. Also stale-proof: the append-only vector can list vertices that have LEFT B; the
    // ancestry walk answers CURRENT membership. O(depth) vs O(size).
    { Sz x = v, g = 0;
      while (x != NIL() && g++ <= mNumNodes) { if (x == B) return true; x = mBParent[x]; }
      return false; }
#else
    for (Sz j = 0; j < mNodeLeaves[B].size(); ++j) if (mNodeLeaves[B][j] == v) return true;
    return false;
#endif
  }
#if defined(MATE_AUDIT) || defined(DIVDBG)
  // is vertex v a leaf of blossom node B?
  bool leafOf(Sz v, Sz B) const {
    if (B < mN) return v == B;
    for (Sz j = 0; j < mNodeLeaves[B].size(); ++j) if (mNodeLeaves[B][j] == v) return true;
    return false;
  }
  // returns # of asymmetric mates (mMate[mMate[v]] != v), sets *firstBad. 0 = clean.
  Sz mateAsym(Sz* firstBad) const {
    Sz asym=0; if(firstBad) *firstBad=NIL();
    for(Sz v=0;v<mN;++v){Sz mv=mMate[v]; if(mv!=NIL()&&(mv>=mN||mMate[mv]!=v)){++asym; if(firstBad&&*firstBad==NIL())*firstBad=v;}}
    return asym;
  }
#endif
#ifdef DIVDBG
  // DIAGNOSTIC (L_Scope livelock): verify every live blossom's cycle-edge arrays are ALIGNED --
  // eu[i] must lift (via mBParent) to cyc[i], ev[i] to cyc[(i+1)%k]. A misalignment makes a later
  // rotateToBase silently no-op (ROTATE-MISS) -> stale mBBase -> dead tree looks live -> livelock.
  // Returns the first offending blossom or NIL.
  mutable std::vector<char> mParWriter;   // DIVDBG: last writer site of mParent[x]
  void parw(Sz x, char t) const { if (x == NIL()) return; if (mParWriter.size() <= x) mParWriter.resize(x + 1, '?'); mParWriter[x] = t; }
  Sz findCycEdgeMisalign(Sz* badIdx) const {
    for (Sz B = mN; B < mNumNodes; ++B) {
      const std::vector<Sz>& cyc = mBCycle[B];
      if (cyc.empty()) continue;
      // only check LIVE blossoms: B must still own its children (mBParent[child]==B)
      bool live = true;
      for (Sz i = 0; i < cyc.size() && live; ++i) if (mBParent[cyc[i]] != B) live = false;
      if (!live) continue;
      Sz k = cyc.size();
      for (Sz i = 0; i < k; ++i) {
        Sz a = mBCycEU[B][i], b = mBCycEV[B][i];
        Sz ca = a; Sz g = 0;
        while (ca != NIL() && mBParent[ca] != B && g++ <= mNumNodes) ca = mBParent[ca];
        Sz cb = b; g = 0;
        while (cb != NIL() && mBParent[cb] != B && g++ <= mNumNodes) cb = mBParent[cb];
        if (ca != cyc[i] || cb != cyc[(i + 1) % k]) { if (badIdx) *badIdx = i; return B; }
      }
    }
    return NIL();
  }
  // DIAGNOSTIC: does any outer node's treeParent chain form a CYCLE? Returns the first offending node or
  // NIL. Used to pinpoint which structural op (shrink/dissolve/relink) corrupts the tree.
  Sz findTreeCycle() const {
    for (Sz s = 0; s < mNumNodes; ++s) {
      if (mBParent[s] != NIL()) continue;            // only outer nodes are tree nodes
      Sz x = s, steps = 0;
      while (x != NIL()) { Sz p = treeParent(x); if (p == NIL()) break; x = p;
        if (++steps > mNumNodes + 2) return s; }
    }
    return NIL();
  }
#endif
  // CONTRACTED-GRAPH tree parent (§5.1): the effective tree parent of outer node
  // X is the OUTER node of X's stored mParent. This resolves through any blossom
  // that has since absorbed the stored parent (or X itself), so the tree stays
  // consistent under shrink WITHOUT a relabel -- X's children automatically see B
  // as their parent once their stored parent is absorbed into B. Returns NIL at a
  // root, or if the resolved parent is X itself (absorbed into the same blossom).
  Sz treeParent(Sz X) const {
    Sz p = mParent[X];
    if (p == NIL()) return NIL();
#ifdef OUTERTAG
    { extern long long gOTtreeParent; gOTtreeParent -= gS5OuterCalls; }
#endif
    // NOTE (2026-07-10, n=176 livelock forensics): resolving the parent via the attach VERTEX
    // (outerNode(mParEV[X])) instead of the node id was tried and REVERTED -- it fixes the
    // stale-node case but desyncs from the label/mate threading elsewhere (ASan: the flip walk
    // reaches a root mid-path, mParEU=NIL -> mMate[NIL] OOB, n=70). The stale node id must be
    // repaired at the SOURCE (the dissolution child-rewire), not reinterpreted here.
    Sz po = outerNode(p);
#ifdef OUTERTAG
    { extern long long gOTtreeParent; gOTtreeParent += gS5OuterCalls; }
#endif
    return (po == X) ? NIL() : po;
  }

  LL w2max(Sz u, Sz v) const;
  LL slackZ(Sz u, Sz v) const;      // yz(e) - 2*w(e)
  LL slackZW(Sz u, Sz v, LL w2) const;  // same, 2w already in hand (no w2max rescan)

  void ensureCap(Sz needNode);
  Sz  newBlossom();
  void shrink(Sz E, Sz W, Sz u, Sz w);
  void rotateToBase(Sz B, Sz baseVtx);
  void flipUpToRoot(Sz N);
  void augment(Sz E, Sz W, Sz u, Sz w);
  void dissolve(Sz B);
  Sz  lcaOuter(Sz a, Sz b);
  void relabel();                   // reset labels, root free outer nodes EVEN
  bool searchOnce();                // one EdmondsSearch (scan-driven, Stage 3.1)

  // ---- Stage 3.2: event-driven search (BucketPQ) ----
  BucketPQ mPQ;
  bool treeLive(Sz N) const;        // tree of N still has an unmatched root base
  void scheduleNodeEdges(Sz E, LL tNow, LL hBase = -1);  // schedule events for EVEN node E's edges
  void reseedAllEven(LL tNow, int site);  // full-EVEN reseed funnel (site = attribution tag, L_Scope Stage 0)
#ifdef PS5_LEAF_LIST
  // M_Scope: the intrusive leaf list (maintenance layer; referees only under SHADOW).
  // Segments are keyed per child node by (head,tail): shrink = O(k) splice, dissolve = O(k) cut.
  std::vector<Sz> mLfNext, mLfPrev;          // per vertex
  std::vector<Sz> mLfHead, mLfTail, mLfCnt;  // per node
  void lfRebuildFromVectors() {              // init/restore: chain each TOP-LEVEL node's vector;
    Sz cap = mNodeLeaves.size();             // interior nodes get front/back/size (their vectors are
    mLfNext.assign(mN, NIL()); mLfPrev.assign(mN, NIL());   // contiguous subranges of the parent's).
    mLfHead.assign(cap, NIL()); mLfTail.assign(cap, NIL()); mLfCnt.assign(cap, 0);
    // The append-only vectors can contain DUPLICATES on carried state (the persist loop's
    // mid-scale restore) -- chaining a duplicate writes its next-pointer twice = a list CYCLE
    // (the n=28 retire spin). DEDUP on first occurrence; the list is the TRUE member sequence.
    std::vector<char> seen(mN, 0);
    for (Sz b = 0; b < mNumNodes; ++b) {
      if (mBParent[b] != NIL()) continue;                       // pass 1: TOP-LEVEL chains
      const std::vector<Sz>& lv = mNodeLeaves[b];
      if (lv.empty()) { if (b < mN) { mLfHead[b]=b; mLfTail[b]=b; mLfCnt[b]=1; } continue; }
      Sz h = NIL(), t = NIL(), cnt = 0;
      for (Sz i = 0; i < lv.size(); ++i) { Sz v = lv[i];
        if (v >= mN || seen[v]) continue; seen[v] = 1;
        if (h == NIL()) { h = v; t = v; }
        else { mLfNext[t] = v; mLfPrev[v] = t; t = v; }
        ++cnt;
      }
      mLfHead[b]=h; mLfTail[b]=t; mLfCnt[b]=cnt;
    }
    for (Sz b = mN; b < mNumNodes; ++b) {
      if (mBParent[b] == NIL()) continue;                       // pass 2: INTERIOR segment keys
      const std::vector<Sz>& lv = mNodeLeaves[b];
      if (lv.empty()) continue;
      // first/last DISTINCT-FIRST-OCCURRENCE members within the parent chain: scan for the first
      // and last v that are actually ON the chain (their links were set in pass 1).
      Sz h = NIL(), t = NIL(), cnt = 0;
      std::vector<char> mine(mN, 0);
      for (Sz i = 0; i < lv.size(); ++i) if (lv[i] < mN) mine[lv[i]] = 1;
      // walk the TOP-LEVEL owner's chain once and mark boundaries (O(owner size))
      Sz top = b; while (top < mBParent.size() && mBParent[top] != NIL()) top = mBParent[top];
      for (Sz v = mLfHead[top]; v != NIL(); v = (v == mLfTail[top] ? NIL() : mLfNext[v])) {
        if (mine[v]) { if (h == NIL()) h = v; t = v; ++cnt; }
      }
      mLfHead[b]=h; mLfTail[b]=t; mLfCnt[b]=cnt;
    }
    for (Sz b = 0; b < mNumNodes; ++b)
      if (mBParent[b] == NIL()) lfCheckNode(b, "rebuild");
  }
  // Uniform leaf iteration for converted readers: binds v; continue/break work; one brace opened.
#define MBX_LF_FOR(v, B) for (Sz v = ((B) < mLfHead.size() ? mLfHead[(B)] : NIL()), v##Nx = NIL(), v##End = ((B) < mLfTail.size() ? mLfTail[(B)] : NIL()); v != NIL(); v = v##Nx) if ((v##Nx = (v == v##End ? NIL() : mLfNext[v])), true)
  void lfSweepAll(const char* tag) {
#ifndef PS5_LEAF_SHADOW
    (void)tag; return;   // vector-comparing referees are SHADOW-only
#endif
    for (Sz b = 0; b < mNumNodes; ++b) {
      if (mBParent[b] != NIL()) continue;
      // skip DEAD blossom ids (dissolved: keys nulled; their stale vectors are not a live contract)
      if (b >= mN && (b >= mLfHead.size() || mLfHead[b] == NIL())) continue;
      lfCheckNode(b, tag);
    }
  }
  void lfCheckNode(Sz B, const char* tag) {  // referee: list contents == vector contents (multiset)
#ifndef PS5_LEAF_SHADOW
    (void)B; (void)tag; return;
#endif
    if (B >= mLfHead.size()) return;
    std::vector<Sz> a; Sz g = 0;
    for (Sz v = mLfHead[B]; v != NIL() && g++ <= mN; v = (v == mLfTail[B] ? NIL() : mLfNext[v])) a.push_back(v);
    std::vector<Sz> c = mNodeLeaves[B];
    std::sort(a.begin(), a.end()); std::sort(c.begin(), c.end());
    c.erase(std::unique(c.begin(), c.end()), c.end());   // vectors are append-only stale multisets
    if (a != c) {
      std::fprintf(stderr, "[LEAF-SHADOW MISMATCH %s B=%lu listN=%lu vecN=%lu]\n",
        tag, (unsigned long)B, (unsigned long)a.size(), (unsigned long)c.size());
      std::abort();
    }
  }
#endif
#ifdef PS5_VTX_BEST
  // *** MEASURED NEGATIVE (2026-07-10) -- KEEP OFF. THE MINIMAL CHERRY-PICK IS UNSOUND. ***
  // Hypothesis: suppress any insert that does not improve the target vertex's best scheduled fire
  // time (Gabow Appx D rule 1: one queue entry per vertex), relying on validate-at-pop + the
  // overrun-guard to rescue suppressed edges whose predecessor went stale. RESULT: 17/348, 71/594.
  // WHY: in THIS engine the best edge's SOURCE goes stale constantly (labels flip, trees are
  // consumed mid-phase, blossoms dissolve), so suppressed second-best edges are MISSED far too
  // often to repair. Gabow can hold one entry per vertex ONLY because his keys are dual-adjustment-
  // INVARIANT and outer nodes never leave S: Appendix D's rules are a PACKAGE -- per-vertex
  // aggregation (rule 1) is unsound without the invariant key frame (rule 2) and once-outer-
  // always-outer. Stage 1 must build rules 1+2 TOGETHER (frozen-Y keys + SLACK(v) with a real
  // invalidation discipline), not cherry-pick. Kept as the recorded negative.
  std::vector<LL> mVBestFire; std::vector<LL> mVBestEp; std::vector<char> mVBestCls; LL mVBestEpoch;
  // v3: suppression is keyed by (epoch, TARGET CLASS): a FREE target counts down at 1/tick, an
  // EVEN target at 1/2-tick -- fire times across classes are NOT comparable, and a class flip
  // (e.g. FREE absorbed into an EVEN blossom) must void the record (the 345/348 residual).
  bool mVbEnabled;   // REGIME GATE (2026-07-13, the islands bisect): the Fig-9 filter's repair
                     // discipline is lethal on DEGENERATE-TIE forests (uniform weights, many tiny
                     // trees: islands 1M DNF>600s with it, 34.6s without) while winning 15-25% on
                     // distinct/tie-dense. Decided per search by a 5-slot distinct-weight scan
                     // (vbRegimeDecide): OFF iff the weight multiset has <=4 distinct values --
                     // exactly the degenerate regime, where it measures neutral-or-lethal.
  void vbRegimeDecide() {
    LL seen[5]; int k = 0;
    for (Sz u = 0; u < mN && k <= 4; ++u)
      for (Sz j = 0; j < mAdjW[u].size() && k <= 4; ++j) {
        LL v = mAdjW[u][j]; bool nw = true;
        for (int t = 0; t < k; ++t) if (seen[t] == v) { nw = false; break; }
        if (nw) { if (k < 5) seen[k] = v; ++k; }
      }
    mVbEnabled = (k > 4);
  }
  bool vbSuppress(Sz w, LL fire, char cls) {
    if (!mVbEnabled) return false;   // filter off: never suppress (and record nothing)
    if (mVBestFire.size() <= w) { mVBestFire.resize(mN + 8, 0); mVBestEp.resize(mN + 8, -1); mVBestCls.resize(mN + 8, 0); }
    if (mVBestEp[w] == mVBestEpoch && mVBestCls[w] == cls && mVBestFire[w] <= fire) {
#ifdef PS5_EVTRACE
      std::fprintf(stderr,"[TR] VB-SUP w=%lu fire=%lld cls=%d (best=%lld)\n",(unsigned long)w,(long long)fire,(int)cls,(long long)mVBestFire[w]);
#endif
#ifdef S5DBG
      ++gS5VbSup;
#endif
      return true;
    }
    mVBestEp[w] = mVBestEpoch; mVBestFire[w] = fire; mVBestCls[w] = cls; return false;
  }
  // v1 REPAIR DISCIPLINE (the missing half that made v0 unsound): whenever a popped entry for
  // target x is DISCARDED (stale source, absorbed, slack-not-there), x's suppressed alternatives
  // must be re-offered -- rescan x's adjacency for the earliest eligible edge INTO x and schedule
  // it (Gabow's key-recomputation analog, O(deg(x)) per discard). Without this, a vertex whose
  // best entry dies loses ALL coverage (the 17/348 failure).
  void vbRepair(Sz x, LL tNow) {
    if (!mVbEnabled) return;
#ifdef S5DBG
    ++gS5VbRepairs; if (x < mN) gS5VbRepairScans += (long long)mAdjV[x].size();
#endif
    if (x >= mN) return;
    Sz X = ((ScalingS5Engine*)this)->outerNode(x);
    int lx = mLabel[X];
    if (lx != LBL_FREE && !treeLive(X)) lx = LBL_FREE;
    if (lx != LBL_EVEN && lx != LBL_FREE) { mVBestEp[x >= mVBestEp.size() ? 0 : x] = -1; return; }
    LL tgt = (mCrit == CRIT1) ? 0 : -2; LL best = -1; Sz bu = NIL();
    Sz ufx = ((ScalingS5Engine*)this)->ufFind(x);
    for (Sz k = 0; k < mAdjV[x].size(); ++k) {
      Sz u2 = mAdjV[x][k];
      if (mIgnoreDummies && u2 >= mNReal) continue;
      if (ufx == ((ScalingS5Engine*)this)->ufFind(u2)) continue;
      Sz U = ((ScalingS5Engine*)this)->outerNode(u2);
      if (U == X) continue;
      if (mLabel[U] != LBL_EVEN || !treeLive(U)) continue;
      LL s = ((ScalingS5Engine*)this)->slackZW(u2, x, 2 * mAdjW[x][k]);
      LL d = (lx == LBL_EVEN) ? (s - tgt) / 2 : (s - tgt);
      if (mAllowOvertight && d < 0) d = 0;
      if (!mAllowOvertight && d < 0) continue;
      if (best < 0 || tNow + d < best) { best = tNow + d; bu = u2; }
    }
    if (mVBestFire.size() <= x) { mVBestFire.resize(mN + 8, 0); mVBestEp.resize(mN + 8, -1); }
    if (best >= 0) {
      if (mVBestCls.size() <= x) mVBestCls.resize(mN + 8, 0);
      mVBestEp[x] = mVBestEpoch; mVBestFire[x] = best; mVBestCls[x] = (lx == LBL_EVEN) ? 1 : 0;
      mPQ.schedule(best, (long long)bu, (long long)x, (lx == LBL_EVEN) ? 1 : 0);
    } else mVBestEp[x] = -1;   // no eligible edge now; future scans re-offer freely
  }
#ifdef PS5_VB_DEFER_REPAIR
  // *** NEGATIVE RESULT v0 (2026-07-12, keep OFF -- UNSOUND: 292/348, weights light by 1e1-1e4 at
  // scale). The deferral window lets same-tick STRUCTURE mutations (shrink/augment/markPathUsed)
  // change x's eligibility between the death and the flush -> the flush repairs against the new
  // state and misses edges the immediate repair would have offered. The COST MODEL validated:
  // 262144 anomaly 386->124s (3x) -- but wrong. A sound v1 needs repair-at-death SEMANTICS with
  // deferred WORK (freeze (x, lx, tree ctx) at death), or the ancestor-hint walk cache instead
  // (read-only, no correctness surface). ***
  // DEFERRED REPAIR (fix #2, 2026-07-12): immediate repairs CHAIN within a tick (the repaired
  // best can die again at the next same-tick pop -> another O(deg) rescan; measured 8.2M repairs
  // / 67M repair-scans at 204800 distinct). Defer: a dead best just queues x; the flush runs at
  // the loop head ONLY when the clock is about to advance (or the PQ empties), i.e. before any
  // decision that depends on x's coverage. Within the tick, x simply has no best -- new inserts
  // for x are not suppressed (ep=-1), and a second same-tick death cannot re-trigger (vbWasBest
  // is false once cleared) -> repeat deaths dedup to ONE rescan per (x, tick). Alternatives fire
  // at >= the dead best's time = within-or-after this tick, and the flush precedes the advance,
  // so nothing is overrun.
  std::vector<Sz> mVbRepList; std::vector<LL> mVbRepStamp; LL mVbRepTok;
  void vbRepairDefer(Sz x) {
    if (x >= mN) return;
    if (mVbRepStamp.size() < mN) mVbRepStamp.assign(mN, -1);
    if (mVbRepStamp[x] == mVbRepTok) return;
    mVbRepStamp[x] = mVbRepTok; mVbRepList.push_back(x);
  }
  void vbRepairFlush(LL tNow) {
    for (Sz i = 0; i < mVbRepList.size(); ++i) vbRepair(mVbRepList[i], tNow);
    mVbRepList.clear(); ++mVbRepTok;
  }
#endif
#endif
                                          // (hBase: horizon window base for PS5_EVT_HORIZON marker
                                          // rescans; -1 = tNow. Ignored without the flag.)
  void applyDualTo(LL newTime, LL* tNow); // advance materialized duals to a time
  bool searchEvent();               // one event-driven EdmondsSearch (no rescan)
#ifdef PS5_SFM
  // CONSISTENCY GATE (no behavior change): compute the next-event fire-time (min over all currently
  // schedulable EVEN-incident edges) two ways -- a brute O(m) scan (ground truth = what the reseed's PQ
  // would pop) vs the split-findmin findmin -- and assert they agree. Validates the query semantics
  // against the engine's LIVE state before split-findmin replaces the reseed. Returns the brute min
  // fire-time (or -1 if none); records (mSfmChecks, mSfmMismatch).
  LL sfmBruteNextFire(LL tNow) const;     // ground-truth min fire-time by full scan
  LL sfmFindminNextFire(LL tNow);         // same via split-findmin (rebuilt over EVEN-incident edges)
  void sfmInit();                         // init structure at search start (all keys INF), enable feed
  LL sfmFindminLive();                    // next-event from the LIVE incrementally-maintained structure
  void sfmCheckAfterReseed(LL tNow);      // run both, compare, count
  Sz mSfmChecks; Sz mSfmMismatch;
#endif
  bool mUseEvent;                   // solve() dispatch: event engine vs scan engine
  bool mEventStuck;                 // searchEvent hit an M==E it can't extend -> solve finishes via scan

  // ---- §5.2 LAZY TIMESTAMP DUALS ----
  // Within one search call, dual adjustments are recorded as t_now advancing (no
  // per-step y/z writes). y(v)/z(B) are computed ON DEMAND from the timestamps, so
  // a dual adjustment is O(1) and -- crucially -- z(B) is EXACT regardless of how
  // B's inner/outer role changed, which fixes the event-engine dissolution race.
  //   t_now      : current time within the call (sum of deltas applied)
  //   mY0/mZ0    : y/z at call start (the materialized base)
  //   mTIn[B]    : time B last became an INNER root blossom (NIL if not inner-root)
  //   mTOut[B]   : time B last became part of an OUTER root blossom (NIL otherwise)
  //   mDelta[B]  : accumulated inner-increments before the current interval
  // A node's role (for the formula) is read from mLabel of its OUTER node:
  //   B is INNER-root  iff outerNode(B)==B and mLabel[B]==LBL_ODD
  //   B is in an OUTER iff mLabel[outerNode(B)]==LBL_EVEN
  LL mTNow;
  bool mLazyActive;                 // slackZ reads lazy duals when true
  std::vector<LL> mY0, mZ0;         // materialized base at search start
  // PER-VERTEX y accumulation (ported from ScalingFastMWPM, valgrind-clean): mYAcc[v]
  // = y-change from v's CLOSED role-intervals; open interval from mVRole[v]/mVStamp[v].
  // Threads y across role/containment changes (a vertex absorbed into a blossom that
  // flips role keeps an exact y). Funnel all role changes through setVtxRole.
  std::vector<LL> mYAcc; std::vector<int> mVRole; std::vector<LL> mVStamp;
  // PER-NODE z accumulation: mZAcc[B] closed-interval z-change; open from mNRole/mNStamp.
  std::vector<LL> mZAcc; std::vector<int> mNRole; std::vector<LL> mNStamp;
  // lazySyncRoles DIRTY-SET (2026-07-02): the full O(n) vertex re-sync is ~91% wasted no-ops (measured).
  // markLabelDirty(X) records a node whose label changed -> lazySyncRoles syncs only those nodes' leaves.
  // mSyncDirtyAll forces the safe full scan (set by any site that can't cheaply enumerate what changed).
  // *** NEGATIVE RESULT (default-OFF, keep PS5_SYNC_DIRTY off): CORRECT (LEMON-optimal to n=25600, 348/348)
  // but SLOWER on blossom-heavy inputs (n=12800 W=1000: 1.75s -> 3.0s). Root flaw: marking a NODE and
  // iterating ALL its leaves is O(blossom size) per mark; a big blossom (maxBloss 6213) re-marked every
  // grow costs K*size, which EXCEEDS the O(n) full scan. The truly-changed vertices are only ~9% but you
  // can't cheaply know WHICH without the scan. Incr% was 91-98% (dirty path IS taken) yet still loses.
  // A correct fix would track individual role-FLIPPING vertices, not nodes -- deferred/likely not worth it.
  // Default keeps mSyncDirtyAll=true; PS5_SYNC_DIRTY enables the (slower) incremental path for the record. ***
  std::vector<Sz> mSyncDirty; bool mSyncDirtyAll;
  void markLabelDirty(Sz X) {
#ifdef PS5_SYNC_DIRTY
    if (!mSyncDirtyAll) mSyncDirty.push_back(X);
#endif
#ifdef PS5_FLIP_DIRTY
    if (!mFlipAll) mFlipNodes.push_back(X);
#endif
#if !defined(PS5_SYNC_DIRTY) && !defined(PS5_FLIP_DIRTY)
    (void)X;
#endif
  }
#ifdef PS5_FLIP_DIRTY
  // EVENT-DRIVEN per-node role-flip list (requires PS5_NODE_Y). Every mLabel write that could change a
  // node's role pushes it here via setLabel(); lazySyncRoles processes ONLY these (O(#flips) via
  // setNodeRole, which under node-y settles members lazily), replacing the O(n) per-advance scan -- the
  // persistent-clock enabler. mFlipAll forces a full sync when a site can't cheaply enumerate (safety).
  std::vector<Sz> mFlipNodes; bool mFlipAll;
  void setLabel(Sz X, int r) { if (mLabel[X] != r) { mLabel[X] = r; if (!mFlipAll) mFlipNodes.push_back(X); } }
  void flipMarkAll() { mFlipAll = true; }
#else
  void setLabel(Sz X, int r) { mLabel[X] = r; }
#endif
#ifdef PS5_INVKEY
#if !(defined(PS5_PERSIST_CLOCK) && defined(PS5_PERSIST_FOREST))
#error PS5_INVKEY requires PS5_PERSIST_CLOCK + PS5_PERSIST_FOREST (absolute fire times + mReoffer)
#endif
  // *** NEGATIVE RESULT (2026-07-11, keep OFF -- UNSOUND: 228/348, 522-557/594, speed probes
  // stall). The Stage-0 "100% repeat" at relabel boundaries measured ROOT persistence only:
  // relabel() TEARS DOWN the forest (every top-level label -> FREE, only free roots re-seeded
  // EVEN), so the post-boundary EVEN set is just the surviving roots. The torn trees change
  // every queued event's TARGET rate class, and roots' own events were consumed mid-phase --
  // a source-side label diff can neither validate the survivors nor restore the consumed.
  // VERDICT: at a teardown boundary the existing reseed is already near-minimal (it scans the
  // new roots, not all of n). The faithful resolution of the re-based-frame deviation is to
  // ELIMINATE teardowns (PS5_PERSIST_SHRINK viability + relink decline coverage), not to diff
  // them. See docs/N_Scope_InvariantKeys.md epilogue. ***
  // Original (falsified) design note: keep the PQ, diff the labels, reschedule only changes.
  std::vector<long long> mIkEp;   // mIkEp[i]==mIkTok  <=>  i was top-level EVEN pre-boundary
  long long mIkTok;
  void ikSnapshot() {             // call BEFORE the boundary's structure mutations
    ++mIkTok;
    if (mIkEp.size() < mZ.size()) mIkEp.resize(mZ.size(), 0);
    for (Sz i = 0; i < mNumNodes; ++i)
      if (mBParent[i] == NIL() && mLabel[i] == LBL_EVEN) mIkEp[i] = mIkTok;
  }
  void ikReschedChanged(LL tNow) {  // call AFTER relabel+lazyBegin: newly-EVEN nodes only
    if (mIkEp.size() < mZ.size()) mIkEp.resize(mZ.size(), 0);
    mPQ.syncNow(tNow);              // the cursor may have peeked ahead of the true clock
    for (Sz i = 0; i < mNumNodes; ++i)
      if (mBParent[i] == NIL() && mLabel[i] == LBL_EVEN && mIkEp[i] != mIkTok)
        scheduleNodeEdges(i, tNow);
    // Events discarded this phase ONLY for mUsed phase-disjointness are still eligible; the
    // full reseed used to restore them implicitly -- the kept-PQ boundary must drain mReoffer
    // (same contract as the persistent pure-augment boundary).
    for (Sz ri = 0; ri < mReoffer.size(); ++ri) {
      const BucketPQ::Event& rev = mReoffer[ri];
      mPQ.schedule(rev.time, rev.a, rev.b, rev.kind);
    }
    mReoffer.clear();
  }
#endif
  std::vector<LL> mLcaSeen;         // lcaOuter scratch (instance, not static). LL: the token
  LL mLcaTok;                       // increments per lcaOuter call (1.3B+/run measured) -- a
                                    // 32-bit token would WRAP -> stale-token collisions (IDX32).
#ifdef PS5_ANCESTOR_HINT
  // O_Scope ancestor hints (docs/O_Scope_AncestorHint.md): per-vertex last-known top-level
  // ancestor + validity. mAhDead[B] set at dissolve (ids not reused within a search: the
  // free-list push is disabled); mAhTok bumps at forest resets (where ids DO get reused).
  mutable std::vector<Sz> mAhHint; mutable std::vector<LL> mAhStamp;
  std::vector<char> mAhDead; LL mAhTok;
  void ahEnsure() {
    if (mAhStamp.size() < mN) { mAhHint.resize(mN, NIL()); mAhStamp.resize(mN, -1); }
    if (mAhDead.size() < mZ.size()) mAhDead.resize(mZ.size(), 0);
  }
#endif
  std::vector<Sz> mFreeBlossoms;    // reclaimed blossom ids (free-list)
  LL yLazy(Sz v) const;             // current y(v) from timestamps
  LL zLazy(Sz B) const;             // current z(B) from timestamps
  void lazyBegin();                 // snapshot mY/mZ -> base; seed roles; t_now=0
  void lazyMaterialize();           // write yLazy/zLazy back into mY/mZ; end of call
  void setVtxRole(Sz v, int role);  // settle vertex open interval -> mYAcc, re-stamp
  void setNodeRole(Sz B, int role); // settle node z open interval -> mZAcc, re-stamp
  void lazySyncRoles();             // sync lazy roles to current materialized labels
  void lazyVerifyStep();            // assert yLazy/zLazy == materialized (verify mode)
  bool mLazyVerify;                 // debug: assert yLazy/zLazy == materialized in scan
  bool mLazyEvent;                  // S5-C: searchEvent drives on lazy duals (mLazyActive)
  bool mNoRelabelShrink;            // S5-C step 2: incremental no-relabel shrink
  bool mNoRelabelAug;               // S5-C step 3: per-level maximal multi-augment (MV/HK phase)
  bool mNoRelabelDissolve;          // S5-C step 4: in-place blossom dissolution (no relabel)
  int  mCrit;                       // eligibility criterion (CRIT1 default; CRIT2/CRIT3 for FVR). Event target = (CRIT1?0:-2).
  bool mIgnoreDummies;               // D: scheduleNodeEdges skips edge targets >= mNReal (dummy slots)
  bool mAllowOvertight;              // D: carried-dual (Property-2) regime -- slack<tgt is eligible now
  Sz   mNReal;                       // D bridge: # vertices that must be ROOTED+MATCHED (reals). mN stays
                                     // the vtx/blossom boundary (= loop blossomBase 2*mN0); the dummy
                                     // slots [mNReal,mN) are vertices but NOT rooted (no FVR edges yet).
                                     // Default mNReal==mN (standalone solve: all vertices are real).
  bool relinkOddCycleS5(Sz B);      // re-thread the tree through B's cycle (no-relabel dissolution)
  bool relinkEligibleS5(Sz B);      // non-mutating: can B be re-threaded in place? (pre-check)
  bool relinkStructEligible(Sz B, Sz depth);  // structural-only: odd cycle + members bare/z=0-recursive (for nested)
  bool relinkRecDissolve(Sz B, std::vector<Sz>* relinkedEven);  // recursively relink+dissolve a now-root z=0 blossom + its z=0 sub-blossoms
  void freeConsumedTree(Sz root);   // settle a consumed tree's roles to FREE (no-relabel augment)
#if defined(PS5_PERSIST_FOREST) && defined(PS5_PERSIST_CLOCK)
  void retireDeadTrees();           // STAGE 2: batch-settle ALL dead trees at a persistent boundary
  std::vector<BucketPQ::Event> mReoffer;  // STAGE 2: events discarded ONLY by mUsed phase-disjointness
                                    // (still eligible!) -- re-scheduled at the next persistent boundary.
                                    // Without this the edge is permanently lost (n=16 rep=11 trace: the
                                    // (3,2) formation popped at its exact fire time, was mUsed-blocked,
                                    // went overtight, and strict mode never schedules d<0).
#endif
  std::vector<char> mUsed;          // per-PHASE: vertex consumed by an augmenting path this level
  void markPathUsed(Sz startEven);  // mark all leaves on the tree path startEven..root used
  bool pathClean(Sz startEven) const;  // true iff no leaf on the path startEven..root is mUsed
public:
#ifdef PS5_STAGE0
  // Persistent-clock STAGE 0 (docs/H_Scope_PersistentClock.md): quantify redundant forest
  // regrowth. At each phase-end teardown, snapshot the forest (mParent/mLabel); during the next
  // phase, a grow that re-creates the IDENTICAL attachment (same parent, same labels) is
  // REDUNDANT -- work the persistent forest would skip. High redundancy validates Stage 2.
  std::vector<Sz> mS0Parent; std::vector<int> mS0Lbl;
  long long mS0Grows, mS0Redundant, mS0Phases;
  void s0Reset() { mS0Grows=mS0Redundant=mS0Phases=0; mS0Parent.clear(); mS0Lbl.clear(); }
  void s0Snapshot() { ++mS0Phases; mS0Parent = mParent;
    mS0Lbl.assign(mLabel.size(), 0); for (Sz i=0;i<mLabel.size();++i) mS0Lbl[i]=mLabel[i]; }
  void s0Print() const { std::fprintf(stderr,
    "[STAGE0] phases=%lld grows=%lld redundant=%lld (%.1f%%)\n",
    mS0Phases,mS0Grows,mS0Redundant,mS0Grows?100.0*mS0Redundant/mS0Grows:0.0); }
#endif
#ifdef PS5_QSTAT
  // Event-queue waste instrumentation (design probe for the reseed cost): distribution of insert
  // fire-distances (fireT - tNow at insert) vs per-phase tick advance. If inserts sit far past the
  // typical phase-end advance, a defer/horizon scheme skips them; if they cluster near tNow, the
  // waste is phase-reset churn instead. Totals across ALL searches of a solve (reset in init()).
  long long mQIns, mQPops, mQAct, mQPhases;
  long long mQInsHist[18];   // log2 buckets of insert distance d
  long long mQAdvHist[18];   // log2 buckets of phase-end advance (augTime)
  void qstatReset() { mQIns=mQPops=mQAct=mQPhases=0; for(int i=0;i<18;++i){mQInsHist[i]=0;mQAdvHist[i]=0;} }
  static int qstatBucket(LL d) { int b=0; while (d>0 && b<17) { d>>=1; ++b; } return b; }
  void qstatPrint() const {
    std::fprintf(stderr,"[QSTAT] inserts=%lld pops=%lld acted=%lld phases=%lld\n",mQIns,mQPops,mQAct,mQPhases);
    std::fprintf(stderr,"[QSTAT] insert-dist:");
    for (int i=0;i<18;++i) if (mQInsHist[i]) std::fprintf(stderr," [%s%lld]=%lld",i?"<=":"=",(long long)(i?((LL)1<<(i-1))*2-1:0),mQInsHist[i]);
    std::fprintf(stderr,"\n[QSTAT] phase-adv :");
    for (int i=0;i<18;++i) if (mQAdvHist[i]) std::fprintf(stderr," [%s%lld]=%lld",i?"<=":"=",(long long)(i?((LL)1<<(i-1))*2-1:0),mQAdvHist[i]);
    std::fprintf(stderr,"\n");
  }
#endif
  // final materialized duals (for post-solve slack / complementary-slackness diagnostics)
  LL finalY(Sz v) const { return (v < mY.size()) ? mY[v] : 0; }
  LL finalSlack(Sz u, Sz v) const { return ((ScalingS5Engine*)this)->slackZW(u, v, 2 * edgeW(u, v)); }
  LL edgeW(Sz u, Sz v) const { for (Sz k = 0; k < mAdjV[u].size(); ++k) if (mAdjV[u][k] == v) return mAdjW[u][k]; return -1; }
  void setUseEvent(bool e) { mUseEvent = e; }
  void setLazyVerify(bool e) { mLazyVerify = e; }
  void setLazyEvent(bool e) { mLazyEvent = e; }   // S5-C: event engine drives on lazy duals
  void setNoRelabelShrink(bool e) { mNoRelabelShrink = e; }  // S5-C step 2: incremental shrink
  void setNoRelabelAug(bool e) { mNoRelabelAug = e; }        // S5-C step 3: persist across augments
  void setNoRelabelDissolve(bool e) { mNoRelabelDissolve = e; }  // S5-C step 4: in-place dissolution
  void setCrit(int c) { mCrit = c; }   // D: eligibility criterion for the per-scale search (CRIT1/2/3)
#ifdef S5DBG
  // HONEST STRUCTURAL PROBE (2026-07-07, user challenged "is it really giant blossoms?"): dump the
  // REAL live top-level blossom sizes -- counted independently, and check the biggest for DUPLICATE
  // vertices in its leaf list (which would inflate gS5MaxBlossom without a true giant blossom).
  void dbgBlossomCensus() const {
    Sz nTop=0, nBloss=0, biggest=0, bigB=NIL();
    for (Sz b = 0; b < mNumNodes; ++b) {
      if (mBParent[b] != NIL()) continue;              // not top-level
      ++nTop;
      if (b >= mN) { ++nBloss; Sz sz = mNodeLeaves[b].size();
        if (sz > biggest) { biggest = sz; bigB = b; } }
    }
    std::fprintf(stderr,"[CENSUS] topLevelNodes=%lu liveBlossoms=%lu biggestBlossomLeaves=%lu\n",
      (unsigned long)nTop,(unsigned long)nBloss,(unsigned long)biggest);
    if (bigB != NIL()) {
      // count DISTINCT vertices in the biggest blossom's leaf list (dup detection)
      const std::vector<Sz>& lv = mNodeLeaves[bigB];
      std::vector<char> seen(mN, 0); Sz distinct=0, dups=0;
      for (Sz i=0;i<lv.size();++i){ Sz v=lv[i]; if(v<mN){ if(seen[v]) ++dups; else {seen[v]=1; ++distinct;} } }
      std::fprintf(stderr,"[CENSUS] biggest blossom %lu: leafList=%lu distinctVerts=%lu dups=%lu (%.1f%% of n)\n",
        (unsigned long)bigB,(unsigned long)lv.size(),(unsigned long)distinct,(unsigned long)dups,
        100.0*(double)distinct/(double)mN);
    }
  }
  // Diagnostic: dump the hot-loop counters (edge-scan volume vs m tells super-linear vs linear).
  void dbgDumpCounters(unsigned long m) const {
    std::fprintf(stderr,
      "[S5DBG] edgeScans=%lld internalScans=%lld (%.1f%%) schedCalls=%lld nodeReSched=%lld inserts=%lld disRelabel=%lld phaseRelabel=%lld persistPhase=%lld maxBlossom=%lld  edgeScans/m=%.1f\n",
      gS5EdgeScans, gS5InternalScans,
      gS5EdgeScans ? 100.0*(double)gS5InternalScans/(double)gS5EdgeScans : 0.0,
      gS5SchedCalls, gS5NodeReSched, gS5SchedInserts,
      gS5DissolveRelabels, gS5PhaseRelabels, gS5PersistPhases, gS5MaxBlossom,
      m ? (double)gS5EdgeScans/(double)m : 0.0);
    std::fprintf(stderr,"[S5DBG] maxBlossomLeaves=%lld maxBlossomDISTINCT=%lld (of n) -- if distinct<<leaves, the 'giant blossom' is leaf-list DUPLICATION not a real blossom\n",
      gS5MaxBlossom, gS5MaxDistinct);
    std::fprintf(stderr,"[S5DBG] PREDICTION: replayedEvents=%lld  w-LABEL-CHANGED=%lld (%.2f%%) -- if SMALL, incremental touch-changed-only fix is O(changes) not O(size)\n",
      gS5ReplayTotal, gS5ReplayWchanged, gS5ReplayTotal ? 100.0*(double)gS5ReplayWchanged/(double)gS5ReplayTotal : 0.0);
    std::fprintf(stderr,"[S5DBG] STAGE-A: IncrSkip=%lld (big-node scans that skipped re-emission, events still live)\n", gS5IncrSkip);
    std::fprintf(stderr,"[S5DBG] lazySyncRoles: calls=%lld (%.2f per vtx) changed=%lld  (%.3f%% of calls*n did real role-flips = the reducible fraction)\n",
      gS5SyncCalls, mN ? (double)gS5SyncCalls/(double)mN : 0.0, gS5SyncChanged,
      (gS5SyncCalls>0 && mN>0) ? 100.0*(double)gS5SyncChanged/((double)gS5SyncCalls*(double)mN) : 0.0);
    std::fprintf(stderr,
      "[S5DBG] bigNodeCalls(>%lld leaves)=%lld  bigNodeLeafSum=%lld  avgBigLeaves=%.0f  viaMarker=%lld viaReseed=%lld\n",
      gS5BigNodeThresh, gS5BigNodeCalls, gS5BigNodeLeafSum,
      gS5BigNodeCalls ? (double)gS5BigNodeLeafSum/(double)gS5BigNodeCalls : 0.0,
      gS5BigViaMarker, gS5BigViaReseed);
    std::fprintf(stderr,
      "[S5DBG] bigRedundant=%lld bigNecessary=%lld  MEMOIZABLE=%.1f%%\n",
      gS5BigRedundant, gS5BigNecessary,
      (gS5BigRedundant+gS5BigNecessary) ? 100.0*(double)gS5BigRedundant/(double)(gS5BigRedundant+gS5BigNecessary) : 0.0);
  }
#endif
  // STEP-6 LEVELING SUPPORT (Fig 5 step 6, docs/Scaling_Paper_Exact_Spec.md L75-81):
  // Root-restricted search: when set, relabel() roots ONLY the listed free vertices (the spec's
  // PQSearch(F) with F = a specific root set). Empty set / clear = root all free reals (default).
  void setRootSet(const std::vector<Sz>& roots) {
    mRootRestrict = true;
    if (mRootAllow.size() < mN) mRootAllow.assign(mN, 0); else std::fill(mRootAllow.begin(), mRootAllow.end(), 0);
    for (Sz i = 0; i < roots.size(); ++i) if (roots[i] < mN) mRootAllow[roots[i]] = 1;
  }
  void clearRootSet() { mRootRestrict = false; }
  // Enumerate the MAXIMAL SMALL blossoms (|leaves| < tauSz): a small top-level blossom is emitted
  // whole; a LARGE one is recursed so its outermost small sub-blossoms are emitted (the spec's
  // "maximal old small blossom" after step-3 liquidated the large shells). Call BEFORE the
  // per-scale reset wipes the blossom structure.
  void collectMaximalSmall(Sz tauSz, std::vector<std::vector<Sz> >* out) const {
    out->clear();
    for (Sz b = mN; b < mNumNodes; ++b)
      if (mBParent[b] == NIL()) collectMaximalSmallRec(b, tauSz, out);
  }
private:
  void collectMaximalSmallRec(Sz b, Sz tauSz, std::vector<std::vector<Sz> >* out) const {
    if (b < mN) return;                            // bare vertex: no blossom to level
    if (mNodeLeaves[b].size() < tauSz) { out->push_back(mNodeLeaves[b]); return; }
    for (Sz k = 0; k < mBCycle[b].size(); ++k) collectMaximalSmallRec(mBCycle[b][k], tauSz, out);
  }
  std::vector<char> mRootAllow;      // step-6: allowed roots when mRootRestrict
  bool mRootRestrict = false;
public:
  void setReals(Sz r) { mNReal = r; }  // D bridge: only the first r vertices must be rooted+matched
  void setIgnoreDummies(bool e) { mIgnoreDummies = e; }  // skip edges to vertices >= mNReal (dummy slots)
  void setAllowOvertight(bool e) { mAllowOvertight = e; }  // D: Property-2 regime (carried duals)
};

// ---------------------------------------------------------------------------
// slack helpers
// ---------------------------------------------------------------------------
inline ScalingS5Engine::LL
ScalingS5Engine::w2max(Sz u, Sz v) const {
  LL best = 0; bool found = false;
#ifdef S5DBG
  ++gS5SlackZCalls; gS5SlackZSteps += (long long)mAdjV[u].size();
#endif
  for (Sz i = 0; i < mAdjV[u].size(); ++i)
    if (mAdjV[u][i] == v && (!found || mAdjW[u][i] > best)) { best = mAdjW[u][i]; found = true; }
  if (!found) return -(LL)4 * (mMaxW + 1) * (LL)(mN + 1);  // not an edge: huge -slack
  return 2 * best;
}

// slackZ(u,v) = y(u)+y(v) + sum z(B enclosing both) - 2w(u,v). Sum z over the
// common suffix of u's and v's blossom-parent chains (O(depth)). Reads duals via
// the lazy timestamps when mLazyActive (event engine), else materialized.
inline ScalingS5Engine::LL
ScalingS5Engine::slackZ(Sz u, Sz v) const { return slackZW(u, v, w2max(u, v)); }

// slackZ with 2*w(u,v) already in hand (the hot caller scheduleNodeEdges has the weight
// as mAdjW[u][k]) -- avoids the O(degree) w2max adjacency rescan (PROFILED as the #1 cost,
// ~1.7B steps = 8x the edge-scans at n=6400). Mirrors ScalingFastMWPM's slackZW win.
inline ScalingS5Engine::LL
ScalingS5Engine::slackZW(Sz u, Sz v, LL w2) const {
  LL yu = mLazyActive ? yLazy(u) : mY[u];
  LL yv = mLazyActive ? yLazy(v) : mY[v];
  LL s = yu + yv - w2;
  Sz a = mBParent[u], b = mBParent[v];
  if (a == NIL() || b == NIL()) return s;
#ifdef PS5_OUTERCACHE
  // FAST PATH: the z-sum term is z over blossoms containing BOTH u and v = the common ancestors from
  // their LCA up to the shared root. If u,v are in DIFFERENT outer blossoms (the common case: cross-tree
  // grow/formation edges), they share NO blossom -> z-sum is 0 -> skip the O(depth) LCA walk entirely.
  // outerNode is O(1) amortized via the cache, so this short-circuit is cheap. (Was O(depth) x millions.)
  if (outerNode(u) != outerNode(v)) return s;
#endif
  Sz da = 0, db = 0;
  for (Sz x = a; x != NIL(); x = mBParent[x]) ++da;
  for (Sz x = b; x != NIL(); x = mBParent[x]) ++db;
  while (da > db) { a = mBParent[a]; --da; }
  while (db > da) { b = mBParent[b]; --db; }
  while (a != b) { a = mBParent[a]; b = mBParent[b]; }
  for (Sz x = a; x != NIL(); x = mBParent[x]) s += (mLazyActive ? zLazy(x) : mZ[x]);
  return s;
}

// ---------------------------------------------------------------------------
// §5.2 lazy-timestamp dual readers + maintenance (PER-VERTEX/PER-NODE accumulation,
// ported from ScalingFastMWPM.h's valgrind-clean scheme). y(v) = y0 + closed-interval
// accumulation (mYAcc) + open-interval term from v's CURRENT role/stamp. Role changes
// funnel through setVtxRole/setNodeRole (settle the open interval into Acc, re-stamp),
// so duals stay EXACT across role/containment changes -- including dissolution -- with
// NO rebuild. A dual adjustment is just mTNow++ (O(1)).
// ---------------------------------------------------------------------------
inline ScalingS5Engine::LL
ScalingS5Engine::yLazy(Sz v) const {
#ifdef PS5_NODE_Y
  // FAITHFUL §5.2 (2026-07-08, docs/J_Scope): a member's in/out accrual = its TOP-LEVEL entity's
  // role (single-level -- mYAcc[v] already holds everything frozen from PRIOR top-levels, settled
  // at each membership change in shrink/dissolve). NO per-advance O(n) role scan. Bare vertex uses
  // its own mVRole; blossom member uses the blossom's mNRole/mNStamp. Rate ODD +1 / EVEN -1.
  Sz On = ((ScalingS5Engine*)this)->outerNode(v);
  LL y = mY0[v] + mYAcc[v];
  int r; LL stamp;
  if (On == v) { r = mVRole[v]; stamp = mVStamp[v]; }
  else         { r = mNRole[On]; stamp = mNStamp[On]; }
  if (r == LBL_ODD)       y += (mTNow - stamp);
  else if (r == LBL_EVEN) y -= (mTNow - stamp);
  return y;
#else
  LL y = mY0[v] + mYAcc[v];
  if (mVRole[v] == LBL_ODD)  y += (mTNow - mVStamp[v]);
  else if (mVRole[v] == LBL_EVEN) y -= (mTNow - mVStamp[v]);
  return y;
#endif
}
inline ScalingS5Engine::LL
ScalingS5Engine::zLazy(Sz B) const {
  if (B < mN) return 0;
  LL z = mZ0[B] + mZAcc[B];
  // (Tested PS5_ZLAZY_ROOTGUARD = freeze non-root z: INERT -- mNRole is already FREE on non-roots, so that
  // was not the drift. The drift's exact location needs the persistent shadow-z detector, not more guesses.)
  if (mNRole[B] == LBL_EVEN) z += 2 * (mTNow - mNStamp[B]);
  else if (mNRole[B] == LBL_ODD) z -= 2 * (mTNow - mNStamp[B]);
  return z;
}
// snapshot materialized duals as base; seed each vertex/node's open-interval role from
// its current label; t_now=0.
inline void ScalingS5Engine::lazyBegin() {
#ifdef S5DBG
  { extern long long gS5LazyBeginCalls; ++gS5LazyBeginCalls; }
#endif
#ifdef PS5_PERSIST_CLOCK
  // STAGE 1 (docs/H_Scope_PersistentClock.md): MONOTONE clock -- the re-snapshot no longer rewinds
  // mTNow. The snapshot is value-preserving at the boundary instant (yLazy is continuous across
  // materialize+begin: mY0 = just-materialized mY, stamps = mTNow), so dual VALUES are unchanged;
  // what changes is that event fire-times (absolute clock values) become survivable across phase
  // boundaries (Stage 2). searchEvent entry resets mTNow=0 (fresh search per scale = faithful);
  // the teardown sites' `tNow = mTNow` pick up the persistent value.
#else
  mTNow = 0;
#endif
  // NOTE (2026-07-02): tried skipping the mY0/mVRole/mVStamp .assign zero-fills here (they ARE fully
  // overwritten by the loops below) -- MEASURED INERT: lazyBegin's cost is the LOOPS (mY0[v]=mY[v],
  // mVRole[v]=mLabel[outerNode(v)] with the per-vertex outerNode walk, mZ0[b]=mZ[b]), not the fills
  // (callgrind: 781M instr both ways, total -1.1%, wall unchanged). The earlier "lazyBegin ~13%" reading
  // aggregated memset/_M_fill_assign across ALL callers, not this phase-reset. Reverted. The real per-phase
  // cost to attack is lazySyncRoles (24%) doing an O(n) role re-sync per time-advance -- needs a dirty-set.
  mY0.assign(mN, 0); mYAcc.assign(mN, 0); mVRole.assign(mN, LBL_FREE); mVStamp.assign(mN, mTNow);
  for (Sz v = 0; v < mN; ++v) mY0[v] = mY[v];
#ifdef SHADOW_DUAL
  mShadowY.assign(mN, 0); for (Sz v = 0; v < mN; ++v) mShadowY[v] = mY[v];  // seed = materialized base
#endif
  // size to mZ.size() (the current node CAPACITY, grown by ensureCap), NOT mNumNodes:
  // a no-relabel shrink creates new blossom nodes mid-search whose lazy z-fields must
  // exist (z0=0/Acc=0/FREE = correct state at creation). (Was the ASan heap-overflow.)
  Sz cap = mZ.size();
  mZ0.assign(cap, 0); mZAcc.assign(cap, 0); mNRole.assign(cap, LBL_FREE); mNStamp.assign(cap, mTNow);
  for (Sz b = 0; b < mNumNodes; ++b) mZ0[b] = mZ[b];
  for (Sz v = 0; v < mN; ++v) { Sz B = outerNode(v); mVRole[v] = mLabel[B]; mVStamp[v] = mTNow; }
  for (Sz b = mN; b < mNumNodes; ++b)
    if (mBParent[b] == NIL()) { mNRole[b] = mLabel[b]; mNStamp[b] = mTNow; }
#ifdef PS5_SYNC_DIRTY
  // lazyBegin just did a FULL sync (mVRole set for every v) -> the incremental tracking starts clean.
  mSyncDirty.clear(); mSyncDirtyAll = false;
#endif
#ifdef PS5_FLIP_DIRTY
  mFlipNodes.clear(); mFlipAll = false;   // full sync done; from here track flips incrementally
#endif
}
// fold lazy values back into materialized mY/mZ (call end). Read BEFORE clearing.
inline void ScalingS5Engine::lazyMaterialize() {
  std::vector<LL> ny(mN);
  for (Sz v = 0; v < mN; ++v) ny[v] = yLazy(v);
  std::vector<LL> nz(mNumNodes, 0);
  for (Sz b = mN; b < mNumNodes; ++b) nz[b] = zLazy(b);
  for (Sz v = 0; v < mN; ++v) mY[v] = ny[v];
  for (Sz b = mN; b < mNumNodes; ++b) mZ[b] = nz[b];
#ifdef PS5_ZSHADOW
  // verify the MATERIALIZED duals (what finalSlack/dissolution read) match the persistent shadow.
  if (mZShadowOk) {
    for (Sz v = 0; v < mN; ++v) if (mY[v] != mYShadow[v] && gZShadowB < 0) {
      std::fprintf(stderr,"[MATERIALIZE-Y-DIVERGE v=%lu mat=%lld shadow=%lld @mTNow=%lld]\n",(unsigned long)v,(long long)mY[v],(long long)mYShadow[v],(long long)mTNow); gZShadowB=mNumNodes+v; }
    for (Sz b = mN; b < mNumNodes; ++b) if (mBParent[b]==NIL() && mZ[b] != mZShadow[b] && gZShadowB < 0) {
      std::fprintf(stderr,"[MATERIALIZE-Z-DIVERGE B=%lu mat=%lld shadow=%lld @mTNow=%lld]\n",(unsigned long)b,(long long)mZ[b],(long long)mZShadow[b],(long long)mTNow); gZShadowB=b; }
  }
#endif
}
// settle a VERTEX's open interval into mYAcc, then set new role + stamp.
inline void ScalingS5Engine::setVtxRole(Sz v, int role) {
#ifdef PS5_NODE_Y
  // FAITHFUL §5.2: a vertex that is currently a BLOSSOM MEMBER accrues y through its outer blossom's
  // role (setNodeRole settles it), NOT through its own mVRole. A stray setVtxRole on a member (e.g. an
  // augment/relink path that re-touches the vertex while its blossom still lists it) would DOUBLE-count
  // the same interval (the n=16 rep=1 bug). Keep mVRole/mVStamp current for when it later becomes bare,
  // but do NOT accrue -- the blossom owns the accrual. Bare vertices (mBParent==NIL) accrue as before.
  if (mBParent[v] != NIL()) {
#ifdef VTXWATCH
    if((long long)v==(VTXWATCH)) std::fprintf(stderr,"[VTXROLE-MEMBER-SKIP v=%lu %d->%d @tNow=%lld inBlossom=%ld]\n",(unsigned long)v,mVRole[v],role,(long long)mTNow,(long)mBParent[v]);
#endif
    mVRole[v] = role; mVStamp[v] = mTNow; return;
  }
#endif
#ifdef VTXWATCH
  if((long long)v==(VTXWATCH)) std::fprintf(stderr,"[VTXROLE v=%lu %d->%d @tNow=%lld stamp=%lld yAcc=%lld interval=%lld]\n",
      (unsigned long)v,mVRole[v],role,(long long)mTNow,(long long)mVStamp[v],(long long)mYAcc[v],(long long)(mTNow-mVStamp[v]));
#endif
  if (mVRole[v] == LBL_ODD)  mYAcc[v] += (mTNow - mVStamp[v]);
  else if (mVRole[v] == LBL_EVEN) mYAcc[v] -= (mTNow - mVStamp[v]);
  mVRole[v] = role; mVStamp[v] = mTNow;
}
// settle a NODE's z open interval into mZAcc, then set new role + stamp.
inline void ScalingS5Engine::setNodeRole(Sz B, int role) {
  if (B < mN) return;
#ifdef PS5_NODE_Y
  // FAITHFUL §5.2: a member's y-rate is driven by its OUTER BLOSSOM's role (yLazy reads mNRole[On]).
  // So whenever B's role changes we must ALSO settle every member's open y-interval at B's OLD role
  // into mYAcc, before restamping -- exactly what setVtxRole does for a bare vertex, funnelled through
  // the blossom. This makes EVERY setNodeRole caller (shrink's cycle-clear, lazySyncRoles, dissolve)
  // correct uniformly. Only touches members when the ROLE actually changes; a stable blossom is O(1).
  if (mNRole[B] != role && mNRole[B] != LBL_FREE) {
    LL span = mTNow - mNStamp[B];
    if (span != 0) {
      // mNodeLeaves[B] can list vertices that have SINCE left B (append-only at shrink; a leaf can be
      // freed by a nested dissolve without B's list being pruned). Only settle a leaf that CURRENTLY
      // resolves to B, else we double-count against its new top-level's settle (the n=16 rep=1 bug).
      // Guard runs only on a ROLE CHANGE (rare), so no per-advance cost is reintroduced.
      int d = (mNRole[B] == LBL_ODD) ? +1 : (mNRole[B] == LBL_EVEN ? -1 : 0);
#ifdef PS5_LEAF_SPLICE
      if (d) MBX_LF_FOR(v, B) { if (v >= mN) continue;
#else
      std::vector<Sz>& lv = mNodeLeaves[B];
      if (d) for (Sz j = 0; j < lv.size(); ++j) { Sz v = lv[j]; if (v >= mN) continue;
#endif
#ifdef VTXWATCH
        if((long long)v==(VTXWATCH)) std::fprintf(stderr,"[NODE-SETTLE v=%lu B=%lu oldRole=%d->%d span=%lld @mTNow=%lld | mBParent_v=%ld outer_v=%lu]\n",(unsigned long)v,(unsigned long)B,mNRole[B],role,(long long)span,(long long)mTNow,(long)mBParent[v],(unsigned long)outerNode(v));
#endif
        mYAcc[v] += d * span;
      }
    }
  }
#endif
  if (mNRole[B] == LBL_EVEN) mZAcc[B] += 2 * (mTNow - mNStamp[B]);
  else if (mNRole[B] == LBL_ODD) mZAcc[B] -= 2 * (mTNow - mNStamp[B]);
  mNRole[B] = role; mNStamp[B] = mTNow;
}
// sync every vertex's/node's lazy role to its CURRENT materialized label (funnel through
// setVtx/NodeRole so the prior open interval is settled into Acc). Called right before a
// dual step in the scan engine: roles are stable from after-grow to the step, so settling
// them now and then advancing mTNow reproduces the materialized y/z increment exactly.
inline void ScalingS5Engine::lazySyncRoles() {
#ifdef S5DBG
  gS5SyncWork += (long long)(mN + (mNumNodes - mN));
  { extern long long gS5SyncCalls; ++gS5SyncCalls; }
#endif
#ifdef OUTERTAG
  { extern long long gOTsync, gOTsyncCalls; gOTsync -= gS5OuterCalls; ++gOTsyncCalls; }
#endif
#ifdef PS5_SYNC_DIRTY
#ifdef S5DBG
  { extern long long gS5SyncFull, gS5SyncIncr; if (mSyncDirtyAll) ++gS5SyncFull; else ++gS5SyncIncr; }
#endif
  if (!mSyncDirtyAll) {
    // INCREMENTAL: only nodes whose label changed since the last sync are dirty. Sync each dirty node's
    // leaves to that node's CURRENT outer-node label (a leaf's role = mLabel of its current OUTER node, so
    // resolve via outerNode of the leaf -- the dirty node may itself now be inside a bigger blossom). The
    // node's own z-role is settled in the node loop below. Idempotent: the mVRole[v]!=r guard makes any
    // stale/duplicate dirty entry a no-op. Clear after processing.
    for (Sz di = 0; di < mSyncDirty.size(); ++di) {
      Sz D = mSyncDirty[di];
      if (D >= mNodeLeaves.size()) continue;
      const std::vector<Sz>& lv = mNodeLeaves[D];
      for (Sz j = 0; j < lv.size(); ++j) {
        Sz v = lv[j]; if (v >= mN) continue;
        Sz V = outerNode(v); int r = mLabel[V];
        if (mVRole[v] != r) {
#ifdef S5DBG
          { extern long long gS5SyncChanged; ++gS5SyncChanged; }
#endif
          setVtxRole(v, r);
        }
      }
    }
    mSyncDirty.clear();
    // node z-roles: still cheap (nodes << n typically); keep the full node loop (correct + not the hot cost).
    for (Sz b = mN; b < mNumNodes; ++b) {
      int r = (mBParent[b] == NIL()) ? mLabel[b] : LBL_FREE;
      if (mNRole[b] != r) setNodeRole(b, r);
    }
    return;
  }
#endif
#ifdef PS5_FLIP_DIRTY
  // EVENT-DRIVEN INCREMENTAL (the persistent-clock enabler): process ONLY the nodes whose label flipped
  // since the last sync (pushed via setLabel), not the O(n) scan. Each flipped TOP-LEVEL entity is synced
  // via setVtxRole (bare) / setNodeRole (blossom) -- which under node-offset y settles members lazily,
  // O(size) ONLY on the real flip (the unavoidable dual work), and ZERO scan cost for the ~99% unchanged.
  // A flipped node that is now a MEMBER (mBParent!=NIL) is FREE-frozen -> setNodeRole(.,FREE). mFlipAll
  // (safety) falls back to the full scan below.
  if (!mFlipAll) {
    for (Sz fi = 0; fi < mFlipNodes.size(); ++fi) {
      Sz X = mFlipNodes[fi];
      if (X < mN) {                                  // bare vertex node
        if (mBParent[X] != NIL()) continue;          // it became a blossom member -> handled via its blossom
        int r = mLabel[X];
        if (mVRole[X] != r) setVtxRole(X, r);
      } else {                                       // blossom node
        int r = (mBParent[X] == NIL()) ? mLabel[X] : LBL_FREE;
        if (mNRole[X] != r) setNodeRole(X, r);
      }
    }
    mFlipNodes.clear();
    return;
  }
  mFlipNodes.clear();   // full scan below covers everything; reset the list + re-arm incremental after
  mFlipAll = false;
#endif
#ifdef PS5_NODE_Y
  // FAITHFUL §5.2: a MEMBER's y-rate comes from its OUTER BLOSSOM's role (read in yLazy), not its own
  // mVRole. So we only sync the role of TOP-LEVEL entities: bare vertices (loop 1) and blossom nodes
  // (loop 2). When a blossom's role changes we FIRST settle every member's open y-interval at the OLD
  // role (freeze into mYAcc) THEN restamp -- exactly what setVtxRole does for a bare vertex, done via
  // the blossom so it is O(members-of-nodes-that-changed) not O(n). This is the whole cliff fix: on a
  // stable giant blossom (role unchanged) this touches ZERO of its members.
  for (Sz v = 0; v < mN; ++v) {
    if (mBParent[v] != NIL()) continue;          // member: role read through its blossom, skip
    int r = mLabel[v];
    if (mVRole[v] != r) {
#ifdef S5DBG
      { extern long long gS5SyncChanged; ++gS5SyncChanged; }
#endif
      setVtxRole(v, r);
    }
  }
  for (Sz b = mN; b < mNumNodes; ++b) {
    int r = (mBParent[b] == NIL()) ? mLabel[b] : LBL_FREE;   // non-root blossom is frozen (FREE)
    if (mNRole[b] != r) setNodeRole(b, r);   // setNodeRole settles this node's members' y AND z, then restamps
  }
  return;
#endif
  for (Sz v = 0; v < mN; ++v) { Sz V = outerNode(v); int r = mLabel[V]; if (mVRole[v] != r) {
#ifdef S5DBG
    { extern long long gS5SyncChanged; ++gS5SyncChanged; }
#endif
    setVtxRole(v, r); } }
#ifdef OUTERTAG
  { extern long long gOTsync; gOTsync += gS5OuterCalls; }
#endif
  for (Sz b = mN; b < mNumNodes; ++b) {
    int r = (mBParent[b] == NIL()) ? mLabel[b] : LBL_FREE;   // non-root z is frozen
    if (mNRole[b] != r) setNodeRole(b, r);
  }
}
// NO-RELABEL augment: settle to FREE every outer node in the tree rooted at `root` (its
// lazy y/z freeze; it leaves the active forest). Called BEFORE the flip so mRoot is still
// valid. Stale PQ events of these nodes die on pop (treeLive false once base is matched).
inline void ScalingS5Engine::freeConsumedTree(Sz root) {
  if (root == NIL()) return;
#ifdef PS5_SYNC_DIRTY
  mSyncDirtyAll = true;   // frees a whole tree's outer nodes -> many leaves change; force full sync next.
#endif
  for (Sz i = 0; i < mNumNodes; ++i) {
    if (mBParent[i] != NIL()) continue;          // only outer nodes
    if (mRoot[i] != root) continue;
    // settle every leaf vertex's y role + the node's z role to FREE (freeze duals).
#ifdef PS5_LEAF_SPLICE
    MBX_LF_FOR(v9, i) setVtxRole(v9, LBL_FREE);
#else
    std::vector<Sz>& lv = mNodeLeaves[i];
    for (Sz j = 0; j < lv.size(); ++j) setVtxRole(lv[j], LBL_FREE);
#endif
    if (i >= mN) setNodeRole(i, LBL_FREE);
    mLabel[i] = LBL_FREE; mParent[i] = NIL(); mParEU[i] = NIL(); mParEV[i] = NIL(); mRoot[i] = NIL();
  }
}
#if defined(PS5_PERSIST_FOREST) && defined(PS5_PERSIST_CLOCK)
// STAGE 2 BATCH RETIREMENT (docs/H_Scope_PersistentClock.md): at a persistent phase boundary,
// settle EVERY dead tree (root's base got matched this phase) in ONE O(numNodes)+O(dead edges)
// pass -- the per-augment freeConsumedTree was an O(n)-scan-per-call quadratic (84% of the run).
// Correct because time does not advance WITHIN a phase: roles only matter when the clock moves,
// and this runs at the boundary, before any advance. For each dead node: settle leaf y-roles and
// the node z-role to FREE, then RESCAN its boundary edges from the surviving EVEN side at the new
// FREE rate (old events fire early and die at pop-validation; without fresh events the TRUE later
// fire time would be missed = clock overrun). Dead-to-dead edges are skipped (no event needed);
// a later-in-pass dead target still labeled EVEN is skipped via treeLive.
inline void ScalingS5Engine::retireDeadTrees() {
  LL fcTgt = (mCrit == CRIT1) ? 0 : -2;
  for (Sz i = 0; i < mNumNodes; ++i) {
    if (mBParent[i] != NIL()) continue;
    if (mRoot[i] == NIL() || treeLive(i)) continue;
#ifdef PS5_EVTRACE
    std::fprintf(stderr,"[TR] RETIRE node=%lu (root=%ld lbl=%d) mTNow=%lld\n",
      (unsigned long)i,(long)mRoot[i],mLabel[i],(long long)mTNow);
#endif
#ifdef PS5_LEAF_SPLICE
    MBX_LF_FOR(v9, i) setVtxRole(v9, LBL_FREE);
    if (i >= mN) setNodeRole(i, LBL_FREE);
    MBX_LF_FOR(v, i) {
      if (v >= mN) continue;
#else
    std::vector<Sz>& lv = mNodeLeaves[i];
    for (Sz j = 0; j < lv.size(); ++j) setVtxRole(lv[j], LBL_FREE);
    if (i >= mN) setNodeRole(i, LBL_FREE);
    for (Sz j = 0; j < lv.size(); ++j) {
      Sz v = lv[j]; if (v >= mN) continue;
#endif
      Sz ufv = ufFind(v);
      for (Sz k = 0; k < mAdjV[v].size(); ++k) {
        Sz u2 = mAdjV[v][k];
        if (mIgnoreDummies && u2 >= mNReal) continue;
        if (ufFind(u2) == ufv) continue;               // internal edge
        Sz W2 = outerNode(u2);
        if (mLabel[W2] != LBL_EVEN || !treeLive(W2)) continue;
        LL s = slackZW(u2, v, 2 * mAdjW[v][k]);
        LL d = s - fcTgt;
        if (mAllowOvertight && d < 0) d = 0;
        if (mAllowOvertight || d >= 0) mPQ.schedule(mTNow + d, (long long)u2, (long long)v, 0);
      }
    }
    mLabel[i] = LBL_FREE; mParent[i] = NIL(); mParEU[i] = NIL(); mParEV[i] = NIL(); mRoot[i] = NIL();
  }
}
#endif
// mark every leaf vertex on the tree path startEven..root as USED this phase (per-level
// disjointness: no later augment THIS LEVEL may reuse a consumed vertex -> progress
// guarantee + correctness of multi-augment within one dual level). Walk via treeParent.
inline void ScalingS5Engine::markPathUsed(Sz startEven) {
  for (Sz x = startEven; ; ) {
#ifdef PS5_LEAF_SPLICE
    MBX_LF_FOR(v9, x) mUsed[v9] = 1;
#else
    std::vector<Sz>& lv = mNodeLeaves[x];
    for (Sz j = 0; j < lv.size(); ++j) mUsed[lv[j]] = 1;
#endif
    Sz p = treeParent(x);
    if (p == NIL()) break;
    x = p;
  }
}
// True full-path disjointness (the real HK phase invariant): a multi-augment this phase is only
// correct if its ENTIRE augmenting path is free of vertices consumed by an earlier augment this
// phase. The endpoint-only mUsed check (mUsed[u]||mUsed[w]) is insufficient: a path whose endpoints
// are fresh can still traverse an INTERIOR node an earlier augment re-mated -> the tree-path then
// encodes a stale matched edge -> a slightly-suboptimal augment (the off-by-~one-edge bug at scale).
// If any path node is used, defer this augment to after the phase relabel (return false -> skip).
inline bool ScalingS5Engine::pathClean(Sz startEven) const {
  for (Sz x = startEven; ; ) {
#ifdef PS5_LEAF_SPLICE
    MBX_LF_FOR(v9, x) if (mUsed[v9]) return false;
#else
    const std::vector<Sz>& lv = mNodeLeaves[x];
    for (Sz j = 0; j < lv.size(); ++j) if (mUsed[lv[j]]) return false;
#endif
    Sz p = treeParent(x);
    if (p == NIL()) break;
    x = p;
  }
  return true;
}
// VERIFY (debug): assert lazy y/z exactly equal the materialized duals. Used in the
// scan engine (which keeps mY/mZ materialized) to prove the lazy scheme before the
// event engine relies on it. Aborts with a diff on mismatch.
inline void ScalingS5Engine::lazyVerifyStep() {
  for (Sz v = 0; v < mN; ++v)
    if (yLazy(v) != mY[v]) {
      std::fprintf(stderr, "LAZYVERIFY y MISMATCH v=%lu lazy=%lld mat=%lld t=%lld\n",
                   (unsigned long)v, (long long)yLazy(v), (long long)mY[v], (long long)mTNow);
      std::abort();
    }
  for (Sz b = mN; b < mNumNodes; ++b)
    if (mBParent[b] == NIL() && zLazy(b) != mZ[b]) {
      std::fprintf(stderr, "LAZYVERIFY z MISMATCH B=%lu lazy=%lld mat=%lld t=%lld\n",
                   (unsigned long)b, (long long)zLazy(b), (long long)mZ[b], (long long)mTNow);
      std::abort();
    }
}

// ---------------------------------------------------------------------------
// node allocation + capacity
// ---------------------------------------------------------------------------
inline void ScalingS5Engine::ensureCap(Sz needNode) {
  if (needNode < mZ.size()) return;
  Sz ncap = 2 * needNode + 16;
  mZ.resize(ncap, 0);
  mBParent.resize(ncap, NIL());
  mBBase.resize(ncap, NIL());
  mBCycle.resize(ncap);
  mBCycEU.resize(ncap);
  mBCycEV.resize(ncap);
  mNodeLeaves.resize(ncap);
  mLabel.resize(ncap, LBL_FREE);
  mParent.resize(ncap, NIL());
  mParEU.resize(ncap, NIL());
  mParEV.resize(ncap, NIL());
  mRoot.resize(ncap, NIL());
  // keep the lazy z-arrays in lockstep so a mid-search no-relabel shrink (new blossom
  // past the current cap) can write mZ0[B]/etc. A new node's lazy z is z0=0/Acc=0/FREE.
  if (!mZ0.empty()) {
    mZ0.resize(ncap, 0); mZAcc.resize(ncap, 0); mNRole.resize(ncap, LBL_FREE); mNStamp.resize(ncap, 0);
  }
}

inline ScalingS5Engine::Sz ScalingS5Engine::newBlossom() {
  Sz B;
  if (!mFreeBlossoms.empty()) {            // reclaim a dissolved id (bounds memory:
    B = mFreeBlossoms.back();              // O(live blossoms) ids, not O(ever-formed)
    mFreeBlossoms.pop_back();
  } else {
    B = mNumNodes++;
    ensureCap(B);
  }
#ifdef PS5_ANCESTOR_HINT
  ahEnsure(); if (B < mAhDead.size()) mAhDead[B] = 0;   // (re)issued id starts alive
#endif
  mZ[B] = 0; mBParent[B] = NIL();
#ifdef PS5_LEAF_LIST
  // fresh id (possibly reclaimed): reset ITS list keys only -- the segment CUT belongs to the
  // real dissolve() (a misplaced cut here ran against the PREVIOUS incarnation's stale mBCycle
  // = the n=28 cross-contamination).
  if (B < mLfHead.size()) { mLfHead[B] = NIL(); mLfTail[B] = NIL(); mLfCnt[B] = 0; }
#endif
  mBCycle[B].clear(); mBCycEU[B].clear(); mBCycEV[B].clear(); mNodeLeaves[B].clear();
  mLabel[B] = LBL_FREE; mParent[B] = NIL(); mParEU[B] = NIL(); mParEV[B] = NIL(); mRoot[B] = NIL();
#ifdef PS5_NODE_Y
  // A reclaimed id keeps its PREVIOUS life's lazy role/stamp (newBlossom didn't reset them). Under
  // node-offset y that stale ODD/EVEN role makes setNodeRole(B,...) settle members at a bogus rate.
  // Reset to a clean FREE base here; the caller sets the real role via setNodeRole after building B.
  // Guard on !empty(): the lazy vectors are only allocated once the lazy path begins (lazyBegin);
  // before that they are size 0 and B would be OOB (the testS5LazyEvent segfault). ensureCap keeps
  // them sized to capacity when the reclaim path grows past the old cap.
  if (B < mNRole.size()) { mNRole[B] = LBL_FREE; mNStamp[B] = mTNow; mZ0[B] = 0; mZAcc[B] = 0; }
#endif
  return B;
}

// ---------------------------------------------------------------------------
// LCA of outer nodes in the alternating forest (via mParent). NIL if diff trees.
// ---------------------------------------------------------------------------
inline ScalingS5Engine::Sz ScalingS5Engine::lcaOuter(Sz a, Sz b) {
#ifdef OUTERTAG
  { extern long long gOTlca, gOTlcaCalls; gOTlca -= gS5OuterCalls; ++gOTlcaCalls; }
#endif
  std::vector<LL>& seen = mLcaSeen; LL& tok = mLcaTok;   // instance state (was static)
  if (seen.size() < mNumNodes) seen.assign(mNumNodes, 0);
  ++tok;
#ifdef S5DBG
  Sz g1 = 0, g2 = 0;
#endif
  // navigate the CONTRACTED tree via treeParent (resolves through absorbed nodes).
  for (Sz x = a; ; ) { seen[x] = tok; Sz p = treeParent(x); if (p == NIL()) break; x = p;
#ifdef S5DBG
    if (++g1 > mNumNodes + 5) { std::fprintf(stderr,"  [lcaOuter CYCLE a-chain from %lu]\n",(unsigned long)a); std::fflush(stderr); return NIL(); }
#endif
  }
  for (Sz y = b; ; ) { if (seen[y] == tok) {
#ifdef OUTERTAG
      { extern long long gOTlca; gOTlca += gS5OuterCalls; }
#endif
      return y; } Sz p = treeParent(y); if (p == NIL()) break; y = p;
#ifdef S5DBG
    if (++g2 > mNumNodes + 5) { std::fprintf(stderr,"  [lcaOuter CYCLE b-chain from %lu]\n",(unsigned long)b); std::fflush(stderr); return NIL(); }
#endif
  }
#ifdef OUTERTAG
  { extern long long gOTlca; gOTlca += gS5OuterCalls; }
#endif
  return NIL();
}

// ---------------------------------------------------------------------------
// rotateToBase: rematch around blossom B so the sub-node carrying baseVtx becomes
// the cycle base (index 0), flipping the even-length segment. (Ported from the
// proven reference engine, adapted to this engine's members.)
// ---------------------------------------------------------------------------
inline void ScalingS5Engine::rotateToBase(Sz B, Sz baseVtx) {
  if (B < mN) return;
  std::vector<Sz>& cyc = mBCycle[B];
  std::vector<Sz>& eu  = mBCycEU[B];
  std::vector<Sz>& ev  = mBCycEV[B];
  Sz k = cyc.size();
#ifdef PS5_ROTATE_XFER
  // NEGATIVE RESULT (2026-07-02): M==E root cause IS traced here -- the internal flip below re-pairs the OLD
  // base vertex to an internal partner, ORPHANING the old base's EXTERNAL matched edge (extMate->oldBase is
  // left dangling asymmetric); that dangling half-edge is the asymmetric mate the later grow step trips over
  // as M==E (LEMON+MATE_AUDIT, n=54 seed=30 W=5: mate[49]=37 but mate[37]=13). BUT neither repair works here:
  //   - TRANSFER extMate to the new base (PS5_ROTATE_XFER_MOVE) -> matching CYCLE + infinite search loop (366s hang).
  //   - CLEAR the dangling half -> corrupts the matching (gate produces no/garbage output).
  // => rotateToBase is the WRONG place: the augment/intact-pole callers already re-route the base's external
  //    edge themselves, so repairing it here double-writes. The fix must be at the CONSUMER (the grow step
  //    that reads the stale mate) or by making the flip external-match-aware per caller-context. Left gated
  //    OFF (default build unaffected). Do NOT re-enable without solving the caller-context problem.
  Sz oldBaseVtx = mBBase[B];
  Sz extMate = (oldBaseVtx < mN) ? mMate[oldBaseVtx] : NIL();
  bool extOutside = (extMate != NIL() && !leafOfNode(extMate, B));
#endif
#ifdef DIVDBG
  // entry/exit alignment self-check for THIS blossom only (first violation wins).
  struct _RotChk {
    ScalingS5Engine* e; Sz B; bool entryOk;
    bool aligned() const {
      const std::vector<Sz>& c = e->mBCycle[B]; Sz kk = c.size();
      if (kk == 0) return true;
      for (Sz i = 0; i < kk; ++i) {
        Sz a = e->mBCycEU[B][i], b = e->mBCycEV[B][i];
        Sz ca = a; Sz g = 0;
        while (ca != NIL() && e->mBParent[ca] != B && g++ <= e->mNumNodes) ca = e->mBParent[ca];
        Sz cb = b; g = 0;
        while (cb != NIL() && e->mBParent[cb] != B && g++ <= e->mNumNodes) cb = e->mBParent[cb];
        if (ca != c[i] || cb != c[(i + 1) % kk]) return false;
      }
      return true;
    }
    _RotChk(ScalingS5Engine* e_, Sz B_) : e(e_), B(B_) { entryOk = aligned(); }
    ~_RotChk() {
      static bool s_rp = false;
      if (!s_rp && entryOk && !aligned()) { s_rp = true;
        std::fprintf(stderr,"[ROT-BROKE-ALIGN B=%lu k=%lu base=%lu]\n",
          (unsigned long)B,(unsigned long)e->mBCycle[B].size(),(unsigned long)e->mBBase[B]);
        std::fflush(stderr); }
      static bool s_re = false;
      if (!s_re && !entryOk) { s_re = true;
        std::fprintf(stderr,"[ROT-ENTRY-ALREADY-BAD B=%lu k=%lu]\n",
          (unsigned long)B,(unsigned long)e->mBCycle[B].size()); std::fflush(stderr); }
    }
  } _rotchk(this, B);
#endif
  Sz target = baseVtx;
  while (target != NIL() && mBParent[target] != B) target = mBParent[target];
  Sz pos = 0; while (pos < k && cyc[pos] != target) ++pos;
  if (pos >= k) {
#ifdef DIVDBG
    std::fprintf(stderr,"[ROTATE-MISS B=%lu baseVtx=%lu target=%ld curBase=%lu k=%lu bp(base)=%ld outer(base)=%lu uf(base)=%lu]\n",
      (unsigned long)B,(unsigned long)baseVtx,(long)target,(unsigned long)mBBase[B],(unsigned long)k,
      (long)(baseVtx < mBParent.size() ? (long)mBParent[baseVtx] : -9),
      (unsigned long)((ScalingS5Engine*)this)->outerNode(baseVtx),
      (unsigned long)(baseVtx < mN ? ((ScalingS5Engine*)this)->ufFind(baseVtx) : (Sz)0));
    std::fflush(stderr);
#ifdef DIVDBG_ABORT
    std::abort();
#endif
#endif
    return;
  }
  if (pos % 2 == 1) {
    std::vector<Sz> rc(k), reu(k), rev2(k);
    rc[0] = cyc[0];
    for (Sz j = 1; j < k; ++j) rc[j] = cyc[k - j];
    for (Sz j = 0; j < k; ++j) {
      Sz oe = (k - 1 - j) % k;
      reu[j] = ev[oe];
      rev2[j] = eu[oe];
    }
    cyc = rc; eu = reu; ev = rev2;
    pos = 0; while (pos < k && cyc[pos] != target) ++pos;
  }
  if (pos != 0) {
    for (Sz i = 0; i + 1 <= pos; i += 2) {
      Sz a = eu[i], b = ev[i];
      if (cyc[i]   >= mN) rotateToBase(cyc[i],   a);
      if (cyc[i+1] >= mN) rotateToBase(cyc[i+1], b);
#ifdef MATE_AUDIT
      // watch for the SPECIFIC steal: pairing a<->b where a (or b) had an old mate that lies OUTSIDE this
      // blossom's leaves -> that outside vertex's mate is never cleared -> dangling. Report those only.
      { Sz oa=mMate[a], ob=mMate[b];
        bool aOutside = (oa!=NIL() && oa!=b && !leafOf(oa,B));
        bool bOutside = (ob!=NIL() && ob!=a && !leafOf(ob,B));
        static int r=0;
        if(r<8 && (aOutside||bOutside)){ ++r;
          std::fprintf(stderr,"[REMATCH-OUTSIDE-STEAL B=%lu base=%lu a=%lu b=%lu | mate[a]=%ld(outside=%d) mate[b]=%ld(outside=%d) | a==base?%d b==base?%d]\n",
            (unsigned long)B,(unsigned long)mBBase[B],(unsigned long)a,(unsigned long)b,(long)oa,(int)aOutside,(long)ob,(int)bOutside,
            (int)(a==mBBase[B]),(int)(b==mBBase[B])); } }
#endif
      MG_CHK("rot-rematch", a, b);
      mMate[a] = b; mMate[b] = a;
    }
    std::vector<Sz> nc(k), neu(k), nev(k);
    for (Sz i = 0; i < k; ++i)
      { nc[i]=cyc[(pos+i)%k]; neu[i]=eu[(pos+i)%k]; nev[i]=ev[(pos+i)%k]; }
    cyc = nc; eu = neu; ev = nev;
  }
  if (cyc[0] >= mN) rotateToBase(cyc[0], baseVtx);
  mBBase[B] = mBBase[cyc[0]];
#ifdef PS5_ROTATE_XFER
  // Heal the dangling external half-edge: if the internal flip stole the old base from its OUTSIDE partner
  // extMate (mMate[oldBaseVtx] no longer == extMate), then extMate->oldBaseVtx is now a dangling asymmetric
  // mate. TRANSFER re-points it to the new base -- but that created a matching CYCLE + infinite search loop
  // (hang) because the augment/intact-pole callers already re-route the base's external edge themselves, so
  // the new base ends up doubly-matched. Instead CLEAR the dangling half (extMate becomes free); the normal
  // search re-matches it with no cycle risk. (Transfer kept under PS5_ROTATE_XFER_MOVE for reference.)
  if (extOutside && oldBaseVtx != baseVtx && mMate[oldBaseVtx] != extMate) {
#ifdef PS5_ROTATE_XFER_MOVE
    if (baseVtx < mN) { mMate[baseVtx] = extMate; mMate[extMate] = baseVtx; }
#else
    if (extMate < mN) mMate[extMate] = NIL();
#endif
  }
#endif
}

// flip matched<->unmatched along the tree path from EVEN node startEven to its
// root, navigating the CONTRACTED tree via treeParent (so absorbed parents are
// resolved to their current outer node -- no relabel needed after shrink).
inline void ScalingS5Engine::flipUpToRoot(Sz startEven) {
  Sz E = startEven;
  for (;;) {
    Sz O = treeParent(E);
    if (O == NIL()) break;
    Sz E2 = treeParent(O);
    Sz uo = mParEU[O];
    Sz ue = mParEV[O];
    if (O  >= mN) rotateToBase(O,  uo);
    if (E2 != NIL() && E2 >= mN) rotateToBase(E2, ue);
#ifdef MATEGUARD
    if (uo >= mN || ue >= mN) {
      std::fprintf(stderr,"[MG-FLIP-CTX start=%lu E=%lu O=%lu E2=%ld | lbl(O)=%d par(O)=%ld root(O)=%ld base(O)=%ld bp(O)=%ld | O>=mN:%d]\n",
        (unsigned long)startEven,(unsigned long)E,(unsigned long)O,(long)E2,
        (O<mLabel.size()?mLabel[O]:-9),(long)(O<mParent.size()?(long)mParent[O]:-9),
        (long)(O<mRoot.size()?(long)mRoot[O]:-9),(long)(O<mBBase.size()?(long)mBBase[O]:-9),
        (long)(O<mBParent.size()?(long)mBParent[O]:-9),(int)(O>=mN));
      // full chain from startEven with labels + parent-writer tags (parw is DIVDBG-gated; '?' if absent)
      { Sz x = startEven;
        for (int st = 0; st < 12 && x != NIL(); ++st) {
          std::fprintf(stderr,"  [MG-CHAIN %d] x=%lu lbl=%d par=%ld parEU=%ld parEV=%ld w=%c bp=%ld\n",
            st,(unsigned long)x,(x<mLabel.size()?mLabel[x]:-9),
            (long)(x<mParent.size()?(long)mParent[x]:-9),
            (long)(x<mParEU.size()?(long)mParEU[x]:-9),(long)(x<mParEV.size()?(long)mParEV[x]:-9),
            (char)(mParWriter.size()>x?mParWriter[x]:'?'),
            (long)(x<mBParent.size()?(long)mBParent[x]:-9));
          x = treeParent(x);
        } }
      std::fflush(stderr);
    }
#endif
    MG_CHK("flipUpToRoot", uo, ue);
    mMate[uo] = ue; mMate[ue] = uo;
    if (E2 == NIL()) break;
    E = E2;
  }
}

// ---------------------------------------------------------------------------
// shrink: EVEN outer nodes E,W (same tree) via eligible edge (uOrig in E, vOrig in
// W). Build blossom B over the odd cycle E..lca..W; z(B)=0. Maintain union-find:
// union all cycle members into B's base vertex. (Ported + UF-augmented.)
// ---------------------------------------------------------------------------
inline void ScalingS5Engine::shrink(Sz E, Sz W, Sz uOrig, Sz vOrig) {
#ifdef TRAJ
  { LL wt=-1; for(Sz k=0;k<mAdjV[uOrig].size();++k) if(mAdjV[uOrig][k]==vOrig){wt=mAdjW[uOrig][k];break;}
    std::fprintf(stderr,"[TRAJ SHRINK (%lu,%lu) w=%lld]\n",(unsigned long)uOrig,(unsigned long)vOrig,(long long)wt); }
#endif
  Sz L = lcaOuter(E, W);
  if (L == NIL()) return;
  Sz B = newBlossom();
  std::vector<Sz> eSide, eU, eV, wSide, wU, wV;
  // Walk E->L and W->L via treeParent (CONTRACTED outer-node navigation). Raw mParent
  // can point at a now-absorbed inner sub-node (after a prior no-relabel shrink), so the
  // raw walk would never reach the outer node L and run off the end (the ASan overflow).
  // treeParent resolves each stored parent to its current OUTER node, matching lcaOuter.
  { Sz x = E; while (x != L) { eSide.push_back(x); eU.push_back(mParEU[x]); eV.push_back(mParEV[x]); x = treeParent(x); } }
  { Sz x = W; while (x != L) { wSide.push_back(x); wU.push_back(mParEU[x]); wV.push_back(mParEV[x]); x = treeParent(x); } }
  mBCycle[B].push_back(L);
  for (Sz i = eSide.size(); i-- > 0; ) {
    mBCycEU[B].push_back(eV[i]); mBCycEV[B].push_back(eU[i]);
    mBCycle[B].push_back(eSide[i]);
  }
  mBCycEU[B].push_back(uOrig); mBCycEV[B].push_back(vOrig);
  for (Sz i = 0; i < wSide.size(); ++i) {
    if (i == 0) { mBCycle[B].push_back(wSide[i]); }
    else { mBCycEU[B].push_back(wU[i-1]); mBCycEV[B].push_back(wV[i-1]); mBCycle[B].push_back(wSide[i]); }
  }
  if (!wSide.empty()) { Sz last = wSide.size()-1; mBCycEU[B].push_back(wU[last]); mBCycEV[B].push_back(wV[last]); }
  else { mBCycEU[B].push_back(mParEU[W]); mBCycEV[B].push_back(mParEV[W]); }
  mBBase[B] = mBBase[L];
#ifdef PS5_NODE_Y
  // FAITHFUL SETTLE-AT-MEMBERSHIP (§5.2, docs/J_Scope): each absorbed vertex/sub-blossom was accruing
  // y under its OLD role until NOW. Freeze that accrual into mYAcc before mBParent redirects to B, so
  // afterward each leaf reads through B fresh. ENGINE-AGNOSTIC: settle each cycle member's accrual at
  // its OWN current role here (bare vertex inline; sub-blossom via setNodeRole(sub,FREE), which settles
  // that sub's members' y + z). This does NOT rely on the lazy-event caller's cycle-clear -- setNodeRole
  // is idempotent (no-op if sub is already FREE), so a cycle-clear that already ran costs nothing, and
  // the SCAN engine (no cycle-clear) is now settled too (the LazyVerify mismatch). O(cycle+members)
  // piggybacks on the leaf-copy. B's own role/stamp are set fresh by the caller (setNodeRole(B,EVEN)).
  // Gate on mLazyActive||mLazyVerify: scan verify maintains lazy timestamps in parallel (mLazyVerify).
  if (mLazyActive || mLazyVerify) {
    for (Sz ci = 0; ci < mBCycle[B].size(); ++ci) {
      Sz sub = mBCycle[B][ci];
      if (sub < mN) {                                  // bare vertex sub == leaf itself
        if (mVRole[sub] == LBL_ODD)       mYAcc[sub] += (mTNow - mVStamp[sub]);
        else if (mVRole[sub] == LBL_EVEN) mYAcc[sub] -= (mTNow - mVStamp[sub]);
        mVRole[sub] = LBL_FREE; mVStamp[sub] = mTNow;  // absorbed: its own role no longer accrues
      } else {
        setNodeRole(sub, LBL_FREE);                    // settles sub's members' y + z; idempotent if already FREE
      }
    }
  }
#endif
  for (Sz i = 0; i < mBCycle[B].size(); ++i) mBParent[mBCycle[B][i]] = B;
#ifdef PS5_CQ_EVENTS
  // STAGE 1b: fuse the cycle sub-blossoms' membership sublists into B's (one rep leaf per cycle member).
  // RECORD the split pivots NOW (each rep's sublist rightmost vtx, in merge order) -- replayed verbatim at
  // dissolve, so the split is immune to rotateToBase reordering mBCycle between shrink and dissolve.
  if (!mCqOrder.empty()) {
    if (mCqPivots.size() <= B) mCqPivots.resize(B + 1);
    mCqPivots[B].clear();
    std::vector<Sz> reps; reps.reserve(mBCycle[B].size());
    for (Sz i = 0; i < mBCycle[B].size(); ++i) {
      Sz sub = mBCycle[B][i];
      // rep MUST be a vertex genuinely in sub's sublist. mNodeLeaves[sub][0] is STALE (append-only) ->
      // findQueue on it can return a FOREIGN sublist and over-merge (the n=120 shrink bug). Use the
      // AUTHORITATIVE base vertex mBBase[sub] (always a real leaf of sub); bare vertex = itself.
      Sz rep = (sub < mN) ? sub : mBBase[sub];
      if (rep != NIL() && rep < mN) reps.push_back(rep);
    }
#ifdef PS5_CQ_REINIT_DISSOLVE
    (void)reps;   // DIAGNOSTIC: rebuild instead of incremental merge (same isolation as dissolve).
    cqInitMembership();
#else
    if (reps.size() > 1) mSFM.merge(reps);
#endif
#ifdef PS5_CQ_CHECK_EVERY_OP
    cqCheckMembership("shrink");
#endif
  }
#endif
#ifdef PS5_OUTERCACHE
  outerCacheBump();   // mBParent changed -> invalidate outerNode cache
#endif
#ifdef S5DBG
  ++gS5ShrinkCalls;
#endif
  // CLIFF FIX #3 (2026-07-02): rebuild mNodeLeaves[B] with ONE up-front reserve (pre-summed size) so the
  // append does zero re-allocation / memmove (callgrind: memmove+malloc were ~13% of the run, the last
  // super-linear cliff term after the relink fixes -- shrinkCopy grew ~7.6x/doubling via realloc churn).
#ifndef PS5_LEAF_NOVEC
  { Sz tot = 0; for (Sz i = 0; i < mBCycle[B].size(); ++i) tot += mNodeLeaves[mBCycle[B][i]].size();
    mNodeLeaves[B].reserve(tot); }
  for (Sz i = 0; i < mBCycle[B].size(); ++i) {
    std::vector<Sz>& cl = mNodeLeaves[mBCycle[B][i]];
#ifdef S5DBG
    gS5ShrinkLeafCopy += (long long)cl.size();
#endif
    mNodeLeaves[B].insert(mNodeLeaves[B].end(), cl.begin(), cl.end());
  }
#endif
#ifdef PS5_LEAF_LIST
  { // splice children's segments in the SAME cycle order as the vector copy above
    Sz h = NIL(), t = NIL(), cnt = 0;
    for (Sz i = 0; i < mBCycle[B].size(); ++i) { Sz c = mBCycle[B][i];
      if (c >= mLfHead.size() || mLfHead[c] == NIL()) continue;
      { // PRE-SPLICE segment validation: c's (head,tail) must be a clean standalone segment --
        // head has no prev, tail has no next (c was top-level), and the walk head->tail must
        // terminate with exactly mLfCnt[c] steps. A violation = the engine absorbed a child
        // that was NOT truly top-level (its vertices still live inside another node's chain).
        Sz g = 0, w2 = mLfHead[c]; bool ok = (mLfPrev[mLfHead[c]] == NIL());
        while (w2 != NIL() && g <= mN) { ++g; if (w2 == mLfTail[c]) break; w2 = mLfNext[w2]; }
        if (!ok || w2 != mLfTail[c] || g != mLfCnt[c] || mLfNext[mLfTail[c]] != NIL())
          std::fprintf(stderr,
            "[LF-SPLICE-BADCHILD B=%lu child=%lu walked=%lu cnt=%lu headPrev=%ld tailNext=%ld bp=%ld]\n",
            (unsigned long)B,(unsigned long)c,(unsigned long)g,(unsigned long)mLfCnt[c],
            (long)mLfPrev[mLfHead[c]],(long)mLfNext[mLfTail[c]],(long)mBParent[c]);
      }
      if (h == NIL()) { h = mLfHead[c]; t = mLfTail[c]; }
      else { mLfNext[t] = mLfHead[c]; mLfPrev[mLfHead[c]] = t; t = mLfTail[c]; }
      cnt += mLfCnt[c];
    }
    if (B >= mLfHead.size()) { mLfHead.resize(B+8,NIL()); mLfTail.resize(B+8,NIL()); mLfCnt.resize(B+8,0); }
    mLfHead[B] = h; mLfTail[B] = t; mLfCnt[B] = cnt;
#ifdef LF_OPTRACE
    { std::fprintf(stderr,"[LF-OP shrink B=%lu cyc=",(unsigned long)B);
      for (Sz i=0;i<mBCycle[B].size();++i) std::fprintf(stderr,"%lu ",(unsigned long)mBCycle[B][i]);
      std::fprintf(stderr,"]");
      std::fprintf(stderr,"%c",10); }
#endif
    lfCheckNode(B, "shrink");
    lfSweepAll("shrink-global");
  }
#endif
#ifdef S5DBG
  { extern long long gS5MaxBlossom, gS5SubBlossomShrinks, gS5MaxDistinct;
    if ((long long)mNodeLeaves[B].size() > gS5MaxBlossom) {
      gS5MaxBlossom = (long long)mNodeLeaves[B].size();
      // DISTINCT check at the peak: is this really that many distinct vertices, or leaf-list dups?
      std::vector<char> seen(mN, 0); long long distinct=0;
      for (Sz i=0;i<mNodeLeaves[B].size();++i){ Sz v=mNodeLeaves[B][i]; if(v<mN&&!seen[v]){seen[v]=1;++distinct;} }
      gS5MaxDistinct = distinct;
    }
    for (Sz i = 0; i < mBCycle[B].size(); ++i) if (mBCycle[B][i] >= mN) { ++gS5SubBlossomShrinks; break; } }
#endif
  // union-find: merge all leaves of B into B's base vertex representative. NOTE: the per-leaf ufFind here
  // ALSO path-compresses every leaf->baseRep (not just merges) -- unioning only the k member reps instead
  // was measured SLOWER + changed the search path (more shrinks), so the per-leaf sweep stays. The residual
  // cliff is the Theta(n^2) leaf VOLUME (blossom grows to ~n via ~n shrinks); needs a lazy leaf structure.
  Sz baseRep = ufFind(mBBase[B]);
#ifdef PS5_LEAF_SPLICE
  { Sz _e = mLfTail[B];
    for (Sz v = mLfHead[B], nx = NIL(); v != NIL(); v = nx) {
      nx = (v == _e) ? NIL() : mLfNext[v];
      mUf[ufFind(v)] = baseRep;
    } }
#else
  for (Sz i = 0; i < mNodeLeaves[B].size(); ++i) mUf[ufFind(mNodeLeaves[B][i])] = baseRep;
#endif
  mUf[baseRep] = baseRep;
  // B inherits L's tree role (EVEN, same parent/parent-edge/root as L).
  mLabel[B] = LBL_EVEN;
  markLabelDirty(B);   // shrink: all of B's leaves now take B's EVEN role (covers the absorbed ODD ones)
  Sz lp = mParent[L];
  Sz lpOuter = (lp == NIL()) ? NIL() : outerNode(lp);
  if (lpOuter == B) lpOuter = NIL();
  mParent[B] = lpOuter;
  mParEU[B] = mParEU[L]; mParEV[B] = mParEV[L];
#ifdef DIVDBG
  parw(B,'S');
#endif
  mRoot[B] = (lpOuter == NIL()) ? B : mRoot[L];
#ifdef DIVDBG
  // BORN-BAD check on THIS blossom, no live filter: lift each recorded endpoint; flag misalignment
  // at construction (stale mParEU/mParEV from a prior re-thread would show here).
  { Sz kk = mBCycle[B].size(); static bool s_bb = false;
    for (Sz i = 0; i < kk && !s_bb; ++i) {
      Sz a = mBCycEU[B][i], b = mBCycEV[B][i];
      Sz ca = a; Sz g = 0;
      while (ca != NIL() && mBParent[ca] != B && g++ <= mNumNodes) ca = mBParent[ca];
      Sz cb = b; g = 0;
      while (cb != NIL() && mBParent[cb] != B && g++ <= mNumNodes) cb = mBParent[cb];
      if (ca != mBCycle[B][i] || cb != mBCycle[B][(i + 1) % kk]) { s_bb = true;
        std::fprintf(stderr,"[SHRINK-BORN-BAD B=%lu k=%lu idx=%lu cyc[i]=%lu eu=%lu(lift %ld) ev=%lu(lift %ld) E=%lu W=%lu u=%lu v=%lu]\n",
          (unsigned long)B,(unsigned long)kk,(unsigned long)i,(unsigned long)mBCycle[B][i],
          (unsigned long)a,(long)ca,(unsigned long)b,(long)cb,
          (unsigned long)E,(unsigned long)W,(unsigned long)uOrig,(unsigned long)vOrig);
        std::fflush(stderr);
        for (Sz j = 0; j < kk; ++j) { Sz cj = mBCycle[B][j];
          Sz pv = mParEV[cj]; Sz lc = pv; Sz g2 = 0;
          while (lc != NIL() && mBParent[lc] != B && g2++ <= mNumNodes) lc = mBParent[lc];
          std::fprintf(stderr,"  cyc[%lu]=%lu writer=%c rawPar=%ld parEV=%ld lift(parEV)=%ld\n",
            (unsigned long)j,(unsigned long)cj,(char)(mParWriter.size()>cj?mParWriter[cj]:63),
            (long)mParent[cj],(long)pv,(long)lc); }
        std::fflush(stderr); } } }
  { Sz bi9=0; Sz mB9 = findCycEdgeMisalign(&bi9); static bool s_r9=false;
    if (mB9 != NIL() && !s_r9) { s_r9=true; std::fprintf(stderr,"[MISALIGN-END-SHRINK B=%lu idx=%lu k=%lu]\n",
      (unsigned long)mB9,(unsigned long)bi9,(unsigned long)mBCycle[mB9].size()); std::fflush(stderr); } }
#endif
}

// S5-C STEP 4: re-thread the alternating tree through ODD blossom B's cycle in place, so a
// live-ODD dissolution needs NO relabel. Edmonds odd-blossom expansion: B entered the tree at
// its ENTRY sub-node (the one holding mParEU[B], reached from B's parent by an unmatched edge);
// the base sub-node (mBCycle[B][0]) is matched downward to B's child. The even-length arc from
// entry to base alternates ODD(entry),EVEN,...,ODD(base), linked by the cycle's matched edges;
// off-arc sub-nodes pair up matched-internally and leave the tree (FREE). Sets mLabel/mParent/
// mParEU/mParEV/mRoot for the on-arc nodes; mLazySyncRoles (caller) then settles the lazy duals.
// Returns false (caller falls back to relabel) if not the clean expandable case. Mirrors the
// proven ScalingMatching::relinkOddCycle (n=180 nested-z0 guard included).
// Non-mutating pre-check: would relinkOddCycleS5(B) succeed? Same bail conditions, no writes.
// Used for the all-or-nothing batch decision (a partial in-place batch corrupts state).
#ifdef RELINKPROBE
extern long long gRPcalls,gRPbadK,gRPnested,gRPnoEntry,gRPmultiChild,gRPok,gRPcycSum,gRPnestSubSum;
#endif
// Structural-only eligibility for a NESTED sub-blossom (no parent/entry checks -- those are done at
// relink time once the parent edge is known). A blossom's cycle is re-threadable iff it is odd-length
// >=3 and every member is either a bare vertex or a z=0 blossom that is itself struct-eligible.
inline bool ScalingS5Engine::relinkStructEligible(Sz B, Sz depth) {
#ifdef SE_PROBE
  extern long long gSEdepth, gSEbadk, gSEdeepfail, gSEok, gSEmaxdepth, gSEcalls;
  ++gSEcalls; if ((long long)depth > gSEmaxdepth) gSEmaxdepth = (long long)depth;
#endif
#ifdef PS5_DEEP_INTACT
  // MEMO: each blossom's struct-eligibility is invariant within one dissolution batch -- check once.
  if (mSEMemo.size() > B && mSEMemoStamp[B] == mSEMemoTok) return mSEMemo[B] == 1;
#endif
#ifdef PS5_DEEP_INTACT
  if (depth > 4000000) {                           // effectively-uncapped (spin-safety only; memo makes
                                                   // deep recursion cheap -- each blossom checked once)
#else
  if (depth > 8) {                                 // depth guard (data: shallow, ~1 sub-blossom)
#endif
#ifdef SE_PROBE
    ++gSEdepth;
#endif
    return false; }
  const std::vector<Sz>& cyc = mBCycle[B];
  Sz k = cyc.size();
#ifdef SE_PROBE
  if (k < 3 || (k % 2) == 0) { ++gSEbadk; return false; }
#endif
#ifdef RELINKPROBE
  if (k < 3 || (k % 2) == 0) { ++gRPnoEntry; return false; }   // REUSE gRPnoEntry = struct-fail: bad k
#else
  if (k < 3 || (k % 2) == 0) return false;
#endif
  for (Sz i = 0; i < k; ++i) {
    Sz sn = cyc[i];
    if (sn >= mN) {
#ifdef PS5_NESTED_INTACT
      if (zLazy(sn) != 0) continue;            // intact sub-blossom: rotated to base on-arc -> OK
#else
      if (zLazy(sn) != 0) return false;        // intact sub-blossom -> unsafe
#endif
      if (!relinkStructEligible(sn, depth + 1)) {
#ifdef SE_PROBE
        ++gSEdeepfail;
#endif
#ifdef PS5_DEEP_INTACT
        if (mSEMemo.size() > B) { mSEMemo[B] = 2; mSEMemoStamp[B] = mSEMemoTok; }
#endif
        return false; }
    }
  }
#ifdef SE_PROBE
  ++gSEok;
#endif
#ifdef PS5_DEEP_INTACT
  if (mSEMemo.size() > B) { mSEMemo[B] = 1; mSEMemoStamp[B] = mSEMemoTok; }
#endif
  return true;
}
inline bool ScalingS5Engine::relinkEligibleS5(Sz B) {
  const std::vector<Sz>& cyc = mBCycle[B];
  Sz k = cyc.size();
#ifdef RELINKPROBE
  ++gRPcalls; gRPcycSum += k;
  if (k < 3 || (k % 2) == 0) { ++gRPbadK; return false; }
#else
  if (k < 3 || (k % 2) == 0) return false;
#endif
  // NESTED case (measured: 99% of uniform-weight bails; avg 1 sub-blossom/cycle, shallow). A cycle
  // member that is itself a blossom is OK to re-thread IF it is z=0 (dissolving this batch too) AND its
  // own cycle is STRUCTURALLY re-threadable (odd length, members likewise). A z!=0 sub-blossom (stays
  // intact) is the dangerous case the old guard avoided -> still bail. The entry/child checks for the
  // nested blossom are done at RELINK time (relinkOddCycleS5), AFTER B's relink sets its parent edge
  // -- they can't be checked here (sn's parent is still B-internal). structEligible = the depth-aware
  // "all members bare or z=0-recursively-structural" precheck.
  for (Sz i = 0; i < k; ++i) {
    Sz sn = cyc[i];
    if (sn >= mN) {
#ifdef PS5_NESTED_INTACT
      if (zLazy(sn) != 0) continue;     // intact sub-blossom: handled by on-arc rotateToBase at relink
      if (!relinkStructEligible(sn, 0)) {
#ifdef RELINKPROBE
        ++gRPnested;
#endif
        return false;
      }
#elif defined(PS5_NESTED_RELINK)
      if (zLazy(sn) != 0) {
#ifdef RELINKPROBE
        ++gRPnested; ++gRPnestSubSum;   // gRPnestSubSum REUSED = count of z!=0 (INTACT) sub-blossom bails
#endif
        return false;
      }
      if (!relinkStructEligible(sn, 0)) {
#ifdef RELINKPROBE
        ++gRPnested;                    // z=0 but struct-ineligible
#endif
        return false;
      }
#else
      // DEFAULT: nested sub-blossom on the cycle -> bail to relabel (the proven-safe path). The nested
      // in-place relink (PS5_NESTED_RELINK) is WIP -- correct target (99% of uniform-wt bails) but the
      // re-thread corrupts the tree on some cases (n=38 rep8 HANG) = the documented correctness trap.
#ifdef RELINKPROBE
      ++gRPnested;
#endif
      return false;
#endif
    }
  }
  Sz par = mParent[B], peU = mParEU[B];
  Sz entryIdx = NIL();
#ifdef PS5_LEAF_SPLICE
  // paper-semantics containment: LIFT peU to B's direct child (O(depth)), then index it -- replaces
  // the O(size*k) member scans; also stale-superset-proof.
  { Sz x = (par != NIL()) ? peU : NIL(); Sz g = 0;
    while (x != NIL() && mBParent[x] != B && g++ <= mNumNodes) x = mBParent[x];
    if (x != NIL() && mBParent[x] == B)
      for (Sz i = 0; i < k; ++i) if (cyc[i] == x) { entryIdx = i; break; }
  }
#else
  for (Sz i = 0; i < k; ++i) {
    const std::vector<Sz>& lv = mNodeLeaves[cyc[i]];
    for (Sz j = 0; j < lv.size(); ++j) if (lv[j] == peU) { entryIdx = i; break; }
    if (entryIdx != NIL()) break;
  }
#endif
  if (par != NIL() && entryIdx == NIL()) {
#ifdef RELINKPROBE
    ++gRPnoEntry;
#endif
    return false; }
  bool haveChild = false;
  for (Sz c = 0; c < mNumNodes; ++c) {
    if (mBParent[c] != NIL() || c == B || mParent[c] == NIL()) continue;
    if (outerNode(mParent[c]) != B) continue;
    Sz aU = mParEV[c]; bool found = false;
#ifdef PS5_LEAF_SPLICE
    { Sz x = aU; Sz g = 0;
      while (x != NIL() && mBParent[x] != B && g++ <= mNumNodes) x = mBParent[x];
      if (x != NIL() && mBParent[x] == B)
        for (Sz i = 0; i < k && !found; ++i) if (cyc[i] == x) found = true;
    }
#else
    for (Sz i = 0; i < k && !found; ++i) {
      const std::vector<Sz>& lv = mNodeLeaves[cyc[i]];
      for (Sz j = 0; j < lv.size(); ++j) if (lv[j] == aU) { found = true; break; }
    }
#endif
    if (!found) {
#ifdef RELINKPROBE
      ++gRPnoEntry;
#endif
      return false; }
    if (haveChild) {
#ifdef RELINKPROBE
      ++gRPmultiChild;
#endif
      return false; }
    haveChild = true;
  }
#ifdef RELINKPROBE
  ++gRPok;
#endif
  return true;
}

// Recursively relink + dissolve a z=0 blossom S that just became a ROOT (its parent B dissolved).
// S already has a label + parent edge set by B's relink. Handles S being ODD (inner) -- relink its
// cycle in place via relinkOddCycleS5, collect newly-EVEN sub-nodes, dissolve S, and recurse into any
// z=0 sub-blossom of S. Returns false (caller falls back to relabel) for cases not yet handled
// (S EVEN/FREE) -- safe. relinkedEven accumulates sub-nodes needing edge scheduling.
inline bool ScalingS5Engine::relinkRecDissolve(Sz B, std::vector<Sz>* relinkedEven) {
  static int s_depth = 0;
#ifdef PS5_DEEP_INTACT
  if (s_depth > 100000) return false;            // effectively-uncapped (spin-safety only)
#else
  if (s_depth > 16) return false;                // hard depth guard (data: shallow ~1) -- prevents spin
#endif
  if (B < mN) return true;                       // bare vertex: nothing to dissolve
  if (mBParent[B] != NIL()) return false;        // not a root yet -> can't dissolve in place; bail
  if (zLazy(B) != 0) return true;                // not dissolving
  if (mLabel[B] != LBL_ODD) return false;        // only ODD (inner) z=0 blossoms dissolve; else bail
  std::vector<Sz> kids = mBCycle[B];             // snapshot (relink/dissolve mutate structure)
  ++s_depth;
  bool ok = relinkOddCycleS5(B);
  --s_depth;
  if (!ok) return false;
  for (Sz i = 0; i < kids.size(); ++i)
    if (mLabel[kids[i]] == LBL_EVEN) relinkedEven->push_back(kids[i]);
  mLabel[B] = LBL_FREE;
  dissolve(B);
  // recurse: any cycle member that is a z=0 blossom is now a root with its label set -- dissolve it too.
  for (Sz i = 0; i < kids.size(); ++i) {
    Sz sn = kids[i];
    if (sn >= mN && mBParent[sn] == NIL() && zLazy(sn) <= 0) {
      ++s_depth; bool r = relinkRecDissolve(sn, relinkedEven); --s_depth;
      if (!r) return false;
    }
  }
  return true;
}

inline bool ScalingS5Engine::relinkOddCycleS5(Sz B) {
#ifdef PS5_SYNC_DIRTY
  // no-relabel dissolution re-threads MANY sub-node labels (on-arc ODD/EVEN, off-arc FREE) + re-parents
  // the child; enumerating exactly which leaves change is fragile. Force the safe FULL sync next time
  // (relink is rare per phase -- a few dissolveRelabels -- so this costs little of the dirty-set win).
  mSyncDirtyAll = true;
#endif
  const std::vector<Sz>& cyc = mBCycle[B];
  Sz k = cyc.size();
  if (k < 3 || (k % 2) == 0) return false;
#if defined(PS5_NESTED_INTACT)
  // INTACT-NESTED: both z=0 (re-threaded recursively by the caller) AND z!=0 (intact, rotated to base
  // on-arc below) sub-blossoms are now handled -> no bail on sub-blossoms here.
#elif defined(PS5_NESTED_RELINK)
  // NESTED (WIP): a z=0 sub-blossom on the cycle is re-threaded recursively by the caller after B
  // dissolves. Only a z!=0 (intact) sub-blossom is unsafe -> bail.
  for (Sz i = 0; i < k; ++i) if (cyc[i] >= mN && zLazy(cyc[i]) != 0) return false;
#else
  for (Sz i = 0; i < k; ++i) if (cyc[i] >= mN) return false;   // DEFAULT: any sub-blossom -> bail (proven-safe)
#endif
  Sz root = mRoot[B], par = mParent[B];
  Sz peU = mParEU[B], peV = mParEV[B];   // B's tree edge: peU in B, peV in parent
  // TREE-HEADLESS DECLINE (2026-07-12, the INTACT 594 heap bug): an ODD blossom with NO tree
  // parent (possible under persistent boundaries -- no relabel cleans ODD roots) has no entry
  // arc; the old par==NIL -> entryIdx=0 threading built a headless ODD top (mParent=NIL,
  // mParEU=NIL) that a later flipUpToRoot walks into: mMate[NIL] = the left-redzone heap
  // overflow (MATEGUARD n=34 rep=2: O=30 lbl=ODD par=-1). Decline -> proven relabel fallback.
  if (par == NIL()) return false;
  // entry = cycle index whose sub-node contains peU.
#ifdef S5DBG
  ++gS5RelinkCalls;
#endif
  // CLIFF FIX (2026-07-02): the old code scanned mNodeLeaves of EVERY cycle member (the whole
  // blossom) to locate the ONE cycle index containing a given vertex -- O(k * blossomLeaves) per
  // relink, near-quadratic in blossom size -> the n>=25600 uniform cliff (callgrind + cliffProbe
  // confirmed: relinkLeaf/relinkChild grow ~14x per doubling). Replace with cycIdxOf(v): the cycle
  // member is the child of B on v's mBParent chain, found in O(nesting-depth) (small) via a tiny
  // cycIdxOfNode[subnode]=i table (O(k)). No leaf scan; no all-mNumNodes child sweep.
  static std::vector<Sz> cycIdxOfNode;   // subnode -> its index in cyc (NIL if not a cycle member)
  if (cycIdxOfNode.size() < mNumNodes) cycIdxOfNode.assign(mNumNodes, NIL());
  for (Sz i = 0; i < k; ++i) cycIdxOfNode[cyc[i]] = i;
  // cycIdxOf(v): walk v's mBParent chain; the node whose mBParent == B is the containing cycle member.
  // (v is a leaf => v<mN or a sub-blossom; its chain passes through exactly one direct child of B.)
  #define CYC_IDX_OF(vtx, outIdx) do { \
      Sz _x = (vtx); Sz _res = NIL(); \
      if (cycIdxOfNode[_x] != NIL()) _res = cycIdxOfNode[_x]; \
      else { while (_x != NIL() && mBParent[_x] != B && mBParent[_x] != NIL()) _x = mBParent[_x]; \
             if (_x != NIL() && mBParent[_x] == B && cycIdxOfNode[_x] != NIL()) _res = cycIdxOfNode[_x]; } \
      (outIdx) = _res; } while(0)
  Sz entryIdx = NIL();
  CYC_IDX_OF(peU, entryIdx);
  if (par != NIL() && entryIdx == NIL()) { for (Sz i = 0; i < k; ++i) cycIdxOfNode[cyc[i]] = NIL(); return false; }
  if (par == NIL()) entryIdx = 0;
  // child sub-tree: a node c with parent resolving to B, attach vtx mParEV[c] inside B.
  Sz childAttachIdx = entryIdx; bool haveChild = false; Sz childNode = NIL();
#if !defined(PS5_NO_FASTCHILD)
  // CLIFF FIX #2 (2026-07-02): the old child search swept ALL mNumNodes (+ outerNode per node) EVERY
  // relink = O(n * numNodes) total, the OTHER half of the n>=25600 cliff (cliffProbe relinkChild ~14x
  // per doubling). B is an ODD tree blossom matched DOWNWARD to its unique child through B's base
  // (mBBase[B]): that base vertex's external mate lies in the child node. So the child + its attach
  // vertex are found in O(1). PS5_CHILDCHECK cross-checks against the old sweep (proves equivalence).
  {
    Sz baseVtx = mBBase[B];
    Sz childMate = (baseVtx < mN) ? mMate[baseVtx] : NIL();
    // child exists iff base is matched externally (to a vertex outside B) and that vertex is in-tree.
    if (childMate != NIL() && (childMate >= mN || outerNode(childMate) != B)) {
      Sz cN = outerNode(childMate);
      // confirm cN is genuinely B's tree child (its resolved parent is B) -- else no clean child.
      if (cN != B && cN != NIL() && mParent[cN] != NIL() && outerNode(mParent[cN]) == B) {
        Sz aU = mParEV[cN];
        Sz attach = NIL();
        CYC_IDX_OF(aU, attach);
        if (attach != NIL()) { childAttachIdx = attach; haveChild = true; childNode = cN; }
      }
    }
  }
#ifdef PS5_CHILDCHECK
  { Sz cChk = NIL(), aChk = NIL(); bool multi = false, badAttach = false;
    for (Sz c = 0; c < mNumNodes; ++c) {
      if (mBParent[c] != NIL() || c == B || mParent[c] == NIL()) continue;
      if (outerNode(mParent[c]) != B) continue;
      Sz aU = mParEV[c]; Sz at = NIL(); CYC_IDX_OF(aU, at);
      if (at == NIL()) { badAttach = true; break; }
      if (cChk != NIL()) { multi = true; break; }
      cChk = c; aChk = at; }
    bool sweepClean = !multi && !badAttach;
    bool fastMatches = sweepClean ? ((cChk == NIL()) ? (childNode == NIL())
                                     : (childNode == cChk && childAttachIdx == aChk))
                                  : true /*sweep would return false; fast may too or diverge -- log*/;
    if (!fastMatches || (!sweepClean)) std::fprintf(stderr,
      "[CHILDCHECK B=%lu fast(node=%ld attach=%ld) sweep(node=%ld attach=%ld multi=%d bad=%d)]\n",
      (unsigned long)B,(long)childNode,(long)childAttachIdx,(long)cChk,(long)aChk,(int)multi,(int)badAttach);
  }
#endif
#else
#ifdef S5DBG
  gS5RelinkChildScan += (long long)mNumNodes;
#endif
  for (Sz c = 0; c < mNumNodes; ++c) {
    if (mBParent[c] != NIL() || c == B || mParent[c] == NIL()) continue;
    if (outerNode(mParent[c]) != B) continue;
    Sz aU = mParEV[c];
    Sz attach = NIL();
    CYC_IDX_OF(aU, attach);
    if (attach == NIL()) { for (Sz i = 0; i < k; ++i) cycIdxOfNode[cyc[i]] = NIL(); return false; }
    if (haveChild) { for (Sz i = 0; i < k; ++i) cycIdxOfNode[cyc[i]] = NIL(); return false; }
    childAttachIdx = attach; haveChild = true; childNode = c;
  }
#endif
  #undef CYC_IDX_OF
  for (Sz i = 0; i < k; ++i) cycIdxOfNode[cyc[i]] = NIL();   // clear the tiny table for reuse
  Sz fwd = (childAttachIdx + k - entryIdx) % k;
  Sz dir, arcLen;
  if ((fwd % 2) == 0) { dir = 1; arcLen = fwd; }
  else { dir = k - 1; arcLen = k - fwd; }
  std::vector<char> onArc(k, 0);
  { Sz i = entryIdx; for (Sz s = 0; s <= arcLen; ++s) { onArc[i] = 1; i = (i + dir) % k; } }
  Sz prev = NIL();
  { Sz i = entryIdx;
    for (Sz s = 0; s <= arcLen; ++s) {
      Sz sn = cyc[i];
      mLabel[sn] = (s % 2 == 0) ? LBL_ODD : LBL_EVEN;
      markLabelDirty(sn);   // relink (dissolution re-thread): sn's role changes -> flip-dirty must catch it
      mRoot[sn] = root;
      if (s == 0) { MG_CHK("relink-entry", peU, peV); mParent[sn] = par; mParEU[sn] = peU; mParEV[sn] = peV;
#ifdef DIVDBG
        parw(sn,'R');
#endif
      }
      else {
        Sz eIdx = (dir == 1) ? ((i + k - 1) % k) : i;
        Sz eu = mBCycEU[B][eIdx], ev = mBCycEV[B][eIdx];
        // ORIENTATION FIX (2026-07-10, the n=176 livelock root): cycle edge j has eu[j] in cyc[j],
        // ev[j] in cyc[j+1]. FORWARD (dir==1) connects prev=cyc[i-1] -> sn=cyc[i] via edge i-1:
        // sn-side endpoint is ev, parent-side is eu. BACKWARD (dir==k-1) connects prev=cyc[i+1] ->
        // sn=cyc[i] via edge i: sn-side is eu, parent-side is ev. The old unconditional
        // {mParEU=ev, mParEV=eu} was wrong for BACKWARD arcs -> self-referential parent edges
        // (mParEV[sn]=sn's own vertex) -> the NEXT shrink over these nodes records misaligned
        // cycle-edge arrays -> rotateToBase silently no-ops -> stale mBBase -> consumed tree
        // reads treeLive -> grow re-parents into its own ancestor chain -> LIVELOCK.
        MG_CHK("relink-arc", eu, ev);
        if (dir == 1) { mParent[sn] = prev; mParEU[sn] = ev; mParEV[sn] = eu; }
        else          { mParent[sn] = prev; mParEU[sn] = eu; mParEV[sn] = ev; }
#ifdef DIVDBG
        parw(sn,'r');
#endif
      }
#ifdef PS5_NESTED_INTACT
      // INTACT on-arc sub-blossom: it stays contracted (z!=0) but the arc enters it at an INTERNAL
      // vertex (mParEU[sn]) -- so its base must be ROTATED to that vertex for the threading + the
      // internal matching to stay consistent (exactly what augment's rotateToBase does when an
      // augmenting path passes through an intact blossom). Without this the n=38-rep8 corruption/hang.
      if (sn >= mN && zLazy(sn) != 0) {
        Sz enterVtx = (s == 0) ? peU : mParEU[sn];   // the arc-entry vertex inside sn
#ifdef INTACT_TRACE
        // compute the EXIT vertex: the sn-side endpoint of the cycle edge to the next arc node (s+1).
        // For s<arcLen: exit is on the outgoing cycle edge. For s==arcLen (base): exits to childNode.
        { Sz exitV = NIL();
          if (s < arcLen) { Sz eIdxOut = (dir == 1) ? i : ((i + k - 1) % k);
            Sz oeu = mBCycEU[B][eIdxOut], oev = mBCycEV[B][eIdxOut];
            // pick the endpoint that is a leaf of sn
            exitV = leafOfNode(oeu, sn) ? oeu : oev; }
          std::fprintf(stderr,"[INTACT-ARC sn=%lu s=%lu/%lu label=%d entry=%lu exit=%ld oldBase=%lu | entry==exit?%d]\n",
            (unsigned long)sn,(unsigned long)s,(unsigned long)arcLen,mLabel[sn],(unsigned long)enterVtx,(long)exitV,(unsigned long)mBBase[sn],(int)(exitV==(long)enterVtx)); }
#endif
#ifdef PS5_INTACT_XFER
        // PARTIAL FIX (M==E root, docs/F_Plan): transfer the old base's EXTERNAL match to the new base.
        // ELIMINATES the corruption + M==E (M==E=no), but yields a VALID-but-SUBOPTIMAL matching (16663 vs
        // 16779 on n=200) -- the transfer is structurally right but not the exact optimal on-arc flip. The
        // exact fix needs the full alternating-flip derivation (s=0 vs s>0, ODD base matched DOWNWARD to
        // child vs EVEN matched UP to parent). Gated OFF; scan-restart stays default (correct+OPTIMAL).
        Sz oldBase = mBBase[sn];
        Sz X = (oldBase < mN) ? mMate[oldBase] : NIL();
        bool xExternal = (X != NIL() && !leafOfNode(X, sn));
        rotateToBase(sn, enterVtx);
        if (xExternal) { mMate[X] = enterVtx; mMate[enterVtx] = X; }
#elif defined(PS5_INTACT_NOROTATE)
        // EXPERIMENT: thread the tree but do NOT rotate the intact blossom's matching. Only valid if
        // enterVtx is already sn's base (no rotation needed). If not, this leaves the tree/matching
        // inconsistent -> test whether the search still converges correctly (relabel may fix it).
        (void)enterVtx;
#elif defined(PS5_INTACT_POLE)
        // POLE FIX (docs/F_Plan): an on-arc intact blossom is TRAVERSED. Its base must sit at the MATCHED
        // pole so its internal alternation is consistent + its external matched edge is preserved:
        //   EVEN on-arc: matched to PARENT (entry side) -> base = entry vertex.
        //   ODD  on-arc: matched DOWNWARD to child (exit side) -> base = EXIT vertex.
        // rotateToBase to that matched pole; its external partner (the matched-edge other end) is then
        // consistent with the arc. The UNMATCHED pole (the other end) needs no external match (it is the
        // exposed/entered side). This replaces the wrong rotate-to-entry-always.
        {
          // Compute sn's two arc poles (vertices inside sn) + their external arc partners:
          //   ENTRY pole = mParEU[sn], partner = mParEV[sn] (the parent-side vertex).
          //   EXIT  pole = the sn-side endpoint of the OUTGOING arc edge, partner = the other endpoint.
          //     s<arcLen: outgoing = the cycle edge to the next arc node.
          //     s==arcLen (arc base): outgoing = the child-attach edge to childNode (partner=mParEU[child]).
          Sz entryPole = enterVtx, entryPartner = (s == 0) ? peV : mParEV[sn];
          Sz exitPole = NIL(), exitPartner = NIL();
          if (s < arcLen) {
            Sz eIdxOut = (dir == 1) ? i : ((i + k - 1) % k);
            Sz oeu = mBCycEU[B][eIdxOut], oev = mBCycEV[B][eIdxOut];
            if (leafOfNode(oeu, sn)) { exitPole = oeu; exitPartner = oev; } else { exitPole = oev; exitPartner = oeu; }
          } else if (haveChild && childNode != NIL()) {
            exitPole = mParEV[childNode];   // vertex inside sn where the child attaches
            exitPartner = mParEU[childNode];
          }
          // MATCHED pole by parity: EVEN on-arc is matched to its PARENT (entry edge matched); ODD is
          // matched DOWNWARD to its child/exit (exit edge matched). Rotate the blossom base to the MATCHED
          // pole (it becomes exposed there), then set its external match to the arc partner.
          bool even = (mLabel[sn] == LBL_EVEN);
          Sz matchedPole  = even ? entryPole : (exitPole != NIL() ? exitPole : entryPole);
          Sz matchedPartn = even ? entryPartner : (exitPartner != NIL() ? exitPartner : entryPartner);
          rotateToBase(sn, matchedPole);
          if (matchedPartn != NIL()) { mMate[matchedPole] = matchedPartn; mMate[matchedPartn] = matchedPole; }
          // NOTE: the ODD-arcLen0-no-child + recursively-dissolved nested sub-blossom cases (n=60 W100 s2,
          // n=150 W5 s3) are the last ~0.8% -- pole logic here does not cover them; caught by dual-cert.
        }
#else
        rotateToBase(sn, enterVtx);
        // KNOWN BUG (docs/F_Plan): rotateToBase re-matches sn's OLD base internally, dropping the base's
        // EXTERNAL match -> asymmetry -> M==E on distinct weights. PS5_INTACT_XFER above transfers the
        // external match (kills M==E) but is suboptimal; the exact-optimal on-arc flip is derived+pending.
        // Scan-restart on M==E keeps results correct+OPTIMAL (just slower) until then.
#endif
      }
#endif
      prev = sn; i = (i + dir) % k;
    }
  }
  for (Sz i = 0; i < k; ++i) if (!onArc[i]) {
    Sz sn = cyc[i];
    mLabel[sn] = LBL_FREE; mParent[sn] = NIL(); mParEU[sn] = NIL(); mParEV[sn] = NIL(); mRoot[sn] = NIL();
  }
  // Re-point the child sub-tree's parent from the (now-dissolved) B to the base sub-node it
  // actually hangs off (cyc[childAttachIdx]). Without this, childNode->mParent still points into
  // dissolved B, and any tree walk (treeParent/outerNode) from C loops forever (the n=18 hang).
  // The child keeps its own mParEU/mParEV (the real attach edge); only the parent NODE changes.
  if (haveChild && childNode != NIL()) { mParent[childNode] = cyc[childAttachIdx];
#ifdef DIVDBG
    parw(childNode,'C');
#endif
  }
#ifdef RELINKDBG
  // SELF-CHECK the alternation: each on-arc node (s>=1) is reached from its parent via an edge whose
  // matched-status must agree with its label (EVEN<-matched, ODD<-unmatched). The base (ODD) is
  // matched DOWNWARD (to child), not to its arc-parent. Print any violation.
  { Sz i2 = entryIdx;
    for (Sz s = 0; s <= arcLen; ++s) {
      Sz sn = cyc[i2];
      if (s >= 1) {
        Sz eu = mParEU[sn], ev = mParEV[sn];
        bool matched = (mMate[eu] == ev);
        bool wantMatched = (mLabel[sn] == LBL_EVEN);   // EVEN child reached via matched edge
        if (matched != wantMatched)
          std::fprintf(stderr,"[ALT-VIOLATION B=%lu s=%lu sn=%lu lbl=%d edge(%lu,%lu) matched=%d want=%d arcLen=%lu dir=%lu]\n",
            (unsigned long)B,(unsigned long)s,(unsigned long)sn,mLabel[sn],(unsigned long)eu,(unsigned long)ev,(int)matched,(int)wantMatched,(unsigned long)arcLen,(unsigned long)dir);
      }
      i2 = (i2 + dir) % k;
    }
  }
#endif
  return true;
#ifdef DIVDBG
  { Sz bi9=0; Sz mB9 = findCycEdgeMisalign(&bi9); static bool s_r9=false;
    if (mB9 != NIL() && !s_r9) { s_r9=true; std::fprintf(stderr,"[MISALIGN-END-RELINK B=%lu idx=%lu k=%lu]\n",
      (unsigned long)mB9,(unsigned long)bi9,(unsigned long)mBCycle[mB9].size()); std::fflush(stderr); } }
#endif
}

// dissolve inner root blossom B (z==0): free its sub-nodes; rebuild union-find
// for the freed leaves (each freed leaf's outer node is re-derived).
inline void ScalingS5Engine::dissolve(Sz B) {
#ifdef LF_OPTRACE
  std::fprintf(stderr,"[LF-OP dissolve B=%lu]%c",(unsigned long)B,10);
#endif
#ifdef TRAJ
  { std::fprintf(stderr,"[TRAJ DISS B=%lu base=%lu cyclelen=%lu lazy=%d | cycle-mates:",(unsigned long)B,(unsigned long)mBBase[B],(unsigned long)mBCycle[B].size(),(int)mLazyActive);
    for(Sz i=0;i<mBCycle[B].size();++i){Sz sn=mBCycle[B][i]; Sz bs=(sn<mN)?sn:mBBase[sn]; std::fprintf(stderr," %lu(m%ld,l%d)",(unsigned long)bs,(long)(bs<mN?mMate[bs]:-1),mLabel[sn]);}
    std::fprintf(stderr,"]\n"); }
#endif
#ifdef PS5_NODE_Y
  // FAITHFUL SETTLE-AT-DISSOLVE (§5.2, docs/J_Scope, mirror of the shrink settle): B's leaves were
  // accruing under B's role (mNRole[B]). B is about to vanish as their top-level; freeze that accrual
  // into mYAcc NOW, then hand each sub-blossom/bare-vertex a fresh stamp = mTNow so the post-dissolve
  // relink+lazySyncRoles starts their new role's interval from zero. O(leaves of B) = piggybacks on the
  // leaf re-split loop below. Without it, the child's stamp would double-count B's interval.
  if (mLazyActive || mLazyVerify) {   // scan verify mode maintains lazy timestamps too (see shrink)
    int br = mNRole[B]; LL bs = mNStamp[B]; LL span = mTNow - bs;
    int d = (br == LBL_ODD) ? +1 : (br == LBL_EVEN ? -1 : 0);
#ifdef PS5_LEAF_SPLICE
    if (d && span) MBX_LF_FOR(v, B) { if (v >= mN) continue;
#else
    std::vector<Sz>& bl = mNodeLeaves[B];
    if (d && span) for (Sz j = 0; j < bl.size(); ++j) { Sz v = bl[j]; if (v >= mN) continue;
#endif
#ifdef VTXWATCH
      if((long long)v==(VTXWATCH)) std::fprintf(stderr,"[DISS-SETTLE v=%lu B=%lu br=%d span=%lld yAccBefore=%lld @mTNow=%lld]\n",(unsigned long)v,(unsigned long)B,br,(long long)span,(long long)mYAcc[v],(long long)mTNow);
#endif
      mYAcc[v] += d * span;
    }
    // each surviving sub-entity restarts its own interval at now (role set by the coming relink/sync).
    for (Sz i = 0; i < mBCycle[B].size(); ++i) { Sz sub = mBCycle[B][i];
      if (sub < mN) { mVRole[sub] = LBL_FREE; mVStamp[sub] = mTNow; }
      else          { mNRole[sub] = LBL_FREE; mNStamp[sub] = mTNow; }
    }
    // B itself is gone: settle its OWN z-interval then mark FREE@now, so a later setNodeRole(B,FREE)
    // (lazySyncRoles/relink) is a no-op and does NOT re-settle B's departed members OR its z for the same
    // interval (the dissolve double-count: DISS-SETTLE + NODE-SETTLE both firing span=79 -- n=16 rep=11).
    // Settle z here (setNodeRole would ALSO re-settle members, which we've already done) so zLazy stays exact.
    if (br == LBL_EVEN)      mZAcc[B] += 2 * span;
    else if (br == LBL_ODD)  mZAcc[B] -= 2 * span;
    mNRole[B] = LBL_FREE; mNStamp[B] = mTNow;
  }
#endif
#ifdef PS5_LEAF_LIST
  { // sever inter-child links AFTER the NODE_Y settle (which walks B's list) and right BEFORE
    // the child release. Placing this at ENTRY nulled the list the settle needed -> zero settles
    // -> the 580 v=6 under-accrual + the 594 n=180 stuck search.
    for (Sz i = 0; i < mBCycle[B].size(); ++i) { Sz c = mBCycle[B][i];
      if (c >= mLfHead.size() || mLfHead[c] == NIL()) continue;
      mLfPrev[mLfHead[c]] = NIL(); mLfNext[mLfTail[c]] = NIL();
      lfCheckNode(c, "dissolve-child");
    }
    if (B < mLfHead.size()) { mLfHead[B] = NIL(); mLfTail[B] = NIL(); mLfCnt[B] = 0; }
    lfSweepAll("dissolve-global");
  }
#endif
  for (Sz i = 0; i < mBCycle[B].size(); ++i) mBParent[mBCycle[B][i]] = NIL();
#ifdef PS5_OUTERCACHE
  outerCacheBump();   // mBParent changed -> invalidate outerNode cache (moved BEFORE the CQ split so
                      // outerNode(v) below reflects the just-completed dissolve)
#endif
#ifdef PS5_CQ_EVENTS
  if (!mCqOrder.empty()) {
#ifdef PS5_CQ_REINIT_DISSOLVE
    cqInitMembership();   // diagnostic reference: rebuild (correct, but O(n)/op -- not shippable)
#else
    // CURRENT-ORDER SPLIT (robust to any merge-order permutation): B's merged sublist now spans several
    // sub-blossoms; each cycle member's mBParent is cleared, so a vertex's outerNode is its NEW top-level.
    // Walk B's sublist in CQ order and split after any leaf whose NEXT leaf has a DIFFERENT outerNode --
    // peels each maximal same-node run into its own sublist. Reads the CQ's ACTUAL order, so it does not
    // depend on merge order or mBCycle. O(|B| log n) at the (rare) dissolve.
    if (!mBCycle[B].empty()) {
      // anyLeaf must be a vertex genuinely in B's (soon-to-split) sublist -- use the AUTHORITATIVE base of
      // the first cycle member (mNodeLeaves is stale). B's sublist = the merged run of all cycle members.
      Sz sub0 = mBCycle[B][0];
      Sz anyLeaf = (sub0 < mN) ? sub0 : mBBase[sub0];
      if (anyLeaf != NIL() && anyLeaf < mN) {
        std::vector<Sz> verts; mSFM.sublistVerts(anyLeaf, verts);
        for (Sz i = 0; i + 1 < verts.size(); ++i) {
          if (((ScalingS5Engine*)this)->outerNode(verts[i]) != ((ScalingS5Engine*)this)->outerNode(verts[i + 1]))
            mSFM.split(verts[i]);   // boundary between two sub-blossoms
        }
      }
    }
#endif
#ifdef PS5_CQ_CHECK_EVERY_OP
    cqCheckMembership("dissolve");
#endif
  }
#endif
  // union-find: re-split. Re-point each leaf to the uf-rep of its (now smaller)
  // outer node. Simplest correct: reset uf for all leaves of B to their own
  // outer node's base, recursively.
  for (Sz i = 0; i < mBCycle[B].size(); ++i) {
    Sz child = mBCycle[B][i];
    Sz rep = mBBase[child];
#ifdef PS5_NESTED_RELINK
    // NESTED FIX: if `child` is itself a LIVE blossom (stays after B dissolves), do NOT flatten its
    // leaves' uf-rep to mBBase[child) -- that breaks the uf<->mBParent consistency (a leaf deep in the
    // child would have uf-rep=child.base but outerNode(leaf)=child, and a later shrink mis-threads ->
    // the 6<->46 parent cycle). Only re-point the child's OWN base to itself; leaves keep resolving
    // within the child via mBParent. Bare-vertex children flatten as before.
    // GUARD (2026-07-01): this is ONLY correct for the NO-RELABEL (lazy) path. In the SCAN engine
    // (mLazyActive==false, relabel() rebuilds the forest every iteration), skipping the flatten leaves the
    // uf inconsistent with the fresh relabel -> the scan engine converges SUBOPTIMALLY (n=60 W100 s2:
    // 2405 vs 2412). The scan path MUST flatten as before. So skip only when lazy.
#if defined(PS5_NO_DISS_UF_FLATTEN)
    // OLD GUARD (default OFF now): SKIP re-pointing a live blossom child's leaves. This LEFT them at their
    // STALE rep from the dissolved parent B -> ufFind(leaf) != child's base = the uf<->mBParent DESYNC that
    // made the search treat internal blossom edges as live cross-blossom + take a wrong augment (n=100 W50
    // s3, ~all the desync-class residuals). The "6<->46 cycle" this guard was added for is no longer an
    // issue post-fixes (broad sweep 299/300, 348/348 gate, zero regressions with the flatten below).
    if (mLazyActive && child >= mN) { mUf[rep] = rep; continue; }
#endif
    // FAITHFULNESS FIX (2026-07-02, DEFAULT-ON): after B dissolves, a live blossom `child` becomes the outer
    // node for its leaves -> every leaf's ufFind MUST resolve to mBBase[child]. Flatten them (same as bare
    // children) so uf agrees with the mBParent tree. Fixes the nested-dissolution uf desync -> the residual
    // suboptimality. Disable with -DPS5_NO_DISS_UF_FLATTEN.
#endif
    std::vector<Sz>& cl = mNodeLeaves[child];
    for (Sz j = 0; j < cl.size(); ++j) mUf[cl[j]] = rep;
    mUf[rep] = rep;
  }
  // NOTE: id reclamation (free-list) is DISABLED -- reusing a dissolved id while
  // other nodes' mParent/mBParent still point to it causes dangling refs (segfault)
  // in the no-relabel path. Safe reclamation needs structural-pointer invalidation
  // (the §5.1 #1b work). With relabel-on-shrink the engine is correct without it.
  // mFreeBlossoms.push_back(B);   // re-enable only with full pointer invalidation
#ifdef PS5_ANCESTOR_HINT
  if (B < mAhDead.size()) mAhDead[B] = 1;   // hints through B are now invalid (id stays dead)
#endif
#ifdef DIVDBG
  { Sz bi9=0; Sz mB9 = findCycEdgeMisalign(&bi9); static bool s_r9=false;
    if (mB9 != NIL() && !s_r9) { s_r9=true; std::fprintf(stderr,"[MISALIGN-END-DISSOLVE B=%lu idx=%lu k=%lu]\n",
      (unsigned long)mB9,(unsigned long)bi9,(unsigned long)mBCycle[mB9].size()); std::fflush(stderr); } }
#endif
}

// augment across cross-tree eligible edge (uOrig in E, vOrig in W)
inline void ScalingS5Engine::augment(Sz E, Sz W, Sz uOrig, Sz vOrig) {
#ifdef OUTERTAG
  { extern long long gOTaug, gOTaugCalls; gOTaug -= gS5OuterCalls; ++gOTaugCalls; }
#endif
#ifdef MATEGUARD
  if (mLabel[E] != LBL_EVEN || mLabel[W] != LBL_EVEN)
    std::fprintf(stderr,"[MG-AUG-ENTRY E=%lu lbl=%d W=%lu lbl=%d u=%lu w=%lu | outer(u)=%lu outer(w)=%lu]\n",
      (unsigned long)E,mLabel[E],(unsigned long)W,mLabel[W],(unsigned long)uOrig,(unsigned long)vOrig,
      (unsigned long)outerNode(uOrig),(unsigned long)outerNode(vOrig));
#endif
#ifdef TRAJ
  { LL wt=-1; for(Sz k=0;k<mAdjV[uOrig].size();++k) if(mAdjV[uOrig][k]==vOrig){wt=mAdjW[uOrig][k];break;}
    std::fprintf(stderr,"[TRAJ AUG (%lu,%lu) w=%lld]\n",(unsigned long)uOrig,(unsigned long)vOrig,(long long)wt); }
#endif
  if (E >= mN) rotateToBase(E, uOrig);
  if (W >= mN) rotateToBase(W, vOrig);
  flipUpToRoot(E);
  // DEAD-TREE-COERCION GUARD (2026-07-12, the INTACT 594 heap bug): the event loop treats a
  // dead-tree ODD target as FREE (lblW coercion) -- but that is a LOCAL relabel; mLabel[W] stays
  // ODD and W may carry a STALE tree chain (the INTACT relink re-threads dead remnants, 'r'
  // writer tags). flipUpToRoot(W) on such a W walks the stale chain into a foreign live root and
  // pairs mMate[mParEU[root]=NIL] = the left-redzone heap overflow (MATEGUARD n=34 rep=2:
  // E=7 lbl=EVEN, W=68 lbl=ODD coerced). Only a genuinely-EVEN W has a live tree to flip; a
  // truly-FREE W's flip was always a no-op (mParent==NIL).
  if (mLabel[W] == LBL_EVEN) flipUpToRoot(W);
  MG_CHK("augment", uOrig, vOrig);
  mMate[uOrig] = vOrig; mMate[vOrig] = uOrig;
#ifdef OUTERTAG
  { extern long long gOTaug; gOTaug += gS5OuterCalls; }
#endif
#ifdef DIVDBG
  { Sz bi9=0; Sz mB9 = findCycEdgeMisalign(&bi9); static bool s_r9=false;
    if (mB9 != NIL() && !s_r9) { s_r9=true; std::fprintf(stderr,"[MISALIGN-END-AUGMENT B=%lu idx=%lu k=%lu]\n",
      (unsigned long)mB9,(unsigned long)bi9,(unsigned long)mBCycle[mB9].size()); std::fflush(stderr); } }
#endif
}

// reset labels; root every FREE outer node EVEN. (Repair malformed blossoms first:
// a free vertex must be a blossom's base.)
inline void ScalingS5Engine::relabel() {
#ifdef OUTERTAG
  { extern long long gOTrelabel, gOTrelabelCalls; gOTrelabel -= gS5OuterCalls; ++gOTrelabelCalls; }
#endif
  for (Sz v = 0; v < mNReal; ++v)              // mNReal: skip dummy slots [mNReal,mN) (not rooted)
    if (mMate[v] == NIL()) {
      Sz N = outerNode(v);
      while (N >= mN && mBBase[N] != v) { dissolve(N); N = outerNode(v); }
    }
  for (Sz i = 0; i < mNumNodes; ++i) {
    if (mBParent[i] != NIL()) continue;
    mLabel[i] = LBL_FREE; mParent[i] = NIL();
    mParEU[i] = NIL(); mParEV[i] = NIL(); mRoot[i] = NIL();
  }
#ifdef PS5_FLIP_DIRTY
  // relabel() bulk-rewrites EVERY root label WITHOUT per-node marking. The event engine always
  // follows with lazyBegin (full sync); the SCAN engine does not -> under flip-dirty the roles
  // went stale after every relabel (the 580 LAZYVERIFY y-mismatch, v=1 t=555). Arm the designed
  // fallback: the next lazySyncRoles does one FULL sync, then re-arms incremental.
  mFlipAll = true;
#endif
  for (Sz v = 0; v < mNReal; ++v)              // root only reals; dummy slots have no FVR edges
    if (mMate[v] == NIL()) {
      if (mRootRestrict && !(v < mRootAllow.size() && mRootAllow[v])) continue;  // step-6 PQSearch(F)
      Sz N = outerNode(v);
      if (mLabel[N] == LBL_FREE) { mLabel[N] = LBL_EVEN; mRoot[N] = N; }
    }
#ifdef OUTERTAG
  { extern long long gOTrelabel; gOTrelabel += gS5OuterCalls; }
#endif
}

// One EdmondsSearch (Stage 3.1: scan-driven; union-find answers same-blossom).
// Returns true if it augmented >=1 path; false if no perfect matching reachable.
inline bool ScalingS5Engine::searchOnce() {
  bool augmentedAny = false;
  for (;;) {
    relabel();
    bool anyRoot = false;
    for (Sz i = 0; i < mNumNodes; ++i)
      if (mBParent[i] == NIL() && mLabel[i] == LBL_EVEN) { anyRoot = true; break; }
    if (!anyRoot) return augmentedAny;

    // ---- GROW to fixpoint ----
    bool didShrink = false;
    std::vector<Sz> wq;
    for (Sz i = 0; i < mNumNodes; ++i)
      if (mBParent[i] == NIL() && mLabel[i] == LBL_EVEN) wq.push_back(i);
    Sz wh = 0;
    while (wh < wq.size() && !didShrink) {
      Sz E = wq[wh++];
      if (mBParent[E] != NIL() || mLabel[E] != LBL_EVEN) continue;
#ifdef PS5_LEAF_SPLICE
      Sz _lfE = (E < mLfTail.size()) ? mLfTail[E] : NIL();
      for (Sz u = (E < mLfHead.size() ? mLfHead[E] : NIL()), _lfNx = NIL(); u != NIL() && !didShrink; u = _lfNx) {
        _lfNx = (u == _lfE) ? NIL() : mLfNext[u];
#else
      std::vector<Sz>& Ev = mNodeLeaves[E];
      for (Sz ui = 0; ui < Ev.size() && !didShrink; ++ui) {
        Sz u = Ev[ui];
#endif
        for (Sz k = 0; k < mAdjV[u].size() && !didShrink; ++k) {
          Sz w = mAdjV[u][k];
          if (ufFind(u) == ufFind(w)) continue;      // same outer blossom (UF, O(a))
          Sz W = outerNode(w);
          if (W == E) continue;
          if (slackZ(u, w) != 0) continue;           // Criterion 1: tight only
          if (mLabel[W] == LBL_EVEN) {
            if (lcaOuter(E, W) != NIL()) { shrink(E, W, u, w); didShrink = true; break; }
            else { augment(E, W, u, w); augmentedAny = true; return true; }
          } else if (mLabel[W] == LBL_ODD) {
            continue;
          } else {
            Sz bw = mBBase[W]; Sz m = mMate[bw];
            if (m == NIL()) { augment(E, W, u, w); augmentedAny = true; return true; }
            Sz M = outerNode(m);
            MG_CHK("grow-scan", m, bw);
            mLabel[W] = LBL_ODD; mParent[W] = E; mParEU[W] = w; mParEV[W] = u; mRoot[W] = mRoot[E];
            mLabel[M] = LBL_EVEN; mParent[M] = W; mParEU[M] = m; mParEV[M] = bw; mRoot[M] = mRoot[E];
#ifdef DIVDBG
            parw(W,'g'); parw(M,'g');
#endif
            markLabelDirty(W); markLabelDirty(M);   // scan-engine grow: W,M flip -> flip-dirty must catch
            wq.push_back(M);
          }
        }
      }
    }
    if (didShrink) continue;

    // ---- DUAL ADJUSTMENT (delta-jump) ----
    LL best = -1;
    for (Sz u = 0; u < mN; ++u) {
      Sz U = outerNode(u);
      if (mLabel[U] != LBL_EVEN) continue;
      for (Sz k = 0; k < mAdjV[u].size(); ++k) {
        Sz w = mAdjV[u][k];
        if (ufFind(u) == ufFind(w)) continue;
        Sz W = outerNode(w);
        if (W == U) continue;
        LL s = slackZ(u, w);
        int lw = mLabel[W];
        if (lw == LBL_EVEN) { LL d = s / 2; if (d >= 0 && (best < 0 || d < best)) best = d; }
        else if (lw == LBL_FREE) { LL d = s; if (d >= 0 && (best < 0 || d < best)) best = d; }
      }
    }
    for (Sz b = mN; b < mNumNodes; ++b)
      if (mBParent[b] == NIL() && mLabel[b] == LBL_ODD) {
        LL d = mZ[b] / 2; if (d >= 0 && (best < 0 || d < best)) best = d;
      }
    if (best < 0) return augmentedAny;        // no event possible -> stop

    // pre-dissolve zero-z ODD blossoms (z-clamp insight from the reference engine)
    bool dissolvedAny = false;
    for (Sz b = mN; b < mNumNodes; ++b)
      if (mBParent[b] == NIL() && mLabel[b] == LBL_ODD && mZ[b] == 0) { dissolve(b); dissolvedAny = true; }
    if (dissolvedAny) continue;

    LL delta = (best == 0) ? 1 : best;
    if (mLazyVerify) { lazySyncRoles(); mTNow += delta; }   // mirror the dual step lazily
    for (Sz v = 0; v < mN; ++v) {
      Sz V = outerNode(v);
      if (mLabel[V] == LBL_EVEN) mY[v] -= delta;
      else if (mLabel[V] == LBL_ODD) mY[v] += delta;
    }
    for (Sz b = mN; b < mNumNodes; ++b)
      if (mBParent[b] == NIL()) {
        if (mLabel[b] == LBL_EVEN) mZ[b] += 2 * delta;
        else if (mLabel[b] == LBL_ODD) mZ[b] -= 2 * delta;
      }
    if (mLazyVerify) lazyVerifyStep();   // lazy y/z must match materialized post-step
    for (Sz b = mN; b < mNumNodes; ++b)
      if (mBParent[b] == NIL() && mLabel[b] == LBL_ODD && mZ[b] == 0) dissolve(b);
  }
}

// ---------------------------------------------------------------------------
// Stage 3.2: event-driven search via BucketPQ. The forest persists across ALL
// events (grow/shrink/augment/dissolution) -- NO relabel between events. The
// union-find same-outer-blossom test (ufFind) is what makes deferred events
// safely validatable on pop (the piece the earlier bolted-on PQ lacked). Each
// event for edge (u,w) is re-validated on pop: u still in a live EVEN node, w not
// in u's outer blossom, edge tight NOW. Stale events are discarded.
// ---------------------------------------------------------------------------

// a tree is LIVE iff its root's base vertex is still unmatched (after an augment
// consumes a tree, its root's base is matched -> the whole tree's events die).
inline bool ScalingS5Engine::treeLive(Sz N) const {
  Sz r = mRoot[N];
  if (r == NIL()) return false;
  // mRoot is always a tree ROOT, which is always an OUTER node (roots are never inside a
  // blossom), so outerNode(r)==r -- skip the redundant O(depth) outerNode walk (hot path:
  // called per edge in scheduleNodeEdges).
  return mMate[mBBase[r]] == NIL();
}

// ---------------------------------------------------------------------------
// FULL-EVEN RESEED FUNNEL (L_Scope_InvariantKeyFrame Stage 0): every "reschedule ALL
// EVEN roots' edges" loop in the engine goes through here so the reseed volume can be
// attributed per call SITE. Behavior-identical to the loops it replaced. Site ids:
//   0 searchEvent entry (once/search, legit)   1 phase-boundary relabel
//   2 NRD dual-fix (in-place dissolution)      3 relink-decline fallback
//   4..14 shrink/augment handler teardowns (one id per source site; see call sites)
// ---------------------------------------------------------------------------
#ifdef S5DBG
static long long gS5ReseedCalls[16] = {0};   // loop executions per site
static long long gS5ReseedSched[16] = {0};   // scheduleNodeEdges calls issued per site
static long long gS5ReseedScans[16] = {0};   // edge scans issued per site (gS5EdgeScans delta)
static double    gS5ReseedSecs[16]  = {0};   // wall seconds per site
// N_Scope Stage 0 (invariant keys): how much of each reseed re-derives SURVIVING state?
// A rescanned node "repeats" if it was already top-level EVEN at the previous reseed (label-only
// proxy -- leaf-set churn between reseeds is not detected; see docs/N_Scope_InvariantKeys.md).
static long long gS5ReseedRepeat[16] = {0};  // rescanned nodes already EVEN last reseed
static long long gS5ReseedRepScan[16] = {0}; // edge scans spent on those repeat nodes
static long long gS5ReseedEpoch = 0;         // reseed sequence number (all sites share one stream)
#endif
inline void ScalingS5Engine::reseedAllEven(LL tNow, int site) {
#ifdef PS5_VTX_BEST
  ++mVBestEpoch;   // every full reseed follows a PQ teardown -> per-vertex bests reset O(1)
#endif
#ifdef S5DBG
  struct timespec _t0, _t1; clock_gettime(CLOCK_MONOTONIC, &_t0);
  long long sched = 0, scans0 = gS5EdgeScans;
  // N_Scope Stage 0: repeat detection (was this node top-level EVEN at the previous reseed too?)
  long long rsEpoch = ++gS5ReseedEpoch, repeat = 0, repScan = 0;
  if (mRsEp.size() < mZ.size()) mRsEp.resize(mZ.size(), -9);
#endif
  for (Sz i = 0; i < mNumNodes; ++i)
    if (mBParent[i] == NIL() && mLabel[i] == LBL_EVEN) {
#ifdef S5DBG
      ++sched;
      bool rep = (mRsEp[i] == rsEpoch - 1); mRsEp[i] = rsEpoch;
      long long preScan = gS5EdgeScans;
#endif
      scheduleNodeEdges(i, tNow);
#ifdef S5DBG
      if (rep) { ++repeat; repScan += gS5EdgeScans - preScan; }
#endif
    }
#ifdef S5DBG
  clock_gettime(CLOCK_MONOTONIC, &_t1);
  if (site >= 0 && site < 16) {
    ++gS5ReseedCalls[site]; gS5ReseedSched[site] += sched;
    gS5ReseedScans[site] += gS5EdgeScans - scans0;
    gS5ReseedSecs[site] += (_t1.tv_sec-_t0.tv_sec)+(_t1.tv_nsec-_t0.tv_nsec)/1e9;
    gS5ReseedRepeat[site] += repeat; gS5ReseedRepScan[site] += repScan;
  }
#else
  (void)site;
#endif
}

// schedule events for every edge out of EVEN node E, at the absolute time the
// edge becomes tight under continued dual adjustment.
inline void ScalingS5Engine::scheduleNodeEdges(Sz E, LL tNow, LL hBase) {
#ifdef OUTERTAG
  { extern long long gOTsched; gOTsched -= gS5OuterCalls; }
#endif
#ifdef PS5_EVT_HORIZON
  // HORIZON DEFER (2026-07-05, from the QSTAT measurement: 92% of inserts are never popped --
  // phases typically end within <=3 ticks while inserts spread out to 100K+ ticks, and far
  // inserts pay the multiset overflow path). Insert only events firing within H ticks of the
  // window base (tNow, or the marker time on a marker rescan); for the rest keep ONE per-node
  // marker (kind 3) at their MIN fire time. Popping a marker rescans the node with the window
  // re-based at the marker time, so the min-fire edge is always inserted exactly -> the clock
  // can never overrun a deferred edge (marker time <= every deferred fire time).
  // SIZE-SCALED HORIZON (2026-07-07, the n=204800 giant-blossom cliff): a marker rescan re-walks
  // the WHOLE node to insert one 8-tick window; a giant blossom's boundary spans a huge tick range,
  // so it needs range/H marker pops, EACH an O(leaves x deg) walk -> 20567 pops x 5437 leaves =
  // ~1B edge-scans = the cliff (viaMarker=20567 vs viaReseed=14, S5DBG). FIX: widen the window in
  // proportion to node size so ONE walk covers ~leaves-worth of events and the pop count drops by
  // ~leaves. The extra inserts are cheap relative to the walk we already pay. Small nodes (the
  // 92%-waste case the horizon targets) keep the tight H window.
  LL qsH = (LL)(PS5_EVT_HORIZON);
#ifndef PS5_HORIZON_SIZEMUL
#define PS5_HORIZON_SIZEMUL 16   // H*|B|/16: measured optimum (n=204800 giant blossom 48->19s; flat 4..64)
#endif
  // NOTE (2026-07-08): a regime-adaptive "no-defer above NODEFER leaves" was tried and REVERTED --
  // qsH=2^60 dumps all far events into the BucketPQ multiset overflow (inserts 44M->120M), strictly
  // WORSE than either the size-scaled window OR full horizon-OFF. The horizon defer BACKFIRES on
  // GIANT blossoms (marker-pop rescans: pathological n=150000 = 104s vs horizon-OFF 30.5s = 3.3x)
  // but plain PS5_NO_EVT_HORIZON is the clean win there; no leaf-count threshold beats it (the
  // "giant" boundary is instance-dependent). See scaling-real-bench-landscape memory.
  { LL L = (LL)mNodeLeaves[E].size();
    if (L > 1) qsH *= L / (LL)(PS5_HORIZON_SIZEMUL) + 1; }  // window ~ H*|B|/k: balance pops vs inserts
  LL qsBase = (hBase < 0) ? tNow : hBase;
#ifdef PS5_MARKER_ADAPT
  // this node has been marker-thrashing -> stop deferring it: schedule its events into a LARGE window
  // that stays INSIDE the bucket array (mWindowCap ~ 1<<20), so no marker re-pop AND no multiset flood.
  if (mNoDeferNode.size() > E && mNoDeferNode[E]) { LL room = ((LL)1 << 20) - 8 - qsBase; if (room > qsH) qsH = room; }
#endif
  LL qsDefMin = -1;
#endif
#ifdef S5DBG
  ++gS5SchedCalls;
  // L_Scope Stage 0: split calls + edge scans by trigger -- marker rescan (hBase>=0) vs plain.
  // RAII so the scan delta is booked on EVERY exit path (cache-skip returns included).
  struct _MkAcct { bool on; _MkAcct(bool o):on(o){ if(on){ extern long long gS5SchedMarkerCalls, gS5SchedMarkerScanBase; ++gS5SchedMarkerCalls; gS5SchedMarkerScanBase -= gS5EdgeScans; } }
    ~_MkAcct(){ if(on){ extern long long gS5SchedMarkerScanBase; gS5SchedMarkerScanBase += gS5EdgeScans; } } } _mk(hBase >= 0);
  // node-level redundancy: was THIS node already scheduled at the same epoch+tNow?
  if (mDbgSchedEpoch.size() > E && mDbgSchedEpoch[E] == mDbgEpoch && mDbgSchedT[E] == tNow)
    ++gS5NodeReSched;
  if (mDbgSchedEpoch.size() > E) { mDbgSchedEpoch[E] = mDbgEpoch; mDbgSchedT[E] = tNow; }
#endif
#ifdef PS5_LEAF_SHADOW
  lfCheckNode(E, "sched");
#endif
  std::vector<Sz>& Ev = mNodeLeaves[E];
#ifdef PS5_BIGNODE_CACHE
  // BIG-NODE EVENT CACHE (2026-07-07, the n>=300k giant-blossom cliff -- 86.5% of big-node scans
  // are provably redundant, MEMOIZABLE probe). Cache a big node's scheduled events as (u,w,kind,
  // offset=fireT-tNow). On a rescan where the node's STRUCTURE EPOCH and LEAF COUNT are unchanged,
  // replay events at tNow+offset instead of the O(leaves x deg) walk. SAFETY: replayed events are
  // re-validated at POP exactly like any scheduled event (label/treeLive/ufFind/slack-at-target),
  // so a STALE cached event (neighbor role changed) is harmlessly discarded -- the engine already
  // tolerates stale queue entries (lazy deletion). The invalidation key is E-side; if a W-side
  // change ever caused a MISSED event the overrun-guard + gates catch it (verified: 594 + same-
  // graph GMG). Big nodes only (small nodes never repeat enough to pay the cache bookkeeping).
  const Sz kBigCacheThresh = 256;
  bool bcBig = (Ev.size() > kBigCacheThresh);
  if (bcBig) {
    if (mBcEpoch.size() <= E) { mBcEpoch.assign(mNumNodes+8, -1); mBcLeaves.assign(mNumNodes+8, (Sz)-1); mBcEvents.resize(mNumNodes+8); }
    if (mBcEpoch[E] == mDbgEpoch && mBcLeaves[E] == Ev.size()) {
#ifdef PS5_INCR_EVENTS
      // STAGE A: if the queue was NOT reset/cleared since this node's events were scheduled
      // (generation unchanged), the events are STILL LIVE in the queue -- do NOT re-emit them.
      // A w that became newly eligible since then schedules the SAME boundary edge from ITS OWN
      // side when scheduleNodeEdges(w) runs (adjacency is symmetric) -- HYPOTHESIS under test by
      // the same-graph-GMG gate @204800. Stale events die at pop. Only the marker is (re)set so
      // the far-window events get a fresh rescan trigger.
      if (mBcGen.size() > E && mBcGen[E] == mPQ.generation()) {
#ifdef S5DBG
        ++gS5BigCacheHit; ++gS5IncrSkip;
#endif
        // TRUE O(1) SKIP: same generation => this node's events (AND its marker) are STILL LIVE in
        // the queue from the prior emission. Re-emitting is pure waste (the 0.04%-change measurement).
        // A newly-eligible neighbor w schedules the boundary edge from ITS side (adjacency symmetry,
        // GMG-exact @204800). So we do NOTHING -- no walk, no reschedule, no marker. O(1).
        return;
      }
#endif
      // CACHE HIT: replay. Offsets are (fire - buildTNow); re-base to current tNow.
      std::vector<BcEvent>& ce = mBcEvents[E];
      for (Sz i = 0; i < ce.size(); ++i) {
        LL fire = tNow + ce[i].off;
#ifdef S5DBG
        // PREDICTION TEST (2026-07-07): how many replayed events have a w whose OUTER LABEL CHANGED
        // since caching? If small + not scaling with blossom size -> an incremental "touch only
        // changed-w events" fix is O(changes) not O(size) = the real fix. If large -> events churn,
        // no handle scheme helps (Piece 3 / CQ needed). Reasoned-from-structure, then measured.
        { Sz w = (Sz)ce[i].w; int nowLbl = mLabel[outerNode(w)];
          if (nowLbl != LBL_FREE && !treeLive(outerNode(w))) nowLbl = LBL_FREE;
          ++gS5ReplayTotal; if (nowLbl != ce[i].wlbl) ++gS5ReplayWchanged; }
#endif
#ifdef PS5_EVT_HORIZON
        if (fire > qsBase + qsH) { if (qsDefMin < 0 || fire < qsDefMin) qsDefMin = fire; continue; }
#endif
        mPQ.schedule(fire, ce[i].u, ce[i].w, ce[i].kind);
      }
#ifdef PS5_EVT_HORIZON
      if (qsDefMin >= 0) mPQ.schedule(qsDefMin, (long long)E, -1, 3);
#endif
#ifdef S5DBG
      ++gS5BigCacheHit;
#endif
      return;   // skipped the whole walk
    }
    // CACHE MISS: rebuild. Collect events into mBcEvents[E] as we walk (below), then stamp key.
    mBcEvents[E].clear();
  }
#endif
#ifdef S5DBG
  if ((LL)Ev.size() > gS5BigNodeThresh) { ++gS5BigNodeCalls; gS5BigNodeLeafSum += (long long)Ev.size();
    if (hBase >= 0) ++gS5BigViaMarker; else ++gS5BigViaReseed;
    // MEMOIZABILITY PROBE: is this big-node scan REDUNDANT? A scan is memoizable if since the last
    // scan of THIS node (a) the structure epoch is unchanged (no shrink/dissolve/relabel touched it)
    // AND (b) its leaf count is unchanged. If so, the events are recomputable from the prior scan +
    // clock delta -> a cache would skip the O(leaves x deg) walk. Counts redundant vs necessary.
    if (mBigLastEpoch.size() <= E) { mBigLastEpoch.assign(mNumNodes+8, -1); mBigLastLeaves.assign(mNumNodes+8, (Sz)-1); }
    bool structSame = (mBigLastEpoch[E] == mDbgEpoch && mBigLastLeaves[E] == Ev.size());
    if (structSame) ++gS5BigRedundant; else ++gS5BigNecessary;
    mBigLastEpoch[E] = mDbgEpoch; mBigLastLeaves[E] = Ev.size();
  }
#endif
#ifdef PS5_LEAF_SPLICE
  // Stage 2: the HOT boundary walk reads the spliced list (vector still maintained as the net).
  Sz _lfE = (E < mLfTail.size()) ? mLfTail[E] : NIL();
  for (Sz u = (E < mLfHead.size() ? mLfHead[E] : NIL()), _lfNx = NIL(); u != NIL(); u = _lfNx) {
    _lfNx = (u == _lfE) ? NIL() : mLfNext[u];
    Sz ufu = ufFind(u);
#else
  for (Sz ui = 0; ui < Ev.size(); ++ui) {
    Sz u = Ev[ui];
    Sz ufu = ufFind(u);                              // loop-invariant: hoist out of edge loop
#endif
#ifdef S5DBG
    gS5EdgeScans += (long long)mAdjV[u].size();
#endif
    for (Sz k = 0; k < mAdjV[u].size(); ++k) {
      Sz w = mAdjV[u][k];
      if (mIgnoreDummies && w >= mNReal) continue;   // D: skip dummy-slot targets (parked separately)
      if (ufu == ufFind(w)) {                         // same outer blossom (u-side hoisted)
#ifdef S5DBG
        ++gS5InternalScans;
#endif
        continue;
      }
      Sz W = outerNode(w);
#ifdef PS5_OUTER_SAMEBLOSS
      // (Tested: authoritative outerNode same-blossom check -- prevents scheduling a phantom FORMATION on an
      // INTERNAL edge when ufFind desyncs from mBParent. SAFE (348/348, no LEMON regression) but INERT on the
      // 2 residuals -> the desync's harm is NOT the phantom scheduling (main loop discards it via E==W at pop);
      // the residual is in the augment PATH choice. Kept gated OFF -- adds a walk per edge for no gain.)
      if (outerNode(u) == W) continue;
#endif
      int lw = mLabel[W];
      if (lw != LBL_FREE && !treeLive(W)) lw = LBL_FREE;
      if (lw != LBL_EVEN && lw != LBL_FREE) continue; // EVEN-ODD: no event -> skip slack calc
      LL s = slackZW(u, w, 2 * mAdjW[u][k]);          // weight in hand -> no w2max rescan
#ifdef PS5_EVTRACE
      std::fprintf(stderr,"[TR] scan node=%lu edge(%lu,%lu) lw=%d s=%lld tNow=%lld\n",
        (unsigned long)E,(unsigned long)u,(unsigned long)w,lw,(long long)s,(long long)tNow);
#endif
      // Criterion target: an edge is eligible when slack reaches `tgt` (CRIT1: 0; CRIT2/3: -2, since
      // grow/formation edges are UNMATCHED -- matched edges are followed, not scheduled). Event fires
      // at tNow + (s-tgt)/step. Matches the loop engine's unmTarget = (crit==CRIT1)?0:-2.
      LL tgt = (mCrit == CRIT1) ? 0 : -2;
      // NEGATIVE-SLACK FIX (D3, the carried-dual / Property-2 regime): if slack is already AT or PAST
      // the eligibility target (d<=0), the edge is eligible NOW -- schedule at tNow (not dropped).
      // The old engine (full-domination Property 1, y=maxw start) never saw d<0; but the SCALING
      // algorithm carries Property-2 duals (near-domination, slack can be <0 after rescale), so a
      // negative-slack edge must be an IMMEDIATE event, else it is invisible -> suboptimal matching.
      if (lw == LBL_EVEN) {
        LL d = (s - tgt) / 2;                          // EVEN-EVEN: 2/step
        if (mAllowOvertight && d < 0) d = 0;           // over-tight -> eligible now (Property-2 regime)
#ifdef S5SCHEDDBG
        { static LL maxd=0; if (d>maxd){maxd=d; std::fprintf(stderr,"[SCHEDDBG] new max d=%lld (EVEN-EVEN) slack=%lld tNow=%lld\n",(long long)d,(long long)s,(long long)tNow);} }
#endif
        if (mAllowOvertight || d >= 0) {
#ifdef PS5_BIGNODE_CACHE
          if (bcBig) { BcEvent be; be.u=(long long)u; be.w=(long long)w; be.kind=1; be.off=(d<0?0:d); be.wlbl=lw; mBcEvents[E].push_back(be); }
#endif
#ifdef PS5_EVT_HORIZON
          if (tNow + (d < 0 ? 0 : d) > qsBase + qsH) {
            LL fire = tNow + (d < 0 ? 0 : d);
            if (qsDefMin < 0 || fire < qsDefMin) qsDefMin = fire;
            continue;
          }
#endif
#ifdef PS5_VTX_BEST
          if (vbSuppress(w, tNow + d, 1)) continue;
#endif
          mPQ.schedule(tNow + d, (long long)u, (long long)w, 1);
#ifdef PS5_QSTAT
          ++mQIns; ++mQInsHist[qstatBucket(d < 0 ? 0 : d)];
#endif
#ifdef PS5_SFM_INCR
          if (mSFMon) mSFM.decreasekey(u, (Val)(tNow + d));   // incremental key(u) = min fire-time
#endif
#ifdef S5DBG
          ++gS5SchedInserts;
#endif
        }
      } else if (lw == LBL_FREE) {
        LL d = s - tgt;                                // EVEN-free: 1/step
        if (mAllowOvertight && d < 0) d = 0;           // over-tight -> eligible now
        if (mAllowOvertight || d >= 0) {
#ifdef PS5_BIGNODE_CACHE
          if (bcBig) { BcEvent be; be.u=(long long)u; be.w=(long long)w; be.kind=0; be.off=(d<0?0:d); be.wlbl=lw; mBcEvents[E].push_back(be); }
#endif
#ifdef PS5_EVT_HORIZON
          if (tNow + (d < 0 ? 0 : d) > qsBase + qsH) {
            LL fire = tNow + (d < 0 ? 0 : d);
            if (qsDefMin < 0 || fire < qsDefMin) qsDefMin = fire;
            continue;
          }
#endif
#ifdef PS5_VTX_BEST
          if (vbSuppress(w, tNow + d, 0)) continue;
#endif
          mPQ.schedule(tNow + d, (long long)u, (long long)w, 0);
#ifdef PS5_QSTAT
          ++mQIns; ++mQInsHist[qstatBucket(d < 0 ? 0 : d)];
#endif
#ifdef PS5_SFM_INCR
          if (mSFMon) mSFM.decreasekey(u, (Val)(tNow + d));   // incremental key(u) = min fire-time
#endif
#ifdef S5DBG
          ++gS5SchedInserts;
#endif
        }
      }
      // EVEN-ODD: yz unchanged by adjustment; no event.
    }
  }
#ifdef PS5_EVT_HORIZON
  if (qsDefMin >= 0) {
    mPQ.schedule(qsDefMin, (long long)E, -1, 3);   // kind 3 = "rescan node E" at its deferred min
#ifdef PS5_EVTRACE
    std::fprintf(stderr,"[TR] sched MARKER node=%lu at t=%lld (base=%lld tNow=%lld)\n",
      (unsigned long)E,(long long)qsDefMin,(long long)qsBase,(long long)tNow);
#endif
#ifdef PS5_QSTAT
    ++mQIns; ++mQInsHist[qstatBucket(qsDefMin - tNow)];
#endif
  }
#endif
#ifdef PS5_BIGNODE_CACHE
  if (bcBig) { mBcEpoch[E] = mDbgEpoch; mBcLeaves[E] = Ev.size();   // stamp: cache now valid for this (epoch,leaves)
#ifdef PS5_INCR_EVENTS
    if (mBcGen.size() <= E) mBcGen.assign(mNumNodes+8, -1);
    mBcGen[E] = mPQ.generation();   // events just scheduled into the queue at THIS generation
#endif
  }
#endif
#ifdef OUTERTAG
  { extern long long gOTsched; gOTsched += gS5OuterCalls; }
#endif
}

#ifdef PS5_SFM
// Ground-truth next-event fire-time: min over all currently schedulable EVEN-incident edges of the SAME
// fire-time scheduleNodeEdges would compute. This is exactly what the reseed builds + the PQ then pops.
inline ScalingS5Engine::LL ScalingS5Engine::sfmBruteNextFire(LL tNow) const {
  LL best = -1;
  LL tgt = (mCrit == CRIT1) ? 0 : -2;
  for (Sz E = 0; E < mNumNodes; ++E) {
    if (mBParent[E] != NIL() || mLabel[E] != LBL_EVEN) continue;
    const std::vector<Sz>& Ev = mNodeLeaves[E];
    for (Sz ui = 0; ui < Ev.size(); ++ui) {
      Sz u = Ev[ui]; Sz ufu = ((ScalingS5Engine*)this)->ufFind(u);
      for (Sz k = 0; k < mAdjV[u].size(); ++k) {
        Sz w = mAdjV[u][k];
        if (mIgnoreDummies && w >= mNReal) continue;
        if (ufu == ((ScalingS5Engine*)this)->ufFind(w)) continue;
        Sz W = outerNode(w);
        int lw = mLabel[W];
        if (lw != LBL_FREE && !treeLive(W)) lw = LBL_FREE;
        if (lw != LBL_EVEN && lw != LBL_FREE) continue;
        LL s = ((ScalingS5Engine*)this)->slackZW(u, w, 2 * mAdjW[u][k]);
        LL d;
        if (lw == LBL_EVEN) { d = (s - tgt) / 2; if (mAllowOvertight && d < 0) d = 0; }
        else                { d = (s - tgt);     if (mAllowOvertight && d < 0) d = 0; }
        if (!mAllowOvertight && d < 0) continue;
        LL fire = tNow + d;
        if (best < 0 || fire < best) best = fire;
      }
    }
  }
  return best;
}

// Same next-event fire-time computed via split-findmin: build a SplitFindmin whose key(u) = the min
// fire-time of an EVEN-incident edge at u (over inner/any vertices), then findmin over the whole vertex
// set. (Rebuilt fresh here -- this validates the QUERY semantics; incremental maintenance is the next
// step. The whole-set findmin = min over all vertices = the global next event.)
inline ScalingS5Engine::LL ScalingS5Engine::sfmFindminNextFire(LL tNow) {
  std::vector<Sz> order; order.reserve(mN);
  for (Sz v = 0; v < mN; ++v) order.push_back(v);
  std::vector<char> end(mN, 0); if (mN) end[mN - 1] = 1;     // one sublist: the whole set
  mSFM.init(order, end);
  LL tgt = (mCrit == CRIT1) ? 0 : -2;
  // key(u) = min fire-time of an EVEN-incident edge whose EVEN endpoint's leaf is u.
  for (Sz E = 0; E < mNumNodes; ++E) {
    if (mBParent[E] != NIL() || mLabel[E] != LBL_EVEN) continue;
    std::vector<Sz>& Ev = mNodeLeaves[E];
    for (Sz ui = 0; ui < Ev.size(); ++ui) {
      Sz u = Ev[ui]; Sz ufu = ufFind(u);
      for (Sz k = 0; k < mAdjV[u].size(); ++k) {
        Sz w = mAdjV[u][k];
        if (mIgnoreDummies && w >= mNReal) continue;
        if (ufu == ufFind(w)) continue;
        Sz W = outerNode(w);
        int lw = mLabel[W];
        if (lw != LBL_FREE && !treeLive(W)) lw = LBL_FREE;
        if (lw != LBL_EVEN && lw != LBL_FREE) continue;
        LL s = slackZW(u, w, 2 * mAdjW[u][k]);
        LL d;
        if (lw == LBL_EVEN) { d = (s - tgt) / 2; if (mAllowOvertight && d < 0) d = 0; }
        else                { d = (s - tgt);     if (mAllowOvertight && d < 0) d = 0; }
        if (!mAllowOvertight && d < 0) continue;
        mSFM.decreasekey(u, (Val)(tNow + d));     // key(u) tracks the min fire-time at u
      }
    }
  }
  Val k; Sz v; mSFM.findmin((Sz)0, &k, &v);       // findmin over the whole vertex set
  if (k >= SplitFindmin::INF) return -1;
  return (LL)k;
}

// init the SFM over all vertices as ONE sublist, all keys INF; mSFMon enables the scheduleNodeEdges feed.
inline void ScalingS5Engine::sfmInit() {
  std::vector<Sz> order; order.reserve(mN);
  for (Sz v = 0; v < mN; ++v) order.push_back(v);
  std::vector<char> end(mN, 0); if (mN) end[mN - 1] = 1;
  mSFM.init(order, end);
  mSFMon = true;
}

// next-event fire-time from the LIVE incrementally-maintained structure (no rebuild): findmin over the
// whole vertex set. (Under PS5_SFM_INCR the keys were fed by scheduleNodeEdges' decreasekey.)
inline ScalingS5Engine::LL ScalingS5Engine::sfmFindminLive() {
  if (mN == 0) return -1;
  Val k; Sz v; mSFM.findmin((Sz)0, &k, &v);
  if (k >= SplitFindmin::INF) return -1;
  return (LL)k;
}

inline void ScalingS5Engine::sfmCheckAfterReseed(LL tNow) {
  LL a = sfmBruteNextFire(tNow);
#ifdef PS5_SFM_INCR
  LL b = sfmFindminLive();        // incremental: query the live structure (the REAL test)
#else
  LL b = sfmFindminNextFire(tNow);// rebuild-each-time (query-semantics-only validation)
#endif
  ++mSfmChecks;
  if (a != b) { ++mSfmMismatch;
    if (mSfmMismatch <= 10) std::fprintf(stderr, "[SFM MISMATCH] tNow=%lld brute=%lld findmin=%lld\n",
                                         (long long)tNow, (long long)a, (long long)b);
  }
}
#endif

#ifdef PS5_CQ_EVENTS
// K_Scope Stage 1a: build ONE ConcatenableQueue queue per EVEN top-level blossom, holding the fire-time
// of every currently-schedulable boundary edge on that blossom's EVEN side (EXACTLY the edges + fire-time
// scheduleNodeEdges/sfmBruteNextFire compute). Then assert the GLOBAL min over all blossom queues equals
// sfmBruteNextFire. This validates: (1) the per-blossom decomposition covers all events, (2) CQ
// minPriority/minElem is correct, (3) fire-time keying matches the brute -- the prerequisites for Stage 2
// (findmin-driven selection) and Stage 1b (incremental merge/split maintenance). REBUILT each call = a
// pure shadow (O(m log n) extra, NOT in the shipping path -- flag-gated, referee only).
inline void ScalingS5Engine::cqCheckAfterReseed(LL tNow) {
  LL tgt = (mCrit == CRIT1) ? 0 : -2;
  mCQ.clear();
  LL globalMin = -1;   // min fire-time across all blossom queues
  for (Sz E = 0; E < mNumNodes; ++E) {
    if (mBParent[E] != NIL() || mLabel[E] != LBL_EVEN) continue;
    ConcatenableQueue::Handle q = cNullItm;
    const std::vector<Sz>& Ev = mNodeLeaves[E];
    for (Sz ui = 0; ui < Ev.size(); ++ui) {
      Sz u = Ev[ui]; Sz ufu = ((ScalingS5Engine*)this)->ufFind(u);
      for (Sz k = 0; k < mAdjV[u].size(); ++k) {
        Sz w = mAdjV[u][k];
        if (mIgnoreDummies && w >= mNReal) continue;
        if (ufu == ((ScalingS5Engine*)this)->ufFind(w)) continue;
        Sz W = outerNode(w);
        int lw = mLabel[W];
        if (lw != LBL_FREE && !treeLive(W)) lw = LBL_FREE;
        if (lw != LBL_EVEN && lw != LBL_FREE) continue;
        LL s = ((ScalingS5Engine*)this)->slackZW(u, w, 2 * mAdjW[u][k]);
        LL d = (lw == LBL_EVEN) ? (s - tgt) / 2 : (s - tgt);
        if (mAllowOvertight && d < 0) d = 0;
        if (!mAllowOvertight && d < 0) continue;
        LL fire = tNow + d;
        // one leaf per boundary edge; payload = w (diagnostic), priority = fire-time.
        mCQ.insertRight(q, (Val)fire, (Size)w);
      }
    }
    if (q != cNullItm) {
      LL bm = (LL)mCQ.minPriority(q);
      if (globalMin < 0 || bm < globalMin) globalMin = bm;
    }
  }
  LL brute = sfmBruteNextFire(tNow);   // ground truth (PS5_SFM must also be defined)
  ++mCqChecks;
  if (brute != globalMin) { ++mCqMismatch;
    if (mCqMismatch <= 10) std::fprintf(stderr, "[CQ MISMATCH] tNow=%lld brute=%lld cqGlobalMin=%lld\n",
                                        (long long)tNow, (long long)brute, (long long)globalMin);
  }
}

// STAGE 1b: init mSFM sublists = the vertices of each top-level node, laid out so each node's leaves are
// contiguous (init's invariant). Order = concatenation of mNodeLeaves[top-level nodes]; listEnd marks each
// node's last leaf. After this, sameList(u,v) iff u,v share a top-level node.
inline void ScalingS5Engine::cqInitMembership() {
  mCqOrder.clear(); mCqOrder.reserve(mN);
  std::vector<char> end; end.reserve(mN);
  // Group vertices by their AUTHORITATIVE current top-level node (outerNode), NOT mNodeLeaves -- which is
  // append-only and STALE (a vertex can appear in mNodeLeaves[b] after it has left b; that stale membership
  // is the referee mismatch). Bucket by outerNode, then emit each bucket contiguously.
  std::vector<std::vector<Sz> > byNode(mNumNodes);
  for (Sz v = 0; v < mN; ++v) { Sz on = ((ScalingS5Engine*)this)->outerNode(v); if (on < mNumNodes) byNode[on].push_back(v); }
  for (Sz b = 0; b < mNumNodes; ++b) {
    if (byNode[b].empty()) continue;
    for (Sz j = 0; j < byNode[b].size(); ++j) { mCqOrder.push_back(byNode[b][j]); end.push_back(0); }
    end.back() = 1;                                      // last leaf of this node ends its sublist
  }
  if (!mCqOrder.empty()) mSFM.init(mCqOrder, end);
}

// STAGE 1b referee: for a sample of adjacent-in-order vertex pairs AND cross-node pairs, assert
// mSFM.sameList(u,v) == (outerNode(u)==outerNode(v)). Fires on the FIRST inconsistency (a merge/split bug).
inline void ScalingS5Engine::cqCheckMembership(const char* where) {
  if (mCqOrder.empty()) return;
  Sz bad = 0;
  for (Sz i = 0; i + 1 < mCqOrder.size(); ++i) {
    Sz u = mCqOrder[i], v = mCqOrder[i + 1];
    bool wantSame = (((ScalingS5Engine*)this)->outerNode(u) == ((ScalingS5Engine*)this)->outerNode(v));
    bool gotSame = mSFM.sameList(u, v);
    if (wantSame != gotSame) {
      ++bad; ++mCqMismatch;
      if (mCqMismatch <= 10) std::fprintf(stderr, "[CQ-MEMBER MISMATCH @%s u=%lu v=%lu want=%d got=%d outU=%lu outV=%lu]\n",
        where, (unsigned long)u, (unsigned long)v, (int)wantSame, (int)gotSame,
        (unsigned long)((ScalingS5Engine*)this)->outerNode(u), (unsigned long)((ScalingS5Engine*)this)->outerNode(v));
      if (bad > 5) return;
    }
  }
  ++mCqChecks;
}
#endif

inline void ScalingS5Engine::applyDualTo(LL newTime, LL* tNow) {
  LL delta = newTime - *tNow;
  // time only ADVANCES: never move tNow (or mTNow) backward. A newTime <= tNow is a no-op (the event is
  // already due now). Moving backward drives the negative-z divergence (seed=0). Belt-and-suspenders with
  // the disT clamp at the call site.
  if (delta <= 0) return;
  if (mLazyActive) {
    // LAZY: a dual adjustment is just advancing mTNow -- y(v)/z(B) are computed on
    // demand from the (verified) timestamps. Roles were synced at the last structural
    // change; advancing time accrues the open-interval term automatically. O(1).
#ifdef SHADOW_DUAL
    // ground-truth materialized y (by CURRENT OUTER LABEL) maintained in parallel; compare to yLazy after
    // the advance -> first divergence = the drift bug. mShadowY seeded at lazyBegin from y0.
    if (!mShadowY.empty()) {
      for (Sz v = 0; v < mN; ++v) { Sz V = outerNode(v);
        if (mLabel[V] == LBL_EVEN) mShadowY[v] -= delta; else if (mLabel[V] == LBL_ODD) mShadowY[v] += delta; }
      mTNow = newTime;
      for (Sz v = 0; v < mN; ++v) if (yLazy(v) != mShadowY[v]) {
        extern long long gShadowV; if (gShadowV < 0 || (Sz)gShadowV == v) {
          std::fprintf(stderr,"[SHADOW-DIVERGE v=%lu yLazy=%lld shadow=%lld role=%d stamp=%lld yAcc=%lld y0=%lld @mTNow=%lld]\n",
            (unsigned long)v,(long long)yLazy(v),(long long)mShadowY[v],mVRole[v],(long long)mVStamp[v],(long long)mYAcc[v],(long long)mY0[v],(long long)mTNow);
          gShadowV = v; break; } }
      *tNow = newTime; return;
    }
#endif
#ifdef PS5_ZSHADOW
    if (mZShadowOk) {
      // bump ROOT blossom z + every vertex y by the exact rule, then compare to zLazy/yLazy.
      for (Sz b = mN; b < mNumNodes; ++b) if (mBParent[b] == NIL()) {
        zShadowEnsure(b);
        if (mLabel[b] == LBL_EVEN) mZShadow[b] += 2 * delta;
        else if (mLabel[b] == LBL_ODD) mZShadow[b] -= 2 * delta;
      }
      for (Sz v = 0; v < mN; ++v) { Sz V = outerNode(v);
        if (mLabel[V] == LBL_EVEN) mYShadow[v] -= delta; else if (mLabel[V] == LBL_ODD) mYShadow[v] += delta; }
      mTNow = newTime;
      for (Sz v = 0; v < mN; ++v) if (yLazy(v) != mYShadow[v]) {
        if (gZShadowB < 0) { std::fprintf(stderr,"[YSHADOW-DIVERGE v=%lu yLazy=%lld shadow=%lld role=%d outerLabel=%d stamp=%lld yAcc=%lld y0=%lld @mTNow=%lld]\n",
          (unsigned long)v,(long long)yLazy(v),(long long)mYShadow[v],mVRole[v],mLabel[outerNode(v)],(long long)mVStamp[v],(long long)mYAcc[v],(long long)mY0[v],(long long)mTNow); gZShadowB = mNumNodes+v; break; } }
      for (Sz b = mN; b < mNumNodes; ++b) if (mBParent[b] == NIL()) {
        zShadowEnsure(b);
        if (zLazy(b) != mZShadow[b]) {
          if (gZShadowB < 0 || (Sz)gZShadowB == b) {
            std::fprintf(stderr,"[ZSHADOW-DIVERGE B=%lu zLazy=%lld shadow=%lld role=%d label=%d stamp=%lld zAcc=%lld z0=%lld @mTNow=%lld]\n",
              (unsigned long)b,(long long)zLazy(b),(long long)mZShadow[b],mNRole[b],mLabel[b],(long long)mNStamp[b],(long long)mZAcc[b],(long long)mZ0[b],(long long)mTNow);
            gZShadowB = b; }
        }
      }
      *tNow = newTime; return;
    }
#endif
    mTNow = newTime; *tNow = newTime; return;
  }
  for (Sz v = 0; v < mN; ++v) {
    Sz V = outerNode(v);
    if (mLabel[V] == LBL_EVEN) mY[v] -= delta;
    else if (mLabel[V] == LBL_ODD) mY[v] += delta;
  }
  for (Sz b = mN; b < mNumNodes; ++b)
    if (mBParent[b] == NIL()) {
      if (mLabel[b] == LBL_EVEN) mZ[b] += 2 * delta;
      else if (mLabel[b] == LBL_ODD) mZ[b] -= 2 * delta;
    }
  *tNow = newTime;
}

inline bool ScalingS5Engine::searchEvent() {
#if defined(S5DBG) || defined(PS5_BIGNODE_CACHE)
  ++mDbgEpoch;   // new search = new epoch (cache-invalidating; also S5DBG redundancy detector)
#endif
#if defined(PS5_VTX_BEST) && defined(PS5_VB_DEFER_REPAIR)
  mVbRepList.clear(); ++mVbRepTok;   // pending repairs from a returned search are moot (entry reseed)
#endif
#ifdef PS5_VTX_BEST
  vbRegimeDecide();   // per-search O(m) 5-slot scan; scale weights change between searches
#endif
#ifdef S5DBG
  ++gS5Searches;
  if (mDbgSchedEpoch.size() < mN) { mDbgSchedEpoch.assign(mN, -1); mDbgSchedT.assign(mN, -1); }
#endif
#ifdef PS5_OUTERCACHE
  outerCacheEnsure();
#endif
  relabel();
  bool anyRoot = false;
  for (Sz i = 0; i < mNumNodes; ++i)
    if (mBParent[i] == NIL() && mLabel[i] == LBL_EVEN) { anyRoot = true; break; }
  if (!anyRoot) return false;

  LL tNow = 0;
  // LAZY mode (S5-C): drive the event search on the verified §5.2 lazy duals. A dual
  // adjustment is mTNow++; y(v)/z(B) read on demand. lazyBegin snapshots the base +
  // seeds roles from the just-relabeled labels. Before each time advance we lazySyncRoles
  // so structural label changes (grow/shrink/dissolve since the last advance) are settled.
  mLazyActive = mLazyEvent;
#ifdef PS5_PERSIST_CLOCK
  mTNow = 0;   // fresh search: the clock is monotone WITHIN a search (per scale = faithful)
#endif
  if (mLazyActive) lazyBegin();
#ifdef PS5_ZSHADOW
  // seed the persistent shadow z from the correct base at search start; never re-seed after this.
  { gZShadowB = -1; mZShadow.assign(mZ.size(), 0); mYShadow.assign(mN, 0);
    for (Sz b = mN; b < mNumNodes; ++b) mZShadow[b] = zLazy(b);
    for (Sz v = 0; v < mN; ++v) mYShadow[v] = yLazy(v); mZShadowOk = true; }
#endif
  if (mNoRelabelAug) mUsed.assign(mN, 0);   // per-phase disjointness marks
  Sz phaseAugs = 0;                          // augmentations batched by the CURRENT phase
#ifndef PS5_NO_PHASE_AUGCAP
  Sz phaseAugSqrt = 1; while (phaseAugSqrt * phaseAugSqrt < mN) ++phaseAugSqrt;  // integer sqrt, once
  Sz mPhaseAugCap = 64 * phaseAugSqrt + 64;  // see the mega-phase cap note below. 64*sqrt(n) chosen
                                             // by sweep (MBX_AUGCAP {4,16,64}rt(n) vs uncapped): any
                                             // finite cap kills the mega-phase (DNF -> ~3s); away from
                                             // the pathology effects are single-seed reorder noise, and
                                             // 64rt(n) is behavior-neutral within noise on BOTH healthy
                                             // regimes (tie 22.4 vs 24.7 uncapped; dist 84.0 vs 85.8).
                                             // The cap is a SAFETY BOUND, not a tuning knob.
  { const char* e = std::getenv("MBX_AUGCAP"); if (e) mPhaseAugCap = (Sz)atoll(e); }  // sweep hook
#else
  Sz mPhaseAugCap = (Sz)-1; (void)mPhaseAugCap;
#endif
  mPQ.clear(); mPQ.init(1);
#if defined(PS5_PERSIST_FOREST) && defined(PS5_PERSIST_CLOCK)
  mReoffer.clear();
#endif
#ifdef PS5_SFM_INCR
  sfmInit();   // structure live; scheduleNodeEdges below feeds key(u) incrementally
#endif
#ifdef PS5_CQ_EVENTS
  cqInitMembership();               // STAGE 1b: live per-node membership sublists (split at dissolve, ...)
#ifdef PS5_CQ_CHECK_EVERY_OP
  cqCheckMembership("search-init"); // must be consistent at start (no shrinks yet)
#endif
#endif
  reseedAllEven(tNow, 0);
#ifdef PS5_DISHEAP
  disHeapRebuild(tNow);   // seed the dissolution heap from the initial ODD roots
#endif
  // NOTE: dissolution is NOT scheduled as a PQ event (that races against the
  // materialized z when a blossom changes inner<->outer role -- the n=8 hang). It
  // is handled the PROVEN scan way: before each PQ event we compute the next
  // dissolution time (min over ODD root blossoms of tNow + z/2, O(#blossoms)); if a
  // dissolution is due no later than the next PQ event, we dissolve + relabel +
  // rebuild the PQ (dissolutions are O(n) total). This keeps grow/formation on the
  // fast event path while dissolution stays exact.

  LL guard = 0, guardCap = 8 * (LL)mN * (LL)mN + 8 * (LL)mN * (mMaxW + 1) + 1000;
  bool augmentedAny = false;   // no-relabel-aug: track augments since we continue past them
  bool phaseAug = false;       // augmented at the CURRENT dual level (since last relabel)
  LL augTime = -1;             // dual level (tNow) of this phase's augments
  bool phaseHadShrink = false; // a shrink happened at THIS dual level (since last relabel). If so, same-T
                               // augments can thread reshaped structure -> must serialize (end phase at ==T).
  Sz phasesDone = 0;           // completed MV multi-augment phases (for mPhaseBudget)
  mLastPhaseCount = 0;
  for (;;) {
    if (++guard > guardCap) break;
#ifdef DIVDBG
    { static Sz s_lastGuard = 0; Sz cy = findTreeCycle();
      if (cy != NIL()) { std::fprintf(stderr,"[CYCLE-DETECTED node=%lu at guard=%lld (loop top)]\n",
        (unsigned long)cy,(long long)guard); std::fflush(stderr); return false; }
      { Sz bi = 0; Sz mB = findCycEdgeMisalign(&bi);
        static bool s_misReported = false;
        if (mB != NIL() && !s_misReported) { s_misReported = true;
          std::fprintf(stderr,"[CYC-EDGE-MISALIGN B=%lu idx=%lu k=%lu at guard=%lld tNow-ish | cyc[i]=%lu eu[i]=%lu ev[i]=%lu]\n",
            (unsigned long)mB,(unsigned long)bi,(unsigned long)mBCycle[mB].size(),(long long)guard,
            (unsigned long)mBCycle[mB][bi],(unsigned long)mBCycEU[mB][bi],(unsigned long)mBCycEV[mB][bi]);
          std::fflush(stderr); } }
      s_lastGuard = (Sz)guard; (void)s_lastGuard; }
#endif
#ifdef MATE_AUDIT
    // report the FIRST guard iteration where the matching became asymmetric (do NOT return -- let the
    // per-op ASYM-AFTER-* checks pinpoint the culprit op first).
    { Sz fb; Sz asym=mateAsym(&fb);
      static bool wasClean=true;
      if(asym>0 && wasClean){ wasClean=false; std::fprintf(stderr,"[MATE-FIRST-ASYM guard=%lld asym=%lu firstBad=%ld]\n",
        (long long)guard,(unsigned long)asym,(long)fb); std::fflush(stderr); } }
#endif
#ifdef PS5_VB_DEFER_REPAIR
    if (!mVbRepList.empty()) {
      bool hp0 = mPQ.advanceToNext();
      if (!hp0 || mPQ.now() > tNow) { vbRepairFlush(tNow); continue; }  // tick done -> restore coverage, re-evaluate
    }
#endif
    bool havePQ = mPQ.advanceToNext();
    LL pqT = havePQ ? mPQ.now() : -1;
#ifdef DIVDBG
    // Divergence detector: once past a small cap, periodically dump state. If tNow / #free / #even-roots
    // stay ~constant across dumps -> INFINITE LOOP (not progressing); if tNow climbs -> slow-but-moving.
    if (guard > 50000 && (guard % 20000) == 0) {
      Sz nfree = 0, neven = 0;
      for (Sz v = 0; v < mNReal; ++v) if (mMate[v] == NIL()) ++nfree;
      for (Sz i = 0; i < mNumNodes; ++i) if (mBParent[i] == NIL() && mLabel[i] == LBL_EVEN) ++neven;
      std::fprintf(stderr,"[DIV guard=%lld tNow=%lld pqT=%lld havePQ=%d free=%lu evenRoots=%lu numNodes=%lu phaseAug=%d augTime=%lld]\n",
        (long long)guard,(long long)tNow,(long long)pqT,(int)havePQ,(unsigned long)nfree,(unsigned long)neven,(unsigned long)mNumNodes,(int)phaseAug,(long long)augTime);
      if (guard > 250000) { std::fprintf(stderr,"[DIV BAILOUT at guard=%lld -- stopping to inspect]\n",(long long)guard); break; }
    }
#endif
#ifdef RELINKDBG
    if (gS5NrdSeen) { static long long lc=0; ++lc; if(lc%5000==0) std::fprintf(stderr,"[LOOP#%lld guard=%lld tNow=%lld havePQ=%d pqT=%lld phaseAug=%d augTime=%lld numNodes=%lu]\n",(long long)lc,(long long)guard,(long long)tNow,(int)havePQ,(long long)pqT,(int)phaseAug,(long long)augTime,(unsigned long)mNumNodes); }
#endif
    // PER-LEVEL MULTI-AUGMENT PHASE BOUNDARY (the MV/HK phase = the sqrt(n) mechanism):
    // we multi-augment all vertex-disjoint paths AT THE SAME dual level (augTime). The
    // instant the next event would advance time past augTime (a later PQ event OR a
    // dissolution), the phase ends: relabel + rebuild the PQ ONCE (the per-phase O(m)),
    // clear the disjointness marks, and restart at the new level. O(sqrt n) phases total.
    // SUBOPTIMALITY FIX (2026-07-02, DEFAULT-ON): normally the phase ends when time advances PAST the augment
    // level (pqT > augTime), batching all same-T augments. But if this level did a SHRINK, a same-T augment
    // can thread the reshaped forest -> stale roles -> negative slack -> suboptimal (trace n=80 W3 s3). So when
    // phaseHadShrink, end the phase at pqT >= augTime (serialize same-T ops). Pure-augment levels (no shrink)
    // keep batching -> the fast greedy phase is unaffected. Disable with -DPS5_NO_SHRINK_SERIALIZE. Also
    // PS5_PHASE_SAMET forces the unconditional >= (the proven-correct-but-slow baseline, 5/6).
#if defined(PS5_PHASE_SAMET)
    bool endAtSameT = true;
#elif !defined(PS5_NO_SHRINK_SERIALIZE)
    bool endAtSameT = phaseHadShrink;   // DEFAULT: end phase at same-T only if this level did a shrink/diss
#ifndef PS5_NO_PHASE_AUGCAP
    // ONE-TICK MEGA-PHASE CAP (2026-07-13, the (51200,W=16,seed0) >900x tail): at high tie
    // density an entire instance can fire at ONE dual level, so pqT never exceeds augTime and
    // the phase batches ~n/2 augmentations with per-augment cost that grows as the phase ages
    // (quadratic batch; DIVDBG: tNow=2 for the whole run). Ending a phase early is ALWAYS sound
    // (budgets/overrun/shrink-serialize already do it), and a same-T end takes the PERSISTENT
    // boundary -- so the cap just breaks the mega-phase into linear chunks. Cap = 4*sqrt(n)+64:
    // far above the natural augs/phase (~20 distinct, low hundreds tie-dense), so ordinary
    // phases never hit it. Opt out: -DPS5_NO_PHASE_AUGCAP.
    if (phaseAugs >= mPhaseAugCap) endAtSameT = true;
#endif
#else
    bool endAtSameT = false;
#endif
    if (mNoRelabelAug && phaseAug && ((!havePQ) || pqT > augTime || (endAtSameT && pqT == augTime))) {
#if defined(PS5_PERSIST_FOREST) && defined(PS5_PERSIST_CLOCK)
      // STAGE 2 (docs/H_Scope_PersistentClock.md): PERSISTENT boundary. On a pure-augment phase
      // the forest, the event queue, and the clock all SURVIVE: consumed trees were retired (and
      // their boundary edges rescheduled at the new FREE rate) inside freeConsumedTree at augment
      // time; surviving trees' events keep their absolute fire times under the monotone clock.
      // Only the phase-scoped bookkeeping resets. Shrink phases keep the proven teardown below
      // (the shrink-serialize correctness class needs the relabel cleaner) unless
      // PS5_PERSIST_SHRINK extends persistence to them (the shrink hazard is same-TICK augment
      // threading WITHIN a phase; at the boundary the level is over -- gate-tested).
#ifdef PS5_PERSIST_SHRINK
      if (true) {
#else
      if (!phaseHadShrink) {
#endif
#ifdef PS5_QSTAT
        ++mQPhases; ++mQAdvHist[qstatBucket(augTime < 0 ? 0 : augTime)];
#endif
        ++phasesDone; mLastPhaseCount = phasesDone;
        if (mPhaseBudget != 0 && phasesDone >= mPhaseBudget) {
          if (mLazyActive) lazyMaterialize();
          mLazyActive = false; return augmentedAny;
        }
#ifdef PS5_EVTRACE
        std::fprintf(stderr,"[TR] BOUNDARY tNow=%lld augTime=%lld reoffers=%lu\n",
          (long long)tNow,(long long)augTime,(unsigned long)mReoffer.size());
#endif
        // The PQ cursor may have PEEKED ahead (advanceToNext ran to find pqT > augTime); inserts
        // below are at the TRUE clock (tNow) and must not be clamped forward past their fire time.
        mPQ.syncNow(tNow);
        retireDeadTrees();   // batch: settle dead trees + reschedule their boundary edges
#ifdef DIVDBG
        { Sz bi2=0; Sz mB2 = findCycEdgeMisalign(&bi2); static bool s_r=false;
            if (mB2 != NIL() && !s_r) { s_r=true; std::fprintf(stderr,"[MISALIGN-AFTER-RETIRE B=%lu idx=%lu k=%lu]\n",
              (unsigned long)mB2,(unsigned long)bi2,(unsigned long)mBCycle[mB2].size()); std::fflush(stderr); } }
        { Sz cy = findTreeCycle();
          if (cy != NIL()) { std::fprintf(stderr,"[CYCLE-AFTER-RETIRE at node=%lu tNow=%lld]\n",
            (unsigned long)cy,(long long)tNow); std::fflush(stderr); return false; } }
#endif
        mUsed.assign(mN, 0);
        // Re-offer events that were mUsed-blocked this phase (still eligible; pop re-validates).
        for (Sz ri = 0; ri < mReoffer.size(); ++ri) {
          const BucketPQ::Event& rev = mReoffer[ri];
          mPQ.schedule(rev.time, rev.a, rev.b, rev.kind);
        }
        mReoffer.clear();
        phaseAug = false; augTime = -1; phaseHadShrink = false; phaseAugs = 0;
#ifdef S5DBG
        { extern long long gS5PersistPhases; ++gS5PersistPhases; }
#endif
        continue;
      }
#endif
      // time would advance past this phase's augment level -> END the phase: relabel +
      // rebuild ONCE, clear marks, restart at the new level. (The per-phase O(m).)
#ifdef S5DBG
      ++gS5PhaseRelabels;
#endif
#ifdef PS5_QSTAT
      ++mQPhases; ++mQAdvHist[qstatBucket(augTime < 0 ? 0 : augTime)];
#endif
#ifdef PS5_STAGE0
      s0Snapshot();   // end-of-phase forest = the baseline the next phase regrows against
#endif
      if (mLazyActive) lazyMaterialize();
      // PHASE BUDGET (the tau bound): one completed MV phase. If we have hit the budget, STOP here
      // -- duals are materialized (mY/mZ current), the matching is valid, and any still-free reals
      // stay free for Phase B to park on dummies. This is the paper-faithful "tau SearchOne's per
      // scale" cap (Fig 5 step 7): bounding the NUMBER OF DUAL-LEVEL PHASES, not searchEvent calls.
      ++phasesDone; mLastPhaseCount = phasesDone;
      if (mPhaseBudget != 0 && phasesDone >= mPhaseBudget) { mLazyActive = false; return augmentedAny; }
#ifdef PS5_INVKEY
      ikSnapshot();   // pre-relabel EVEN set = what the surviving PQ's events cover
#endif
      relabel();
      if (mLazyActive) lazyBegin();
      mUsed.assign(mN, 0);
      phaseAug = false; augTime = -1; phaseHadShrink = false; phaseAugs = 0; tNow = mTNow;
#if defined(S5DBG) || defined(PS5_BIGNODE_CACHE)
      ++mDbgEpoch;   // phase-boundary relabel+reseed = new epoch
#endif
#ifdef PS5_INVKEY
      // INVARIANT-KEY boundary (N_Scope Stage 2): keep the PQ (fire times are absolute under the
      // monotone clock; lazyBegin was value-preserving), reschedule only newly-EVEN nodes. VTX_BEST
      // epoch NOT bumped: the queued bests survive with the queue.
      ikReschedChanged(tNow);
#else
      mPQ.reset();   // empty buckets but KEEP capacity (clear()+init(1) reallocated the
                     // whole bucket array every relabel -- the ~44% malloc/resize storm).
      reseedAllEven(tNow, 1);   // phase-boundary relabel reseed
#endif
#ifdef PS5_DISHEAP
      disHeapRebuild(tNow);   // phase boundary: rebuild dissolution heap at the new dual level
#endif
      continue;
    }
    // next dissolution time over ODD root blossoms (materialized z is exact here:
    // we recompute at the moment of use, and ODD roots' z only drops as we advance).
#ifdef S5DBG
    ++gS5Iters;
#ifndef PS5_DISHEAP
    gS5DisScanWork += (long long)(mNumNodes - mN);   // only the real scan path costs this
#endif
#endif
    LL disT = -1; Sz disB = NIL();
#ifdef PS5_DISHEAP
    // O(log n) amortized via the dissolution min-heap (lazy-validated) instead of the O(blossoms) scan.
    { LL t; Sz b = disHeapNext(tNow, &t); if (b != NIL()) { disT = t; disB = b; } }
    // (DISHEAP_SHADOW cross-check confirmed heap==scan on the min dissolution; removed.)
#else
    for (Sz b = mN; b < mNumNodes; ++b)
      if (mBParent[b] == NIL() && mLabel[b] == LBL_ODD) {
        LL zb = mLazyActive ? zLazy(b) : mZ[b];      // lazy z is exact across role changes
        LL t = tNow + zb / 2;
        if (disT < 0 || t < disT) { disT = t; disB = b; }
      }
#endif
    // NEGATIVE-Z CLAMP (the seed=0 divergence fix): a dissolution time disT = tNow + z/2 can be LESS than
    // tNow when z has gone negative (an over-tight ODD blossom whose z dropped past 0 without dissolving,
    // in the carried-dual / Property-2 regime). Firing at a PAST time drives applyDualTo BACKWARD (tNow
    // spirals to -infinity, free stuck, guard-cap failure). A due dissolution fires NOW, never in the
    // past: clamp disT up to tNow.
    if (disB != NIL() && disT < tNow) disT = tNow;
    if (!havePQ && disB == NIL()) break;            // no event at all -> done
    // pick the earliest of (next PQ event, next dissolution).
    bool doDissolve = (disB != NIL()) && (!havePQ || disT <= pqT);
    // TICK BUDGET (paper Fig 5 step 7 semantics): stop this search once the scale's total dual
    // adjustment would exceed the budget. Events AT the current time (zero advance) still process --
    // the paper's SearchOne also completes its maximal augmenting set at the current duals before
    // adjusting. Exiting mid-phase is safe for the scaling-loop caller: duals are materialized here
    // and the driver rebuilds the forest per scale (resetForestKeepDuals).
    if (mTickBudget > 0) {
      LL nxt = doDissolve ? disT : pqT;
      if (nxt > tNow) {
        LL adv = nxt - tNow;
        if (mTicksUsed + adv > mTickBudget) {
          // The budget is a number of dual-adjustment units, not a permission to process only
          // complete future events.  SearchOne must still advance through the remaining units
          // when the next tight edge lies farther away, then return at that intermediate dual
          // level.  The old early return consumed zero ticks in this case, so Figure-6 stage 1
          // frequently did no dual work at all.
          LL remain = mTickBudget - mTicksUsed;
          if (remain > 0) {
            if (mLazyActive) lazySyncRoles();
            applyDualTo(tNow + remain, &tNow);
            mTicksUsed += remain;
          }
          if (mLazyActive) { lazyMaterialize(); mLazyActive = false; }
          return augmentedAny;
        }
        mTicksUsed += adv;
      }
    }
#ifdef DIVDBG
    { static long long dc=0; if(doDissolve && ((++dc)%20000)==0)
        std::fprintf(stderr,"[DISSOLVE-SPIN #%lld disB=%lu disT=%lld tNow=%lld pqT=%lld zLazy=%lld label=%d bparent=%ld]\n",
          dc,(unsigned long)disB,(long long)disT,(long long)tNow,(long long)pqT,(long long)zLazy(disB),mLabel[disB],(long)mBParent[disB]); }
#endif
    if (doDissolve) {
      if (mLazyActive && disT > tNow) lazySyncRoles();   // settle ONLY when time advances
      applyDualTo(disT, &tNow);
      // S5-C STEP 4: NO-RELABEL in-place dissolution. Re-thread the tree through each z=0 ODD
      // root blossom's cycle (relinkOddCycleS5) instead of materialize+relabel+rebuild-PQ.
      // lazySyncRoles settles the lazy duals for the new labels; only the freed sub-nodes' edges
      // are (re)scheduled, not all m. Falls back to the relabel path if any relink declines.
      if (mNoRelabelDissolve && mLazyActive) {
#ifdef MATE_AUDIT
        { Sz fb; Sz a=mateAsym(&fb); static int r=0; if(a>0 && r++<3)
          std::fprintf(stderr,"[ASYM-BEFORE-DISSOLVE-BATCH asym=%lu firstBad=%ld]\n",(unsigned long)a,(long)fb); }
#endif
        // ALL-OR-NOTHING: pre-check EVERY z=0 ODD root is relink-eligible BEFORE mutating. A
        // partial in-place batch (relink+dissolve some, then a later one declines -> fall to
        // relabel) leaves the half-expanded ones inconsistent (the n=30/32 suboptimal bug). Only
        // commit the in-place path if all qualify; else relabel the whole batch (proven path).
#ifdef PS5_DEEP_INTACT
        // fresh memo epoch for this dissolution batch (struct-eligibility is invariant within it).
        if (mSEMemo.size() < mNumNodes) { mSEMemo.resize(mNumNodes, 0); mSEMemoStamp.resize(mNumNodes, -1); }
        ++mSEMemoTok;
#endif
        bool allEligible = true;
        for (Sz b = mN; b < mNumNodes && allEligible; ++b)
          if (mBParent[b] == NIL() && mLabel[b] == LBL_ODD && zLazy(b) <= 0)
            if (!relinkEligibleS5(b)) allEligible = false;
        std::vector<Sz> relinkedNodes;     // sub-nodes that became EVEN (need edge scheduling)
#ifdef PS5_DISHEAP
        std::vector<Sz> newOddRoots;       // sub-nodes that became ODD roots (new dissolution candidates)
#endif
        bool allInPlace = allEligible;
        if (allEligible)
          for (Sz b = mN; b < mNumNodes; ++b)
            if (mBParent[b] == NIL() && mLabel[b] == LBL_ODD && zLazy(b) <= 0) {
              std::vector<Sz> kids = mBCycle[b];   // snapshot before relink/dissolve mutate
#ifdef DISS_TRACE
              gDissKids = kids;
              { std::fprintf(stderr,"[DISS b=%lu z=%lld cycle=", (unsigned long)b,(long long)zLazy(b));
                for(Sz i=0;i<kids.size();++i) std::fprintf(stderr,"%lu(base%lu,mate%ld) ",(unsigned long)kids[i],(unsigned long)mBBase[kids[i]],(long)(mBBase[kids[i]]<mN?mMate[mBBase[kids[i]]]:-1));
                std::fprintf(stderr,"| tNow-region mates:"); for(Sz i=0;i<kids.size();++i){Sz bs=mBBase[kids[i]]; if(bs<mN)std::fprintf(stderr," %lu-%ld",(unsigned long)bs,(long)mMate[bs]);} std::fprintf(stderr,"]\n"); }
#endif
              bool ok = relinkOddCycleS5(b);
              if (!ok) { allInPlace = false; break; }   // should not happen (pre-checked) -- safety
              for (Sz i = 0; i < kids.size(); ++i)
                if (mLabel[kids[i]] == LBL_EVEN) relinkedNodes.push_back(kids[i]);
#ifdef PS5_DISHEAP
              for (Sz i = 0; i < kids.size(); ++i)
                if (kids[i] >= mN && mBParent[kids[i]] == NIL() && mLabel[kids[i]] == LBL_ODD)
                  newOddRoots.push_back(kids[i]);
#endif
              mLabel[b] = LBL_FREE;
              dissolve(b);
#ifdef PS5_NESTED_RELINK
              // NESTED (WIP): recursively relink+dissolve any z=0 sub-blossom now exposed as a root.
              for (Sz i = 0; i < kids.size() && allInPlace; ++i) {
                Sz sn = kids[i];
                if (sn >= mN && mBParent[sn] == NIL() && zLazy(sn) <= 0) {
                  if (!relinkRecDissolve(sn, &relinkedNodes)) { allInPlace = false; break; }
                }
              }
              if (!allInPlace) break;
#else
              (void)kids;
#endif
            }
        if (allInPlace) {
#ifdef DIVDBG
          { Sz bi2=0; Sz mB2 = findCycEdgeMisalign(&bi2); static bool s_r=false;
            if (mB2 != NIL() && !s_r) { s_r=true; std::fprintf(stderr,"[MISALIGN-AFTER-NRDISSOLVE B=%lu idx=%lu k=%lu]\n",
              (unsigned long)mB2,(unsigned long)bi2,(unsigned long)mBCycle[mB2].size()); std::fflush(stderr); } }
#endif
#ifdef PS5_UF_CHECK
          ufCheck("dissolve");
#endif
#ifdef PS5_UF_RESYNC
          ufResyncFromTree();   // heal the uf<->mBParent desync introduced by the nested dissolution
#endif
          lazySyncRoles();                 // settle lazy duals for all the relink label changes
#ifndef PS5_NRD_NODUALFIX
          // DUAL-FIX (correctness, default-on): after the in-place dissolution, REFRESH the lazy-dual base
          // + tNow + PQ consistently (like the relabel path's materialize+lazyBegin+reseed, but WITHOUT the
          // re-root, so the no-relabel tree is kept). The incremental lazy-dual settling across dissolutions
          // DRIFTS -> dual-INFEASIBLE y/z (n=20 seed2: edge 7-10 slack=-46) -> the search stops at a
          // SUBOPTIMAL perfect matching. Re-snapshotting restores complementary slackness. MEASURED: fixes
          // ~28% of distinct-weight suboptimality (LEMON 162/164 optimal), does NOT reintroduce the uniform
          // cliff (dissolutions batch on uniform; n=12800 stays 0.73s, beats GMG), 348 gate green.
          // (Disable via PS5_NRD_NODUALFIX to A/B the drift.)
          lazyMaterialize();
#ifdef PS5_NRD_RELABEL
          // FLAVOR-2 TEST: the no-relabel dissolution frees off-arc cycle nodes (LBL_FREE) that the RELABEL
          // path would re-root into OTHER trees -> augments through them are lost -> suboptimal (n=80 W5 s0).
          // relabel() re-roots the whole forest, matching the relabel path's reachability. Tests whether the
          // off-arc re-integration is what flavor 2 needs (at the cost of the no-relabel dissolution's point).
          relabel();
#endif
          lazyBegin();                         // fresh consistent dual base, mTNow=0
          tNow = mTNow; mPQ.reset();
          reseedAllEven(tNow, 2);   // NRD dual-fix reseed (after successful in-place dissolution)
#ifdef PS5_DISHEAP
          disHeapRebuild(tNow);
#endif
          continue;   // skip the normal relinkedNodes-only reschedule below
#endif
#ifdef PS5_DISHEAP
          for (Sz i = 0; i < newOddRoots.size(); ++i)
            if (mBParent[newOddRoots[i]] == NIL() && mLabel[newOddRoots[i]] == LBL_ODD)
              disHeapSchedule(newOddRoots[i], tNow);   // new ODD roots -> future dissolution candidates
#endif
#ifdef NRD_PROBE
          { extern long long gNRDsucc, gNRDrelinked; ++gNRDsucc; gNRDrelinked += (long long)relinkedNodes.size(); }
#endif
          for (Sz i = 0; i < relinkedNodes.size(); ++i)
            scheduleNodeEdges(relinkedNodes[i], tNow);   // only the newly-EVEN nodes' edges
#ifdef PS5_DISS_RESCHED_ALL
          // RESIDUAL FIX TEST (150 W5 s3, a dual-overrun): a no-relabel dissolution can leave an EVEN node's
          // boundary edge unscheduled (its slack changed but it was not in relinkedNodes) -> the dual clock
          // overruns that edge -> z overgrows -> neg slack. Reschedule ALL current EVEN nodes' edges so no
          // boundary event is missing. (O(m) per dissolution -- dissolutions are rare, ~few per solve.)
          reseedAllEven(tNow, 4);
#endif
#ifdef DISS_TRACE
          { std::fprintf(stderr,"[POST-DISS relinkedEVEN={");
            for(Sz i=0;i<relinkedNodes.size();++i) std::fprintf(stderr,"%lu ",(unsigned long)relinkedNodes[i]);
            std::fprintf(stderr,"} freed-sub-node labels:");
            for(Sz i=0;i<gDissKids.size();++i){Sz sn=gDissKids[i]; std::fprintf(stderr," %lu(l%d,par%ld,base%lu,mate%ld)",(unsigned long)sn,mLabel[sn],(long)mParent[sn],(unsigned long)(sn>=mN?mBBase[sn]:sn),(long)(mBBase[sn]<mN?mMate[mBBase[sn]]:-1));}
            std::fprintf(stderr,"]\n"); }
#endif
#ifdef RELINKDBG
          { gS5NrdSeen=true; static long long dc=0; if(++dc<=60) std::fprintf(stderr,"[NRD] tNow=%lld relinked %lu even-nodes, disT=%lld havePQ=%d pqT=%lld\n",(long long)tNow,(unsigned long)relinkedNodes.size(),(long long)disT,(int)havePQ,(long long)pqT); }
#endif
#if defined(S5DBG) || defined(PS5_BIGNODE_CACHE)
          ++mDbgEpoch;
#endif
#ifdef DIVDBG
          { Sz cy = findTreeCycle();
            if (cy != NIL()) { std::fprintf(stderr,"[CYCLE-AFTER-DISSOLVE at node=%lu disB=%lu]\n",
              (unsigned long)cy,(unsigned long)disB); std::fflush(stderr); return false; } }
#endif
#ifdef MATE_AUDIT
          { Sz fb; Sz a=mateAsym(&fb); static int r=0; if(a>0 && r++<2) {
            Sz p = (fb<mN)?mMate[fb]:NIL();
            std::fprintf(stderr,"[ASYM-AFTER-DISSOLVE disB=%lu asym=%lu firstBad=%ld mate[%ld]=%ld mate[%ld]=%ld]\n",
              (unsigned long)disB,(unsigned long)a,(long)fb,(long)fb,(long)p,(long)p,(long)(p!=NIL()&&p<mN?mMate[p]:-1));
            // dump the freed cycle(s): which vertex re-paired fb spuriously
            for(Sz bb=mN; bb<mNumNodes; ++bb) if(!mBCycle[bb].empty()){ bool has=false;
              for(Sz i=0;i<mBCycle[bb].size();++i){Sz sn=mBCycle[bb][i]; if(sn<mN&&(sn==fb||(p!=NIL()&&sn==(Sz)p))) has=true; if(sn>=mN){for(Sz j=0;j<mNodeLeaves[sn].size();++j)if(mNodeLeaves[sn][j]==fb||(p!=NIL()&&mNodeLeaves[sn][j]==(Sz)p))has=true;}}
              if(has){ std::fprintf(stderr,"   [cycle of node %lu contains fb/p, len=%lu bparent=%ld]\n",(unsigned long)bb,(unsigned long)mBCycle[bb].size(),(long)mBParent[bb]); } }
          } }
#endif
          phaseHadShrink = true;   // no-relabel DISSOLUTION also reshapes the forest -> serialize same-T augs
          continue;                        // NO relabel, NO full reschedule, NO tNow reset
        }
        // a relink declined: fall through to the proven relabel path (correct, just O(m)).
      }
      // S5-C STEP 1 fallback: MATERIALIZE the lazy duals, relabel, and rebuild the PQ.
      if (mLazyActive) lazyMaterialize();            // write lazy y/z -> mY/mZ for dissolve+relabel
#ifdef PS5_INVKEY
      ikSnapshot();   // BEFORE the dissolve loop: the queued events cover the pre-dissolve frame
#endif
      bool any = false;
      for (Sz b = mN; b < mNumNodes; ++b)
        if (mBParent[b] == NIL() && mLabel[b] == LBL_ODD && mZ[b] == 0) { dissolve(b); any = true; }
#ifdef DIVDBG
      { Sz bi2=0; Sz mB2 = findCycEdgeMisalign(&bi2); static bool s_r=false;
            if (mB2 != NIL() && !s_r) { s_r=true; std::fprintf(stderr,"[MISALIGN-AFTER-RELDISSOLVE B=%lu idx=%lu k=%lu]\n",
              (unsigned long)mB2,(unsigned long)bi2,(unsigned long)mBCycle[mB2].size()); std::fflush(stderr); } }
#endif
#ifdef S5DBG
      ++gS5DissolveRelabels;
#endif
      relabel();
      if (mLazyActive) lazyBegin();                  // re-snapshot base + seed roles post-relabel
      tNow = mTNow;
#if defined(S5DBG) || defined(PS5_BIGNODE_CACHE)
      ++mDbgEpoch;   // dissolution relabel+reseed = new epoch
#endif
#ifdef PS5_INVKEY
      ikReschedChanged(tNow);   // kept-PQ boundary: reschedule only newly-EVEN (Stage-0: 100% repeat here)
#else
      mPQ.reset();   // KEEP capacity (was clear()+init(1) = realloc storm)
      // The post-dissolution RESEED (re-schedule all EVEN nodes' edges) -- MEASURED to be ~100% of the
      // uniform-weight cliff (n=12800: 37.1s -> 0.03s if skipped; reseed-ON doubling 16x/7x super-linear,
      // reseed-OFF 1.7x near-linear). This is what the split-findmin findmin-driven scheduling replaces.
      // (2026-07-02 callgrind: scheduleNodeEdges ~19.5% of a distinct-weight run; but split-findmin does
      // NOT close the complexity + does NOT fix the residuals -- see docs/G_Scope_SplitFindminDriver.md.)
      reseedAllEven(tNow, 3);   // relink-decline fallback relabel reseed
#endif
#ifdef PS5_DISHEAP
      disHeapRebuild(tNow);   // relabel fallback: rebuild the dissolution heap
#endif
#ifdef PS5_SFM
      // CONSISTENCY GATE: split-findmin findmin must agree with the reseed's next-event (no behavior
      // change -- the reseed above still runs; this only validates the query before it replaces it).
      sfmCheckAfterReseed(tNow);
#endif
#ifdef PS5_CQ_EVENTS
      // K_Scope Stage 1a REFEREE: per-blossom CQ global-min must equal the brute next-fire (shadow only).
      cqCheckAfterReseed(tNow);
#endif
      if (!any) {
        if (!mPQ.advanceToNext()) break;
      }
      continue;
    }

    LL T = pqT;
    BucketPQ::Event ev;
#ifdef DIVDBG
    { static long long popOk=0, popFail=0; bool got = mPQ.popAtNow(&ev);
      if (got) ++popOk; else ++popFail;
      if (((popOk+popFail) % 20000)==0) std::fprintf(stderr,"[PQ popOk=%lld popFail=%lld pqT=%lld tNow=%lld ev=(%ld,%ld)]\n",
        popOk,popFail,(long long)pqT,(long long)tNow,(long)ev.a,(long)ev.b);
      if (!got) continue; }
#else
    if (!mPQ.popAtNow(&ev)) continue;
#endif
#ifdef PS5_QSTAT
    ++mQPops;
#endif
#ifdef PS5_EVTRACE
    std::fprintf(stderr,"[TR] pop T=%lld kind=%d a=%ld b=%ld tNow=%lld\n",
      (long long)T,ev.kind,(long)ev.a,(long)ev.b,(long long)tNow);
#endif
#ifdef PS5_EVT_HORIZON
    // Marker (kind 3): rescan the node's edges with the window re-based at the marker time. Runs
    // BEFORE the overrun guard / mUsed / applyDualTo -- markers carry no endpoint pair (ev.b=-1),
    // must not advance the clock, and a stale marker (node absorbed/label changed) simply drops:
    // whatever changed the node already rescheduled the covering node's edges (the existing
    // structural-change invariant).
    if (ev.kind == 3) {
      // Walk to the node's CURRENT top-level: no-relabel shrink absorbs nodes WITHOUT rescanning
      // (the stale member events in the PQ are the coverage, validated at pop via outerNode) --
      // so a marker on an absorbed member must re-aim at the containing blossom, not drop (the
      // n=4 rep=0 microscope: B={0,1,3} formed, member markers dropped, the (3,2) augment event
      // was never re-created). A non-EVEN top drops the marker, exactly like a stale event would
      // be discarded at pop: every transition back to EVEN comes with its own scan.
      Sz En = (Sz)ev.a;
      if (En < mNumNodes) {
        Sz On = En;
        while (On < mNumNodes && mBParent[On] != NIL()) On = mBParent[On];
        if (mLabel[On] == LBL_EVEN) {
#ifdef PS5_MARKER_ADAPT
          // SELF-TUNING no-defer (2026-07-08): a node whose marker keeps popping is in the giant-
          // blossom RESCAN-thrash regime (pathological n=150000: one blossom's marker popped 29155x,
          // each an O(leaves) rescan = the cliff). Count per-node marker pops; once a node exceeds
          // kThresh, flag it so scheduleNodeEdges stops deferring it (schedules its events into a
          // LARGE-BUT-BOUNDED window -- NOT 2^60, which floods the multiset overflow). Benign nodes
          // pop 0-1x and never trip it -> the horizon's win on the ring+3n family is preserved.
          if (mMarkerPops.size() <= On) mMarkerPops.assign(mNumNodes+8, 0);
          if (++mMarkerPops[On] > (unsigned)(PS5_MARKER_ADAPT)) {
            if (mNoDeferNode.size() <= On) mNoDeferNode.assign(mNumNodes+8, 0);
            mNoDeferNode[On] = 1;
          }
#endif
          scheduleNodeEdges(On, tNow, T);
        }
#ifdef PS5_EVTRACE
        else std::fprintf(stderr,"[TR] MARKER-DROP node=%lu top=%lu T=%lld label=%d\n",
          (unsigned long)En,(unsigned long)On,(long long)T,mLabel[On]);
#endif
      }   // En < mNumNodes
      continue;
    }
#endif
#ifdef PS5_OVERRUN_LOG
    // OVERRUN DETECTOR: compute the true next-fire (min over all EVEN-incident external edges of the time
    // they reach the target). If the popped event time T EXCEEDS it, the PQ is MISSING that edge's event ->
    // log the culprit edge + its type. This names WHICH event scheduling is incomplete (the real bug).
    if (T > tNow) {
      LL bf = -1; Sz be_u = NIL(), be_w = NIL(); int be_lw = -1; LL tgt = (mCrit == CRIT1) ? 0 : -2;
      for (Sz En = 0; En < mNumNodes; ++En) {
        if (mBParent[En] != NIL() || mLabel[En] != LBL_EVEN) continue;
        std::vector<Sz>& Ev = mNodeLeaves[En];
        for (Sz ui = 0; ui < Ev.size(); ++ui) { Sz uu = Ev[ui]; Sz ufu = ufFind(uu);
          for (Sz kk = 0; kk < mAdjV[uu].size(); ++kk) { Sz ww = mAdjV[uu][kk];
            if (mIgnoreDummies && ww >= mNReal) continue;
            if (ufu == ufFind(ww)) continue;
            Sz Wn = outerNode(ww); int lw = mLabel[Wn];
            if (lw != LBL_FREE && !treeLive(Wn)) lw = LBL_FREE;
            if (lw != LBL_EVEN && lw != LBL_FREE) continue;
            LL s = slackZW(uu, ww, 2 * mAdjW[uu][kk]);
            LL d = (lw == LBL_EVEN) ? (s - tgt) / 2 : (s - tgt);
            if (!mAllowOvertight && d < 0) continue; if (mAllowOvertight && d < 0) d = 0;
            LL fire = tNow + d; if (bf < 0 || fire < bf) { bf = fire; be_u = uu; be_w = ww; be_lw = lw; } } } }
      if (bf >= 0 && T > bf) { static int r=0; if (r++ < 12) {
        Sz oU = outerNode(be_u), oW = outerNode(be_w);
        std::fprintf(stderr,"   [UF-CHECK ufFind(%lu)=%lu ufFind(%lu)=%lu same=%d | outerNode says same=%d]\n",
          (unsigned long)be_u,(unsigned long)ufFind(be_u),(unsigned long)be_w,(unsigned long)ufFind(be_w),
          (int)(ufFind(be_u)==ufFind(be_w)),(int)(oU==oW));
        std::fprintf(stderr,"[OVERRUN T=%lld > trueFire=%lld by %lld | MISSING edge (%lu,%lu) type=%s | outer(%lu)=%lu[lbl=%d root=%ld] outer(%lu)=%lu[lbl=%d root=%ld] slack=%lld tNow=%lld]\n",
          (long long)T,(long long)bf,(long long)(T-bf),(unsigned long)be_u,(unsigned long)be_w,
          be_lw==LBL_EVEN?"FORMATION":"GROW",(unsigned long)be_u,(unsigned long)oU,mLabel[oU],(long)mRoot[oU],
          (unsigned long)be_w,(unsigned long)oW,mLabel[oW],(long)mRoot[oW],
          (long long)slackZW(be_u,be_w,2*edgeW(be_u,be_w)),(long long)tNow); } }
    }
#endif
#ifdef PS5_CLAMP_DUAL
    // DIAGNOSTIC (NOT a fix -- confirmed the residual mechanism): clamp the dual advance to the ground-truth
    // next-fire (inline O(m) scan) + re-pop the too-early event. RESULT: fixes n=150 W5 s3 (== LEMON) ->
    // confirms 150 is a dual-clock OVERRUN (the PQ is missing a boundary-edge event so the clock advances
    // past when that edge went tight -> blossom overgrows its z -> negative slack). But n=100 W50 s3 then
    // MISSES a perfect matching (S5=-1) -> 100 is a DIFFERENT bug, not an overrun. So the 2 residuals are
    // distinct: 150 = missing-event overrun (real fix = schedule that event); 100 = separate. Gated OFF.
    { LL bf = -1; LL tgt = (mCrit == CRIT1) ? 0 : -2;
      for (Sz En = 0; En < mNumNodes; ++En) {
        if (mBParent[En] != NIL() || mLabel[En] != LBL_EVEN) continue;
        std::vector<Sz>& Ev = mNodeLeaves[En];
        for (Sz ui = 0; ui < Ev.size(); ++ui) { Sz uu = Ev[ui]; Sz ufu = ufFind(uu);
          for (Sz kk = 0; kk < mAdjV[uu].size(); ++kk) { Sz ww = mAdjV[uu][kk];
            if (mIgnoreDummies && ww >= mNReal) continue;
            if (ufu == ufFind(ww)) continue;
            Sz Wn = outerNode(ww); int lw = mLabel[Wn];
            if (lw != LBL_FREE && !treeLive(Wn)) lw = LBL_FREE;
            if (lw != LBL_EVEN && lw != LBL_FREE) continue;
            LL s = slackZW(uu, ww, 2 * mAdjW[uu][kk]);
            LL d = (lw == LBL_EVEN) ? (s - tgt) / 2 : (s - tgt);
            if (!mAllowOvertight && d < 0) continue; if (mAllowOvertight && d < 0) d = 0;
            LL fire = tNow + d; if (bf < 0 || fire < bf) bf = fire;
          } } }
      if (bf >= 0 && T > bf) {
        if (mLazyActive && bf > tNow) lazySyncRoles();
        applyDualTo(bf, &tNow);
        mPQ.schedule(T, (long long)ev.a, (long long)ev.b, ev.kind);
        continue;
      }
    }
#endif
#if !defined(PS5_NO_OVERRUN_GUARD)
    // SUBOPTIMALITY FIX (2026-07-02, DEFAULT-ON): after a same-level augment reshaped the forest, an edge can
    // become eligible NOW (slack already at target) without a PQ event -> advancing to the popped T would
    // OVERRUN it (Property-1.3 stop signal missed) -> negative slack -> suboptimal (the multi-augment residual
    // class, e.g. n=150 W5 s3; the duals are faithful -- proven -- so this missing-event overrun is the bug).
    // Only when we ALREADY augmented this phase AND time would advance (T>tNow) do the O(m) true-next-fire
    // scan (rare -- not per pop). If T overruns it, END the phase (relabel once) instead of overrunning.
    // O(m) only on post-augment time-advancing pops -> keeps multi-augment speed. Disable -DPS5_NO_OVERRUN_GUARD.
    if (mLazyActive && mNoRelabelAug && phaseAug && T > tNow) {
      LL bf = -1; LL tgt = (mCrit == CRIT1) ? 0 : -2;
      for (Sz En = 0; En < mNumNodes && (bf < 0 || bf > tNow); ++En) {
        if (mBParent[En] != NIL() || mLabel[En] != LBL_EVEN) continue;
        std::vector<Sz>& Ev = mNodeLeaves[En];
        for (Sz ui = 0; ui < Ev.size() && (bf < 0 || bf > tNow); ++ui) { Sz uu = Ev[ui]; Sz ufu = ufFind(uu);
          for (Sz kk = 0; kk < mAdjV[uu].size(); ++kk) { Sz ww = mAdjV[uu][kk];
            if (mIgnoreDummies && ww >= mNReal) continue;
            if (ufu == ufFind(ww)) continue;
            Sz Wn = outerNode(ww); if (Wn == En) continue; int lw = mLabel[Wn];
            if (lw != LBL_FREE && !treeLive(Wn)) lw = LBL_FREE;
            if (lw != LBL_EVEN && lw != LBL_FREE) continue;
            LL s = slackZW(uu, ww, 2 * mAdjW[uu][kk]);
            LL d = (lw == LBL_EVEN) ? (s - tgt) / 2 : (s - tgt);
            if (d < 0) continue;
            LL fire = tNow + d; if (bf < 0 || fire < bf) bf = fire; } } }
      if (bf >= 0 && T > bf) {
        // an event is due at bf < T -> end the phase (rebuild forest from the matching) instead of overrun.
        if (mLazyActive) lazyMaterialize();
        relabel();
        if (mLazyActive) lazyBegin();
        mUsed.assign(mN, 0); phaseAug = false; augTime = -1; phaseHadShrink = false; phaseAugs = 0; tNow = mTNow; mPQ.reset();
        reseedAllEven(tNow, 5);
#ifdef PS5_DISHEAP
        disHeapRebuild(tNow);
#endif
        continue;
      }
    }
#endif
    if (mLazyActive && T > tNow) lazySyncRoles();      // settle ONLY when time advances
    applyDualTo(T, &tNow);

    Sz u = (Sz)ev.a, w = (Sz)ev.b;
    if (mNoRelabelAug && (mUsed[u] || mUsed[w])) {          // consumed by an augment this phase
#if defined(PS5_PERSIST_FOREST) && defined(PS5_PERSIST_CLOCK)
      mReoffer.push_back(ev);   // still eligible -- re-offer after the boundary clears mUsed
#endif
#ifdef PS5_EVTRACE
      std::fprintf(stderr,"[TR] DISC-USED (%lu,%lu) T=%lld\n",(unsigned long)u,(unsigned long)w,(long long)T);
#endif
      continue;
    }
#ifdef PS5_VTX_BEST
    // entry (.,w) has left the queue. If it WAS w's covering best (recorded fire == T at
    // this epoch), clear the record and let the discard sites REPAIR w; a stale duplicate
    // carried no coverage -- skip both (the 52M-repairs/30% profile: repairing on EVERY
    // discard re-inserted millions of already-covered targets).
    bool vbWasBest = (w < mVBestEp.size() && mVBestEp[w] == mVBestEpoch && mVBestFire[w] == T);
    if (vbWasBest) mVBestEp[w] = -1;
#endif
    if (ufFind(u) == ufFind(w)) {
#ifdef PS5_EVTRACE
      std::fprintf(stderr,"[TR] DISC-UF (%lu,%lu) T=%lld\n",(unsigned long)u,(unsigned long)w,(long long)T);
#endif
#ifdef PS5_VTX_BEST
      if (vbWasBest) {
#ifdef PS5_VB_DEFER_REPAIR
        vbRepairDefer(w);
#else
        vbRepair(w, tNow);
#endif
      }
#endif
      continue;              // now same outer blossom
    }
    Sz E = outerNode(u), W = outerNode(w);
    if (mLabel[E] != LBL_EVEN || !treeLive(E)) {
#ifdef PS5_EVTRACE
      std::fprintf(stderr,"[TR] DISC-ELBL (%lu,%lu) T=%lld E=%lu lbl=%d live=%d\n",
        (unsigned long)u,(unsigned long)w,(long long)T,(unsigned long)E,mLabel[E],(int)treeLive(E));
#endif
#ifdef PS5_VTX_BEST
      if (vbWasBest) {
#ifdef PS5_VB_DEFER_REPAIR
        vbRepairDefer(w);
#else
        vbRepair(w, tNow);
#endif
      }
#endif
      continue;
    }
    if (E == W) {
#ifdef PS5_VTX_BEST
      if (vbWasBest) {
#ifdef PS5_VB_DEFER_REPAIR
        vbRepairDefer(w);
#else
        vbRepair(w, tNow);
#endif
      }
#endif
      continue;
    }
    // Eligibility at pop: a correctly-scheduled grow/formation event is due when slack reaches the
    // criterion target (CRIT1: 0; CRIT2/3: -2, unmatched edges). Stale events (slack not yet there)
    // are discarded. (Matched-edge eligibility for CRIT2 is handled by the grow following matched
    // edges, not by these scheduled events.)
    // Eligible when slack == target (CRIT1: 0, CRIT2/3: -2). The over-tight (<tgt) handling is gated:
    // mAllowOvertight lets the carried-dual (Property-2) regime treat slack<tgt as eligible-now; the
    // standalone full-domination engine keeps strict ==tgt (348 gate).
    { LL tgt = (mCrit == CRIT1) ? 0 : -2; LL sl = slackZ(u, w);
      if (mAllowOvertight ? (sl > tgt) : (sl != tgt)) {
#ifdef PS5_EVTRACE
        std::fprintf(stderr,"[TR] DISC-SLACK (%lu,%lu) T=%lld sl=%lld\n",(unsigned long)u,(unsigned long)w,(long long)T,(long long)sl);
#endif
#ifdef PS5_VTX_BEST
        if (vbWasBest) {
#ifdef PS5_VB_DEFER_REPAIR
        vbRepairDefer(w);
#else
        vbRepair(w, tNow);
#endif
      }
#endif
        continue;
      } }
#ifdef PS5_QSTAT
    ++mQAct;   // event survived all staleness/eligibility filters -> does real work
#endif

    int lblW = mLabel[W];
    if (lblW != LBL_FREE && !treeLive(W)) lblW = LBL_FREE;
#ifdef PS5_EVTRACE
    std::fprintf(stderr,"[TR] ACT (%lu,%lu) T=%lld lblW=%d tNow=%lld\n",
      (unsigned long)u,(unsigned long)w,(long long)T,lblW,(long long)tNow);
#endif

    if (lblW == LBL_EVEN) {
      if (lcaOuter(E, W) != NIL()) {
#ifdef PS5_SHRINK_AFTER_AUG_RELABEL
        // SCOPED FIX (trace-derived): the corruption is a SHRINK at the SAME dual time T after an augment
        // already happened this T -- the shrink reshapes the forest the augment reshaped -> stale roles ->
        // negative slack. Disjoint same-T AUGMENTS are safe (vertex-disjoint via mUsed); only a same-T SHRINK
        // after an augment is unsafe. So: if we already augmented at THIS T, END the phase (relabel) BEFORE
        // this shrink -- restart clean. Keeps multi-augment (many disjoint augments per T) but serializes
        // augment-then-shrink. Much narrower than pqT>=augTime (which relabels on every same-T event).
        if (mNoRelabelAug && phaseAug && augTime == tNow) {
          if (mLazyActive) lazyMaterialize();
          relabel();
          if (mLazyActive) lazyBegin();
          mUsed.assign(mN, 0); phaseAug = false; augTime = -1; phaseAugs = 0; tNow = mTNow; mPQ.reset();
          reseedAllEven(tNow, 6);
#ifdef PS5_DISHEAP
          disHeapRebuild(tNow);
#endif
          continue;
        }
#endif
        if (mLazyActive && mNoRelabelShrink) {
          // S5-C STEP 2: INCREMENTAL no-relabel shrink (the persistent-forest path).
          // Keep the PQ + the rest of the forest. Collect the formerly-ODD cycle nodes
          // (they become EVEN inside B; their edges were never scheduled) BEFORE shrink,
          // and settle every absorbed node's lazy role to FREE (so its z freezes -- it is
          // no longer a root). After shrink, B is EVEN with z0=0 fresh; setNodeRole(B,EVEN)
          // starts its z growing now; schedule ONLY the newlyEven nodes' edges. Formerly-
          // EVEN cycle nodes' events stay in the PQ and are discarded on pop (ufFind same).
          Sz Lb = lcaOuter(E, W);
          std::vector<Sz> newlyEven;
          for (Sz x = E; x != Lb; x = treeParent(x)) {
            if (mLabel[x] == LBL_ODD) newlyEven.push_back(x);
            if (x >= mN) setNodeRole(x, LBL_FREE);
          }
          for (Sz x = W; x != Lb; x = treeParent(x)) {
            if (mLabel[x] == LBL_ODD) newlyEven.push_back(x);
            if (x >= mN) setNodeRole(x, LBL_FREE);
          }
          if (Lb >= mN) setNodeRole(Lb, LBL_FREE);
          shrink(E, W, u, w);                     // sets mLabel[B]=EVEN, reparents cycle
          Sz B = outerNode(u);
#ifdef PS5_UF_CHECK
          ufCheck("shrink");
#endif
          mZ0[B] = 0; mZAcc[B] = 0; mNRole[B] = LBL_FREE; mNStamp[B] = mTNow;  // fresh z base
          setNodeRole(B, LBL_EVEN);               // B's z grows from now
#ifdef DIVDBG
          { Sz bi2=0; Sz mB2 = findCycEdgeMisalign(&bi2); static bool s_r=false;
            if (mB2 != NIL() && !s_r) { s_r=true; std::fprintf(stderr,"[MISALIGN-AFTER-NRSHRINK B=%lu idx=%lu k=%lu]\n",
              (unsigned long)mB2,(unsigned long)bi2,(unsigned long)mBCycle[mB2].size()); std::fflush(stderr); } }
#endif
#ifdef PS5_ZSHADOW
          if (mZShadowOk) { zShadowEnsure(B); mZShadow[B] = 0; }   // new blossom starts at z=0
#endif
          // NOTE (docs/F_Plan): a shrink-time lazySyncRoles() does NOT fix the shrink-drift (n=40 W5 s4
          // still subopt) -- only a full lazyMaterialize+lazyBegin rebase does, i.e. the error is already
          // in yAcc/y0 before this point, not a settle-timing issue at the shrink. Exact-cheap fix pending.
          for (Sz i = 0; i < newlyEven.size(); ++i) scheduleNodeEdges(newlyEven[i], tNow);
#ifdef PS5_VTX_BEST
          // v5 (replaces v4's O(boundary) whole-blossom rescan -- the residual super-linear term):
          // absorbed formerly-EVEN members' pre-existing entries stay EXACTLY valid across the
          // shrink (EVEN->EVEN rate unchanged; external edges carry no z(B) term), and formerly-
          // ODD members get the newlyEven rescan above. The ONLY coverage lost is the two
          // endpoints of the CONSUMED formation event (clear-on-pop wiped their records and their
          // runner-ups were suppressed). Repair exactly those two: O(deg), not O(boundary).
          #ifdef PS5_VB_DEFER_REPAIR
          vbRepairDefer(u); vbRepairDefer(w);
#else
          vbRepair(u, tNow); vbRepair(w, tNow);
#endif
#endif
#ifdef PS5_SHRINK_DUALFIX
          // RESULT (2026-07-02): refreshing the dual base after EVERY shrink+dissolution STILL fails n=100 W50
          // s3 + n=150 W5 s3 -> the drift is NOT "refresh too rarely"; the REFRESH ITSELF is lossy because
          // lazyMaterialize folds the ALREADY-DRIFTED zLazy into mZ. Proves the lazy z/y ARITHMETIC is wrong,
          // not the refresh cadence. Only the scan engine's from-scratch recompute is correct. Kept gated OFF.
          lazyMaterialize(); relabel(); lazyBegin(); tNow = mTNow; mPQ.reset();
          reseedAllEven(tNow, 7);
#ifdef PS5_DISHEAP
          disHeapRebuild(tNow);
#endif
          continue;
#endif
#ifdef PS5_RESCHED_BLOSSOM
          // ATTEMPT (INERT, gated OFF): reschedule ALL of B's boundary edges after shrink, in case B overgrows
          // its z past a boundary edge (trace n=100 W50 s3: blossom 117 z=2 -> edge 70-79 slack -10). Did NOT
          // fix the 2 residuals -- the z is already WRONG when B forms in the dense same-T dissolution->shrink
          // cascade, not a post-formation scheduling gap. Rescheduling can't repair a pre-existing bad z.
          scheduleNodeEdges(B, tNow);
#endif
#ifdef PS5_SHRINK_RESCHED_ALL
          // RESIDUAL FIX TEST (150 = dual overrun): after a no-relabel shrink, reschedule ALL EVEN nodes so
          // no boundary event is missing (the clamp diagnosis showed 150 overruns a boundary edge).
          reseedAllEven(tNow, 8);
#endif
#ifdef DIVDBG
          { Sz cy = findTreeCycle();
            if (cy != NIL()) { std::fprintf(stderr,"[CYCLE-AFTER-SHRINK B=%lu at node=%lu E=%lu W=%lu Lb=%lu]\n",
              (unsigned long)B,(unsigned long)cy,(unsigned long)E,(unsigned long)W,(unsigned long)Lb); std::fflush(stderr); return false; } }
#endif
#ifdef MATE_AUDIT
          { Sz fb; Sz a=mateAsym(&fb); static int r=0; if(a>0 && r++<3)
            std::fprintf(stderr,"[ASYM-AFTER-SHRINK B=%lu E=%lu W=%lu u=%lu w=%lu asym=%lu firstBad=%ld]\n",
              (unsigned long)B,(unsigned long)E,(unsigned long)W,(unsigned long)u,(unsigned long)w,(unsigned long)a,(long)fb); }
#endif
          phaseHadShrink = true;   // a shrink at this dual level -> a later same-T augment may thread the
          continue;                // reshaped forest -> serialize same-T augments (checked at the phase boundary).
        }
        // STEP 1 fallback (correct, relabels): materialize, shrink+relabel, re-lazyBegin.
        if (mLazyActive) lazyMaterialize();
        shrink(E, W, u, w);
#ifdef DIVDBG
        { Sz bi2=0; Sz mB2 = findCycEdgeMisalign(&bi2); static bool s_r=false;
            if (mB2 != NIL() && !s_r) { s_r=true; std::fprintf(stderr,"[MISALIGN-AFTER-RELSHRINK B=%lu idx=%lu k=%lu]\n",
              (unsigned long)mB2,(unsigned long)bi2,(unsigned long)mBCycle[mB2].size()); std::fflush(stderr); } }
#endif
        relabel();
        if (mLazyActive) lazyBegin();
        tNow = mTNow; mPQ.reset();   // KEEP capacity (was clear()+init(1) = realloc storm)
        reseedAllEven(tNow, 9);
#ifdef PS5_DISHEAP
        disHeapRebuild(tNow);   // post-augment relabel: rebuild the dissolution heap
#endif
        continue;
      } else {
        if (mLazyActive && mNoRelabelAug) {
#ifdef PS5_BLOSSOM_AUG_SERIALIZE
          // ATTEMPT (INERT on the residual, gated OFF): relabel before an augment through a blossom after a
          // same-level shrink. Did NOT fix n=100 W50 s3 / n=150 W5 s3 -- the residual is a DEEPER layer: the
          // blossom z-dual is already built WRONG by the same-T shrink CASCADE (7 shrinks build E=111), so
          // relabeling before the augment does not repair the pre-existing bad z. Needs a fix in the same-T
          // shrink z-arithmetic itself, upstream of the augment. Left off.
          if (phaseHadShrink && (E >= mN || W >= mN)) {
            if (mLazyActive) lazyMaterialize();
            relabel();
            if (mLazyActive) lazyBegin();
            mUsed.assign(mN, 0); phaseAug = false; augTime = -1; phaseHadShrink = false; phaseAugs = 0; tNow = mTNow; mPQ.reset();
            reseedAllEven(tNow, 10);
#ifdef PS5_DISHEAP
            disHeapRebuild(tNow);
#endif
            continue;
          }
#endif
#ifdef PS5_RELABEL_COMPLEX_AUG
          // SCOPED FIX (trace-derived): a 2nd+ same-T augment through SHARED structure (matched vertices or a
          // blossom) corrupts the forest -> stale roles -> negative slack. Disjoint FREE-FREE augments (both
          // endpoints unmatched, singleton nodes) are SAFE (the fast greedy phase). So relabel before a same-T
          // augment ONLY if it touches shared structure -> keeps the greedy phase fast, serializes complex augs.
          if (phaseAug && augTime == tNow &&
              (E >= mN || W >= mN || mMate[mBBase[E]] != NIL() || mMate[mBBase[W]] != NIL())) {
            if (mLazyActive) lazyMaterialize();
            relabel();
            if (mLazyActive) lazyBegin();
            mUsed.assign(mN, 0); phaseAug = false; augTime = -1; tNow = mTNow; mPQ.reset();
            reseedAllEven(tNow, 11);
#ifdef PS5_DISHEAP
            disHeapRebuild(tNow);
#endif
            continue;
          }
#endif
          // S5-C STEP 3: per-level MULTI-AUGMENT (MV/HK phase). Mark the whole augmenting
          // path USED (vertex-disjointness this level: no later augment THIS LEVEL reuses
          // a consumed vertex -> progress + correctness), augment, set the phase level,
          // and CONTINUE at the SAME dual level. The phase boundary above relabels ONCE
          // when time must advance. Re-routing through just-matched edges is handled
          // BETWEEN phases (the relabel rebuilds the forest); within a level disjoint
          // augments are correct. This is what gives O(sqrt n) phases/scale.
#ifdef PS5_RELABEL_AUG2
          // ATTEMPT (BROKEN, gated OFF): relabel before a 2nd+ same-level augment (phaseAug true). The
          // `continue` after relabel corrupts the nrAug path (gate 342/348, misses/suboptimal) -- the popped
          // edge/PQ state is inconsistent with the freshly-relabeled forest. Reverted. The correct-but-slow
          // PS5_RELABEL_AUG (relabel AFTER every augment) works; the challenge is doing it cheaply.
          if (phaseAug) {
            if (mLazyActive) lazyMaterialize();
            relabel();
            if (mLazyActive) lazyBegin();
            mUsed.assign(mN, 0); phaseAug = false; augTime = -1;
            continue;
          }
#endif
#ifdef PS5_PATHCLEAN
          // (INERT on the tested cases -- pathClean was already true; the bug is NOT interior-vertex reuse.
          // Kept off. The real cause is multiple same-level augments corrupting the carried duals, see below.)
          if (!pathClean(E) || !pathClean(W)) continue;
#endif
#ifdef PS5_SINGLE_AUG
          // DIAGNOSTIC: force ONE augment per phase (relabel after each) while staying on the event path,
          // to test whether MULTIPLE same-level augments are what corrupt the duals -> suboptimal (LEMON).
          markPathUsed(E); markPathUsed(W);
#ifdef PS5_EVTRACE
          std::fprintf(stderr,"[TR] AUG2 (%lu,%lu) E=%lu W=%lu tNow=%lld\n",
            (unsigned long)u,(unsigned long)w,(unsigned long)E,(unsigned long)W,(long long)tNow);
#endif
          if (mLazyActive) lazyMaterialize();
          augment(E, W, u, w);
          augmentedAny = true; return true;
#endif
          markPathUsed(E); markPathUsed(W);
#ifdef PS5_EVTRACE
          std::fprintf(stderr,"[TR] AUG2 (%lu,%lu) E=%lu W=%lu tNow=%lld\n",
            (unsigned long)u,(unsigned long)w,(unsigned long)E,(unsigned long)W,(long long)tNow);
#endif
#ifdef PS5_AUGTRACE
          // trace same-level augment ordering: print E,W,u,w, augTime, and whether this is a 2nd+ aug at the
          // SAME dual level (phaseAug already true) -- the suspect for the stale-forest corruption.
          std::fprintf(stderr,"[AUG tNow=%lld augTime=%lld %s SAME-T?%d pathCleanE=%d pathCleanW=%d E=%lu W=%lu (u=%lu,w=%lu) mate[u]=%ld mate[w]=%ld]\n",
            (long long)tNow, (long long)augTime, phaseAug?"2nd+":"first", (int)(phaseAug && augTime==tNow),
            (int)pathClean(E),(int)pathClean(W),
            (unsigned long)E,(unsigned long)W, (unsigned long)u,(unsigned long)w,(long)mMate[u],(long)mMate[w]);
#endif
#if defined(PS5_FREE_CONSUMED)
          // (INERT: settling the consumed tree's roles to FREE after augment does NOT fix the suboptimality
          // -- the corruption is in the carried lazy dual BASE (yAcc/y0), not the consumed roles. Kept off.
          // STAGE 2 does the equivalent BATCHED at the phase boundary: retireDeadTrees() -- per-augment
          // freeConsumedTree was 84% of the run, an O(n)-scan-per-call quadratic.)
          Sz rootE = mRoot[E], rootW = mRoot[W];
#endif
#ifdef PS5_SETTLE_TWO
          Sz rE = mRoot[E], rW = mRoot[W];   // capture consumed-tree roots BEFORE augment clears them
#endif
#ifdef PS5_SETTLE_PATH
          // SURGICAL: capture EXACTLY the augmenting-path nodes (E..rootE, W..rootW via treeParent) BEFORE
          // augment, then settle ONLY those nodes' leaf roles to FREE after. O(path), not O(all-trees).
          // Tests whether the stale roles are confined to the PATH (cheap fix) vs whole coupled forest.
          std::vector<Sz> pathNodes;
          { for (Sz x = E; ; ) { pathNodes.push_back(x); Sz p = treeParent(x); if (p==NIL()) break; x = p; }
            for (Sz x = W; ; ) { pathNodes.push_back(x); Sz p = treeParent(x); if (p==NIL()) break; x = p; } }
#endif
          augment(E, W, u, w);
#ifdef DIVDBG
          { Sz bi2=0; Sz mB2 = findCycEdgeMisalign(&bi2); static bool s_r=false;
            if (mB2 != NIL() && !s_r) { s_r=true; std::fprintf(stderr,"[MISALIGN-AFTER-AUG2 B=%lu idx=%lu k=%lu E=%lu W=%lu]\n",
              (unsigned long)mB2,(unsigned long)bi2,(unsigned long)mBCycle[mB2].size(),(unsigned long)E,(unsigned long)W); std::fflush(stderr); } }
          { Sz cy = findTreeCycle();
            if (cy != NIL()) { std::fprintf(stderr,"[CYCLE-AFTER-AUG2 at node=%lu E=%lu W=%lu u=%lu w=%lu tNow=%lld]\n",
              (unsigned long)cy,(unsigned long)E,(unsigned long)W,(unsigned long)u,(unsigned long)w,(long long)tNow); std::fflush(stderr); return false; } }
#endif
#ifdef PS5_SETTLE_PATH
          for (Sz pi = 0; pi < pathNodes.size(); ++pi) { Sz nd = pathNodes[pi];
            if (nd >= mNodeLeaves.size()) continue;
            std::vector<Sz>& lv = mNodeLeaves[nd];
            for (Sz j = 0; j < lv.size(); ++j) setVtxRole(lv[j], LBL_FREE);
            if (nd >= mN && mBParent[nd]==NIL()) setNodeRole(nd, LBL_FREE);
            if (mBParent[nd]==NIL()) { mLabel[nd]=LBL_FREE; mParent[nd]=NIL(); mParEU[nd]=NIL(); mParEV[nd]=NIL(); mRoot[nd]=NIL(); }
          }
#endif
#ifdef PS5_SLACKMON
          // after each same-level augment, scan ALL edges for the FIRST negative-slack (infeasible dual) and
          // report which augment created it -> pinpoints the corrupting operation. finalSlack = y+y+z-2w.
          { static bool firstNeg=false; if(!firstNeg){ for(Sz uu=0; uu<mN && !firstNeg; ++uu)
              for(Sz kk=0; kk<mAdjV[uu].size(); ++kk){ Sz ww=mAdjV[uu][kk]; if(ww<uu) continue;
                if(((ScalingS5Engine*)this)->ufFind(uu)==((ScalingS5Engine*)this)->ufFind(ww)) continue; // skip INTERNAL blossom edges (slack>=0 only required for EXTERNAL edges)
                LL sl=((ScalingS5Engine*)this)->finalSlack(uu,ww);
                if(sl<0){ ScalingS5Engine* self=(ScalingS5Engine*)this;
                  Sz ou=self->outerNode(uu), ow=self->outerNode(ww);
                  std::fprintf(stderr,"[FIRST-NEG-SLACK after AUG(E=%lu,W=%lu,u=%lu,w=%lu): edge %lu-%lu slack=%lld w=%lld | y[%lu]=%lld y[%lu]=%lld | outer[%lu]=%lu(z=%lld) outer[%lu]=%lu(z=%lld) sameBloss=%d]\n",
                  (unsigned long)E,(unsigned long)W,(unsigned long)u,(unsigned long)w,(unsigned long)uu,(unsigned long)ww,(long long)sl,(long long)mAdjW[uu][kk],
                  (unsigned long)uu,(long long)self->finalY(uu),(unsigned long)ww,(long long)self->finalY(ww),
                  (unsigned long)uu,(unsigned long)ou,(long long)(ou>=mN?self->zLazy(ou):0),
                  (unsigned long)ww,(unsigned long)ow,(long long)(ow>=mN?self->zLazy(ow):0),(int)(self->ufFind(uu)==self->ufFind(ww))); firstNeg=true; break; } } } }
#endif
#if defined(PS5_FREE_CONSUMED)
          freeConsumedTree(rootE); freeConsumedTree(rootW);
#endif
#ifdef PS5_REBASE_AUG
          // DIAGNOSTIC: rebase the lazy duals after each same-level augment (materialize -> lazyBegin), WITHOUT
          // the full O(m) edge reschedule. Tests whether the cheap dual rebase alone fixes the suboptimality
          // (isolating it from the expensive relabel/reschedule that makes PS5_SINGLE_AUG slow).
          if (mLazyActive) { lazyMaterialize(); lazyBegin(); }
#endif
#ifdef PS5_SETTLE_AUG
          // The corruption is stale EVEN/ODD ROLES on consumed nodes -> wrong y at materialize. Settling ALL
          // in-tree nodes FREE after the augment fixes correctness (relabel does the same, and does NOT touch
          // y) -- but tearing down EVERY tree forces a full regrow (edgeScan ~13x, 20x slow). SETTLE_TWO only
          // settles the 2 CONSUMED trees (rootE,rootW), leaving the other same-level trees intact -> cheap.
#ifdef PS5_SETTLE_TWO
          { Sz roots[2] = { rE, rW };
            for (int ri = 0; ri < 2; ++ri) { Sz root = roots[ri]; if (root == NIL()) continue;
              for (Sz i = 0; i < mNumNodes; ++i) {
                if (mBParent[i] != NIL()) continue;
                if (mRoot[i] != root) continue;
                std::vector<Sz>& lv = mNodeLeaves[i];
                for (Sz j = 0; j < lv.size(); ++j) setVtxRole(lv[j], LBL_FREE);
                if (i >= mN) setNodeRole(i, LBL_FREE);
                mLabel[i] = LBL_FREE; mParent[i] = NIL(); mParEU[i] = NIL(); mParEV[i] = NIL(); mRoot[i] = NIL();
              } } }
#else
          for (Sz i = 0; i < mNumNodes; ++i) {
            if (mBParent[i] != NIL()) continue;
            if (mLabel[i] == LBL_FREE) continue;
            std::vector<Sz>& lv = mNodeLeaves[i];
            for (Sz j = 0; j < lv.size(); ++j) setVtxRole(lv[j], LBL_FREE);
            if (i >= mN) setNodeRole(i, LBL_FREE);
            mLabel[i] = LBL_FREE; mParent[i] = NIL(); mParEU[i] = NIL(); mParEV[i] = NIL(); mRoot[i] = NIL();
          }
          phaseAug = false; augTime = -1;
#endif
#endif
#ifdef PS5_RELABEL_AUG
          // DIAGNOSTIC: structural RELABEL after each same-level augment (rebuild the forest from the matching)
          // but KEEP the PQ (no O(m) reschedule). Tests whether the STRUCTURE rebuild is the fix (vs duals).
          if (mLazyActive) lazyMaterialize();
          relabel();
          if (mLazyActive) lazyBegin();
          mUsed.assign(mN, 0); phaseAug = false; augTime = -1;
#endif
#ifdef MATE_AUDIT
          { Sz fb; Sz a=mateAsym(&fb); static int r=0; if(a>0 && r++<3)
            std::fprintf(stderr,"[ASYM-AFTER-EE-AUG E=%lu W=%lu u=%lu w=%lu asym=%lu firstBad=%ld]\n",
              (unsigned long)E,(unsigned long)W,(unsigned long)u,(unsigned long)w,(unsigned long)a,(long)fb); }
#endif
#ifdef S5DBG
          ++gS5Augs;
#endif
#ifdef PS5_SERIALIZE_2ND_AUG
          // (INERT for perf, OFF: relabel after a 2nd+ same-T complex augment fixes n=150 + gate but is
          // ~10x slow -- per-augment relabel+reseed is O(augments*m); the complexity gate is nearly always
          // true post-greedy. The phase-boundary phaseHadShrink approach is the fast one; extend THAT.)
          if (phaseAug && augTime == tNow &&
              (E >= mN || W >= mN || (mBBase[E] < mN && mMate[mBBase[E]] != NIL())
                                  || (mBBase[W] < mN && mMate[mBBase[W]] != NIL()))) {
            if (mLazyActive) lazyMaterialize();
            relabel();
            if (mLazyActive) lazyBegin();
            mUsed.assign(mN, 0); phaseAug = false; augTime = -1; phaseHadShrink = false; tNow = mTNow; mPQ.reset();
            reseedAllEven(tNow, 12);
#ifdef PS5_DISHEAP
            disHeapRebuild(tNow);
#endif
            augmentedAny = true;
            continue;
          }
#endif
#ifdef PS5_COMPLEX_AUG_FLAG
          // (INERT for perf, OFF: flagging phaseHadShrink when a 2nd+ same-T augment threads matched/blossom
          // structure catches n=150 but is still ~10x slow -- the condition is nearly always true post-greedy,
          // so it serializes almost every same-T event. The harmful-vs-safe distinction needs PATH analysis,
          // not a cheap structural predicate. Kept OFF.)
          if (phaseAug && augTime == tNow &&
              (E >= mN || W >= mN || (mBBase[E] < mN && mMate[mBBase[E]] != NIL())
                                  || (mBBase[W] < mN && mMate[mBBase[W]] != NIL())))
            phaseHadShrink = true;
#endif
          augmentedAny = true; phaseAug = true; augTime = tNow; ++phaseAugs;
          continue;
        }
        // HALT after the first augment (relabel path): the outer solve loop relabels.
        if (mLazyActive) { lazyMaterialize(); mLazyActive = false; }   // mY/mZ current for next search
        augment(E, W, u, w);
        return true;
      }
    } else if (lblW == LBL_FREE) {
      Sz bw = mBBase[W]; Sz m = mMate[bw];
#if !defined(PS5_NO_STALE_MATE)
      // M==E ROOT-CAUSE FIX (2026-07-02, DEFAULT-ON): W's base bw may carry a STALE mate m -- a prior
      // rotateToBase flip re-matched m internally to its own blossom, ORPHANING the bw->m half-edge
      // (mMate[m] != bw). Then bw is TRULY free; treating m as a live mate drives the bogus grow-cycle
      // (M==E) + the O(nm) scan restart = the DISTINCT-weight cliff (n=25600 W=1000: 14.4s->2.5s, restarts
      // 1->0, LEMON-optimal). Detect the asymmetry and treat bw as free (m==NIL path = a clean augment).
      // Traced n=54 seed=30 W=5: bw=49, m=37, but mMate[37]=13 != 49 -> stale. CONSUMER-side fix (the two
      // rotateToBase-side repairs both failed: transfer hangs, clear corrupts). 348/348 brute gate + LEMON
      // optimal to n=25600. Disable with -DPS5_NO_STALE_MATE.
      if (m != NIL() && (m >= mN || mMate[m] != bw)) m = NIL();
#endif
      if (m == NIL()) {
        if (mLazyActive && mNoRelabelAug) {
          markPathUsed(E);                       // W is a free node (its own path end)
          mUsed[w] = 1;
#ifdef PS5_EVTRACE
          std::fprintf(stderr,"[TR] AUGF (%lu,%lu) E=%lu W=%lu tNow=%lld\n",
            (unsigned long)u,(unsigned long)w,(unsigned long)E,(unsigned long)W,(long long)tNow);
#endif
          augment(E, W, u, w);
#ifdef MEE_DUMP
          { Sz asym=0,fb=NIL(); for(Sz v=0;v<mN;++v){Sz mv=mMate[v]; if(mv!=NIL()&&(mv>=mN||mMate[mv]!=v)){++asym;if(fb==NIL())fb=v;}}
            static int rep=0; if(asym>0 && rep++<5) std::fprintf(stderr,"[AUG-CORRUPT after augment(E=%lu,W=%lu,u=%lu,w=%lu): asym=%lu firstBad=%ld]\n",
              (unsigned long)E,(unsigned long)W,(unsigned long)u,(unsigned long)w,(unsigned long)asym,(long)fb); }
#endif
#ifdef S5DBG
          ++gS5Augs;
#endif
#ifdef PS5_SERIALIZE_2ND_AUG
          // (INERT for perf, OFF: relabel after a 2nd+ same-T complex augment fixes n=150 + gate but is
          // ~10x slow -- per-augment relabel+reseed is O(augments*m); the complexity gate is nearly always
          // true post-greedy. The phase-boundary phaseHadShrink approach is the fast one; extend THAT.)
          if (phaseAug && augTime == tNow &&
              (E >= mN || W >= mN || (mBBase[E] < mN && mMate[mBBase[E]] != NIL())
                                  || (mBBase[W] < mN && mMate[mBBase[W]] != NIL()))) {
            if (mLazyActive) lazyMaterialize();
            relabel();
            if (mLazyActive) lazyBegin();
            mUsed.assign(mN, 0); phaseAug = false; augTime = -1; phaseHadShrink = false; tNow = mTNow; mPQ.reset();
            reseedAllEven(tNow, 13);
#ifdef PS5_DISHEAP
            disHeapRebuild(tNow);
#endif
            augmentedAny = true;
            continue;
          }
#endif
#ifdef PS5_COMPLEX_AUG_FLAG
          // (INERT for perf, OFF: flagging phaseHadShrink when a 2nd+ same-T augment threads matched/blossom
          // structure catches n=150 but is still ~10x slow -- the condition is nearly always true post-greedy,
          // so it serializes almost every same-T event. The harmful-vs-safe distinction needs PATH analysis,
          // not a cheap structural predicate. Kept OFF.)
          if (phaseAug && augTime == tNow &&
              (E >= mN || W >= mN || (mBBase[E] < mN && mMate[mBBase[E]] != NIL())
                                  || (mBBase[W] < mN && mMate[mBBase[W]] != NIL())))
            phaseHadShrink = true;
#endif
          augmentedAny = true; phaseAug = true; augTime = tNow; ++phaseAugs;
          continue;
        }
        if (mLazyActive) { lazyMaterialize(); mLazyActive = false; }   // mY/mZ current for next search
        augment(E, W, u, w); return true;
      } else {
        Sz M = outerNode(m);
        // GROW-CYCLE CASE (M==E): the free node W's matched partner m resolves to the SAME node E we grow
        // FROM. Growing would set mParent[M=E]=W, overwriting E's parent -> a treeParent 2-cycle (the
        // n=200/W=200 divergence). Neither an in-place relabel-and-continue NOR a mid-way handoff resolves
        // it in the event-PQ order (both re-hit the same M==E deferral / lose optimality). The correct fix
        // = re-solve this instance with searchOnce, the SAME engine's exhaustive SCAN-order search (NOT a
        // foreign engine; same duals/blossoms/Property-1 invariants), whose post-relabel full re-scan
        // takes the augmenting path the PQ order keeps deferring. Rare (~5% of distinct-weight graphs, ~1
        // trigger). Correct + optimal (verified n=200: wt=16779==GMG).
        if (M == E) {
#ifdef MEE_COUNT
          { extern long long gMEEtriggers; ++gMEEtriggers; }
#endif
#ifdef MEE_DUMP
          { static int dumped=0; if(dumped++<3) {
            // matching-symmetry audit: count v with mate[mate[v]] != v (corruption).
            Sz asym=0; Sz firstA=NIL();
            for (Sz v=0; v<mN; ++v) { Sz mv=mMate[v]; if(mv!=NIL() && (mv>=mN || mMate[mv]!=v)) { ++asym; if(firstA==NIL()) firstA=v; } }
            std::fprintf(stderr,"[MEE-DUMP u=%lu w=%lu m=%lu | mate[u]=%ld mate[w]=%ld u==m?%d | ASYMMETRIC-mates=%lu firstBad=%ld]\n",
              (unsigned long)u,(unsigned long)w,(unsigned long)m,(long)mMate[u],(long)mMate[w],(int)(u==m),
              (unsigned long)asym,(long)firstA); std::fflush(stderr); } }
#endif
#ifdef PS5_MEE_SKIP
          // EXPERIMENT: with the negative-z spin fixed, does simply SKIPPING the M==E edge (the scan
          // engine treats this grow as a transient no-op cleaned by relabel) now stay correct? The event
          // stays on the fast path -- no scan restart. Discards this edge; the augment is found via others.
          continue;
#elif defined(PS5_MEE_RELABEL)
          // EXPERIMENT (re-try now that negative-z is fixed): in-EVENT relabel-and-continue. Materialize,
          // relabel (re-derive forest from the matching, dissolving the transient M==E), reseed the PQ, and
          // continue the SAME event loop. Stays on the fast path (no scan restart). Earlier this MISSED the
          // match -- but that was before the negative-z spin fix; re-evaluate.
          if (mLazyActive) lazyMaterialize();
          relabel();
          if (mLazyActive) lazyBegin();
          tNow = mTNow; mPQ.reset();
#if defined(S5DBG) || defined(PS5_BIGNODE_CACHE)
          ++mDbgEpoch;
#endif
          reseedAllEven(tNow, 14);
#ifdef PS5_DISHEAP
          disHeapRebuild(tNow);
#endif
          phaseAug = false; augTime = -1;
          continue;
#else
          mEventStuck = true; if (mLazyActive) { lazyMaterialize(); mLazyActive = false; }
          return augmentedAny;
#endif
        }
#ifdef PS5_STAGE0
        ++mS0Grows;
        if (W < mS0Parent.size() && M < mS0Parent.size() &&
            mS0Parent[W] == E && mS0Lbl[W] == LBL_ODD &&
            mS0Parent[M] == W && mS0Lbl[M] == LBL_EVEN)
          ++mS0Redundant;
#endif
#ifdef PS5_EVTRACE
        std::fprintf(stderr,"[TR] GROW (%lu,%lu): W=%lu->ODD M=%lu->EVEN (m=%lu bw=%lu) tNow=%lld\n",
          (unsigned long)u,(unsigned long)w,(unsigned long)W,(unsigned long)M,(unsigned long)m,(unsigned long)bw,(long long)tNow);
#endif
#ifdef DIVDBG
        // PRE-GROW cycle forensics: if W or M is already on E's ancestor chain, this grow will
        // create a parent cycle -- dump the chain with per-node label/root/liveness to expose
        // WHICH bookkeeping (mRoot/treeLive vs mParent) is stale.
        { bool bad = false; Sz gsteps = 0;
          for (Sz x = E; x != NIL() && gsteps <= mNumNodes + 5; x = treeParent(x), ++gsteps)
            if (x == W || x == M) { bad = true; break; }
          if (bad) {
            std::fprintf(stderr,"[GROW-WOULD-CYCLE E=%lu W=%lu M=%lu tNow=%lld | E-chain:",
              (unsigned long)E,(unsigned long)W,(unsigned long)M,(long long)tNow);
            Sz st2 = 0;
            for (Sz x = E; x != NIL() && st2 <= mNumNodes + 5; x = treeParent(x), ++st2)
              std::fprintf(stderr," %lu(lbl%d,root%ld,live%d)",(unsigned long)x,mLabel[x],
                (long)(x < mRoot.size() ? (long)mRoot[x] : -9),(int)treeLive(x));
            std::fprintf(stderr," | W: lbl%d root%ld live%d | M: lbl%d root%ld live%d",
              mLabel[W],(long)mRoot[W],(int)treeLive(W),mLabel[M],(long)mRoot[M],(int)treeLive(M));
            { Sz bw2 = mBBase[W], m2 = (bw2 < mN ? mMate[bw2] : NIL());
              Sz bM = mBBase[M];
              std::fprintf(stderr," | bases: bw=%lu m=%ld mMate[m]=%ld | mBBase[M]=%lu mMate[base(M)]=%ld m-leafOf-M=%d bParent(M)=%ld]\n",
                (unsigned long)bw2,(long)m2,(long)(m2!=NIL()&&m2<mN?(long)mMate[m2]:-9),
                (unsigned long)bM,(long)(bM<mN?(long)mMate[bM]:-9),
                (int)(m2!=NIL()&&M>=mN&&M<mNodeLeaves.size()?[&]{for(Sz j=0;j<mNodeLeaves[M].size();++j)if(mNodeLeaves[M][j]==m2)return 1;return 0;}():-1),(long)mBParent[M]); }
            std::fflush(stderr); return false; } }
#endif
        MG_CHK("grow-event", m, bw);
        mLabel[W] = LBL_ODD; mParent[W] = E; mParEU[W] = w; mParEV[W] = u; mRoot[W] = mRoot[E];
        mLabel[M] = LBL_EVEN; mParent[M] = W; mParEU[M] = m; mParEV[M] = bw; mRoot[M] = mRoot[E];
#ifdef DIVDBG
        parw(W,'G'); parw(M,'G');
#endif
        markLabelDirty(W); markLabelDirty(M);   // event grow: only W,M change label -> dirty their leaves
#ifdef DIVDBG
        { Sz cy = findTreeCycle();
          if (cy != NIL()) { std::fprintf(stderr,"[CYCLE-AFTER-GROW at node=%lu E=%lu W=%lu M=%lu tNow=%lld]\n",
            (unsigned long)cy,(unsigned long)E,(unsigned long)W,(unsigned long)M,(long long)tNow); std::fflush(stderr); return false; } }
#endif
#ifdef PS5_GROW_SERIALIZE
        // TEST (INERT for perf): serialize after a GROW-following-augment too. FIXES n=150 W5 s3 (the last
        // residual) + 348/348 gate, BUT ~14x SLOWER on distinct (28s@25600) -- grows are frequent, so this
        // serializes nearly every augment, killing multi-augment batching (same trap as PHASE_SAMET). The
        // grow predicate is too broad; the harmful case is narrower (a same-T augment that THREADS the grown
        // structure). Kept OFF; the correct narrow trigger is the next target.
        if (phaseAug) phaseHadShrink = true;
#endif
#ifdef PS5_DISHEAP
        if (W >= mN && mBParent[W] == NIL()) disHeapSchedule(W, tNow);   // W is now an ODD root
#endif
        scheduleNodeEdges(M, tNow);   // new EVEN node's edges (dissolution handled by scan)
      }
    }
  }
  if (mLazyActive) { lazyMaterialize(); mLazyActive = false; }   // keep mY/mZ current
  // no-relabel-aug: we continued past augments, so report whether ANY happened (the outer
  // solve loop needs this -- else it sees "false" and wrongly concludes no perfect matching).
  return augmentedAny;
}

inline ScalingS5Engine::LL
ScalingS5Engine::solve(std::vector<Sz>* mateOut) {
  if (mN == 0) { if (mateOut) mateOut->clear(); return 0; }
  if (mN % 2 != 0) return -1;
  Sz cap = 2 * mN + 1;
  bool savedUseEvent = mUseEvent, savedNRS = mNoRelabelShrink, savedNRA = mNoRelabelAug,
       savedNRD = mNoRelabelDissolve, savedLazyE = mLazyEvent;
 restart:                                   // M==E re-enters here in PURE SCAN mode (identical to a fresh
  mNumNodes = mN;                           // scan solve -> LEMON-optimal, unlike the in-place reset).
#ifdef PS5_ANCESTOR_HINT
  ++mAhTok; ahEnsure();   // fresh solve: ids reused -> all hints invalid
#endif
  mMate.assign(mN, NIL());
  mY.assign(mN, mMaxW);
  mZ.assign(cap, 0);
  mBParent.assign(cap, NIL());
  mBBase.assign(cap, NIL());
  mBCycle.assign(cap, std::vector<Sz>());
  mBCycEU.assign(cap, std::vector<Sz>());
  mBCycEV.assign(cap, std::vector<Sz>());
  mNodeLeaves.assign(cap, std::vector<Sz>());
  mLabel.assign(cap, LBL_FREE);
  mParent.assign(cap, NIL());
  mParEU.assign(cap, NIL());
  mParEV.assign(cap, NIL());
  mRoot.assign(cap, NIL());
  mUf.assign(mN, 0);
  mFreeBlossoms.clear();
  for (Sz v = 0; v < mN; ++v) { mBBase[v] = v; mNodeLeaves[v].assign(1, v); mUf[v] = v; }
#ifdef PS5_LEAF_LIST
    lfRebuildFromVectors();
#endif
  if (mLazyVerify) {
    // size lazy arrays to cap (covers blossoms created later: a not-yet-existing node
    // has z0=0/Acc=0/role=FREE, which is exactly its state at creation) + seed base.
    mTNow = 0;
    mY0.assign(mN, mMaxW); mYAcc.assign(mN, 0); mVRole.assign(mN, LBL_FREE); mVStamp.assign(mN, 0);
    mZ0.assign(cap, 0); mZAcc.assign(cap, 0); mNRole.assign(cap, LBL_FREE); mNStamp.assign(cap, 0);
  }

  mEventStuck = false;
  for (Sz guard = 0; guard <= 2 * mN + 2; ++guard) {
    bool perfect = true;
    for (Sz v = 0; v < mN; ++v) if (mMate[v] == NIL()) { perfect = false; break; }
    if (perfect) break;
#ifdef DIVDBG
    { Sz nf = 0; for (Sz v = 0; v < mN; ++v) if (mMate[v] == NIL()) ++nf;
      std::fprintf(stderr,"[SOLVE iter=%lu free=%lu useEvent=%d]\n",(unsigned long)guard,(unsigned long)nf,(int)mUseEvent); }
#endif
    bool ok = mUseEvent ? searchEvent() : searchOnce();
    if (mEventStuck) {
#ifdef S5DBG
      { extern long long gS5ScanRestarts; ++gS5ScanRestarts; }
#endif
      // M==E: re-solve in pure scan mode (goto restart). NOTE (docs/F_Plan): the scan engine searchOnce is
      // ORDER-DEPENDENT and does NOT guarantee max-weight -- it can land on a suboptimal perfect matching
      // (LEMON: 162825 vs 162840 on n=400 seed=2, edge-order-sensitive). This is a CORE-ENGINE correctness
      // gap, separate from M==E; the scan-restart inherits it. Documented; fix pending.
      mEventStuck = false;
      mUseEvent = false; mNoRelabelShrink = false; mNoRelabelAug = false;
      mNoRelabelDissolve = false; mLazyEvent = false; mLazyActive = false;
      goto restart;
    }
    if (!ok) { mUseEvent=savedUseEvent; mNoRelabelShrink=savedNRS; mNoRelabelAug=savedNRA;
               mNoRelabelDissolve=savedNRD; mLazyEvent=savedLazyE; return -1; }
  }
  for (Sz v = 0; v < mN; ++v) if (mMate[v] == NIL()) { mUseEvent=savedUseEvent; mNoRelabelShrink=savedNRS;
    mNoRelabelAug=savedNRA; mNoRelabelDissolve=savedNRD; mLazyEvent=savedLazyE; return -1; }
#ifdef PS5_DUAL_CERT
  // OPTIMALITY CERTIFICATE (safety net for the pole-fix / no-relabel dissolution): a max-weight PERFECT
  // matching REQUIRES dual feasibility -- every edge slack y[u]+y[v]-2w >= 0. The no-relabel path can leave
  // INFEASIBLE duals on rare configs (pole-fix uncovered nested cases: ~6/224, e.g. n=60 W100 s2 edge 21-30
  // slack=-22) -> a suboptimal (but perfect) matching. Detect it in O(m): if ANY edge has slack < 0, the
  // event result is not certified optimal -> re-solve via the pure SCAN engine (always optimal). This makes
  // the fast pole-fix path DEFAULT-SAFE: correct always, fast on the ~97% it certifies, scan-fallback only
  // on the rare uncertified ones. (Runs only if the event engine was used.)
  if (mUseEvent) {
    bool feasible = true;
    for (Sz u = 0; u < mN && feasible; ++u)
      for (Sz k = 0; k < mAdjV[u].size(); ++k) {
        Sz v = mAdjV[u][k]; if (v < u) continue;
        if (mY[u] + mY[v] - 2 * mAdjW[u][k] < 0) { feasible = false; break; }
      }
    if (!feasible) {
      mUseEvent = false; mNoRelabelShrink = false; mNoRelabelAug = false;
      mNoRelabelDissolve = false; mLazyEvent = false; mLazyActive = false;
      goto restart;   // uncertified -> pure-scan re-solve
    }
  }
#endif
  mUseEvent=savedUseEvent; mNoRelabelShrink=savedNRS; mNoRelabelAug=savedNRA;
  mNoRelabelDissolve=savedNRD; mLazyEvent=savedLazyE;   // restore (engine reusable)

  LL total = 0;
  if (mateOut) mateOut->assign(mN, NIL());
  for (Sz v = 0; v < mN; ++v) {
    Sz m = mMate[v];
    if (m != NIL() && v < m) {
      LL w2 = w2max(v, m);
      if (w2 < 0) return -1;                  // matched on a non-edge
      total += w2 / 2;
      if (mateOut) { (*mateOut)[v] = m; (*mateOut)[m] = v; }
    }
  }
  return total;
}

} // namespace Matchbox

#endif // SCALING_S5_ENGINE_H
