#ifndef MATCHING_ENGINE_H
#define MATCHING_ENGINE_H

#define EPS 0.0000001

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm> // for fill, find, sort.
#include <functional> // for std::greater.
#include <numeric> // for accumulate.
#include <cassert>
#include <cmath>  
#include <limits>
#include <queue>
#include <stack>
#include <map>
#include <set>
#include <omp.h>
#include <ctime>
#include <cstdlib>
#ifdef MBX_SHARED_D3_TIE_ORDER
#include <stdexcept>
#endif

namespace Matchbox {
inline long double profNow() {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return (long double)ts.tv_sec + (long double)ts.tv_nsec * 1e-9L;
}
}
// Hot-path fine-grained timers (per-event/per-edge profNow/clock_gettime syscalls).
// Compiling with -DMBX_NO_HOTPROF removes them so we can measure how much of the
// runtime is the profiling instrumentation itself vs real algorithm work, AND so the
// optimized/benchmark binaries don't pay the clock-syscall tax. profNow() is ~20-50ns
// each; with several firing per event-loop iteration (and per pushD3/make_edge) this is
// seconds of overhead at the 100M+ iterations big Gabow runs reach. NOTE: ALL hot-loop
// timers must go through MBX_HOTNOW() (not bare profNow()) or -DMBX_NO_HOTPROF won't
// fully remove them; the per-region PROFLOOP breakdown reads zero under the flag, by
// design (use a normal build for the diagnostic breakdown).
#ifdef MBX_NO_HOTPROF
#  define MBX_HOTNOW() (0.0L)
#else
#  define MBX_HOTNOW() profNow()
#endif

#include "Functional.h" // for ValGreater.
#include "Types.h"
#include "Graph.h"
#include "Utility.h" // for ResizeVector.
#include "VecItmQue.h"
#include "LstItmQue.h"


#include "VecItmStk.h"
#include "LstItmStk.h"
#include "LstItmIdxdQue.h"
#include "LstItmAndValIdxdMinPriQue.h"
#include "BhpItmAndValIdxdMinPriQue.h"
#include "LstItmAndValIdxdMaxPriQue.h"
#include "BhpItmAndValIdxdMaxPriQue.h"
#include "Timer.h"
#include "ConcatenableQueue.h"
#include "DynamicNCA.h"
#include "BlossomMerging.h"
#include "QueuedRowPopulation.h"
#include "FrontierMembership.h"
#include "LargestReuseAudit.h"
#include "TwinCleanupAudit.h"
#include "MonotoneSupportAudit.h"
#include "dominant_membership_merge.h"
#include "MembershipReuseAudit.h"
#if defined(MBX_GABOW_FRONTIER_AUDIT) || defined(MBX_GABOW_FRONTIER_VERIFY)
#include "GabowFrontierAudit.h"
#endif
#ifdef MBX_CLASSIC_EVENT_TRACE
#include "classic_event_trace.h"
#endif

namespace Matchbox {

class FastMatchingEngineV2 {
public:
  // Backend selects the D3 (shrink/augment) data structure:
  //   kGMG1986  - Galil-Micali-Gabow 1986: global heap + per-edge stale-pop.
  //   kGabow2018 - Gabow 2018: supporting forest with dynamic NCA, back-edge
  //                reduction in pushD3, BlossomMerging packet structure for
  //                the blossom-merging problem (Section 3.1 of the paper).
  //                F2: forest only (placeholder mPQD3 still in use).
  // Default is kGMG1986 so existing callers see identical behavior.
  enum Backend { kGMG1986 = 0, kGabow2018 = 1 };

  FastMatchingEngineV2() : mBackend(kGMG1986) {}
  explicit FastMatchingEngineV2(Backend b) : mBackend(b) {}

  void compute(const Graph& graph,
               std::vector<Size>* mateVec,
               Size* card,
               Val* edgWght);
  bool verifyCertificate(const std::vector<Size>& mateVec,
                         Size* card,
                         Val* edgWght,
                         Val* dualObj,
                         Val* minSlack) const;

private:
#ifdef MBX_GABOW_FRONTIER_TEST_ACCESS
  friend struct FrontierReuseTest;
#endif
#ifdef MBX_CLASSIC_EVENT_TRACE
  ClassicEventTrace::Sink mEventTrace;
  std::vector<ClassicEventTrace::Signature> mEventMembers;
  Val traceEffectiveY(Size v, Size b) const {
    // The accepted candidate already supplied b. Do not call topLevel here:
    // even its const overload writes caches and diagnostic counters.
    if (b==cNullItm || b>=mBlsm.size()) return mVtx[v].pot;
    if (mBlsm[b].status==BS_EVEN) return mVtx[v].pot-(mDeltaSum-mBlsm[b].offset);
    if (mBlsm[b].status==BS_ODD) return mVtx[v].pot+(mDeltaSum-mBlsm[b].offset);
    return mVtx[v].pot;
  }
  void traceAction(ClassicEventTrace::Action action, Size u, Size v,
                   Size blossom=cNullItm, Size aux1=cNullItm, Size aux2=cNullItm,
                   Size tree1=cNullItm, Size tree2=cNullItm, int variant=0) {
    if (!mEventTrace.enabled() || mAbortRun) return;
    const bool have = blossom < mBlsm.size() && blossom < mEventMembers.size();
    mEventTrace.commit(action,u,v,have?mBlsm[blossom].base:cNullItm,
      have?mEventMembers[blossom]:ClassicEventTrace::Signature(),aux1,aux2,tree1,tree2,variant);
  }
#endif
  /*
   * Pure GMG-style Edmonds engine. Correctness comes from the persistent
   * primal-dual event loop below; the slow Edmonds oracle, output repair,
   * global heap compaction, and profiling scaffolding are deliberately absent.
   */
  enum BlossomStatus {
    BS_EVEN    = 0,
    BS_MATCHED = 1,
    BS_ODD     = 2,
    BS_UNMATCHED = 3
  };

  typedef ::FrontierMembership<Size> Membership;
  typedef ::FrontierMembershipView<Size> MembershipView;
  struct BlossomData {
    BlossomStatus status;
    Size parent;
    Size base;
    Size tree;
    Size pred;
    Size next;
    Val pot;
    Val offset;
    std::vector<Size> children;
    std::vector<Size> childEdges;
    std::vector<ConcatenableQueue::Handle> splitLeaves;
    Membership vertices;
    // vertex -> sub-blossom index, for blossomContainsVertex / childIndexContaining.
    // Was std::map<Size,Size>, but it is rebuilt from scratch on EVERY shrink and the
    // per-node red-black-tree malloc/free dominated the shrink path (~4.4s of 5.1s,
    // ~25% allocator churn on blossom-heavy graphs). Replaced with a SORTED vector of
    // (vtx,childIdx) pairs: build = push_back + one sort; lookup = lower_bound (O(log));
    // clear() keeps capacity so rebuilds reuse the buffer (no per-shrink malloc storm).
    // Safe because there is no per-key erase/update anywhere -- only bulk build, clear,
    // and find. Sorted on the vtx key.
    std::vector<std::pair<Size, Size> > childOfVtx;
    ConcatenableQueue::Handle cqRoot;
    ConcatenableQueue::Handle cqLastLeaf;
    Size materializedBase;
    Size d1MinVtx;
    bool d1MinValid;
    Val d2Key;
    Size d2Vtx;
    Size d2Edge;
    Size d2Tree;
    bool d2Queued;
    Size generation;
    Size topCache;
    Size topCacheGeneration;
    Size listedTree;
    bool alive;
  };

  struct VertexData {
    Val pot;
    ConcatenableQueue::Handle cqLeaf;
    std::map<Size, std::pair<Val, Size> > bestByTree;
    Val d1Key;
    Size d1Version;
    bool d1Queued;
    Val d2Key;
    Size d2Edge;
    Size d2Blsm;
    Size d2Tree;
    Size d2Version;
    bool d2Queued;
    Size lastEvenScanTree;
    Size lastEvenScanTop;
  };

  struct HeapEntryD1 {
    Val key;
    Size vtx;
    Size version;
    bool operator>(const HeapEntryD1& o) const { return key > o.key; }
  };

  struct HeapEntryD2 {
    Val key;
    Size matchedVtx;
    Size edgeIdxInMatched;
    Size matchedBlsm;
    Size tree;
    Size version;
    bool operator>(const HeapEntryD2& o) const { return key > o.key; }
  };

  struct D2Less {
    bool operator()(const HeapEntryD2& a, const HeapEntryD2& b) const
    {
      if (a.key < b.key) return true;
      if (b.key < a.key) return false;
      if (a.matchedVtx < b.matchedVtx) return true;
      if (b.matchedVtx < a.matchedVtx) return false;
      if (a.edgeIdxInMatched < b.edgeIdxInMatched) return true;
      if (b.edgeIdxInMatched < a.edgeIdxInMatched) return false;
      if (a.matchedBlsm < b.matchedBlsm) return true;
      if (b.matchedBlsm < a.matchedBlsm) return false;
      return a.tree < b.tree;
    }
  };

  struct HeapEntryD3 {
    Val key;
    Size u;
    Size edgeIdx;
#ifdef MBX_GMG_D3_OWNERSHIP_EXPERIMENT
    unsigned long long publication;
#endif
    bool operator>(const HeapEntryD3& o) const {
#ifdef MBX_SHARED_D3_TIE_ORDER
      // pushD3Impl publishes canonical u<v in sorted normalized rows. Thus
      // (u,edgeIdx) has exactly the order of the checked BM original-edge rank.
      // Costs are finite and compared exactly; EPS remains an acceptance rule.
      if (key > o.key) return true;
      if (o.key > key) return false;
      if (u != o.u) return u > o.u;
      return edgeIdx > o.edgeIdx;
#else
      return key > o.key;
#endif
    }
  };

  struct HeapEntryD4 {
    Val key;
    Size blsm;
    bool operator>(const HeapEntryD4& o) const { return key > o.key; }
  };

  typedef std::priority_queue<HeapEntryD1, std::vector<HeapEntryD1>,
                              std::greater<HeapEntryD1> > PQD1;
#ifndef MBX_D2_LAZY
  // DEFAULT: RB-tree PQD2 with exact-entry erase. A lazy-deletion binary-heap variant was
  // A/B-tested 2026-07-06 (MBX_D2_LAZY below): NEUTRAL on distinct weights (the gprof 12.5%
  // equal_range share was instrumentation inflation, not real -O2 cost) and it HANGS on
  // uniform W=1 (>600s where this version takes 2s at n=12800): the physical erase maintains
  // the one-live-entry-per-blossom invariant that the D2 discard/re-push branches rely on
  // under heavy blossom churn. Keep the set.
  class PQD2 {
  public:
    bool empty() const { return mSet.empty(); }
    Size size() const { return mSet.size(); }
    const HeapEntryD2& top() const { return *mSet.begin(); }
    void pop()
    {
      if (!mSet.empty()) mSet.erase(mSet.begin());
    }
    void push(const HeapEntryD2& e) { mSet.insert(e); }
    void erase(const HeapEntryD2& e) { mSet.erase(e); }
    void clearStorage()
    {
      mSet.clear();
    }
  private:
    std::set<HeapEntryD2, D2Less> mSet;
  };
#else
  // EXPERIMENT (MBX_D2_LAZY, do NOT default): binary heap + lazy deletion. Correct + neutral on
  // distinct weights (exact vs LEMON to n=102400) but HANGS on uniform W=1 (see note above).
  class PQD2 : public std::priority_queue<HeapEntryD2, std::vector<HeapEntryD2>,
                                          std::greater<HeapEntryD2> > {
  public:
    void erase(const HeapEntryD2&) {}   // lazy: pop-validation discards stale entries
    void clearStorage()
    {
      std::vector<HeapEntryD2> empty;
      this->c.swap(empty);
    }
  };
#endif
  class PQD3 : public std::priority_queue<HeapEntryD3,
                                          std::vector<HeapEntryD3>,
                                          std::greater<HeapEntryD3> > {
  public:
    void clearStorage()
    {
      std::vector<HeapEntryD3> empty;
      this->c.swap(empty);
    }
  };
  typedef std::priority_queue<HeapEntryD4, std::vector<HeapEntryD4>,
                              std::greater<HeapEntryD4> > PQD4;

  // CSR-backed replacement for a jagged std::vector<std::vector<T> > that is
  // indexed [u][i] with i in [0, rowLen(u)). Stores all rows in ONE flat
  // std::vector + a row-offset array, eliminating the per-row heap allocation
  // (and its header/rounding/fragmentation tax) that dominated Gabow's memory.
  // Exposes just the surface the engine uses: arr.size(), arr.assign(rowLens,
  // fill), arr[u] -> Row, Row.size(), Row[i] (read/write), Row.assign(n, fill).
  // C++98-compatible. The row lengths must match the corresponding mAdj rows;
  // callers build them from mAdj before assign().
  template <class T>
  class FlatJagged {
  public:
    FlatJagged() {}
    // Build offsets from per-row lengths and fill all cells with `fill`.
    void assign(const std::vector<Size>& rowLens, const T& fill) {
      mOff.assign(rowLens.size() + 1, 0);
      for (Size r = 0; r < rowLens.size(); ++r)
        mOff[r + 1] = mOff[r] + rowLens[r];
      mData.assign(mOff.empty() ? 0 : mOff[mOff.size() - 1], fill);
    }
    void clear() { mData.clear(); mOff.clear(); }
    Size size() const { return mOff.empty() ? 0 : mOff.size() - 1; }
    bool empty() const { return size() == 0; }
    class Row {
    public:
      Row(T* base, Size len) : mBase(base), mLen(len) {}
      Size size() const { return mLen; }
      T& operator[](Size i) { return mBase[i]; }
      const T& operator[](Size i) const { return mBase[i]; }
      // Re-fill this row in place (lengths are fixed by the CSR layout, so the
      // requested n must equal the row's length; we assert that in debug).
      void assign(Size n, const T& fill) {
        (void)n;
        for (Size i = 0; i < mLen; ++i) mBase[i] = fill;
      }
    private:
      T* mBase;
      Size mLen;
    };
    Row operator[](Size u) {
      return Row(mData.empty() ? (T*)0 : &mData[mOff[u]], mOff[u + 1] - mOff[u]);
    }
    // Read-only row view, for access from const methods.
    class ConstRow {
    public:
      ConstRow(const T* base, Size len) : mBase(base), mLen(len) {}
      Size size() const { return mLen; }
      const T& operator[](Size i) const { return mBase[i]; }
    private:
      const T* mBase;
      Size mLen;
    };
    ConstRow operator[](Size u) const {
      return ConstRow(mData.empty() ? (const T*)0 : &mData[mOff[u]],
                      mOff[u + 1] - mOff[u]);
    }
  private:
    std::vector<T> mData;
    std::vector<Size> mOff;  // size N+1; row u is [mOff[u], mOff[u+1])
  };

  Size mNumVtxs;
  Size mUnmatched;
  Val  mDeltaSum;
  Size mBlossomCap;
  bool mAbortRun;

  // CSR-flattened (FlatJagged): adjacency + weights, fixed length per row after the
  // init build. Collapses 2N tiny per-vertex heap allocations into 2 flat arrays.
  // Hot-loop access mAdj[u][i] / mAdj[u].size() works via FlatJagged::Row / ConstRow.
  FlatJagged<Size> mAdj;
  FlatJagged<Val>  mWgt;
  // NOTE: mOrigWgt was a verbatim copy of mWgt used to recover pre-tie-break
  // weights. With cardTieEps hardwired to 0.0, mWgt is never perturbed and is
  // never mutated by the algorithm (it reweights via dual potentials, not mWgt),
  // so mWgt IS the original weight. The redundant ~8 B/edge copy was removed;
  // the former mOrigWgt[u][i] reads now use mWgt[u][i]. If a nonzero cardTieEps
  // is ever reintroduced, restore a separate pre-eps array for those 3 readers.
  // CSR-flattened (FlatJagged): per-edge, fixed length per row after init. Collapses
  // N tiny per-vertex heap allocations into one flat array each (the vector<vector>
  // allocator tax was the dominant per-edge memory cost). Access mRevIdx[u][i] /
  // mRevIdx[u].size() works via FlatJagged::Row.
  FlatJagged<Size> mRevIdx;
  FlatJagged<char> mD3Queued;
#ifdef MBX_GMG_D3_OWNERSHIP_EXPERIMENT
  // Only the current publication may consume or clear the canonical slot.
  // Superseded physical records remain lazy, but cannot become owners again.
  struct GmgD3Publication {
    Val key;
    unsigned long long version;
    GmgD3Publication() : key(0), version(0) {}
  };
  FlatJagged<GmgD3Publication> mGmgD3Publication;
#endif
  // Profiling counters (Gabow path only).
  Size mProfCleanupCount;
  long double mProfCleanupTime;
  Size mProfCleanupEdgesRemoved;
  Size mProfResetCount;
  long double mProfResetTime;
  Size mProfPushD3Count;
  Size mProfPushD3Early;
  Size mProfMakeEdgeCount;
  Size mProfMergeCount;
  Size mProfDisposeCount;
  Size mProfMinEdgeCount;
  Size mProfStalePopCount;
  Size mProfNcaCount;
  long double mProfPushD3Time;          // total pushD3 wallclock (incl make_edge)
  long double mProfMinEdgeTime;
  long double mProfMergeTime;
  long double mProfDisposeTime;
  long double mProfMakeEdgeTime;        // make_edge call sites in pushD3
  long double mProfNcaTime;             // mAuxNCA.nca() call sites
  long double mProfStaleCheckTime;      // consumer Bv/status check + remove
  Size mProfStaleCostRepush;            // pops hitting stale-cost re-push path
  Size mProfShrinks;                    // blossom shrinks (both backends)
  Size mProfAugments;                   // augmentations (both backends)
  Size mProfEvictCalls;                 // matched-blossom D3 retirement calls
  long double mProfEvictTime;           // time in eager eviction
  bool mAnyShrinkYet;                   // true once any blossom shrink occurred
  long double mProfWallStart;           // wall-clock anchor for periodic dumps
  long double mProfDestroyTreeTime;     // destroyTree calls in stale branch
  long double mProfD2CandTime;          // addMatchedVertexD2Candidate in stale branch
  // Breakdown of the previously-"unaccounted" event-loop time (both backends).
  long double mProfD124ScanTime;        // scanning mPQD1/D2/D4 heap-tops per iter
  long double mProfActExtendTime;       // extendOnArc (grow) dispatch
  long double mProfActShrinkTime;       // shrinkOnEdge (blossom) dispatch
  long double mProfShrinkHookTime;      // Gabow-only mBM.merge+discard+repush block inside shrink
  unsigned long long mProfShrinkDiscarded; // total mMergeDiscarded re-pushes across shrinks
  long double mProfShrinkRefreshTime;   // refreshEvenBlossomEdgesSubset inside shrink (both backends)
  long double mProfShrinkVtxBuildTime;  // B.vertices + childOfVtx(map) build loop (both backends)
  long double mProfShrinkCqMergeTime;   // mCQ.merge inside shrink (both backends)
  unsigned long long mProfEvictSlotsScanned;  // adj slots scanned in evict (O(degree) waste signal)
  unsigned long long mProfEvictLiveRemovals;  // actual live-handle removes in evict
  long double mProfActAugmentTime;      // augmentOnArc/Edge dispatch
  long double mProfActUnmatchTime;      // unmatchNode (d1) dispatch
  unsigned long long mProfLoopIters;    // total event-loop iterations
  Size mProfPushD3CallsAll;             // pushD3 calls, both backends (volume compare)
  // Stale-pop categorization (which reason dominates pdb1HYS churn).
  Size mProfStaleSameBlsm;
  Size mProfStaleOdd;
  Size mProfStaleMatched;
  Size mProfStaleInactiveTree;
  Size mProfStaleOther;
  // Re-push multiplicity: per-edge make_edge count + aggregates. High average
  // = re-push churn = the suspected dense-graph (pdb1HYS) inefficiency.
  std::vector<std::vector<Size> > mPushMult;       // [u][edgeIdx] -> push count
  unsigned long long mProfTotalEdgePushes;         // sum over all edges
  unsigned long long mProfDistinctEdgesPushed;     // edges pushed >= 1 time
  unsigned long long mProfMaxPushMult;             // max single-edge push count

#ifdef MBX_GABOW_SEARCH_AUDIT
  // ==== Per-augmentation-interval audit (compile-gated; absent normally) ====
  // The engine is a CONTINUOUS multi-tree loop: there is no per-search object
  // lifetime. The operational reporting boundary is the AUGMENTATION EVENT:
  // each successful augmentOnArc/augmentOnEdge ends one interval (destroying
  // the consumed tree(s); all other trees persist logically). audEmitSearch()
  // snapshots every cumulative counter there and reports the delta since the
  // previous augmentation, so each IAUD line is one augmentation interval.
  // GABOW_SEARCH_AUDIT_EMIT=1 is retained as a compatibility environment name;
  // it prints per-interval lines. IAUDSUM and run-level peaks print at the end.
  enum { kAudN = 47 };
  unsigned long long mAudAddLeafFirst;    // first-time installs of an original vertex
  unsigned long long mAudAddLeafReentry;  // re-entry installs (vertex was in T before)
  unsigned long long mAudRetires;         // retireGabowForestVertex (occurrence++)
  unsigned long long mAudNcaAncestor;     // nca(u,v) returned u or v (zero-span side)
  unsigned long long mAudSingleRec;       // reduction with one nondegenerate BM record
  unsigned long long mAudDoubleRec;       // reduction with two nondegenerate BM records
  unsigned long long mAudR8Reduced;       // original edges reduced through NCA/root
  unsigned long long mAudR8OneSide;       // one nondegenerate side (ancestor case)
  unsigned long long mAudR8TwoSide;       // two nondegenerate path sides
  unsigned long long mAudR8DummyRoot;     // different-tree reductions through rho
  unsigned long long mAudR8Ancestor;      // same-tree reductions with a zero-length side
  unsigned long long mAudSplits;          // splitBlossom (dissolve) events
  unsigned long long mAudTreeAdds;        // addTreeMember calls (outer-occurrence proxy)
  unsigned long long mAudMergeUpSteps;    // shrink-hook merge steps where depth INCREASED
  unsigned long long mAudMergeDownSteps;  // ... where depth decreased or held
  unsigned long long mAudMergeChains;     // shrink-hook merge chains (one per shrink w/ >=1 merge)
  unsigned long long mAudMergeChainMax;   // longest single merge chain
  unsigned long long mAudMergeOrderLogged;// examples printed (capped)
  unsigned long long mAudSearches;        // AUGMENTATION INTERVALS emitted (not paper searches)
  std::vector<unsigned char> mAudSeen;    // vertex ever installed in T (first-vs-reentry)
  unsigned long long mAudPrev[kAudN];
  unsigned long long mAudSum[kAudN];
  unsigned long long mAudMax[kAudN];
  bool mAudEmitPerSearch;
  // Interval-local edge multiplicity (separate from the run-global mPushMult):
  // epoch-stamped per-edge counts; O(1)/push, allocated once (audit builds only).
  std::vector<std::vector<unsigned> > mAudIvStamp;   // last interval that pushed [u][i]
  std::vector<std::vector<unsigned> > mAudIvCount;   // pushes of [u][i] in that interval
  unsigned long long mAudIvDupPushes;     // submissions beyond the first in one interval
  unsigned long long mAudIvMaxMult;       // max submissions of one edge in one interval
  // ---- Stage-1 paired-chassis counters (BOTH backends; same gating) ----
  mutable unsigned long long mAudTopLevelCalls;
  unsigned long long mAudEvD1, mAudEvD2, mAudEvD3, mAudEvD4;   // consumed events
  unsigned long long mAudGrows;                        // extendOnArc dispatches
  unsigned long long mAudD2CandAdds;                   // addMatchedVertexD2Candidate
  unsigned long long mAudD2BestPushes;                 // pushD2BestForBlossom
  unsigned long long mAudD2Drops;                      // dropD2TreeCandidate
  unsigned long long mAudGmgD3Pushes;                  // kGMG1986 physical mPQD3 pushes
  unsigned long long mAudGmgD3Stale;                   // kGMG1986 stale D3 pops discarded
  static const char* audName(int i);
  void audCollect(unsigned long long* out) const;
  void audNoteAddLeaf(Size w);
  void audNoteIntervalPush(Size u, Size edgeIdx);
  void audEmitSearch(char kind);
  void audPrintSummary() const;
  void audPrintChassis() const;           // CHASSIS line, printed for BOTH backends
#endif // MBX_GABOW_SEARCH_AUDIT

  mutable std::vector<BlossomData> mBlsm;
  std::vector<VertexData>  mVtx;
  mutable std::vector<Size> mTopBlsmOfVtx;
  mutable std::vector<Size> mTopCacheOfVtx;
  mutable std::vector<Size> mTopCacheGenerationOfVtx;
  std::vector<Size>        mMate;
  std::vector<Size>        mTreeParent;
  std::vector<std::vector<Size> > mTreeMembers;
  // Diagnostic accounting for retained capacities across all historical tree
  // roots. clear() preserves these buffers, so their sum can greatly exceed
  // the O(n) number of currently live tree members.
  unsigned long long       mProfTreeMemberCapacity;
  Size                     mProfTreeMemberMaxRowCapacity;
  std::vector<Size>        mBestD2Stamp;
  std::vector<Size>        mBestD2Edge;
  std::vector<Size>        mBestD2Touched;
  std::vector<Val>         mBestD2Key;
  Size                     mBestD2CurrentStamp;
  std::vector<Size>        mShrinkStamp;
  std::vector<Size>        mShrinkLeftIndex;
  Size                     mShrinkCurrentStamp;
  Size                     mMaterializeDepth;

  PQD1 mPQD1;
  PQD2 mPQD2;
  PQD3 mPQD3;
  PQD4 mPQD4;

  std::vector<Size> mFreeBlsmIds;
  ConcatenableQueue mCQ;

  // ---- Gabow 2018 (kGabow2018) backend members ----
  // mBackend selects which D3 backend is live in the current compute() run.
  // mAuxNCA is the supporting forest T per Section 3 / Definition 3.1.
  // Maintained by add_leaf only (no cuts). Initial state: every vertex is
  // a leaf attached to a dummy root (mAuxRoot).
  // mInT[v] tracks whether v has been added to T (avoids double-addLeaf).
  // mBM is the BlossomMerging packet structure (Section 3.1). Allocated but
  // not yet consulted in F2 (placeholder mPQD3 still drives D3 events).
  Backend mBackend;
  DynamicNCA mAuxNCA;
  std::vector<char> mInT;
  Size mAuxRoot;
  BlossomMerging mBM;
  // Stride for packing (u, edgeIdx) into a Size payload for mBM:
  // payload = u * mAdjStride + edgeIdx. Set to max adjacency size at init.
  Size mAdjStride;
  // For each pushed-into-mBM edge (u, edgeIdx), the two back-edge handles
  // returned by mBM.make_edge. mBmHandleU[u][edgeIdx] is for back-edge
  // (u, a); mBmHandleV[u][edgeIdx] is for (v, a) where v = mAdj[u][edgeIdx].
  // cNullItm if not currently pushed.
  FlatJagged<BlossomMerging::EdgeHandle> mBmHandleU;
  FlatJagged<BlossomMerging::EdgeHandle> mBmHandleV;
  // Generations captured at push time, paired with the handles above, so a
  // later remove_edge_if can detect a freed/reused slot and no-op safely.
  FlatJagged<Size> mBmGenU;
  FlatJagged<Size> mBmGenV;
  // Mapping from matchbox blossom ID to mBM blossom ID. For atomic vertex
  // blossoms (blsmId < mNumVtxs and free EVEN), the mapping points at the
  // mBM singleton blossom. For super-blossoms and non-free vertices, this
  // is cNullItm (not tracked in mBM). Sized to mBlossomCap.
  std::vector<BlossomMerging::BlossomId> mMatchboxToBm;
  // GABOW_LOOPDIAG=1 enables the time-based event-loop progress dump.
  bool mGabowLoopDiag;
  // GABOW_DEEPT=1 (R_Scope prototype): faithful supporting-forest membership.
  // Vertices LEAVE T when they leave a live tree (whole-tree retirement via the
  // destroyTree conversions + dissolve-off-path leaves), so the existing grow
  // hooks re-install them at REAL depths. Default off pending the decision
  // gate; off = the permanent-star approximation (disclosed deviation).
  bool mGabowDeepForest;
  bool mGabowDeepSplitClear;   // bisect sub-flag (GABOW_DEEPT=1 on, =2 off)
  bool mGabowD3Log;            // GABOW_D3LOG=1: trace D3 pushes/fires (forensics)
  // Each vertex's supporting-forest occurrence is bumped when it leaves T.
  // Two 32-bit endpoint occurrences are packed into each live BM edge record,
  // avoiding an O(m) side table while making deep-forest lazy records safe.
  // GABOW_FORESTCHECK=1 enables reporting of the repair counters.
  bool mGabowForestCheck;
  std::vector<unsigned int> mForestGeneration;
  // Structural validation: same-tree NCAs that resolved to a non-root vertex.
  // Under the permanent star this is always 0; a faithful forest must show >0
  // on grow-heavy inputs.
  unsigned long long mProfNcaNonRoot;
  // Deep-forest install-coverage gaps repaired defensively at pushD3 (root
  // install, star semantics for that edge). Must be ~0; >0 flags a missing
  // install hook.
  unsigned long long mProfForestFallback;
  unsigned long long mProfForestStaleQueued;
  unsigned long long mProfForestStaleFired;
  unsigned long long mProfForestStaleKept;   // Stage-3: occurrence-stale records KEPT (price current)
  unsigned long long mProfForestStaleRebuilt;
  // Compile-gated MBX_GABOW_NOREPUSH_EXPERIMENT + GABOW_NOREPUSH=1
  // (breaks correctness): skip merge-discard re-publication to measure only an
  // upper-bound performance ceiling.  Unreachable in production builds.
  bool mGabowNoRepushExpt;
#ifdef MBX_GABOW_EDGE_WATCH
  Size mWatchEdgeU, mWatchEdgeV;
  Val mWatchLastSlack;
  bool mWatchEdgeEnabled, mWatchHaveSlack;
  void watchEdgeState(const char* where, Size eventCount, int eventKind);
#endif
  // Scratch buffer reused inside splitBlossom's Gabow hook to stash each
  // sub-blossom's atomic vertex list across the dispose_blossom + make_blossom
  // round (dispose clears mVtxBlossom mid-loop, so we capture before disposing).
  std::vector<std::vector<Size> > mSplitChildVtxsScratch;
  // Scratch buffer for payloads of edges discarded by I[B] dedup during
  // mBM.merge. After each merge, the engine iterates this and clears
  // mD3Queued so refresh can re-push the affected matchbox edges.
  std::vector<Size> mMergeDiscarded;

  void init(const Graph& graph);
  void run();
  void extractMatching(std::vector<Size>* mateVec, Size* card, Val* edgWght) const;
  void extractBlossomMatching(Size blossom, Size base, Size matching,
                              std::vector<Size>* mateVec) const;
  void setExtractMate(Size u, Size v, std::vector<Size>* mateVec) const;
  bool scoreMateVector(const std::vector<Size>& mateVec, Size* card,
                       Val* edgWght) const;
  void collectPositiveDualAncestors(
    Size vtx,
    const std::vector<Val>& zEff,
    std::vector<std::pair<Size, Val> >* ancestors) const;
  Val commonPositiveDual(Size u, Size v,
                         const std::vector<Val>& zEff) const;
  bool activeTree(Size tree) const;
  Size topLevel(Size vtx) const;
  void start();

  Val effY(Size v) const;
  Val effZ(Size b) const;
  Val edgeSlack(Size u, Size edgeIdx) const;
  Val d2Key(Size matchedVtx, Size edgeIdxInMatched) const;
  Val d3Key(Size u, Size edgeIdx) const;

  void matchedToEven(Size blossom, Size tree);
  void matchedToOdd(Size blossom);
  void evenToMatched(Size blossom, Size tree);
  void oddToMatched(Size blossom);
  void oddToEven(Size blossom, Size tree);
  void bakeEvenToMatchedForShrink(Size blossom);
  void bakeOddToMatchedForShrink(Size blossom);

  void alternatePath(Size even, Size tree);
  void destroyTree(Size tree);
  void unmatchNode(Size vtx);
  void extendOnArc(Size u, Size v, Size edgeIdxInU);
  void augmentOnArc(Size u, Size v, Size edgeIdxInU);
  void augmentOnEdge(Size u, Size v, Size edgeIdxInU);
  void shrinkOnEdge(Size u, Size v, Size edgeIdxInU, Size tree);
  void splitBlossom(Size blossom);
  bool restoreGabowFrontierMembership(Size parent,const std::vector<Size>& children,
      const std::vector<ConcatenableQueue::Handle>& roots,Membership* savedParent);
  void verifyLargestBmPartition(const std::vector<Size>& children,Size oldBm,bool after) const;
  void verifyGabowFrontierStable(bool force=false) const;

  void setGabowBlossomActive(Size blossom, bool active);
  void retireGabowForestVertex(Size v);
  void retireGabowD3ForMatchedBlossom(Size blossom);
  Size gabowForestOccurrenceTag(Size u, Size v) const;

  // R7: non-mutating walk of the EVEN-length alternating arc of an odd/outer
  // blossom's child cycle from the entry child to the base child. Emits the
  // ordered (child, entryAtom) sequence (same index walk as
  // augmentThroughBlossom, without flipping anything). Returns false on a
  // malformed cycle (caller falls back to the flat install and counts it).
  bool collectPathChildren(Size blossom, Size entryVtx,
                           std::vector<std::pair<Size,Size> >& out) const;
  // R7 grow-geometry install: chain the path children's entry atoms
  // (prev <- entry_1 <- entry_2 <- ... <- base side), attach each child's
  // remaining atoms under its entry atom, and return the LAST chained atom
  // (the base-side attachment point for the mate blossom). Falls back to the
  // flat-under-anchor install when the cycle is malformed.
  typedef std::vector<std::pair<Size,Size> > GabowDeferredTwins;
  void stageGabowTwins(Size u,Size i,GabowDeferredTwins& deferred);
  void removeStagedGabowTwins(const GabowDeferredTwins& deferred);
  std::vector<Size> mSupportHead,mSupportNext,mSupportTreeOf,mSupportSingleton;
  std::vector<unsigned int> mSupportOwnerStamp;
  unsigned int mSupportOwnerEpoch;
  void supportPath(Size b,Size x,std::vector<Size>& out);
  void supportAdd(Size v,Size parent,bool singleton,int kind);
  Size supportInstallPath(Size b,Size entry,Size anchor);
  void supportInner(Size b);
  void supportOuter(Size b);
  Size supportMergePath(Size v,Size ancestor);
  void supportBlossom(Size b,Size u,Size v);
  void supportRepublish();
  void supportRetireTree(Size tree);
  void supportVerifyStable(const char* where) const;
  Size installBlossomChainInForest(Size blossom, Size entryVtx, Size anchor);
  // Add Bv's atomic vertices to mAuxNCA under u and register Bv as an alive
  // mBM blossom (so back-edges targeting these vertices are valid).
  void addBvAtomicsToGabowForest(Size Bv, Size u, Size v);

  bool trivial(Size blsm) const { return blsm < mNumVtxs; }
  Size acquireBlsmId();
  void releaseBlsmId(Size id);
  void abortRun(const char* reason);

  template<class F>
  void forEachVertex(Size b, F& f) const;

  Val vertexBestSlack(Size v) const;
  void updateVertexCQLeaf(Size v);
  void clearVertexD2(Size v);
  void eraseD2Entry(Size v);
  void eraseD2BlossomEntry(Size blossom);
  void addMatchedVertexD2Candidate(Size matchedVtx, Size edgeIdxInMatched);
  void refreshMatchedVertexD2(Size v);
  void refreshMatchedVertexD2Impl(Size v, bool pushRoot);
  void pushD1(Size v);
  void pushD1ForBlossom(Size blossom);
  void pushD2Best(Size matchedVtx);
  void pushD2BestForBlossom(Size blossom);
  void dropD2TreeCandidate(Size matchedVtx, Size tree);
  void pushD3(Size u, Size edgeIdx);
  void pushD3Impl(Size u, Size edgeIdx);
  void pushD4(Size blossom);
  void addTreeMember(Size blossom, Size tree);
  void clearMate(Size vtx);
  void setMatePair(Size u, Size v);
  void flipPathFromEven(Size evenVtx);
  void flipPathFromTopEven(Size evenBlossom, Size entryVtx);
  void augmentThroughBlossom(Size blossom, Size entryVtx);
  void materializeBlossomForCurrentBoundary(Size blossom);
  void materializeBlossomMatching(Size blossom, Size baseVtx);
  void rotateBlossomToChild(Size blossom, Size childPos);
  void collectTreePath(Size blossom, std::vector<Size>* path) const;
  Size orientedTreeEdge(Size fromBlossom, Size toBlossom) const;
  Size reverseEncodedEdge(Size edge) const;
  Size edgeSrc(Size edge) const;
  Size edgeDst(Size edge) const;
  Size childIndexContaining(Size blossom, Size vtx) const;
  bool blossomContainsVertex(Size blossom, Size vtx) const;
  void collectBlossomVertices(Size blossom, std::vector<Size>* vtxs) const;
  void setTopLevelForBlossom(Size blossom, Size top);
  void clearBlossomD2(Size blossom);
  void refreshEvenBlossomEdges(Size blossom);
  void refreshEvenBlossomEdgesSubset(Size blossom,
                                     MembershipView vtxs);
  void refreshMatchedBlossomD2(Size blossom);
  void invalidateMaterializedAncestors(Size vtx);
  Val originalEdgeWeight(Size u, Size v) const;
  Size findEdgeIdx(Size u, Size v) const;
  Size reverseEdgeIdx(Size u, Size edgeIdx) const;
  static Size encodeEdge(Size u, Size v, Size numVtxs)
  { return u * numVtxs + v; }
  static Val infKey()
  { return std::numeric_limits<Val>::max() / 16.0; }
  // Append storage: preserve offsets of every original engine member.
  QueuedRowPopulation<Size> mD3QueuedPopulation;
  void setD3QueuedFlag(Size u, Size i, char value);
};

inline void FastMatchingEngineV2::setD3QueuedFlag(Size u, Size i, char value)
{
  if (mBackend == kGabow2018) mD3QueuedPopulation.set(mD3Queued, u, i, value);
  else mD3Queued[u][i] = value;
}

#include "MonotoneSupport.inc"

inline void FastMatchingEngineV2::compute(const Graph& graph,
                                          std::vector<Size>* mateVec,
                                          Size* card,
                                          Val* edgWght)
{
#ifdef MBX_GABOW_FRONTIER_AUDIT
  FrontierReuseAudit::Session frontierSession((int)mBackend);
#endif

#ifdef MBX_LARGEST_REUSE_AUDIT
  LargestReuseAudit::Session largestSession((int)mBackend);
#endif
  if (mateVec == NULL || card == NULL || edgWght == NULL) return;
#ifdef MBX_MEMBERSHIP_REUSE_AUDIT
  MembershipReuseAudit::Session membershipReuseSession((int)mBackend);
#endif
  init(graph);
#ifdef MBX_TWIN_CLEANUP_AUDIT
  TwinCleanupAudit::Session twinCleanupSession((int)mBackend);
#endif
#ifdef MBX_CLASSIC_EVENT_TRACE
  mEventTrace.begin(mNumVtxs,mBackend==kGabow2018?"gabow":"edmonds",mAdj,mWgt);
  mEventMembers.clear();
  if (mEventTrace.enabled()) {
    mEventMembers.resize(mBlsm.size());
    for (Size v=0;v<mNumVtxs;++v) mEventMembers[v]=ClassicEventTrace::Signature::atom(v);
  }
#endif
#ifdef MBX_MONOTONE_SUPPORT_AUDIT
  MonotoneSupportAudit::Session monotoneSession((int)mBackend);
#endif
  if(mBackend==kGabow2018)supportVerifyStable("initial");
  run();
  if(mBackend==kGabow2018&&!mAbortRun)supportVerifyStable("final");
  if (mAbortRun) {
    mateVec->assign(mNumVtxs, cNullItm);
    *card = 0;
    *edgWght = 0;
#ifdef MBX_CLASSIC_EVENT_TRACE
    mEventTrace.finish(true,mUnmatched,mDeltaSum,*card,*edgWght);
#endif
    return;
  }
  Size directUnmatched = 0;
  for (Size v = 0; v < mMate.size(); ++v) {
    if (mMate[v] == cNullItm) ++directUnmatched;
  }
  (void)directUnmatched;
  extractMatching(mateVec, card, edgWght);
#ifdef MBX_CLASSIC_EVENT_TRACE
  mEventTrace.finish(false,mUnmatched,mDeltaSum,*card,*edgWght);
#endif
#ifndef MBX_NO_DIAGNOSTIC_OUTPUT
  // Structural stats (both backends): shrink count, augment count, and the
  // shrink/augment ratio = "blossom-cascade intensity". High ratio predicts
  // a Gabow win (it's the regime that floods GMG'86's global PQ with useless
  // edges, the Section 2.3 pathology). Low ratio predicts a GMG win.
  {
    double ratio = (mProfAugments > 0)
      ? (double)mProfShrinks / (double)mProfAugments : 0.0;
    std::fprintf(stderr,
      "STRUCT %s: shrinks=%lu augments=%lu shrink/aug=%.3f pushD3calls=%lu\n",
      (mBackend == kGabow2018) ? "gabow" : "gmg",
      (unsigned long)mProfShrinks, (unsigned long)mProfAugments, ratio,
      (unsigned long)mProfPushD3CallsAll);
    // Event-loop time breakdown (both backends): localizes the previously-
    // "unaccounted" cost. iters = event-loop iterations (~min_edge calls).
    {
      long double total = mProfWallStart > 0 ? (profNow() - mProfWallStart) : 0;
      std::fprintf(stderr,
        "PROFLOOP %s: iters=%llu  D124scan=%.2Lfs extend=%.2Lfs shrink=%.2Lfs "
        "augment=%.2Lfs unmatch=%.2Lfs  minEdge=%.2Lfs evict=%.2Lfs "
        "pushD3=%.2Lfs staleChk=%.2Lfs  wall~=%.2Lfs\n",
        (mBackend == kGabow2018) ? "gabow" : "gmg",
        (unsigned long long)mProfLoopIters,
        mProfD124ScanTime, mProfActExtendTime, mProfActShrinkTime,
        mProfActAugmentTime, mProfActUnmatchTime,
        mProfMinEdgeTime, mProfEvictTime, mProfPushD3Time, mProfStaleCheckTime,
        total);
      std::fprintf(stderr,
        "  SHRINKDETAIL %s: refreshInShrink=%.2Lfs vtxBuild(+childOfVtxMap)=%.2Lfs cqMerge=%.2Lfs",
        (mBackend == kGabow2018) ? "gabow" : "gmg", mProfShrinkRefreshTime,
        mProfShrinkVtxBuildTime, mProfShrinkCqMergeTime);
      if (mBackend == kGabow2018)
        std::fprintf(stderr, "  mBMhook=%.2Lfs discardedRepushes=%llu",
          mProfShrinkHookTime, mProfShrinkDiscarded);
      std::fprintf(stderr, "\n");
    }
    if (mBackend == kGabow2018) {
      {
        unsigned long long dirEdges = 0;  // directed edge slots = 2m
        for (Size vv = 0; vv < mAdj.size(); ++vv) dirEdges += mAdj[vv].size();
        std::fprintf(stderr,
          "MBMEDGES gabow: peakAliveEdges=%llu  (2m=%llu, ratio peak/2m=%.2f)\n",
          bmPeakAliveEdges(), dirEdges,
          dirEdges > 0 ? (double)bmPeakAliveEdges() / (double)dirEdges : 0.0);
      }
      if (mGabowLoopDiag) {
        double avgMult = (mProfDistinctEdgesPushed > 0)
          ? (double)mProfTotalEdgePushes / (double)mProfDistinctEdgesPushed : 0.0;
        std::fprintf(stderr,
          "REPUSH gabow: totalEdgePushes=%llu distinctEdges=%llu "
          "avgMultiplicity=%.2f maxMultiplicity=%llu\n",
          mProfTotalEdgePushes, mProfDistinctEdgesPushed,
          avgMult, mProfMaxPushMult);
      } else {
        std::fprintf(stderr,
          "REPUSH gabow: (disabled; set GABOW_LOOPDIAG=1 to measure)\n");
      }
    }
  }
  if (mBackend == kGabow2018) {
    long double measured = mProfPushD3Time + mProfMinEdgeTime +
                            mProfMergeTime + mProfDisposeTime +
                            mProfStaleCheckTime;
    std::fprintf(stderr,
      "PROF gabow:\n"
      "  pushD3=%lu (%.3Lfs total, early=%lu)\n"
      "    inner: make_edge=%lu (%.3Lfs)  nca=%lu (%.3Lfs)\n"
      "  merge=%lu (%.3Lfs)  dispose=%lu (%.3Lfs)\n"
      "  min_edge=%lu (%.3Lfs)  stalePops=%lu (staleCheck=%.3Lfs)\n"
      "  matched-retire=%lu (%.3Lfs) slots=%llu live-handles=%llu\n"
      "  total measured = %.3Lfs\n",
      (unsigned long)mProfPushD3Count, mProfPushD3Time,
      (unsigned long)mProfPushD3Early,
      (unsigned long)mProfMakeEdgeCount, mProfMakeEdgeTime,
      (unsigned long)mProfNcaCount, mProfNcaTime,
      (unsigned long)mProfMergeCount, mProfMergeTime,
      (unsigned long)mProfDisposeCount, mProfDisposeTime,
      (unsigned long)mProfMinEdgeCount, mProfMinEdgeTime,
      (unsigned long)mProfStalePopCount, mProfStaleCheckTime,
      (unsigned long)mProfEvictCalls, mProfEvictTime,
      mProfEvictSlotsScanned, mProfEvictLiveRemovals,
      measured + mProfEvictTime);
    std::fprintf(stderr,
      "BMREFRESH gabow: calls=%llu packetSteps=%llu maxPackets=%llu\n",
      bmRefreshCount(), bmRefreshPacketSteps(), bmMaxBlossomPackets());
    // R_Scope structural validation: a faithful forest must show maxDepth>1 and
    // nonRootNca>0 on grow-heavy inputs; the permanent star shows 1 and 0.
    std::fprintf(stderr,
      "FOREST gabow: deepT=%s maxDepth=%lu nonRootNca=%llu (of %lu same-tree ncas) fallbackInstalls=%llu staleQueued=%llu staleKept=%llu staleFired=%llu rebuiltValid=%llu\n",
      mGabowDeepForest ? "ON" : "off",
      (unsigned long)mAuxNCA.maxDepth(), mProfNcaNonRoot,
      (unsigned long)mProfNcaCount, mProfForestFallback,
      mProfForestStaleQueued, mProfForestStaleKept, mProfForestStaleFired,
      mProfForestStaleRebuilt);
  }
#endif
#ifdef MBX_GABOW_FRONTIER_VERIFY
  verifyGabowFrontierStable(true);
#endif
#ifdef MBX_GABOW_SEARCH_AUDIT
  if (mBackend == kGabow2018) {
    audEmitSearch('f');   // trailing partial interval (after the last augment)
    audPrintSummary();
  }
  audPrintChassis();      // paired shared-chassis accounting, BOTH backends
#endif
}

inline void FastMatchingEngineV2::run()
{
  start();
}

#ifdef MBX_GABOW_SEARCH_AUDIT
inline const char* FastMatchingEngineV2::audName(int i)
{
  static const char* n[kAudN] = {
    "loopIters",     "pushD3",        "pushDupEarly",  "makeEdge",
    "recSingle",     "recDouble",     "nca",           "ncaNonRoot",
    "ncaAncestor",   "dncaNcaSteps",  "dncaAddSteps",  "addLeafFirst",
    "addLeafReentry","retires",       "fallbackInstall","staleQueued",
    "bmMerges",      "minEdge",       "stalePops",     "disposes",
    "shrinks",       "splits",        "treeAdds",      "repushDiscard",
    "totalPushes",   "distinctPushed","bmToPacketZ",   "bmToLoose",
    "bmLooseToPkt",  "bmPktToLoose",  "bmReclass",     "bmDiscSame",
    "bmDiscDedup",   "mergeUpSteps",  "mergeDownSteps","ncaTimeUs",
    "makeEdgeTimeUs","minEdgeTimeUs", "mergeTimeUs",   "pushD3TimeUs",
    "ivDupPushes",   "grows",          "r8Reduced",     "r8OneSide",
    "r8TwoSide",     "r8DummyRoot",    "r8Ancestor" };
  return n[i];
}

inline void FastMatchingEngineV2::audCollect(unsigned long long* o) const
{
  o[0]  = mProfLoopIters;                  o[1]  = mProfPushD3Count;
  o[2]  = mProfPushD3Early;                o[3]  = mProfMakeEdgeCount;
  o[4]  = mAudSingleRec;                   o[5]  = mAudDoubleRec;
  o[6]  = mProfNcaCount;                   o[7]  = mProfNcaNonRoot;
  o[8]  = mAudNcaAncestor;                 o[9]  = dncaNcaSteps();
  o[10] = dncaAddLeafSteps();              o[11] = mAudAddLeafFirst;
  o[12] = mAudAddLeafReentry;              o[13] = mAudRetires;
  o[14] = mProfForestFallback;             o[15] = mProfForestStaleQueued;
  o[16] = mProfMergeCount;                 o[17] = mProfMinEdgeCount;
  o[18] = mProfStalePopCount;              o[19] = mProfDisposeCount;
  o[20] = mProfShrinks;                    o[21] = mAudSplits;
  o[22] = mAudTreeAdds;                    o[23] = mProfShrinkDiscarded;
  o[24] = mProfTotalEdgePushes;            o[25] = mProfDistinctEdgesPushed;
  o[26] = bmR22ToPacketZ();                o[27] = bmR22ToLoose();
  o[28] = bmLooseToPacket();               o[29] = bmPacketToLoose();
  o[30] = bmReclassCount();                o[31] = bmDiscardSameBlsm();
  o[32] = bmDiscardDedup();                o[33] = mAudMergeUpSteps;
  o[34] = mAudMergeDownSteps;              o[35] = (unsigned long long)(mProfNcaTime * 1e6L);
  o[36] = (unsigned long long)(mProfMakeEdgeTime * 1e6L);
  o[37] = (unsigned long long)(mProfMinEdgeTime * 1e6L);
  o[38] = (unsigned long long)(mProfMergeTime * 1e6L);
  o[39] = (unsigned long long)(mProfPushD3Time * 1e6L);
  o[40] = mAudIvDupPushes;
  o[41] = mAudGrows;
  o[42] = mAudR8Reduced;
  o[43] = mAudR8OneSide;
  o[44] = mAudR8TwoSide;
  o[45] = mAudR8DummyRoot;
  o[46] = mAudR8Ancestor;
}

inline void FastMatchingEngineV2::audNoteIntervalPush(Size u, Size edgeIdx)
{
  if (u >= mAudIvStamp.size() || edgeIdx >= mAudIvStamp[u].size()) return;
  const unsigned iv = (unsigned)(mAudSearches + 1);   // current (unemitted) interval
  if (mAudIvStamp[u][edgeIdx] != iv) {
    mAudIvStamp[u][edgeIdx] = iv;
    mAudIvCount[u][edgeIdx] = 1;
  } else {
    ++mAudIvCount[u][edgeIdx];
    ++mAudIvDupPushes;
    if (mAudIvCount[u][edgeIdx] > mAudIvMaxMult) mAudIvMaxMult = mAudIvCount[u][edgeIdx];
  }
}

inline void FastMatchingEngineV2::audNoteAddLeaf(Size w)
{
  if (w >= mAudSeen.size()) return;
  if (mAudSeen[w]) ++mAudAddLeafReentry;
  else { mAudSeen[w] = 1; ++mAudAddLeafFirst; }
}

inline void FastMatchingEngineV2::audEmitSearch(char kind)
{
  unsigned long long cur[kAudN];
  audCollect(cur);
  unsigned long long d[kAudN];
  for (int i = 0; i < kAudN; ++i) {
    d[i] = cur[i] - mAudPrev[i];
    mAudSum[i] += d[i];
    if (d[i] > mAudMax[i]) mAudMax[i] = d[i];
    mAudPrev[i] = cur[i];
  }
  ++mAudSearches;
  if (mAudEmitPerSearch) {
    // IAUD = one AUGMENTATION INTERVAL (kind '0' = the initialization flood,
    // 'a'/'e' = arc/edge augmentation, 'f' = trailing partial interval).
    // These are NOT Gabow's independent per-search budgets: surviving trees
    // and records persist across intervals in the continuous engine.
    std::fprintf(stderr, "IAUD,%llu,%c", mAudSearches, kind);
    for (int i = 0; i < kAudN; ++i)
      if (d[i]) std::fprintf(stderr, ",%s=%llu", audName(i), d[i]);
    std::fprintf(stderr, "\n");
  }
}

inline void FastMatchingEngineV2::audPrintSummary() const
{
  std::fprintf(stderr, "IAUDSUM intervals=%llu n=%lu (metric total perIntervalMax mean; "
               "interval '0' = initialization flood, reported separately in IAUD)\n",
               mAudSearches, (unsigned long)mNumVtxs);
  for (int i = 0; i < kAudN; ++i)
    std::fprintf(stderr, "IAUDSUM %-15s %12llu %10llu %10.1f\n", audName(i),
      mAudSum[i], mAudMax[i],
      mAudSearches ? (double)mAudSum[i] / (double)mAudSearches : 0.0);
  std::fprintf(stderr,
    "IAUDPEAK maxForestDepth=%lu dncaMaxLift=%llu bmPeakEdges=%llu bmLoosePeak=%llu "
    "maxLoose=%llu maxPkt=%llu maxPushMultRun=%llu maxPushMultInterval=%llu "
    "mergeChainMax=%llu mergeChains=%llu\n",
    (unsigned long)mAuxNCA.maxDepth(), dncaMaxLift(), bmPeakAliveEdges(),
    bmLooseLivePeak(), bmMaxLooseLen(), bmMaxPacketLen(), mProfMaxPushMult,
    mAudIvMaxMult, mAudMergeChainMax, mAudMergeChains);
  // Work-order ratios (denominators as available; 0-guarded).
  unsigned long long mkE = mAudSum[3], push = mAudSum[1], ncaQ = mAudSum[6],
    mrg = mAudSum[16], minE = mAudSum[17], tAdd = mAudSum[22],
    tot = mAudSum[24], dst = mAudSum[25];
  std::fprintf(stderr,
    "IAUDRATIO makeEdge/pushD3=%.3f physIns/distinctEdge=%.3f "
    "repushRunGlobal=%.3f ivDup/pushD3=%.3f nca/pushD3=%.3f "
    "merges/treeAdds=%.4f minEdge/treeAdds=%.4f liftsPerNca=%.2f\n",
    push ? (double)mkE / push : 0.0,
    dst ? (double)mkE / dst : 0.0,
    dst ? (double)(tot - dst) / dst : 0.0,
    push ? (double)mAudSum[40] / push : 0.0,
    push ? (double)ncaQ / push : 0.0,
    tAdd ? (double)mrg / tAdd : 0.0,
    tAdd ? (double)minE / tAdd : 0.0,
    ncaQ ? (double)mAudSum[9] / ncaQ : 0.0);
}

inline void FastMatchingEngineV2::audPrintChassis() const
{
  // Stage-1 paired shared-chassis accounting: printed for BOTH backends so
  // GMG-vs-Gabow differences in SHARED work are directly comparable.
  std::fprintf(stderr,
    "CHASSIS backend=%s n=%lu loopIters=%llu evD1=%llu evD2=%llu evD3=%llu evD4=%llu "
    "grows=%llu shrinks=%lu splits=%llu augments=%lu topLevel=%llu "
    "d2CandAdds=%llu d2BestPushes=%llu d2Drops=%llu "
    "d3Logical=%lu d3PhysGabow=%lu d3PhysGmg=%llu d3StaleGabow=%lu d3StaleGmg=%llu\n",
    (mBackend == kGabow2018) ? "gabow" : "gmg", (unsigned long)mNumVtxs,
    mProfLoopIters, mAudEvD1, mAudEvD2, mAudEvD3, mAudEvD4,
    mAudGrows, (unsigned long)mProfShrinks, mAudSplits, (unsigned long)mProfAugments,
    mAudTopLevelCalls, mAudD2CandAdds, mAudD2BestPushes, mAudD2Drops,
    (unsigned long)mProfPushD3CallsAll, (unsigned long)mProfMakeEdgeCount,
    mAudGmgD3Pushes, (unsigned long)mProfStalePopCount, mAudGmgD3Stale);
}
#endif // MBX_GABOW_SEARCH_AUDIT


inline Val FastMatchingEngineV2::effY(Size v) const
{
  Size b = topLevel(v);
  if (b == cNullItm || b >= mBlsm.size()) return mVtx[v].pot;
  if (mBlsm[b].status == BS_EVEN)
    return mVtx[v].pot - (mDeltaSum - mBlsm[b].offset);
  if (mBlsm[b].status == BS_ODD)
    return mVtx[v].pot + (mDeltaSum - mBlsm[b].offset);
  return mVtx[v].pot;
}

inline Val FastMatchingEngineV2::effZ(Size b) const
{
  if (mBlsm[b].status == BS_EVEN)
    return mBlsm[b].pot + 2.0 * (mDeltaSum - mBlsm[b].offset);
  if (mBlsm[b].status == BS_ODD)
    return mBlsm[b].pot - 2.0 * (mDeltaSum - mBlsm[b].offset);
  return mBlsm[b].pot;
}

inline Val FastMatchingEngineV2::edgeSlack(Size u, Size edgeIdx) const
{
  Size v = mAdj[u][edgeIdx];
  return effY(u) + effY(v) - mWgt[u][edgeIdx];
}

inline Val FastMatchingEngineV2::d2Key(Size matchedVtx, Size edgeIdxInMatched) const
{
  Val key = mDeltaSum + edgeSlack(matchedVtx, edgeIdxInMatched);
  return key < mDeltaSum ? mDeltaSum : key;
}

inline Val FastMatchingEngineV2::d3Key(Size u, Size edgeIdx) const
{
  Val key = mDeltaSum + edgeSlack(u, edgeIdx) / 2.0;
  return key < mDeltaSum ? mDeltaSum : key;
}

inline void FastMatchingEngineV2::init(const Graph& graph)
{
  Size numEdgs = 0;
  graph.GetSize(&mNumVtxs, &numEdgs);
  mUnmatched = mNumVtxs;
  mDeltaSum = 0;
  mAbortRun = false;

  // Build adjacency into temporary jagged vectors (self-loops filtered, so per-row
  // length isn't known until after the scan), then flatten into the CSR FlatJagged.
  std::vector<std::vector<Size> > tmpAdj(mNumVtxs);
  std::vector<std::vector<Val> >  tmpWgt(mNumVtxs);
  for (Size v = 0; v < mNumVtxs; ++v) {
    Size deg = 0;
    std::vector<Size> rawAdj;
    std::vector<Val> rawWgt;
    graph.GetAdj(v, &deg, &rawAdj, &rawWgt);
    tmpAdj[v].reserve(rawAdj.size());
    tmpWgt[v].reserve(rawWgt.size());
    for (Size i = 0; i < rawAdj.size(); ++i) {
      if (rawAdj[i] == v) continue;
      tmpAdj[v].push_back(rawAdj[i]);
      tmpWgt[v].push_back(rawWgt[i]);
    }
  }
  std::vector<Size> rowLens(mNumVtxs, 0);
  for (Size u = 0; u < mNumVtxs; ++u) rowLens[u] = tmpAdj[u].size();
  mAdj.assign(rowLens, (Size)0);
  mWgt.assign(rowLens, (Val)0);
  for (Size u = 0; u < mNumVtxs; ++u) {
    FlatJagged<Size>::Row ar = mAdj[u];
    FlatJagged<Val>::Row  wr = mWgt[u];
    for (Size i = 0; i < tmpAdj[u].size(); ++i) { ar[i] = tmpAdj[u][i]; wr[i] = tmpWgt[u][i]; }
  }
  // Graph::rFormAdjs supplies sorted, duplicate-free, symmetrically stored
  // adjacency rows.  Pair each undirected edge once by binary-searching the
  // reciprocal sorted row, then fill both reverse slots together.  The former
  // implementation built and sorted a 2m-entry (encoded-edge,slot) table;
  // besides its O(m log m) work and temporary memory, encodeEdge(u,v,n) could
  // overflow for sufficiently large Size products.  The graph invariant also
  // means there is no parallel-slot ambiguity here: input parallel edges have
  // already been reduced exactly to their maximum-weight representative.
  mRevIdx.assign(rowLens, cNullItm);
  Size pairedEdgs = 0;
  for (Size u = 0; u < mNumVtxs; ++u) {
    for (Size i = 0; i < mAdj[u].size(); ++i) {
      const Size v = mAdj[u][i];
      if (v >= mNumVtxs || v <= u) continue;
      const FlatJagged<Size>::ConstRow vr =
        static_cast<const FlatJagged<Size>&>(mAdj)[v];
      Size lo = 0;
      Size hi = vr.size();
      while (lo < hi) {
        const Size mid = lo + (hi - lo) / 2;
        if (vr[mid] < u) lo = mid + 1;
        else hi = mid;
      }
      if (lo < vr.size() && vr[lo] == u) {
        mRevIdx[u][i] = lo;
        mRevIdx[v][lo] = i;
        ++pairedEdgs;
#ifndef NDEBUG
        assert(mWgt[u][i] == mWgt[v][lo]);
#endif
      }
    }
  }
#ifndef NDEBUG
  assert(pairedEdgs == numEdgs);
  for (Size u = 0; u < mNumVtxs; ++u)
    for (Size i = 0; i < mAdj[u].size(); ++i) {
      const Size v = mAdj[u][i];
      assert(v < mNumVtxs && mRevIdx[u][i] < mAdj[v].size());
      assert(mAdj[v][mRevIdx[u][i]] == u);
      assert(mRevIdx[v][mRevIdx[u][i]] == i);
    }
#else
  (void)pairedEdgs;
  (void)numEdgs;
#endif
  mD3Queued.assign(rowLens, (char)0);
  if (mBackend == kGabow2018) mD3QueuedPopulation.resetZeroRows(mNumVtxs);
#ifdef MBX_GMG_D3_OWNERSHIP_EXPERIMENT
  if (mBackend == kGMG1986)
    mGmgD3Publication.assign(rowLens, GmgD3Publication());
  else mGmgD3Publication.clear();
#endif
  // Re-push multiplicity table (parallel to mD3Queued). This per-edge array
  // costs ~16 B/edge and is pure instrumentation, so it is only allocated when
  // diagnostics are enabled (GABOW_LOOPDIAG=1). When empty, the pushD3
  // increment is skipped (its size() guard fails) and the REPUSH report is
  // omitted. Allocated below, after mGabowLoopDiag is read from the env.
  mPushMult.clear();
  // Reset profiling counters.
  mProfResetCount = 0; mProfResetTime = 0;
  mProfPushD3Count = 0; mProfPushD3Early = 0; mProfPushD3Time = 0;
  mProfMakeEdgeCount = 0; mProfMakeEdgeTime = 0;
  mProfMergeCount = 0; mProfMergeTime = 0;
  mProfDisposeCount = 0; mProfDisposeTime = 0;
  mProfMinEdgeCount = 0; mProfMinEdgeTime = 0;
  mProfStalePopCount = 0; mProfStaleCheckTime = 0;
  mProfNcaCount = 0; mProfNcaTime = 0;
  mProfStaleCostRepush = 0;
  mProfShrinks = 0; mProfAugments = 0;
  mProfEvictCalls = 0; mProfEvictTime = 0;
  mAnyShrinkYet = false;
  mProfDestroyTreeTime = 0; mProfD2CandTime = 0;
  mProfD124ScanTime = 0; mProfActExtendTime = 0; mProfActShrinkTime = 0;
  mProfShrinkHookTime = 0; mProfShrinkDiscarded = 0; mProfShrinkRefreshTime = 0;
  mProfShrinkVtxBuildTime = 0; mProfShrinkCqMergeTime = 0;
  mProfEvictSlotsScanned = 0; mProfEvictLiveRemovals = 0;
  mProfActAugmentTime = 0; mProfActUnmatchTime = 0; mProfLoopIters = 0;
  mProfPushD3CallsAll = 0;
  mProfStaleSameBlsm = 0; mProfStaleOdd = 0; mProfStaleMatched = 0;
  mProfStaleInactiveTree = 0; mProfStaleOther = 0;
  mProfTotalEdgePushes = 0; mProfDistinctEdgesPushed = 0; mProfMaxPushMult = 0;
  mProfCleanupCount = 0; mProfCleanupTime = 0; mProfCleanupEdgesRemoved = 0;
  mProfTreeMemberCapacity = 0;
  mProfTreeMemberMaxRowCapacity = 0;
#ifdef MBX_GABOW_SEARCH_AUDIT
  mAudAddLeafFirst = mAudAddLeafReentry = mAudRetires = 0;
  mAudNcaAncestor = mAudSingleRec = mAudDoubleRec = 0;
  mAudR8Reduced = mAudR8OneSide = mAudR8TwoSide = 0;
  mAudR8DummyRoot = mAudR8Ancestor = 0;
  mAudSplits = mAudTreeAdds = 0;
  mAudMergeUpSteps = mAudMergeDownSteps = mAudMergeChains = 0;
  mAudMergeChainMax = mAudMergeOrderLogged = mAudSearches = 0;
  mAudIvDupPushes = mAudIvMaxMult = 0;
  mAudTopLevelCalls = 0;
  mAudEvD1 = mAudEvD2 = mAudEvD3 = mAudEvD4 = mAudGrows = 0;
  mAudD2CandAdds = mAudD2BestPushes = mAudD2Drops = 0;
  mAudGmgD3Pushes = mAudGmgD3Stale = 0;
  mAudIvStamp.assign(mNumVtxs, std::vector<unsigned>());
  mAudIvCount.assign(mNumVtxs, std::vector<unsigned>());
  for (Size v = 0; v < mNumVtxs; ++v) {
    mAudIvStamp[v].assign(mAdj[v].size(), 0U);
    mAudIvCount[v].assign(mAdj[v].size(), 0U);
  }
  mAudSeen.assign(mNumVtxs + 1, 0);
  for (int ai = 0; ai < kAudN; ++ai) { mAudPrev[ai] = mAudSum[ai] = mAudMax[ai] = 0; }
  {
    const char* ev = std::getenv("GABOW_SEARCH_AUDIT_EMIT");
    mAudEmitPerSearch = (ev != NULL && ev[0] == '1');
  }
  // Process-global structure counters restart with this engine instance so
  // per-search deltas are attributable (audit builds only).
  dncaNcaCalls() = dncaNcaSteps() = dncaAddLeafCalls() = dncaAddLeafSteps() = 0;
  dncaMaxLift() = 0;
  bmR22ToPacketZ() = bmR22ToLoose() = bmLooseToPacket() = bmPacketToLoose() = 0;
  bmReclassCount() = bmDiscardSameBlsm() = bmDiscardDedup() = 0;
  bmLooseLive() = bmLooseLivePeak() = bmMaxLooseLen() = 0;
#endif

  // Set up the Gabow 2018 backend members. Allocated unconditionally so
  // members are in a defined state; only consulted when mBackend == kGabow2018.
  //
  // Vertices 0..mNumVtxs-1 are real vertices; mAuxRoot == mNumVtxs is the
  // synthetic dummy root r per paper Section 3 ("Each supporting tree is
  // rooted at a child of r"). Capacity = mNumVtxs + 1.
  //
  // F2 simplification: Matchbox starts with every vertex as a free EVEN
  // search-tree root, so every vertex is installed as a child of rho. This
  // permanent-star adaptation avoids per-augment rebuilding, but it also
  // trivializes NCA depth after the first augmentation; it must not be used as
  // evidence that a faithful dynamic supporting forest has negligible cost.
  mAuxRoot = mNumVtxs;
  mAuxNCA.reset(mNumVtxs + 1);
  mInT.assign(mNumVtxs + 1, 0);
  // mBM also needs room for mAuxRoot (back-edges have a = mAuxRoot under our
  // depth-1 init). +1 for the dummy root vertex slot.
  mBM.reset(mNumVtxs + 1);
  // Gabow generally uses lazy invalidation.  The one narrow eager policy is
  // matched-blossom retirement: when augmentation destroys a whole search
  // tree, queued D3 records incident to its now-matched blossoms are removed
  // through their exact generation-guarded handles.  This is not BM-owner-wide
  // deactivation and does not run for internal label changes during shrink.
  // GABOW_LOOPDIAG=1 enables a time-based event-loop progress dump (every ~10s
  // wall) to stderr. Off by default so the per-event profNow() clock syscall
  // doesn't tax the hot loop on normal timed runs. Used to localize bottlenecks
  // (e.g. the hugebubbles init-flood) where the D3 consumer rarely fires.
  {
    const char* ev = std::getenv("GABOW_LOOPDIAG");
    mGabowLoopDiag = (ev != NULL && ev[0] == '1');
  }
  {
#ifdef MBX_GABOW_NOREPUSH_EXPERIMENT
    const char* ev = std::getenv("GABOW_NOREPUSH");
    mGabowNoRepushExpt = (ev != NULL && ev[0] == '1');
#else
    // The no-republish path is known incorrect.  Production binaries must not
    // change mathematical behavior because of an inherited environment value.
    mGabowNoRepushExpt = false;
#endif
  }
#ifdef MBX_GABOW_EDGE_WATCH
  {
    const char* eu = std::getenv("GABOW_EDGE_U");
    const char* ev = std::getenv("GABOW_EDGE_V");
    mWatchEdgeEnabled = eu != NULL && ev != NULL;
    mWatchEdgeU = mWatchEdgeEnabled ? (Size)std::strtoul(eu, NULL, 10) : cNullItm;
    mWatchEdgeV = mWatchEdgeEnabled ? (Size)std::strtoul(ev, NULL, 10) : cNullItm;
    mWatchLastSlack = 0;
    mWatchHaveSlack = false;
  }
#endif
  {
#ifdef MBX_GABOW_STAR_EXPERIMENT
    // Experiment-only A/B: GABOW_DEEPT=0 forces the RETIRED star geometry
    // (unfaithful; measurement builds only). Absent in production builds.
    const char* ev = std::getenv("GABOW_DEEPT");
    mGabowDeepForest = !(ev != NULL && ev[0] == '0');
    mGabowDeepSplitClear = mGabowDeepForest;
#else
    // R7 (adopted): the faithful supporting forest IS the production Gabow
    // geometry. No user-facing star/deep policy choice remains; the old star
    // mode lives in git history and behind MBX_GABOW_STAR_EXPERIMENT.
    mGabowDeepForest = true;
    mGabowDeepSplitClear = true;
#endif
  }
  {
    const char* ev = std::getenv("GABOW_FORESTCHECK");
    mGabowForestCheck = (ev != NULL && ev[0] == '1');
  }
  if (mBackend == kGabow2018 && mGabowDeepForest)
    mForestGeneration.assign(mNumVtxs, 1U);
  else
    mForestGeneration.clear();
  mProfNcaNonRoot = 0;
  mProfForestFallback = 0;
  mProfForestStaleQueued = 0;
  mProfForestStaleFired = 0;
  mProfForestStaleKept = 0;
  mProfForestStaleRebuilt = 0;
  {
    const char* ev = std::getenv("GABOW_D3LOG");
    mGabowD3Log = (ev != NULL && ev[0] == '1');
  }
  // Allocate the re-push multiplicity table only when diagnostics are on
  // (saves ~16 B/edge in normal runs). Gabow-only; never touches the GMG path.
  if (mBackend == kGabow2018 && mGabowLoopDiag) {
    mPushMult.assign(mNumVtxs, std::vector<Size>());
    for (Size v = 0; v < mNumVtxs; ++v)
      mPushMult[v].assign(mAdj[v].size(), 0);
  }
  // Stride for payload packing.
  mAdjStride = 1;
  for (Size v = 0; v < mNumVtxs; ++v) {
    if (mAdj[v].size() > mAdjStride) mAdjStride = mAdj[v].size();
  }
  if (mAdjStride == 0) mAdjStride = 1;
#ifdef MBX_SHARED_D3_TIE_ORDER
  // Every packed original-edge ordinal is below cNullItm. Check before any
  // multiplication, not after an overflow. This metadata is fixed per graph.
  if (mNumVtxs && mAdjStride > cNullItm / mNumVtxs)
    throw std::overflow_error("shared D3 original-edge ordinal overflow");
#endif
  // Per-(u, edgeIdx) handle tracking, sized to match mAdj. CSR-flat (one heap
  // block each, not mNumVtxs tiny vectors). Row lengths come from mAdj.
  {
    std::vector<Size> rowLens(mNumVtxs);
    for (Size v = 0; v < mNumVtxs; ++v) rowLens[v] = mAdj[v].size();
    mBmHandleU.assign(rowLens, (BlossomMerging::EdgeHandle)cNullItm);
    mBmHandleV.assign(rowLens, (BlossomMerging::EdgeHandle)cNullItm);
    mBmGenU.assign(rowLens, (Size)0);
    mBmGenV.assign(rowLens, (Size)0);
  }
  // mMatchboxToBm mapping. Initial identity for vertices (matchbox blsmId =
  // mBM blossomId = v for v < mNumVtxs); root slot handled below.
  // Sized later when mBlossomCap is known.
  if (mBackend == kGabow2018) {
    // Root r has no parent.
    mAuxNCA.addLeaf(mAuxRoot, cNullItm);
    mInT[mAuxRoot] = 1;
    mBM.set_depth(mAuxRoot, 0);
    for (Size v = 0; v < mNumVtxs; ++v) {
      mAuxNCA.addLeaf(v, mAuxRoot);
      mInT[v] = 1;
      mBM.set_depth(v, 1);
    }
#ifdef MBX_GABOW_SEARCH_AUDIT
    // The engine PRE-INSTALLS every vertex as a depth-1 star child (F2). Mark
    // them seen: every later search-driven addLeaf is therefore a RE-INSTALL
    // (addLeafReentry); addLeafFirst stays 0 by construction under this model.
    for (Size v = 0; v <= mNumVtxs; ++v) mAudSeen[v] = 1;
#endif
    // Each vertex starts as its own singleton blossom in mBM. Root is also
    // a singleton blossom so back-edges with a = root are valid.
    std::vector<Size> singleton(1);
    for (Size v = 0; v < mNumVtxs; ++v) {
      singleton[0] = v;
      mBM.make_blossom(singleton);
    }
    singleton[0] = mAuxRoot;
    mBM.make_blossom(singleton);
    mSupportHead.assign(mNumVtxs,cNullItm);
    mSupportNext.assign(mNumVtxs+1,cNullItm);
    mSupportTreeOf.assign(mNumVtxs+1,cNullItm);
    mSupportOwnerStamp.assign(mNumVtxs+1,0);mSupportOwnerEpoch=0;
    mSupportSingleton.assign(1,0);
    for(Size v=0;v<mNumVtxs;++v){mSupportHead[v]=v;mSupportTreeOf[v]=v;}

  }

  mBlossomCap = mNumVtxs + (mNumVtxs / 2) + 1;
  // A new graph discards every old membership. Free old per-slot buffers
  // before vector assignment can retain a large prior frontier as a tiny/dead
  // header in the new graph. Fresh benchmark engines have no old slots here.
  if(mBackend==kGabow2018) {
    for(Size b=0;b<mBlsm.size();++b) {
      mBlsm[b].vertices.release();
      std::vector<std::pair<Size,Size> >().swap(mBlsm[b].childOfVtx);
    }
    std::vector<std::vector<Size> >().swap(mSplitChildVtxsScratch);
  }
  mBlsm.assign(mBlossomCap, BlossomData());
  // mMatchboxToBm mapping: matchbox blsmId for trivial blsms 0..mNumVtxs-1
  // maps to mBM blossomId of the same index (set at init below); super-blsm
  // IDs assigned later by acquireBlsmId get filled in shrinkOnEdge.
  mMatchboxToBm.assign(mBlossomCap, cNullItm);
  if (mBackend == kGabow2018) {
    for (Size v = 0; v < mNumVtxs; ++v) mMatchboxToBm[v] = v;
  }
  mVtx.assign(mNumVtxs, VertexData());
  mTopBlsmOfVtx.assign(mNumVtxs, cNullItm);
  mTopCacheOfVtx.assign(mNumVtxs, cNullItm);
  mTopCacheGenerationOfVtx.assign(mNumVtxs, 0);
  mMate.assign(mNumVtxs, cNullItm);
  mTreeParent.assign(mNumVtxs, cNullItm);
  mTreeMembers.assign(mNumVtxs, std::vector<Size>());
  mBestD2Stamp.assign(mNumVtxs, 0);
  mBestD2Edge.assign(mNumVtxs, cNullItm);
  mBestD2Touched.clear();
  mBestD2Key.assign(mNumVtxs, infKey());
  mBestD2CurrentStamp = 0;
  mShrinkStamp.assign(mBlossomCap, 0);
  mShrinkLeftIndex.assign(mBlossomCap, cNullItm);
  mShrinkCurrentStamp = 0;
  mMaterializeDepth = 0;
  mFreeBlsmIds.clear();
  mCQ.clear();
  mPQD1 = PQD1();
  mPQD2 = PQD2();
  mPQD3 = PQD3();
  mPQD4 = PQD4();

  const Val cardTieEps = 0.0;
  for (Size u = 0; u < mNumVtxs; ++u) {
    for (Size i = 0; i < mWgt[u].size(); ++i) mWgt[u][i] += cardTieEps;
  }

  for (Size v = 0; v < mNumVtxs; ++v) {
    Val maxIncident = 0;
    for (Size i = 0; i < mWgt[v].size(); ++i) {
      if (mWgt[v][i] > maxIncident) maxIncident = mWgt[v][i];
    }
    mVtx[v].pot = maxIncident / 2.0;
    mVtx[v].bestByTree.clear();
    mVtx[v].d1Key = infKey();
    mVtx[v].d1Version = 0;
    mVtx[v].d1Queued = false;
    mVtx[v].d2Key = infKey();
    mVtx[v].d2Edge = cNullItm;
    mVtx[v].d2Blsm = cNullItm;
    mVtx[v].d2Tree = cNullItm;
    mVtx[v].d2Version = 0;
    mVtx[v].d2Queued = false;
    mVtx[v].lastEvenScanTree = cNullItm;
    mVtx[v].lastEvenScanTop = cNullItm;

    BlossomData& B = mBlsm[v];
    B.status = BS_EVEN;
    B.parent = cNullItm;
    B.base = v;
    B.tree = v;
    B.pred = cNullItm;
    B.next = cNullItm;
    B.pot = 0;
    B.offset = 0;
    B.children.clear();
    B.childEdges.clear();
    B.splitLeaves.clear();
    B.vertices.clear();
    B.vertices.push_back(v);
    B.childOfVtx.clear();
    B.cqRoot = cNullItm;
    B.cqLastLeaf = cNullItm;
    B.materializedBase = cNullItm;
    B.d1MinVtx = v;
    B.d1MinValid = true;
    B.d2Key = infKey();
    B.d2Vtx = cNullItm;
    B.d2Edge = cNullItm;
    B.d2Tree = cNullItm;
    B.d2Queued = false;
    B.generation = 1;
    B.topCache = v;
    B.topCacheGeneration = B.generation;
    B.listedTree = cNullItm;
    B.alive = true;

    mVtx[v].cqLeaf = mCQ.insertRight(B.cqRoot, infKey(), v);
    B.cqLastLeaf = mVtx[v].cqLeaf;
    mTopBlsmOfVtx[v] = v;
    mTopCacheOfVtx[v] = v;
    mTopCacheGenerationOfVtx[v] = B.generation;
    mTreeMembers[v].push_back(v);
    pushD1(v);
  }

  for (Size u = 0; u < mNumVtxs; ++u) {
    for (Size i = 0; i < mAdj[u].size(); ++i) {
      if (mAdj[u][i] > u) pushD3(u, i);
    }
  }
#ifdef MBX_GABOW_SEARCH_AUDIT
  // Close the INITIALIZATION interval ('0') so the init flood is reported
  // separately and never conflated with the first augmentation interval.
  if (mBackend == kGabow2018) audEmitSearch('0');
#endif
}

inline void FastMatchingEngineV2::extractMatching(std::vector<Size>* mateVec,
                                                  Size* card,
                                                  Val* edgWght) const
{
  std::vector<Size> recursiveMate(mNumVtxs, cNullItm);
  std::vector<Size> directMate(mNumVtxs, cNullItm);
  for (Size v = 0; v < mNumVtxs; ++v) directMate[v] = mMate[v];

  std::set<Size> roots;
  for (Size v = 0; v < mNumVtxs; ++v) {
    Size b = topLevel(v);
    if (b != cNullItm && b < mBlsm.size() && mBlsm[b].alive)
      roots.insert(b);
  }

  std::set<Size>::const_iterator rootIt;
  for (rootIt = roots.begin(); rootIt != roots.end(); ++rootIt) {
    Size b = *rootIt;
    if (trivial(b)) {
      Size matching = mBlsm[b].next;
      if (mMate[b] != cNullItm) {
        setExtractMate(b, mMate[b], &recursiveMate);
      } else if (matching != cNullItm && edgeSrc(matching) == b) {
        setExtractMate(b, edgeDst(matching), &recursiveMate);
      }
      continue;
    }
    Size matching = mBlsm[b].next;
    Size base = matching == cNullItm ? mBlsm[b].base : edgeSrc(matching);
    extractBlossomMatching(b, base, matching, &recursiveMate);
  }

  Size recCard = 0, dirCard = 0;
  Val recWgt = 0, dirWgt = 0;
  bool recValid = scoreMateVector(recursiveMate, &recCard, &recWgt);
  bool dirValid = scoreMateVector(directMate, &dirCard, &dirWgt);

  const std::vector<Size>* bestMate = &recursiveMate;
  *card = recCard;
  *edgWght = recWgt;

  if (!recValid || (dirValid && dirWgt > recWgt + EPS)) {
    bestMate = &directMate;
    *card = dirCard;
    *edgWght = dirWgt;
  }

  ResizeVector<Size>(mateVec, mNumVtxs);
  for (Size v = 0; v < mNumVtxs; ++v) (*mateVec)[v] = (*bestMate)[v];
}

inline bool FastMatchingEngineV2::scoreMateVector(
  const std::vector<Size>& mateVec, Size* card, Val* edgWght) const
{
  if (card == NULL || edgWght == NULL || mateVec.size() != mNumVtxs)
    return false;

  *card = 0;
  *edgWght = 0;
  bool valid = true;
  for (Size v = 0; v < mNumVtxs; ++v) {
    Size mate = mateVec[v];
    if (mate == cNullItm) continue;
    if (mate >= mNumVtxs || mateVec[mate] != v) {
      valid = false;
      continue;
    }
    if (mate > v) {
      bool found = false;
      for (Size i = 0; i < mAdj[v].size(); ++i) {
        if (mAdj[v][i] == mate) {
          ++(*card);
          *edgWght += mWgt[v][i];  // mWgt == original (cardTieEps=0, never mutated)
          found = true;
          break;
        }
      }
      if (!found) valid = false;
    }
  }
  return valid;
}

inline bool FastMatchingEngineV2::verifyCertificate(
  const std::vector<Size>& mateVec,
  Size* card,
  Val* edgWght,
  Val* dualObj,
  Val* minSlack) const
{
  if (card != NULL) *card = 0;
  if (edgWght != NULL) *edgWght = 0;
  if (dualObj != NULL) *dualObj = 0;
  if (minSlack != NULL) *minSlack = 0;
  if (mateVec.size() != mNumVtxs) return false;

  const Val tol = 1.0e-5;
  const Val objTol = 1.0e-4;
  bool valid = !mAbortRun;

  Size certCard = 0;
  Val certWeight = 0;
  if (!scoreMateVector(mateVec, &certCard, &certWeight)) valid = false;
  if (card != NULL) *card = certCard;
  if (edgWght != NULL) *edgWght = certWeight;

  std::vector<Val> yEff(mNumVtxs, 0);
  Val dual = 0;
  for (Size v = 0; v < mNumVtxs; ++v) {
    yEff[v] = effY(v);
    if (yEff[v] < -tol) valid = false;
    dual += yEff[v];
  }

  std::vector<Val> zEff(mBlsm.size(), 0);
  for (Size b = mNumVtxs; b < mBlsm.size(); ++b) {
    if (!mBlsm[b].alive || trivial(b)) continue;
    Val z = effZ(b);
    if (z < -tol) valid = false;
    if (z <= tol) continue;

    std::vector<Size> vtxs;
    collectBlossomVertices(b, &vtxs);
    std::sort(vtxs.begin(), vtxs.end());
    vtxs.erase(std::unique(vtxs.begin(), vtxs.end()), vtxs.end());
    if (vtxs.size() < 3 || (vtxs.size() % 2) == 0) {
      valid = false;
      continue;
    }

    zEff[b] = z;
    Size expectedInside = (vtxs.size() - 1) / 2;
    dual += static_cast<Val>(expectedInside) * z;

    Size internalMatched = 0;
    for (Size i = 0; i < vtxs.size(); ++i) {
      Size v = vtxs[i];
      Size mate = mateVec[v];
      if (mate != cNullItm && mate < mNumVtxs && v < mate &&
          std::binary_search(vtxs.begin(), vtxs.end(), mate))
        ++internalMatched;

      bool foundAncestor = false;
      Size cur = v;
      Size guard = 0;
      while (cur != cNullItm && cur < mBlsm.size() &&
             guard <= mBlsm.size()) {
        if (cur == b) {
          foundAncestor = true;
          break;
        }
        cur = mBlsm[cur].parent;
        ++guard;
      }
      if (!foundAncestor) valid = false;
    }
    if (internalMatched != expectedInside) valid = false;
  }

  Val edgeMinSlack = 0;
  bool sawEdge = false;
  for (Size u = 0; u < mNumVtxs; ++u) {
    for (Size i = 0; i < mAdj[u].size(); ++i) {
      Size v = mAdj[u][i];
      if (v >= mNumVtxs) {
        valid = false;
        continue;
      }
      Val slack = yEff[u] + yEff[v] + commonPositiveDual(u, v, zEff) -
                  mWgt[u][i];
      if (!sawEdge || slack < edgeMinSlack) {
        edgeMinSlack = slack;
        sawEdge = true;
      }
      if (slack < -tol) valid = false;
      if (mateVec[u] == v) {
        Val tight = slack < 0 ? -slack : slack;
        if (tight > objTol) valid = false;
      }
    }
  }
  if (minSlack != NULL) *minSlack = sawEdge ? edgeMinSlack : 0;
  if (dualObj != NULL) *dualObj = dual;

  Val diff = certWeight - dual;
  if (diff < 0) diff = -diff;
  if (diff > objTol) valid = false;
  return valid;
}

inline void FastMatchingEngineV2::collectPositiveDualAncestors(
  Size vtx,
  const std::vector<Val>& zEff,
  std::vector<std::pair<Size, Val> >* ancestors) const
{
  if (ancestors == NULL) return;
  ancestors->clear();
  if (vtx >= mNumVtxs) return;

  Size b = vtx;
  Size guard = 0;
  while (b != cNullItm && b < mBlsm.size() && guard <= mBlsm.size()) {
    if (!trivial(b) && mBlsm[b].alive &&
        b < zEff.size() && zEff[b] > 1.0e-5)
      ancestors->push_back(std::make_pair(b, zEff[b]));
    b = mBlsm[b].parent;
    ++guard;
  }
}

inline Val FastMatchingEngineV2::commonPositiveDual(
  Size u, Size v, const std::vector<Val>& zEff) const
{
  std::vector<std::pair<Size, Val> > au, av;
  collectPositiveDualAncestors(u, zEff, &au);
  collectPositiveDualAncestors(v, zEff, &av);

  Val sum = 0;
  for (Size i = 0; i < au.size(); ++i) {
    for (Size j = 0; j < av.size(); ++j) {
      if (au[i].first == av[j].first) {
        sum += au[i].second;
        break;
      }
    }
  }
  return sum;
}

#ifdef MBX_GABOW_EDGE_WATCH
inline void FastMatchingEngineV2::watchEdgeState(
  const char* where, Size eventCount, int eventKind)
{
  if (!mWatchEdgeEnabled || mWatchEdgeU >= mNumVtxs ||
      mWatchEdgeV >= mNumVtxs) return;
  const Size idx = findEdgeIdx(mWatchEdgeU, mWatchEdgeV);
  if (idx == cNullItm) return;
  Val commonZ = 0;
  for (Size a = mWatchEdgeU, ga = 0;
       a != cNullItm && a < mBlsm.size() && ga <= mBlsm.size();
       a = mBlsm[a].parent, ++ga) {
    if (trivial(a) || !mBlsm[a].alive) continue;
    for (Size b = mWatchEdgeV, gb = 0;
         b != cNullItm && b < mBlsm.size() && gb <= mBlsm.size();
         b = mBlsm[b].parent, ++gb) {
      if (a == b) {
        const Val z = effZ(a);
        if (z > 1.0e-5) commonZ += z;
        break;
      }
    }
  }
  const Val slack = effY(mWatchEdgeU) + effY(mWatchEdgeV) + commonZ -
                    mWgt[mWatchEdgeU][idx];
  if (!mWatchHaveSlack || std::fabs(slack - mWatchLastSlack) > EPS) {
    const Size bu = topLevel(mWatchEdgeU), bv = topLevel(mWatchEdgeV);
    Size hu = cNullItm, hv = cNullItm, gu = cNullItm, gv = cNullItm;
    Val cu = -std::numeric_limits<Val>::max();
    Val cv = -std::numeric_limits<Val>::max();
    bool liveU = false, liveV = false;
    int queued = 0;
    if (mWatchEdgeU < mBmHandleU.size() && idx < mBmHandleU[mWatchEdgeU].size()) {
      hu = mBmHandleU[mWatchEdgeU][idx]; hv = mBmHandleV[mWatchEdgeU][idx];
      gu = mBmGenU[mWatchEdgeU][idx]; gv = mBmGenV[mWatchEdgeU][idx];
      liveU = mBM.edge_cost_if(hu, gu, &cu);
      liveV = mBM.edge_cost_if(hv, gv, &cv);
      queued = mD3Queued[mWatchEdgeU][idx] != 0;
    }
    std::fprintf(stderr,
      "EDGEWATCH %s event=%lu kind=%d delta=%.17g slack=%.17g yu=%.17g "
      "yv=%.17g commonZ=%.17g d3key=%.17g queued=%d "
      "hU=%lu/%lu/live%d/c%.17g hV=%lu/%lu/live%d/c%.17g "
      "top=%lu/%lu status=%d/%d\n",
      where, (unsigned long)eventCount, eventKind, (double)mDeltaSum,
      (double)slack, (double)effY(mWatchEdgeU), (double)effY(mWatchEdgeV),
      (double)commonZ, (double)d3Key(mWatchEdgeU, idx), queued,
      (unsigned long)hu, (unsigned long)gu, liveU ? 1 : 0, (double)cu,
      (unsigned long)hv, (unsigned long)gv, liveV ? 1 : 0, (double)cv,
      (unsigned long)bu, (unsigned long)bv,
      bu < mBlsm.size() ? (int)mBlsm[bu].status : -1,
      bv < mBlsm.size() ? (int)mBlsm[bv].status : -1);
    if (slack <= 0.5) {
      mBM.debug_edge_state(stderr, "U", hu, gu);
      mBM.debug_edge_state(stderr, "V", hv, gv);
    }
    mWatchLastSlack = slack;
    mWatchHaveSlack = true;
  }
}
#endif

inline bool FastMatchingEngineV2::activeTree(Size tree) const
{
  if (tree == cNullItm || tree >= mNumVtxs) return false;
  Size B = topLevel(tree);
  if (B == cNullItm || B >= mBlsm.size()) return false;
  return mMate[tree] == cNullItm &&
         mBlsm[B].alive &&
         mBlsm[B].tree == tree &&
         mBlsm[B].status == BS_EVEN;
}

inline Size FastMatchingEngineV2::topLevel(Size vtx) const
{
#ifdef MBX_GABOW_SEARCH_AUDIT
  ++mAudTopLevelCalls;
#endif
  if (vtx >= mNumVtxs) return cNullItm;
  Size cached = mTopCacheOfVtx[vtx];
  if (cached != cNullItm && cached < mBlsm.size() &&
      mBlsm[cached].alive &&
      mBlsm[cached].parent == cNullItm &&
      mTopCacheGenerationOfVtx[vtx] == mBlsm[cached].generation) {
    return cached;
  }

  Size b = mTopBlsmOfVtx[vtx];
  if (b == cNullItm || b >= mBlsm.size() || !mBlsm[b].alive)
    return cNullItm;
  std::vector<Size> path;
  Size guard = 0;
  while (b != cNullItm && b < mBlsm.size() && guard <= mBlsm.size()) {
    if (!mBlsm[b].alive) return cNullItm;
    path.push_back(b);

    Size blossomCached = mBlsm[b].topCache;
    if (blossomCached != cNullItm &&
        blossomCached < mBlsm.size() &&
        mBlsm[blossomCached].alive &&
        mBlsm[blossomCached].parent == cNullItm &&
        mBlsm[b].topCacheGeneration == mBlsm[blossomCached].generation) {
      mTopCacheOfVtx[vtx] = blossomCached;
      mTopCacheGenerationOfVtx[vtx] = mBlsm[blossomCached].generation;
      for (Size i = 0; i < path.size(); ++i) {
        Size p = path[i];
        mBlsm[p].topCache = blossomCached;
        mBlsm[p].topCacheGeneration = mBlsm[blossomCached].generation;
      }
      return blossomCached;
    }

    if (mBlsm[b].parent == cNullItm) {
      mTopCacheOfVtx[vtx] = b;
      mTopCacheGenerationOfVtx[vtx] = mBlsm[b].generation;
      for (Size i = 0; i < path.size(); ++i) {
        Size p = path[i];
        mBlsm[p].topCache = b;
        mBlsm[p].topCacheGeneration = mBlsm[b].generation;
      }
      return b;
    }
    b = mBlsm[b].parent;
    ++guard;
  }
  return cNullItm;
}

inline void FastMatchingEngineV2::setExtractMate(
  Size u, Size v, std::vector<Size>* mateVec) const
{
  if (mateVec == NULL || u >= mNumVtxs) return;
  if (v >= mNumVtxs) {
    (*mateVec)[u] = cNullItm;
    return;
  }
  Size oldU = (*mateVec)[u];
  Size oldV = (*mateVec)[v];
  if (oldU != cNullItm && oldU < mNumVtxs && oldU != v &&
      (*mateVec)[oldU] == u) {
    (*mateVec)[oldU] = cNullItm;
  }
  if (oldV != cNullItm && oldV < mNumVtxs && oldV != u &&
      (*mateVec)[oldV] == v) {
    (*mateVec)[oldV] = cNullItm;
  }
  (*mateVec)[u] = v;
  (*mateVec)[v] = u;
}

inline void FastMatchingEngineV2::extractBlossomMatching(
  Size blossom, Size base, Size matching, std::vector<Size>* mateVec) const
{
  if (mateVec == NULL) return;
  // The nesting depth can approach n on a large final blossom. Preserve the
  // recursive visit order with explicit frames so extraction uses heap space
  // proportional to depth rather than exhausting the process stack.
  struct Frame {
    Size blossom, base, matching, k, ib, step, phase, edge;
    Frame(Size b, Size v, Size m)
      : blossom(b), base(v), matching(m), k(0), ib(0), step(0),
        phase(0), edge(cNullItm) {}
  };
  std::vector<Frame> frames;
  frames.push_back(Frame(blossom, base, matching));
  while (!frames.empty()) {
    Frame& f = frames.back();
    if (f.phase == 0) {
      if (f.blossom == cNullItm || f.blossom >= mBlsm.size()) {
        frames.pop_back();
        continue;
      }
      if (trivial(f.blossom)) {
        if (f.matching != cNullItm)
          setExtractMate(f.base, edgeDst(f.matching), mateVec);
        frames.pop_back();
        continue;
      }
      f.k = mBlsm[f.blossom].children.size();
      if (f.k == 0 || mBlsm[f.blossom].childEdges.size() != f.k) {
        frames.pop_back();
        continue;
      }
      f.ib = childIndexContaining(f.blossom, f.base);
      if (f.ib == cNullItm || f.ib >= f.k) f.ib = 0;
      f.step = 1;
      f.phase = 1;
      const Size child = mBlsm[f.blossom].children[f.ib];
      frames.push_back(Frame(child, f.base, f.matching));
      continue;
    }
    if (f.phase == 1) {
      if (f.step >= f.k) {
        frames.pop_back();
        continue;
      }
      const Size si = (f.ib + f.step) % f.k;
      const Size ti = (f.ib + f.step + 1) % f.k;
      const Size child = mBlsm[f.blossom].children[ti];
      f.edge = mBlsm[child].next;
      if (f.edge == cNullItm)
        f.edge = mBlsm[f.blossom].childEdges[si];
      if (f.edge == cNullItm) {
        f.step += 2;
        continue;
      }
      const Size first = mBlsm[f.blossom].children[si];
      const Size y = edgeDst(f.edge);
      const Size reverse = reverseEncodedEdge(f.edge);
      f.phase = 2;
      frames.push_back(Frame(first, y, reverse));
      continue;
    }
    if (f.phase == 2) {
      const Size ti = (f.ib + f.step + 1) % f.k;
      const Size second = mBlsm[f.blossom].children[ti];
      const Size x = edgeSrc(f.edge);
      const Size edge = f.edge;
      f.phase = 3;
      frames.push_back(Frame(second, x, edge));
      continue;
    }
    f.step += 2;
    f.phase = 1;
  }
}

inline void FastMatchingEngineV2::start()
{
  mProfWallStart = profNow();
  const Val INF = infKey();
  Size edgeRefs = 0;
  for (Size er = 0; er < mAdj.size(); ++er) edgeRefs += mAdj[er].size();
  Size eventCount = 0;
  Size heapScanCount = 0;
  Size eventCap;
  Size heapScanCap;
  Size capUnit = edgeRefs + mNumVtxs + 1000;
  Size capScale = (mNumVtxs == 0) ? Size(1) : mNumVtxs;
  Size capLimit = cNullItm / 80;
  if (capUnit != 0 && capScale > capLimit / capUnit) {
    eventCap = cNullItm;
    heapScanCap = cNullItm;
  } else {
    eventCap = 20 * capUnit * capScale + 1000;
    heapScanCap = 40 * capUnit * capScale + 1000;
  }

  long double mProfLoopLastDump = mProfWallStart;
  unsigned long long mProfLoopLastEventCount = 0;
  while (mUnmatched > 0 && !mAbortRun) {
    if (++eventCount > eventCap) {
      abortRun("event-cap");
      break;
    }
#ifdef MBX_GABOW_EDGE_WATCH
    watchEdgeState("before", eventCount, 0);
#endif

    // DIAG: time-based event-loop progress dump. Fires every ~10s of wall for
    // EITHER backend when GABOW_LOOPDIAG=1, so we can compare GMG vs Gabow's
    // event-loop breakdown (D124 scan vs extend/shrink/augment dispatch) and
    // see where the previously-"unaccounted" 75% actually goes.
    if (mGabowLoopDiag) {
      long double wnow = profNow();
      if (wnow - mProfLoopLastDump >= 10.0L) {
        unsigned long long bestByTreeEntries = 0;
        for (Size dv = 0; dv < mVtx.size(); ++dv)
          bestByTreeEntries += mVtx[dv].bestByTree.size();
        unsigned long long capChildren = 0, capChildEdges = 0, capSplitLeaves = 0;
        unsigned long long capVertices = 0, capChildOfVtxPairs = 0;
        unsigned long long liveCapVertices = 0, liveCapChildOfVtxPairs = 0;
        Size blossomMaxVerticesCapacity = 0;
        for (Size db = 0; db < mBlsm.size(); ++db) {
          const BlossomData& D = mBlsm[db];
          capChildren += D.children.capacity();
          capChildEdges += D.childEdges.capacity();
          capSplitLeaves += D.splitLeaves.capacity();
          capVertices += D.vertices.capacity();
          capChildOfVtxPairs += D.childOfVtx.capacity();
          if (D.alive) {
            liveCapVertices += D.vertices.capacity();
            liveCapChildOfVtxPairs += D.childOfVtx.capacity();
          }
          if (D.vertices.capacity() > blossomMaxVerticesCapacity)
            blossomMaxVerticesCapacity = D.vertices.capacity();
        }
        long double accounted = mProfPushD3Time + mProfMinEdgeTime +
                                mProfMergeTime + mProfDisposeTime +
                                mProfStaleCheckTime + mProfEvictTime +
                                mProfD124ScanTime + mProfActExtendTime +
                                mProfActShrinkTime + mProfActAugmentTime +
                                mProfActUnmatchTime;
        std::fprintf(stderr,
          "  [LOOP %s wall=%.1Lfs iters=%llu augments=%lu minEdgeCalls=%lu]\n"
          "      pushD3=%.2Lfs(mkEdge=%.2Lfs nca=%.2Lfs) minEdge=%.2Lfs "
          "evict=%.2Lfs staleChk=%.2Lfs merge=%.2Lfs\n"
          "      D124scan=%.2Lfs extend=%.2Lfs shrink=%.2Lfs augment=%.2Lfs "
          "unmatch=%.2Lfs\n"
          "      accounted=%.2Lfs unaccounted=%.2Lfs  maxTdepth=%lu"
          " treeMemberCap=%llu(%.1fMiB) maxRowCap=%lu\n"
          "      queues: D1=%lu D2=%lu D3=%lu D4=%lu CQpool=%lu"
          " bestByTree=%llu blossomCapsMiB: child=%.1f edges=%.1f split=%.1f"
          " vertices=%.1f(live %.1f) childMap=%.1f(live %.1f) maxVerticesCap=%lu\n",
          (mBackend == kGabow2018) ? "gabow" : "gmg",
          wnow - mProfWallStart, (unsigned long long)mProfLoopIters,
          (unsigned long)mProfAugments, (unsigned long)mProfMinEdgeCount,
          mProfPushD3Time, mProfMakeEdgeTime, mProfNcaTime, mProfMinEdgeTime,
          mProfEvictTime, mProfStaleCheckTime, mProfMergeTime,
          mProfD124ScanTime, mProfActExtendTime, mProfActShrinkTime,
          mProfActAugmentTime, mProfActUnmatchTime,
          accounted, (wnow - mProfWallStart) - accounted,
          (unsigned long)mAuxNCA.maxDepth(), mProfTreeMemberCapacity,
          (double)(mProfTreeMemberCapacity * sizeof(Size)) / (1024.0 * 1024.0),
          (unsigned long)mProfTreeMemberMaxRowCapacity,
          (unsigned long)mPQD1.size(), (unsigned long)mPQD2.size(),
          (unsigned long)mPQD3.size(), (unsigned long)mPQD4.size(),
          (unsigned long)mCQ.poolSize(), bestByTreeEntries,
          (double)(capChildren * sizeof(Size)) / (1024.0 * 1024.0),
          (double)(capChildEdges * sizeof(Size)) / (1024.0 * 1024.0),
          (double)(capSplitLeaves * sizeof(ConcatenableQueue::Handle)) / (1024.0 * 1024.0),
          (double)(capVertices * sizeof(Size)) / (1024.0 * 1024.0),
          (double)(liveCapVertices * sizeof(Size)) / (1024.0 * 1024.0),
          (double)(capChildOfVtxPairs * sizeof(std::pair<Size, Size>)) / (1024.0 * 1024.0),
          (double)(liveCapChildOfVtxPairs * sizeof(std::pair<Size, Size>)) / (1024.0 * 1024.0),
          (unsigned long)blossomMaxVerticesCapacity);
        std::fflush(stderr);
        mProfLoopLastDump = wnow;
        mProfLoopLastEventCount = eventCount;
      }
    }

    Val d1 = INF, d2 = INF, d3 = INF, d4 = INF;
#ifdef MBX_CLASSIC_EVENT_TRACE
    Val traceValidated[5] = {0,INF,INF,INF,INF};
    Val traceYu[5]={},traceYv[5]={},traceWeight[5]={};
    Size traceBlossom[5]={cNullItm,cNullItm,cNullItm,cNullItm,cNullItm};
    Size traceTree1[5]={cNullItm,cNullItm,cNullItm,cNullItm,cNullItm};
    Size traceTree2[5]={cNullItm,cNullItm,cNullItm,cNullItm,cNullItm};
#endif
    Size d1Vtx = cNullItm;
    Size d2Matched = cNullItm, d2EdgeIdx = cNullItm;
    Size d3U = cNullItm, d3EdgeIdx = cNullItm;
    BlossomMerging::EdgeHandle d3BmHandle = cNullItm;
    Size d3BmGeneration = cNullItm;
    Size d4Blsm = cNullItm;
    ++mProfLoopIters;
    long double profD124T0 = MBX_HOTNOW();

    while (!mPQD1.empty()) {
      if (++heapScanCount > heapScanCap) {
        abortRun("heap-d1-cap");
        return;
      }
      const HeapEntryD1& top = mPQD1.top();
      Size v = top.vtx;
      if (v >= mNumVtxs) {
        mPQD1.pop();
        continue;
      }
      if (top.version != mVtx[v].d1Version) {
        mPQD1.pop();
        continue;
      }
      Size B = topLevel(v);
      if (B == cNullItm || B >= mBlsm.size() ||
          mBlsm[B].status != BS_EVEN) {
        mPQD1.pop();
        mVtx[v].d1Queued = false;
        continue;
      }
      Val cur = mDeltaSum + effY(v);
      if (cur < top.key - EPS || cur > top.key + EPS) {
        mPQD1.pop();
        mVtx[v].d1Queued = false;
        pushD1(v);
        continue;
      }
      d1 = top.key;
#ifdef MBX_CLASSIC_EVENT_TRACE
      traceValidated[1] = cur;
      traceBlossom[1]=B;traceTree1[1]=mBlsm[B].tree;
#endif
      d1Vtx = v;
      break;
    }

    while (!mPQD2.empty()) {
      if (++heapScanCount > heapScanCap) {
        abortRun("heap-d2-cap");
        return;
      }
      // Copy before any pop(): PQD2 is backed by std::set, so pop erases the
      // node and invalidates a reference returned by top(). Several stale-entry
      // branches still need top.tree/matchedBlsm after removing that node.
      HeapEntryD2 top = mPQD2.top();
      Size v = top.matchedVtx;
      Size i = top.edgeIdxInMatched;
      if (v >= mNumVtxs || i >= mAdj[v].size()) {
        mPQD2.pop();
        continue;
      }
      Size u = mAdj[v][i];
      Size Bv = topLevel(v);
      Size Bu = topLevel(u);
      Size evenTree = (Bu != cNullItm && Bu < mBlsm.size()) ?
                      mBlsm[Bu].tree : cNullItm;
      if (top.matchedBlsm >= mBlsm.size() ||
          !mBlsm[top.matchedBlsm].d2Queued ||
          mBlsm[top.matchedBlsm].d2Vtx != top.matchedVtx ||
          mBlsm[top.matchedBlsm].d2Edge != top.edgeIdxInMatched ||
          mBlsm[top.matchedBlsm].d2Tree != top.tree ||
          top.key < mBlsm[top.matchedBlsm].d2Key - EPS ||
          top.key > mBlsm[top.matchedBlsm].d2Key + EPS) {
        mPQD2.pop();
        if (Bv != cNullItm && Bv < mBlsm.size() &&
            mBlsm[Bv].status == BS_MATCHED)
          pushD2BestForBlossom(Bv);
        continue;
      }
      if (Bv == cNullItm || Bv != top.matchedBlsm ||
          mBlsm[Bv].status != BS_MATCHED) {
        mPQD2.pop();
        eraseD2BlossomEntry(top.matchedBlsm);
        if (Bv != cNullItm && Bv < mBlsm.size() &&
            mBlsm[Bv].status == BS_MATCHED)
          pushD2BestForBlossom(Bv);
        continue;
      }
      if (Bu == cNullItm || Bv == Bu ||
          mBlsm[Bu].status != BS_EVEN ||
          !activeTree(evenTree)) {
        if (evenTree != cNullItm && evenTree < mNumVtxs &&
            !activeTree(evenTree)) destroyTree(evenTree);
        mPQD2.pop();
        eraseD2BlossomEntry(Bv);
        dropD2TreeCandidate(v, top.tree);
        continue;
      }
      if (top.tree != evenTree) {
        mPQD2.pop();
        eraseD2BlossomEntry(Bv);
        dropD2TreeCandidate(v, top.tree);
        addMatchedVertexD2Candidate(v, i);
        continue;
      }
      Val cur = d2Key(v, i);
      if (cur < top.key - EPS || cur > top.key + EPS) {
        mPQD2.pop();
        eraseD2BlossomEntry(Bv);
        dropD2TreeCandidate(v, top.tree);
        addMatchedVertexD2Candidate(v, i);
        continue;
      }
      d2 = top.key;
#ifdef MBX_CLASSIC_EVENT_TRACE
      traceValidated[2] = cur;
      traceBlossom[2]=Bv;traceTree1[2]=evenTree;
      if (mEventTrace.enabled()) { traceYu[2]=traceEffectiveY(u,Bu);traceYv[2]=traceEffectiveY(v,Bv);traceWeight[2]=mWgt[v][i]; }
#endif
      d2Matched = v;
      d2EdgeIdx = i;
      break;
    }
    mProfD124ScanTime += MBX_HOTNOW() - profD124T0;  // d1+d2 scan (d3 profiled separately)

    // D3 consumer. kGMG1986 uses mPQD3; kGabow2018 uses mBM.min_edge with
    // payload decoding to recover the original (u, edgeIdx).
    if (mBackend == kGabow2018) {
      // Gabow path: pop from mBM.
      while (true) {
        BlossomMerging::VertexId bmV = cNullItm, bmW = cNullItm;
        Val bmCost = 0;
        Size bmPayload = cNullItm;
        bool gotEdge;
        BlossomMerging::EdgeHandle bmHandle = cNullItm;
        { long double t0 = MBX_HOTNOW();
          gotEdge = mBM.min_edge(bmV, bmW, bmCost, &bmPayload, &bmHandle);
          ++mProfMinEdgeCount;
          mProfMinEdgeTime += MBX_HOTNOW() - t0; }
#ifndef MBX_NO_DIAGNOSTIC_OUTPUT
        if ((mProfMinEdgeCount % 100000ULL) == 0) {
          long double wallNow = profNow();
          long double wallElapsed = wallNow - mProfWallStart;
          // Cumulative time per bucket (seconds). pushD3Time already INCLUDES
          // make_edge+nca (they're called inside it), so we print pushD3 as a
          // whole plus its sub-parts, and separately the consumer-side buckets.
          long double accounted = mProfPushD3Time + mProfMinEdgeTime +
                                  mProfMergeTime + mProfDisposeTime +
                                  mProfStaleCheckTime + mProfEvictTime;
          std::fprintf(stderr,
            "  [minEdge=%lu wall=%.1Lfs]\n"
            "      pushD3=%.2Lfs (makeEdge=%.2Lfs nca=%.2Lfs)  minEdge=%.2Lfs\n"
            "      merge=%.2Lfs dispose=%.2Lfs staleChk=%.2Lfs evict=%.2Lfs\n"
            "      (staleChk split: destroyTree=%.2Lfs d2cand=%.2Lfs packetScan=%llu maxPktLen=%llu)\n"
            "      accounted=%.2Lfs  unaccounted=%.2Lfs  (stalePops=%lu merge=%lu)\n"
            "      staleWHY: sameBlsm=%lu odd=%lu matched=%lu inactiveTree=%lu other=%lu\n",
            (unsigned long)mProfMinEdgeCount, wallElapsed,
            mProfPushD3Time, mProfMakeEdgeTime, mProfNcaTime, mProfMinEdgeTime,
            mProfMergeTime, mProfDisposeTime, mProfStaleCheckTime, mProfEvictTime,
            mProfDestroyTreeTime, mProfD2CandTime, bmPacketScanCounter(),
            bmMaxPacketLen(),
            accounted, wallElapsed - accounted,
            (unsigned long)mProfStalePopCount, (unsigned long)mProfMergeCount,
            (unsigned long)mProfStaleSameBlsm, (unsigned long)mProfStaleOdd,
            (unsigned long)mProfStaleMatched, (unsigned long)mProfStaleInactiveTree,
            (unsigned long)mProfStaleOther);
          std::fflush(stderr);
        }
#endif
        if (!gotEdge) break;
        if (++heapScanCount > heapScanCap) {
          abortRun("heap-d3-cap"); return;
        }
        Size u = bmPayload / mAdjStride;
        Size i = bmPayload % mAdjStride;
        if (u >= mNumVtxs || i >= mAdj[u].size()) {
          // Payload malformed; defensively remove both handles and continue.
          if (u < mBmHandleU.size() && i < mBmHandleU[u].size()) {
            mBM.remove_edge(mBmHandleU[u][i]);
            mBM.remove_edge(mBmHandleV[u][i]);
            mBmHandleU[u][i] = cNullItm;
            mBmHandleV[u][i] = cNullItm;
            if (u < mD3Queued.size() && i < mD3Queued[u].size())
              mD3QueuedPopulation.set(mD3Queued, u, i, (char)0);
          }
          continue;
        }
        long double profStaleT0 = MBX_HOTNOW();
        Size v = mAdj[u][i];
        Size firedTag = 0;
        const Size currentTag = gabowForestOccurrenceTag(u, v);
        const Size firedGen = mBM.edge_generation(bmHandle);
        const bool firedTagLive =
          mBM.edge_occurrence_tag_if(bmHandle, firedGen, &firedTag);
        const bool occurrenceStale =
          mGabowDeepForest && (!firedTagLive || firedTag != currentTag);
        if (occurrenceStale) {
          ++mProfForestStaleFired;
          if (mGabowD3Log)
            std::fprintf(stderr,
              "FOREST stale-fired %lu %lu oldTag=%lu nowTag=%lu\n",
              (unsigned long)u, (unsigned long)v,
              (unsigned long)firedTag, (unsigned long)currentTag);
        }
        Size Bu = topLevel(u);
        Size Bv = topLevel(v);
        Size treeU = (Bu != cNullItm && Bu < mBlsm.size()) ?
                     mBlsm[Bu].tree : cNullItm;
        Size treeV = (Bv != cNullItm && Bv < mBlsm.size()) ?
                     mBlsm[Bv].tree : cNullItm;
        bool consumerStale =
          (Bu == cNullItm || Bv == cNullItm || Bu == Bv ||
           !(mBlsm[Bu].status == BS_EVEN &&
             mBlsm[Bv].status == BS_EVEN &&
             activeTree(treeU) && activeTree(treeV)));
        if (consumerStale) {
          ++mProfStalePopCount;
          // Categorize WHY this pop is stale (mutually-exclusive, first match):
          //  a) same blossom (endpoints merged)
          //  b) an endpoint is ODD (a T/inner blossom)  <- p.q.2 deactivation target
          //  c) an endpoint is MATCHED (left the tree)
          //  d) tree inactive
          //  e) other (cNullItm etc.)
          if (Bu == Bv) {
            ++mProfStaleSameBlsm;
          } else if ((Bu != cNullItm && mBlsm[Bu].status == BS_ODD) ||
                     (Bv != cNullItm && mBlsm[Bv].status == BS_ODD)) {
            ++mProfStaleOdd;
          } else if ((Bu != cNullItm && mBlsm[Bu].status == BS_MATCHED) ||
                     (Bv != cNullItm && mBlsm[Bv].status == BS_MATCHED)) {
            ++mProfStaleMatched;
          } else if ((treeU != cNullItm && !activeTree(treeU)) ||
                     (treeV != cNullItm && !activeTree(treeV))) {
            ++mProfStaleInactiveTree;
          } else {
            ++mProfStaleOther;
          }
          { long double dtT0 = MBX_HOTNOW();
          if (treeU != cNullItm && treeU < mNumVtxs &&
              !activeTree(treeU)) destroyTree(treeU);
          if (treeV != cNullItm && treeV < mNumVtxs &&
              !activeTree(treeV)) destroyTree(treeV);
            mProfDestroyTreeTime += MBX_HOTNOW() - dtT0; }
          { long double d2T0 = MBX_HOTNOW();
          if (Bu != cNullItm && Bv != cNullItm && Bu != Bv &&
              Bu < mBlsm.size() && Bv < mBlsm.size()) {
            if (mBlsm[Bu].status == BS_MATCHED &&
                mBlsm[Bv].status == BS_EVEN && activeTree(treeV)) {
              addMatchedVertexD2Candidate(u, i);
            } else if (mBlsm[Bu].status == BS_EVEN &&
                       mBlsm[Bv].status == BS_MATCHED &&
                       activeTree(treeU)) {
              Size rev = reverseEdgeIdx(u, i);
              if (rev != cNullItm) addMatchedVertexD2Candidate(v, rev);
            }
          }
            mProfD2CandTime += MBX_HOTNOW() - d2T0; }
          // The matchbox edge (u,v) is consumer-stale (an endpoint left the
          // EVEN-EVEN active state). BOTH its back-edges (u,a) and (v,a) are
          // therefore stale; remove both so neither resurfaces. Removing the
          // mate is safe here (unlike the cost-stale path, which re-pushes):
          // a consumer-stale edge will be freshly re-pushed by refresh if/when
          // its endpoints become EVEN-EVEN again.
          // destroyTree() above may have eagerly retired this same record.
          // Use the extraction-time generation so an already-freed (or ever
          // reused) slot is a guaranteed no-op.
          mBM.remove_edge_if(bmHandle, firedGen);
          if (u < mBmHandleU.size() && i < mBmHandleU[u].size()) {
            // Remove the mate back-edge too, guarded by generation so a stale
            // handle (slot freed/reused) is a safe no-op.
            if (mBmHandleU[u][i] != cNullItm && mBmHandleU[u][i] != bmHandle)
              mBM.remove_edge_if(mBmHandleU[u][i], mBmGenU[u][i]);
            if (mBmHandleV[u][i] != cNullItm && mBmHandleV[u][i] != bmHandle)
              mBM.remove_edge_if(mBmHandleV[u][i], mBmGenV[u][i]);
            mBmHandleU[u][i] = cNullItm;
            mBmHandleV[u][i] = cNullItm;
            if (u < mD3Queued.size() && i < mD3Queued[u].size())
              mD3QueuedPopulation.set(mD3Queued, u, i, (char)0);
          }
          mProfStaleCheckTime += MBX_HOTNOW() - profStaleT0;
          continue;
        }
        mProfStaleCheckTime += MBX_HOTNOW() - profStaleT0;
        Val cur = d3Key(u, i);
        if (cur < bmCost - EPS || cur > bmCost + EPS) {
          ++mProfStaleCostRepush;
          // Occurrence currency does not imply key currency: both back-edge
          // records were created with the same old slack. Retire the exact
          // returned record and its tracked mate, then rebuild once.
          mBM.remove_edge(bmHandle);
          if (u < mBmHandleU.size() && i < mBmHandleU[u].size()) {
            mBM.remove_edge_if(mBmHandleU[u][i], mBmGenU[u][i]);
            mBM.remove_edge_if(mBmHandleV[u][i], mBmGenV[u][i]);
            mBmHandleU[u][i] = cNullItm;
            mBmHandleV[u][i] = cNullItm;
            if (u < mD3Queued.size() && i < mD3Queued[u].size())
              mD3QueuedPopulation.set(mD3Queued, u, i, (char)0);
          }
          pushD3(u, i);
          continue;
        }
        if (occurrenceStale) {
          // Endpoint state and key are current, but the reduction still names
          // an old supporting-forest occurrence. Rebuild only this genuinely
          // live candidate. Status-stale records were already discarded above,
          // avoiding the former rebuild-before-status churn without allowing an
          // obsolete NCA/rank to drive a shrink or augmentation.
          ++mProfForestStaleRebuilt;
          bool removedReturned = false;
          if (u < mBmHandleU.size() && i < mBmHandleU[u].size()) {
            const bool removedU =
              mBM.remove_edge_if(mBmHandleU[u][i], mBmGenU[u][i]);
            const bool removedV =
              mBM.remove_edge_if(mBmHandleV[u][i], mBmGenV[u][i]);
            removedReturned =
              (mBmHandleU[u][i] == bmHandle && removedU) ||
              (mBmHandleV[u][i] == bmHandle && removedV);
            mBmHandleU[u][i] = cNullItm;
            mBmHandleV[u][i] = cNullItm;
          }
          if (!removedReturned) mBM.remove_edge(bmHandle);
          mD3QueuedPopulation.set(mD3Queued, u, i, (char)0);
          mProfStaleCheckTime += MBX_HOTNOW() - profStaleT0;
          pushD3(u, i);
          continue;
        }
        d3 = bmCost;
#ifdef MBX_CLASSIC_EVENT_TRACE
        traceValidated[3] = cur;
        traceTree1[3]=treeU;traceTree2[3]=treeV;
        if (mEventTrace.enabled()) { traceYu[3]=traceEffectiveY(u,Bu);traceYv[3]=traceEffectiveY(v,Bv);traceWeight[3]=mWgt[u][i]; }
#endif
        d3U = u;
        d3EdgeIdx = i;
        d3BmHandle = bmHandle;
        d3BmGeneration = firedGen;
        break;
      }
      goto d3_done;
    }
    while (!mPQD3.empty()) {
      if (++heapScanCount > heapScanCap) {
        abortRun("heap-d3-cap");
        return;
      }
      const HeapEntryD3& top = mPQD3.top();
      Size u = top.u;
      Size i = top.edgeIdx;
      if (u >= mNumVtxs || i >= mAdj[u].size()) {
        mPQD3.pop();
        continue;
      }
      Size v = mAdj[u][i];
#ifdef MBX_GMG_D3_OWNERSHIP_EXPERIMENT
      const bool ownsD3Slot = mD3Queued[u][i] != 0 &&
        top.publication == mGmgD3Publication[u][i].version;
#endif
      Size Bu = topLevel(u);
      Size Bv = topLevel(v);
      Size treeU = (Bu != cNullItm && Bu < mBlsm.size()) ?
                   mBlsm[Bu].tree : cNullItm;
      Size treeV = (Bv != cNullItm && Bv < mBlsm.size()) ?
                   mBlsm[Bv].tree : cNullItm;
      if (Bu == cNullItm || Bv == cNullItm || Bu == Bv ||
          !(mBlsm[Bu].status == BS_EVEN &&
            mBlsm[Bv].status == BS_EVEN &&
            activeTree(treeU) && activeTree(treeV))) {
        if (treeU != cNullItm && treeU < mNumVtxs &&
            !activeTree(treeU)) destroyTree(treeU);
        if (treeV != cNullItm && treeV < mNumVtxs &&
            !activeTree(treeV)) destroyTree(treeV);
        if (Bu != cNullItm && Bv != cNullItm && Bu != Bv &&
            Bu < mBlsm.size() && Bv < mBlsm.size()) {
          if (mBlsm[Bu].status == BS_MATCHED &&
              mBlsm[Bv].status == BS_EVEN && activeTree(treeV)) {
            addMatchedVertexD2Candidate(u, i);
          } else if (mBlsm[Bu].status == BS_EVEN &&
                     mBlsm[Bv].status == BS_MATCHED &&
                     activeTree(treeU)) {
            Size rev = reverseEdgeIdx(u, i);
            if (rev != cNullItm) addMatchedVertexD2Candidate(v, rev);
          }
        }
        mPQD3.pop();
#ifdef MBX_GABOW_SEARCH_AUDIT
        ++mAudGmgD3Stale;
#endif
#ifdef MBX_GMG_D3_OWNERSHIP_EXPERIMENT
        if (ownsD3Slot)
#endif
        if (u < mD3Queued.size() && i < mD3Queued[u].size())
          mD3Queued[u][i] = 0;
        continue;
      }
#ifdef MBX_GMG_D3_OWNERSHIP_EXPERIMENT
      if (!ownsD3Slot) {
        // Keep the state-stale cleanup/D2 discovery above. A superseded
        // price is never allowed to fire or clear a newer publication.
        mPQD3.pop();
#ifdef MBX_GABOW_SEARCH_AUDIT
        ++mAudGmgD3Stale;
#endif
        // Refresh if the live edge has repriced again; unchanged current
        // publications are suppressed by the publisher, including lower keys.
        pushD3(u, i);
        continue;
      }
#endif
      Val cur = d3Key(u, i);
      if (cur < top.key - EPS || cur > top.key + EPS) {
        mPQD3.pop();
#ifdef MBX_GABOW_SEARCH_AUDIT
        ++mAudGmgD3Stale;
#endif
        if (u < mD3Queued.size() && i < mD3Queued[u].size())
          mD3Queued[u][i] = 0;
        pushD3(u, i);
        continue;
      }
      d3 = top.key;
#ifdef MBX_CLASSIC_EVENT_TRACE
      traceValidated[3] = cur;
      traceTree1[3]=treeU;traceTree2[3]=treeV;
      if (mEventTrace.enabled()) { traceYu[3]=traceEffectiveY(u,Bu);traceYv[3]=traceEffectiveY(v,Bv);traceWeight[3]=mWgt[u][i]; }
#endif
      d3U = u;
      d3EdgeIdx = i;
      break;
    }
    d3_done: ;

    long double profD4T0 = MBX_HOTNOW();
    while (!mPQD4.empty()) {
      if (++heapScanCount > heapScanCap) {
        abortRun("heap-d4-cap");
        return;
      }
      const HeapEntryD4& top = mPQD4.top();
      Size B = top.blsm;
      if (B >= mBlsm.size() || !mBlsm[B].alive ||
          mBlsm[B].parent != cNullItm ||
          mBlsm[B].status != BS_ODD) {
        mPQD4.pop();
        continue;
      }
      Val cur = mDeltaSum + effZ(B) / 2.0;
      if (cur < top.key - EPS || cur > top.key + EPS) {
        mPQD4.pop();
        continue;
      }
      d4 = top.key;
#ifdef MBX_CLASSIC_EVENT_TRACE
      traceValidated[4] = cur;
      traceBlossom[4]=B;traceTree1[4]=mBlsm[B].tree;
#endif
      d4Blsm = B;
      break;
    }
    mProfD124ScanTime += MBX_HOTNOW() - profD4T0;  // + d4 scan

    Val mn = INF;
    int which = 0;
    if (d1 < mn - EPS) { mn = d1; which = 1; }
    if (d2 < mn - EPS) { mn = d2; which = 2; }
    if (d3 < mn - EPS) { mn = d3; which = 3; }
    if (d4 < mn - EPS) { mn = d4; which = 4; }


    if (mn >= INF - EPS) break;
#ifdef MBX_CLASSIC_EVENT_TRACE
    const Val traceQueueKey = mn;
    const Val traceOldClock = mDeltaSum;
#endif
    if (mn < mDeltaSum) mn = mDeltaSum;
    mDeltaSum = mn;
#ifdef MBX_CLASSIC_EVENT_TRACE
    if (mEventTrace.enabled()) {
      Size traceU=cNullItm,traceV=cNullItm;
      if (which==1) traceU=d1Vtx;
      else if (which==2) { traceU=mAdj[d2Matched][d2EdgeIdx];traceV=d2Matched; }
      else if (which==3) { traceU=d3U;traceV=mAdj[d3U][d3EdgeIdx]; }
      else if (which==4) traceU=mBlsm[d4Blsm].base;
      mEventTrace.selected(which,traceQueueKey,traceValidated[which],traceOldClock,mDeltaSum,traceU,traceV,
        traceYu[which],traceYv[which],traceWeight[which]);
      traceAction(ClassicEventTrace::Dispatch,traceU,traceV,traceBlossom[which],cNullItm,cNullItm,traceTree1[which],traceTree2[which]);
    }
#endif

    if (which == 1) {
      long double t0 = MBX_HOTNOW();
#ifdef MBX_GABOW_SEARCH_AUDIT
      ++mAudEvD1;
#endif
      mPQD1.pop();
      if (d1Vtx < mNumVtxs) mVtx[d1Vtx].d1Queued = false;
      unmatchNode(d1Vtx);
      --mUnmatched;
      mProfActUnmatchTime += MBX_HOTNOW() - t0;
    } else if (which == 2) {
#ifdef MBX_GABOW_SEARCH_AUDIT
      ++mAudEvD2;
#endif
      mPQD2.pop();
      Size matched = d2Matched;
      Size matchedBEntry = topLevel(matched);
      if (matchedBEntry != cNullItm && matchedBEntry < mBlsm.size())
        eraseD2BlossomEntry(matchedBEntry);
      Size even = mAdj[matched][d2EdgeIdx];
      Size idxEven = findEdgeIdx(even, matched);
      Size matchedB = topLevel(matched);
      bool hasMatchedNext = (matchedB != cNullItm && matchedB < mBlsm.size() &&
                             mBlsm[matchedB].next != cNullItm);
      if (!hasMatchedNext && matchedB != cNullItm && matchedB < mBlsm.size() &&
          trivial(matchedB) && mMate[matched] != cNullItm) {
        hasMatchedNext = true;
      }
      bool d2Augment = (matchedB != cNullItm && matchedB < mBlsm.size() &&
                        !hasMatchedNext);
      if (d2Augment) {
        long double t0 = MBX_HOTNOW();
        augmentOnArc(even, matched, idxEven);
        mProfActAugmentTime += MBX_HOTNOW() - t0;
      }
      else {
        long double t0 = MBX_HOTNOW();
#ifdef MBX_GABOW_SEARCH_AUDIT
        ++mAudGrows;
#endif
        extendOnArc(even, matched, idxEven);
        mProfActExtendTime += MBX_HOTNOW() - t0;
      }
    } else if (which == 3) {
#ifdef MBX_GABOW_SEARCH_AUDIT
      ++mAudEvD3;
#endif
      Size u = d3U;
      Size edgeIdx = d3EdgeIdx;
      if (mBackend == kGabow2018) {
        // Remove the exact returned record plus its tracked mate. Generation
        // guards prevent a recycled handle from deleting another edge.
        mBM.remove_edge_if(d3BmHandle, d3BmGeneration);
        if (u < mBmHandleU.size() && edgeIdx < mBmHandleU[u].size()) {
          mBM.remove_edge_if(mBmHandleU[u][edgeIdx], mBmGenU[u][edgeIdx]);
          mBM.remove_edge_if(mBmHandleV[u][edgeIdx], mBmGenV[u][edgeIdx]);
          mBmHandleU[u][edgeIdx] = cNullItm;
          mBmHandleV[u][edgeIdx] = cNullItm;
        }
      } else {
        mPQD3.pop();
      }
      if (u < mD3Queued.size() && edgeIdx < mD3Queued[u].size())
        setD3QueuedFlag(u, edgeIdx, (char)0);
      Size v = mAdj[u][edgeIdx];
      Size Bu = topLevel(u);
      Size Bv = topLevel(v);
      if (mGabowD3Log)
        std::fprintf(stderr, "F %lu %lu k=%.0f %s\n",
          (unsigned long)(u < v ? u : v), (unsigned long)(u < v ? v : u),
          (double)d3, (mBlsm[Bu].tree != mBlsm[Bv].tree) ? "aug" : "shr");
      if (mBlsm[Bu].tree != mBlsm[Bv].tree) {
        long double t0 = MBX_HOTNOW();
        augmentOnEdge(u, v, edgeIdx);
        mProfActAugmentTime += MBX_HOTNOW() - t0;
      } else {
        long double t0 = MBX_HOTNOW();
        shrinkOnEdge(u, v, edgeIdx, mBlsm[Bu].tree);
        mProfActShrinkTime += MBX_HOTNOW() - t0;
      }
    } else {
      long double t0 = MBX_HOTNOW();
#ifdef MBX_GABOW_SEARCH_AUDIT
      ++mAudEvD4;
#endif
      mPQD4.pop();
      splitBlossom(d4Blsm);
      mProfActShrinkTime += MBX_HOTNOW() - t0;
    }
#ifdef MBX_MONOTONE_SUPPORT_VERIFY
    if(mBackend==kGabow2018&&(mNumVtxs<=32||eventCount<=3||eventCount%16384==0))
      supportVerifyStable("event");
#endif
#ifdef MBX_GABOW_EDGE_WATCH
    watchEdgeState("after", eventCount, which);
#endif
  }

}

inline void FastMatchingEngineV2::setGabowBlossomActive(Size blossom, bool active)
{
  if (mBackend != kGabow2018 || blossom >= mBlsm.size()) return;
  // A reduced back-edge is owned by the blossom containing its Table-1
  // charged endpoint, which can be an internal NCA rather than either original
  // endpoint.  Therefore a matchbox label transition for that owner does NOT
  // imply that all original edges stored there are inactive.  Deactivating the
  // BM header hid live EVEN-EVEN candidates from min_edge (and allowed the dual
  // search to advance past their tight time).  Keep headers visible and rely on
  // the D3 consumer's original-endpoint status/tree validation to discard stale
  // records safely.
  if (!active) return;
  Size atom = cNullItm;
  if (trivial(blossom)) atom = blossom;
  else if (!mBlsm[blossom].vertices.empty()) atom = mBlsm[blossom].vertices[0];
  if (atom != cNullItm && atom < mNumVtxs)
    mBM.set_vertex_blossom_active(atom, active);
}

inline void FastMatchingEngineV2::retireGabowForestVertex(Size v)
{
  if (mBackend != kGabow2018 || !mGabowDeepForest ||
      v >= mNumVtxs || v >= mInT.size() || !mInT[v]) return;
  mInT[v] = 0;
#ifdef MBX_GABOW_SEARCH_AUDIT
  ++mAudRetires;
#endif
  if (v < mForestGeneration.size()) {
    ++mForestGeneration[v];
    if (mForestGeneration[v] == 0) ++mForestGeneration[v];
  }
}

inline void FastMatchingEngineV2::retireGabowD3ForMatchedBlossom(Size blossom)
{
  if (mBackend != kGabow2018 || blossom >= mBlsm.size()) return;

  ++mProfEvictCalls;

  std::vector<Size> collected;
  MembershipView vertices;
  const bool atomic = trivial(blossom);
  if (!atomic && !mBlsm[blossom].vertices.empty()) {
    vertices = &mBlsm[blossom].vertices;
  } else if (!atomic) {
    collectBlossomVertices(blossom, &collected);
    vertices = &collected;
  }

  // D3 bookkeeping is stored only at the canonical (lower-numbered)
  // endpoint.  Scan just the original edges incident to this blossom, map
  // each slot to that canonical record, and retire its two exact BM records.
  // Generation checks make duplicate visits (including internal edges) safe.
  // No BM blossom/header is disabled: unrelated records owned by the same BM
  // blossom remain visible to min_edge().
  const Size numVertices = atomic ? 1 : vertices->size();
  for (Size vi = 0; vi < numVertices; ++vi) {
    const Size v = atomic ? blossom : (*vertices)[vi];
    if (v >= mNumVtxs) continue;
    for (Size i = 0; i < mAdj[v].size(); ++i) {
      ++mProfEvictSlotsScanned;
      const Size w = mAdj[v][i];
      Size u = v;
      Size edgeIdx = i;
      if (w < v) {
        u = w;
        edgeIdx = reverseEdgeIdx(v, i);
        if (edgeIdx == cNullItm) continue;
      }
      if (u >= mD3Queued.size() || edgeIdx >= mD3Queued[u].size() ||
          mD3Queued[u][edgeIdx] == 0)
        continue;

      if (u < mBmHandleU.size() && edgeIdx < mBmHandleU[u].size()) {
        if (mBM.remove_edge_if_batched(
              mBmHandleU[u][edgeIdx], mBmGenU[u][edgeIdx]))
          ++mProfEvictLiveRemovals;
        if (mBM.remove_edge_if_batched(
              mBmHandleV[u][edgeIdx], mBmGenV[u][edgeIdx]))
          ++mProfEvictLiveRemovals;
        mBmHandleU[u][edgeIdx] = cNullItm;
        mBmHandleV[u][edgeIdx] = cNullItm;
      }
      mD3QueuedPopulation.set(mD3Queued, u, edgeIdx, (char)0);
    }
  }
}

inline Size FastMatchingEngineV2::gabowForestOccurrenceTag(Size u, Size v) const
{
  if (mBackend != kGabow2018 || !mGabowDeepForest) return 0;
  assert(sizeof(Size) >= 8);
  assert(u < mForestGeneration.size() && v < mForestGeneration.size());
  const unsigned long long hi = (unsigned long long)mForestGeneration[u];
  const unsigned long long lo = (unsigned long long)mForestGeneration[v];
  return (Size)((hi << 32) | lo);
}

inline void FastMatchingEngineV2::matchedToEven(Size blossom, Size tree)
{
  eraseD2BlossomEntry(blossom);
  setGabowBlossomActive(blossom, true);
  if (!trivial(blossom)) {
    mBlsm[blossom].status = BS_EVEN;
    mBlsm[blossom].tree = tree;
    mBlsm[blossom].offset = mDeltaSum;
    if(mBackend==kGabow2018)supportOuter(blossom);
    addTreeMember(blossom, tree);
    refreshEvenBlossomEdges(blossom);
    pushD1ForBlossom(blossom);
    return;
  }

  Size v = blossom;
  mBlsm[blossom].status = BS_EVEN;
  mBlsm[blossom].tree = tree;
  mBlsm[blossom].offset = mDeltaSum;
  if(mBackend==kGabow2018)supportOuter(blossom);
  addTreeMember(blossom, tree);
  clearVertexD2(v);

  pushD1(v);

  for (Size i = 0; i < mAdj[v].size(); ++i) {
    Size w = mAdj[v][i];
    Size Bw = topLevel(w);
    if (Bw == blossom) continue;
    if (mBlsm[Bw].status == BS_EVEN) {
      pushD3(v, i);
    }
    else if (mBlsm[Bw].status == BS_MATCHED) {
      Size rev = reverseEdgeIdx(v, i);
      if (rev != cNullItm) addMatchedVertexD2Candidate(w, rev);
    }
  }
}

inline void FastMatchingEngineV2::matchedToOdd(Size blossom)
{
  eraseD2BlossomEntry(blossom);
  setGabowBlossomActive(blossom, false);
  if (!trivial(blossom)) {
    mBlsm[blossom].status = BS_ODD;
    mBlsm[blossom].offset = mDeltaSum;
    if(mBackend==kGabow2018)supportInner(blossom);
    pushD4(blossom);
    return;
  }

  Size v = blossom;
  clearVertexD2(v);
  mBlsm[blossom].status = BS_ODD;
  mBlsm[blossom].offset = mDeltaSum;
  if(mBackend==kGabow2018)supportInner(blossom);
}

inline void FastMatchingEngineV2::evenToMatched(Size blossom, Size /*tree*/)
{
  setGabowBlossomActive(blossom, false);
  Val deltaSinceEven = mDeltaSum - mBlsm[blossom].offset;
  if (!trivial(blossom)) {
    std::vector<Size> vtxs;
    MembershipView pvtxs = &mBlsm[blossom].vertices;
    if (pvtxs->empty()) {
      collectBlossomVertices(blossom, &vtxs);
      pvtxs = &vtxs;
    }
    for (Size i = 0; i < pvtxs->size(); ++i) {
      Size v = (*pvtxs)[i];
      mVtx[v].pot -= deltaSinceEven;
      mVtx[v].lastEvenScanTree = cNullItm;
      mVtx[v].lastEvenScanTop = cNullItm;
      // R_Scope deep forest: leaving a live tree retires the T occurrence, so
      // a later grow re-installs at the REAL depth (sole caller: destroyTree).
      retireGabowForestVertex(v);
    }
    mBlsm[blossom].d1MinValid = false;
    mBlsm[blossom].pot += 2.0 * deltaSinceEven;
  } else {
    Size v = blossom;
    mVtx[v].pot -= deltaSinceEven;
    mVtx[v].lastEvenScanTree = cNullItm;
    mVtx[v].lastEvenScanTop = cNullItm;
    retireGabowForestVertex(v);
  }

  mBlsm[blossom].status = BS_MATCHED;
  mBlsm[blossom].tree = cNullItm;
  mBlsm[blossom].offset = 0;

  if (trivial(blossom)) {
    refreshMatchedVertexD2(blossom);
  } else {
    refreshMatchedBlossomD2(blossom);
  }
}

inline void FastMatchingEngineV2::oddToMatched(Size blossom)
{
  setGabowBlossomActive(blossom, false);
  Val deltaSinceOdd = mDeltaSum - mBlsm[blossom].offset;
  if (!trivial(blossom)) {
    std::vector<Size> vtxs;
    MembershipView pvtxs = &mBlsm[blossom].vertices;
    if (pvtxs->empty()) {
      collectBlossomVertices(blossom, &vtxs);
      pvtxs = &vtxs;
    }
    for (Size i = 0; i < pvtxs->size(); ++i) {
      Size v = (*pvtxs)[i];
      mVtx[v].pot += deltaSinceOdd;
      mVtx[v].lastEvenScanTree = cNullItm;
      mVtx[v].lastEvenScanTop = cNullItm;
      // R_Scope deep forest: retire the T occurrence (sole caller: destroyTree).
      retireGabowForestVertex(v);
    }
    mBlsm[blossom].d1MinValid = false;
    mBlsm[blossom].pot -= 2.0 * deltaSinceOdd;
  } else {
    Size v = blossom;
    mVtx[v].pot += deltaSinceOdd;
    mVtx[v].lastEvenScanTree = cNullItm;
    mVtx[v].lastEvenScanTop = cNullItm;
    retireGabowForestVertex(v);
  }

  mBlsm[blossom].status = BS_MATCHED;
  mBlsm[blossom].tree = cNullItm;
  mBlsm[blossom].offset = 0;

  if (trivial(blossom)) refreshMatchedVertexD2(blossom);
  else {
    refreshMatchedBlossomD2(blossom);
  }
}

inline void FastMatchingEngineV2::oddToEven(Size blossom, Size tree)
{
  Val deltaSinceOdd = mDeltaSum - mBlsm[blossom].offset;
  if (!trivial(blossom)) {
    std::vector<Size> vtxs;
    MembershipView pvtxs = &mBlsm[blossom].vertices;
    if (pvtxs->empty()) {
      collectBlossomVertices(blossom, &vtxs);
      pvtxs = &vtxs;
    }
    for (Size i = 0; i < pvtxs->size(); ++i) mVtx[(*pvtxs)[i]].pot += deltaSinceOdd;
    mBlsm[blossom].d1MinValid = false;
    mBlsm[blossom].pot -= 2.0 * deltaSinceOdd;
  } else {
    mVtx[blossom].pot += deltaSinceOdd;
  }
  matchedToEven(blossom, tree);
}

inline void FastMatchingEngineV2::bakeEvenToMatchedForShrink(Size blossom)
{
  Val deltaSinceEven = mDeltaSum - mBlsm[blossom].offset;
  if (!trivial(blossom)) {
    std::vector<Size> vtxs;
    MembershipView pvtxs = &mBlsm[blossom].vertices;
    if (pvtxs->empty()) {
      collectBlossomVertices(blossom, &vtxs);
      pvtxs = &vtxs;
    }
    for (Size i = 0; i < pvtxs->size(); ++i) {
      Size v = (*pvtxs)[i];
      mVtx[v].pot -= deltaSinceEven;
      mVtx[v].lastEvenScanTree = cNullItm;
      mVtx[v].lastEvenScanTop = cNullItm;
    }
    mBlsm[blossom].d1MinValid = false;
    mBlsm[blossom].pot += 2.0 * deltaSinceEven;
  } else {
    mVtx[blossom].pot -= deltaSinceEven;
    mVtx[blossom].lastEvenScanTree = cNullItm;
    mVtx[blossom].lastEvenScanTop = cNullItm;
  }

  mBlsm[blossom].status = BS_MATCHED;
  mBlsm[blossom].tree = cNullItm;
  mBlsm[blossom].offset = 0;
}

inline void FastMatchingEngineV2::bakeOddToMatchedForShrink(Size blossom)
{
  Val deltaSinceOdd = mDeltaSum - mBlsm[blossom].offset;
  if (!trivial(blossom)) {
    std::vector<Size> vtxs;
    MembershipView pvtxs = &mBlsm[blossom].vertices;
    if (pvtxs->empty()) {
      collectBlossomVertices(blossom, &vtxs);
      pvtxs = &vtxs;
    }
    for (Size i = 0; i < pvtxs->size(); ++i) {
      Size v = (*pvtxs)[i];
      mVtx[v].pot += deltaSinceOdd;
      mVtx[v].lastEvenScanTree = cNullItm;
      mVtx[v].lastEvenScanTop = cNullItm;
    }
    mBlsm[blossom].d1MinValid = false;
    mBlsm[blossom].pot -= 2.0 * deltaSinceOdd;
  } else {
    mVtx[blossom].pot += deltaSinceOdd;
    mVtx[blossom].lastEvenScanTree = cNullItm;
    mVtx[blossom].lastEvenScanTop = cNullItm;
  }

  mBlsm[blossom].status = BS_MATCHED;
  mBlsm[blossom].tree = cNullItm;
  mBlsm[blossom].offset = 0;
}

inline void FastMatchingEngineV2::alternatePath(Size even, Size /*tree*/)
{
  Size entry = trivial(even) ? even : mBlsm[even].base;
  flipPathFromTopEven(even, entry);
}

inline void FastMatchingEngineV2::destroyTree(Size tree)
{
  if (tree >= mTreeMembers.size()) return;

  std::vector<Size> members = mTreeMembers[tree];
  if (mBackend == kGabow2018) {
    const long double t0 = MBX_HOTNOW();
    mBM.begin_edge_removal_batch();
    for (Size i = 0; i < members.size(); ++i) {
      const Size b = members[i];
      if (b >= mBlsm.size() || !mBlsm[b].alive ||
          mBlsm[b].parent != cNullItm || mBlsm[b].tree != tree)
        continue;
      if (mBlsm[b].status == BS_EVEN || mBlsm[b].status == BS_ODD)
        retireGabowD3ForMatchedBlossom(b);
    }
    mBM.end_edge_removal_batch();
    mProfEvictTime += MBX_HOTNOW() - t0;
  }
  for (Size i = 0; i < members.size(); ++i) {
    Size b = members[i];
    if (b < mBlsm.size() && mBlsm[b].listedTree == tree)
      mBlsm[b].listedTree = cNullItm;
    if (b >= mBlsm.size() || !mBlsm[b].alive || mBlsm[b].parent != cNullItm)
      continue;
    if (mBlsm[b].tree != tree) continue;

    if (mBlsm[b].status == BS_EVEN) evenToMatched(b, tree);
    else if (mBlsm[b].status == BS_ODD) oddToMatched(b);
    else mBlsm[b].tree = cNullItm;

    mBlsm[b].pred = cNullItm;
    if (trivial(b)) mTreeParent[b] = cNullItm;
  }

  mTreeMembers[tree].clear();
  if(mBackend==kGabow2018)supportRetireTree(tree);
}

inline void FastMatchingEngineV2::unmatchNode(Size vtx)
{
  Size B = topLevel(vtx);
  Size tree = mBlsm[B].tree;
  flipPathFromTopEven(B, vtx);
  destroyTree(tree);
  B = topLevel(vtx);
  clearMate(vtx);
  mBlsm[B].base = vtx;
  mBlsm[B].next = cNullItm;
#ifdef MBX_CLASSIC_EVENT_TRACE
  traceAction(ClassicEventTrace::Settle,vtx,cNullItm,B,cNullItm,cNullItm,tree);
#endif
}

inline bool FastMatchingEngineV2::collectPathChildren(Size blossom, Size entryVtx,
    std::vector<std::pair<Size,Size> >& out) const
{
  out.clear();
  if (blossom == cNullItm || blossom >= mBlsm.size() || trivial(blossom)) return false;
  Size entry = childIndexContaining(blossom, entryVtx);
  Size k = mBlsm[blossom].children.size();
  if (entry == cNullItm || k == 0 || mBlsm[blossom].childEdges.size() != k) return false;
  Size baseIdx = childIndexContaining(blossom, mBlsm[blossom].base);
  if (baseIdx == cNullItm || baseIdx >= k) baseIdx = 0;
  out.push_back(std::make_pair(mBlsm[blossom].children[entry], entryVtx));
  Size cur = entry, guard = 0;
  while (cur != baseIdx) {
    if (++guard > k) { out.clear(); return false; }
    Size rel = (cur + k - baseIdx) % k;
    if ((rel % 2) == 0) {
      Size byIdx = (cur + k - 1) % k;
      Size bxIdx = (cur + k - 2) % k;
      Size edge = mBlsm[blossom].childEdges[bxIdx];
      if (edge == cNullItm) { out.clear(); return false; }
      Size x = edgeSrc(edge), y = edgeDst(edge);
      if (x >= mNumVtxs || y >= mNumVtxs) { out.clear(); return false; }
      out.push_back(std::make_pair(mBlsm[blossom].children[byIdx], y));
      out.push_back(std::make_pair(mBlsm[blossom].children[bxIdx], x));
      cur = bxIdx;
    } else {
      Size edgeIdx = (cur + 1) % k;
      Size bxIdx = edgeIdx;
      Size byIdx = (edgeIdx + 1) % k;
      Size edge = mBlsm[blossom].childEdges[edgeIdx];
      if (edge == cNullItm) { out.clear(); return false; }
      Size x = edgeSrc(edge), y = edgeDst(edge);
      if (x >= mNumVtxs || y >= mNumVtxs) { out.clear(); return false; }
      out.push_back(std::make_pair(mBlsm[blossom].children[bxIdx], x));
      out.push_back(std::make_pair(mBlsm[blossom].children[byIdx], y));
      cur = byIdx;
    }
  }
  return true;
}

inline Size FastMatchingEngineV2::installBlossomChainInForest(Size blossom,
    Size entryVtx,Size anchor)
{
  return supportInstallPath(blossom,entryVtx,anchor);
}

inline void FastMatchingEngineV2::extendOnArc(Size u, Size v, Size edgeIdxInU)
{
  Size Bu = topLevel(u);
  Size Bv = topLevel(v);
  if (Bu == cNullItm || Bv == cNullItm) return;
  if (mBlsm[Bu].status != BS_EVEN || mBlsm[Bv].status != BS_MATCHED) return;

  Size next = mBlsm[Bv].next;
  if (trivial(Bv) && v == Bv && mMate[v] != cNullItm) {
    next = encodeEdge(v, mMate[v], mNumVtxs);
    mBlsm[Bv].next = next;
  }
  if (next == cNullItm) {
    augmentOnArc(u, v, cNullItm);
    return;
  }

  Size baseV = edgeSrc(next);
  Size mv = edgeDst(next);
  if (baseV >= mNumVtxs || mv >= mNumVtxs ||
      topLevel(baseV) != Bv) {
    if (mv < mNumVtxs && topLevel(mv) == Bv &&
        baseV < mNumVtxs) {
      next = reverseEncodedEdge(next);
      mBlsm[Bv].next = next;
      baseV = edgeSrc(next);
      mv = edgeDst(next);
    } else {
      abortRun("extend-next-not-in-blossom");
      return;
    }
  }
  Size Bm = topLevel(mv);
  if (Bm != cNullItm && Bm != Bv && Bm < mBlsm.size() &&
      mBlsm[Bm].status == BS_EVEN &&
      mBlsm[Bm].tree == mBlsm[Bu].tree) {
    Size tree = mBlsm[Bm].tree;
    mBlsm[Bv].tree = tree;
    mBlsm[Bv].pred = reverseEncodedEdge(next);
    mBlsm[Bv].next = next;
    addTreeMember(Bv, tree);
    matchedToOdd(Bv);
#ifdef MBX_CLASSIC_EVENT_TRACE
    traceAction(ClassicEventTrace::Grow,u,v,Bv,baseV,mv,tree,cNullItm,1);
#endif
    return;
  }
  if (Bm != cNullItm && Bm != Bv && Bm < mBlsm.size() &&
      mBlsm[Bm].status == BS_EVEN &&
      mBlsm[Bm].tree != mBlsm[Bu].tree) {
    Size treeM = mBlsm[Bm].tree;
    mBlsm[Bv].tree = treeM;
    mBlsm[Bv].pred = reverseEncodedEdge(next);
    mBlsm[Bv].next = next;
    addTreeMember(Bv, treeM);
    matchedToOdd(Bv);
#ifdef MBX_CLASSIC_EVENT_TRACE
    traceAction(ClassicEventTrace::Grow,u,v,Bv,baseV,mv,treeM,cNullItm,2);
#endif
    return;
  }
  if (Bm != cNullItm && Bm != Bv && Bm < mBlsm.size() &&
      mBlsm[Bm].status == BS_ODD &&
      mBlsm[Bm].tree == mBlsm[Bu].tree) {
    Size tree = mBlsm[Bm].tree;
    mBlsm[Bv].tree = tree;
    mBlsm[Bv].pred = reverseEncodedEdge(next);
    mBlsm[Bv].next = next;
    matchedToEven(Bv, tree);
#ifdef MBX_CLASSIC_EVENT_TRACE
    traceAction(ClassicEventTrace::Grow,u,v,Bv,baseV,mv,tree,cNullItm,3);
#endif
    Size idx = edgeIdxInU;
    if (idx == cNullItm) idx = findEdgeIdx(u, v);
    if (idx != cNullItm) shrinkOnEdge(u, v, idx, tree);
    return;
  }
  if (Bm == cNullItm || Bm == Bv || mBlsm[Bm].status != BS_MATCHED) {
    abortRun("extend-mate-blossom-invalid");
    return;
  }

  Size tree = mBlsm[Bu].tree;
  mTreeParent[v] = u;

  // Supporting-forest admission occurs in matchedToOdd/Even below.
  mBlsm[Bv].tree = tree;
  mBlsm[Bv].pred = encodeEdge(u, v, mNumVtxs);
  mBlsm[Bv].next = next;
  addTreeMember(Bv, tree);
  matchedToOdd(Bv);

  mTreeParent[mv] = baseV;
  mBlsm[Bm].pred = next;
  matchedToEven(Bm, tree);
#ifdef MBX_CLASSIC_EVENT_TRACE
  traceAction(ClassicEventTrace::Grow,u,v,Bv,baseV,mv,tree,cNullItm,0);
#endif
}

inline void FastMatchingEngineV2::augmentOnArc(Size u, Size v, Size /*edgeIdxInU*/)
{
  Size Bu = topLevel(u);
  Size Bv = topLevel(v);
  if (Bu == cNullItm || Bv == cNullItm) return;
  if (mBlsm[Bu].status != BS_EVEN || mBlsm[Bv].status != BS_MATCHED) return;

  Size tree = mBlsm[Bu].tree;
  flipPathFromTopEven(Bu, u);
  augmentThroughBlossom(Bv, v);
  setMatePair(u, v);
  Bu = topLevel(u);
  Bv = topLevel(v);
  if (Bu != cNullItm && Bu < mBlsm.size())
    mBlsm[Bu].next = encodeEdge(u, v, mNumVtxs);
  if (Bv != cNullItm && Bv < mBlsm.size())
    mBlsm[Bv].next = encodeEdge(v, u, mNumVtxs);
  destroyTree(tree);
  refreshMatchedVertexD2(v);
  --mUnmatched;
  ++mProfAugments;
#ifdef MBX_CLASSIC_EVENT_TRACE
  traceAction(ClassicEventTrace::AugmentArc,u,v,cNullItm,cNullItm,cNullItm,tree);
#endif
#ifdef MBX_GABOW_SEARCH_AUDIT
  if (mBackend == kGabow2018) audEmitSearch('a');   // arc augment = one tree consumed
#endif
}

inline void FastMatchingEngineV2::augmentOnEdge(Size u, Size v, Size /*edgeIdxInU*/)
{
  Size Bu = topLevel(u);
  Size Bv = topLevel(v);
  if (Bu == cNullItm || Bv == cNullItm) return;
  BlossomStatus statusU = mBlsm[Bu].status;
  BlossomStatus statusV = mBlsm[Bv].status;
  Size treeU = mBlsm[Bu].tree;
  Size treeV = mBlsm[Bv].tree;

  if (statusU == BS_EVEN && statusV == BS_EVEN) {
    flipPathFromTopEven(Bu, u);
    flipPathFromTopEven(Bv, v);
    setMatePair(u, v);
    Bu = topLevel(u);
    Bv = topLevel(v);
    if (Bu != cNullItm && Bu < mBlsm.size())
      mBlsm[Bu].next = encodeEdge(u, v, mNumVtxs);
    if (Bv != cNullItm && Bv < mBlsm.size())
      mBlsm[Bv].next = encodeEdge(v, u, mNumVtxs);

    destroyTree(treeU);
    destroyTree(treeV);
    mUnmatched -= 2;
    ++mProfAugments;
#ifdef MBX_CLASSIC_EVENT_TRACE
    traceAction(ClassicEventTrace::AugmentEdge,u,v,cNullItm,cNullItm,cNullItm,treeU,treeV);
#endif
#ifdef MBX_GABOW_SEARCH_AUDIT
    if (mBackend == kGabow2018) audEmitSearch('e');  // edge augment = two trees consumed
#endif
    return;
  }

  (void)statusU;
  (void)statusV;
  (void)treeU;
  (void)treeV;
}

inline void FastMatchingEngineV2::shrinkOnEdge(Size u, Size v,
                                               Size /*edgeIdxInU*/, Size tree)
{
  Size Bu = topLevel(u);
  Size Bv = topLevel(v);
  if (Bu == cNullItm || Bv == cNullItm || Bu == Bv) return;

  std::vector<Size> leftPath, rightPath;
  collectTreePath(Bu, &leftPath);
  collectTreePath(Bv, &rightPath);

  Size nca = cNullItm;
  Size leftNca = cNullItm;
  Size rightNca = cNullItm;
  ++mShrinkCurrentStamp;
  if (mShrinkCurrentStamp == cNullItm) {
    std::fill(mShrinkStamp.begin(), mShrinkStamp.end(), 0);
    std::fill(mShrinkLeftIndex.begin(), mShrinkLeftIndex.end(), cNullItm);
    mShrinkCurrentStamp = 1;
  }
  if (mShrinkStamp.size() < mBlsm.size()) {
    mShrinkStamp.resize(mBlsm.size(), 0);
    mShrinkLeftIndex.resize(mBlsm.size(), cNullItm);
  }
  for (Size i = 0; i < leftPath.size(); ++i) {
    Size b = leftPath[i];
    if (b < mShrinkStamp.size()) {
      mShrinkStamp[b] = mShrinkCurrentStamp;
      mShrinkLeftIndex[b] = i;
    }
  }
  for (Size i = 0; i < rightPath.size(); ++i) {
    Size b = rightPath[i];
    if (b < mShrinkStamp.size() &&
        mShrinkStamp[b] == mShrinkCurrentStamp) {
      nca = b;
      leftNca = mShrinkLeftIndex[b];
      rightNca = i;
      break;
    }
  }
  if (nca == cNullItm) {
    Size edgeRefs = 0;
    for (Size i = 0; i < mAdj.size(); ++i) edgeRefs += mAdj[i].size();
    if (edgeRefs / 2 <= 50000) abortRun("shrink-no-nca");
    return;
  }

  if (leftNca == cNullItm || rightNca == cNullItm) return;

  std::vector<Size> subblossoms;
  std::vector<Size> childNext;
  Size prev = encodeEdge(u, v, mNumVtxs);
  for (Size i = 0; i < leftNca; i += 2) {
    if (i + 1 >= leftPath.size()) {
      abortRun("shrink-left-path-malformed");
      return;
    }
    subblossoms.push_back(leftPath[i]);
    childNext.push_back(prev);

    subblossoms.push_back(leftPath[i + 1]);
    childNext.push_back(mBlsm[leftPath[i + 1]].next);
    prev = reverseEncodedEdge(mBlsm[leftPath[i + 1]].pred);
  }

  subblossoms.push_back(nca);
  childNext.push_back(prev);

  for (int i = static_cast<int>(rightNca) - 2; i >= 0; i -= 2) {
    Size odd = rightPath[static_cast<Size>(i + 1)];
    Size even = rightPath[static_cast<Size>(i)];
    subblossoms.push_back(odd);
    childNext.push_back(mBlsm[odd].pred);
    subblossoms.push_back(even);
    childNext.push_back(mBlsm[even].next);
  }

  if (subblossoms.size() < 3 || (subblossoms.size() % 2) == 0) {
    return;
  }

  std::vector<Size> cycleEdges;
  for (Size i = 0; i < subblossoms.size(); ++i) {
    Size from = subblossoms[i];
    Size to = subblossoms[(i + 1) % subblossoms.size()];
    Size edge = cNullItm;
    if (from == Bu && to == Bv) edge = encodeEdge(u, v, mNumVtxs);
    else if (from == Bv && to == Bu) edge = encodeEdge(v, u, mNumVtxs);
    else edge = orientedTreeEdge(from, to);
    if (edge == cNullItm) {
      return;
    }
    cycleEdges.push_back(edge);
  }
  if (childNext.size() != subblossoms.size()) {
    abortRun("shrink-child-next-size");
    return;
  }

  Size surfacePred = mBlsm[nca].pred;
  Size surfaceBase = mBlsm[nca].base;

  std::vector<ConcatenableQueue::Handle> subRoots;
  std::vector<ConcatenableQueue::Handle> splitLeaves;
  for (Size i = 0; i < subblossoms.size(); ++i) {
    Size child = subblossoms[i];
    if (child >= mBlsm.size() || !mBlsm[child].alive ||
        mBlsm[child].parent != cNullItm) {
      abortRun("shrink-child-invalid");
      return;
    }
    if (mBlsm[child].cqRoot == cNullItm ||
        mBlsm[child].cqRoot >= mCQ.poolSize()) {
      Size edgeRefs = 0;
      for (Size j = 0; j < mAdj.size(); ++j) edgeRefs += mAdj[j].size();
      abortRun("shrink-cq-root-invalid");
      return;
    }
    if (mBlsm[child].cqLastLeaf == cNullItm ||
        mBlsm[child].cqLastLeaf >= mCQ.poolSize()) {
      abortRun("shrink-cq-last-invalid");
      return;
    }
    subRoots.push_back(mBlsm[child].cqRoot);
    if (i + 1 < subblossoms.size()) {
      splitLeaves.push_back(mBlsm[child].cqLastLeaf);
    }
  }

  std::vector<Size> shrinkSurfaceScanVtxs;
  for (Size i = 0; i < subblossoms.size(); ++i) {
    Size child = subblossoms[i];
    BlossomStatus oldStatus = mBlsm[child].status;
    if (oldStatus != BS_EVEN) {
      std::vector<Size> childVtxs;
      MembershipView pChildVtxs = &mBlsm[child].vertices;
      if (pChildVtxs->empty()) {
        collectBlossomVertices(child, &childVtxs);
        pChildVtxs = &childVtxs;
      }
      shrinkSurfaceScanVtxs.insert(shrinkSurfaceScanVtxs.end(),
                                   pChildVtxs->begin(), pChildVtxs->end());
    }
    if (oldStatus == BS_ODD) bakeOddToMatchedForShrink(child);
    else if (oldStatus == BS_EVEN) bakeEvenToMatchedForShrink(child);
    else mBlsm[child].status = BS_MATCHED;
    mBlsm[child].status = BS_MATCHED;
    mBlsm[child].tree = cNullItm;
    mBlsm[child].next = childNext[i];
  }

  ++mProfShrinks;
  mAnyShrinkYet = true;
  Size surface = acquireBlsmId();
  if (surface >= mBlsm.size()) return;
  BlossomData& B = mBlsm[surface];
  B.status = BS_EVEN;
  B.parent = cNullItm;
  B.base = surfaceBase;
  B.tree = tree;
  B.pred = surfacePred;
  B.next = surfacePred;
  B.pot = 0;
  B.offset = mDeltaSum;
  B.children = subblossoms;
  B.childEdges = cycleEdges;
  B.splitLeaves = splitLeaves;
  B.vertices.clear();
  B.childOfVtx.clear();
  B.d1MinVtx = cNullItm;
  B.d1MinValid = false;
  B.d2Key = infKey();
  B.d2Vtx = cNullItm;
  B.d2Edge = cNullItm;
  B.d2Tree = cNullItm;
  B.d2Queued = false;
  long double profVtxBuildT0 = MBX_HOTNOW();
  for (Size i = 0; i < subblossoms.size(); ++i) {
    std::vector<Size> childVtxs;
    MembershipView pChildVtxs = &mBlsm[subblossoms[i]].vertices;
    if (pChildVtxs->empty()) {
      collectBlossomVertices(subblossoms[i], &childVtxs);
      pChildVtxs = &childVtxs;
    }
    for (Size j = 0; j < pChildVtxs->size(); ++j) {
      Size cv = (*pChildVtxs)[j];
      B.vertices.push_back(cv);
      B.childOfVtx.push_back(std::make_pair(cv, i));
      if (B.d1MinVtx == cNullItm ||
          mVtx[cv].pot < mVtx[B.d1MinVtx].pot) {
        B.d1MinVtx = cv;
        B.d1MinValid = true;
      }
    }
  }
  // Sort childOfVtx by vtx key once, so lookups use lower_bound. (Each vtx appears
  // in exactly one sub-blossom, so keys are unique -- no dedup needed.)
#ifdef MBX_MEMBERSHIP_REUSE_VERIFY
  std::vector<std::pair<Size, Size> > referenceMembership(B.childOfVtx);
  std::sort(referenceMembership.begin(), referenceMembership.end());
#endif
  // Keep B.vertices, pair construction and minimum-potential traversal intact.
  // Only replace the final full sort with reuse of a complete sorted child map.
  MR_ADD(contractions,1);MR_ADD(vertices,B.childOfVtx.size());
  MR_ADD(reusedSurfaces,B.generation>1);
  Size dominantIndex=cNullItm,dominantStart=0,dominantCount=0;
  for(Size at=0;at<B.childOfVtx.size();) {
    const Size index=B.childOfVtx[at].second;
    Size end=at+1;
    while(end<B.childOfVtx.size()&&B.childOfVtx[end].second==index)++end;
    const std::vector<std::pair<Size,Size> >& cached=mBlsm[subblossoms[index]].childOfVtx;
    if(!trivial(subblossoms[index]))MR_ADD(nestedChildren,1);
    if(cached.empty())MR_ADD(emptyCaches,1);
    if(cached.size()==end-at&&end-at>dominantCount) {
      dominantIndex=index;dominantStart=at;dominantCount=end-at;
    }
    at=end;
  }
  bool mergedMembership=false;
  if(dominantIndex!=cNullItm&&dominantCount>=B.childOfVtx.size()-dominantCount) {
    mergedMembership=mbxDominantMembershipMerge(B.childOfVtx,
      mBlsm[subblossoms[dominantIndex]].childOfVtx,
      dominantIndex,dominantStart,dominantCount);
  }
  if(mergedMembership) {
    MR_ADD(shortcuts,1);MR_ADD(donorVertices,dominantCount);
    MR_ADD(remainderVertices,B.childOfVtx.size()-dominantCount);
  } else {
    MR_ADD(fallbacks,1);std::sort(B.childOfVtx.begin(),B.childOfVtx.end());
  }
#ifdef MBX_MEMBERSHIP_REUSE_VERIFY
  if(B.childOfVtx!=referenceMembership)
    throw std::runtime_error("membership-reuse full-sort differential mismatch");
  MR_ADD(verified,1);
#endif
#ifdef GMG_MAXBLOSSOM_PROBE
  { extern long long gGmgMaxBlossom; if ((long long)B.vertices.size() > gGmgMaxBlossom) gGmgMaxBlossom = (long long)B.vertices.size(); }
#endif
  mProfShrinkVtxBuildTime += MBX_HOTNOW() - profVtxBuildT0;
  long double profCqMergeT0 = MBX_HOTNOW();
  B.cqRoot = mCQ.merge(subRoots);
  mProfShrinkCqMergeTime += MBX_HOTNOW() - profCqMergeT0;
  B.cqLastLeaf = mBlsm[subblossoms.back()].cqLastLeaf;
  B.materializedBase = cNullItm;
  B.listedTree = cNullItm;
  B.topCache = surface;
  B.topCacheGeneration = B.generation;
  B.alive = true;

  for (Size i = 0; i < subblossoms.size(); ++i) {
    Size child = subblossoms[i];
    eraseD2BlossomEntry(child);
    mBlsm[child].parent = surface;
    mBlsm[child].topCache = surface;
    mBlsm[child].topCacheGeneration = B.generation;
    (void)child;
  }

  // Monotone support: add missing inner atoms and merge actual parent paths.
  // No full-inner BM owner is manufactured or later split during expansion.
  if(mBackend==kGabow2018)supportBlossom(surface,u,v);

  addTreeMember(surface, tree);
  { long double t0 = MBX_HOTNOW();
    refreshEvenBlossomEdgesSubset(surface, shrinkSurfaceScanVtxs);
    mProfShrinkRefreshTime += MBX_HOTNOW() - t0; }
  pushD1ForBlossom(surface);

  // A child is now represented by the parent's flattened membership cache.
#ifdef MBX_CLASSIC_EVENT_TRACE
  if (mEventTrace.enabled()) {
    if (mEventMembers.size()<mBlsm.size())mEventMembers.resize(mBlsm.size());
    ClassicEventTrace::Signature signature;
    bool traceNested=false;
    for (Size i=0;i<subblossoms.size();++i) {
      signature.add(mEventMembers[subblossoms[i]]);
      if (mEventMembers[subblossoms[i]].count>1)traceNested=true;
    }
    mEventTrace.shrinkLifecycle(mEventMembers[surface].count>0,traceNested,signature.count);
    mEventMembers[surface]=signature;
    traceAction(ClassicEventTrace::Shrink,u,v,surface,cNullItm,cNullItm,tree);
  }
#endif
  // Keeping the child's own flattened vertices and vertex-to-child map would
  // retain one full copy at every level of a nested blossom chain. The
  // hierarchy (children/parent) remains intact for expansion and for the
  // uncached membership queries below. Trivial blossoms keep their one vertex.
  for (Size i = 0; i < subblossoms.size(); ++i) {
    Size child = subblossoms[i];
    if (trivial(child)) continue;
    mBlsm[child].vertices.release();
    std::vector<std::pair<Size, Size> >().swap(mBlsm[child].childOfVtx);
  }
#ifdef MBX_GABOW_FRONTIER_VERIFY
  verifyGabowFrontierStable();
#endif
}


inline void FastMatchingEngineV2::addBvAtomicsToGabowForest(Size Bv,Size u,Size v)
{
  supportInstallPath(Bv,v,u);
}

// Gabow-only frontier cache restoration. CQ and membership use the same ordered
// child concatenation; do not use this to retain flattened nested ancestors.
inline bool FastMatchingEngineV2::restoreGabowFrontierMembership(
    Size parent,const std::vector<Size>& children,
    const std::vector<ConcatenableQueue::Handle>& roots,
    Membership* savedParent)
{
  if(mBackend!=kGabow2018)return true;
  if(savedParent==NULL || parent>=mBlsm.size() || trivial(parent) || !mBlsm[parent].alive || children.empty() || children.size()!=roots.size()) {
    abortRun("frontier-invalid-header");return false;
  }
  savedParent->clear();
  Size total=0;
  for(Size i=0;i<children.size();++i) {
    if(children[i]>=mBlsm.size() || !mBlsm[children[i]].alive || roots[i]==cNullItm || roots[i]>=mCQ.poolSize()) {
      abortRun("frontier-invalid-child-root");return false;
    }
    Size count=mCQ.size(roots[i]);
    if(count==0 || total>mNumVtxs || count>mNumVtxs-total) {
      abortRun("frontier-invalid-range-size");return false;
    }
    total+=count;
  }
  // Rebuild only if the parent cache is genuinely absent, while its hierarchy
  // still exists. Normal candidate splits transfer a pre-existing parent list.
  const bool cached=!mBlsm[parent].vertices.empty();
  if(!cached){std::vector<Size> fallback;collectBlossomVertices(parent,&fallback);savedParent->assign(fallback.begin(),fallback.end());}
  const Membership& members=cached?mBlsm[parent].vertices:*savedParent;
  if(total!=members.size()) {abortRun("frontier-parent-size-mismatch");return false;}
  Size offset=0;
  for(Size i=0;i<children.size();++i) {
    const Size child=children[i],count=mCQ.size(roots[i]);
    if(count>members.size()-offset) {abortRun("frontier-range-overflow");return false;}
    if(trivial(child) && (count!=1 || members[offset]!=child || mBlsm[child].vertices.size()!=1 || mBlsm[child].vertices[0]!=child)) {
      abortRun("frontier-invalid-singleton");return false;
    }
#ifdef MBX_GABOW_FRONTIER_VERIFY
    // Independently traverse hierarchy, explicitly ignoring nontrivial caches.
    std::vector<Size> expected,pending(1,child);
    while(!pending.empty()) {
      Size b=pending.back();pending.pop_back();
      if(b>=mBlsm.size())throw std::runtime_error("frontier verification: bad hierarchy");
      if(trivial(b))expected.push_back(b);
      else for(Size j=mBlsm[b].children.size();j>0;--j)pending.push_back(mBlsm[b].children[j-1]);
    }
    if(expected.size()!=count || !std::equal(expected.begin(),expected.end(),members.begin()+offset))
      throw std::runtime_error("frontier verification: ordered hierarchy mismatch");
    std::vector<ConcatenableQueue::Handle> leaves;
    mCQ.collectLeaves(roots[i],leaves);
    if(leaves.size()!=count)throw std::runtime_error("frontier verification: CQ size mismatch");
    for(Size j=0;j<count;++j)if(mCQ.payload(leaves[j])!=members[offset+j])
      throw std::runtime_error("frontier verification: ordered CQ mismatch");
    ++FrontierReuseAudit::state().orderedChecks;++FrontierReuseAudit::state().cqChecks;
#endif
    offset+=count;
  }
  if(offset!=members.size()) {abortRun("frontier-unconsumed-parent");return false;}
  // No mutation of child caches occurs until the complete partition is valid.
  if(cached)savedParent->swap(mBlsm[parent].vertices);
#ifdef MBX_GABOW_FRONTIER_AUDIT
  FrontierReuseAudit::State& observed=FrontierReuseAudit::state();
  ++observed.calls;if(cached)++observed.parentHits;else ++observed.parentFallbacks;
  observed.parentLeaves+=savedParent->size();observed.children+=children.size();
  FrontierReuseAudit::maximum(observed.maxParentSize,savedParent->size());
  FrontierReuseAudit::maximum(observed.maxParentCapacity,savedParent->capacity());
#endif

  Size largest=0,largestOffset=0;
  offset=0;
  for(Size i=0;i<children.size();++i) {
    const Size count=mCQ.size(roots[i]);
    if(count>mCQ.size(roots[largest])) {largest=i;largestOffset=offset;}
    offset+=count;
  }
  LR_ADD(calls,1); LR_ADD(parentLeaves,savedParent->size());
  LR_ADD(largestLeaves,mCQ.size(roots[largest]));
  offset=0;
  for(Size i=0;i<children.size();++i) {
    const Size child=children[i],count=mCQ.size(roots[i]);
    if(!trivial(child)) {
      if(i!=largest) {
        mBlsm[child].vertices.assign(savedParent->begin()+offset,savedParent->begin()+offset+count);
        LR_ADD(smallerCopied,count);
      }
#ifdef MBX_GABOW_FRONTIER_AUDIT
      ++observed.restored;observed.restoredLeaves+=count;
#endif
    }
#ifdef MBX_GABOW_FRONTIER_AUDIT
    else ++observed.singletons;
#endif
    offset+=count;
  }
  if(!trivial(children[largest])) {
    std::size_t copied=0;
    if(!mBlsm[children[largest]].vertices.takeRange(*savedParent,largestOffset,mCQ.size(roots[largest]),&copied))
      throw std::runtime_error("largest membership range transfer failed");
    LR_ADD(compactedCopied,copied); LR_ADD(largestMoved,1);
  } else savedParent->release();
#ifdef MBX_LARGEST_REUSE_VERIFY
  for(Size i=0;i<children.size();++i)
    if(!mBlsm[children[i]].vertices.storageWithinBound())
      throw std::runtime_error("largest exposed membership exceeds backing bound");
#endif
  return true;
}


inline void FastMatchingEngineV2::verifyLargestBmPartition(
    const std::vector<Size>& children, Size oldBm, bool after) const
{
#ifdef MBX_LARGEST_REUSE_VERIFY
  // This independent checker is intentionally expensive and never timed as
  // clean solver work. Epoch stamps avoid per-expansion O(n) clearing.
  static std::vector<unsigned long long> seen;
  static unsigned long long epoch=0;
  if(seen.size()!=mNumVtxs) {seen.assign(mNumVtxs,0);epoch=0;}
  if(++epoch==0) {std::fill(seen.begin(),seen.end(),0);epoch=1;}
  Size total=0;
  if(!after) {
    for(Size i=0;i<children.size();++i)for(Size j=0;j<mBlsm[children[i]].vertices.size();++j) {
      const Size w=mBlsm[children[i]].vertices[j];
      if(w>=mNumVtxs || seen[w]==epoch || mBM.mVtxBlossom[w]!=oldBm)
        throw std::runtime_error("largest pre partition duplicate/foreign atom");
      seen[w]=epoch;++total;
    }
    if(total!=mBM.blossom_size(oldBm))throw std::runtime_error("largest pre missing atoms");
    const Size head=mBM.mBlossoms[oldBm].vtxListHead;
    Size cur=head;
    for(Size j=0;j<total;++j) {
      if(cur>=mNumVtxs || seen[cur]!=epoch || mBM.mVtxBlossom[cur]!=oldBm)
        throw std::runtime_error("largest pre ring ownership");
      // Clear visited stamp: a repeated short cycle cannot pass the size test.
      seen[cur]=0;
      const Size next=mBM.mVtxNext[cur],prev=mBM.mVtxPrev[cur];
      if(next>=mNumVtxs || prev>=mNumVtxs || mBM.mVtxPrev[next]!=cur || mBM.mVtxNext[prev]!=cur)
        throw std::runtime_error("largest pre broken ring");
      cur=next;
    }
    if(cur!=head)throw std::runtime_error("largest pre ring size");
    LR_ADD(preChecks,1);
  } else {
    for(Size i=0;i<children.size();++i) {
      const Size child=children[i],id=mMatchboxToBm[child];
      MembershipView members=mBlsm[child].vertices;
      if(!mBM.blossom_alive(id) || mBM.blossom_size(id)!=members.size())
        throw std::runtime_error("largest post child header");
      Size expectedRoot=members[0];
      for(Size j=0;j<members.size();++j) {
        const Size w=members[j];
        if(w>=mNumVtxs || seen[w]==epoch || mBM.mVtxBlossom[w]!=id || mMatchboxToBm[w]!=id)
          throw std::runtime_error("largest post child membership/map");
        seen[w]=epoch;
        if(mBM.mDepth[w]<mBM.mDepth[expectedRoot])expectedRoot=w;
      }
      if(mBM.mBlossoms[id].rootVertex!=expectedRoot || mBM.mBlossoms[id].rank!=mBM.logFloor(members.size()))
        throw std::runtime_error("largest post child root/rank");
      Size head=mBM.mBlossoms[id].vtxListHead,cur=head;
      for(Size j=0;j<members.size();++j) {
        if(cur>=mNumVtxs || seen[cur]!=epoch || mBM.mVtxBlossom[cur]!=id)
          throw std::runtime_error("largest post ring membership/size");
        seen[cur]=0;
        const Size next=mBM.mVtxNext[cur],prev=mBM.mVtxPrev[cur];
        if(next>=mNumVtxs || prev>=mNumVtxs || mBM.mVtxPrev[next]!=cur || mBM.mVtxNext[prev]!=cur)
          throw std::runtime_error("largest post broken ring");
        cur=next;
      }
      if(cur!=head)throw std::runtime_error("largest post ring closure");
    }
    LR_ADD(postChecks,1);
  }
#else
  (void)children;(void)oldBm;(void)after;
#endif
}

inline void FastMatchingEngineV2::verifyGabowFrontierStable(bool force) const
{
#ifdef MBX_GABOW_FRONTIER_VERIFY
  if(mBackend!=kGabow2018 || mAbortRun || (!force && mNumVtxs>4096))return;
  unsigned long long liveSize=0,liveCapacity=0,deadCapacity=0,scratchCapacity=0,liveMapCapacity=0,deadMapCapacity=0;
  std::vector<unsigned char> frontierOwner(mNumVtxs,0);
  for(Size b=0;b<mBlsm.size();++b) {
    const BlossomData& B=mBlsm[b];
    if(B.alive) {
      if(!B.vertices.storageWithinBound())throw std::runtime_error("largest stable backing bound");
      liveSize+=B.vertices.size();liveCapacity+=B.vertices.capacity();
      liveMapCapacity+=B.childOfVtx.capacity();
      if(!trivial(b) && B.parent!=cNullItm && (!B.vertices.empty() || B.vertices.capacity()!=0 || !B.childOfVtx.empty() || B.childOfVtx.capacity()!=0))
        throw std::runtime_error("frontier verification: retained nested cache");
      if(!trivial(b) && B.parent==cNullItm && B.vertices.empty())
        throw std::runtime_error("frontier verification: missing exposed cache");
      if(!trivial(b) && B.parent==cNullItm)for(Size j=0;j<B.vertices.size();++j) {
        Size v=B.vertices[j];
        if(v>=mNumVtxs || frontierOwner[v])throw std::runtime_error("frontier verification: overlapping exposed membership");
        frontierOwner[v]=1;
      }
    } else {
      deadCapacity+=B.vertices.capacity();
      deadMapCapacity+=B.childOfVtx.capacity();
      if(B.vertices.capacity()>4096)throw std::runtime_error("frontier verification: oversized dead cache");
      if(B.childOfVtx.capacity()>4096)throw std::runtime_error("frontier verification: oversized dead membership map");
    }
  }
  for(Size i=0;i<mSplitChildVtxsScratch.size();++i)scratchCapacity+=mSplitChildVtxsScratch[i].capacity();
  if(scratchCapacity!=0)throw std::runtime_error("frontier verification: legacy scratch populated");
  if(liveSize>2*(unsigned long long)mNumVtxs)throw std::runtime_error("frontier verification: non-disjoint live cache sizes");
  FrontierReuseAudit::State& s=FrontierReuseAudit::state();++s.stableChecks;
  FrontierReuseAudit::maximum(s.maxLiveSize,liveSize);FrontierReuseAudit::maximum(s.maxLiveCapacity,liveCapacity);
  FrontierReuseAudit::maximum(s.maxDeadCapacity,deadCapacity);FrontierReuseAudit::maximum(s.maxScratchCapacity,scratchCapacity);
  FrontierReuseAudit::maximum(s.maxLiveMapCapacity,liveMapCapacity);FrontierReuseAudit::maximum(s.maxDeadMapCapacity,deadMapCapacity);
#else
  (void)force;
#endif
}


// Records may be owned outside the parent currently being partitioned. Keep
// their exact generations while handles are cleared, but do not run BM removal
// until the partition's coverage/ring invariants have been completely restored.
inline void FastMatchingEngineV2::stageGabowTwins(
    Size u, Size i, GabowDeferredTwins& deferred)
{
  if(mBmHandleU[u][i]!=cNullItm) {
    deferred.push_back(std::make_pair(mBmHandleU[u][i],mBmGenU[u][i]));
    TC_ADD(saved,1);
  }
  if(mBmHandleV[u][i]!=cNullItm) {
    deferred.push_back(std::make_pair(mBmHandleV[u][i],mBmGenV[u][i]));
    TC_ADD(saved,1);
  }
  TC_PEAK(deferred.size());
}

inline void FastMatchingEngineV2::removeStagedGabowTwins(
    const GabowDeferredTwins& deferred)
{
  if(deferred.empty())return;
  // Only called after finish_largest_child_partition, before any republication.
  // No broad owner deactivation, no payload-based deletion, no stale raw handle.
  TC_ADD(batches,1);
  mBM.begin_edge_removal_batch();
  for(Size k=0;k<deferred.size();++k) {
    TC_ADD(attempted,1);
    if(mBM.remove_edge_if_batched(deferred[k].first,deferred[k].second)) {
      TC_ADD(removed,1);
    } else {
      TC_ADD(stale,1);
    }
  }
  mBM.end_edge_removal_batch();
#ifdef MBX_TWIN_CLEANUP_VERIFY
  for(Size k=0;k<deferred.size();++k) {
    Size tag=0;
    if(mBM.edge_occurrence_tag_if(deferred[k].first,deferred[k].second,&tag))
      throw std::runtime_error("staged Gabow twin survived exact cleanup");
    TC_ADD(verified,1);
  }
#endif
}

inline void FastMatchingEngineV2::splitBlossom(Size blossom)
{
  if (blossom >= mBlsm.size() || trivial(blossom) ||
      !mBlsm[blossom].alive || mBlsm[blossom].parent != cNullItm ||
      mBlsm[blossom].status != BS_ODD) {
    return;
  }
#ifdef MBX_GABOW_SEARCH_AUDIT
  ++mAudSplits;
#endif

  std::vector<Size> children = mBlsm[blossom].children;
#ifdef MBX_CLASSIC_EVENT_TRACE
  const Size traceBase=mBlsm[blossom].base;
  const ClassicEventTrace::Signature traceMembers=mEventTrace.enabled()?mEventMembers[blossom]:ClassicEventTrace::Signature();
#endif
  std::vector<Size> childEdges = mBlsm[blossom].childEdges;
  Size tree = mBlsm[blossom].tree;
  Size pred = mBlsm[blossom].pred;
  Size next = mBlsm[blossom].next;
  ConcatenableQueue::Handle mergedRoot = mBlsm[blossom].cqRoot;
  ConcatenableQueue::Handle mergedLastLeaf = mBlsm[blossom].cqLastLeaf;
  std::vector<ConcatenableQueue::Handle> splitLeaves =
    mBlsm[blossom].splitLeaves;
  if (children.empty() || childEdges.size() != children.size() ||
      pred == cNullItm || next == cNullItm ||
      edgeSrc(pred) >= mNumVtxs || edgeDst(pred) >= mNumVtxs ||
      edgeSrc(next) >= mNumVtxs || edgeDst(next) >= mNumVtxs) {
    abortRun("split-invalid-header");
    return;
  }

  Size ib = childIndexContaining(blossom, edgeDst(pred));
  Size id = childIndexContaining(blossom, edgeSrc(next));
  if (ib == cNullItm || id == cNullItm ||
      ib >= children.size() || id >= children.size()) {
    abortRun("split-entry-index-invalid");
    return;
  }
  Val parentOddDelta = mDeltaSum - mBlsm[blossom].offset;
  if (parentOddDelta != 0) {
    std::vector<Size> parentVtxs;
    MembershipView pParentVtxs = &mBlsm[blossom].vertices;
    if (pParentVtxs->empty()) {
      collectBlossomVertices(blossom, &parentVtxs);
      pParentVtxs = &parentVtxs;
    }
    for (Size i = 0; i < pParentVtxs->size(); ++i)
      mVtx[(*pParentVtxs)[i]].pot += parentOddDelta;
  }

  std::vector<ConcatenableQueue::Handle> splitRoots;
  if (mergedRoot != cNullItm && children.size() > 1) {
    if (mergedRoot >= mCQ.poolSize()) {
      abortRun("split-root-invalid");
      return;
    }
    for (Size i = 0; i < splitLeaves.size(); ++i) {
      if (splitLeaves[i] >= mCQ.poolSize()) {
        abortRun("split-leaf-invalid");
        return;
      }
    }
    mCQ.split(mergedRoot, splitLeaves, splitRoots);
  } else {
    splitRoots.assign(children.size(), mergedRoot);
  }
  if (splitRoots.size() != children.size()) {
    abortRun("split-root-count-mismatch");
    return;
  }

  if(mBackend==kGabow2018) {
    Membership gabowSavedParentMembership;
    if(!restoreGabowFrontierMembership(blossom,children,splitRoots,&gabowSavedParentMembership))return;
  } // Child lists own their data; release the detached parent before BM work.

  for (Size i = 0; i < children.size(); ++i) {
    Size child = children[i];
    if (child >= mBlsm.size() || !mBlsm[child].alive) {
      abortRun("split-child-invalid");
      return;
    }
    mBlsm[child].parent = cNullItm;
    mBlsm[child].cqRoot = splitRoots[i];
    mBlsm[child].cqLastLeaf =
      (i < splitLeaves.size()) ? splitLeaves[i] : mergedLastLeaf;
    mBlsm[child].status = BS_MATCHED;
    mBlsm[child].tree = cNullItm;
    mBlsm[child].offset = 0;
    mBlsm[child].materializedBase = cNullItm;
    mBlsm[child].d1MinValid = false;
    mBlsm[child].d2Key = infKey();
    mBlsm[child].d2Vtx = cNullItm;
    mBlsm[child].d2Edge = cNullItm;
    mBlsm[child].d2Tree = cNullItm;
    mBlsm[child].d2Queued = false;
    mBlsm[child].topCache = child;
    mBlsm[child].topCacheGeneration = mBlsm[child].generation;
    setTopLevelForBlossom(child, child);
  }

  releaseBlsmId(blossom);

  // Expansion preserves the supporting forest and existing singleton path
  // owners. matchedToEven adds only newly outer atoms before publication.
  if(mBackend==kGabow2018){mMatchboxToBm[blossom]=cNullItm;MS_ADD(splits,1);}

  bool forward = id > ib ? (((id - ib) % 2) == 0) :
                           (((ib - id) % 2) == 1);
  Size cur = ib;
  Size x = edgeSrc(pred);
  Size y = edgeDst(pred);
  Size guard = 0;
  while (cur != id && !mAbortRun) {
    if (++guard > children.size()) {
      abortRun("split-walk-guard");
      return;
    }
    Size tChild = children[cur];
    if (tChild >= mBlsm.size()) {
      abortRun("split-walk-child-invalid");
      return;
    }

    mBlsm[tChild].tree = tree;
    mBlsm[tChild].pred = encodeEdge(x, y, mNumVtxs);
    if (y < mNumVtxs) mTreeParent[y] = x;
    addTreeMember(tChild, tree);
    matchedToOdd(tChild);

    Size sIdx;
    Size nextT;
    Size matchEdge;
    Size nextEdge;
    if (forward) {
      sIdx = (cur + 1) % children.size();
      nextT = (cur + 2) % children.size();
      matchEdge = childEdges[cur];
      nextEdge = childEdges[sIdx];
      x = edgeSrc(nextEdge);
      y = edgeDst(nextEdge);
      cur = nextT;
    } else {
      sIdx = (cur + children.size() - 1) % children.size();
      nextT = (cur + children.size() - 2) % children.size();
      matchEdge = reverseEncodedEdge(childEdges[sIdx]);
      nextEdge = childEdges[nextT];
      x = edgeDst(nextEdge);
      y = edgeSrc(nextEdge);
      cur = nextT;
    }

    if (matchEdge == cNullItm || nextEdge == cNullItm ||
        sIdx >= children.size()) {
      abortRun("split-walk-edge-invalid");
      return;
    }

    Size sChild = children[sIdx];
    if (sChild >= mBlsm.size()) {
      abortRun("split-walk-even-child-invalid");
      return;
    }
    mBlsm[tChild].next = matchEdge;
    mBlsm[sChild].pred = matchEdge;
    if (edgeDst(matchEdge) < mNumVtxs)
      mTreeParent[edgeDst(matchEdge)] = edgeSrc(matchEdge);
    matchedToEven(sChild, tree);
  }

  if (!mAbortRun && cur == id && children[id] < mBlsm.size()) {
    mBlsm[children[id]].status = BS_ODD;
    mBlsm[children[id]].tree = tree;
    mBlsm[children[id]].pred = encodeEdge(x, y, mNumVtxs);
    mBlsm[children[id]].next = next;
    if (y < mNumVtxs) mTreeParent[y] = x;
    addTreeMember(children[id], tree);
    matchedToOdd(children[id]);
  }

  if (!mAbortRun) {
    for (Size i = 0; i < children.size(); ++i) {
      Size child = children[i];
      if (child < mBlsm.size() && mBlsm[child].alive &&
          mBlsm[child].parent == cNullItm &&
          mBlsm[child].status == BS_MATCHED) {
        refreshMatchedBlossomD2(child);
        // Off-path matched atoms were never inserted into this live forest.
        // Do not retire/reparent any supporting-tree node during expansion.

      }
    }
  }
#ifdef MBX_CLASSIC_EVENT_TRACE
  if (mEventTrace.enabled() && !mAbortRun)
    mEventTrace.commit(ClassicEventTrace::Expand,edgeSrc(pred),edgeDst(pred),traceBase,traceMembers,
      edgeSrc(next),edgeDst(next),tree);
#endif
#ifdef MBX_GABOW_FRONTIER_VERIFY
  verifyGabowFrontierStable();
#endif
}

inline Size FastMatchingEngineV2::acquireBlsmId()
{
  if (!mFreeBlsmIds.empty()) {
    Size id = mFreeBlsmIds.back();
    mFreeBlsmIds.pop_back();
    ++mBlsm[id].generation;
    mBlsm[id].topCache = cNullItm;
    mBlsm[id].topCacheGeneration = 0;
    // R_Scope mirror fix: a recycled id must not inherit its previous life's
    // mBM image (dissolving via a stale image would dispose a live blossom).
    if (id < mMatchboxToBm.size()) mMatchboxToBm[id] = cNullItm;
    return id;
  }
  Size id = mBlossomCap;
  mBlsm.push_back(BlossomData());
  mShrinkStamp.push_back(0);
  mShrinkLeftIndex.push_back(cNullItm);
  // R_Scope mirror fix: keep mMatchboxToBm in lockstep with mBlsm. It was
  // sized once at init (mBlossomCap), so every super-blossom allocated later
  // fell OUTSIDE the map: the shrink hook's guarded write silently no-op'd,
  // the dissolve mirror then found cNullItm ("don't touch mBM") and mBM kept
  // the vertices FUSED while matchbox split them. Invisible under the
  // permanent-star forest (rho-anchored back-edges can never be refused by
  // make_edge's Bv==Bw check); fatal under the deep forest (lost candidates
  // -> valid-but-suboptimal matchings).
  if (mMatchboxToBm.size() < mBlsm.size())
    mMatchboxToBm.resize(mBlsm.size(), cNullItm);
  mBlsm[id].generation = 1;
  mBlsm[id].topCache = cNullItm;
  mBlsm[id].topCacheGeneration = 0;
  ++mBlossomCap;
  return id;
}

inline void FastMatchingEngineV2::releaseBlsmId(Size id)
{
  BlossomData& B = mBlsm[id];
  B.alive = false;
  B.children.clear();
  B.childEdges.clear();
  B.splitLeaves.clear();
  B.vertices.clear();
  B.childOfVtx.clear();
  // Recycled blossom IDs used to retain the largest buffers from every
  // historical blossom. On a 262K seed this left 575 MiB in vertices plus
  // 1.145 GiB in childOfVtx immediately before bad_alloc. Preserve reuse for
  // ordinary small blossoms, but return giant dissolved-blossom buffers to
  // the allocator. This bounds historical retention without changing state.
  const Size retainCap = 4096;
  if (B.children.capacity() > retainCap) std::vector<Size>().swap(B.children);
  if (B.childEdges.capacity() > retainCap) std::vector<Size>().swap(B.childEdges);
  if (B.splitLeaves.capacity() > retainCap)
    std::vector<ConcatenableQueue::Handle>().swap(B.splitLeaves);
  if (B.vertices.capacity() > retainCap) B.vertices.release();
  if (B.childOfVtx.capacity() > retainCap)
    std::vector<std::pair<Size, Size> >().swap(B.childOfVtx);
  B.cqRoot = cNullItm;
  B.cqLastLeaf = cNullItm;
  B.materializedBase = cNullItm;
  B.d1MinVtx = cNullItm;
  B.d1MinValid = false;
  B.d2Key = infKey();
  B.d2Vtx = cNullItm;
  B.d2Edge = cNullItm;
  B.d2Tree = cNullItm;
  B.d2Queued = false;
  ++B.generation;
  B.topCache = cNullItm;
  B.topCacheGeneration = 0;
  B.listedTree = cNullItm;
  mFreeBlsmIds.push_back(id);
}

inline void FastMatchingEngineV2::abortRun(const char* reason)
{
#ifdef MBX_CLASSIC_EVENT_TRACE
  mEventTrace.aborted(reason);
#endif
  (void)reason;
  mAbortRun = true;
}

inline Val FastMatchingEngineV2::vertexBestSlack(Size v) const
{
  Val best = infKey();
  std::map<Size, std::pair<Val, Size> >::const_iterator it;
  for (it = mVtx[v].bestByTree.begin(); it != mVtx[v].bestByTree.end(); ++it) {
    if (it->second.first < best) best = it->second.first;
  }
  return best;
}

inline void FastMatchingEngineV2::updateVertexCQLeaf(Size v)
{
  if (mVtx[v].cqLeaf == cNullItm) return;
  mCQ.setPriority(mVtx[v].cqLeaf, vertexBestSlack(v));
}

inline void FastMatchingEngineV2::clearVertexD2(Size v)
{
  eraseD2Entry(v);
  mVtx[v].bestByTree.clear();
  updateVertexCQLeaf(v);
  pushD2Best(v);
}

inline void FastMatchingEngineV2::eraseD2Entry(Size v)
{
  if (v >= mNumVtxs) return;
  Size B = topLevel(v);
  if (B == cNullItm || B >= mBlsm.size()) return;
  eraseD2BlossomEntry(B);
}

inline void FastMatchingEngineV2::eraseD2BlossomEntry(Size blossom)
{
  if (blossom >= mBlsm.size() || !mBlsm[blossom].d2Queued) return;
  HeapEntryD2 e;
  e.key = mBlsm[blossom].d2Key;
  e.matchedVtx = mBlsm[blossom].d2Vtx;
  e.edgeIdxInMatched = mBlsm[blossom].d2Edge;
  e.matchedBlsm = blossom;
  e.tree = mBlsm[blossom].d2Tree;
  e.version = 0;
  mPQD2.erase(e);
  mBlsm[blossom].d2Queued = false;
  mBlsm[blossom].d2Key = infKey();
  mBlsm[blossom].d2Vtx = cNullItm;
  mBlsm[blossom].d2Edge = cNullItm;
  mBlsm[blossom].d2Tree = cNullItm;
}

inline void FastMatchingEngineV2::addMatchedVertexD2Candidate(
  Size matchedVtx, Size edgeIdxInMatched)
{
#ifdef MBX_GABOW_SEARCH_AUDIT
  ++mAudD2CandAdds;
#endif
  if (matchedVtx >= mNumVtxs ||
      edgeIdxInMatched >= mAdj[matchedVtx].size())
    return;

  Size B = topLevel(matchedVtx);
  if (B == cNullItm || B >= mBlsm.size() ||
      mBlsm[B].status != BS_MATCHED) {
    clearVertexD2(matchedVtx);
    return;
  }

  Size even = mAdj[matchedVtx][edgeIdxInMatched];
  if (mMate[matchedVtx] == even) return;
  Size Be = topLevel(even);
  if (Be == cNullItm || Be == B || Be >= mBlsm.size() ||
      mBlsm[Be].status != BS_EVEN || !activeTree(mBlsm[Be].tree))
    return;

  Size tree = mBlsm[Be].tree;
  Val key = d2Key(matchedVtx, edgeIdxInMatched);
  std::map<Size, std::pair<Val, Size> >::iterator it =
    mVtx[matchedVtx].bestByTree.find(tree);
  if (it != mVtx[matchedVtx].bestByTree.end() &&
      it->second.first <= key + EPS)
    return;

  eraseD2Entry(matchedVtx);
  mVtx[matchedVtx].bestByTree[tree] =
    std::make_pair(key, edgeIdxInMatched);
  updateVertexCQLeaf(matchedVtx);
  pushD2Best(matchedVtx);
}

inline void FastMatchingEngineV2::refreshMatchedVertexD2(Size v)
{
  refreshMatchedVertexD2Impl(v, true);
}

inline void FastMatchingEngineV2::refreshMatchedVertexD2Impl(Size v,
                                                             bool pushRoot)
{
  if (v >= mNumVtxs) return;
  Size B = topLevel(v);
  if (B == cNullItm || mBlsm[B].status != BS_MATCHED) {
    clearVertexD2(v);
    return;
  }
  if (pushRoot) eraseD2BlossomEntry(B);

  std::map<Size, std::pair<Val, Size> > newBest;
  for (Size i = 0; i < mAdj[v].size(); ++i) {
    Size u = mAdj[v][i];
    if (mMate[v] == u) continue;
    Size Bu = topLevel(u);
    if (Bu == cNullItm || Bu == B || mBlsm[Bu].status != BS_EVEN) continue;
    Size tree = mBlsm[Bu].tree;
    Val key = d2Key(v, i);
    std::map<Size, std::pair<Val, Size> >::iterator it =
      newBest.find(tree);
    if (it == newBest.end() || key < it->second.first) {
      newBest[tree] = std::make_pair(key, i);
    }
  }

  bool same = (newBest.size() == mVtx[v].bestByTree.size());
  if (same) {
    std::map<Size, std::pair<Val, Size> >::const_iterator a =
      newBest.begin();
    std::map<Size, std::pair<Val, Size> >::const_iterator b =
      mVtx[v].bestByTree.begin();
    for (; a != newBest.end(); ++a, ++b) {
      if (a->first != b->first ||
          a->second.second != b->second.second ||
          a->second.first < b->second.first - EPS ||
          a->second.first > b->second.first + EPS) {
        same = false;
        break;
      }
    }
  }

  if (!same) {
    mVtx[v].bestByTree.swap(newBest);
  } else {
    mVtx[v].bestByTree.swap(newBest);
  }

  updateVertexCQLeaf(v);
  if (pushRoot) pushD2BestForBlossom(B);
}

inline void FastMatchingEngineV2::pushD1(Size v)
{
  // D1 keys are invariant while a vertex remains in an EVEN top blossom:
  // mDeltaSum + effY(v) collapses to the baked vertex pot plus the blossom
  // offset. Reinsertions are therefore charged to label/potential baking
  // transitions, and stale versions are discarded lazily.
  if (v >= mNumVtxs) return;
  Val key = mDeltaSum + effY(v);
  if (mVtx[v].d1Queued &&
      key >= mVtx[v].d1Key - EPS &&
      key <= mVtx[v].d1Key + EPS) {
    return;
  }
  ++mVtx[v].d1Version;
  mVtx[v].d1Key = key;
  mVtx[v].d1Queued = true;

  HeapEntryD1 e;
  e.key = key;
  e.vtx = v;
  e.version = mVtx[v].d1Version;
  mPQD1.push(e);
}

inline void FastMatchingEngineV2::pushD1ForBlossom(Size blossom)
{
  if (blossom >= mBlsm.size()) return;
  if (trivial(blossom)) {
    pushD1(blossom);
    return;
  }

  std::vector<Size> vtxs;
  MembershipView pvtxs = &mBlsm[blossom].vertices;
  if (pvtxs->empty()) {
    collectBlossomVertices(blossom, &vtxs);
    pvtxs = &vtxs;
  }
  // D1's nontrivial-blossom scan is amortized over blossom status changes:
  // d1MinValid is invalidated only when vertex potentials are baked during
  // grow/shrink/split transitions, not on every event-loop validation.
  if (!mBlsm[blossom].d1MinValid ||
      mBlsm[blossom].d1MinVtx == cNullItm ||
      mBlsm[blossom].d1MinVtx >= mNumVtxs) {
    Size best = cNullItm;
    for (Size i = 0; i < pvtxs->size(); ++i) {
      Size v = (*pvtxs)[i];
      if (best == cNullItm || mVtx[v].pot < mVtx[best].pot) best = v;
    }
    mBlsm[blossom].d1MinVtx = best;
    mBlsm[blossom].d1MinValid = (best != cNullItm);
  } else {
  }
  if (mBlsm[blossom].d1MinValid) pushD1(mBlsm[blossom].d1MinVtx);
}

inline void FastMatchingEngineV2::pushD2Best(Size matchedVtx)
{
  if (matchedVtx >= mNumVtxs) return;
  Size B = topLevel(matchedVtx);
  if (B == cNullItm || B >= mBlsm.size()) return;
  pushD2BestForBlossom(B);
}

inline void FastMatchingEngineV2::pushD2BestForBlossom(Size blossom)
{
#ifdef MBX_GABOW_SEARCH_AUDIT
  ++mAudD2BestPushes;
#endif
  // D2 is aggregated at the current matched top blossom: vertex leaves hold
  // per-tree best arcs, and the CQ root contributes one heap event.
  if (blossom >= mBlsm.size() || !mBlsm[blossom].alive ||
      mBlsm[blossom].parent != cNullItm ||
      mBlsm[blossom].status != BS_MATCHED ||
      mBlsm[blossom].cqRoot == cNullItm) {
    if (blossom < mBlsm.size()) eraseD2BlossomEntry(blossom);
    return;
  }

  ConcatenableQueue::Handle leaf = mCQ.minElem(mBlsm[blossom].cqRoot);
  if (leaf == cNullItm || leaf >= mCQ.poolSize()) {
    eraseD2BlossomEntry(blossom);
    return;
  }
  Size matchedVtx = mCQ.payload(leaf);
  if (matchedVtx >= mNumVtxs || topLevel(matchedVtx) != blossom) {
    eraseD2BlossomEntry(blossom);
    return;
  }

  Val best = infKey();
  Size bestEdge = cNullItm;
  Size bestTree = cNullItm;
  std::map<Size, std::pair<Val, Size> >::const_iterator it;
  for (it = mVtx[matchedVtx].bestByTree.begin();
       it != mVtx[matchedVtx].bestByTree.end(); ++it) {
    if (it->second.first < best) {
      best = it->second.first;
      bestEdge = it->second.second;
      bestTree = it->first;
    }
  }
  if (bestEdge == cNullItm) {
    eraseD2BlossomEntry(blossom);
    return;
  }

  if (mBlsm[blossom].d2Queued &&
      mBlsm[blossom].d2Vtx == matchedVtx &&
      mBlsm[blossom].d2Edge == bestEdge &&
      mBlsm[blossom].d2Tree == bestTree &&
      best >= mBlsm[blossom].d2Key - EPS &&
      best <= mBlsm[blossom].d2Key + EPS)
    return;

  eraseD2BlossomEntry(blossom);

  HeapEntryD2 e;
  e.key = best;
  e.matchedVtx = matchedVtx;
  e.edgeIdxInMatched = bestEdge;
  e.matchedBlsm = blossom;
  e.tree = bestTree;
  e.version = 0;
  mPQD2.push(e);
  mBlsm[blossom].d2Queued = true;
  mBlsm[blossom].d2Key = e.key;
  mBlsm[blossom].d2Vtx = e.matchedVtx;
  mBlsm[blossom].d2Edge = e.edgeIdxInMatched;
  mBlsm[blossom].d2Tree = e.tree;
}

inline void FastMatchingEngineV2::dropD2TreeCandidate(Size matchedVtx,
                                                      Size tree)
{
#ifdef MBX_GABOW_SEARCH_AUDIT
  ++mAudD2Drops;
#endif
  if (matchedVtx >= mNumVtxs) return;
  eraseD2Entry(matchedVtx);
  if (tree != cNullItm) mVtx[matchedVtx].bestByTree.erase(tree);
  updateVertexCQLeaf(matchedVtx);
  pushD2Best(matchedVtx);
}

inline void FastMatchingEngineV2::pushD3(Size u, Size edgeIdx)
{
  ++mProfPushD3CallsAll;  // counts calls for BOTH backends (compare push volume)
  if (mBackend != kGabow2018) { pushD3Impl(u, edgeIdx); return; }
  long double t0 = MBX_HOTNOW();
  ++mProfPushD3Count;
  pushD3Impl(u, edgeIdx);
  mProfPushD3Time += MBX_HOTNOW() - t0;
}

inline void FastMatchingEngineV2::pushD3Impl(Size u, Size edgeIdx)
{
  // D3 events are lazy: an edge is queued when endpoint labels make it an
  // EVEN-EVEN candidate, and stale entries are discarded/requeued only after
  // endpoint tree/blossom state or key changes.
  if (u >= mNumVtxs || edgeIdx >= mAdj[u].size()) return;
  Size v = mAdj[u][edgeIdx];
  if (v < u) {
    Size oldU = u;
    Size idx = reverseEdgeIdx(u, edgeIdx);
    if (idx == cNullItm) return;
    u = v;
    edgeIdx = idx;
    v = oldU;
  }
  Size Bu = topLevel(u);
  Size Bv = topLevel(v);
  if (Bu == cNullItm || Bv == cNullItm || Bu == Bv) return;
  BlossomStatus su = mBlsm[Bu].status;
  BlossomStatus sv = mBlsm[Bv].status;
  if (!(su == BS_EVEN && sv == BS_EVEN)) {
    return;
  }

  Val key = d3Key(u, edgeIdx);
#ifdef MBX_SHARED_D3_TIE_ORDER
  if (!std::isfinite(key))
    throw std::invalid_argument("shared D3 candidate cost must be finite");
#endif
  const Size occurrenceTag = gabowForestOccurrenceTag(u, v);

#ifdef MBX_GMG_D3_OWNERSHIP_EXPERIMENT
  if (mBackend == kGMG1986 && mD3Queued[u][edgeIdx] != 0) {
    if (mGmgD3Publication[u][edgeIdx].key == key) return;
    // Changed prices replace logical ownership without consulting Gabow
    // handles. A previous physical record stays available only for cleanup.
    mD3Queued[u][edgeIdx] = 0;
  }
#endif
  if (u < mD3Queued.size() && edgeIdx < mD3Queued[u].size() &&
      mD3Queued[u][edgeIdx] != 0) {
    bool sameOccurrence = true;
    Size queuedTag = 0;
    if (mGabowDeepForest) {
      sameOccurrence = false;
      if (u < mBmHandleU.size() && edgeIdx < mBmHandleU[u].size()) {
        if (mBM.edge_occurrence_tag_if(mBmHandleU[u][edgeIdx],
                                       mBmGenU[u][edgeIdx], &queuedTag) ||
            mBM.edge_occurrence_tag_if(mBmHandleV[u][edgeIdx],
                                       mBmGenV[u][edgeIdx], &queuedTag))
          sameOccurrence = (queuedTag == occurrenceTag);
      }
    }
    if (!sameOccurrence) {
      // Stage-3 BOUNDED CHURN REPAIR: an occurrence-stale tag does not by
      // itself invalidate the record as a SCHEDULER entry. The D3 consumer
      // validates against CURRENT matchbox state (topLevel/status/trees) and
      // re-checks the price at fire time (stale-cost repush), and BM
      // classification is allowed to be lazily stale until R0. So a record is
      // reconstructed ONLY when its stored price no longer equals the current
      // key; otherwise the physical records are KEPT (saving remove_edge x2 +
      // nca + make_edge x2 per re-entry touch). Original-edge identity,
      // event-time correctness and container invariants are unchanged; no
      // O(n)/O(m) work is added per augmentation; old occurrences do not
      // accumulate (still <= 2 records per adjacency slot).
      Val keptCost;
      bool priceCurrent = false;
      if (u < mBmHandleU.size() && edgeIdx < mBmHandleU[u].size()) {
        if (mBM.edge_cost_if(mBmHandleU[u][edgeIdx], mBmGenU[u][edgeIdx], &keptCost) ||
            mBM.edge_cost_if(mBmHandleV[u][edgeIdx], mBmGenV[u][edgeIdx], &keptCost))
          priceCurrent = (keptCost == key);
      }
      if (priceCurrent) {
        ++mProfForestStaleKept;
        if (mBackend == kGabow2018) ++mProfPushD3Early;
        return;                 // records stay queued; mD3Queued remains 1
      }
      ++mProfForestStaleQueued;
      if (mGabowD3Log)
        std::fprintf(stderr,
          "FOREST stale-queued %lu %lu oldTag=%lu nowTag=%lu\n",
          (unsigned long)u, (unsigned long)v,
          (unsigned long)queuedTag, (unsigned long)occurrenceTag);
      // Price-stale: retire both old back-edge records, clear the dedup bit,
      // and rebuild below using the endpoints' current ancestor chain.
      if (u < mBmHandleU.size() && edgeIdx < mBmHandleU[u].size()) {
        mBM.remove_edge_if(mBmHandleU[u][edgeIdx], mBmGenU[u][edgeIdx]);
        mBM.remove_edge_if(mBmHandleV[u][edgeIdx], mBmGenV[u][edgeIdx]);
        mBmHandleU[u][edgeIdx] = cNullItm;
        mBmHandleV[u][edgeIdx] = cNullItm;
      }
      setD3QueuedFlag(u, edgeIdx, (char)0);
    } else {
      if (mBackend == kGabow2018) ++mProfPushD3Early;
      return;
    }
  }

  if (mBackend == kGabow2018) {
    if(u>=mInT.size()||v>=mInT.size()||!mInT[u]||!mInT[v]||
       mBM.blossom_of(u)==cNullItm||mBM.blossom_of(v)==cNullItm)
      throw std::runtime_error("outer candidate is absent from monotone support");
    if (u >= mInT.size() || v >= mInT.size() || !mInT[u] || !mInT[v]) {
      // Still not coverable (out-of-range ids) — fall through to legacy push.
    } else {
      // Re-push multiplicity instrumentation: count how many times THIS
      // matchbox edge (u, edgeIdx) is make_edge'd over the whole run. The
      // paper's phased model would push each edge ~once per phase; our
      // continuous loop re-pushes across implicit phase boundaries. A high
      // average multiplicity here is the suspected reason dense+blossom-heavy
      // graphs (pdb1HYS) don't show Gabow's expected asymptotic win.
      if (u < mPushMult.size() && edgeIdx < mPushMult[u].size()) {
        if (mPushMult[u][edgeIdx] == 0) ++mProfDistinctEdgesPushed;
        ++mPushMult[u][edgeIdx];
        ++mProfTotalEdgePushes;
        if ((unsigned long long)mPushMult[u][edgeIdx] > mProfMaxPushMult)
          mProfMaxPushMult = mPushMult[u][edgeIdx];
      }
#ifdef MBX_GABOW_SEARCH_AUDIT
      audNoteIntervalPush(u, edgeIdx);   // interval-LOCAL multiplicity (separate)
#endif
      const bool crossTree =
        (mBlsm[Bu].tree != cNullItm && mBlsm[Bv].tree != cNullItm &&
         mBlsm[Bu].tree != mBlsm[Bv].tree);
      Size a = mAuxRoot;
      if (!crossTree) {
        long double t0 = MBX_HOTNOW();
        a = mAuxNCA.nca(u, v);
        ++mProfNcaCount;
        if (a != mAuxRoot) ++mProfNcaNonRoot;   // structural: real meeting point
#ifdef MBX_GABOW_SEARCH_AUDIT
        if (a == u || a == v) {
          ++mAudNcaAncestor;
          ++mAudR8Ancestor;
        }
#endif
        mProfNcaTime += MBX_HOTNOW() - t0;
      }
      if (mGabowD3Log)
        std::fprintf(stderr, "P %lu %lu a=%ld x=%d k=%.0f bm=%ld/%ld mbx=%ld/%ld\n",
          (unsigned long)(u < v ? u : v), (unsigned long)(u < v ? v : u),
          (long)a, (int)crossTree, (double)key,
          (long)mBM.blossom_of(u < v ? u : v), (long)mBM.blossom_of(u < v ? v : u),
          (long)topLevel(u < v ? u : v), (long)topLevel(u < v ? v : u));
      Size payload = u * mAdjStride + edgeIdx;
      BlossomMerging::EdgeHandle hU = cNullItm;
      BlossomMerging::EdgeHandle hV = cNullItm;
      if (u >= mNumVtxs + 1 || a >= mNumVtxs + 1 || v >= mNumVtxs + 1) {
        std::fprintf(stderr,
          "pushD3 oob: u=%lu v=%lu a=%lu N=%lu\n",
          (unsigned long)u, (unsigned long)v, (unsigned long)a,
          (unsigned long)(mNumVtxs + 1));
        std::abort();
      }
      // R8 (adopted): the paper's reduction of an original edge e = uv through
      // its supporting-forest NCA a produces TWO path-side back-edge records
      // (u,a) and (v,a) -- for SAME-TREE endpoints a = nca(u,v), for
      // DIFFERENT-TREE endpoints a = the dummy root rho, consistently. Both
      // records carry the same original-edge payload; the D3 consumer removes
      // both on consumption and mD3Queued guarantees at most one solver event
      // per original edge. A ZERO-LENGTH side (a == u or a == v: a is an
      // ancestor endpoint) is safely omitted because BM's make_edge requires
      // a strictly positive span (depth(v) > depth(w)); the surviving side
      // retains the edge's key/payload. (The former one-cross-tree-record
      // policy was an optimization, retired per the R8 work order.)
      if (u != a) {
        long double t0 = MBX_HOTNOW();
#ifdef MBX_SHARED_D3_TIE_ORDER
        hU = mBM.make_edge(u, a, key, payload, occurrenceTag, payload);
#else
        hU = mBM.make_edge(u, a, key, payload, occurrenceTag);
#endif
        ++mProfMakeEdgeCount;
        mProfMakeEdgeTime += MBX_HOTNOW() - t0;
      }
      if (v != a) {
        long double t0 = MBX_HOTNOW();
#ifdef MBX_SHARED_D3_TIE_ORDER
        hV = mBM.make_edge(v, a, key, payload, occurrenceTag, payload);
#else
        hV = mBM.make_edge(v, a, key, payload, occurrenceTag);
#endif
        ++mProfMakeEdgeCount;
        mProfMakeEdgeTime += MBX_HOTNOW() - t0;
      }
#ifdef MBX_GABOW_SEARCH_AUDIT
      {
        ++mAudR8Reduced;
        if (crossTree) ++mAudR8DummyRoot;
        const int recs = (u != a ? 1 : 0) + (v != a ? 1 : 0);
        if (recs >= 2) {
          ++mAudDoubleRec;
          ++mAudR8TwoSide;
        } else {
          ++mAudSingleRec;
          ++mAudR8OneSide;
        }
      }
#endif
      if (u < mBmHandleU.size() && edgeIdx < mBmHandleU[u].size()) {
        mBmHandleU[u][edgeIdx] = hU;
        mBmHandleV[u][edgeIdx] = hV;
        mBmGenU[u][edgeIdx] = mBM.edge_generation(hU);
        mBmGenV[u][edgeIdx] = mBM.edge_generation(hV);
      }
      if (u < mD3Queued.size() && edgeIdx < mD3Queued[u].size())
        mD3QueuedPopulation.set(mD3Queued, u, edgeIdx, (char)1);
      return;
    }
  }


  HeapEntryD3 e;
  e.key = key;
  e.u = u;
  e.edgeIdx = edgeIdx;
#ifdef MBX_GMG_D3_OWNERSHIP_EXPERIMENT
  e.publication = 0;
  if (mBackend == kGMG1986) {
    if (mGmgD3Publication[u][edgeIdx].version ==
        std::numeric_limits<unsigned long long>::max()) {
      abortRun("gmg-d3-publication-overflow");
      return;
    }
    e.publication = mGmgD3Publication[u][edgeIdx].version + 1;
  }
#endif
  mPQD3.push(e);
#ifdef MBX_GMG_D3_OWNERSHIP_EXPERIMENT
  if (mBackend == kGMG1986) {
    mGmgD3Publication[u][edgeIdx].key = key;
    mGmgD3Publication[u][edgeIdx].version = e.publication;
  }
#endif
#ifdef MBX_GABOW_SEARCH_AUDIT
  ++mAudGmgD3Pushes;
#endif
  if (u < mD3Queued.size() && edgeIdx < mD3Queued[u].size())
    setD3QueuedFlag(u, edgeIdx, (char)1);
}

inline void FastMatchingEngineV2::pushD4(Size blossom)
{
  // D4 has one live event for each top-level nontrivial ODD blossom exposure.
  // Its key is stable while the blossom remains ODD, so stale pops are charged
  // to split/parent/status changes rather than to every dual update.
  if (blossom >= mBlsm.size() || trivial(blossom) ||
      !mBlsm[blossom].alive || mBlsm[blossom].parent != cNullItm ||
      mBlsm[blossom].status != BS_ODD)
    return;
  HeapEntryD4 e;
  e.key = mDeltaSum + effZ(blossom) / 2.0;
  e.blsm = blossom;
  mPQD4.push(e);
}

inline void FastMatchingEngineV2::addTreeMember(Size blossom, Size tree)
{
  if (tree == cNullItm) return;
  if (tree >= mTreeMembers.size()) mTreeMembers.resize(tree + 1);
  if (mBlsm[blossom].tree != tree) mBlsm[blossom].tree = tree;
  if (mBlsm[blossom].listedTree == tree) return;
  mBlsm[blossom].listedTree = tree;
#ifdef MBX_GABOW_SEARCH_AUDIT
  ++mAudTreeAdds;
#endif
  const Size oldCapacity = mGabowLoopDiag ? mTreeMembers[tree].capacity() : 0;
  mTreeMembers[tree].push_back(blossom);
  const Size newCapacity = mGabowLoopDiag ? mTreeMembers[tree].capacity() : 0;
  if (mGabowLoopDiag && newCapacity > oldCapacity) {
    mProfTreeMemberCapacity += newCapacity - oldCapacity;
    if (newCapacity > mProfTreeMemberMaxRowCapacity)
      mProfTreeMemberMaxRowCapacity = newCapacity;
  }
}

inline void FastMatchingEngineV2::clearMate(Size vtx)
{
  if (vtx >= mNumVtxs) return;
  Size old = mMate[vtx];
  if (old != cNullItm && old < mNumVtxs && mMate[old] == vtx)
    mMate[old] = cNullItm;
  mMate[vtx] = cNullItm;
}

inline void FastMatchingEngineV2::setMatePair(Size u, Size v)
{
  if (u >= mNumVtxs || v >= mNumVtxs || u == v) return;
  Size oldU = mMate[u];
  Size oldV = mMate[v];
  if (mMaterializeDepth == 0) {
    invalidateMaterializedAncestors(u);
    invalidateMaterializedAncestors(v);
    if (oldU != cNullItm && oldU < mNumVtxs) invalidateMaterializedAncestors(oldU);
    if (oldV != cNullItm && oldV < mNumVtxs) invalidateMaterializedAncestors(oldV);
  }
  if (oldU != cNullItm && oldU < mNumVtxs && oldU != v &&
      mMate[oldU] == u) {
    mMate[oldU] = cNullItm;
    Size BoldU = topLevel(oldU);
    if (BoldU == oldU && BoldU < mBlsm.size()) mBlsm[BoldU].next = cNullItm;
  }
  if (oldV != cNullItm && oldV < mNumVtxs && oldV != u &&
      mMate[oldV] == v) {
    mMate[oldV] = cNullItm;
    Size BoldV = topLevel(oldV);
    if (BoldV == oldV && BoldV < mBlsm.size()) mBlsm[BoldV].next = cNullItm;
  }
  mMate[u] = v;
  mMate[v] = u;
  Size Bu = topLevel(u);
  Size Bv = topLevel(v);
  if (Bu == u && Bu < mBlsm.size()) mBlsm[Bu].next = encodeEdge(u, v, mNumVtxs);
  if (Bv == v && Bv < mBlsm.size()) mBlsm[Bv].next = encodeEdge(v, u, mNumVtxs);
}

inline Size FastMatchingEngineV2::edgeSrc(Size edge) const
{
  if (edge == cNullItm || mNumVtxs == 0) return cNullItm;
  return edge / mNumVtxs;
}

inline Size FastMatchingEngineV2::edgeDst(Size edge) const
{
  if (edge == cNullItm || mNumVtxs == 0) return cNullItm;
  return edge % mNumVtxs;
}

inline Size FastMatchingEngineV2::reverseEncodedEdge(Size edge) const
{
  if (edge == cNullItm) return cNullItm;
  return encodeEdge(edgeDst(edge), edgeSrc(edge), mNumVtxs);
}

inline bool FastMatchingEngineV2::blossomContainsVertex(Size blossom,
                                                        Size vtx) const
{
  if (blossom == cNullItm || blossom >= mBlsm.size() || vtx >= mNumVtxs)
    return false;
  if (trivial(blossom)) return blossom == vtx;
  if (!mBlsm[blossom].childOfVtx.empty()) {
    const std::vector<std::pair<Size, Size> >& cov = mBlsm[blossom].childOfVtx;
    std::vector<std::pair<Size, Size> >::const_iterator it =
      std::lower_bound(cov.begin(), cov.end(),
                       std::make_pair(vtx, (Size)0));
    return it != cov.end() && it->first == vtx;
  }
  // The leaf-to-root parent chain answers uncached membership without
  // scanning every vertex or recursively descending a deep blossom tree.
  Size current = vtx;
  for (Size steps = 0; steps < mBlsm.size() && current < mBlsm.size();
       ++steps) {
    if (current == blossom) return true;
    current = mBlsm[current].parent;
    if (current == cNullItm) break;
  }
  return false;
}

inline Size FastMatchingEngineV2::childIndexContaining(Size blossom,
                                                       Size vtx) const
{
  if (blossom == cNullItm || blossom >= mBlsm.size() || trivial(blossom))
    return cNullItm;
  if (!mBlsm[blossom].childOfVtx.empty()) {
    const std::vector<std::pair<Size, Size> >& cov = mBlsm[blossom].childOfVtx;
    std::vector<std::pair<Size, Size> >::const_iterator it =
      std::lower_bound(cov.begin(), cov.end(),
                       std::make_pair(vtx, (Size)0));
    if (it != cov.end() && it->first == vtx) return it->second;
    return cNullItm;
  }
  Size current = vtx;
  for (Size steps = 0; steps < mBlsm.size() && current < mBlsm.size();
       ++steps) {
    Size parent = mBlsm[current].parent;
    if (parent == blossom) {
      const std::vector<Size>& children = mBlsm[blossom].children;
      for (Size i = 0; i < children.size(); ++i)
        if (children[i] == current) return i;
      return cNullItm;
    }
    if (parent == cNullItm) break;
    current = parent;
  }
  return cNullItm;
}

inline void FastMatchingEngineV2::collectBlossomVertices(
  Size blossom, std::vector<Size>* vtxs) const
{
  if (vtxs == NULL || blossom == cNullItm || blossom >= mBlsm.size()) return;
  std::vector<Size> pending(1, blossom);
  while (!pending.empty()) {
    Size current = pending.back();
    pending.pop_back();
    if (current >= mBlsm.size()) continue;
    const BlossomData& B = mBlsm[current];
    if (!B.vertices.empty()) {
      vtxs->insert(vtxs->end(), B.vertices.begin(), B.vertices.end());
    } else if (trivial(current)) {
      vtxs->push_back(current);
    } else {
      // Reverse push order preserves the original recursive child order.
      for (Size i = B.children.size(); i > 0; --i)
        pending.push_back(B.children[i - 1]);
    }
  }
}

inline void FastMatchingEngineV2::setTopLevelForBlossom(Size blossom, Size top)
{
  (void)blossom;
  (void)top;
}

inline void FastMatchingEngineV2::clearBlossomD2(Size blossom)
{
  std::vector<Size> vtxs;
  MembershipView pvtxs = &mBlsm[blossom].vertices;
  if (pvtxs->empty()) {
    collectBlossomVertices(blossom, &vtxs);
    pvtxs = &vtxs;
  }
  for (Size i = 0; i < pvtxs->size(); ++i) clearVertexD2((*pvtxs)[i]);
}

inline void FastMatchingEngineV2::refreshMatchedBlossomD2(Size blossom)
{
  if (blossom >= mBlsm.size()) return;
  materializeBlossomForCurrentBoundary(blossom);
  eraseD2BlossomEntry(blossom);
  std::vector<Size> vtxs;
  MembershipView pvtxs = &mBlsm[blossom].vertices;
  if (pvtxs->empty()) {
    collectBlossomVertices(blossom, &vtxs);
    pvtxs = &vtxs;
  }
  // This whole-blossom scan is charged to exposing a top blossom as MATCHED
  // (or to split-created matched children). Per-vertex best arcs are stored in
  // CQ leaves; after this scan the blossom contributes one D2 root event.
  for (Size i = 0; i < pvtxs->size(); ++i)
    refreshMatchedVertexD2Impl((*pvtxs)[i], false);
  pushD2BestForBlossom(blossom);
}

inline void FastMatchingEngineV2::invalidateMaterializedAncestors(Size vtx)
{
  if (vtx >= mNumVtxs) return;
  Size b = topLevel(vtx);
  Size guard = 0;
  while (b != cNullItm && b < mBlsm.size() && guard <= mBlsm.size()) {
    mBlsm[b].materializedBase = cNullItm;
    b = mBlsm[b].parent;
    ++guard;
  }
}

inline Val FastMatchingEngineV2::originalEdgeWeight(Size u, Size v) const
{
  Size edgeIdx = findEdgeIdx(u, v);
  if (edgeIdx == cNullItm) return 0;
  return mWgt[u][edgeIdx];  // mWgt == original (cardTieEps=0, never mutated)
}

inline void FastMatchingEngineV2::refreshEvenBlossomEdges(Size blossom)
{
  std::vector<Size> vtxs;
  MembershipView pvtxs = &mBlsm[blossom].vertices;
  if (pvtxs->empty()) {
    collectBlossomVertices(blossom, &vtxs);
    pvtxs = &vtxs;
  }
  refreshEvenBlossomEdgesSubset(blossom, *pvtxs);
}

inline void FastMatchingEngineV2::refreshEvenBlossomEdgesSubset(
  Size blossom, MembershipView vtxs)
{
  Size tree = mBlsm[blossom].tree;
  ++mBestD2CurrentStamp;
  if (mBestD2CurrentStamp == cNullItm) {
    std::fill(mBestD2Stamp.begin(), mBestD2Stamp.end(), 0);
    mBestD2CurrentStamp = 1;
  }
  mBestD2Touched.clear();

  for (Size j = 0; j < vtxs.size(); ++j) {
    Size v = vtxs[j];
    if (mVtx[v].lastEvenScanTree == tree &&
        mVtx[v].lastEvenScanTop == blossom)
      continue;
    mVtx[v].lastEvenScanTree = tree;
    mVtx[v].lastEvenScanTop = blossom;
    clearVertexD2(v);
    for (Size i = 0; i < mAdj[v].size(); ++i) {
      Size w = mAdj[v][i];
      Size Bw = topLevel(w);
      if (Bw == cNullItm || Bw == blossom) continue;
      if (mBlsm[Bw].status == BS_EVEN) {
        pushD3(v, i);
      }
      else if (mBlsm[Bw].status == BS_MATCHED) {
        Size rev = reverseEdgeIdx(v, i);
        if (rev != cNullItm && mMate[w] != v) {
          Val key = d2Key(w, rev);
          if (mBestD2Stamp[w] != mBestD2CurrentStamp) {
            mBestD2Stamp[w] = mBestD2CurrentStamp;
            mBestD2Key[w] = key;
            mBestD2Edge[w] = rev;
            mBestD2Touched.push_back(w);
          } else if (key < mBestD2Key[w]) {
            mBestD2Key[w] = key;
            mBestD2Edge[w] = rev;
          }
        }
      }
    }
  }

  for (Size i = 0; i < mBestD2Touched.size(); ++i) {
    Size matchedVtx = mBestD2Touched[i];
    Size edgeIdx = mBestD2Edge[matchedVtx];
    if (matchedVtx >= mNumVtxs || edgeIdx >= mAdj[matchedVtx].size())
      continue;
    Size B = topLevel(matchedVtx);
    if (B == cNullItm || B >= mBlsm.size() ||
        mBlsm[B].status != BS_MATCHED)
      continue;

    std::map<Size, std::pair<Val, Size> >::iterator it =
      mVtx[matchedVtx].bestByTree.find(tree);
    if (it != mVtx[matchedVtx].bestByTree.end() &&
        it->second.first <= mBestD2Key[matchedVtx] + EPS)
      continue;

    eraseD2Entry(matchedVtx);
    mVtx[matchedVtx].bestByTree[tree] =
      std::make_pair(mBestD2Key[matchedVtx], edgeIdx);
    updateVertexCQLeaf(matchedVtx);
    pushD2Best(matchedVtx);
  }
}

inline void FastMatchingEngineV2::collectTreePath(
  Size blossom, std::vector<Size>* path) const
{
  if (path == NULL) return;
  path->clear();
  Size cur = blossom;
  while (cur != cNullItm && cur < mBlsm.size()) {
    path->push_back(cur);
    Size pred = mBlsm[cur].pred;
    if (pred == cNullItm) break;
    Size src = edgeSrc(pred);
    if (src >= mNumVtxs) break;
    Size parent = topLevel(src);
    if (parent == cNullItm || parent == cur) break;
    cur = parent;
  }
}

inline Size FastMatchingEngineV2::orientedTreeEdge(Size fromBlossom,
                                                   Size toBlossom) const
{
  if (fromBlossom == cNullItm || toBlossom == cNullItm ||
      fromBlossom >= mBlsm.size() || toBlossom >= mBlsm.size())
    return cNullItm;

  Size candidates[2];
  candidates[0] = mBlsm[toBlossom].pred;
  candidates[1] = mBlsm[fromBlossom].pred;
  for (int i = 0; i < 2; ++i) {
    Size e = candidates[i];
    if (e == cNullItm) continue;
    Size s = edgeSrc(e);
    Size t = edgeDst(e);
    if (s >= mNumVtxs || t >= mNumVtxs) continue;
    Size Bs = topLevel(s);
    Size Bt = topLevel(t);
    if (Bs == fromBlossom && Bt == toBlossom) return e;
    if (Bs == toBlossom && Bt == fromBlossom) return reverseEncodedEdge(e);
  }
  return cNullItm;
}

inline void FastMatchingEngineV2::rotateBlossomToChild(Size blossom,
                                                       Size childPos)
{
  if (blossom >= mBlsm.size() || trivial(blossom) ||
      childPos == 0 || childPos >= mBlsm[blossom].children.size()) {
    return;
  }

  Size child = mBlsm[blossom].children[childPos];
  if (child < mBlsm.size()) mBlsm[blossom].base = mBlsm[child].base;
}

inline void FastMatchingEngineV2::augmentThroughBlossom(Size blossom,
                                                        Size entryVtx)
{
  if (blossom == cNullItm || blossom >= mBlsm.size() || trivial(blossom))
    return;
  Size entry = childIndexContaining(blossom, entryVtx);
  if (entry == cNullItm) {
    Size edgeRefs = 0;
    for (Size i = 0; i < mAdj.size(); ++i) edgeRefs += mAdj[i].size();
    if (edgeRefs / 2 <= 50000) abortRun("augment-through-entry-missing");
    return;
  }

  Size k = mBlsm[blossom].children.size();
  if (k == 0 || mBlsm[blossom].childEdges.size() != k) return;

  Size baseIdx = childIndexContaining(blossom, mBlsm[blossom].base);
  if (baseIdx == cNullItm || baseIdx >= k) baseIdx = 0;

  if (entry != baseIdx) mBlsm[blossom].materializedBase = cNullItm;
  Size cur = entry;
  while (cur != baseIdx && !mAbortRun) {
    Size bxIdx, byIdx, edge;
    Size rel = (cur + k - baseIdx) % k;
    if ((rel % 2) == 0) {
      byIdx = (cur + k - 1) % k;
      bxIdx = (cur + k - 2) % k;
      edge = mBlsm[blossom].childEdges[bxIdx];
      cur = bxIdx;
    } else {
      Size edgeIdx = (cur + 1) % k;
      bxIdx = edgeIdx;
      byIdx = (edgeIdx + 1) % k;
      edge = mBlsm[blossom].childEdges[edgeIdx];
      cur = byIdx;
    }

    Size x = edgeSrc(edge);
    Size y = edgeDst(edge);
    if (x >= mNumVtxs || y >= mNumVtxs) {
      abortRun("augment-through-edge-invalid");
      return;
    }
    Size bx = mBlsm[blossom].children[bxIdx];
    Size by = mBlsm[blossom].children[byIdx];
    augmentThroughBlossom(bx, x);
    augmentThroughBlossom(by, y);
    setMatePair(x, y);
    if (bx < mBlsm.size()) mBlsm[bx].next = edge;
    if (by < mBlsm.size()) mBlsm[by].next = reverseEncodedEdge(edge);
  }

  rotateBlossomToChild(blossom, entry);
  materializeBlossomMatching(blossom, entryVtx);
}

inline void FastMatchingEngineV2::materializeBlossomForCurrentBoundary(
  Size blossom)
{
  if (blossom == cNullItm || blossom >= mBlsm.size() || trivial(blossom))
    return;

  Size base = mBlsm[blossom].base;
  if (mBlsm[blossom].materializedBase == base) {
    return;
  }
  materializeBlossomMatching(blossom, base);
}

inline void FastMatchingEngineV2::materializeBlossomMatching(Size blossom,
                                                             Size baseVtx)
{
  // Memoization and child visit order match the recursive version. Explicit
  // frames are necessary because a deep nested blossom can exceed the native
  // call stack even when its membership buffers fit comfortably in memory.
  struct Frame {
    Size blossom, baseVtx, k, ib, step, phase, edge;
    Frame(Size b, Size v)
      : blossom(b), baseVtx(v), k(0), ib(0), step(0), phase(0),
        edge(cNullItm) {}
  };
  std::vector<Frame> frames;
  frames.push_back(Frame(blossom, baseVtx));
  while (!frames.empty()) {
    Frame& f = frames.back();
    if (f.phase == 0) {
      if (f.blossom == cNullItm || f.blossom >= mBlsm.size() ||
          trivial(f.blossom) ||
          mBlsm[f.blossom].materializedBase == f.baseVtx) {
        frames.pop_back();
        continue;
      }
      f.k = mBlsm[f.blossom].children.size();
      if (f.k == 0 || mBlsm[f.blossom].childEdges.size() != f.k) {
        frames.pop_back();
        continue;
      }
      ++mMaterializeDepth;
      f.ib = childIndexContaining(f.blossom, f.baseVtx);
      if (f.ib == cNullItm || f.ib >= f.k) f.ib = 0;
      f.step = 1;
      f.phase = 1;
      const Size child = mBlsm[f.blossom].children[f.ib];
      frames.push_back(Frame(child, f.baseVtx));
      continue;
    }
    if (f.phase == 1) {
      if (f.step >= f.k) {
        const Size baseChild = mBlsm[f.blossom].children[f.ib];
        if (baseChild < mBlsm.size())
          mBlsm[f.blossom].base = mBlsm[baseChild].base;
        mBlsm[f.blossom].materializedBase = f.baseVtx;
        --mMaterializeDepth;
        frames.pop_back();
        continue;
      }
      const Size si = (f.ib + f.step) % f.k;
      f.edge = mBlsm[f.blossom].childEdges[si];
      if (f.edge == cNullItm || edgeSrc(f.edge) >= mNumVtxs ||
          edgeDst(f.edge) >= mNumVtxs) {
        f.step += 2;
        continue;
      }
      const Size first = mBlsm[f.blossom].children[si];
      const Size x = edgeSrc(f.edge);
      f.phase = 2;
      frames.push_back(Frame(first, x));
      continue;
    }
    if (f.phase == 2) {
      const Size ti = (f.ib + f.step + 1) % f.k;
      const Size second = mBlsm[f.blossom].children[ti];
      const Size y = edgeDst(f.edge);
      f.phase = 3;
      frames.push_back(Frame(second, y));
      continue;
    }
    const Size si = (f.ib + f.step) % f.k;
    const Size ti = (f.ib + f.step + 1) % f.k;
    const Size bs = mBlsm[f.blossom].children[si];
    const Size bt = mBlsm[f.blossom].children[ti];
    const Size x = edgeSrc(f.edge);
    const Size y = edgeDst(f.edge);
    setMatePair(x, y);
    if (bs < mBlsm.size()) mBlsm[bs].next = f.edge;
    if (bt < mBlsm.size()) mBlsm[bt].next = reverseEncodedEdge(f.edge);
    f.step += 2;
    f.phase = 1;
  }
}

inline void FastMatchingEngineV2::flipPathFromTopEven(Size evenBlossom,
                                                      Size entryVtx)
{
  Size evenB = evenBlossom;
  Size entry = entryVtx;
  while (evenB != cNullItm && evenB < mBlsm.size()) {
    augmentThroughBlossom(evenB, entry);
    Size evenPred = mBlsm[evenB].pred;
    if (evenPred == cNullItm) break;

    Size oddV = edgeSrc(evenPred);
    if (oddV >= mNumVtxs) break;
    Size oddB = topLevel(oddV);
    if (oddB == cNullItm || oddB >= mBlsm.size()) break;

    Size oddPred = mBlsm[oddB].pred;
    if (oddPred == cNullItm) break;
    Size prevEvenV = edgeSrc(oddPred);
    Size oddEntryV = edgeDst(oddPred);
    if (prevEvenV >= mNumVtxs || oddEntryV >= mNumVtxs) break;

    augmentThroughBlossom(oddB, oddEntryV);
    setMatePair(prevEvenV, oddEntryV);
    mBlsm[oddB].next = reverseEncodedEdge(oddPred);
    Size prevEvenB = topLevel(prevEvenV);
    if (prevEvenB != cNullItm && prevEvenB < mBlsm.size())
      mBlsm[prevEvenB].next = oddPred;

    evenB = prevEvenB;
    entry = prevEvenV;
  }
}

inline void FastMatchingEngineV2::flipPathFromEven(Size evenVtx)
{
  if (evenVtx >= mNumVtxs) return;
  Size B = topLevel(evenVtx);
  flipPathFromTopEven(B, evenVtx);
}

inline Size FastMatchingEngineV2::findEdgeIdx(Size u, Size v) const
{
  if (u >= mNumVtxs) return cNullItm;
  for (Size i = 0; i < mAdj[u].size(); ++i) {
    if (mAdj[u][i] == v) return i;
  }
  return cNullItm;
}

inline Size FastMatchingEngineV2::reverseEdgeIdx(Size u, Size edgeIdx) const
{
  if (u >= mRevIdx.size() || edgeIdx >= mRevIdx[u].size())
    return cNullItm;
  return mRevIdx[u][edgeIdx];
}

template<class F>
void FastMatchingEngineV2::forEachVertex(Size b, F& f) const
{
  if (trivial(b)) {
    f(b);
    return;
  }
  for (Size i = 0; i < mBlsm[b].children.size(); ++i)
    forEachVertex(mBlsm[b].children[i], f);
}

// ============================================================
// GeneralizedPQ1: A priority queue supporting O(1) uniform
// decrease of all elements' priorities.
//
// Operations (all O(log n) except decreaseAll which is O(1)):
//   insert(id, priority)  — add element with given priority
//   deleteMin()           — remove and return the element with minimum priority
//   findMin()             — return (id, effective_priority) of minimum element
//   decreaseAll(delta)    — subtract delta from all elements' priorities
//   empty() / size()      — query state
//
// Implementation: binary min-heap storing (stored_priority, id).
// A global offset tracks cumulative decreases.
// stored_priority = real_priority + offset at time of insertion.
// effective_priority = stored_priority - current_offset.
//
// Pre-C++11 compatible.
// ============================================================
class GeneralizedPQ1 {
public:
  GeneralizedPQ1() : mOffset(0) {}

  bool empty() const { return mHeap.empty(); }
  Size size() const { return mHeap.size(); }

  // Insert element with given effective priority
  void insert(Size id, Val priority)
  {
    // Store priority + offset so that effective = stored - offset = priority
    mHeap.push_back(std::make_pair(priority + mOffset, id));
    siftUp(mHeap.size() - 1);
  }

  // Return (id, effective_priority) of minimum element
  std::pair<Size, Val> findMin() const
  {
    return std::make_pair(mHeap[0].second, mHeap[0].first - mOffset);
  }

  // Remove and return (id, effective_priority) of minimum element
  std::pair<Size, Val> deleteMin()
  {
    std::pair<Size, Val> result = findMin();
    // Move last element to root and sift down
    mHeap[0] = mHeap.back();
    mHeap.pop_back();
    if (!mHeap.empty())
      siftDown(0);
    return result;
  }

  // Decrease all elements' effective priorities by delta.
  // O(1) operation.
  void decreaseAll(Val delta)
  {
    mOffset += delta;
  }

  // Clear the priority queue
  void clear()
  {
    mHeap.clear();
    mOffset = 0;
  }

private:
  std::vector<std::pair<Val, Size> > mHeap; // (stored_priority, id)
  Val mOffset;

  void siftUp(Size idx)
  {
    while (idx > 0)
    {
      Size parent = (idx - 1) / 2;
      if (mHeap[idx].first < mHeap[parent].first)
      {
        std::swap(mHeap[idx], mHeap[parent]);
        idx = parent;
      }
      else break;
    }
  }

  void siftDown(Size idx)
  {
    Size n = mHeap.size();
    while (true)
    {
      Size smallest = idx;
      Size left = 2 * idx + 1;
      Size right = 2 * idx + 2;
      if (left < n && mHeap[left].first < mHeap[smallest].first)
        smallest = left;
      if (right < n && mHeap[right].first < mHeap[smallest].first)
        smallest = right;
      if (smallest != idx)
      {
        std::swap(mHeap[idx], mHeap[smallest]);
        idx = smallest;
      }
      else break;
    }
  }
};

class MatchingEngine {
public:
  /// Constructor.
  MatchingEngine():
      mInlz(false), mQueAndStkType(eVecQueAndStk), mIdxdQueType(eLstIdxdQue),
      mIdxdPriQueType(eBhpIdxdPriQue), mCardGraphSrchType(eSglSrcBfs),
      mEdgWghtGraphSrchType(eSglSrcBfs), mVtxWghtAlgType(eDrct),
      mPrcmptdMaxCard(false), mPrecision(cMaxPrecision), mFullPrint(false),
      mStatsPrint(false) {}

  /// Destructor.
  ~MatchingEngine() {}

  ///
  void ComputeMxmlMatching(const Graph& graph,
                           std::vector<Size>* mateVec,
                           Size* card) const;

  ///
  void ComputeMaxCardMatching(const Graph& graph,
                              std::vector<Size>* mateVec,
                              Size* card) const;

  ///
  void ComputeAprxMaxEdgWghtMatching4(const Graph& graph,
                                      std::vector<Size>* mateVec,
                                      Size* card,
                                      Val* edgWght) const;
  ///
  void ComputeMaxEdgWghtMatchingEdmonds(const Graph& graph,
                                       std::vector<Size>* mateVec,
                                       Size* card,
                                       Val* edgWght) const;
  ///
  void ComputeMaxEdgWghtMatchingEdmondsFastV2(const Graph& graph,
                                             std::vector<Size>* mateVec,
                                             Size* card,
                                             Val* edgWght) const;
  ///
  void ComputeMaxEdgWghtMatchingGabow(const Graph& graph,
                                      std::vector<Size>* mateVec,
                                      Size* card,
                                      Val* edgWght) const;
  ///
  void ComputeMaxVtxWghtMatching(const Graph& graph,
                                      std::vector<Size>* mateVec,
                                      Size* card,
                                      Val* vtxWght, Size augPathLenBnd) const;

  ///
  void ComputeHalfVtxWghtMatching(const Graph& graph,
                                      std::vector<Size>* mateVec,
                                      Size* card,
                                      Val* vtxWght) const;
  ///
  void ComputeTwoThirdVtxWghtMatching(const Graph& graph,
                                      std::vector<Size>* mateVec,
                                      Size* card,
                                      Val* vtxWght) const;
  /////
  void SuitorVtxWghtMatching(const Graph& graph,
     							std::vector<Size>* mateVec, 
								 Size* card, 
								 Val* vtxWght) const;
  ///
  void ParSuitorVtxWghtMatching(const Graph& graph,
     								std::vector<Size>* mateVec, 
									 Size* card, 
									 Val* vtxWght) const;
  ///
  void PGVtxWghtMatching(const Graph& graph,
     						std::vector<Size>* mateVec,
							Size* card, 
							Val* vtxWght) const;  
	///
  void PGDPVtxWght(const Graph& graph,
     					std::vector<Size>* mateVec, 
						Size* card, 
						Val* vtxWght) const;
	////
	void LocallyDominantVtxWght(const Graph& graph, 
									std::vector<Size>* mateVec, 
									Size* card, 
									Val* vtxWght)  const;
	///
	void ScaleOneMinusEpsVtxWght(const Graph& graph,
     							std::vector<Size>* mateVec, 
								Size* card, 
								Val* vtxWght, 
								Val ep) const;
     ////
	void ComputeHalfEdgWghtMatching(const Graph& graph, 
									std::vector<Size>* mateVec, 
									Size* card, 
									Val* edgWght) const;
  //////////
  Err GetMatchedEdgWght(const Graph& graph,
                        const std::vector<Size>& mateVec,
                        Size card,
                        Val *edgWght) const;

  ///
  Err GetMatchedVtxWght(const Graph& graph,
                        const std::vector<Size>& mateVec,
                        Size card,
                        Val *vtxWght,
                        Val *edgWght) const;

  ///
  void SetInlz(bool inlz) { mInlz = inlz; }

  ///
  Err SetQueAndStkType(QueAndStkType queAndStkType) {
    switch (queAndStkType) {
    case eVecQueAndStk:
    case eLstQueAndStk:
      mQueAndStkType = queAndStkType;
      return eErrNone;
    default:
      return eErrInvType;
    }
  }

  ///
  Err SetIdxdQueType(IdxdQueType idxdQueType) {
    switch (idxdQueType) {
    case eLstIdxdQue:
      mIdxdQueType = idxdQueType;
      return eErrNone;
    default:
      return eErrInvType;
    }
  }

  ///
  Err SetIdxdPriQueType(IdxdPriQueType idxdPriQueType) {
    switch (idxdPriQueType) {
    case eLstIdxdPriQue:
    case eBhpIdxdPriQue:
      mIdxdPriQueType = idxdPriQueType;
      return eErrNone;
    case eFhpIdxdPriQue:
      // TODO: add fibonacci heap indexed priority queue data structure.
    default:
      return eErrInvType;
    }
  }

  ///
  Err SetCardGraphSrchType(GraphSrchType cardGraphSrchType) {
    switch (cardGraphSrchType) {
    case eSglSrcBfs:
      mCardGraphSrchType = cardGraphSrchType;
      return eErrNone;
    default:
      // TODO: add other cardinality computations, especially multiple path.
      return eErrInvType;
    }
  }

  ///
  Err SetEdgWghtGraphSrchType(GraphSrchType edgWghtGraphSrchType) {
    switch (edgWghtGraphSrchType) {
    case eSglSrcBfs:
      mEdgWghtGraphSrchType = edgWghtGraphSrchType;
      return eErrNone;
    case eMplSrcSglPath:
      // TODO: add multiple source edge-weight computations.
    default:
      // TODO: reformulate the comment below in the context of more algorithms.
      // multiple path edge-weight computations are not implemented since they
      // are generally not expected to improve the performance, given the fact
      // that it is generally unlikely to find a significant number of
      // augmenting paths in the tight graph during a particular iteration.
      return eErrInvType;
    }
  }

  ///
  Err SetVtxWghtAlgType(VtxWghtAlgType vtxWghtAlgType) {
    switch (vtxWghtAlgType) {
    case eDrct:
    case eSpcl:
      mVtxWghtAlgType = vtxWghtAlgType;
      return eErrNone;
    default:
      return eErrInvType;
    }
  }

  ///
  void SetPrcmptdMaxCard(bool prcmptdMaxCard) {
    mPrcmptdMaxCard = prcmptdMaxCard;
  }

  ///
  Err SetPrecision(Size precision) {
    if (precision > cMaxPrecision) {
      return eErrInvPrecision;
    }
    mPrecision = precision;
    return eErrNone;
  }

  ///
  void SetFullPrint(bool fullPrint) { mFullPrint = fullPrint; }

  ///
  void SetStatsPrint(bool statsPrint) { mStatsPrint = statsPrint; }

  ///
  bool CheckMatching(const Graph& graph,
                     const std::vector<Size>& mateVec,
                     Size card) const;

  ///
  void PrintMatching(const std::vector<Size>& mateVec,
                     Size card) const;

  ///
  void PrintMatching(const std::vector<Size>& mateVec,
                     Size card,
                     Val wght) const;

private:
  // controls the greedy initialization of a matching: if false then do not
  // perform greedy initialization, if true then perform greedy initialization.
  bool mInlz;
  // type of queue and stack used for vertex processing.
  QueAndStkType mQueAndStkType;
  // type of indexed queue used for vertex and weight processing.
  IdxdQueType mIdxdQueType;
  // type of indexed priority queue used for vertex and weight processing.
  IdxdPriQueType mIdxdPriQueType;
  // type of graph search used for the computation of a maximum cardinality
  // matching.
  GraphSrchType mCardGraphSrchType;
  // type of graph search used for the computation of a maximum edge-weight
  // matching.
  GraphSrchType mEdgWghtGraphSrchType;
  // type of algorithm used for the computation of a maximum vertex-weight
  // matching.
  VtxWghtAlgType mVtxWghtAlgType;
  // TODO: add comment that explains why mutable.
  // controls the precomputation of a maximum cardinality matching for the
  // speculative computation of a maximum vertex-weight matching.
  mutable bool mPrcmptdMaxCard;
  // precision used for printing floating point numbers.
  Size mPrecision;
  // controls the amount of printing: if false then print less information,
  // if true then print more information.
  bool mFullPrint;
  // controls the printing of computation statistics: if false then do not
  // print computation statistics, if true then print computation statistics.
  bool mStatsPrint;

  // forbid the copy constructor and the assignment operator.
  MatchingEngine(const MatchingEngine&);
  MatchingEngine& operator=(const MatchingEngine&);

  ///
  void rInlzGrdyForCard(const Graph& graph,
                        Size* mateArr,
                        Size* card) const;
  ///
  void rAugment(Size* mateArr,
                Size* ptr1Arr,
                Size* ptr2Arr,
                Size sLast,
                Size tLast) const;
  ///       
  void rAugment(Size* mateArr,
                Size* ptr1Arr,
                Size* ptr2Arr,
                Size sLast,
                Size tLast,
				Val* auglen) const;
				////////      
  void fix(Size* mateArr,
                Size* ptr1Arr,
                Size* ptr2Arr,
                Size sLast,
                Size tLast) const;
  ///
  Size rFind(Size* linkArr,
             Size u) const;

  ///
  void rUnion(Size* linkArr,
              Size* rankArr,
              Size u,
              Size v) const;

  ///
  void rLink(Size* linkArr,
             Size* rankArr,
             Size u,
             Size v) const;

  ///
  /// without disjoint set support.
  template<class ItmQue>
  void rProcessBlsm1(const Graph& graph,
                     Size* mateArr,
                     Size* ptr1Arr,
                     Size* ptr2Arr,
                     Size* blsmArr,
                     Size* rprsArr,
                     Stt* sttArr,
                     ItmQue& prcsbQue,
                     Size s1,
                     Size s2,
                     Size blsm) const;

  ///
  /// witho disjoint set support.
  template<class ItmQue>
  void rProcessBlsm2(const Graph& graph,
                     Size* mateArr,
                     Size* ptr1Arr,
                     Size* ptr2Arr,
                     Size* blsmArr,
                     Size* rprsArr,
                     Size* linkArr,
                     Size* rankArr,
                     Stt* sttArr,
                     ItmQue& prcsbQue,
                     Size s1,
                     Size s2,
                     Size blsm) const;
///
///used by speculative MVM  
template<class ItmQue>
void rProcessBlsm6(const Graph& graph, 
				Size* mateArr,
    			Size* ptr1Arr, 
				Size* ptr2Arr, 
				Size* blsmArr, 
				Size* rprsArr, 
				Size* linkArr,
    			Size* rankArr, 
				Stt* sttArr, 
				ItmQue& prcsbQue, 
				Size s1, 
				Size s2, 
				Size blsm, 
				Val* minVtxWght, 
				Size* sLast, 
				Size* tLast) const;
	////
  ///
  /// implicitly breadth-first.
  /// without disjoint set support.
  ///
  /// used in:
  /// - maximum cardinality computations.
  template<class ItmQue>
  void rFindAugPathCardSglSrc1(const Graph& graph,
                               Size* mateArr,
                               Size* ptr1Arr,
                               Size* ptr2Arr,
                               Size* blsmArr,
                               Size* rprsArr,
                               Stt* sttArr,
                               ItmQue& prcsbQue,
                               ItmQue& prcsdQue,
                               Size sFirst,
                               Size* sLast,
                               Size* tLast) const;

  ///
  /// implicitly breadth-first.
  /// with disjoint set support.
  ///
  /// used in:
  /// - maximum cardinality computations.
  template<class ItmQue>
  void rFindAugPathCardSglSrc2(const Graph& graph,
                               Size* mateArr,
                               Size* ptr1Arr,
                               Size* ptr2Arr,
                               Size* blsmArr,
                               Size* rprsArr,
                               Size* linkArr,
                               Size* rankArr,
                               Stt* sttArr,
                               ItmQue& prcsbQue,
                               ItmQue& prcsdQue,
                               Size sFirst,
                               Size* sLast,
                               Size* tLast) const;

  ///
  ///used by MVM  
	template<class ItmQue>
	void rFindAugPathCardSglSrc3(const Graph& graph, 
								Size* mateArr,
    							Size* ptr1Arr, 
								Size* ptr2Arr, 
								Size* blsmArr, 
								Size* rprsArr, 
								Size* linkArr,
    							Size* rankArr, 
								Stt* sttArr, 
								ItmQue& prcsbQue, 
								ItmQue& prcsdQue,
    							Size sFirst, 
								Size* sLast, 
								Size* tLast,
								Val w2)const;
    ///
	///
	///used in restricted augmenting paths MVM  
  	template<class ItmQue>
	void rFindAugPathCardSglSrc1Bnd(const Graph& graph, 
									Size* mateArr,
    								Size* ptr1Arr, 
									Size* ptr2Arr, 
									Size* rprsArr, 
									Stt* sttArr, 
									Size* sLenArr,
    								ItmQue& prcsbQue, 
									ItmQue& prcsdQue, 
									Size sFirst, 
									Size* sLast, 
									Size* tLast,
									Size augPathLenBnd)const;
  ///
  /// the core of the single source maximum cardinality matching computation.
  /// without disjoint set support.
  template<class ItmQue>
  void rComputeMaxCardMatchingSglSrc1(const Graph& graph,
                                      Size* mateArr,
                                      Size* card) const;

  ///
  /// the core of the single source maximum cardinality matching computation.
  /// with disjoint set support.
  template<class ItmQue>
  void rComputeMaxCardMatchingSglSrc2(const Graph& graph,
                                      Size* mateArr,
                                      Size* card) const;
    ///  
	/// 
	/// the core of the single source maximum cardinality matching computation.                              
	template<class ItmQue>
  	void rComputeMaxCardMatchingSglSrc3(const Graph& graph,
                                      Size* mateArr,
                                      Size* card,
									  std::pair<Size,Val>* sortedArr) const;
  ///
  /// the core of the single source approximate maximum weighted matching computation.
  template<class Que>
  void rComputeAprxMaxEdgWghtMatching4(const Graph& graph,
                                       Size* mateArr,
                                       Size* card,
                                       Val* edgWght) const;
  ///
  /// the core of the Edmonds primal-dual exact maximum edge-weight matching computation.
  template<class ItmQue, class ItmStk, class ItmIdxdQue>
  void rComputeMaxEdgWghtMatchingEdmonds(const Graph& graph,
                                        Size* mateArr,
                                        Size* card,
                                        Val* edgWght) const;
  ///
  ///Find increasing paths it used in Specluative MVM                                     
  template<class ItmQue> void rFindIncPathCardSglSrc3(const Graph& graph, 
  													Size* mateArr,
    												Size* ptr1Arr, 
													Size* ptr2Arr, 
													Size* blsmArr, 
													Size* rprsArr, 
													Size* linkArr,
    												Size* rankArr, 
													Stt* sttArr, 
													ItmQue& prcsbQue, 
													ItmQue& prcsdQue,
    												Size sFirst, 
													Size* sLast, 
													Size* tLast) const;
	////
	/// the core of the single source maximum vertex weighted matching computation.
  template<class ItmQue, class ItmStk>
  void rComputeMaxVtxWght(const Graph& graph,
    					  Size* mateArr, 
						  Size* card, 
						  Val* vtxWght, 
						  Size augPathLenBnd) const;
    ///
    /// the core of the single source greedy 1/2-approx. vertex weighted matching computation.
	template<class ItmQue>
  	void rComputeHalfVtxWghtMatching(const Graph& graph,
    					  			Size* mateArr, 
						  			Size* card, 
						  			Val* vtxWght) const;
    ///
    /// the core of the single source greedy 2/3-approx. vertex weighted matching computation.
   	template<class ItmQue>
  	void rComputeTwoThirdVtxWghtMatching(const Graph& graph,
    					  				Size* mateArr, 
						  				Size* card, 
						  				Val* vtxWght) const;
    ///
    /// the core of the single source Suitor 1/2-approx. vertex weighted matching computation.
    template<class ItmQue>
    void rSuitorVtxWght(const Graph& graph,
    					   Size* mateArr, 
						   Size* card, 
						   Val* vtxWght) const;
  ///
 	template<class ItmQue>
 	/// the core of the single source locally dominant 1/2-approx. vertex weighted matching computation.
	void rLocallyDominantVtxWght(const Graph& graph,
    								Size* mateArr, 
									Size* card, 
									Val* vtxWght) const;
  	///
  	/// the core of the shared memory parallel Suitor 1/2-approx. vertex weighted matching computation.
  	template<class ItmQue>
	void rParSuitorVtxWght(const Graph& graph,
    						Size* mateArr, 
							Size* card, 
							Val* vtxWght) const;
  	///
  	/// the core of the single source path growing 1/2-approx. vertex weighted matching computation.
  	template<class ItmQue>
   	void rPGPVtxWght(const Graph& graph,
    					Size* mateArr, 
						Size* card, 
						Val* vtxWght) const;
    ////
    /// the core of the single source path growing and DP 1/2-approx. vertex weighted matching computation.
	template<class ItmQue>
	void rPGDPVtxWght(const Graph& graph,
    					Size* mateArr, 
						Size* card, 
						Val* vtxWght) const;
	///
	///Scaling approx. algorithm functions
	/// the core of multiple source Scaling 1-eps -approx. vertex weighted matching computation.
	template<class ItmQue, class ItmStk, class ItmIdxdQue>
	void rScaleOneMinusEpsVtxWght(const Graph& graph, 
								Size* mateArr, 
								Size* card, 
								Val* vtxWght, 
								Val ep) const;
    ///
    // augment augemnting path 
    // used in Scaling approximation algorithm
    template<class ItmIdxdQue,class ItmQue>
  	void rAugment(Size* mateArr, 
	  			Size* ptrArr, 
				ItmIdxdQue& expsdQue, 
				Size sLast,
                Size tLast,
				Size* blsmArr, 
				Size* rprsArr, 
				Size* blsmParent, 
				std::vector<Size>* blsmChildren,
				std::vector<Size>* blsmL,
				Size num, 
				Size* label,
				Size* blsmin) const;
	///
	// augment through blossoms
    // used in Scaling approximation algorithm
	template<class ItmQue>
	void augmentBlsm(Size b,
					Size s,
					Size* mateArr,
					Size* blsmArr, 
					Size* rprsArr, 
					Size* blsmParent, 
					std::vector<Size>* blsmChildren,
					std::vector<Size>* blsmL,
					Size num) const;
	///
	// find which blossom a vertex belong to
    // used in Scaling approximation algorithm
	template<class ItmQue>	
	Size rFindblsm(Size* p,
					Size* blsm, 
					Size* rprs, 
					Size v) const;
	///
	// dissolve a blossom
    // used in Scaling approximation algorithm
	template<class ItmQue>
	void rDissolveBlsm3(const Graph& graph, 
						Size* blsmArr,
						Size* blsmParentArr, 
						std::vector<Size>* blsmChildren, 
						Size b) const;
	///
	inline bool rAcquireBlsmId(Size* blsm,
						Size numVtxs,
						Size maxBlsm,
						std::vector<Size>* blsmChildren,
						std::vector<Size>* blsmL,
						Size* blsmin,
						std::stack<Size>* freeBlsmStk,
						Size* b) const;
	///
	// dissolve a blossom and recursively reset blsmArr for ALL contained vertices.
	// Also clears reusable blossom state and pushes freed IDs onto freeBlsmStk.
	template<class ItmQue>
	void rDissolveBlsmRecursive(const Graph& graph,
						Size* blsmArr,
						Size* blsmParentArr,
						std::vector<Size>* blsmChildren,
						std::vector<Size>* blsmL,
						Size* rprsArr,
						Size* blsmin,
						Size b,
						std::stack<Size>& freeBlsmStk) const;
	///
	// create and process a blossom
    // used in Scaling approximation algorithm
	template<class ItmQue,class ItmStk>
	void rProcessBlsm3(const Graph& graph, 
						Size* mateArr,
    					Size* ptr1Arr, 
						Size* blsmArr, 
						Size* rprsArr, 
						Size* blsmParent, 
						std::vector<Size>* blsmChildren,
						std::vector<Size>* blsmL,
						Size* label,
    					Stt* sttArr, 
						ItmStk& prcsbStk, 
						Size s1, 
						Size s2, 
						Size* blsm, 
						Val* dualArr,
						Size* blsmin, 
						Size* treeNumArr,
						Size maxBlsm,
						std::stack<Size>* freeBlsmStk,
						Size* newBlsmId) const;
    ///
    // find a maximal vertex disjoint augmenting paths 
    // used in Scaling approximation algorithm
	template<class ItmQue, class ItmStk, class ItmIdxdQue>
	void rFindMxmlSetAugPathsCard(const Graph& graph,
    							Size* mateArr, 
								Size* ptrArr, 
								Stt* sttArr,
    							Size* idxArr, 
								ItmStk& prcsbStk, 
								ItmQue& prcsdQue, 
								ItmIdxdQue& expsdQue,
    							ItmQue& sLastQue, 
								ItmQue& tLastQue, 
								Size* blsmArr,
								Size* rprsArr, 
								Size* blsmParent,
								std::vector<Size>* blsmChildren,
								std::vector<Size>* blsmL,
								Size* label, 
								Val* dualArr, 
								Size* blsm, 
								Val delta, 
								Size* blsmin, 
								Size* treeNumArr, 
								Size iter, 
								Val ep) const;
	///
	/// the core of the single source greedy 1/2-approx. edge weighted matching computation.
	template<class ItmQue>
	void rComputeHalfEdgWghtMatching(const Graph& graph, 
									Size* mateArr, 
									Size* card, 
									Val* edgWght) const;
	/// 
	Size logbase2(Val x) const;

	template<class ItmQue, class ItmStk>
	void rProcessBlsmEdmonds(const Graph& graph, Size* mateArr,
	    Size* ptrArr, Size* blsmArr, Size* rprsArr, Size* blsmParent,
	    std::vector<Size>* blsmChildren, std::vector<Size>* blsmL,
	    Size* label, Stt* sttArr, ItmStk& prcsbStk,
	    Size s1, Size s2, Size* blsm, Val* dualArr,
	    Size* blsmin, Size* treeNumArr, Size maxBlsm,
	    std::stack<Size>& freeBlsmStk, Size* newBlsmId) const;

	template<class ItmIdxdQue, class ItmQue>
	void rAugmentEdmonds(Size* mateArr, Size* ptrArr,
	    ItmIdxdQue& expsdQue, Size sLast, Size tLast,
	    Size* blsmArr, Size* rprsArr, Size* blsmParent,
	    std::vector<Size>* blsmChildren, std::vector<Size>* blsmL,
	    Size num, Size* label, Size* blsmin) const;

	template<class ItmIdxdQue, class ItmQue>
	void rAugmentOneSide(Size* mateArr, Size* ptrArr,
	    Size sStart,
	    Size* blsmArr, Size* rprsArr, Size* blsmParent,
	    std::vector<Size>* blsmChildren, std::vector<Size>* blsmL,
	    Size num, Size* label, Size* blsmin) const;

	/// Expand a T-blossom in-tree (event-driven D4).
	/// Re-labels blossom cycle children alternately T/S from entry child to base.
	template<class ItmQue, class ItmStk>
	void rSplitBlossom(const Graph& graph,
	    Size b, Size* mateArr, Size* ptrArr,
	    Size* blsmArr, Size* rprsArr, Size* blsmParentArr,
	    std::vector<Size>* blsmChildrenArr, std::vector<Size>* blsmL,
	    Size* label, Stt* sttArr, ItmStk& prcsbStk,
	    Val* dualArr, Size* blsmin, Size* treeNumArr,
	    std::stack<Size>& freeBlsmStk,
	    Size numVtxs) const;

	///////////
}; // class MatchingEngine
  ///

	

template<class ItmQue>
void MatchingEngine::rComputeHalfEdgWghtMatching(
    const Graph& graph, Size* mateArr,Size* card, Val* edgWght) const {
  	
	Size numVtxs =graph.mNumVtxs;
	Size numEdgs =graph.mNumEdgs;
	const std::vector<Size>* vtxVecArr = (numVtxs == 0) ? NULL : &graph.mVtxVecVec[0];
    const std::vector<Val>* edgWghtVecArr = (numVtxs == 0) ? NULL : &graph.mEdgWghtVecVec[0];
  	std::vector<std::pair<std::pair<Size,Size>, Val> > sExpsdLst2;
  	ResizeVector<std::pair<std::pair<Size,Size>, Val> >(&sExpsdLst2, numEdgs);
  	std::pair<std::pair<Size,Size>,Val> * sExpsdLst2Arr=  (numEdgs == 0) ? NULL : &sExpsdLst2[0];
	if (mateArr != NULL) {
      std::fill(&mateArr[0], &mateArr[numVtxs], cNullItm);
    }
	
  	*card = 0; 
    *edgWght=0.0;
	
  	Size edgeCounter =0;
  	// list of edges and their weights
  	for (Size sFirst = 0; sFirst < numVtxs; ++sFirst) {
  	   	Size numEdgs = vtxVecArr[sFirst].size();  	   
    	const Size* vtxArr = (numEdgs == 0) ? NULL : &vtxVecArr[sFirst][0];
    	const Val* edgWghtArr = (numEdgs == 0) ? NULL : &edgWghtVecArr[sFirst][0];
   	 	for (Size i = 0; i < numEdgs; ++i) {
      		Size t = vtxArr[i];
	  		Val wgt = edgWghtArr[i];
	  		if(sFirst>t)
	  		 continue;
	  		std::pair<Size,Size> p (sFirst,t);
      		sExpsdLst2Arr[edgeCounter]=std::pair<std::pair<Size,Size>, Val>(p, wgt);
      		edgeCounter++;
		}
     
  	}
  	// sort the edges
	std::sort(sExpsdLst2.begin(),sExpsdLst2.end(),ValGreaterPairs<std::pair<std::pair<Size,Size>, Val> >());

	std::pair<Size,Size> sedge;

    for (Size current = 0; current < numEdgs; ++current) {
    	sedge =  sExpsdLst2Arr[current].first;
    	Size s = sedge.first;
    	Size t = sedge.second;
    	if(mateArr[s]!=cNullItm || mateArr[t]!=cNullItm)
    	 continue;
    	else
    	{
    		// match an edge and udate the cardinality the weight
    		mateArr[t]=s;
			mateArr[s]=t;
			++(*card);
			*edgWght += sExpsdLst2Arr[current].second;			
		}
	}
}



template<class ItmQue>
void MatchingEngine::rProcessBlsm1(const Graph& graph, Size* mateArr,
    Size* ptr1Arr, Size* ptr2Arr, Size* blsmArr, Size* rprsArr,
    Stt* sttArr, ItmQue& prcsbQue, Size s1, Size s2, Size blsm) const {
  // find vertex rprs, the representative of the blossom. search alternatively
  // along both paths toward the base of the blossom, instead of searching
  // along one path first and then along the other one, in order for the number
  // of steps in the search to be proportional to the number of vertices in the
  // blossom rather than being bounded by the number of vertices in the graph.
  // for a graph G = (V,E) the former approach would determine a complexity of
  // O(|V|) per augmenting path search, as opposed to the O(|V|^2) complexity
  // determined by the latter. the number of blossoms discovered during an
  // augmenting path search is O(|V|). the number of vertices "shrunk" into
  // blossoms during an augmenting path search is also O(|V|).
  Size rprs(cNullItm);
  Size tmpRprsArr[2];
  tmpRprsArr[0] = rprsArr[s1];
  tmpRprsArr[1] = rprsArr[s2];
  blsmArr[tmpRprsArr[0]] = blsm;
  blsmArr[tmpRprsArr[1]] = blsm;
  bool done = false;
  do {
    for (Size k = 0; k < 2; ++k) {
      Size r = tmpRprsArr[k];
      Size s = ptr1Arr[r];
      if (s != cNullItm) {
        r = rprsArr[s];
        if (blsmArr[r] == blsm) {
          rprs = r;
          done = true;
          break;
        }
        blsmArr[r] = blsm;
        tmpRprsArr[k] = r;
      }
    }
  } while (done == false);
  // make every T-vertex in the blossom an S-vertex and make vertex rprs the
  // representative of such a T-vertex that becomes an S-vertex. this time
  // separately along both paths toward the base of the blossom since it is
  // easy now to prevent the procedure from continuing beyond the base. for a
  // graph G = (V,E), for the reasons stated above, the complexity of this task
  // is O(|V|) per augmenting path search.
  tmpRprsArr[0] = rprsArr[s1];
  tmpRprsArr[1] = rprsArr[s2];
  for (Size k = 0; k < 2; ++k) {
    Size r = tmpRprsArr[k];
    while (r != rprs) {
      Size t = mateArr[r];
      ptr1Arr[t] = s1;
      ptr2Arr[t] = s2;
      rprsArr[t] = rprs;
      prcsbQue.Push(t);
      sttArr[t] = eSttBfsPrcsb;
      Size s = ptr1Arr[r];
      r = rprsArr[s];
    }
  }
  // update the representative of every S-vertex (original S-vertex or T-vertex
  // that later became S-vertex) whose representative is supposed to become
  // rprs, the representative of the new blossom. this is the most
  // expensive part of the blossom processing task because it visits every
  // vertex in the graph. since, for a graph G = (V,E), the number of blossoms
  // discovered during an augmenting path search is O(|V|), this determines a
  // complexity of O(|V|^2) per augmenting path search.
  for (Size s = 0; s < graph.mNumVtxs; ++s) {
    Size r = rprsArr[s];
    // check if rprs is already the representative of s or if r was involved
    // in the task of determining rprs.
    if ((r == rprs) || (blsmArr[r] != blsm)) {
      continue;
    }
    Size t = mateArr[r];
    // check if r is actually inside the new blossom. it could have been
    // involved in the task of determining rprs, but it could still be
    // outside the blossom (beyond the blossom's base).
    if ((t == cNullItm) || (rprsArr[t] != rprs)) {
      continue;
    }
    rprsArr[s] = rprs;
  }
}

template<class ItmQue>
void MatchingEngine::rProcessBlsm2(const Graph& graph, Size* mateArr,
    Size* ptr1Arr, Size* ptr2Arr, Size* blsmArr, Size* rprsArr, Size* linkArr,
    Size* rankArr, Stt* sttArr, ItmQue& prcsbQue, Size s1, Size s2, Size blsm)
    const {
  // find vertex rprs, the representative of the blossom. search alternatively
  // along both paths toward the base of the blossom, instead of searching
  // along one path first and then along the other one, in order for the number
  // of steps in the search to be proportional to the number of vertices in the
  // blossom rather than being bounded by the number of vertices in the graph.
  // for a graph G = (V,E) the former approach would determine a complexity of
  // O(|V|) per augmenting path search, as opposed to the O(|V|^2) complexity
  // determined by the latter. the number of blossoms discovered during an
  // augmenting path search is O(|V|). the number of vertices "shrunk" into
  // blossoms during an augmenting path search is also O(|V|).
  Size rprs(cNullItm);
  Size tmpRprsArr[2];
  tmpRprsArr[0] = rprsArr[rFind(linkArr, s1)];
  tmpRprsArr[1] = rprsArr[rFind(linkArr, s2)];
  blsmArr[tmpRprsArr[0]] = blsm;
  blsmArr[tmpRprsArr[1]] = blsm;
  bool done = false;
  do {
    for (Size k = 0; k < 2; ++k) {
      Size r = tmpRprsArr[k];
      Size s = ptr1Arr[r];
      if (s != cNullItm) {
        r = rprsArr[rFind(linkArr, s)];
        if (blsmArr[r] == blsm) {
          rprs = r;
          done = true;
          break;
        }
        blsmArr[r] = blsm;
        tmpRprsArr[k] = r;
      }
    }
  } while (done == false);
  // make every T-vertex in the blossom an S-vertex and make vertex rprs the
  // representative of such a T-vertex that becomes an S-vertex. this time
  // separately along both paths toward the base of the blossom since it is
  // easy now to prevent the procedure from continuing beyond the base. for a
  // graph G = (V,E), for the reasons stated above, the complexity of this task
  // is O(|V|) per augmenting path search.
  tmpRprsArr[0] = rprsArr[rFind(linkArr, s1)];
  tmpRprsArr[1] = rprsArr[rFind(linkArr, s2)];
  for (Size k = 0; k < 2; ++k) {
    Size r = tmpRprsArr[k];
    while (r != rprs) {
      rUnion(linkArr, rankArr, r, rprs);
      Size t = mateArr[r];
      ptr1Arr[t] = s1;
      ptr2Arr[t] = s2;
      rUnion(linkArr, rankArr, t, rprs);
      	prcsbQue.Push(t);
      	sttArr[t] = eSttBfsPrcsb;
      rprsArr[rFind(linkArr, rprs)] = rprs;
      Size s = ptr1Arr[r];
      r = rprsArr[rFind(linkArr, s)];
    }
  }
}

template<class ItmQue>
void MatchingEngine::rProcessBlsm6(const Graph& graph, Size* mateArr,
    Size* ptr1Arr, Size* ptr2Arr, Size* blsmArr, Size* rprsArr, Size* linkArr,
    Size* rankArr, Stt* sttArr, ItmQue& prcsbQue, Size s1, Size s2, Size blsm, Val* minVtxWght, Size* sLast, Size* tLast)
    const {
  // find vertex rprs, the representative of the blossom. search alternatively
  // along both paths toward the base of the blossom, instead of searching
  // along one path first and then along the other one, in order for the number
  // of steps in the search to be proportional to the number of vertices in the
  // blossom rather than being bounded by the number of vertices in the graph.
  // for a graph G = (V,E) the former approach would determine a complexity of
  // O(|V|) per augmenting path search, as opposed to the O(|V|^2) complexity
  // determined by the latter. the number of blossoms discovered during an
  // augmenting path search is O(|V|). the number of vertices "shrunk" into
  // blossoms during an augmenting path search is also O(|V|).
  // Also return the lightest vertex in the blossom which is used by the speculative MVM
  const Val* vtxWghtArr = (graph.mNumVtxs == 0) ? NULL : &graph.mVtxWghtVec[0];
  Size rprs(cNullItm);
  Size tmpRprsArr[2];
  tmpRprsArr[0] = rprsArr[rFind(linkArr, s1)];
  tmpRprsArr[1] = rprsArr[rFind(linkArr, s2)];
  blsmArr[tmpRprsArr[0]] = blsm;
  blsmArr[tmpRprsArr[1]] = blsm;
  bool done = false;
  do {
    for (Size k = 0; k < 2; ++k) {
      Size r = tmpRprsArr[k];
      Size s = ptr1Arr[r];
      if (s != cNullItm) {
        r = rprsArr[rFind(linkArr, s)];
        if (blsmArr[r] == blsm) {
          rprs = r;
          done = true;
          break;
        }
        blsmArr[r] = blsm;
        tmpRprsArr[k] = r;
      }
    }
  } while (done == false);
  // make every T-vertex in the blossom an S-vertex and make vertex rprs the
  // representative of such a T-vertex that becomes an S-vertex. this time
  // separately along both paths toward the base of the blossom since it is
  // easy now to prevent the procedure from continuing beyond the base. for a
  // graph G = (V,E), for the reasons stated above, the complexity of this task
  // is O(|V|) per augmenting path search.
  tmpRprsArr[0] = rprsArr[rFind(linkArr, s1)];
  tmpRprsArr[1] = rprsArr[rFind(linkArr, s2)];
  for (Size k = 0; k < 2; ++k) {
    Size r = tmpRprsArr[k];
    while (r != rprs) {
      rUnion(linkArr, rankArr, r, rprs);
      Size t = mateArr[r];
      ptr1Arr[t] = s1;
      ptr2Arr[t] = s2;
      rUnion(linkArr, rankArr, t, rprs);
      	prcsbQue.Push(t);
      	sttArr[t] = eSttBfsPrcsb;
      	if (*minVtxWght > vtxWghtArr[t]) {
          *sLast = r;
          *tLast = r;
          *minVtxWght = vtxWghtArr[t];
        }
      
      rprsArr[rFind(linkArr, rprs)] = rprs;
      Size s = ptr1Arr[r];
      r = rprsArr[rFind(linkArr, s)];
    }
  }
}

template<class ItmQue>
void MatchingEngine::rFindAugPathCardSglSrc1(const Graph& graph, Size* mateArr,
    Size* ptr1Arr, Size* ptr2Arr, Size* blsmArr, Size* rprsArr, Stt* sttArr,
    ItmQue& prcsbQue, ItmQue& prcsdQue, Size sFirst, Size* sLast, Size* tLast)
    const {
  *sLast = cNullItm;
  *tLast = cNullItm;
  const std::vector<Size>*
    vtxVecArr = (graph.mNumVtxs == 0) ? NULL : &graph.mVtxVecVec[0];
  Size blsm = 0;
  prcsbQue.Push(sFirst);
  sttArr[sFirst] = eSttBfsPrcsb;
  while (prcsbQue.Empty() == false) {
    Size s = prcsbQue.Front();
    prcsbQue.Pop();
    prcsdQue.Push(s);
    sttArr[s] = eSttBfsPrcsd;
    Size sNumEdgs = vtxVecArr[s].size();
    const Size* sVtxArr = (sNumEdgs == 0) ? NULL : &vtxVecArr[s][0];
    for (Size i = 0; i < sNumEdgs; ++i) {
      Size t = sVtxArr[i];
      assert(s != t);
      // if t is an outer vertex and s and t are in different blossoms then creat a blossom
      if (sttArr[t] == eSttBfsPrcsb) {
        if (rprsArr[s] != rprsArr[t]) {
          rProcessBlsm1(graph, mateArr, ptr1Arr, ptr2Arr, blsmArr, rprsArr,
                        sttArr, prcsbQue, s, t, blsm);
        }
        ++blsm;
        continue;
      } else if (sttArr[t] == eSttBfsPrcsd) {
        continue;
      }
      Size ss = mateArr[t];
      if ((ss != cNullItm) && (sttArr[ss] != eSttIdle)) {
        continue;
      }
      if (ss == cNullItm) {
      	// an unmatched vertex is found
        *sLast = s;
        *tLast = t;
        return;
      }
      if (sttArr[ss] != eSttIdle) {
        continue;
      }
      ptr1Arr[ss] = s;
      prcsbQue.Push(ss);
      sttArr[ss] = eSttBfsPrcsb;
    }
  }
}

template<class ItmQue>
void MatchingEngine::rFindAugPathCardSglSrc2(const Graph& graph, Size* mateArr,
    Size* ptr1Arr, Size* ptr2Arr, Size* blsmArr, Size* rprsArr, Size* linkArr,
    Size* rankArr, Stt* sttArr, ItmQue& prcsbQue, ItmQue& prcsdQue,
    Size sFirst, Size* sLast, Size* tLast) const {
  *sLast = cNullItm;
  *tLast = cNullItm;
  const std::vector<Size>*
    vtxVecArr = (graph.mNumVtxs == 0) ? NULL : &graph.mVtxVecVec[0];
  Size blsm = graph.mNumVtxs;
  prcsbQue.Push(sFirst);
  sttArr[sFirst] = eSttBfsPrcsb;
  while (prcsbQue.Empty() == false) {
    Size s = prcsbQue.Front();
    prcsbQue.Pop();
    prcsdQue.Push(s);
    sttArr[s] = eSttBfsPrcsd;
    Size sNumEdgs = vtxVecArr[s].size();
    const Size* sVtxArr = (sNumEdgs == 0) ? NULL : &vtxVecArr[s][0];
    for (Size i = 0; i < sNumEdgs; ++i) {
      Size t = sVtxArr[i];
      assert(s != t);
      // if t is an outer vertex and s and t are in different blossoms then creat a blossom
      if (sttArr[t] != eSttIdle) {
      	Size b1 =rprsArr[rFind(linkArr, s)] ;
      	Size b2 = rprsArr[rFind(linkArr, t)];
        if (b1!= b2) {
          rProcessBlsm2(graph, mateArr, ptr1Arr, ptr2Arr, blsmArr, rprsArr,
                        linkArr, rankArr, sttArr, prcsbQue, s, t, blsm);
        }
        ++blsm;
        continue;
      } else if (sttArr[t] == eSttBfsPrcsd) {
        continue;
      }
      Size ss = mateArr[t];
      
      if ((ss != cNullItm) && (sttArr[ss] != eSttIdle)) {
        continue;
      }
      if (ss == cNullItm) {
      	// an unmatched vertex is found
        *sLast = s;
        *tLast = t;
       
        return;
      }
      if (sttArr[ss] != eSttIdle) {
        continue;
      }
      ptr1Arr[ss] = s;
      prcsbQue.Push(ss);
      sttArr[ss] = eSttBfsPrcsb;
    }
  }
}
//////////////////////////////////////////////

template<class ItmQue>
void MatchingEngine::rFindAugPathCardSglSrc3(const Graph& graph, Size* mateArr,
    Size* ptr1Arr, Size* ptr2Arr, Size* blsmArr, Size* rprsArr, Size* linkArr,
    Size* rankArr, Stt* sttArr, ItmQue& prcsbQue, ItmQue& prcsdQue,
    Size sFirst, Size* sLast, Size* tLast, Val w2)
    const {
  *sLast = cNullItm;
  *tLast = cNullItm;
  const std::vector<Size>* vtxVecArr = (graph.mNumVtxs == 0) ? NULL : &graph.mVtxVecVec[0];
  const Val* vtxWghtArr = (graph.mNumVtxs == 0) ? NULL : &graph.mVtxWghtVec[0];
  Size blsm = graph.mNumVtxs;
  Val maxVtxWght = cNegInfVal;
  prcsbQue.Push(sFirst);
  sttArr[sFirst] = eSttBfsPrcsb;

  while (prcsbQue.Empty() == false) {
    Size s = prcsbQue.Front();
    prcsbQue.Pop();
    prcsdQue.Push(s);
    sttArr[s] = eSttBfsPrcsd;
    Size sNumEdgs = vtxVecArr[s].size();
    const Size* sVtxArr = (sNumEdgs == 0) ? NULL : &vtxVecArr[s][0];
    for (Size i = 0; i < sNumEdgs; ++i) {
      Size t = sVtxArr[i];
      assert(s != t);
      // if t is an outer vertex and s and t are in different blossoms then creat a blossom
      if (sttArr[t] == eSttBfsPrcsd || sttArr[t] == eSttBfsPrcsb) {
        if (rprsArr[rFind(linkArr, s)] != rprsArr[rFind(linkArr, t)]) {
          rProcessBlsm2(graph, mateArr, ptr1Arr, ptr2Arr, blsmArr, rprsArr,
                        linkArr, rankArr, sttArr, prcsbQue, s, t, blsm);                     
                        ++blsm;               
        }

        continue;
      } else if (sttArr[t] == eSttLast) {
        continue;
      }
      Size ss = mateArr[t];
      if ((ss != cNullItm) && (sttArr[ss] != eSttIdle)) {
        continue;
      }
      if (ss == cNullItm) {
      if (maxVtxWght < vtxWghtArr[t]) {
      	  // a heavier unmatched vertex is found
          *sLast = s;
          *tLast = t;
          maxVtxWght = vtxWghtArr[t];
          if(fabs(w2 - maxVtxWght) <= 0.00001)
           return;
        }
        continue;
      }
      if (sttArr[ss] != eSttIdle) {
        continue;
      }
      
      ptr1Arr[ss] = s;
      prcsbQue.Push(ss);
      sttArr[ss] = eSttBfsPrcsb;
    }
  }
}
      
template<class ItmQue>
void MatchingEngine::rFindIncPathCardSglSrc3(const Graph& graph, Size* mateArr,
    Size* ptr1Arr, Size* ptr2Arr, Size* blsmArr, Size* rprsArr, Size* linkArr,
    Size* rankArr, Stt* sttArr, ItmQue& prcsbQue, ItmQue& prcsdQue,
    Size sFirst, Size* sLast, Size* tLast)
    const {
  *sLast = cNullItm;
  *tLast = cNullItm;
  const std::vector<Size>*
    vtxVecArr = (graph.mNumVtxs == 0) ? NULL : &graph.mVtxVecVec[0];
  const Val* vtxWghtArr = (graph.mNumVtxs == 0) ? NULL : &graph.mVtxWghtVec[0];
  Size blsm = graph.mNumVtxs;
  Val minVtxWght= vtxWghtArr[sFirst];
  prcsbQue.Push(sFirst);
  sttArr[sFirst] = eSttBfsPrcsb;

  while (prcsbQue.Empty() == false) {
    Size s = prcsbQue.Front();
    prcsbQue.Pop();
    prcsdQue.Push(s);
    sttArr[s] = eSttBfsPrcsd;
    Size sNumEdgs = vtxVecArr[s].size();
    const Size* sVtxArr = (sNumEdgs == 0) ? NULL : &vtxVecArr[s][0];
	
	
    for (Size i = 0; i < sNumEdgs; ++i) {
      Size t = sVtxArr[i];
      assert(s != t);
      // if t is an outer vertex and s and t are in different blossoms then creat a blossom
      if (sttArr[t] == eSttBfsPrcsd || sttArr[t] == eSttBfsPrcsb) {
        if (rprsArr[rFind(linkArr, s)] != rprsArr[rFind(linkArr, t)]) {
          rProcessBlsm6<ItmQue>(graph, mateArr, ptr1Arr, ptr2Arr, blsmArr, rprsArr,
                        linkArr, rankArr, sttArr, prcsbQue, s, t, blsm, &minVtxWght,sLast,tLast);
                        ++blsm;
        }
        
        continue;
      } else if (sttArr[t] == eSttBfsPrcsd) {
        continue;
      }
      Size ss = mateArr[t];
      if ((ss != cNullItm) && (sttArr[ss] != eSttIdle)) {
        continue;
      }
      if (minVtxWght > vtxWghtArr[ss]) {
      	// a lighter matched vertex is found
          *sLast = s;
          *tLast = t;
          minVtxWght = vtxWghtArr[ss];
          ptr1Arr[ss] = s;
          prcsbQue.Push(ss);
          sttArr[ss] = eSttBfsPrcsb;
          continue;
        }
        
      if (sttArr[ss] != eSttIdle) {
        continue;
      }
      ptr1Arr[ss] = s;
      prcsbQue.Push(ss);
      sttArr[ss] = eSttBfsPrcsb;
    }
  }
}

/////////////// Scaling approximation Algorithm functions////
template<class ItmQue>
Size MatchingEngine::rFindblsm(Size* p ,Size* blsm, Size* rprs, Size v) const 
{
   Size vblsm=blsm[v];
	if(rprs[vblsm] != cNullItm)
	{
		if(p[vblsm] == cNullItm)
		{
			return vblsm;
		}			
		else
		{
			// if the blossom tree was updated then find the root blossom by climbing the blossom tree
			while(p[vblsm] != cNullItm)
			{
			    vblsm = p[vblsm];    
			}
			blsm[v]=vblsm;
			return vblsm;
		}
	}
	else
	{
		// if the blossom has been dissolved then move up the tree until the root is found
		Size newBlsm = v;
		while(p[newBlsm] != cNullItm)
		{
			newBlsm = p[newBlsm];
		}
		blsm[v]=newBlsm;
		return newBlsm;
	}
}

template<class ItmQue>
void MatchingEngine::augmentBlsm(Size b, Size s,
                                 Size* mateArr,
                                 Size* blsmArr,
                                 Size* rprsArr,
                                 Size* blsmParent,
                                 std::vector<Size>* blsmChildren,
                                 std::vector<Size>* blsmL,
                                 Size num) const
{
  std::stack<Size> tempS;
  tempS.push(b);
  tempS.push(s);

  while (!tempS.empty())
  {
    s = tempS.top(); tempS.pop();
    b = tempS.top(); tempS.pop();

    Size shiftB = b - num;
    std::vector<Size>& kids  = blsmChildren[shiftB];
    std::vector<Size>& links = blsmL[shiftB];

    if (kids.empty())
    {
      rprsArr[b] = s;
      continue;
    }

    // Normalize: rotate kids so the current representative-child is at index 0.
    Size baseV = rprsArr[b];
    if (baseV != cNullItm)
    {
      Size baseChild = baseV;
      while ((baseChild != cNullItm) && (blsmParent[baseChild] != b))
      {
        baseChild = blsmParent[baseChild];
      }
      if (baseChild != cNullItm)
      {
        Size baseIdx = 0;
        for (; baseIdx < kids.size(); ++baseIdx)
        {
          if (kids[baseIdx] == baseChild) break;
        }
        if ((baseIdx < kids.size()) && (baseIdx != 0))
        {
          std::rotate(kids.begin(), kids.begin() + baseIdx, kids.end());
          if (!links.empty())
          {
            Size pivotL = (2 * baseIdx) % links.size();
            std::rotate(links.begin(), links.begin() + pivotL, links.end());
          }
        }
      }
    }

    // Find direct child t of blossom b containing vertex s.
    Size t = s;
    while ((t != cNullItm) && (blsmParent[t] != b))
    {
      t = blsmParent[t];
    }
    if (t == cNullItm)
    {
      // Defensive: avoid crashes if structure is corrupted.
      rprsArr[b] = s;
      continue;
    }

    // If t is a nested blossom, recurse into it first.
    if (t >= num)
    {
      tempS.push(t);
      tempS.push(s);
    }

    // Locate t in kids.
    Size idxEntry = 0;
    for (; idxEntry < kids.size(); ++idxEntry)
    {
      if (kids[idxEntry] == t) break;
    }
    if (idxEntry >= kids.size())
    {
      rprsArr[b] = s;
      continue;
    }

    const long k = (long)kids.size();
    const long L = (long)links.size();

    // Direction rule preserved from the original implementation:
    // even index -> walk backwards, odd index -> walk forwards.
    int inc = ((idxEntry % 2) == 0) ? -1 : 1;

    long idx  = (long)idxEntry;

    // Links format from rProcessBlsm3:
    //   Edge e connects child[e] to child[(e+1)%k]:
    //     links[2*e]   = vertex in child[e]
    //     links[2*e+1] = vertex in child[(e+1)%k]

    while (idx != 0)
    {
      // step 1: advance to child1
      idx += inc;
      idx %= k;
      if (idx < 0) idx += k;
      long idx1 = idx;

      // step 2: advance to child2
      idx += inc;
      idx %= k;
      if (idx < 0) idx += k;
      long idx2 = idx;

      // child1 and child2 are adjacent. The edge connecting them:
      //   Forward (inc=+1): child1 at idx1, child2 at idx1+1. Edge = idx1.
      //     child1 = child[idx1] → links[2*idx1]    (lower side)
      //     child2 = child[idx2] → links[2*idx1+1]  (higher side)
      //   Backward (inc=-1): child1 at idx1, child2 at idx1-1. Edge = idx2.
      //     child1 = child[idx1] → links[2*idx2+1]  (higher side)  
      //     child2 = child[idx2] → links[2*idx2]    (lower side)
      long edgeNum = (inc == 1) ? idx1 : idx2;
      if (edgeNum < 0) edgeNum += k;
      edgeNum %= k;

      Size child1 = kids[(Size)idx1];
      Size a1 = child1;
      if (child1 >= num && L > 0)
      {
        // child1 is on the higher side if backward, lower if forward
        long pick = (inc == 1) ? (2 * edgeNum) : (2 * edgeNum + 1);
        pick %= L;
        if (pick < 0) pick += L;
        a1 = links[(Size)pick];
        tempS.push(child1);
        tempS.push(a1);
      }

      Size child2 = kids[(Size)idx2];
      Size a2 = child2;
      if (child2 >= num && L > 0)
      {
        // child2 is on the higher side if forward, lower if backward
        long pick = (inc == 1) ? (2 * edgeNum + 1) : (2 * edgeNum);
        pick %= L;
        if (pick < 0) pick += L;
        a2 = links[(Size)pick];
        tempS.push(child2);
        tempS.push(a2);
      }

      mateArr[a1] = a2;
      mateArr[a2] = a1;
    }

    // Rebase blossom so that it begins at entry t.
    if (idxEntry != 0)
    {
      std::rotate(kids.begin(), kids.begin() + idxEntry, kids.end());
      if (!links.empty())
      {
        Size pivotL = (2 * idxEntry) % links.size();
        std::rotate(links.begin(), links.begin() + pivotL, links.end());
      }
    }

    rprsArr[b] = s;
  }
}




template<class ItmIdxdQue,class ItmQue>
void MatchingEngine::rAugment(Size* mateArr,Size* ptrArr,
    ItmIdxdQue& expsdQue, Size sLast, Size tLast, Size* blsmArr, 
	Size* rprsArr, Size* blsmParent, std::vector<Size>* blsmChildren,std::vector<Size>* blsmL,Size num, Size* label, Size* blsmin) const {
  
  Size s = sLast;
  Size t = tLast;
  // Don't erase from queue - not used properly in Edmonds
  // expsdQue.Erase(t);
  while (s != cNullItm) {
    Size blsmt=t;
    if(blsmParent[t]!=cNullItm)
   		blsmt = rFindblsm<ItmQue>(blsmParent,blsmArr, rprsArr, t);
    label[blsmt]=0;
   	if(blsmt >= num && blsmin[blsmt-num]!= cNullItm)
    {	
     	Size tempT =blsmin[blsmt-num];
     	//augment the blossom
		augmentBlsm<ItmQue>(blsmArr[t],tempT,mateArr,blsmArr, rprsArr, blsmParent, blsmChildren,blsmL,num);
		t= tempT;
    }
    
    Size blsms=s;
    if(blsmParent[s]!=cNullItm)
	 blsms =rFindblsm<ItmQue>(blsmParent,blsmArr, rprsArr, s);
    if(blsms >= num && rprsArr[blsms]!=s)
    {	
		//if s is in a blossom 
    	label[blsms]=0;
		Size base = rprsArr[blsms];
		Size mateBase = mateArr[base];
		//augment the blossom
    	augmentBlsm<ItmQue>(blsms,s,mateArr,blsmArr, rprsArr, blsmParent, blsmChildren,blsmL,num);
    	
    	label[s]=0;
		label[t]=0;
		mateArr[s] = t;
    	mateArr[t] = s;
    	s=base;
    	if(ptrArr[s] == cNullItm || mateBase ==cNullItm )
    	{
    		// expsdQue.Erase(s);
		}
    	s = ptrArr[s];
    	t=mateBase;	
	}
	else
	{
		//if s is not in blossom
		Size tt = mateArr[s];
		label[s]=0;
		label[t]=0;
		mateArr[s] = t;
    	mateArr[t] = s;
    	if(ptrArr[s] == cNullItm)
    	{
    		// expsdQue.Erase(s);
		}
    	s = ptrArr[s];
    	t=tt;
	}
  }
}

template<class ItmQue>
void MatchingEngine::rDissolveBlsm3(const Graph& graph, Size* blsmArr,Size* blsmParent , std::vector<Size>* blsmChildren, Size b) const {

    for(Size i=0; i < blsmChildren[b-graph.mNumVtxs].size();i++)
    {
        Size child = blsmChildren[b-graph.mNumVtxs][i];
        blsmParent[child] = cNullItm;
        if( child < graph.mNumVtxs)
		{
			blsmArr[child]= child;	
		}
	}
}

inline bool MatchingEngine::rAcquireBlsmId(Size* blsm, Size numVtxs, Size maxBlsm,
    std::vector<Size>* blsmChildren, std::vector<Size>* blsmL,
    Size* blsmin, std::stack<Size>* freeBlsmStk, Size* b) const
{
  // Reuse only ids that were dead before the current stage started.
  // The stage-start cache reset ensures no vertex blsmArr still points to them.
  if (freeBlsmStk != NULL && freeBlsmStk->empty() == false)
  {
    *b = freeBlsmStk->top();
    freeBlsmStk->pop();
  }
  else
  {
    *b = *blsm;
    if (*b >= numVtxs + maxBlsm)
    {
      return false;
    }
    ++(*blsm);
  }

  assert(*b >= numVtxs);
  assert(*b < numVtxs + maxBlsm);

  Size shiftB = *b - numVtxs;
  blsmChildren[shiftB].clear();
  blsmL[shiftB].clear();
  blsmin[shiftB] = cNullItm;
  return true;
}

// rDissolveBlsmRecursive: dissolves blossom b, recursively resets blsmArr[v]=v
// for ALL vertices inside (including nested sub-blossoms), clears reusable
// blossom state, and pushes freed IDs onto freeBlsmStk.
// Total work is O(k) where k = number of vertices in the blossom.
template<class ItmQue>
void MatchingEngine::rDissolveBlsmRecursive(const Graph& graph,
    Size* blsmArr, Size* blsmParentArr,
    std::vector<Size>* blsmChildren, std::vector<Size>* blsmL,
    Size* rprsArr, Size* blsmin, Size b,
    std::stack<Size>& freeBlsmStk) const
{
  Size numVtxs = graph.mNumVtxs;
  Size shiftB = b - numVtxs;

  // Use a stack-based traversal to reset all vertices inside b
  std::stack<Size> workStk;
  workStk.push(b);

  while (!workStk.empty())
  {
    Size cur = workStk.top();
    workStk.pop();

    if (cur < numVtxs)
    {
      // Vertex: reset blsmArr
      blsmArr[cur] = cur;
    }
    else
    {
      // Blossom: push all children for processing
      Size shiftCur = cur - numVtxs;
      for (Size i = 0; i < blsmChildren[shiftCur].size(); ++i)
      {
        Size child = blsmChildren[shiftCur][i];
        blsmParentArr[child] = cNullItm;
        workStk.push(child);
      }

      // Clear the blossom's data (safe for reuse)
      blsmChildren[shiftCur].clear();
      blsmL[shiftCur].clear();
      blsmin[shiftCur] = cNullItm;

      // Only recycle non-root blossoms (root is b itself, handled below)
      if (cur != b)
      {
        rprsArr[cur] = cNullItm;
        freeBlsmStk.push(cur);
      }
    }
  }

  // Clear root blossom
  blsmChildren[shiftB].clear();
  blsmL[shiftB].clear();
  blsmin[shiftB] = cNullItm;
  rprsArr[b] = cNullItm;
  freeBlsmStk.push(b);
}

template <class ItmQue,class ItmStk>
void MatchingEngine::rProcessBlsm3(const Graph& graph, Size* mateArr,
    Size* ptr1Arr, Size* blsmArr, Size* rprsArr, Size* blsmParent, std::vector<Size>* blsmChildren,std::vector<Size>* blsmL, Size* label,
    Stt* sttArr, ItmStk& prcsbStk, Size s1, Size s2, Size* blsm, Val* dualArr,
    Size* blsmin, Size* treeNumArr, Size maxBlsm,
    std::stack<Size>* freeBlsmStk, Size* newBlsmId) const {
  
  	std::queue<Size> tempQ; 
  	Size numVtxs = graph.mNumVtxs;
  	Size bs1=s1;
  	Size bs2=s2;
  	if(blsmParent[s1]!=cNullItm)	
  	 bs1 = rFindblsm<ItmQue>(blsmParent,blsmArr, rprsArr, s1);
	if(blsmParent[s2]!=cNullItm)
	 bs2 = rFindblsm<ItmQue>(blsmParent,blsmArr, rprsArr, s2);
    Size base =rprsArr[bs2];      
    // Create blossom.
    Size b = cNullItm;
    if(rAcquireBlsmId(blsm, numVtxs, maxBlsm, blsmChildren, blsmL,
        blsmin, freeBlsmStk, &b) == false) {
      return;
    }
    if (newBlsmId != NULL)
      *newBlsmId = b;
    Size shiftB =b-numVtxs;
    rprsArr[b] = base;
    treeNumArr[b]=treeNumArr[bs2];
    blsmParent[b] = cNullItm;
    blsmParent[bs2] = b;
    blsmChildren[shiftB].push_back(bs2);
    blsmL[shiftB].push_back(s2);
    blsmL[shiftB].push_back(s1);
    // iterate over the vertices in backward order until the base is reached
    while (bs1 != b)
    {	
		Size tempS1 = s1;
    	blsmParent[bs1] = b;
    	blsmChildren[shiftB].push_back(bs1);
		Size mates1;
		if(bs1 < numVtxs)
		{
			blsmL[shiftB].push_back(bs1);
			mates1 = mateArr[bs1];
		}		
		else
		{
			Size rprsbs1 = rprsArr[bs1];
			if(tempS1 != rprsbs1)
			{
				ptr1Arr[tempS1]= ptr1Arr[base];
				s1= rprsbs1;
				tempS1 = rprsbs1;
			}
			
			blsmL[shiftB].push_back(rprsbs1);
			mates1 = mateArr[rprsbs1];
		}
		blsmL[shiftB].push_back(mates1);
		bs1= mates1;
		if(blsmParent[mates1]!=cNullItm)
		 bs1 = rFindblsm<ItmQue>(blsmParent,blsmArr, rprsArr, mates1);
		blsmParent[bs1] = b;
    	blsmChildren[shiftB].push_back(bs1);
    	
		if(bs1 < numVtxs)
		{
			blsmL[shiftB].push_back(mates1);
			if(sttArr[bs1] == eSttIdle || (label[bs1]== 2 || label[bs1]== 0)) 
			{
				prcsbStk.Push(bs1);
      			sttArr[bs1] = eSttDfsPrcsb;
      			ptr1Arr[bs1]= ptr1Arr[base];
			}	
		}		
		else
		{	
			if(blsmin[bs1-numVtxs]==cNullItm || blsmin[bs1-numVtxs]==rprsArr[bs1])
			{
				blsmL[shiftB].push_back(mates1);
				ptr1Arr[mates1]= ptr1Arr[base];
			}	
			else
			{
				ptr1Arr[blsmin[bs1-numVtxs]]= ptr1Arr[base];
				blsmL[shiftB].push_back(blsmin[bs1-numVtxs]);	
			}
	
			if(sttArr[rprsArr[bs1]] == eSttIdle)
			{	
				tempQ.push(bs1-numVtxs);					
				while (tempQ.empty() == false) 
				{     	
					Size tempB = tempQ.front();
					tempQ.pop();
					for(Size i=0; i < blsmChildren[tempB].size();i++)
					{
						Size child = blsmChildren[tempB][i];
						if( child < numVtxs)
						{
							if(sttArr[child] == eSttIdle)
							{
								prcsbStk.Push(child);
								sttArr[child] = eSttDfsPrcsb;
								ptr1Arr[child] = ptr1Arr[base];
							}
						}
						else											
						{
							tempQ.push(child-numVtxs);
						} 
					}
				}
			}
			label[bs1] = 1;
		}	
		s1 = ptr1Arr[s1];
		ptr1Arr[tempS1]= ptr1Arr[base];
		blsmL[shiftB].push_back(s1);
		bs1=s1;
		if(blsmParent[s1]!=cNullItm)
		 bs1 = rFindblsm<ItmQue>(blsmParent,blsmArr, rprsArr, s1);
	}
    label[b] = 1;
    dualArr[b] = 0;
}


template<class ItmQue, class ItmStk, class ItmIdxdQue>
void MatchingEngine::rFindMxmlSetAugPathsCard(const Graph& graph,
    Size* mateArr, Size* ptrArr, Stt* sttArr,
    Size* idxArr, ItmStk& prcsbStk, ItmQue& prcsdQue, ItmIdxdQue& expsdQue,
    ItmQue& sLastQue, ItmQue& tLastQue, Size* blsmArr,Size* rprsArr, Size* blsmParent,
	std::vector<Size>* blsmChildren,std::vector<Size>* blsmL, Size* label, Val* dualArr,
	 Size* blsm, Val delta, Size* blsmin, Size* treeNumArr,Size iter, Val ep) const {

  	Size numVtxs = graph.mNumVtxs; 
   	Val gamma = logbase2(1.0/ep);
  	const std::vector<Size>* vtxVecArr = (graph.mNumVtxs == 0) ? NULL : &graph.mVtxVecVec[0];
  	const Val* vtxWghtArr = (graph.mNumVtxs == 0) ? NULL : &graph.mVtxWghtVec[0];   
  	std::queue<Size> tempQ;
  	Size sFirst = expsdQue.First();
  	while (sFirst != cNullItm) 
  	{
  		if(sttArr[sFirst] != eSttIdle)
  		{
  			sFirst = expsdQue.Next(sFirst);
  			continue;
	  	}
		prcsbStk.Push(sFirst);
    	sttArr[sFirst] = eSttDfsPrcsb;
    	label[sFirst]=1;
        treeNumArr[sFirst]=sFirst;

    	Size blsmf= sFirst;
		if(blsmParent[sFirst]!=cNullItm)
    	 blsmf = rFindblsm<ItmQue>(blsmParent,blsmArr, rprsArr, sFirst);
    	// if sFirst is in a blosso then add all vertices within the blossom to the search stack
		if(blsmf  >= numVtxs)
		{
			label[blsmf]=1;
			treeNumArr[blsmf]=sFirst; 
			tempQ.push(blsmf-numVtxs);					
			while (tempQ.empty() == false) 
			{     	
				Size tempB = tempQ.front();
				tempQ.pop();
				for(Size i=0; i < blsmChildren[tempB].size();i++)
				{
					Size child = blsmChildren[tempB][i];
					if( child < numVtxs)
					{
						if(sttArr[child] == eSttIdle)
						{
							prcsbStk.Push(child);
							sttArr[child] = eSttDfsPrcsb;
						}
					}
					else
					{
						tempQ.push(child-numVtxs);
					} 
				}
			}
		}
    	while (prcsbStk.Empty() == false) 
		{
			
      		Size s = prcsbStk.Top();
      		Size sNumEdgs = vtxVecArr[s].size();
      		const Size* vtxArr = (sNumEdgs == 0) ? NULL : &vtxVecArr[s][0];
      		Size cases = 0;
      		while (idxArr[s] < sNumEdgs) 
			{
        		Size t = vtxArr[idxArr[s]];
        		Val wght = vtxWghtArr[s]+vtxWghtArr[t];
        		Size ss;
        		Size base;
        		/////////
        		if(t == mateArr[s])
        		{
        			++(idxArr[s]);
        			continue;
				}
				
				Size blsmt=t;
				Size blsms =s;
				if(blsmParent[t]!=cNullItm)				
				 blsmt = rFindblsm<ItmQue>(blsmParent,blsmArr, rprsArr, t);
				if(blsmParent[s]!=cNullItm)
				 blsms = rFindblsm<ItmQue>(blsmParent,blsmArr, rprsArr, s);
				
				if(label[blsmt]==1 &&  treeNumArr[blsms] == treeNumArr[blsmt])
				{				
			  		if(fabs((dualArr[s]+dualArr[t] )- ((int)(wght/delta)*delta - delta)) < EPS && blsms != blsmt && sttArr[t] == eSttDfsPrcsd && logbase2(wght)>=iter-gamma)
			  		{
			  			//create a blossom
						rProcessBlsm3<ItmQue, ItmStk>(graph, mateArr,ptrArr,blsmArr, rprsArr, blsmParent,blsmChildren,blsmL,label,sttArr, prcsbStk,t, s, blsm, dualArr,blsmin,treeNumArr, numVtxs, NULL, NULL);
						break;
					}
			  		++(idxArr[s]);
			  		continue;	
				}
				
        		if(blsmt >=numVtxs)
        		{
        			base = rprsArr[blsmt];
					ss = mateArr[base];
					cases=1;
				}
				else
				{
					ss = mateArr[t];
					if(ss !=cNullItm)
					 cases= 2;
					else
					 cases= 3;
				}
        		if (((ss == cNullItm) && (sttArr[t] != eSttIdle)) || label[blsmt]==2 || (treeNumArr[blsmt] !=cNullItm && treeNumArr[blsms] != treeNumArr[blsmt] )) {
	
          				++(idxArr[s]);
        		} 
				else
				{
					// check if the edge (s,t) is an eligible edge
					if(fabs((dualArr[s]+dualArr[t] )- ((int)(wght/delta)*delta - delta)) <= EPS && logbase2(wght)>=iter-gamma) 
					{
					
						if(ss == cNullItm)
						 break;
						Size tempv = t;
						if(cases == 1 )
						{
							tempv = base;
							if(ss == cNullItm)
							{
								break;
							}
						}
						double result = dualArr[ss]+dualArr[tempv] - (floor((vtxWghtArr[ss]+vtxWghtArr[tempv])/delta)*delta );
						bool eligible = false;
						if(result >= 0)
						{    // check if (t,matet) is an eligble edge                                                                   
							if(fabs(result - floorf(result)) <= EPS  &&  fabs(((round(result/delta)*delta)- result)) <= EPS && logbase2(vtxWghtArr[ss]+vtxWghtArr[tempv])>=iter-gamma)
							{
								eligible = true;
							}
						}
						if(eligible)
						{
							break;
						}			
						else
						{
							label[t]=2;
							++(idxArr[s]);
						}	
					}
					else
					 ++(idxArr[s]);	
        		}
    		}
    		if(idxArr[s] < sNumEdgs) 
			{
        		Size t = vtxArr[idxArr[s]];
        		Size blsmt=t;
				Size blsms =s;
				if(blsmParent[t]!=cNullItm)	
        		 blsmt = rFindblsm<ItmQue>(blsmParent,blsmArr, rprsArr, t);
        		if(blsmParent[s]!=cNullItm)
        		 blsms = rFindblsm<ItmQue>(blsmParent,blsmArr, rprsArr, s);
				if(blsms == blsmt)
				{
					++(idxArr[s]);
					continue;
				}
				Size ss = mateArr[t];
				
				if(ss != cNullItm)
				{
					Size blsmss=ss;
					if(blsmParent[ss]!=cNullItm)
					 blsmss = rFindblsm<ItmQue>(blsmParent,blsmArr, rprsArr, ss);
					if(label[blsmt]==0)
					{
						//if t has not been considered					
						if(blsmt >= numVtxs) 
						{
							// if t is on a blossom
							Size base = rprsArr[blsmt];		
							blsmin[blsmt-numVtxs]=t;	
							label[blsmt]=2;
							treeNumArr[blsmt]=treeNumArr[blsms];
							ss = mateArr[base];
							if(ss == cNullItm)
							{
								// if the base of a blossom is unmatched then an augmenting path is found
								if(sttArr[base] == eSttIdle)
								{
									label[blsmt]=1;
									sLastQue.Push(s);
					          		tLastQue.Push(base);
					          		sttArr[base] = eSttLast;
					          		while (prcsbStk.Empty() == false) 
									{
					            		Size sss = prcsbStk.Top();			
					            		prcsbStk.Pop();
					            		prcsdQue.Push(sss);
					            		sttArr[sss] = eSttDfsPrcsd;	
					          		}
					          		break;
								}
								else
								{
									++(idxArr[s]);
						        	continue;
								}
							}
							// other wise move to the base of a blossom
							ptrArr[ss] = s;	
							blsmss=ss;
							if(blsmParent[ss]!=cNullItm)
							 blsmss = rFindblsm<ItmQue>(blsmParent,blsmArr, rprsArr, ss);
							treeNumArr[blsmss]=treeNumArr[blsms];
							if(blsmss >= numVtxs)
							{
								// if mste of t is in a blosso then add all vertices within the blossom to the search stack
								prcsbStk.Push(ss);
								sttArr[ss] = eSttDfsPrcsb;
								label[blsmss]=1;
								tempQ.push(blsmss-numVtxs);
								while (tempQ.empty() == false) 
								{     	
									Size tempB = tempQ.front();
									tempQ.pop();
									for(Size i=0; i < blsmChildren[tempB].size();i++)
									{
									    Size child = blsmChildren[tempB][i];
									    if( child < numVtxs)
										{
									      	if(sttArr[child] == eSttIdle)
									      	{
									      		prcsbStk.Push(child);
									      		sttArr[child] = eSttDfsPrcsb;
									      		ptrArr[child] = s;
											}
										}
										else
										{
											tempQ.push(child-numVtxs);
										} 
									}
								}
							}
							else
							{
								label[ss]=1;
								prcsbStk.Push(ss);
						        sttArr[ss] = eSttDfsPrcsb;
							}
							++(idxArr[s]);	
						}
						else if(blsmss >= numVtxs)							
						{
							// if mate of t is in a blosso then add all vertices within the blossom to the search stack
							label[blsmt]=2;
							treeNumArr[blsmt]=treeNumArr[blsms];
							ptrArr[ss] = s;
							prcsbStk.Push(ss);
							sttArr[ss] = eSttDfsPrcsb;
							label[blsmss]=1;
							treeNumArr[blsmss]=treeNumArr[blsms];
							tempQ.push(blsmss-numVtxs);
							while (tempQ.empty() == false) 
							{     	
								Size tempB = tempQ.front();
								tempQ.pop();
								for(Size i=0; i < blsmChildren[tempB].size();i++)
								{
									Size child = blsmChildren[tempB][i];
									if( child < numVtxs)
									{					
									    if(sttArr[child] == eSttIdle)
									    {
									      	ptrArr[child] = s;
									      	prcsbStk.Push(child);
									    	sttArr[child] = eSttDfsPrcsb;
										}			
									}
									else
									{
										tempQ.push(child-numVtxs);
									} 
								}
							}
							++(idxArr[s]);
						}
						else
						{
							// add the mate of t to the search stack
							label[blsmt]=2;	
							treeNumArr[blsmt]=treeNumArr[blsms];							
							label[ss]=1;
							treeNumArr[ss]=treeNumArr[blsms];
							ptrArr[ss] = s;
						    prcsbStk.Push(ss);
						    sttArr[ss] = eSttDfsPrcsb;
						    ++(idxArr[s]);
						}	
					}
					else
					 ++(idxArr[s]);
				}
				else 
				{
					// an ugmenting path is found
	          		sLastQue.Push(s);
	          		tLastQue.Push(t);
	          		sttArr[t] = eSttLast;
	          		label[blsmt]=1;
	          		if(blsmt>=numVtxs)	          				
	          		 blsmin[blsmt-numVtxs]=cNullItm;
	          		while (prcsbStk.Empty() == false) 
					{
	            		Size sss = prcsbStk.Top();
	            		prcsbStk.Pop();
	            		prcsdQue.Push(sss);
	            		sttArr[sss] = eSttDfsPrcsd;	
	          		}
	          		break;
	        	}		       
      		} 
	  		else 
	  		{
		        prcsbStk.Pop();
		        prcsdQue.Push(s);
		        sttArr[s] = eSttDfsPrcsd;
      		}
    	}
    	sFirst = expsdQue.Next(sFirst);
  	}
}

template<class ItmQue, class ItmStk, class ItmIdxdQue>
void MatchingEngine::rScaleOneMinusEpsVtxWght(const Graph& graph, Size* mateArr, Size* card, Val* vtxWght, Val ep) const {

  Size numVtxs = graph.mNumVtxs; 
  if (mInlz == false) {
    if (mateArr != NULL) {
      std::fill(&mateArr[0], &mateArr[numVtxs], cNullItm);
    }
    *card = 0;
  } else {
    Size iniCard(0);
    rInlzGrdyForCard(graph, mateArr, &iniCard);

    *card = iniCard;
  }
  //to store a pointer from an outer vertex to its parent in the search tree
  std::vector<Size> ptrVec;
  ResizeVector<Size>(&ptrVec,numVtxs);
  Size* ptrArr = (numVtxs == 0) ? NULL : &ptrVec[0];
  if (ptrArr != NULL) {
    std::fill(&ptrArr[0], &ptrArr[numVtxs], cNullItm);
  }
  //to store a blossom id for a vertex
  std::vector<Size> blsmVec;
  ResizeVector<Size>(&blsmVec, numVtxs);
  Size* blsmArr = (numVtxs == 0) ? NULL : &blsmVec[0];
   for (Size v = 0; v < graph.mNumVtxs; ++v) {
    blsmArr[v] = v;
  }
  //to store a base vertex for a vertex or a blossom
  std::vector<Size> rprsVec;
  ResizeVector<Size>(&rprsVec, numVtxs*2);
  Size* rprsArr = (numVtxs == 0) ? NULL : &rprsVec[0];
  for (Size v = 0; v < numVtxs; ++v) {
    rprsArr[v] = v;
  }
  std::fill(&rprsArr[numVtxs], &rprsArr[numVtxs*2], cNullItm);
  for (Size v = numVtxs; v < numVtxs*2; ++v) {
    rprsArr[v] = cNullItm;
  }
  //2
  //to store a label for a vertex (0 for unlabelled vertex 1 for inner vertex 2 for an outer vertex)
  std::vector<Size> labelVec;
  ResizeVector<Size>(&labelVec, numVtxs*2);
  Size* label = (numVtxs == 0) ? NULL : &labelVec[0];
  if (label != NULL) {
    std::fill(&label[0], &label[numVtxs], 0);
  }
  
  //to store a blossom parent in a blossom tree
  std::vector<Size> blsmParentVec;
  ResizeVector<Size>(&blsmParentVec, numVtxs*2);
  Size* blsmParentArr = (numVtxs == 0) ? NULL : &blsmParentVec[0];
  if (blsmParentArr != NULL) {
    std::fill(&blsmParentArr[0], &blsmParentArr[numVtxs*2], cNullItm);
   
  }
	
	//to store a blossom children in a blossom tree
  std::vector<std::vector<Size> > blsmChildrenVec;
  ResizeVector<std::vector<Size> >(&blsmChildrenVec, numVtxs);
  std::vector<Size>* blsmChildrenArr = (numVtxs == 0) ? NULL : &blsmChildrenVec[0];
	
	//to store a blossom edges
  std::vector<std::vector<Size> > blsmLVec;
  ResizeVector<std::vector<Size> >(&blsmLVec, numVtxs);
  std::vector<Size>* blsmL = (numVtxs == 0) ? NULL : &blsmLVec[0];
  
  //to store which vertex an augmenting path is using in a blossom
  std::vector<Size> blsminVec;
  ResizeVector<Size>(&blsminVec, numVtxs);
  Size* blsmin = (numVtxs == 0) ? NULL : &blsminVec[0];
  if (blsmin != NULL) {
    std::fill(&blsmin[0], &blsmin[numVtxs], cNullItm);
  }
  //2
  //to store a tree id for each vertex or a blossom
  std::vector<Size> treeNumVec;
  ResizeVector<Size>(&treeNumVec, numVtxs*2);
  Size* treeNumArr = (numVtxs == 0) ? NULL : &treeNumVec[0];
    if (treeNumArr != NULL) {
    std::fill(&treeNumArr[0], &treeNumArr[numVtxs*2], cNullItm);
  }
  
  const std::vector<Size>*
    vtxVecArr = (numVtxs == 0) ? NULL : &graph.mVtxVecVec[0];
   //const std::vector<Val>* edgWghtVecArr = (graph.mNumVtxs == 0) ? NULL : &graph.mEdgWghtVecVec[0]; 
   const Val* vtxWghtArr = (numVtxs == 0) ? NULL : &graph.mVtxWghtVec[0];    

  Val maxweight = 0.0;
  // find the max edge weight
  for (Size v = 0; v < numVtxs; ++v) {
    Size sNumEdgs = vtxVecArr[v].size();
    const Size* vtxArr = (sNumEdgs == 0) ? NULL : &vtxVecArr[v][0];
	 for (Size i = 0; i < sNumEdgs; ++i) {
	      Size t = vtxArr[i];
		      
	      	if(vtxWghtArr[v]+vtxWghtArr[t] > maxweight)
	      	 maxweight = vtxWghtArr[v]+vtxWghtArr[t];
  }
}
  
 
  Val delta = maxweight*ep;
  //2
  //to store dual values for each vertex of a blossom
  std::vector<Val> dualVec;
  ResizeVector<Val>(&dualVec, numVtxs*2);
  Val* dualArr = (numVtxs== 0) ? NULL : &dualVec[0];
  std::fill(&dualArr[0], &dualArr[numVtxs], maxweight/2.0 - delta/2.0);
  std::fill(&dualArr[numVtxs], &dualArr[numVtxs*2], 0);
  
  //to store a status for each vertex or a blossom
  std::vector<Stt> sttVec;
  ResizeVector<Stt>(&sttVec, numVtxs);
  Stt* sttArr = (numVtxs == 0) ? NULL : &sttVec[0];
  if (sttArr != NULL) {
    std::fill(&sttArr[0], &sttArr[numVtxs], eSttIdle);
}
  std::vector<Size> idxVec;
  ResizeVector<Size>(&idxVec, numVtxs);
  Size* idxArr = (numVtxs == 0) ? NULL : &idxVec[0];
  if (idxArr != NULL) {
    std::fill(&idxArr[0], &idxArr[numVtxs], 0);
  }

  	ItmQue sLastQue(numVtxs);
  	ItmQue tLastQue(numVtxs);
  
  	ItmStk prcsbStk(numVtxs);
	ItmQue prcsdQue(numVtxs);
  	ItmIdxdQue expsdQue(numVtxs);
	
	Size blsm = numVtxs;

	for (Size s = 0; s < numVtxs; ++s) 
	{
      	expsdQue.Push(s);
  	}
	
  	for(Size i = 0; i <= logbase2(maxweight) ; i++)
  	{		
  		bool doneScale=true;
  		
  		do
  		{
  			doneScale=true;
  			
  			//maximal set of vertex disjoint augmenting paths
      		rFindMxmlSetAugPathsCard<ItmQue, ItmStk, ItmIdxdQue>
        	(graph, mateArr, ptrArr, sttArr, idxArr,
         	prcsbStk, prcsdQue, expsdQue, sLastQue, tLastQue,blsmArr,rprsArr,blsmParentArr,blsmChildrenArr,blsmL,label,
			  dualArr, &blsm, delta,blsmin,treeNumArr,i,ep);
			
			// augmenting the set of maximal  augmenting paths
    		while (tLastQue.Empty() == false) 
			{
      			Size sLast = sLastQue.Front();
      			sLastQue.Pop();
      			Size tLast = tLastQue.Front();
      			tLastQue.Pop();
      			sttArr[tLast] = eSttIdle;
      			rAugment<ItmIdxdQue,ItmQue>(mateArr, ptrArr, expsdQue, sLast, tLast, blsmArr, rprsArr, blsmParentArr,
				    blsmChildrenArr,blsmL,numVtxs,label,blsmin);
      			
      			++(*card);
    		}
    		// vertices dual updates
  			for (Size v = 0; v < numVtxs; ++v) 
			{
				Size blsmv=v;
				if(blsmParentArr[v]!=cNullItm)
	    		{
				blsmv = rFindblsm<ItmQue>(blsmParentArr,blsmArr, rprsArr, v);
				}
			    if( blsmv>=graph.mNumVtxs && label[blsmv]==1)
				{
					dualArr[v] = fabs(dualArr[v] - delta/2);
				}
				else if(blsmv>=graph.mNumVtxs && label[ blsmv]==2)
				{
					dualArr[v] += delta/2;
				}
  				else if(label[v]==1 || mateArr[v]==cNullItm) 
    				{
    					dualArr[v] = fabs(dualArr[v] - delta/2);
					}
    			else if(label[v]==2)
    			{
    				dualArr[v] += delta/2;
				}
				label[v]=0;
				treeNumArr[v]=cNullItm;
    				
  			}

            //blossoms dual updates
  			for(Size v = numVtxs; v < blsm; ++v) 
			{
  				if(rprsArr[v]!=cNullItm && blsmParentArr[v]==cNullItm)
  				{
  					if(label[v]==1)
		    			{
		    				dualArr[v] += delta;
						}
		    		else if(label[v]==2)
		    			{
		    				dualArr[v] = fabs(dualArr[v]-delta);
						}
	  			}
	  			label[v]=0;
	  			treeNumArr[v]=cNullItm;
  			}
  		
  			//dissolving inner blossoms if the dual becomes zero
  			for (Size v = numVtxs; v < blsm; ++v)
			{
  				if(rprsArr[v]!=cNullItm && blsmParentArr[v]==cNullItm && fabs(dualArr[v]) <= EPS)
  				{
  						rprsArr[v]=cNullItm;
		    			rDissolveBlsm3<ItmQue>(graph, blsmArr,blsmParentArr, blsmChildrenArr, v);		    			
	  			}
  			}
  
  			//check if we reach the limit of iterations within a scale
  			if(i <logbase2(maxweight) )
  			{
  				Size v = expsdQue.First();
				while (v != cNullItm) 
				{
				  	
					if(dualArr[v] > (maxweight/(pow(2,i+2)))-delta/2)
  					{
  						doneScale=false;
  						break;
					}
					v = expsdQue.Next(v);				
				}
			}
			else
			{
				Size v = expsdQue.First();
				while (v != cNullItm) 
				{
				  	
					if(dualArr[v] > EPS)
  					{
  						doneScale=false;
  						break;
					}
					v = expsdQue.Next(v);				
				}
			}
			while (prcsdQue.Empty() == false) 
			{
      			Size s = prcsdQue.Front();
      			prcsdQue.Pop();
      			ptrArr[s] = cNullItm;
     			sttArr[s] = eSttIdle;
      			idxArr[s] = 0;	
    		}
		}while(!doneScale);
  		
  		// update the duals
  		delta = delta/2.0;
  		for (Size v = 0; v < numVtxs; ++v) 
		{
    		dualArr[v] += delta;
		}
	}// main scale loop
	
	*vtxWght = 0.0;
  for (Size s = 0; s < numVtxs; ++s) {
    if (mateArr[s] != cNullItm) {
      *vtxWght += vtxWghtArr[s];
 
    }
   }
}

template<class ItmQue>
void MatchingEngine::rComputeMaxCardMatchingSglSrc1(const Graph& graph,
    Size* mateArr, Size* card) const {
  if (mInlz == false) {
    if (mateArr != NULL) {
      std::fill(&mateArr[0], &mateArr[graph.mNumVtxs], cNullItm);
    }
    *card = 0;
  } else {
    Size iniCard(0);
    //initialize with maximal cardinality
    rInlzGrdyForCard(graph, mateArr, &iniCard);
    *card = iniCard;
  }
  //to store a pointer from an outer vertex to its parent in the search tree
  std::vector<Size> ptr1Vec;
  ResizeVector<Size>(&ptr1Vec, graph.mNumVtxs);
  Size* ptr1Arr = (graph.mNumVtxs == 0) ? NULL : &ptr1Vec[0];
  if (ptr1Arr != NULL) {
    std::fill(&ptr1Arr[0], &ptr1Arr[graph.mNumVtxs], cNullItm);
  }
  //to store a pointer in from an outer vertex to its parent after its status is changed in a blossom
  std::vector<Size> ptr2Vec;
  ResizeVector<Size>(&ptr2Vec, graph.mNumVtxs);
  Size* ptr2Arr = (graph.mNumVtxs == 0) ? NULL : &ptr2Vec[0];
  if (ptr2Arr != NULL) {
    std::fill(&ptr2Arr[0], &ptr2Arr[graph.mNumVtxs], cNullItm);
  }
  
  //to store a blossom id
  std::vector<Size> blsmVec;
  ResizeVector<Size>(&blsmVec, graph.mNumVtxs);
  Size* blsmArr = (graph.mNumVtxs == 0) ? NULL : &blsmVec[0];
  if (blsmArr != NULL) {
    std::fill(&blsmArr[0], &blsmArr[graph.mNumVtxs], cNullItm);
  }
  //to store a base id of a blossom
  std::vector<Size> rprsVec;
  ResizeVector<Size>(&rprsVec, graph.mNumVtxs);
  Size* rprsArr = (graph.mNumVtxs == 0) ? NULL : &rprsVec[0];
  for (Size v = 0; v < graph.mNumVtxs; ++v) {
    rprsArr[v] = v;
  }
  //to store a status of a vertex
  std::vector<Stt> sttVec;
  ResizeVector<Stt>(&sttVec, graph.mNumVtxs);
  Stt* sttArr = (graph.mNumVtxs == 0) ? NULL : &sttVec[0];
  if (sttArr != NULL) {
    std::fill(&sttArr[0], &sttArr[graph.mNumVtxs], eSttIdle);
  }
  ItmQue prcsbQue(graph.mNumVtxs);
  ItmQue prcsdQue(graph.mNumVtxs);
  for (Size sFirst = 0; sFirst < graph.mNumVtxs; ++sFirst) {
    if (mateArr[sFirst] != cNullItm) {
      continue;
    }
    Size sLast(cNullItm);
    Size tLast(cNullItm);
    // find an augmenting path
    rFindAugPathCardSglSrc1<ItmQue>
      (graph, mateArr, ptr1Arr, ptr2Arr, blsmArr, rprsArr, sttArr, prcsbQue,
       prcsdQue, sFirst, &sLast, &tLast);
    if (tLast != cNullItm) {
    	//augment an augmenting path
      rAugment(mateArr, ptr1Arr, ptr2Arr, sLast, tLast);
      ++(*card);
    }
    //reset working variable for visited vertices
    while (prcsbQue.Empty() == false) {
      Size v = prcsbQue.Front();
      prcsbQue.Pop();
      ptr1Arr[v] = cNullItm;
      ptr2Arr[v] = cNullItm;
      blsmArr[v] = cNullItm;
      rprsArr[v] = v;
      sttArr[v] = eSttIdle;
    }
    while (prcsdQue.Empty() == false) {
      Size v = prcsdQue.Front();
      prcsdQue.Pop();
      ptr1Arr[v] = cNullItm;
      ptr2Arr[v] = cNullItm;
      blsmArr[v] = cNullItm;
      rprsArr[v] = v;
      sttArr[v] = eSttIdle;
    }
  }
}

template<class ItmQue>
void MatchingEngine::rComputeMaxCardMatchingSglSrc2(const Graph& graph,
    Size* mateArr, Size* card) const {
  if (mInlz == false) {
    if (mateArr != NULL) {
      std::fill(&mateArr[0], &mateArr[graph.mNumVtxs], cNullItm);
    }
    *card = 0;
  } else {
    Size iniCard(0);
    //initialize with maximal cardinality
    rInlzGrdyForCard(graph, mateArr, &iniCard);
    *card = iniCard;
  }

  //to store a pointer from an outer vertex to its parent in the search tree
  std::vector<Size> ptr1Vec;
  ResizeVector<Size>(&ptr1Vec, graph.mNumVtxs);
  Size* ptr1Arr = (graph.mNumVtxs == 0) ? NULL : &ptr1Vec[0];
  if (ptr1Arr != NULL) {
    std::fill(&ptr1Arr[0], &ptr1Arr[graph.mNumVtxs], cNullItm);
  }
  //to store a pointer in from an outer vertex to its parent after its status is changed in a blossom
  std::vector<Size> ptr2Vec;
  ResizeVector<Size>(&ptr2Vec, graph.mNumVtxs);
  Size* ptr2Arr = (graph.mNumVtxs == 0) ? NULL : &ptr2Vec[0];
  if (ptr2Arr != NULL) {
    std::fill(&ptr2Arr[0], &ptr2Arr[graph.mNumVtxs], cNullItm);
  }
  //to store a blossom id
  std::vector<Size> blsmVec;
  ResizeVector<Size>(&blsmVec, graph.mNumVtxs);
  Size* blsmArr = (graph.mNumVtxs == 0) ? NULL : &blsmVec[0];
  if (blsmArr != NULL) {
    std::fill(&blsmArr[0], &blsmArr[graph.mNumVtxs], cNullItm);
  }
  //to store a base id of a blossom
  std::vector<Size> rprsVec;
  ResizeVector<Size>(&rprsVec, graph.mNumVtxs);
  Size* rprsArr = (graph.mNumVtxs == 0) ? NULL : &rprsVec[0];
  for (Size v = 0; v < graph.mNumVtxs; ++v) {
    rprsArr[v] = v;
  }
  //to store links between items in union find tree
  std::vector<Size> linkVec;
  ResizeVector<Size>(&linkVec, graph.mNumVtxs);
  Size* linkArr = (graph.mNumVtxs == 0) ? NULL : &linkVec[0];
  for (Size v = 0; v < graph.mNumVtxs; ++v) {
    linkArr[v] = v;
  }
  //to store a rank of an item in union find tree
  std::vector<Size> rankVec;
  ResizeVector<Size>(&rankVec, graph.mNumVtxs);
  Size* rankArr = (graph.mNumVtxs == 0) ? NULL : &rankVec[0];
  if (rankArr != NULL) {
    std::fill(&rankArr[0], &rankArr[graph.mNumVtxs], 0);
  }
  //to store a status of a vertex
  std::vector<Stt> sttVec;
  ResizeVector<Stt>(&sttVec, graph.mNumVtxs);
  Stt* sttArr = (graph.mNumVtxs == 0) ? NULL : &sttVec[0];
  if (sttArr != NULL) {
    std::fill(&sttArr[0], &sttArr[graph.mNumVtxs], eSttIdle);
  }
  ItmQue prcsbQue(graph.mNumVtxs);
  ItmQue prcsdQue(graph.mNumVtxs);
  for (Size sFirst = 0; sFirst < graph.mNumVtxs; ++sFirst) {
    if (mateArr[sFirst] != cNullItm) {
      continue;
    }
    Size sLast(cNullItm);
    Size tLast(cNullItm);
    // find an augmenting path
    rFindAugPathCardSglSrc2<ItmQue>
      (graph, mateArr, ptr1Arr, ptr2Arr, blsmArr, rprsArr, linkArr, rankArr,
       sttArr, prcsbQue, prcsdQue, sFirst, &sLast, &tLast);
    if (tLast != cNullItm) {
    	//augment an augmenting path	
      rAugment(mateArr, ptr1Arr, ptr2Arr, sLast, tLast);
      ++(*card);
      //reset working variable for visited vertices
      while (prcsbQue.Empty() == false) {
      Size v = prcsbQue.Front();
      prcsbQue.Pop();
      ptr1Arr[v] = cNullItm;
      ptr2Arr[v] = cNullItm;
      blsmArr[v] = cNullItm;
      rprsArr[v] = v;
      linkArr[v] = v;
      rankArr[v] = 0;
      sttArr[v] = eSttIdle;
    }
    while (prcsdQue.Empty() == false) {
      Size v = prcsdQue.Front();
      prcsdQue.Pop();
      ptr1Arr[v] = cNullItm;
      ptr2Arr[v] = cNullItm;
      blsmArr[v] = cNullItm;
      rprsArr[v] = v;
      linkArr[v] = v;
      rankArr[v] = 0;
      sttArr[v] = eSttIdle;
    }
    continue;
    }
    //mark visited vertices so they are not visited again in future searches
    while (prcsbQue.Empty() == false) {
      Size v = prcsbQue.Front();
      prcsbQue.Pop();
      ptr1Arr[v] = cNullItm;
      ptr2Arr[v] = cNullItm;
      blsmArr[v] = cNullItm;
      rprsArr[v] = v;
      linkArr[v] = v;
      rankArr[v] = 0;
      sttArr[v] = eSttLast;
    }
    while (prcsdQue.Empty() == false) {
      Size v = prcsdQue.Front();
      prcsdQue.Pop();
      ptr1Arr[v] = cNullItm;
      ptr2Arr[v] = cNullItm;
      blsmArr[v] = cNullItm;
      rprsArr[v] = v;
      linkArr[v] = v;
      rankArr[v] = 0;
      sttArr[v] = eSttLast;
    }
  }	
}


template<class ItmQue>
void MatchingEngine::rComputeMaxCardMatchingSglSrc3(const Graph& graph,
    Size* mateArr, Size* card,std::pair<Size,Val>* sortedArr) const {
  if (mInlz == false) {
    if (mateArr != NULL) {
      std::fill(&mateArr[0], &mateArr[graph.mNumVtxs], cNullItm);
    }
    *card = 0;
  } else {
  	//initialize with maximal cardinality
    Size iniCard(0);
    rInlzGrdyForCard(graph, mateArr, &iniCard);
    *card = iniCard;
  }
  //Timer timer;
  //timer.Start();
  //to store a pointer from an outer vertex to its parent in the search tree
  std::vector<Size> ptr1Vec;
  ResizeVector<Size>(&ptr1Vec, graph.mNumVtxs);
  Size* ptr1Arr = (graph.mNumVtxs == 0) ? NULL : &ptr1Vec[0];
  if (ptr1Arr != NULL) {
    std::fill(&ptr1Arr[0], &ptr1Arr[graph.mNumVtxs], cNullItm);
  }
  //to store a pointer in from an outer vertex to its parent after its status is changed in a blossom
  std::vector<Size> ptr2Vec;
  ResizeVector<Size>(&ptr2Vec, graph.mNumVtxs);
  Size* ptr2Arr = (graph.mNumVtxs == 0) ? NULL : &ptr2Vec[0];
  if (ptr2Arr != NULL) {
    std::fill(&ptr2Arr[0], &ptr2Arr[graph.mNumVtxs], cNullItm);
  }
  //to store a blossom id
  std::vector<Size> blsmVec;
  ResizeVector<Size>(&blsmVec, graph.mNumVtxs);
  Size* blsmArr = (graph.mNumVtxs == 0) ? NULL : &blsmVec[0];
  if (blsmArr != NULL) {
    std::fill(&blsmArr[0], &blsmArr[graph.mNumVtxs], cNullItm);
  }
  //to store a base id of a blossom
  std::vector<Size> rprsVec;
  ResizeVector<Size>(&rprsVec, graph.mNumVtxs);
  Size* rprsArr = (graph.mNumVtxs == 0) ? NULL : &rprsVec[0];
  for (Size v = 0; v < graph.mNumVtxs; ++v) {
    rprsArr[v] = v;
  }
  //to store links between items in union find tree
  std::vector<Size> linkVec;
  ResizeVector<Size>(&linkVec, graph.mNumVtxs);
  Size* linkArr = (graph.mNumVtxs == 0) ? NULL : &linkVec[0];
  for (Size v = 0; v < graph.mNumVtxs; ++v) {
    linkArr[v] = v;
  }
  //to store a rank of an item in union find tree
  std::vector<Size> rankVec;
  ResizeVector<Size>(&rankVec, graph.mNumVtxs);
  Size* rankArr = (graph.mNumVtxs == 0) ? NULL : &rankVec[0];
  if (rankArr != NULL) {
    std::fill(&rankArr[0], &rankArr[graph.mNumVtxs], 0);
  }
  //to store a status of a vertex
  std::vector<Stt> sttVec;
  ResizeVector<Stt>(&sttVec, graph.mNumVtxs);
  Stt* sttArr = (graph.mNumVtxs == 0) ? NULL : &sttVec[0];
  if (sttArr != NULL) {
    std::fill(&sttArr[0], &sttArr[graph.mNumVtxs], eSttIdle);
  }
  ItmQue prcsbQue(graph.mNumVtxs);
  ItmQue prcsdQue(graph.mNumVtxs);
   for (Size current = 0; current <graph.mNumVtxs; ++current) {
    	Size sFirst =  sortedArr[current].first;
    if (mateArr[sFirst] != cNullItm) {
      continue;
    }
    Size sLast(cNullItm);
    Size tLast(cNullItm);
    // find an augmenting path
    rFindAugPathCardSglSrc2<ItmQue>
      (graph, mateArr, ptr1Arr, ptr2Arr, blsmArr, rprsArr, linkArr, rankArr,
       sttArr, prcsbQue, prcsdQue, sFirst, &sLast, &tLast);
    if (tLast != cNullItm) {
    	//augment an augmenting path
      rAugment(mateArr, ptr1Arr, ptr2Arr, sLast, tLast);
      ++(*card);
      //reset visited vertices working variables
      while (prcsbQue.Empty() == false) {
      Size v = prcsbQue.Front();
      prcsbQue.Pop();
      ptr1Arr[v] = cNullItm;
      ptr2Arr[v] = cNullItm;
      blsmArr[v] = cNullItm;
      rprsArr[v] = v;
      linkArr[v] = v;
      rankArr[v] = 0;
      sttArr[v] = eSttIdle;
    }
    while (prcsdQue.Empty() == false) {
      Size v = prcsdQue.Front();
      prcsdQue.Pop();
      ptr1Arr[v] = cNullItm;
      ptr2Arr[v] = cNullItm;
      blsmArr[v] = cNullItm;
      rprsArr[v] = v;
      linkArr[v] = v;
      rankArr[v] = 0;
      sttArr[v] = eSttIdle;
    }
    continue;
    }
    //mark visited vertices so they are not visited again in future searches
    while (prcsbQue.Empty() == false) {
      Size v = prcsbQue.Front();
      prcsbQue.Pop();
      ptr1Arr[v] = cNullItm;
      ptr2Arr[v] = cNullItm;
      blsmArr[v] = cNullItm;
      rprsArr[v] = v;
      linkArr[v] = v;
      rankArr[v] = 0;
      sttArr[v] = eSttLast;
    }
    while (prcsdQue.Empty() == false) {
      Size v = prcsdQue.Front();
      prcsdQue.Pop();
      ptr1Arr[v] = cNullItm;
      ptr2Arr[v] = cNullItm;
      blsmArr[v] = cNullItm;
      rprsArr[v] = v;
      linkArr[v] = v;
      rankArr[v] = 0;
      sttArr[v] = eSttLast;
    }
  }		
}

template<class ItmQue>
void MatchingEngine::rComputeAprxMaxEdgWghtMatching4(const Graph& graph,
    Size* mateArr, Size* card, Val* edgWght) const {
  if (mateArr != NULL) {
    std::fill(&mateArr[0], &mateArr[graph.mNumVtxs], cNullItm);
  }
  *card = 0;
  *edgWght = 0.0;
  const std::vector<Size>*
    vtxVecArr = (graph.mNumVtxs == 0) ? NULL : &graph.mVtxVecVec[0];
  const std::vector<Val>*
    edgWghtVecArr = (graph.mNumVtxs == 0) ? NULL : &graph.mEdgWghtVecVec[0];
  // candidates.
  std::vector<Size> cndtMateVec;
  ResizeVector<Size>(&cndtMateVec, graph.mNumVtxs);
  Size* cndtMateArr = (graph.mNumVtxs == 0) ? NULL : &cndtMateVec[0];
  if (cndtMateArr != NULL) {
    std::fill(&cndtMateArr[0], &cndtMateArr[graph.mNumVtxs], cNullItm);
  }
  // queue of matched vertices to be processed.
  ItmQue matchedQue(graph.mNumVtxs);
  // initialize the queue of matched vertices.
  for (Size v = 0; v < graph.mNumVtxs; ++v) {
    Size w = cNullItm;
    Val maxEdgWght = cNegInfVal;
    Size vNumEdgs = vtxVecArr[v].size();
    const Size* vVtxArr = (vNumEdgs == 0) ? NULL : &vtxVecArr[v][0];
    const Val* vEdgWghtArr = (vNumEdgs == 0) ? NULL : &edgWghtVecArr[v][0];
    for (Size j = 0; j < vNumEdgs; ++j) {
      Size z = vVtxArr[j];
      Val vzEdgWght = vEdgWghtArr[j];
      if (z == v) {
        continue;
      }
      if (maxEdgWght < vzEdgWght) {
        w = z;
        maxEdgWght = vzEdgWght;
      }
    }
    cndtMateArr[v] = w;
    if ((w != cNullItm) && (cndtMateArr[w] == v)) {
      mateArr[v] = w;
      mateArr[w] = v;
      matchedQue.Push(v);
      ++(*card);
      *edgWght += maxEdgWght;
    }
  }
  // only one vertex out of each matched pair is queued, use an array to
  // process both.
  Size matchedArr[2];
  // process the matched vertices from the queue.
  while (matchedQue.Empty() == false) {
    Size u = matchedQue.Front();
    matchedQue.Pop();
    matchedArr[0] = u;
    matchedArr[1] = mateArr[u];
    for (Size k = 0; k < 2; ++k) {
      u = matchedArr[k];
      Size uNumEdgs = vtxVecArr[u].size();
      const Size* uVtxArr = (uNumEdgs == 0) ? NULL : &vtxVecArr[u][0];
      for (Size i = 0; i < uNumEdgs; ++i) {
        Size v = uVtxArr[i];
        if (v == u) {
          continue;
        }
        if (mateArr[v] != cNullItm) {
          continue;
        }
        if (cndtMateArr[v] != u) {
          continue;
        }
        Size w = cNullItm;
        Val maxEdgWght = cNegInfVal;
        Size vNumEdgs = vtxVecArr[v].size();
        const Size* vVtxArr = (vNumEdgs == 0) ? NULL : &vtxVecArr[v][0];
        const Val* vEdgWghtArr = (vNumEdgs == 0) ? NULL : &edgWghtVecArr[v][0];
        for (Size j = 0; j < vNumEdgs; ++j) {
          Size z = vVtxArr[j];
          Val vzEdgWght = vEdgWghtArr[j];
          if (z == v) {
            continue;
          }
          if (mateArr[z] != cNullItm) {
            continue;
          }
          if (maxEdgWght < vzEdgWght) {
            w = z;
            maxEdgWght = vzEdgWght;
          }
        }
        cndtMateArr[v] = w;
        if ((w != cNullItm) && (cndtMateArr[w] == v)) {
          mateArr[v] = w;
          mateArr[w] = v;
          matchedQue.Push(v);
          ++(*card);
          *edgWght += maxEdgWght;
        }
      }
    }
  }
}

//max vertex weighted matching using augmenting paths that reach a heaviest vertex
template<class ItmQue, class ItmStk>
void MatchingEngine::rComputeMaxVtxWght(const Graph& graph,
    Size* mateArr, Size* card, Val* vtxWght, Size augPathLenBnd) const {

    if (mateArr != NULL) {
      std::fill(&mateArr[0], &mateArr[graph.mNumVtxs], cNullItm);
    }
    *card = 0;
	Size numVtxs = graph.mNumVtxs;
	
  //to store a pointer from an outer vertex to its parent in the search tree
  std::vector<Size> ptr1Vec;
  ResizeVector<Size>(&ptr1Vec, graph.mNumVtxs);
  Size* ptr1Arr = (graph.mNumVtxs == 0) ? NULL : &ptr1Vec[0];
  if (ptr1Arr != NULL) {
    std::fill(&ptr1Arr[0], &ptr1Arr[graph.mNumVtxs], cNullItm);
  }
  //to store a pointer in from an outer vertex to its parent after its status is changed in a blossom
  std::vector<Size> ptr2Vec;
  ResizeVector<Size>(&ptr2Vec, graph.mNumVtxs);
  Size* ptr2Arr = (graph.mNumVtxs == 0) ? NULL : &ptr2Vec[0];
  if (ptr2Arr != NULL) {
    std::fill(&ptr2Arr[0], &ptr2Arr[graph.mNumVtxs], cNullItm);
  }
  //to store a blossom id
  std::vector<Size> blsmVec;
  ResizeVector<Size>(&blsmVec, graph.mNumVtxs);
  Size* blsmArr = (graph.mNumVtxs == 0) ? NULL : &blsmVec[0];
  if (blsmArr != NULL) {
     for (Size v = 0; v < graph.mNumVtxs; ++v) {
    blsmArr[v] = v;
  }
  }
  //to store a base id of a blossom
  std::vector<Size> rprsVec;
  ResizeVector<Size>(&rprsVec, graph.mNumVtxs*2);
  Size* rprsArr = (graph.mNumVtxs == 0) ? NULL : &rprsVec[0];
  for (Size v = 0; v < graph.mNumVtxs; ++v) {
    rprsArr[v] = v;
  }
  std::fill(&rprsArr[numVtxs], &rprsArr[numVtxs*2], cNullItm);
  //to store a status of a vertex
  std::vector<Stt> sttVec;
  ResizeVector<Stt>(&sttVec, graph.mNumVtxs);
  Stt* sttArr = (graph.mNumVtxs == 0) ? NULL : &sttVec[0];
  if (sttArr != NULL) {
    std::fill(&sttArr[0], &sttArr[graph.mNumVtxs], eSttIdle);
  }
  //to store links between items in union find tree
  std::vector<Size> linkVec;
  ResizeVector<Size>(&linkVec, graph.mNumVtxs);
  Size* linkArr = (graph.mNumVtxs == 0) ? NULL : &linkVec[0];
  for (Size v = 0; v < graph.mNumVtxs; ++v) {
    linkArr[v] = v;
  }
  //to store a rank of an item in union find tree
  std::vector<Size> rankVec;
  ResizeVector<Size>(&rankVec, graph.mNumVtxs);
  Size* rankArr = (graph.mNumVtxs == 0) ? NULL : &rankVec[0];
  if (rankArr != NULL) {
    std::fill(&rankArr[0], &rankArr[graph.mNumVtxs], 0);
  }
  

  ItmQue prcsbQue(graph.mNumVtxs);
  ItmQue prcsdQue(graph.mNumVtxs);
  
  const Val* vtxWghtArr = (graph.mNumVtxs == 0) ? NULL : &graph.mVtxWghtVec[0];
  
  std::vector<std::pair<Size, Val> > sExpsdLst2;
  ResizeVector<std::pair<Size, Val> >(&sExpsdLst2, numVtxs);
  std::pair<Size,Val> * sExpsdLst2Arr=  (numVtxs == 0) ? NULL : &sExpsdLst2[0]; 

  for (Size first = 0; first < numVtxs; ++first) {
    sExpsdLst2Arr[first]=std::pair<Size, Val>(first, vtxWghtArr[first]); 
  } 
  //sort vertices in non-increasing order of their weights
  std::sort(sExpsdLst2.begin(),sExpsdLst2.end(),ValGreater2<std::pair<Size, Val> >());

   for (Size current = 0; current < numVtxs; ++current) {
    	Size sFirst =  sExpsdLst2Arr[current].first; 
		
    if(mateArr[sFirst]!=cNullItm)
    continue;

    //find the next heaviest vertex
    Val secondWght =sExpsdLst2Arr[current].second;
    
    if(current+1<numVtxs)
    {
    	Size c = current+1;
    	while(c<numVtxs)
    	{
    		if(mateArr[sExpsdLst2Arr[c].first]==cNullItm)
    		{
    			secondWght = sExpsdLst2Arr[c].second;
    			break;
			}
			c++;
		}
	}
    else
    break;
    Size sLast(cNullItm);
    Size tLast(cNullItm);
    
       //find an augmenting path to a heaviest unmatched vertex
      rFindAugPathCardSglSrc3<ItmQue>
        (graph, mateArr, ptr1Arr,  ptr2Arr,blsmArr, rprsArr,linkArr, rankArr,
       sttArr, prcsbQue, prcsdQue, sFirst, &sLast, &tLast,secondWght);
   
    if (tLast != cNullItm) {
       //augment an augmenting path 
      rAugment(mateArr, ptr1Arr, ptr2Arr, sLast, tLast);
      ++(*card);
	  
	//reset visited vertices working variables
	while (prcsbQue.Empty() == false) {
      Size v = prcsbQue.Front();
      prcsbQue.Pop();
      ptr1Arr[v] = cNullItm;
      ptr2Arr[v] = cNullItm;
      blsmArr[v] = cNullItm;
      linkArr[v] = v;
      rankArr[v] = 0;
      rprsArr[v] = v;
      sttArr[v] = eSttIdle;
    } 
    
    while (prcsdQue.Empty() == false) {
      Size v = prcsdQue.Front();
      prcsdQue.Pop();
      ptr1Arr[v] = cNullItm;
      ptr2Arr[v] = cNullItm;
      blsmArr[v] = cNullItm;
      linkArr[v] = v;
      rankArr[v] = 0;
      rprsArr[v] = v;
      sttArr[v] = eSttIdle;
    }
      continue;
    }
    //mark visited vertices so they are not visited again in future searches
    while (prcsbQue.Empty() == false) {
      Size v = prcsbQue.Front();
      prcsbQue.Pop();
      ptr1Arr[v] = cNullItm;
      ptr2Arr[v] = cNullItm;
      blsmArr[v] = cNullItm;
      linkArr[v] = v;
      rankArr[v] = 0;
      rprsArr[v] = v;
      sttArr[v] = eSttLast;
    }

    while (prcsdQue.Empty() == false) {
      Size v = prcsdQue.Front();
      prcsdQue.Pop();
      ptr1Arr[v] = cNullItm;
      ptr2Arr[v] = cNullItm;
      blsmArr[v] = cNullItm;
      linkArr[v] = v;
      rankArr[v] = 0;
      rprsArr[v] = v;
      sttArr[v] = eSttLast;
    }
  
  }
  *vtxWght = 0.0;
  for (Size s = 0; s < numVtxs; ++s) {
    if (mateArr[s] != cNullItm) {
      *vtxWght += vtxWghtArr[s];
    }
   }
}

template<class ItmQue>
void MatchingEngine::rSuitorVtxWght(const Graph& graph,
    Size* mateArr, Size* card, Val* vtxWght) const {
    if (mateArr != NULL) {
      std::fill(&mateArr[0], &mateArr[graph.mNumVtxs], cNullItm);
    }
    *card = 0;
	Size numVtxs = graph.mNumVtxs;
  
  std::vector<Val> sWs;
  ResizeVector<Val>(&sWs, graph.mNumVtxs);
  Val* ws = (graph.mNumVtxs == 0) ? NULL : &sWs[0];
  if (ws != NULL) {
    std::fill(&ws[0], &ws[graph.mNumVtxs], cZeroVal);
  }

////
  const Val* vtxWghtArr = (graph.mNumVtxs == 0) ? NULL : &graph.mVtxWghtVec[0];
  
   const std::vector<Size>*
    vtxVecArr = (graph.mNumVtxs == 0) ? NULL : &graph.mVtxVecVec[0];
  for (Size first = 0; first < numVtxs; ++first)
   {
    Size current = first;
    bool done = false;
    while(!done)
    {
    	Size partner = cNullItm;
    	Val heaviest = 0.0;
    	
    	Size sNumEdgs = vtxVecArr[current].size();
	    const Size* currentVtxArr = (sNumEdgs == 0) ? NULL : &vtxVecArr[current][0];
	    //find a heaviest neighbor with a lower or no proposal
	    for (Size i = 0; i < sNumEdgs; ++i) {
	      Size t = currentVtxArr[i];

	      if (vtxWghtArr[t] + vtxWghtArr[current] >  heaviest || ((fabs(vtxWghtArr[t] + vtxWghtArr[current] -  heaviest)< EPS) && partner > t)){
				if( vtxWghtArr[t] + vtxWghtArr[current] > ws[t])
	            	{
					partner = t;
    				heaviest = vtxWghtArr[t] + vtxWghtArr[current];
					}
					else if (mateArr[t] !=cNullItm && fabs((vtxWghtArr[t] + vtxWghtArr[current] )- ws[t]) < EPS && mateArr[t]>current)
					{
						partner = t;
    					heaviest = vtxWghtArr[t] + vtxWghtArr[current];  
					}
    			
	      }
	    }
	    done = true;
	    if(heaviest > 0.0)
	    {
	    	//change partner and annul y proposal is any
	    	Size y = mateArr[partner];
	    	mateArr[partner]=current;
	    	ws[partner]=heaviest;
	    	if(y != cNullItm)
	    	{
	    		current = y;
	    		done= false;
			}
		}
	}
}

  *vtxWght = 0.0;
  for (Size s = 0; s < numVtxs; ++s) {
    if (mateArr[s] != cNullItm) {
      *vtxWght += vtxWghtArr[s];
      ++(*card);
    }
   }
   (*card)=(*card)/2;
}

template<class ItmQue>
void MatchingEngine::rParSuitorVtxWght(const Graph& graph,
    Size* mateArr, Size* card, Val* vtxWght) const {
	
	Size numVtxss = graph.mNumVtxs;
  	std::vector<Val> sWs;
  	ResizeVector<Val>(&sWs, graph.mNumVtxs);
  	Val* ws = (graph.mNumVtxs == 0) ? NULL : &sWs[0];
  	const Val* vtxWghtArr = (graph.mNumVtxs == 0) ? NULL : &graph.mVtxWghtVec[0];
   	const std::vector<Size>*
    vtxVecArr = (graph.mNumVtxs == 0) ? NULL : &graph.mVtxVecVec[0];
    Size nlocks[numVtxss];
    Val tempvtxWght = 0.0;
	Size tempcard=0;
  	#pragma omp parallel
  	{
  		Size numVtxs = numVtxss;
  		//setting locks for each vertex in parallel
		#pragma omp for schedule(static,256) 
		for(Size i=0;i<numVtxs;i++)  
    	{
        	nlocks[i]=0;
        	ws[i]=cZeroVal;
        	mateArr[i]=cNullItm;
    	}
    	// in parallel find a heaviest partner with a lower or no proposal
  		#pragma omp for schedule(static,256) 
  		for (Size first = 0; first < numVtxs; ++first)
   		{
    		Size current = first;
    		bool done = false;
    		while(!done)
    		{
    			Size partner = cNullItm;
    			Val heaviest = 0;
    	
    			Size sNumEdgs = vtxVecArr[current].size();
	    		const Size* currentVtxArr = (sNumEdgs == 0) ? NULL : &vtxVecArr[current][0];
	    		for (Size i = 0; i < sNumEdgs; ++i) {
	      		Size t = currentVtxArr[i];
	      		//while(__sync_lock_test_and_set(&nlocks[t], 1));
	      		if (vtxWghtArr[t] + vtxWghtArr[current] >  heaviest || ((fabs(vtxWghtArr[t] + vtxWghtArr[current] -  heaviest)< EPS) && partner > t)){
					if( vtxWghtArr[t] + vtxWghtArr[current] > ws[t])
	            	{
						partner = t;
    					heaviest = vtxWghtArr[t] + vtxWghtArr[current];
					}
					else if (mateArr[t] !=cNullItm && fabs((vtxWghtArr[t] + vtxWghtArr[current] )- ws[t]) < EPS && mateArr[t]>current)
					{
						partner = t;
    					heaviest = vtxWghtArr[t] + vtxWghtArr[current];  
					}
    			
	      		}
	      		//__sync_lock_release(&nlocks[t]);
	    	}
	    	done = true;
	    	if(heaviest > 0.0)
	    	{
	    		//lock partner
	    		while(__sync_lock_test_and_set(&nlocks[partner], 1))
				 ;
				Size y = mateArr[partner];
	    		if(y != cNullItm)
	    		{
					if(vtxWghtArr[y] < vtxWghtArr[current] || (vtxWghtArr[y] == vtxWghtArr[current] && current < y))
					{
						mateArr[partner]=current;
	    				ws[partner]=heaviest;
						current = y;
	    				done= false;
					}
					else
					{
						done =false;
					}	
				}
				else
				{
					mateArr[partner]=current;
	    			ws[partner]=heaviest;	
				}
				//release the lock
				__sync_lock_release(&nlocks[partner]);
			}
		}
	}

	#pragma omp for schedule(static,256) nowait reduction(+:tempvtxWght,tempcard) 
	for (Size s = 0; s < numVtxs; ++s)
	{
 		if (mateArr[s] != cNullItm) {
      	tempvtxWght += vtxWghtArr[s];
      	++(tempcard);
    	}
	}
}
  *vtxWght =tempvtxWght;
   *card=tempcard/2;
}

template<class ItmQue>
void MatchingEngine::rPGPVtxWght(const Graph& graph,
    Size* mateArr, Size* card, Val* vtxWght) const {
  
    if (mateArr != NULL) {
      std::fill(&mateArr[0], &mateArr[graph.mNumVtxs], cNullItm);
    }
    *card = 0;
Val vtxWght1=0.0;
Val vtxWght2=0.0;
Size card1=0;
Size card2=0;
Size numVtxs = graph.mNumVtxs;
//paths counter
Size pc=0;
//mates of vertices for M1
std::vector<Size> mateVec1;
  ResizeVector<Size>(&mateVec1, graph.mNumVtxs);
  Size* mateArr1 = (graph.mNumVtxs == 0) ? NULL : &mateVec1[0];
  if (mateArr1 != NULL) {
    std::fill(&mateArr1[0], &mateArr1[graph.mNumVtxs], cNullItm);
  } 
  //mate of vertices for M2
  std::vector<Size> mateVec2;
  ResizeVector<Size>(&mateVec2, graph.mNumVtxs);
  Size* mateArr2 = (graph.mNumVtxs == 0) ? NULL : &mateVec2[0];
  if (mateArr2 != NULL) {
    std::fill(&mateArr2[0], &mateArr2[graph.mNumVtxs], cNullItm);
  } 
  //marking the visitied vertices
  std::vector<Size> markVec;
  ResizeVector<Size>(&markVec, graph.mNumVtxs);
  Size* mark = (graph.mNumVtxs == 0) ? NULL : &markVec[0];
  if (mark != NULL) {
    std::fill(&mark[0], &mark[graph.mNumVtxs], 0);
  }
const std::vector<Size>*
    vtxVecArr = (graph.mNumVtxs == 0) ? NULL : &graph.mVtxVecVec[0];
  std::vector<std::vector<Size> > tmpVtxVecVec;
  ResizeVector<std::vector<Size> >(&tmpVtxVecVec, numVtxs);
  std::vector<Size>* tmpVtxVecArr = (numVtxs == 0) ? NULL : &tmpVtxVecVec[0];
  for (Size first = 0; first < numVtxs; ++first){
    	Size sNumEdgs = vtxVecArr[first].size();
	    for (Size i = 0; i < sNumEdgs; ++i)
	    {
	    	tmpVtxVecArr[first].push_back(0);
		}  
	  }
  const Val* vtxWghtArr = (graph.mNumVtxs == 0) ? NULL : &graph.mVtxWghtVec[0];
  //to alternate between M1 and M2
  bool M1turn =true;
  bool cycle =false;
  //length of the path
  Size lp=0;
  for (Size vfirst = 0; vfirst < numVtxs; ++vfirst){
    Size first = vfirst;
    cycle=false;
    pc++;
    mark[first]=pc;
    lp=0;
    M1turn=true;
    bool done = false;
    // construct a path from first vertex
    while(!done)
    {
    // the path could not be extended
    if(M1turn && mateArr1[first]!=cNullItm)
    break;
    if(!M1turn && mateArr2[first]!=cNullItm)
    break;
    	Size sNumEdgs = vtxVecArr[first].size();
    	Size second=cNullItm;
    	Val heaviest = 0.0;
    	Size eIndex=cNullItm;
	    const Size* currentVtxArr = (sNumEdgs == 0) ? NULL : &vtxVecArr[first][0];
	    // find a heaviest neighbor to extend a path
	    for (Size i = 0; i < sNumEdgs; ++i) {
	      Size tempSecond = currentVtxArr[i];
	      
	      assert(first != tempSecond);
	      // if a vertex is already visited
	      if(mark[tempSecond]==pc)
	      {
	      	if(second==vfirst && lp %2 ==1)
		     cycle=true;
		    else
		    continue;
		  }
		  // already constructed path
		   if(mark[tempSecond]<pc && mark[tempSecond]>0)
		   continue;
		   // save a neighbor if it is heavier than the previous one
	      if (vtxWghtArr[tempSecond]> heaviest && tmpVtxVecArr[first][i] ==0 ) {
	            	second= tempSecond;
    				heaviest = vtxWghtArr[tempSecond];
    				eIndex=i;
	      }
	    }
		
		done =true;
		if(heaviest >0.0)
	    {
	    	mark[second]=pc;
	    	// increment the length
	    	lp++;
	    	if(M1turn)
	    	{
	    		if (mateArr1[second] == cNullItm)
				{
	    		mateArr1[second]=first;
	    		mateArr1[first]=second;
	    		vtxWght1+=vtxWghtArr[first];
	    		vtxWght1+=vtxWghtArr[second];
	    		card1++;
	    		tmpVtxVecArr[first][eIndex]=1;
	    		Size sNumEdgs1 = vtxVecArr[second].size();
	    		 const Size* currentVtxArr1 = (sNumEdgs == 0) ? NULL : &vtxVecArr[second][0];
			    for (Size j = 0; j < sNumEdgs1; ++j) {
			      if(currentVtxArr1[j]==first)
			      {
			      	tmpVtxVecArr[second][j]=1;
			      	break;
				  }
					}
				M1turn=false;
			  }
			}
			else
	    	{
	    	if (mateArr2[second] == cNullItm)
				{
	    		mateArr2[second]=first;
	    		mateArr2[first]=second;
	    		vtxWght2+=vtxWghtArr[first];
	    		vtxWght2+=vtxWghtArr[second];
	    		card2++;
	    		tmpVtxVecArr[first][eIndex]=1;
	    		Size sNumEdgs1 = vtxVecArr[second].size();
	    		 const Size* currentVtxArr1 = (sNumEdgs == 0) ? NULL : &vtxVecArr[second][0];
			    for (Size j = 0; j < sNumEdgs1; ++j) {
			      if(currentVtxArr1[j]==first)
			      {
			      	tmpVtxVecArr[second][j]=1;
			      	break;
				  }
					}
				M1turn=true;
			  }  		
			}
			done= false;
			first=second;
			if(cycle)
			done= true;
		}	
	}
  }

 if(vtxWght1 >= vtxWght2)
 {
 	*vtxWght =vtxWght1;	
 	*card=card1;
 	 for (Size s = 0; s < numVtxs; ++s) {
    	mateArr[s] = mateArr1[s];
   }
 }
 else
 {
 	*vtxWght =vtxWght2;	
 	*card=card2;
 	 for (Size s = 0; s < numVtxs; ++s) {
    	mateArr[s] = mateArr2[s];
   }
 }
}

template<class ItmQue>
void MatchingEngine::rPGDPVtxWght(const Graph& graph,
    Size* mateArr, Size* card, Val* vtxWght) const {
  
    if (mateArr != NULL) {
      std::fill(&mateArr[0], &mateArr[graph.mNumVtxs], cNullItm);
    }
    *card = 0;
  
Size numVtxs = graph.mNumVtxs;
 
    Size midx1=0;
    Size midx2=0;
    Size pidx=0;
    Val W1=0;
    Val W2=0;
    Val temp;
	Size* tp;
	Size pc=0;
//mates of vertices for M1
std::vector<Size> M1Vec;
  ResizeVector<Size>(&M1Vec, 3*graph.mNumVtxs);
  Size* M1 = (graph.mNumVtxs == 0) ? NULL : &M1Vec[0];
  if (M1 != NULL) {
    std::fill(&M1[0], &M1[3*graph.mNumVtxs], 0);
  }
 //mates of vertices for M2
std::vector<Size> M2Vec;
  ResizeVector<Size>(&M2Vec, 3*graph.mNumVtxs);
  Size* M2 = (graph.mNumVtxs == 0) ? NULL : &M2Vec[0];
  if (M2 != NULL) {
    std::fill(&M2[0], &M2[3*graph.mNumVtxs], 0);
  }
//mark visited vertices
std::vector<Size> markVec;
  ResizeVector<Size>(&markVec, graph.mNumVtxs);
  Size* mark = (graph.mNumVtxs == 0) ? NULL : &markVec[0];
  if (mark != NULL) {
    std::fill(&mark[0], &mark[graph.mNumVtxs], 0);
  }
  
const std::vector<Size>*
    vtxVecArr = (graph.mNumVtxs == 0) ? NULL : &graph.mVtxVecVec[0];

  std::vector<std::vector<Size> > tmpVtxVecVec;
  ResizeVector<std::vector<Size> >(&tmpVtxVecVec, numVtxs);
  std::vector<Size>* tmpVtxVecArr = (numVtxs == 0) ? NULL : &tmpVtxVecVec[0];
  
  for (Size first = 0; first < numVtxs; ++first){
    	Size sNumEdgs = vtxVecArr[first].size();
	    for (Size i = 0; i < sNumEdgs; ++i)
	    {
	    	tmpVtxVecArr[first].push_back(0);
		}
	      
	  }

  const Val* vtxWghtArr = (graph.mNumVtxs == 0) ? NULL : &graph.mVtxWghtVec[0];
	 
    bool cycle =true;
    //edge number in apth
    Size eNumPath=0;
    //path ength
    Size lp=0;
	
	for (Size vfirst = 0; vfirst < numVtxs; ++vfirst){
    Size first = vfirst;
    
    bool done = false;
    midx1=0;
    midx2=0;
    W1=0;
    W2=0;
    eNumPath=0;
    pc++;
    cycle=false;
    mark[first]=pc;
    lp=0;
    while(!done)
    {
    	if(mateArr[first]!=cNullItm)
    	break;

    	Size sNumEdgs = vtxVecArr[first].size();
    	Size second=cNullItm;
    	Val heaviest = 0.0;
    	Size eIndex=cNullItm;
	    const Size* currentVtxArr = (sNumEdgs == 0) ? NULL : &vtxVecArr[first][0];
	    for (Size i = 0; i < sNumEdgs; ++i) {
	      Size tempSecond = currentVtxArr[i];
	      
	      assert(first != tempSecond);
	      // check if the vertex is already visited
	      if(mark[tempSecond]==pc)
	      {
	      	if(second==vfirst && lp %2 ==1)
		     cycle=true;
		    else
		    continue;
		  }
		   if(mateArr[tempSecond]!=cNullItm)
		   continue;
		   //the vertex belong to different path
		   if(mark[tempSecond]<pc && mark[tempSecond]>0)
		   continue;
	    
	      if (vtxWghtArr[tempSecond]> heaviest && tmpVtxVecArr[first][i] ==0 ) {
	            	second= tempSecond;
    				heaviest = vtxWghtArr[tempSecond];
    				eIndex=i;
	      }
	    }
		
		done =true;
		if(heaviest >0.0)
	    {
	    	mark[second]=pc;
	    	//increment the path length
	    	eNumPath++;
	    	if(eNumPath==1)
	    	{
	    		// the first edge to be included in a path
				W1=heaviest+vtxWghtArr[first];
            	M1[midx1]=first;
            	M1[midx1+1]=second;
            	M1[midx1+2]=heaviest;
            	midx1+=3;
	    		tmpVtxVecArr[first][eIndex]=1;
	    		Size sNumEdgs1 = vtxVecArr[second].size();
	    		const Size* currentVtxArr1 = (sNumEdgs == 0) ? NULL : &vtxVecArr[second][0];
			    for (Size j = 0; j < sNumEdgs1; ++j) {
			      if(currentVtxArr1[j]==first)
			      {
			      	tmpVtxVecArr[second][j]=1;
			      	break;
				  }
				}
			}
			else
	    	{
	    		// using dynamic programming to select the heaviest matching edges in a path
				if(W2 + heaviest+vtxWghtArr[first] > W1)
	            {
	                temp=W2;
	                W2=W1;
	                W1=temp+heaviest+vtxWghtArr[first];
	
	                M2[midx2]=first;
	                M2[midx2+1]=second;
	                M2[midx2+2]=heaviest+vtxWghtArr[first];
	                midx2+=3;
	
	                tp=M2;
	                M2=M1;
	                M1=tp;
	
	                pidx=midx2;
	                midx2=midx1;
	                midx1=pidx;
	
	            }
	            else
	            {
	                W2=W1;
	                for(Size i=0;i<midx1;i++)
	                 M2[i]=M1[i];
	                midx2=midx1;
	            }
	    	
	    		tmpVtxVecArr[first][eIndex]=1;
	    		Size sNumEdgs1 = vtxVecArr[second].size();
	    		 const Size* currentVtxArr1 = (sNumEdgs == 0) ? NULL : &vtxVecArr[second][0];
			    for (Size j = 0; j < sNumEdgs1; ++j) {
			      if(currentVtxArr1[j]==first)
			      {
			      	tmpVtxVecArr[second][j]=1;
			      	break;
				  }
				} 		
			}
			eNumPath++;
			done= false;
			first=second;
			if(cycle)
			done=true;
		}	
	}
  	
  	//matched vertices found by using dynamic programming
  	Size u=0,v=0;
	if(midx1>0)
        {
            *vtxWght+=W1;
           for(Size j=0;j< midx1;j+=3)
           {
               u=M1[j];
               v=M1[j+1];  
                mateArr[u]=v;
                mateArr[v]=u;
                *card=*card+1;
           }
        }	
  }
}

template<class ItmQue>
void MatchingEngine::rLocallyDominantVtxWght(const Graph& graph,
    Size* mateArr, Size* card, Val* vtxWght) const {
    if (mateArr != NULL) {
      std::fill(&mateArr[0], &mateArr[graph.mNumVtxs], cNullItm);
    }
    *card = 0;
  
Size numVtxs = graph.mNumVtxs;
 // live vertices are considered (they are not matched)
std::vector<Size> aliveVec;
  ResizeVector<Size>(&aliveVec, graph.mNumVtxs);
  Size* alive = (graph.mNumVtxs == 0) ? NULL : &aliveVec[0];
  if (alive != NULL) {
    std::fill(&alive[0], &alive[graph.mNumVtxs], 1);
  }
  
const std::vector<Size>*
    vtxVecArr = (graph.mNumVtxs == 0) ? NULL : &graph.mVtxVecVec[0];

   const Val* vtxWghtArr = (graph.mNumVtxs == 0) ? NULL : &graph.mVtxWghtVec[0];      
   bool done = false;
    
    while(!done)
    {
      
	    for (Size vfirst = 0; vfirst < numVtxs; ++vfirst){
 
    	Val heaviest=0.0;
    	//skip if the vertex is already matched		
    	if(alive[ vfirst]==0 || mateArr[vfirst]!=cNullItm)
    	continue;

    	Size sNumEdgs = vtxVecArr[vfirst].size();
    	Size second = cNullItm;
	    const Size* currentVtxArr = (sNumEdgs == 0) ? NULL : &vtxVecArr[vfirst][0];
	    for (Size i = 0; i < sNumEdgs; ++i) {
	      Size tempSecond = currentVtxArr[i];
	      
	      assert(vfirst != tempSecond);
		   if(alive[tempSecond]==0)
		   continue;
	      if (vtxWghtArr[tempSecond] > heaviest ||(vtxWghtArr[tempSecond]==heaviest && tempSecond>second)) {
	            	second= tempSecond;
    				heaviest = vtxWghtArr[tempSecond];
	      }
	    }
		
		done =true;
		if(heaviest >0.0)
	    { 
	    	mateArr[vfirst]=second;
		}
		else
	    {		
	    	alive[vfirst]=0;
	   	}
		
		}
		
		for (Size vfirst = 0; vfirst < numVtxs; ++vfirst){
			 if(alive[vfirst]==0)
                    continue;
                    
            Size curr =mateArr[vfirst];
            if(mateArr[curr]!= vfirst)
            {
            	//if the edge is not locally dominant then reset mate to be considered in the next iteration
            	mateArr[vfirst]=cNullItm;
            	done=false;
			}
			else
			{
				// match an edge
				*vtxWght+=(vtxWghtArr[vfirst]+vtxWghtArr[curr]);
				*card=*card+1;
				mateArr[vfirst]=curr;
				mateArr[curr]=vfirst;
				alive[vfirst]=0;
				alive[curr]=0;
			}   
		}	
	}
}


template<class ItmQue>
void MatchingEngine::rComputeHalfVtxWghtMatching(
    const Graph& graph, Size* mateArr, Size* card,
    Val* vtxWght) const {
  
	Size numVtxs =graph.mNumVtxs;
	
	if (mateArr != NULL) {
      std::fill(&mateArr[0], &mateArr[graph.mNumVtxs], cNullItm);
    }
	
  	*card = 0; 
   	const std::vector<Size>* vtxVecArr = (graph.mNumVtxs == 0) ? NULL : &graph.mVtxVecVec[0];
    const Val* vtxWghtArr = (graph.mNumVtxs == 0) ? NULL : &graph.mVtxWghtVec[0];
  	std::vector<std::pair<Size, Val> > sExpsdLst2;
  	ResizeVector<std::pair<Size, Val> >(&sExpsdLst2, numVtxs);
  	std::pair<Size,Val> * sExpsdLst2Arr=  (numVtxs == 0) ? NULL : &sExpsdLst2[0]; 

  	for (Size first = 0; first < numVtxs; ++first) {
    	sExpsdLst2Arr[first]=std::pair<Size, Val>(first, vtxWghtArr[first]); 
  	} 
  	//sort vertices in non-increasing order of their weights
  	std::sort(sExpsdLst2.begin(),sExpsdLst2.end(),ValGreater<std::pair<Size, Val> >());

    for (Size current = 0; current < numVtxs; ++current) {
    	Size sFirst =  sExpsdLst2[current].first;
    	if(mateArr[sFirst]!=cNullItm)
    		continue;
    	Val heaviest = 0.0;
    	Size mateT = cNullItm;
    	Size numEdgs = vtxVecArr[sFirst].size();
    	const Size* vtxArr = (numEdgs == 0) ? NULL : &vtxVecArr[sFirst][0];
    	//find  heaviest neighbor
    	for (Size i = 0; i < numEdgs; ++i) {
      		Size t = vtxArr[i];
	  		if(mateArr[t]== cNullItm && vtxWghtArr[t] > heaviest)
	  		{
				heaviest = 	vtxWghtArr[t];
				mateT=t;
	  		}
		}
		if(heaviest > 0.0)
		{
			//match an edge
			mateArr[mateT]=sFirst;
			mateArr[sFirst]=mateT;
			++(*card);	
		}
	
	}
  *vtxWght = 0.0;
  for (Size s = 0; s < numVtxs; ++s) {
    if (mateArr[s] != cNullItm) {
      *vtxWght += vtxWghtArr[s];
    }
  }
}


template<class ItmQue>
void MatchingEngine::rComputeTwoThirdVtxWghtMatching(
    const Graph& graph, Size* mateArr, Size* card,
    Val* vtxWght) const {
 Size numVtxs =graph.mNumVtxs;
 
 if (mateArr != NULL) {
      std::fill(&mateArr[0], &mateArr[graph.mNumVtxs], cNullItm);
    }
    

  *card = 0;
  
   const std::vector<Size>* vtxVecArr = (graph.mNumVtxs == 0) ? NULL : &graph.mVtxVecVec[0];

    const Val* vtxWghtArr = (graph.mNumVtxs == 0) ? NULL : &graph.mVtxWghtVec[0];
  
  // a queue to save matched neighbors
  std::deque<Size> sPrcsbQue;
 
  std::vector<std::pair<Size, Val> > sExpsdLst2;
  ResizeVector<std::pair<Size, Val> >(&sExpsdLst2, numVtxs);
  std::pair<Size,Val> * sExpsdLst2Arr=  (numVtxs == 0) ? NULL : &sExpsdLst2[0]; 

  for (Size first = 0; first < numVtxs; ++first) {
    sExpsdLst2Arr[first]=std::pair<Size, Val>(first, vtxWghtArr[first]); 
  } 
  //sort vertices in non-increasing order
  std::sort(sExpsdLst2.begin(),sExpsdLst2.end(),ValGreater<std::pair<Size, Val> >());
    // for each vertex find an augmenting path of length at most 3 that reach a heaviest unmatched vertex
    for (Size current = 0; current < numVtxs; ++current) {
		Size sFirst =  sExpsdLst2Arr[current].first;    	
		if(mateArr[sFirst]!=cNullItm)
    	continue;
    	Size numEdgs = vtxVecArr[sFirst].size();
    	const Size* vtxArr = (numEdgs == 0) ? NULL : &vtxVecArr[sFirst][0];
		Val heaviest =0;
		Size mateT =cNullItm;
		Size tempMate=cNullItm;
		sPrcsbQue.clear();
    	for (Size i = 0; i < numEdgs; ++i) {
      		Size t = vtxArr[i];
	  		if(mateArr[t]== cNullItm)
	  		{
	  			if(vtxWghtArr[t] > heaviest)
	  			{
	  				heaviest = vtxWghtArr[t];
	  				tempMate=cNullItm;
	  				mateT = t;	
				}
	  		}
	  		else
	  		{
	  			sPrcsbQue.push_back(mateArr[t]);
	  		}
		}
		while (sPrcsbQue.empty() == false) {
		  	Size s = sPrcsbQue.front();
		  	sPrcsbQue.pop_front();
		  	Size numEdgs = vtxVecArr[s].size();
			const Size* vtxArr = (numEdgs == 0) ? NULL : &vtxVecArr[s][0];
			for (Size i = 0; i < numEdgs; ++i) {
				Size t = vtxArr[i];
			  	if(mateArr[t] == cNullItm && t != sFirst)
			  	{
			  		if(vtxWghtArr[t] > heaviest)
			  		{
			  		heaviest = vtxWghtArr[t];
	  				mateT = t;
					tempMate=mateArr[s];
					}
				
			  	}
			}
		}
		if(heaviest > 0)
		{
			//augment an augmenting path
			if(tempMate == cNullItm )
			{
				mateArr[mateT]=sFirst;
				mateArr[sFirst]=mateT;
				++(*card);
			}
			else
			{
				Size tempS = mateArr[tempMate];
				mateArr[mateT]=tempS;
				mateArr[tempS]=mateT;
				mateArr[tempMate]=sFirst;
				mateArr[sFirst]=tempMate;
				++(*card);
			}
		}
		
	}
  	*vtxWght = 0.0;
  	for (Size s = 0; s < numVtxs; ++s) {
    	if (mateArr[s] != cNullItm) {
      	*vtxWght += vtxWghtArr[s];
    	}
  	}
}



// rProcessBlsmEdmonds: finds proper LCA, then delegates to rProcessBlsm3.
// rProcessBlsm3 assumes bs2 is the ancestor. We find the real LCA and
// arrange args so that the ancestor side is s2 in the rProcessBlsm3 call.
template <class ItmQue, class ItmStk>
void MatchingEngine::rProcessBlsmEdmonds(const Graph& graph, Size* mateArr,
    Size* ptrArr, Size* blsmArr, Size* rprsArr, Size* blsmParent,
    std::vector<Size>* blsmChildren, std::vector<Size>* blsmL,
    Size* label, Stt* sttArr, ItmStk& prcsbStk,
    Size s1, Size s2, Size* blsm, Val* dualArr,
    Size* blsmin, Size* treeNumArr, Size maxBlsm,
    std::stack<Size>& freeBlsmStk, Size* newBlsmId) const {

  if (newBlsmId != NULL)
    *newBlsmId = cNullItm;

  Size numVtxs = graph.mNumVtxs;

  Size bs1 = s1;
  Size bs2 = s2;
  if(blsmParent[s1] != cNullItm)
    bs1 = rFindblsm<ItmQue>(blsmParent, blsmArr, rprsArr, s1);
  if(blsmParent[s2] != cNullItm)
    bs2 = rFindblsm<ItmQue>(blsmParent, blsmArr, rprsArr, s2);
  if(bs1 == bs2) return;

  // --- Find LCA via two-sided walk through ptrArr ---
  std::vector<bool> vis(2 * numVtxs, false);
  vis[bs1] = true;
  vis[bs2] = true;
  Size Ba = bs1, Bb = bs2;
  Size LCA = cNullItm;

  while(LCA == cNullItm)
  {
    if(Ba != cNullItm)
    {
      Size baseA = (Ba >= numVtxs) ? rprsArr[Ba] : Ba;
      if(ptrArr[baseA] != cNullItm)
      {
        Ba = ptrArr[baseA];
        if(blsmParent[Ba] != cNullItm)
          Ba = rFindblsm<ItmQue>(blsmParent, blsmArr, rprsArr, Ba);
        if(vis[Ba]) { LCA = Ba; break; }
        vis[Ba] = true;
      }
      else Ba = cNullItm;
    }
    if(Bb != cNullItm && LCA == cNullItm)
    {
      Size baseB = (Bb >= numVtxs) ? rprsArr[Bb] : Bb;
      if(ptrArr[baseB] != cNullItm)
      {
        Bb = ptrArr[baseB];
        if(blsmParent[Bb] != cNullItm)
          Bb = rFindblsm<ItmQue>(blsmParent, blsmArr, rprsArr, Bb);
        if(vis[Bb]) { LCA = Bb; break; }
        vis[Bb] = true;
      }
      else Bb = cNullItm;
    }
    if(Ba == cNullItm && Bb == cNullItm) break;
  }
  if(LCA == cNullItm) return;


  // --- Delegate to rProcessBlsm3 ---
  if(LCA == bs2)
  {
    rProcessBlsm3<ItmQue, ItmStk>(graph, mateArr, ptrArr, blsmArr, rprsArr,
        blsmParent, blsmChildren, blsmL, label, sttArr, prcsbStk,
        s1, s2, blsm, dualArr, blsmin, treeNumArr, maxBlsm, &freeBlsmStk,
        newBlsmId);
  }
  else if(LCA == bs1)
  {
    rProcessBlsm3<ItmQue, ItmStk>(graph, mateArr, ptrArr, blsmArr, rprsArr,
        blsmParent, blsmChildren, blsmL, label, sttArr, prcsbStk,
        s2, s1, blsm, dualArr, blsmin, treeNumArr, maxBlsm, &freeBlsmStk,
        newBlsmId);
  }
  else
  {
    // LCA is a third node. Build manually with two half-walks.
    std::queue<Size> tempQ;
    Size b = cNullItm;
    if(rAcquireBlsmId(blsm, numVtxs, maxBlsm, blsmChildren, blsmL,
        blsmin, &freeBlsmStk, &b) == false) {
      return;
    }
    if (newBlsmId != NULL)
      *newBlsmId = b;
    Size shiftB = b - numVtxs;
    Size base = (LCA >= numVtxs) ? rprsArr[LCA] : LCA;

    rprsArr[b] = base;
    treeNumArr[b] = treeNumArr[LCA];
    blsmParent[b] = cNullItm;
    dualArr[b] = 0;
    label[b] = 1;
    blsmChildren[shiftB].push_back(LCA);
    // DEFERRED: blsmParent[LCA] = b is set after the arm walks below,
    // because the arm walks call rFindblsm which would see LCA as part
    // of the not-yet-constructed blossom b and fail to terminate.

    // Arm from bs1 to LCA (will be reversed)
    std::vector<Size> arm1_ch, arm1_lk;
    {
      Size cur = bs1;
      while(cur != LCA)
      {
        arm1_ch.push_back(cur);
        Size baseCur = (cur >= numVtxs) ? rprsArr[cur] : cur;
        Size mateCur = mateArr[baseCur];
        arm1_lk.push_back(baseCur);
        arm1_lk.push_back(mateCur);

        Size Tblsm = mateCur;
        if(blsmParent[mateCur] != cNullItm)
          Tblsm = rFindblsm<ItmQue>(blsmParent, blsmArr, rprsArr, mateCur);
        if(Tblsm == LCA) break;

        arm1_ch.push_back(Tblsm);
        Size exitT;
        if(Tblsm >= numVtxs && blsmin[Tblsm - numVtxs] != cNullItm)
          exitT = blsmin[Tblsm - numVtxs];
        else
          exitT = mateCur;
        Size nextS = ptrArr[baseCur];
        arm1_lk.push_back(exitT);
        arm1_lk.push_back(nextS);

        Size nextSB = nextS;
        if(blsmParent[nextS] != cNullItm)
          nextSB = rFindblsm<ItmQue>(blsmParent, blsmArr, rprsArr, nextS);
        if(nextSB == LCA) break;
        cur = nextSB;
      }
    }

    // Arm from bs2 to LCA (in order)
    std::vector<Size> arm2_ch, arm2_lk;
    {
      Size cur = bs2;
      while(cur != LCA)
      {
        arm2_ch.push_back(cur);
        Size baseCur = (cur >= numVtxs) ? rprsArr[cur] : cur;
        Size mateCur = mateArr[baseCur];
        arm2_lk.push_back(baseCur);
        arm2_lk.push_back(mateCur);

        Size Tblsm = mateCur;
        if(blsmParent[mateCur] != cNullItm)
          Tblsm = rFindblsm<ItmQue>(blsmParent, blsmArr, rprsArr, mateCur);
        if(Tblsm == LCA) break;

        arm2_ch.push_back(Tblsm);
        Size exitT;
        if(Tblsm >= numVtxs && blsmin[Tblsm - numVtxs] != cNullItm)
          exitT = blsmin[Tblsm - numVtxs];
        else
          exitT = mateCur;
        Size nextS = ptrArr[baseCur];
        arm2_lk.push_back(exitT);
        arm2_lk.push_back(nextS);

        Size nextSB = nextS;
        if(blsmParent[nextS] != cNullItm)
          nextSB = rFindblsm<ItmQue>(blsmParent, blsmArr, rprsArr, nextS);
        if(nextSB == LCA) break;
        cur = nextSB;
      }
    }

    Size p = arm1_ch.size();
    Size q = arm2_ch.size();

    // Now safe to set blsmParent[LCA] (arm walks are done)
    blsmParent[LCA] = b;

    for(int i = (int)p - 1; i >= 0; i--) {
      blsmChildren[shiftB].push_back(arm1_ch[i]);
      blsmParent[arm1_ch[i]] = b;
    }
    for(Size i = 0; i < q; i++) {
      blsmChildren[shiftB].push_back(arm2_ch[i]);
      blsmParent[arm2_ch[i]] = b;
    }

    Size totalKids = blsmChildren[shiftB].size();
    blsmL[shiftB].resize(2 * totalKids, cNullItm);

    if(p > 0) {
      blsmL[shiftB][0] = arm1_lk.back();
      blsmL[shiftB][1] = arm1_lk[arm1_lk.size()-2];
    } else if(q > 0) {
      blsmL[shiftB][0] = s1;
      blsmL[shiftB][1] = s2;
    }

    for(Size i = 0; i + 1 < p; i++) {
      Size ci = 1 + i;
      Size j = p - 2 - i;
      blsmL[shiftB][2*ci]     = arm1_lk[2*j + 1];
      blsmL[shiftB][2*ci + 1] = arm1_lk[2*j];
    }

    {
      Size ci = p;
      blsmL[shiftB][2*ci]   = s1;
      blsmL[shiftB][2*ci+1] = s2;
    }

    for(Size j = 0; j + 1 < q; j++) {
      Size ci = p + 1 + j;
      blsmL[shiftB][2*ci]   = arm2_lk[2*j];
      blsmL[shiftB][2*ci+1] = arm2_lk[2*j + 1];
    }

    {
      Size ci = totalKids - 1;
      if(q > 0) {
        blsmL[shiftB][2*ci]   = arm2_lk[arm2_lk.size()-2];
        blsmL[shiftB][2*ci+1] = arm2_lk.back();
      } else if(p > 0) {
        blsmL[shiftB][2*ci]   = s2;
        blsmL[shiftB][2*ci+1] = s1;
      }
    }

    for(Size ci = 0; ci < totalKids; ci++) {
      Size child = blsmChildren[shiftB][ci];
      label[child] = 1;
      if(child < numVtxs) {
        if(sttArr[child] == eSttIdle || label[child] == 2 || label[child] == 0) {
          prcsbStk.Push(child);
          sttArr[child] = eSttDfsPrcsb;
          ptrArr[child] = ptrArr[base];
        }
      } else {
        if(blsmin[child - numVtxs] == cNullItm || blsmin[child - numVtxs] == rprsArr[child])
          ptrArr[rprsArr[child]] = ptrArr[base];
        else
          ptrArr[blsmin[child - numVtxs]] = ptrArr[base];

        if(sttArr[rprsArr[child]] == eSttIdle) {
          tempQ.push(child - numVtxs);
          while(tempQ.empty() == false) {
            Size tempB = tempQ.front(); tempQ.pop();
            for(Size i = 0; i < blsmChildren[tempB].size(); i++) {
              Size sub = blsmChildren[tempB][i];
              if(sub < numVtxs) {
                if(sttArr[sub] == eSttIdle) {
                  prcsbStk.Push(sub);
                  sttArr[sub] = eSttDfsPrcsb;
                  ptrArr[sub] = ptrArr[base];
                }
              } else tempQ.push(sub - numVtxs);
            }
          }
        }
        label[child] = 1;
      }
    }
  }

  // Debug print for small graphs
}



// rAugmentEdmonds: same as rAugment but fixes nested blossom bug.
// rAugment line 1423 passes blsmArr[t] (lowest-level blossom) to augmentBlsm.
// For nested blossoms this is wrong; we pass blsmt (top-level) instead.
template<class ItmIdxdQue, class ItmQue>
void MatchingEngine::rAugmentEdmonds(Size* mateArr, Size* ptrArr,
    ItmIdxdQue& expsdQue, Size sLast, Size tLast,
    Size* blsmArr, Size* rprsArr, Size* blsmParent,
    std::vector<Size>* blsmChildren, std::vector<Size>* blsmL,
    Size num, Size* label, Size* blsmin) const {

  Size s = sLast;
  Size t = tLast;
  while(s != cNullItm) {
    Size blsmt = t;
    if(blsmParent[t] != cNullItm)
      blsmt = rFindblsm<ItmQue>(blsmParent, blsmArr, rprsArr, t);
    label[blsmt] = 0;
    if(blsmt >= num && blsmin[blsmt - num] != cNullItm)
    {
      Size tempT = blsmin[blsmt - num];
      augmentBlsm<ItmQue>(blsmt, tempT, mateArr, blsmArr, rprsArr,
                           blsmParent, blsmChildren, blsmL, num);
      t = tempT;
    }

    Size blsms = s;
    if(blsmParent[s] != cNullItm)
      blsms = rFindblsm<ItmQue>(blsmParent, blsmArr, rprsArr, s);
    if(blsms >= num && rprsArr[blsms] != s)
    {
      label[blsms] = 0;
      Size base = rprsArr[blsms];
      Size mateBase = mateArr[base];
      augmentBlsm<ItmQue>(blsms, s, mateArr, blsmArr, rprsArr,
                           blsmParent, blsmChildren, blsmL, num);
      label[s] = 0;
      label[t] = 0;
      mateArr[s] = t;
      mateArr[t] = s;
      s = base;
      s = ptrArr[s];
      t = mateBase;
    }
    else
    {
      Size tt = mateArr[s];
      label[s] = 0;
      label[t] = 0;
      mateArr[s] = t;
      mateArr[t] = s;
      s = ptrArr[s];
      t = tt;
    }
  }
}

// rAugmentOneSide: augment one side of a cross-tree augmenting path.
// Walks from sStart up to its tree root via ptrArr, flipping matched edges.
// Does NOT set mateArr[sStart] — the caller sets the cross-tree edge after both sides.
//
// Strategy: Use EXACTLY the same code path as rAugmentEdmonds, but pass a
// dummy "virtual" tLast. We observe that rAugmentEdmonds(sLast, tLast) does:
//   iteration 1: match sLast <-> tLast, advance
//   iteration 2+: match s <-> t, advance
// For one-sided augmentation, we want to skip the matching of sLast but still
// walk the same path. We achieve this by running rAugmentEdmonds logic but
// NOT writing the first mateArr assignment, then writing mateArr for all
// subsequent steps.
template<class ItmIdxdQue, class ItmQue>
void MatchingEngine::rAugmentOneSide(Size* mateArr, Size* ptrArr,
    Size sStart,
    Size* blsmArr, Size* rprsArr, Size* blsmParent,
    std::vector<Size>* blsmChildren, std::vector<Size>* blsmL,
    Size num, Size* label, Size* blsmin) const {

  // Exactly mirrors rAugmentEdmonds but:
  // - First iteration: DON'T set mateArr[s]=t or mateArr[t]=s (skip the cross-tree edge)
  // - All subsequent iterations: normal

  Size s = sStart;
  Size t = cNullItm;  // dummy — won't be used in first iteration
  bool firstIter = true;

  while(s != cNullItm) {
    if(!firstIter) {
      // Handle T-blossom containing t
      Size blsmt = t;
      if(blsmParent[t] != cNullItm)
        blsmt = rFindblsm<ItmQue>(blsmParent, blsmArr, rprsArr, t);
      label[blsmt] = 0;
      if(blsmt >= num && blsmin[blsmt - num] != cNullItm) {
        Size tempT = blsmin[blsmt - num];
        augmentBlsm<ItmQue>(blsmt, tempT, mateArr, blsmArr, rprsArr,
                             blsmParent, blsmChildren, blsmL, num);
        t = tempT;
      }
    }

    // Handle S-blossom containing s
    Size blsms = s;
    if(blsmParent[s] != cNullItm)
      blsms = rFindblsm<ItmQue>(blsmParent, blsmArr, rprsArr, s);
    if(blsms >= num && rprsArr[blsms] != s) {
      label[blsms] = 0;
      Size base = rprsArr[blsms];
      Size mateBase = mateArr[base];
      augmentBlsm<ItmQue>(blsms, s, mateArr, blsmArr, rprsArr,
                           blsmParent, blsmChildren, blsmL, num);
      if(!firstIter) {
        label[s] = 0;
        label[t] = 0;
        mateArr[s] = t;
        mateArr[t] = s;
      }
      s = ptrArr[base];
      t = mateBase;
    }
    else
    {
      Size tt = mateArr[s];
      if(!firstIter) {
        label[s] = 0;
        label[t] = 0;
        mateArr[s] = t;
        mateArr[t] = s;
      }
      s = ptrArr[s];
      t = tt;
    }
    firstIter = false;
  }
}


// rSplitBlossom: expand a T-blossom in-tree (D4 event).
// When a T-blossom b has z=0, we dissolve it and relabel its cycle children.
// The entry child (containing blsmin[b]) keeps T-label.
// Walking the blossom cycle from entry to base, we alternate T/S labels on pairs.
// The base child gets S-label and its vertices are pushed for scanning.
// Children that become S get their vertices pushed onto prcsbStk.
// Children that become T get blsmin set to the link vertex that enters them.
template<class ItmQue, class ItmStk>
void MatchingEngine::rSplitBlossom(const Graph& graph,
    Size b, Size* mateArr, Size* ptrArr,
    Size* blsmArr, Size* rprsArr, Size* blsmParentArr,
    std::vector<Size>* blsmChildrenArr, std::vector<Size>* blsmL,
    Size* label, Stt* sttArr, ItmStk& prcsbStk,
    Val* dualArr, Size* blsmin, Size* treeNumArr,
    std::stack<Size>& freeBlsmStk,
    Size numVtxs) const
{
  Size shiftB = b - numVtxs;
  std::vector<Size>& kids  = blsmChildrenArr[shiftB];
  std::vector<Size>& links = blsmL[shiftB];
  Size k = kids.size();
  if (k == 0) return;

  Size bTreeNum = treeNumArr[b];
  Size entryVtx = blsmin[shiftB]; // vertex through which the tree entered this T-blossom

  // Find which child contains the entry vertex by walking up from entryVtx
  Size entryChildIdx = 0;
  {
    // Walk up parent chain from entryVtx to find the direct child of b
    Size cursor = entryVtx;
    while (cursor != cNullItm && blsmParentArr[cursor] != b)
      cursor = blsmParentArr[cursor];
    // cursor is now the direct child of b containing entryVtx
    // (could be entryVtx itself if it's a direct vertex child)
    if (cursor == cNullItm) cursor = entryVtx; // fallback
    for (Size ci = 0; ci < k; ++ci)
    {
      if (kids[ci] == cursor) { entryChildIdx = ci; break; }
    }
  }

  // Find which child contains the base vertex (rprsArr[b])
  Size baseV = rprsArr[b];
  Size baseChildIdx = 0;
  {
    Size cursor = baseV;
    while (cursor != cNullItm && blsmParentArr[cursor] != b)
      cursor = blsmParentArr[cursor];
    if (cursor == cNullItm) cursor = baseV;
    for (Size ci = 0; ci < k; ++ci)
    {
      if (kids[ci] == cursor) { baseChildIdx = ci; break; }
    }
  }

  // Reset cached top-level blossom ids for every vertex contained in b.
  // After the split, vertices should rediscover their new top-level blossom
  // by climbing blsmParentArr rather than reusing the dissolved shell id.
  {
    std::stack<Size> resetStk;
    for (Size ci = 0; ci < k; ++ci)
      resetStk.push(kids[ci]);

    while (!resetStk.empty())
    {
      Size cur = resetStk.top();
      resetStk.pop();
      if (cur < numVtxs)
      {
        blsmArr[cur] = cur;
      }
      else
      {
        Size shiftCur = cur - numVtxs;
        for (Size i = 0; i < blsmChildrenArr[shiftCur].size(); ++i)
          resetStk.push(blsmChildrenArr[shiftCur][i]);
      }
    }
  }

  // Dissolve: detach children from b
  rprsArr[b] = cNullItm;
  dualArr[b] = 0;
  label[b] = 0;
  treeNumArr[b] = cNullItm;
  blsmin[shiftB] = cNullItm;

  for (Size ci = 0; ci < k; ++ci)
  {
    Size child = kids[ci];
    blsmParentArr[child] = cNullItm;
    if (child < numVtxs)
      blsmArr[child] = child;
  }

  // Now relabel: entry child is T, then alternate T/S in pairs walking toward base.
  // Direction: we pick the direction (CW or CCW) from entry toward base that has even distance.
  // (blossom cycle has odd length k; one direction gives even dist, other gives odd)
  Size distCW = 0;
  {
    Size ci = entryChildIdx;
    while (ci != baseChildIdx) { ci = (ci + 1) % k; ++distCW; }
  }
  Size distCCW = 0;
  {
    Size ci = entryChildIdx;
    while (ci != baseChildIdx) { ci = (ci + k - 1) % k; ++distCCW; }
  }

  // The path from entry to base should have even length so that alternation (T,S,T,S,...) 
  // ends with base getting S-label. In a blossom with odd-length cycle, one direction
  // gives even distance, the other odd. Pick the even one.
  int inc;
  if (distCW % 2 == 0)
    inc = 1;  // clockwise
  else
    inc = -1; // counterclockwise

  std::queue<Size> tempQ;

  // Helper lambda-like: push all vertices in a child (possibly blossom) as S-processable
  // (We can't use lambdas in pre-C++11, so we'll inline the logic)

  // Label the entry child as T
  {
    Size child = kids[entryChildIdx];
    label[child] = 2;
    treeNumArr[child] = bTreeNum;
    if (child >= numVtxs)
      blsmin[child - numVtxs] = entryVtx;
    else
      sttArr[child] = eSttIdle; // T-vertex, not scanned
  }

  // Walk from entry toward base, alternating S/T labels.
  // Pattern: entry=T (done above), then S, T, S, T, ..., base=S
  Size ci = entryChildIdx;
  while (ci != baseChildIdx)
  {
    // Advance to next child — this one becomes S
    Size prevCi = ci;
    ci = (Size)(((int)ci + inc + (int)k) % (int)k);
    if (ci == baseChildIdx) break;

    Size childS = kids[ci]; // S-labeled
    label[childS] = 1;
    treeNumArr[childS] = bTreeNum;

    // Set ptrArr: the mate of the previous T-child's base points here
    {
      Size edgeIdx = (inc == 1) ? prevCi : ci;
      edgeIdx = edgeIdx % k;
      // This S-child was reached from the T-child before it
      // ptrArr for this S-vertex should point to the link vertex in the prev T-child
      Size fromVtx = links[2 * edgeIdx + (inc == 1 ? 0 : 1)]; // vertex in prev (T) child
      if (childS < numVtxs)
      {
        ptrArr[childS] = fromVtx;
      }
      else
      {
        // For S-blossom child, ptrArr of base points back
        ptrArr[rprsArr[childS]] = fromVtx;
      }
    }

    if (childS >= numVtxs)
    {
      tempQ.push(childS - numVtxs);
      while (tempQ.empty() == false)
      {
        Size tempB = tempQ.front(); tempQ.pop();
        for (Size i = 0; i < blsmChildrenArr[tempB].size(); ++i)
        {
          Size sub = blsmChildrenArr[tempB][i];
          if (sub < numVtxs)
          {
            if (sttArr[sub] == eSttIdle)
            {
              label[sub] = 1;
              treeNumArr[sub] = bTreeNum;
              sttArr[sub] = eSttDfsPrcsb;
              prcsbStk.Push(sub);
            }
          }
          else
          {
            label[sub] = 1;
            treeNumArr[sub] = bTreeNum;
            tempQ.push(sub - numVtxs);
          }
        }
      }
    }
    else
    {
      if (sttArr[childS] == eSttIdle)
      {
        sttArr[childS] = eSttDfsPrcsb;
        prcsbStk.Push(childS);
      }
    }

    // Advance to next child — this one becomes T
    Size prevCi2 = ci;
    ci = (Size)(((int)ci + inc + (int)k) % (int)k);
    if (ci == baseChildIdx) break;

    Size childT = kids[ci]; // T-labeled
    label[childT] = 2;
    treeNumArr[childT] = bTreeNum;
    if (childT >= numVtxs)
    {
      Size edgeIdx2 = (inc == 1) ? prevCi2 : ci;
      edgeIdx2 = edgeIdx2 % k;
      Size linkIn = links[2 * edgeIdx2 + (inc == 1 ? 1 : 0)];
      blsmin[childT - numVtxs] = linkIn;
    }
    else
    {
      sttArr[childT] = eSttIdle;
    }

    // Set ptrArr for mate of this T-child's base
    {
      Size edgeIdx2 = (inc == 1) ? prevCi2 : ci;
      edgeIdx2 = edgeIdx2 % k;
      Size fromVtx = links[2 * edgeIdx2 + (inc == 1 ? 0 : 1)]; // vertex in S child
      Size tBase = (childT >= numVtxs) ? rprsArr[childT] : childT;
      Size mateT = mateArr[tBase];
      if (mateT != cNullItm)
        ptrArr[mateT] = fromVtx;
    }
  }

  // Label the base child as S
  {
    Size child = kids[baseChildIdx];
    label[child] = 1;
    treeNumArr[child] = bTreeNum;
    if (child >= numVtxs)
    {
      tempQ.push(child - numVtxs);
      while (tempQ.empty() == false)
      {
        Size tempB = tempQ.front(); tempQ.pop();
        for (Size i = 0; i < blsmChildrenArr[tempB].size(); ++i)
        {
          Size sub = blsmChildrenArr[tempB][i];
          if (sub < numVtxs)
          {
            if (sttArr[sub] == eSttIdle)
            {
              label[sub] = 1;
              treeNumArr[sub] = bTreeNum;
              sttArr[sub] = eSttDfsPrcsb;
              prcsbStk.Push(sub);
            }
          }
          else
          {
            label[sub] = 1;
            treeNumArr[sub] = bTreeNum;
            tempQ.push(sub - numVtxs);
          }
        }
      }
    }
    else
    {
      if (sttArr[child] == eSttIdle)
      {
        sttArr[child] = eSttDfsPrcsb;
        prcsbStk.Push(child);
      }
    }
  }

  kids.clear();
  links.clear();
  freeBlsmStk.push(b);
}


inline void MatchingEngine::ComputeMaxEdgWghtMatchingEdmonds(const Graph& graph,
                                                            std::vector<Size>* mateVec,
                                                            Size* card,
                                                            Val* edgWght) const
{
  if ((mateVec == NULL) || (card == NULL) || (edgWght == NULL))
  {
    return;
  }

  Size numVtxs = graph.mNumVtxs;
  ResizeVector<Size>(mateVec, numVtxs);
  Size* mateArr = (numVtxs == 0) ? NULL : &((*mateVec)[0]);

  rComputeMaxEdgWghtMatchingEdmonds<VecItmQue, VecItmStk, LstItmIdxdQue>(graph,
                                                                                           mateArr,
                                                                                           card,
                                                                                           edgWght);
}

inline void MatchingEngine::ComputeMaxEdgWghtMatchingEdmondsFastV2(
    const Graph& graph,
    std::vector<Size>* mateVec,
    Size* card,
    Val* edgWght) const
{
  FastMatchingEngineV2 engine;
  engine.compute(graph, mateVec, card, edgWght);
}

inline void MatchingEngine::ComputeMaxEdgWghtMatchingGabow(
    const Graph& graph,
    std::vector<Size>* mateVec,
    Size* card,
    Val* edgWght) const
{
  FastMatchingEngineV2 engine(FastMatchingEngineV2::kGabow2018);
  engine.compute(graph, mateVec, card, edgWght);
}


template<class ItmQue, class ItmStk, class ItmIdxdQue>
void MatchingEngine::rComputeMaxEdgWghtMatchingEdmonds(const Graph& graph,
                                                      Size* mateArr,
                                                      Size* card,
                                                      Val* edgWght) const
{
  // Edmonds' primal-dual blossom algorithm for maximum edge-weight matching.
  // The active exact path still uses the validated rebuild-based phase-2
  // logic; experimental incremental event code is kept disabled in-source.

  if ((mateArr == NULL) || (card == NULL) || (edgWght == NULL))
    return;

  Size numVtxs = graph.mNumVtxs;
  *card = 0;
  *edgWght = 0;
  if (numVtxs == 0)
    return;

  const std::vector<std::vector<Size> >& vtxVecArr = graph.mVtxVecVec;
  const std::vector<std::vector<Val> >& edgWghtVecArr = graph.mEdgWghtVecVec;

  // ------------------------
  // Allocate/initialize state
  // ------------------------

  std::vector<Size> ptrVec;
  ResizeVector<Size>(&ptrVec, numVtxs);
  Size* ptrArr = &ptrVec[0];

  std::vector<Size> blsmVec;
  ResizeVector<Size>(&blsmVec, numVtxs);
  Size* blsmArr = &blsmVec[0];

  // Blossom arrays: Edmonds can recycle freed IDs via freeBlsmStk.
  // 20*n remains the current capacity bound.
  const Size maxBlsm = 20 * numVtxs;

  std::vector<Size> rprsVec;
  ResizeVector<Size>(&rprsVec, numVtxs + maxBlsm);
  Size* rprsArr = &rprsVec[0];

  std::vector<Size> labelVec;
  ResizeVector<Size>(&labelVec, numVtxs + maxBlsm);
  Size* label = &labelVec[0];

  std::vector<Size> blsmParentVec;
  ResizeVector<Size>(&blsmParentVec, numVtxs + maxBlsm);
  Size* blsmParentArr = &blsmParentVec[0];

  std::vector<std::vector<Size> > blsmChildrenVec;
  ResizeVector<std::vector<Size> >(&blsmChildrenVec, maxBlsm);
  std::vector<Size>* blsmChildrenArr = &blsmChildrenVec[0];

  std::vector<std::vector<Size> > blsmLVec;
  ResizeVector<std::vector<Size> >(&blsmLVec, maxBlsm);
  std::vector<Size>* blsmL = &blsmLVec[0];

  std::vector<Size> blsminVec;
  ResizeVector<Size>(&blsminVec, maxBlsm);
  Size* blsmin = &blsminVec[0];

  std::vector<Val> dualVec;
  ResizeVector<Val>(&dualVec, numVtxs + maxBlsm);
  Val* dualArr = &dualVec[0];

  std::vector<Stt> sttVec;
  ResizeVector<Stt>(&sttVec, numVtxs);
  Stt* sttArr = &sttVec[0];

  std::vector<Size> idxVec;
  ResizeVector<Size>(&idxVec, numVtxs);
  Size* idxArr = &idxVec[0];

  std::vector<Size> treeNumVec;
  ResizeVector<Size>(&treeNumVec, numVtxs + maxBlsm);
  Size* treeNumArr = &treeNumVec[0];

  // matching empty
  for (Size v = 0; v < numVtxs; ++v)
  {
    mateArr[v] = cNullItm;
  }

  // compute max edge weight (for a feasible initial dual solution)
  Val maxW = 0;
  for (Size u = 0; u < numVtxs; ++u)
  {
    const std::vector<Val>& wVec = edgWghtVecArr[u];
    for (Size i = 0; i < wVec.size(); ++i)
    {
      if (wVec[i] > maxW)
      {
        maxW = wVec[i];
      }
    }
  }
  Val initY = maxW / 2;

  //Initialization of data structures

  for (Size v = 0; v < numVtxs; ++v)
  {
    dualArr[v] = initY; // y[v]
    rprsArr[v] = v;
    blsmArr[v] = v;
    ptrArr[v] = cNullItm;
    sttArr[v] = eSttIdle;
    idxArr[v] = 0;
  }
  for (Size b = numVtxs; b < numVtxs + maxBlsm; ++b)
  {
    dualArr[b] = 0; // z[b]
    rprsArr[b] = cNullItm;
  }

  for (Size i = 0; i < numVtxs + maxBlsm; ++i)
  {
    label[i] = 0;
    blsmParentArr[i] = cNullItm;
    treeNumArr[i] = cNullItm;
  }
  for (Size i = 0; i < maxBlsm; ++i)
  {
    blsmin[i] = cNullItm;
  }



  // Current number of blossoms in use (next available blossom id).
  Size blsm = numVtxs;

  // Slack: vertex duals are updated directly (like 1-eps), so for edges
  // between different top-level blossoms, slack = y[u]+y[v]-w.
  // Blossom duals z[B] only needed for Delta4 (T-blossom expansion).

  ItmIdxdQue expsdQue(numVtxs);
  for (Size v = 0; v < numVtxs; ++v)
    expsdQue.Push(v);

  const Val INF = std::numeric_limits<Val>::max();

  // ================================================================
  // Event-driven main loop.
  //
  // Heap entries: std::pair<Val, Size> where first=key (deltaSum at
  // which the event fires), second=payload.
  //
  // D1 heap: key = dualArr[v] for free S-vertex v (fires when y_eff=0)
  // D2 heap: key = dualArr[u]+dualArr[v]-w for edge u-v, S-to-unlabeled
  //          payload = edge encoding  (u * numVtxs + edgeIdx)
  // D3 heap: key = (dualArr[u]+dualArr[v]-w)/2 for edge u-v, S-to-S
  //          payload = edge encoding
  // D4 heap: key = dualArr[b]/2 for T-blossom b
  //
  // When a vertex becomes S-labeled at time deltaSum=d, we set
  //   dualArr[v] += d  (so effective y = dualArr[v] - deltaSum)
  // When a vertex becomes T-labeled at time deltaSum=d, we set
  //   dualArr[v] -= d  (so effective y = dualArr[v] + deltaSum)
  // When a blossom b is created (S) at time deltaSum=d, we set
  //   dualArr[b] = -2*d  (effective z = dualArr[b] + 2*deltaSum)
  // When a blossom b is T-labeled at time deltaSum=d, we set
  //   dualArr[b] += 2*d  (effective z = dualArr[b] - 2*deltaSum)
  // ================================================================

  // --- Lazy-deletion min-heaps ---
  // std::greater gives min-heap (smallest key at top).

  // Offset state tracking: records whether each vertex's dualArr has been
  // adjusted for S(+deltaSum) or T(-deltaSum) labeling.
  //   0 = not adjusted (unlabeled or raw)
  //   1 = adjusted for S:  dualArr[v] += labelDelta[v]; effective = dualArr[v] - deltaSum
  //   2 = adjusted for T:  dualArr[v] -= labelDelta[v]; effective = dualArr[v] + deltaSum
  std::vector<Size> vtxOffsetState(numVtxs, 0);
  std::vector<Size> offsetVtxs;
  offsetVtxs.reserve(numVtxs);
  std::vector<unsigned char> inOffsetVtx(numVtxs, 0);
  std::vector<Size> sScanVtxs;
  sScanVtxs.reserve(numVtxs);
  std::vector<unsigned char> inSScan(numVtxs, 0);
  std::vector<Size> tLabelVtxs;
  tLabelVtxs.reserve(numVtxs);
  std::vector<unsigned char> inTLabel(numVtxs, 0);
  std::vector<Size> topLevelSBlsms;
  topLevelSBlsms.reserve(maxBlsm);
  std::vector<unsigned char> inTopLevelSBlsm(maxBlsm, 0);
  std::vector<Size> topLevelTBlsms;
  topLevelTBlsms.reserve(maxBlsm);
  std::vector<unsigned char> inTopLevelTBlsm(maxBlsm, 0);
  std::vector<Size> labeledBlsms;
  labeledBlsms.reserve(maxBlsm);
  std::vector<unsigned char> inLabeledBlsm(maxBlsm, 0);
  std::vector<Size> blsminTouched;
  blsminTouched.reserve(maxBlsm);
  std::vector<unsigned char> inBlsminTouched(maxBlsm, 0);

  Val deltaSum = 0;
  // B1: topLevelBlsmOfVtx is maintained incrementally — initialize to
  // identity (no blossoms yet), and update at shrink/dissolve sites.
  std::vector<Size> topLevelBlsmOfVtx(numVtxs);
  for (Size v = 0; v < numVtxs; ++v) topLevelBlsmOfVtx[v] = v;
  std::vector<Size> freeVtxs;
  freeVtxs.reserve(numVtxs);
  for (Size v = 0; v < numVtxs; ++v)
  {
    if (mateArr[v] == cNullItm)
      freeVtxs.push_back(v);
  }

  // Helper: scan edges of a newly S-labeled vertex and insert into heaps.
  // This is called when vertex `s` becomes S.
  // `s` must be a vertex (< numVtxs), already S-labeled with dualArr adjusted.
  // We scan its adjacency list and insert D2/D3 entries.
  // Note: we also check for immediate augmentation / grow / shrink.

  bool augmented = true;
  Size stageNum = 0;
  std::stack<Size> freeBlsmStk; // for blossom ID recycling

  // B4: persist freeBlsmStk across stages. Blossoms dissolved mid-stage
  // (D4 expand sites) append here instead of being recollected via an
  // O(blsm) rebuild scan at the next stage start.
  std::vector<Size> pendingFreeIds;

  // B5+B6+B7: lift event-queue and scratch containers out of the stage loop.
  // Re-using these across stages avoids O(n) per-stage allocation of the
  // heap backing vectors and the 16n-capacity DFS stack / processed queue.
  // At stage start we only need an O(1) Clear()/clear(), and the Phase-2
  // heaps are already drained to empty by the augmentation branch.
  GeneralizedPQ1 pqD1;
  typedef std::pair<Size, Size> EdgeEventRef;
  typedef std::pair<Val, EdgeEventRef> EdgeEvent;
  std::priority_queue<EdgeEvent, std::vector<EdgeEvent>, std::greater<EdgeEvent> > pqD2;
  std::priority_queue<EdgeEvent, std::vector<EdgeEvent>, std::greater<EdgeEvent> > pqD3;
  std::priority_queue<std::pair<Val, Size>, std::vector<std::pair<Val, Size> >,
      std::greater<std::pair<Val, Size> > > pqD4;
  ItmStk prcsbStk(16 * numVtxs);
  ItmQue prcsdQue(16 * numVtxs);
  std::queue<Size> tempQ;

  while (augmented)
  {
    augmented = false;
    ++stageNum;

    // Check if any free vertex exists.
    Size freeCount = freeVtxs.size();
    if (freeCount == 0) break;

#ifdef EDMONDS_VERBOSE
    if (stageNum <= 20 || stageNum % 500 == 0 || freeCount <= 20)
      std::cerr << "[S" << stageNum << "] free=" << freeCount << " blsm=" << blsm << std::endl;
#endif

    // --- Stage start: expand z=0 blossoms, reset labels ---
    // B3: iterate only the blossoms touched in the previous stage
    // (labeledBlsms — still populated until its reset below) instead of
    // scanning the full [numVtxs, blsm) range. Any blossom that hasn't
    // been labeled since its creation has z unchanged from creation;
    // z=0 blossoms that existed at the previous stage's end must have
    // been among the labeled ones, so this covers every real candidate.
    // Dissolved ids are collected in pendingFreeIds; drained into
    // freeBlsmStk after the trim-tail loop below (B4).
    for (Size li = 0; li < labeledBlsms.size(); ++li)
    {
      Size b = labeledBlsms[li];
      if ((rprsArr[b] != cNullItm) && (blsmParentArr[b] == cNullItm) && (dualArr[b] <= EPS))
      {
        dualArr[b] = 0;
        rprsArr[b] = cNullItm;
        rDissolveBlsm3<ItmQue>(graph, blsmArr, blsmParentArr, blsmChildrenArr, b);
        blsmin[b - numVtxs] = cNullItm;
        // B1: after dissolving top-level b, its direct children become
        // top-level. Refresh topLevelBlsmOfVtx only for vertices inside b.
        {
          const std::vector<Size>& kids = blsmChildrenArr[b - numVtxs];
          for (Size ci = 0; ci < kids.size(); ++ci) {
            Size c = kids[ci];
            if (c < numVtxs) {
              topLevelBlsmOfVtx[c] = c;
            } else {
              tempQ.push(c - numVtxs);
              while (!tempQ.empty()) {
                Size bb = tempQ.front(); tempQ.pop();
                const std::vector<Size>& kk = blsmChildrenArr[bb];
                for (Size j = 0; j < kk.size(); ++j) {
                  Size gc = kk[j];
                  if (gc < numVtxs) {
                    topLevelBlsmOfVtx[gc] = c;
                    blsmArr[gc] = c;
                  }
                  else tempQ.push(gc - numVtxs);
                }
              }
            }
          }
        }
        pendingFreeIds.push_back(b);
      }
    }

    // B2: blsmArr is invalidated incrementally at every dissolve site
    // (we set blsmArr[gc] = c for each nested leaf and rDissolveBlsm3
    // sets blsmArr[c] = c for direct leaves), so no stage-start reset
    // is needed. No vertex retains a recycled blossom id in its cache.

    // Trim trailing dead blossom ids so blsm doesn't grow unbounded.
    while (blsm > numVtxs && rprsArr[blsm - 1] == cNullItm)
    {
      Size dead = blsm - 1;
      Size shiftDead = dead - numVtxs;
      dualArr[dead] = 0;
      blsmParentArr[dead] = cNullItm;
      blsmin[shiftDead] = cNullItm;
      blsmChildrenArr[shiftDead].clear();
      blsmL[shiftDead].clear();
      --blsm;
    }

    // Drain pending-free ids into freeBlsmStk (B4). Filter out ids that
    // were trimmed away above, and normalize per-id state for reuse.
    for (Size i = 0; i < pendingFreeIds.size(); ++i)
    {
      Size b = pendingFreeIds[i];
      if (b >= blsm) continue;              // trimmed away
      if (rprsArr[b] != cNullItm) continue; // re-used already
      Size shiftB = b - numVtxs;
      dualArr[b] = 0;
      blsmParentArr[b] = cNullItm;
      blsmin[shiftB] = cNullItm;
      blsmChildrenArr[shiftB].clear();
      blsmL[shiftB].clear();
      freeBlsmStk.push(b);
    }
    pendingFreeIds.clear();

    for (Size i = 0; i < sScanVtxs.size(); ++i)
    {
      Size v = sScanVtxs[i];
      ptrArr[v] = cNullItm;
      sttArr[v] = eSttIdle;
      idxArr[v] = 0;
      label[v] = 0;
      treeNumArr[v] = cNullItm;
      inSScan[v] = 0;
    }
    sScanVtxs.clear();
    for (Size i = 0; i < tLabelVtxs.size(); ++i)
    {
      Size v = tLabelVtxs[i];
      label[v] = 0;
      treeNumArr[v] = cNullItm;
      inTLabel[v] = 0;
    }
    tLabelVtxs.clear();
    for (Size i = 0; i < offsetVtxs.size(); ++i)
    {
      Size v = offsetVtxs[i];
      vtxOffsetState[v] = 0;
      inOffsetVtx[v] = 0;
    }
    offsetVtxs.clear();
    for (Size i = 0; i < topLevelSBlsms.size(); ++i)
      inTopLevelSBlsm[topLevelSBlsms[i] - numVtxs] = 0;
    topLevelSBlsms.clear();
    for (Size i = 0; i < topLevelTBlsms.size(); ++i)
      inTopLevelTBlsm[topLevelTBlsms[i] - numVtxs] = 0;
    topLevelTBlsms.clear();
    for (Size i = 0; i < labeledBlsms.size(); ++i)
    {
      Size b = labeledBlsms[i];
      label[b] = 0;
      treeNumArr[b] = cNullItm;
      inLabeledBlsm[b - numVtxs] = 0;
    }
    labeledBlsms.clear();
    for (Size i = 0; i < blsminTouched.size(); ++i)
    {
      Size b = blsminTouched[i];
      blsmin[b - numVtxs] = cNullItm;
      inBlsminTouched[b - numVtxs] = 0;
    }
    blsminTouched.clear();

    // B1: topLevelBlsmOfVtx is now maintained incrementally at shrink /
    // dissolve sites. No full rebuild needed at stage start.

    // Clear candidate structures for the current stage.
    // Heaps and scratch containers are declared outside the stage loop
    // (see B5+B6+B7 lift above); here we just reset them.
    deltaSum = 0;
    pqD1.clear();
    while (!pqD2.empty()) pqD2.pop();
    while (!pqD3.empty()) pqD3.pop();
    while (!pqD4.empty()) pqD4.pop();

    // Materialize duals: at stage start deltaSum=0, so dualArr values
    // are already the effective (materialized) values. No adjustment needed.

    prcsbStk.Clear();
    prcsdQue.Clear();
    while (!tempQ.empty()) tempQ.pop();

    // Multi-tree initialization: all free vertices become S-roots.
    // At stage start, deltaSum=0, so dualArr[v] += 0 is a no-op.
    // But we must set vtxOffsetState[v] = 1 so later transitions are correct.
    for (Size fi = 0; fi < freeVtxs.size(); ++fi)
    {
      Size fv = freeVtxs[fi];
      if (mateArr[fv] != cNullItm) continue;
      Size Bfv = topLevelBlsmOfVtx[fv];
      if (label[Bfv] != 0) continue;

      if (Bfv >= numVtxs)
      {
        label[Bfv] = 1;
        treeNumArr[Bfv] = fv;
        if (inLabeledBlsm[Bfv - numVtxs] == 0)
        {
          inLabeledBlsm[Bfv - numVtxs] = 1;
          labeledBlsms.push_back(Bfv);
        }
        if (inTopLevelSBlsm[Bfv - numVtxs] == 0)
        {
          inTopLevelSBlsm[Bfv - numVtxs] = 1;
          topLevelSBlsms.push_back(Bfv);
        }
        pqD1.insert(fv, dualArr[fv]);
        tempQ.push(Bfv - numVtxs);
        while (tempQ.empty() == false)
        {
          Size tempB2 = tempQ.front(); tempQ.pop();
          for (Size ii = 0; ii < blsmChildrenArr[tempB2].size(); ++ii)
          {
            Size child = blsmChildrenArr[tempB2][ii];
            if (child < numVtxs)
            {
              if (sttArr[child] == eSttIdle)
              {
                label[child] = 1; treeNumArr[child] = fv;
                if (inSScan[child] == 0)
                {
                  inSScan[child] = 1;
                  sScanVtxs.push_back(child);
                }
                vtxOffsetState[child] = 1; // S-adjusted (deltaSum=0 so no-op on dualArr)
                if (inOffsetVtx[child] == 0)
                {
                  inOffsetVtx[child] = 1;
                  offsetVtxs.push_back(child);
                }
                sttArr[child] = eSttDfsPrcsb;
                prcsbStk.Push(child);
              }
            }
            else
            {
              label[child] = 1; treeNumArr[child] = fv;
              if (inLabeledBlsm[child - numVtxs] == 0)
              {
                inLabeledBlsm[child - numVtxs] = 1;
                labeledBlsms.push_back(child);
              }
              tempQ.push(child - numVtxs);
            }
          }
        }
      }
      else
      {
        label[fv] = 1; treeNumArr[fv] = fv;
        if (inSScan[fv] == 0)
        {
          inSScan[fv] = 1;
          sScanVtxs.push_back(fv);
        }
        vtxOffsetState[fv] = 1; // S-adjusted
        if (inOffsetVtx[fv] == 0)
        {
          inOffsetVtx[fv] = 1;
          offsetVtxs.push_back(fv);
        }
        sttArr[fv] = eSttDfsPrcsb;
        prcsbStk.Push(fv);
        pqD1.insert(fv, dualArr[fv]);
      }
    }

    // Initial S-vertices are already on prcsbStk.
    // Phase 1 of the event loop will scan them, processing tight edges
    // inline and pushing non-tight edges to heaps. No separate initial
    // scan needed.


    Size sLast = cNullItm;
    Size tLast = cNullItm;
    bool foundAugPath = false;
    bool stageDone = false;

    // ================================================================
    // Event loop with two phases per iteration:
    //   Phase 1: DFS scan on tight edges (grow/shrink/augment inline)
    //   Phase 2: If DFS exhausted, advance deltaSum to next heap event
    //
    // This avoids pushing tight edges to heaps (which causes O(nm) heap
    // entries and bad_alloc on large graphs).
    // ================================================================
#ifdef EDMONDS_VERBOSE
    Size innerIter = 0;
#endif
    while (!foundAugPath && !stageDone)
    {
#ifdef EDMONDS_VERBOSE
      ++innerIter;
      if (innerIter > 200 * numVtxs) {
        std::cerr << "[STUCK] stage=" << stageNum << " innerIter=" << innerIter
                  << " deltaSum=" << deltaSum << std::endl;
        stageDone = true;
        break;
      }
#endif
      // ---- Phase 1: DFS on tight edges ----
      // Process all tight edges reachable from S-vertices on prcsbStk.
      // This is essentially the same DFS as the old rebuild-based code,
      // except we also push non-tight edges to heaps.
      while ((prcsbStk.Empty() == false) && (foundAugPath == false))
      {
        Size sv = prcsbStk.Top();
        assert(sv < numVtxs);

        Size blsmsv = topLevelBlsmOfVtx[sv];
        assert(blsmsv < numVtxs + maxBlsm);

        if (label[blsmsv] != 1) {
          prcsbStk.Pop(); sttArr[sv] = eSttDfsPrcsd;
          continue;
        }

        // Adjust dual offset for newly S-labeled vertex (from shrink T→S)
        if (vtxOffsetState[sv] == 2)
        {
          dualArr[sv] += 2.0 * deltaSum;
          vtxOffsetState[sv] = 1;
        }
        else if (vtxOffsetState[sv] == 0)
        {
          dualArr[sv] += deltaSum;
          vtxOffsetState[sv] = 1;
          if (inOffsetVtx[sv] == 0)
          {
            inOffsetVtx[sv] = 1;
            offsetVtxs.push_back(sv);
          }
        }

        Size svTreeNum = treeNumArr[blsmsv];
        bool advanced = false;
        const std::vector<Size>& svNbrVec = vtxVecArr[sv];
        const std::vector<Val>& svWVec = edgWghtVecArr[sv];
        Size& svIdx = idxArr[sv];

        for ( ; svIdx < svNbrVec.size(); ++svIdx)
        {
          Size t2 = svNbrVec[svIdx];
          if (t2 == mateArr[sv]) continue;

          Size blsmt2 = topLevelBlsmOfVtx[t2];
          if (blsmsv == blsmt2) continue;

          Val rawSlack = dualArr[sv] + dualArr[t2] - svWVec[svIdx];
          if (rawSlack < -EPS) rawSlack = 0;

          if (label[blsmt2] == 0)
          {
            // S-to-unlabeled: effective slack = rawSlack - deltaSum
            if (rawSlack <= deltaSum + EPS)
            {
              // TIGHT: grow immediately
              if (blsmt2 >= numVtxs)
              {
                assert(blsmt2 - numVtxs < maxBlsm);
                assert(rprsArr[blsmt2] != cNullItm);
                blsmin[blsmt2 - numVtxs] = t2;
                if (inBlsminTouched[blsmt2 - numVtxs] == 0)
                {
                  inBlsminTouched[blsmt2 - numVtxs] = 1;
                  blsminTouched.push_back(blsmt2);
                }
                label[blsmt2] = 2; treeNumArr[blsmt2] = svTreeNum;
                if (inLabeledBlsm[blsmt2 - numVtxs] == 0)
                {
                  inLabeledBlsm[blsmt2 - numVtxs] = 1;
                  labeledBlsms.push_back(blsmt2);
                }
                if (inTopLevelTBlsm[blsmt2 - numVtxs] == 0)
                {
                  inTopLevelTBlsm[blsmt2 - numVtxs] = 1;
                  topLevelTBlsms.push_back(blsmt2);
                }
                tempQ.push(blsmt2 - numVtxs);
                while (tempQ.empty() == false)
                {
                  Size tempB3 = tempQ.front(); tempQ.pop();
                  for (Size iii = 0; iii < blsmChildrenArr[tempB3].size(); ++iii)
                  {
                    Size child = blsmChildrenArr[tempB3][iii];
                    if (child < numVtxs)
                    {
                      if (inTLabel[child] == 0)
                      {
                        inTLabel[child] = 1;
                        tLabelVtxs.push_back(child);
                      }
                    }
                    else
                    {
                      tempQ.push(child - numVtxs);
                    }
                  }
                }
                dualArr[blsmt2] += 2 * deltaSum;
                label[t2] = 2; treeNumArr[t2] = svTreeNum;
                dualArr[t2] -= deltaSum; vtxOffsetState[t2] = 2;
                if (inOffsetVtx[t2] == 0)
                {
                  inOffsetVtx[t2] = 1;
                  offsetVtxs.push_back(t2);
                }

                assert(rprsArr[blsmt2] < numVtxs);
                Size ss2 = mateArr[rprsArr[blsmt2]];
                if (ss2 == cNullItm)
                {
                  sLast = sv; tLast = t2;
                  foundAugPath = true; break;
                }
                ptrArr[ss2] = sv;
                Size blsmss2 = topLevelBlsmOfVtx[ss2];
                label[blsmss2] = 1; treeNumArr[blsmss2] = svTreeNum;
                if (blsmss2 >= numVtxs)
                {
                  if (inLabeledBlsm[blsmss2 - numVtxs] == 0)
                  {
                    inLabeledBlsm[blsmss2 - numVtxs] = 1;
                    labeledBlsms.push_back(blsmss2);
                  }
                  if (inTopLevelSBlsm[blsmss2 - numVtxs] == 0)
                  {
                    inTopLevelSBlsm[blsmss2 - numVtxs] = 1;
                    topLevelSBlsms.push_back(blsmss2);
                  }
                  tempQ.push(blsmss2 - numVtxs);
                  while (tempQ.empty() == false)
                  {
                    Size tempB = tempQ.front(); tempQ.pop();
                    for (Size i = 0; i < blsmChildrenArr[tempB].size(); ++i)
                    {
                      Size child = blsmChildrenArr[tempB][i];
                      if (child < numVtxs) {
                        if (sttArr[child] == eSttIdle) {
                          label[child] = 1; treeNumArr[child] = svTreeNum;
                          if (inSScan[child] == 0)
                          {
                            inSScan[child] = 1;
                            sScanVtxs.push_back(child);
                          }
                          sttArr[child] = eSttDfsPrcsb; prcsbStk.Push(child);
                        }
                      } else {
                        label[child] = 1; treeNumArr[child] = svTreeNum;
                        if (inLabeledBlsm[child - numVtxs] == 0)
                        {
                          inLabeledBlsm[child - numVtxs] = 1;
                          labeledBlsms.push_back(child);
                        }
                        tempQ.push(child - numVtxs);
                      }
                    }
                  }
                }
                else
                {
                  if (sttArr[ss2] == eSttIdle) {
                    label[ss2] = 1; treeNumArr[ss2] = svTreeNum;
                    if (inSScan[ss2] == 0)
                    {
                      inSScan[ss2] = 1;
                      sScanVtxs.push_back(ss2);
                    }
                    sttArr[ss2] = eSttDfsPrcsb; prcsbStk.Push(ss2);
                  }
                }
              }
              else
              {
                // t2 is a plain unlabeled vertex
                label[t2] = 2; treeNumArr[t2] = svTreeNum;
                dualArr[t2] -= deltaSum; vtxOffsetState[t2] = 2;
                if (inOffsetVtx[t2] == 0)
                {
                  inOffsetVtx[t2] = 1;
                  offsetVtxs.push_back(t2);
                }
                if (inTLabel[t2] == 0)
                {
                  inTLabel[t2] = 1;
                  tLabelVtxs.push_back(t2);
                }

                Size ss2 = mateArr[t2];
                if (ss2 == cNullItm)
                {
                  sLast = sv; tLast = t2;
                  foundAugPath = true; break;
                }
                ptrArr[ss2] = sv;
                Size blsmss2 = topLevelBlsmOfVtx[ss2];
                label[blsmss2] = 1; treeNumArr[blsmss2] = svTreeNum;
                if (blsmss2 >= numVtxs)
                {
                  if (inLabeledBlsm[blsmss2 - numVtxs] == 0)
                  {
                    inLabeledBlsm[blsmss2 - numVtxs] = 1;
                    labeledBlsms.push_back(blsmss2);
                  }
                  if (inTopLevelSBlsm[blsmss2 - numVtxs] == 0)
                  {
                    inTopLevelSBlsm[blsmss2 - numVtxs] = 1;
                    topLevelSBlsms.push_back(blsmss2);
                  }
                  tempQ.push(blsmss2 - numVtxs);
                  while (tempQ.empty() == false)
                  {
                    Size tempB = tempQ.front(); tempQ.pop();
                    for (Size i = 0; i < blsmChildrenArr[tempB].size(); ++i)
                    {
                      Size child = blsmChildrenArr[tempB][i];
                      if (child < numVtxs) {
                        if (sttArr[child] == eSttIdle) {
                          label[child] = 1; treeNumArr[child] = svTreeNum;
                          if (inSScan[child] == 0)
                          {
                            inSScan[child] = 1;
                            sScanVtxs.push_back(child);
                          }
                          sttArr[child] = eSttDfsPrcsb; prcsbStk.Push(child);
                        }
                      } else {
                        label[child] = 1; treeNumArr[child] = svTreeNum;
                        if (inLabeledBlsm[child - numVtxs] == 0)
                        {
                          inLabeledBlsm[child - numVtxs] = 1;
                          labeledBlsms.push_back(child);
                        }
                        tempQ.push(child - numVtxs);
                      }
                    }
                  }
                }
                else
                {
                  if (sttArr[ss2] == eSttIdle) {
                    label[ss2] = 1; treeNumArr[ss2] = svTreeNum;
                    if (inSScan[ss2] == 0)
                    {
                      inSScan[ss2] = 1;
                      sScanVtxs.push_back(ss2);
                    }
                    sttArr[ss2] = eSttDfsPrcsb; prcsbStk.Push(ss2);
                  }
                }
              }
              advanced = true;
              break; // restart DFS from stack top
            }
            else
            {
              // Experimental D2 event queues are disabled on the active exact path.
            }
          }
          else if (label[blsmt2] == 1)
          {
            // S-to-S: effective slack = rawSlack - 2*deltaSum
            if (rawSlack <= 2 * deltaSum + EPS)
            {
              // TIGHT: shrink or augment immediately
              Size t2TreeNum = treeNumArr[blsmt2];
              if (svTreeNum != cNullItm && t2TreeNum != cNullItm && svTreeNum != t2TreeNum)
              {
                sLast = sv; tLast = t2;
                foundAugPath = true; break;
              }
              else
              {
                Size newBlsmId = cNullItm;
                rProcessBlsmEdmonds<ItmQue, ItmStk>(graph,
                    mateArr, ptrArr, blsmArr, rprsArr,
                    blsmParentArr, blsmChildrenArr, blsmL,
                    label, sttArr, prcsbStk,
                    t2, sv, &blsm, dualArr, blsmin, treeNumArr,
                    maxBlsm, freeBlsmStk, &newBlsmId);
                if (newBlsmId != cNullItm)
                {
                  if (inLabeledBlsm[newBlsmId - numVtxs] == 0)
                  {
                    inLabeledBlsm[newBlsmId - numVtxs] = 1;
                    labeledBlsms.push_back(newBlsmId);
                  }
                  if (inTopLevelSBlsm[newBlsmId - numVtxs] == 0)
                  {
                    inTopLevelSBlsm[newBlsmId - numVtxs] = 1;
                    topLevelSBlsms.push_back(newBlsmId);
                  }
                  std::queue<Size> scanQ;
                  scanQ.push(newBlsmId - numVtxs);
                  while (scanQ.empty() == false)
                  {
                    Size tempB = scanQ.front(); scanQ.pop();
                    for (Size ci = 0; ci < blsmChildrenArr[tempB].size(); ++ci)
                    {
                      Size child = blsmChildrenArr[tempB][ci];
                      if (child < numVtxs)
                      {
                        topLevelBlsmOfVtx[child] = newBlsmId;
                        if (inSScan[child] == 0)
                        {
                          inSScan[child] = 1;
                          sScanVtxs.push_back(child);
                        }
                      }
                      else
                      {
                        scanQ.push(child - numVtxs);
                      }
                    }
                  }
                  dualArr[newBlsmId] = -2.0 * deltaSum;
                }
                advanced = true;
                break; // restart DFS — shrink pushed new S-vertices
              }
            }
            else
            {
              // Experimental D3 event queues are disabled on the active exact path.
            }
          }
          // T-labeled: skip
        }

        if (foundAugPath) break;

        if (!advanced)
        {
          prcsbStk.Pop();
          sttArr[sv] = eSttDfsPrcsd;
        }
      } // end Phase 1 DFS

      if (foundAugPath) break;

      bool rebuildForest = false;
      while (false && !foundAugPath && !stageDone && !rebuildForest)
      {
        Val d1 = INF, d2 = INF, d3 = INF, d4 = INF;

        while (!pqD1.empty())
        {
          std::pair<Size, Val> d1Min = pqD1.findMin();
          Size u = d1Min.first;
          Size Bu = u;
          if (blsmParentArr[u] != cNullItm)
            Bu = rFindblsm<ItmQue>(blsmParentArr, blsmArr, rprsArr, u);
          if (mateArr[u] == cNullItm && label[Bu] == 1)
          {
            d1 = d1Min.second;
            if (d1 < deltaSum && d1 > deltaSum - EPS)
              d1 = deltaSum;
            break;
          }
          pqD1.deleteMin();
        }

        while (!pqD2.empty())
        {
          EdgeEvent d2Top = pqD2.top();
          Size u = d2Top.second.first;
          Size edgeIdx = d2Top.second.second;
          if (edgeIdx >= vtxVecArr[u].size())
          {
            pqD2.pop();
            continue;
          }
          Size v = vtxVecArr[u][edgeIdx];
          if (v == mateArr[u])
          {
            pqD2.pop();
            continue;
          }
          Size Bu = u;
          if (blsmParentArr[u] != cNullItm)
            Bu = rFindblsm<ItmQue>(blsmParentArr, blsmArr, rprsArr, u);
          if (label[Bu] != 1)
          {
            pqD2.pop();
            continue;
          }
          Size Bv = v;
          if (blsmParentArr[v] != cNullItm)
            Bv = rFindblsm<ItmQue>(blsmParentArr, blsmArr, rprsArr, v);
          if (Bu == Bv || label[Bv] != 0)
          {
            pqD2.pop();
            continue;
          }
          Val rawSlack = dualArr[u] + dualArr[v] - edgWghtVecArr[u][edgeIdx];
          if (rawSlack < -EPS) rawSlack = 0;
          if (rawSlack < d2Top.first - EPS || rawSlack > d2Top.first + EPS)
          {
            pqD2.pop();
            continue;
          }
          d2 = d2Top.first;
          if (d2 < deltaSum && d2 > deltaSum - EPS)
            d2 = deltaSum;
          break;
        }

        while (!pqD3.empty())
        {
          EdgeEvent d3Top = pqD3.top();
          Size u = d3Top.second.first;
          Size edgeIdx = d3Top.second.second;
          if (edgeIdx >= vtxVecArr[u].size())
          {
            pqD3.pop();
            continue;
          }
          Size v = vtxVecArr[u][edgeIdx];
          if (v == mateArr[u])
          {
            pqD3.pop();
            continue;
          }
          Size Bu = u;
          if (blsmParentArr[u] != cNullItm)
            Bu = rFindblsm<ItmQue>(blsmParentArr, blsmArr, rprsArr, u);
          if (label[Bu] != 1)
          {
            pqD3.pop();
            continue;
          }
          Size Bv = v;
          if (blsmParentArr[v] != cNullItm)
            Bv = rFindblsm<ItmQue>(blsmParentArr, blsmArr, rprsArr, v);
          if (Bu == Bv || label[Bv] != 1)
          {
            pqD3.pop();
            continue;
          }
          Val rawSlack = dualArr[u] + dualArr[v] - edgWghtVecArr[u][edgeIdx];
          if (rawSlack < -EPS) rawSlack = 0;
          Val halfSlack = rawSlack / 2.0;
          if (halfSlack < d3Top.first - EPS || halfSlack > d3Top.first + EPS)
          {
            pqD3.pop();
            continue;
          }
          d3 = d3Top.first;
          if (d3 < deltaSum && d3 > deltaSum - EPS)
            d3 = deltaSum;
          break;
        }

        while (!pqD4.empty())
        {
          std::pair<Val, Size> d4Top = pqD4.top();
          Size b = d4Top.second;
          if (b < numVtxs || rprsArr[b] == cNullItm || blsmParentArr[b] != cNullItm ||
              label[b] != 2)
          {
            pqD4.pop();
            continue;
          }
          Val cand = dualArr[b] / 2.0;
          if (cand < d4Top.first - EPS || cand > d4Top.first + EPS)
          {
            pqD4.pop();
            continue;
          }
          d4 = d4Top.first;
          if (d4 < deltaSum && d4 > deltaSum - EPS)
            d4 = deltaSum;
          break;
        }

        Val nextDelta = d1;
        if (d2 < nextDelta) nextDelta = d2;
        if (d3 < nextDelta) nextDelta = d3;
        if (d4 < nextDelta) nextDelta = d4;

        if (nextDelta >= INF || nextDelta < -EPS)
        {
          stageDone = true;
          break;
        }

        if (nextDelta > deltaSum)
          deltaSum = nextDelta;

        bool processedEvent = false;

        while (!pqD2.empty())
        {
          EdgeEvent d2Top = pqD2.top();
          if (d2Top.first > deltaSum + EPS)
            break;
          pqD2.pop();

          Size suEv = d2Top.second.first;
          Size edgeIdx = d2Top.second.second;
          if (edgeIdx >= vtxVecArr[suEv].size())
            continue;

          Size t2 = vtxVecArr[suEv][edgeIdx];
          if (t2 == mateArr[suEv])
            continue;

          Size Bu = suEv;
          if (blsmParentArr[suEv] != cNullItm)
            Bu = rFindblsm<ItmQue>(blsmParentArr, blsmArr, rprsArr, suEv);
          if (label[Bu] != 1)
            continue;

          Size blsmt2 = t2;
          if (blsmParentArr[t2] != cNullItm)
            blsmt2 = rFindblsm<ItmQue>(blsmParentArr, blsmArr, rprsArr, t2);
          if (Bu == blsmt2 || label[blsmt2] != 0)
            continue;

          Val rawSlack = dualArr[suEv] + dualArr[t2] - edgWghtVecArr[suEv][edgeIdx];
          if (rawSlack < -EPS) rawSlack = 0;
          if (rawSlack > deltaSum + EPS)
            continue;

          Size suTreeNum = treeNumArr[Bu];
          if (blsmt2 >= numVtxs)
          {
            assert(blsmt2 - numVtxs < maxBlsm);
            assert(rprsArr[blsmt2] != cNullItm);
            blsmin[blsmt2 - numVtxs] = t2;
            label[blsmt2] = 2; treeNumArr[blsmt2] = suTreeNum;
            if (inTopLevelTBlsm[blsmt2 - numVtxs] == 0)
            {
              inTopLevelTBlsm[blsmt2 - numVtxs] = 1;
              topLevelTBlsms.push_back(blsmt2);
            }
            tempQ.push(blsmt2 - numVtxs);
            while (tempQ.empty() == false)
            {
              Size tempB3 = tempQ.front(); tempQ.pop();
              for (Size iii = 0; iii < blsmChildrenArr[tempB3].size(); ++iii)
              {
                Size child = blsmChildrenArr[tempB3][iii];
                if (child < numVtxs)
                {
                  if (inTLabel[child] == 0)
                  {
                    inTLabel[child] = 1;
                    tLabelVtxs.push_back(child);
                  }
                }
                else
                {
                  tempQ.push(child - numVtxs);
                }
              }
            }
            dualArr[blsmt2] += 2 * deltaSum;
            pqD4.push(std::make_pair(dualArr[blsmt2] / 2.0, blsmt2));
            label[t2] = 2; treeNumArr[t2] = suTreeNum;
            dualArr[t2] -= deltaSum; vtxOffsetState[t2] = 2;
            if (inOffsetVtx[t2] == 0)
            {
              inOffsetVtx[t2] = 1;
              offsetVtxs.push_back(t2);
            }

            assert(rprsArr[blsmt2] < numVtxs);
            Size ss2 = mateArr[rprsArr[blsmt2]];
            if (ss2 == cNullItm)
            {
              sLast = suEv; tLast = t2;
              foundAugPath = true; break;
            }
            ptrArr[ss2] = suEv;
            Size blsmss2 = ss2;
            if (blsmParentArr[ss2] != cNullItm)
              blsmss2 = rFindblsm<ItmQue>(blsmParentArr, blsmArr, rprsArr, ss2);
            label[blsmss2] = 1; treeNumArr[blsmss2] = suTreeNum;
            if (blsmss2 >= numVtxs)
            {
              if (inTopLevelSBlsm[blsmss2 - numVtxs] == 0)
              {
                inTopLevelSBlsm[blsmss2 - numVtxs] = 1;
                topLevelSBlsms.push_back(blsmss2);
              }
              tempQ.push(blsmss2 - numVtxs);
              while (tempQ.empty() == false)
              {
                Size tempB = tempQ.front(); tempQ.pop();
                for (Size i = 0; i < blsmChildrenArr[tempB].size(); ++i)
                {
                  Size child = blsmChildrenArr[tempB][i];
                  if (child < numVtxs) {
                    if (sttArr[child] == eSttIdle) {
                      label[child] = 1; treeNumArr[child] = suTreeNum;
                      if (inSScan[child] == 0)
                      {
                        inSScan[child] = 1;
                        sScanVtxs.push_back(child);
                      }
                      sttArr[child] = eSttDfsPrcsb; prcsbStk.Push(child);
                    }
                  } else {
                    label[child] = 1; treeNumArr[child] = suTreeNum;
                    tempQ.push(child - numVtxs);
                  }
                }
              }
            }
            else
            {
              if (sttArr[ss2] == eSttIdle) {
                label[ss2] = 1; treeNumArr[ss2] = suTreeNum;
                if (inSScan[ss2] == 0)
                {
                  inSScan[ss2] = 1;
                  sScanVtxs.push_back(ss2);
                }
                sttArr[ss2] = eSttDfsPrcsb; prcsbStk.Push(ss2);
              }
            }
          }
          else
          {
            label[t2] = 2; treeNumArr[t2] = suTreeNum;
            dualArr[t2] -= deltaSum; vtxOffsetState[t2] = 2;
            if (inOffsetVtx[t2] == 0)
            {
              inOffsetVtx[t2] = 1;
              offsetVtxs.push_back(t2);
            }
            if (inTLabel[t2] == 0)
            {
              inTLabel[t2] = 1;
              tLabelVtxs.push_back(t2);
            }

            Size ss2 = mateArr[t2];
            if (ss2 == cNullItm)
            {
              sLast = suEv; tLast = t2;
              foundAugPath = true; break;
            }
            ptrArr[ss2] = suEv;
            Size blsmss2 = ss2;
            if (blsmParentArr[ss2] != cNullItm)
              blsmss2 = rFindblsm<ItmQue>(blsmParentArr, blsmArr, rprsArr, ss2);
            label[blsmss2] = 1; treeNumArr[blsmss2] = suTreeNum;
            if (blsmss2 >= numVtxs)
            {
              if (inTopLevelSBlsm[blsmss2 - numVtxs] == 0)
              {
                inTopLevelSBlsm[blsmss2 - numVtxs] = 1;
                topLevelSBlsms.push_back(blsmss2);
              }
              tempQ.push(blsmss2 - numVtxs);
              while (tempQ.empty() == false)
              {
                Size tempB = tempQ.front(); tempQ.pop();
                for (Size i = 0; i < blsmChildrenArr[tempB].size(); ++i)
                {
                  Size child = blsmChildrenArr[tempB][i];
                  if (child < numVtxs) {
                    if (sttArr[child] == eSttIdle) {
                      label[child] = 1; treeNumArr[child] = suTreeNum;
                      if (inSScan[child] == 0)
                      {
                        inSScan[child] = 1;
                        sScanVtxs.push_back(child);
                      }
                      sttArr[child] = eSttDfsPrcsb; prcsbStk.Push(child);
                    }
                  } else {
                    label[child] = 1; treeNumArr[child] = suTreeNum;
                    tempQ.push(child - numVtxs);
                  }
                }
              }
            }
            else
            {
              if (sttArr[ss2] == eSttIdle) {
                label[ss2] = 1; treeNumArr[ss2] = suTreeNum;
                if (inSScan[ss2] == 0)
                {
                  inSScan[ss2] = 1;
                  sScanVtxs.push_back(ss2);
                }
                sttArr[ss2] = eSttDfsPrcsb; prcsbStk.Push(ss2);
              }
            }
          }

          processedEvent = true;
          break;
        }

        if (foundAugPath) break;

        while (!processedEvent && !pqD3.empty())
        {
          EdgeEvent d3Top = pqD3.top();
          if (d3Top.first > deltaSum + EPS)
            break;
          pqD3.pop();

          Size suEv = d3Top.second.first;
          Size edgeIdx = d3Top.second.second;
          if (edgeIdx >= vtxVecArr[suEv].size())
            continue;

          Size t2 = vtxVecArr[suEv][edgeIdx];
          if (t2 == mateArr[suEv])
            continue;

          Size Bu = suEv;
          if (blsmParentArr[suEv] != cNullItm)
            Bu = rFindblsm<ItmQue>(blsmParentArr, blsmArr, rprsArr, suEv);
          if (label[Bu] != 1)
            continue;

          Size blsmt2 = t2;
          if (blsmParentArr[t2] != cNullItm)
            blsmt2 = rFindblsm<ItmQue>(blsmParentArr, blsmArr, rprsArr, t2);
          if (Bu == blsmt2 || label[blsmt2] != 1)
            continue;

          Val rawSlack = dualArr[suEv] + dualArr[t2] - edgWghtVecArr[suEv][edgeIdx];
          if (rawSlack < -EPS) rawSlack = 0;
          if (rawSlack > 2.0 * deltaSum + EPS)
            continue;

          Size suTreeNum = treeNumArr[Bu];
          Size t2TreeNum = treeNumArr[blsmt2];
          if (suTreeNum != cNullItm && t2TreeNum != cNullItm && suTreeNum != t2TreeNum)
          {
            sLast = suEv; tLast = t2;
            foundAugPath = true; break;
          }

          Size newBlsmId = cNullItm;
          rProcessBlsmEdmonds<ItmQue, ItmStk>(graph,
              mateArr, ptrArr, blsmArr, rprsArr,
              blsmParentArr, blsmChildrenArr, blsmL,
              label, sttArr, prcsbStk,
              t2, suEv, &blsm, dualArr, blsmin, treeNumArr,
              maxBlsm, freeBlsmStk, &newBlsmId);
          if (newBlsmId != cNullItm)
          {
            if (inTopLevelSBlsm[newBlsmId - numVtxs] == 0)
            {
              inTopLevelSBlsm[newBlsmId - numVtxs] = 1;
              topLevelSBlsms.push_back(newBlsmId);
            }
            dualArr[newBlsmId] = -2.0 * deltaSum;
          }
          processedEvent = true;
          break;
        }

        if (foundAugPath) break;

        while (!processedEvent && !pqD4.empty())
        {
          std::pair<Val, Size> d4Top = pqD4.top();
          if (d4Top.first > deltaSum + EPS)
            break;
          pqD4.pop();

          Size b = d4Top.second;
          if (b < numVtxs || rprsArr[b] == cNullItm || blsmParentArr[b] != cNullItm ||
              label[b] != 2)
            continue;

          Val cand = dualArr[b] / 2.0;
          if (cand > deltaSum + EPS)
            continue;

          for (Size oi = 0; oi < offsetVtxs.size(); ++oi)
          {
            Size v = offsetVtxs[oi];
            if (vtxOffsetState[v] == 1)
              dualArr[v] -= deltaSum;
            else if (vtxOffsetState[v] == 2)
              dualArr[v] += deltaSum;
            vtxOffsetState[v] = 0;
            inOffsetVtx[v] = 0;
            if (dualArr[v] < 0 && dualArr[v] > -EPS) dualArr[v] = 0;
          }
          offsetVtxs.clear();
          for (Size sb = 0; sb < topLevelSBlsms.size(); ++sb)
          {
            Size sbId = topLevelSBlsms[sb];
            if (rprsArr[sbId] == cNullItm || blsmParentArr[sbId] != cNullItm ||
                label[sbId] != 1)
              continue;
            dualArr[sbId] += 2.0 * deltaSum;
            if (dualArr[sbId] < 0 && dualArr[sbId] > -EPS) dualArr[sbId] = 0;
          }
          for (Size tb = 0; tb < topLevelTBlsms.size(); ++tb)
          {
            Size tbId = topLevelTBlsms[tb];
            if (rprsArr[tbId] == cNullItm || blsmParentArr[tbId] != cNullItm ||
                label[tbId] != 2)
              continue;
            dualArr[tbId] -= 2.0 * deltaSum;
            if (dualArr[tbId] < 0 && dualArr[tbId] > -EPS) dualArr[tbId] = 0;
          }
          deltaSum = 0;

          bool expandedAny = false;
          for (Size tb = 0; tb < topLevelTBlsms.size(); ++tb)
          {
            Size tbId = topLevelTBlsms[tb];
            if ((rprsArr[tbId] != cNullItm) && (blsmParentArr[tbId] == cNullItm) &&
                (label[tbId] == 2) && (dualArr[tbId] <= EPS))
            {
              dualArr[tbId] = 0;
              label[tbId] = 0; treeNumArr[tbId] = cNullItm;
              Size shiftBd = tbId - numVtxs;
              for (Size ci = 0; ci < blsmChildrenArr[shiftBd].size(); ++ci) {
                Size child = blsmChildrenArr[shiftBd][ci];
                label[child] = 0; treeNumArr[child] = cNullItm;
                if (child < numVtxs) sttArr[child] = eSttIdle;
              }
              rprsArr[tbId] = cNullItm;
              rDissolveBlsm3<ItmQue>(graph, blsmArr, blsmParentArr, blsmChildrenArr, tbId);
              blsmin[tbId - numVtxs] = cNullItm;
              // B1: refresh topLevelBlsmOfVtx for vertices that were in tbId.
              {
                const std::vector<Size>& kids = blsmChildrenArr[tbId - numVtxs];
                for (Size ci = 0; ci < kids.size(); ++ci) {
                  Size c = kids[ci];
                  if (c < numVtxs) {
                    topLevelBlsmOfVtx[c] = c;
                  } else {
                    tempQ.push(c - numVtxs);
                    while (!tempQ.empty()) {
                      Size bb = tempQ.front(); tempQ.pop();
                      const std::vector<Size>& kk = blsmChildrenArr[bb];
                      for (Size j = 0; j < kk.size(); ++j) {
                        Size gc = kk[j];
                        if (gc < numVtxs) {
                          topLevelBlsmOfVtx[gc] = c;
                          blsmArr[gc] = c;
                        }
                        else tempQ.push(gc - numVtxs);
                      }
                    }
                  }
                }
              }
              pendingFreeIds.push_back(tbId);
              expandedAny = true;
            }
          }
          if (expandedAny)
            rebuildForest = true;
          processedEvent = true;
          break;
        }

        if (foundAugPath || stageDone || rebuildForest)
          break;
        if (processedEvent)
          break;

        if (d1 <= deltaSum + EPS)
        {
          stageDone = true;
          break;
        }
      }

      if (foundAugPath || stageDone)
        break;
      if (false && !rebuildForest)
        continue;

      if (true)
      {
      // ---- Phase 2: Materialize duals, delta scan, apply, re-push ----
      // Materialize all duals to raw form, then use the old-style delta scan.
      const int offsetCount = static_cast<int>(offsetVtxs.size());
      #pragma omp parallel for if(offsetCount > 1024) schedule(static)
      for (int oi = 0; oi < offsetCount; ++oi)
      {
        Size v = offsetVtxs[static_cast<Size>(oi)];
        if (vtxOffsetState[v] == 1)
          dualArr[v] -= deltaSum;
        else if (vtxOffsetState[v] == 2)
          dualArr[v] += deltaSum;
        vtxOffsetState[v] = 0;
        inOffsetVtx[v] = 0;
        if (dualArr[v] < 0 && dualArr[v] > -EPS) dualArr[v] = 0;
      }
      offsetVtxs.clear();
      const int topLevelSCount = static_cast<int>(topLevelSBlsms.size());
      #pragma omp parallel for if(topLevelSCount > 256) schedule(static)
      for (int sb = 0; sb < topLevelSCount; ++sb)
      {
        Size b = topLevelSBlsms[static_cast<Size>(sb)];
        if (rprsArr[b] == cNullItm || blsmParentArr[b] != cNullItm || label[b] != 1)
          continue;
        dualArr[b] += 2.0 * deltaSum;
        if (dualArr[b] < 0 && dualArr[b] > -EPS) dualArr[b] = 0;
      }
      const int topLevelTCount = static_cast<int>(topLevelTBlsms.size());
      #pragma omp parallel for if(topLevelTCount > 256) schedule(static)
      for (int tb = 0; tb < topLevelTCount; ++tb)
      {
        Size b = topLevelTBlsms[static_cast<Size>(tb)];
        if (rprsArr[b] == cNullItm || blsmParentArr[b] != cNullItm || label[b] != 2)
          continue;
        dualArr[b] -= 2.0 * deltaSum;
        if (dualArr[b] < 0 && dualArr[b] > -EPS) dualArr[b] = 0;
      }
      deltaSum = 0;

      // Delta scan with raw (materialized) duals — same as old code
      {
        Val d1 = INF, d2 = INF, d3 = INF, d4 = INF;
        const int sScanCount = static_cast<int>(sScanVtxs.size());
        if (!pqD1.empty())
          d1 = pqD1.findMin().second;

        {
          const int topLevelTCount2 = static_cast<int>(topLevelTBlsms.size());
          #pragma omp parallel if(sScanCount > 256 || topLevelTCount2 > 256)
          {
            Val localD2 = INF;
            Val localD3 = INF;
            Val localD4 = INF;

            #pragma omp for schedule(static) nowait
            for (int su = 0; su < sScanCount; ++su)
            {
              Size u = sScanVtxs[static_cast<Size>(su)];
              Size Bu = topLevelBlsmOfVtx[u];
              if (label[Bu] != 1) continue;

              const std::vector<Size>& uNbrVec = vtxVecArr[u];
              const std::vector<Val>& uWVec = edgWghtVecArr[u];
              for (Size i = 0; i < uNbrVec.size(); ++i)
              {
                Size v = uNbrVec[i];
                if (v == mateArr[u]) continue;

                Size Bv = topLevelBlsmOfVtx[v];
                if (Bu == Bv) continue;

                Val slack = dualArr[u] + dualArr[v] - uWVec[i];
                if (slack < -EPS) slack = 0;

                if (label[Bv] == 0)
                {
                  if (slack < localD2) localD2 = slack;
                }
                else if (label[Bv] == 1)
                {
                  Val half = slack / 2.0;
                  if (half < localD3) localD3 = half;
                }
              }
            }

            #pragma omp for schedule(static) nowait
            for (int tb = 0; tb < topLevelTCount2; ++tb)
            {
              Size b = topLevelTBlsms[static_cast<Size>(tb)];
              if ((rprsArr[b] != cNullItm) && (blsmParentArr[b] == cNullItm) && (label[b] == 2))
              {
                Val cand = dualArr[b] / 2.0;
                if (cand < localD4) localD4 = cand;
              }
            }

            #pragma omp critical
            {
              if (localD2 < d2) d2 = localD2;
              if (localD3 < d3) d3 = localD3;
              if (localD4 < d4) d4 = localD4;
            }
          }
        }

        Val Delta = d1;
        if (d2 < Delta) Delta = d2;
        if (d3 < Delta) Delta = d3;
        if (d4 < Delta) Delta = d4;

        if (Delta >= INF || Delta < -EPS)
        {
          stageDone = true;
          break;
        }

        if (Delta <= EPS)
        {
          // Delta=0: expand T-blossoms with z=0, or stage done
          bool expandedAny = false;
          for (Size tb = 0; tb < topLevelTBlsms.size(); ++tb)
          {
            Size b = topLevelTBlsms[tb];
            if ((rprsArr[b] != cNullItm) && (blsmParentArr[b] == cNullItm) &&
                (label[b] == 2) && (dualArr[b] <= EPS))
            {
              dualArr[b] = 0;
              label[b] = 0; treeNumArr[b] = cNullItm;
              Size shiftBd = b - numVtxs;
              for (Size ci = 0; ci < blsmChildrenArr[shiftBd].size(); ++ci) {
                Size child = blsmChildrenArr[shiftBd][ci];
                label[child] = 0; treeNumArr[child] = cNullItm;
                if (child < numVtxs) sttArr[child] = eSttIdle;
              }
              rprsArr[b] = cNullItm;
              rDissolveBlsm3<ItmQue>(graph, blsmArr, blsmParentArr, blsmChildrenArr, b);
              blsmin[b - numVtxs] = cNullItm;
              // B1: refresh topLevelBlsmOfVtx for vertices in b.
              {
                const std::vector<Size>& kids = blsmChildrenArr[b - numVtxs];
                for (Size ci = 0; ci < kids.size(); ++ci) {
                  Size c = kids[ci];
                  if (c < numVtxs) {
                    topLevelBlsmOfVtx[c] = c;
                  } else {
                    tempQ.push(c - numVtxs);
                    while (!tempQ.empty()) {
                      Size bb = tempQ.front(); tempQ.pop();
                      const std::vector<Size>& kk = blsmChildrenArr[bb];
                      for (Size j = 0; j < kk.size(); ++j) {
                        Size gc = kk[j];
                        if (gc < numVtxs) {
                          topLevelBlsmOfVtx[gc] = c;
                          blsmArr[gc] = c;
                        }
                        else tempQ.push(gc - numVtxs);
                      }
                    }
                  }
                }
              }
              pendingFreeIds.push_back(b);
              expandedAny = true;
            }
          }
          if (!expandedAny) { stageDone = true; break; }
          // Fall through to rebuild below
        }

        if (Delta > EPS)
        {
          // Apply dual update with raw duals
          #pragma omp parallel for if(sScanCount > 256) schedule(static)
          for (int su = 0; su < sScanCount; ++su)
          {
            Size u = sScanVtxs[static_cast<Size>(su)];
            Size Bu = topLevelBlsmOfVtx[u];
            if (label[Bu] == 1) {
              dualArr[u] -= Delta;
              if (dualArr[u] < 0 && dualArr[u] > -EPS) dualArr[u] = 0;
            }
          }
          const int tLabelCount = static_cast<int>(tLabelVtxs.size());
          #pragma omp parallel for if(tLabelCount > 256) schedule(static)
          for (int tu = 0; tu < tLabelCount; ++tu)
          {
            Size u = tLabelVtxs[static_cast<Size>(tu)];
            Size Bu = topLevelBlsmOfVtx[u];
            if (label[Bu] == 2) {
              dualArr[u] += Delta;
            }
          }
          #pragma omp parallel for if(topLevelSCount > 256) schedule(static)
          for (int sb = 0; sb < topLevelSCount; ++sb)
          {
            Size b = topLevelSBlsms[static_cast<Size>(sb)];
            if (rprsArr[b] == cNullItm || blsmParentArr[b] != cNullItm || label[b] != 1)
              continue;
            dualArr[b] += 2 * Delta;
          }
          #pragma omp parallel for if(topLevelTCount > 256) schedule(static)
          for (int tb = 0; tb < topLevelTCount; ++tb)
          {
            Size b = topLevelTBlsms[static_cast<Size>(tb)];
            if (rprsArr[b] == cNullItm || blsmParentArr[b] != cNullItm || label[b] != 2)
              continue;
            dualArr[b] -= 2 * Delta;
            if (dualArr[b] < 0 && dualArr[b] > -EPS) dualArr[b] = 0;
          }

          // Expand T-blossoms with z==0
          for (Size tb = 0; tb < topLevelTBlsms.size(); ++tb)
          {
            Size b = topLevelTBlsms[tb];
            if ((rprsArr[b] != cNullItm) && (blsmParentArr[b] == cNullItm) &&
                (label[b] == 2) && (dualArr[b] <= EPS))
            {
              dualArr[b] = 0;
              label[b] = 0; treeNumArr[b] = cNullItm;
              Size shiftBd = b - numVtxs;
              for (Size ci = 0; ci < blsmChildrenArr[shiftBd].size(); ++ci) {
                Size child = blsmChildrenArr[shiftBd][ci];
                label[child] = 0; treeNumArr[child] = cNullItm;
                if (child < numVtxs) sttArr[child] = eSttIdle;
              }
              rprsArr[b] = cNullItm;
              rDissolveBlsm3<ItmQue>(graph, blsmArr, blsmParentArr, blsmChildrenArr, b);
              blsmin[b - numVtxs] = cNullItm;
              // B1: refresh topLevelBlsmOfVtx for vertices in b.
              {
                const std::vector<Size>& kids = blsmChildrenArr[b - numVtxs];
                for (Size ci = 0; ci < kids.size(); ++ci) {
                  Size c = kids[ci];
                  if (c < numVtxs) {
                    topLevelBlsmOfVtx[c] = c;
                  } else {
                    tempQ.push(c - numVtxs);
                    while (!tempQ.empty()) {
                      Size bb = tempQ.front(); tempQ.pop();
                      const std::vector<Size>& kk = blsmChildrenArr[bb];
                      for (Size j = 0; j < kk.size(); ++j) {
                        Size gc = kk[j];
                        if (gc < numVtxs) {
                          topLevelBlsmOfVtx[gc] = c;
                          blsmArr[gc] = c;
                        }
                        else tempQ.push(gc - numVtxs);
                      }
                    }
                  }
                }
              }
              pendingFreeIds.push_back(b);
            }
          }
        }

      }

        // Rebuild: reset scan state, re-init S-roots for Phase 1
        while (!prcsbStk.Empty()) prcsbStk.Pop();
        for (Size i = 0; i < sScanVtxs.size(); ++i)
        {
          Size v = sScanVtxs[i];
          ptrArr[v] = cNullItm;
          sttArr[v] = eSttIdle;
          idxArr[v] = 0;
          label[v] = 0;
          treeNumArr[v] = cNullItm;
          inSScan[v] = 0;
        }
        sScanVtxs.clear();
        for (Size i = 0; i < tLabelVtxs.size(); ++i)
        {
          Size v = tLabelVtxs[i];
          label[v] = 0;
          treeNumArr[v] = cNullItm;
          inTLabel[v] = 0;
        }
        tLabelVtxs.clear();
        for (Size i = 0; i < offsetVtxs.size(); ++i)
        {
          Size v = offsetVtxs[i];
          vtxOffsetState[v] = 0;
          inOffsetVtx[v] = 0;
        }
        offsetVtxs.clear();
        for (Size i = 0; i < topLevelSBlsms.size(); ++i)
          inTopLevelSBlsm[topLevelSBlsms[i] - numVtxs] = 0;
        topLevelSBlsms.clear();
        for (Size i = 0; i < topLevelTBlsms.size(); ++i)
          inTopLevelTBlsm[topLevelTBlsms[i] - numVtxs] = 0;
        topLevelTBlsms.clear();
        for (Size i = 0; i < labeledBlsms.size(); ++i)
        {
          Size b = labeledBlsms[i];
          label[b] = 0;
          treeNumArr[b] = cNullItm;
          inLabeledBlsm[b - numVtxs] = 0;
        }
        labeledBlsms.clear();
        for (Size i = 0; i < blsminTouched.size(); ++i)
        {
          Size b = blsminTouched[i];
          blsmin[b - numVtxs] = cNullItm;
          inBlsminTouched[b - numVtxs] = 0;
        }
        blsminTouched.clear();

        // B1: topLevelBlsmOfVtx is maintained incrementally at each
        // dissolve site; no full rebuild needed here.

        pqD1.clear();
        while (!pqD2.empty()) pqD2.pop();
        while (!pqD3.empty()) pqD3.pop();
        while (!pqD4.empty()) pqD4.pop();

        for (Size fi = 0; fi < freeVtxs.size(); ++fi)
        {
          Size fv = freeVtxs[fi];
          if (mateArr[fv] != cNullItm) continue;
          Size Bfv = topLevelBlsmOfVtx[fv];
          if (label[Bfv] != 0) continue;

          if (Bfv >= numVtxs)
          {
            label[Bfv] = 1; treeNumArr[Bfv] = fv;
            if (inLabeledBlsm[Bfv - numVtxs] == 0)
            {
              inLabeledBlsm[Bfv - numVtxs] = 1;
              labeledBlsms.push_back(Bfv);
            }
            if (inTopLevelSBlsm[Bfv - numVtxs] == 0)
            {
              inTopLevelSBlsm[Bfv - numVtxs] = 1;
              topLevelSBlsms.push_back(Bfv);
            }
            pqD1.insert(fv, dualArr[fv]);
            tempQ.push(Bfv - numVtxs);
            while (tempQ.empty() == false)
            {
              Size tempB2 = tempQ.front(); tempQ.pop();
              for (Size ii = 0; ii < blsmChildrenArr[tempB2].size(); ++ii)
              {
                Size child = blsmChildrenArr[tempB2][ii];
                if (child < numVtxs) {
                  if (sttArr[child] == eSttIdle) {
                    label[child] = 1; treeNumArr[child] = fv;
                    if (inSScan[child] == 0)
                    {
                      inSScan[child] = 1;
                      sScanVtxs.push_back(child);
                    }
                    vtxOffsetState[child] = 1;
                    if (inOffsetVtx[child] == 0)
                    {
                      inOffsetVtx[child] = 1;
                      offsetVtxs.push_back(child);
                    }
                    sttArr[child] = eSttDfsPrcsb; prcsbStk.Push(child);
                  }
                } else {
                  label[child] = 1; treeNumArr[child] = fv;
                  if (inLabeledBlsm[child - numVtxs] == 0)
                  {
                    inLabeledBlsm[child - numVtxs] = 1;
                    labeledBlsms.push_back(child);
                  }
                  tempQ.push(child - numVtxs);
                }
              }
            }
          }
          else
          {
            label[fv] = 1; treeNumArr[fv] = fv;
            pqD1.insert(fv, dualArr[fv]);
            if (inSScan[fv] == 0)
            {
              inSScan[fv] = 1;
              sScanVtxs.push_back(fv);
            }
            vtxOffsetState[fv] = 1;
            if (inOffsetVtx[fv] == 0)
            {
              inOffsetVtx[fv] = 1;
              offsetVtxs.push_back(fv);
            }
            sttArr[fv] = eSttDfsPrcsb; prcsbStk.Push(fv);
          }
        }
      }
    } // end event loop

    if (foundAugPath)
    {
      // Materialize all duals before augmentation.
      //   0 = raw (unlabeled): effective = dualArr[v], no change needed
      //   1 = S-adjusted (dualArr += d_label): effective = dualArr[v] - deltaSum
      //   2 = T-adjusted (dualArr -= d_label): effective = dualArr[v] + deltaSum
      const int offsetCount = static_cast<int>(offsetVtxs.size());
      #pragma omp parallel for if(offsetCount > 1024) schedule(static)
      for (int oi = 0; oi < offsetCount; ++oi)
      {
        Size v = offsetVtxs[static_cast<Size>(oi)];
        if (vtxOffsetState[v] == 1)
        {
          dualArr[v] -= deltaSum;
        }
        else if (vtxOffsetState[v] == 2)
        {
          dualArr[v] += deltaSum;
        }
        vtxOffsetState[v] = 0;
        inOffsetVtx[v] = 0;
        if (dualArr[v] < 0 && dualArr[v] > -EPS) dualArr[v] = 0;
      }
      offsetVtxs.clear();
      // Materialize blossom duals
      const int topLevelSCount = static_cast<int>(topLevelSBlsms.size());
      #pragma omp parallel for if(topLevelSCount > 256) schedule(static)
      for (int sb = 0; sb < topLevelSCount; ++sb)
      {
        Size b = topLevelSBlsms[static_cast<Size>(sb)];
        if (rprsArr[b] == cNullItm || blsmParentArr[b] != cNullItm || label[b] != 1)
          continue;
        dualArr[b] += 2.0 * deltaSum;
        if (dualArr[b] < 0 && dualArr[b] > -EPS) dualArr[b] = 0;
      }
      const int topLevelTCount = static_cast<int>(topLevelTBlsms.size());
      #pragma omp parallel for if(topLevelTCount > 256) schedule(static)
      for (int tb = 0; tb < topLevelTCount; ++tb)
      {
        Size b = topLevelTBlsms[static_cast<Size>(tb)];
        if (rprsArr[b] == cNullItm || blsmParentArr[b] != cNullItm || label[b] != 2)
          continue;
        dualArr[b] -= 2.0 * deltaSum;
        if (dualArr[b] < 0 && dualArr[b] > -EPS) dualArr[b] = 0;
      }

      // Perform augmentation
      Size blsm_sLast = topLevelBlsmOfVtx[sLast];
      Size blsm_tLast = topLevelBlsmOfVtx[tLast];

      bool isCrossTree = (label[blsm_sLast] == 1 && label[blsm_tLast] == 1 &&
                          treeNumArr[blsm_sLast] != treeNumArr[blsm_tLast]);

      if (isCrossTree)
      {
        rAugmentOneSide<ItmIdxdQue, ItmQue>(mateArr, ptrArr, sLast,
                                     blsmArr, rprsArr, blsmParentArr,
                                     blsmChildrenArr, blsmL, numVtxs,
                                     label, blsmin);
        rAugmentOneSide<ItmIdxdQue, ItmQue>(mateArr, ptrArr, tLast,
                                     blsmArr, rprsArr, blsmParentArr,
                                     blsmChildrenArr, blsmL, numVtxs,
                                     label, blsmin);
        mateArr[sLast] = tLast;
        mateArr[tLast] = sLast;
      }
      else
      {
        rAugmentEdmonds<ItmIdxdQue, ItmQue>(mateArr,
                                     ptrArr, expsdQue, sLast, tLast,
                                     blsmArr, rprsArr, blsmParentArr,
                                     blsmChildrenArr, blsmL, numVtxs,
                                     label, blsmin);
      }

      (*card)++;
      Size freeWrite = 0;
      for (Size fi = 0; fi < freeVtxs.size(); ++fi)
      {
        Size fv = freeVtxs[fi];
        if (mateArr[fv] == cNullItm)
          freeVtxs[freeWrite++] = fv;
      }
      freeVtxs.resize(freeWrite);
      augmented = true;
    }
    else
    {
      // Stage ended without augmentation (D1 or INF).
      // Materialize duals using vtxOffsetState.
      const int offsetCount = static_cast<int>(offsetVtxs.size());
      #pragma omp parallel for if(offsetCount > 1024) schedule(static)
      for (int oi = 0; oi < offsetCount; ++oi)
      {
        Size v = offsetVtxs[static_cast<Size>(oi)];
        if (vtxOffsetState[v] == 1)
          dualArr[v] -= deltaSum;
        else if (vtxOffsetState[v] == 2)
          dualArr[v] += deltaSum;
        vtxOffsetState[v] = 0;
        inOffsetVtx[v] = 0;
        if (dualArr[v] < 0 && dualArr[v] > -EPS) dualArr[v] = 0;
      }
      offsetVtxs.clear();
      const int topLevelSCount = static_cast<int>(topLevelSBlsms.size());
      #pragma omp parallel for if(topLevelSCount > 256) schedule(static)
      for (int sb = 0; sb < topLevelSCount; ++sb)
      {
        Size b = topLevelSBlsms[static_cast<Size>(sb)];
        if (rprsArr[b] == cNullItm || blsmParentArr[b] != cNullItm || label[b] != 1)
          continue;
        dualArr[b] += 2.0 * deltaSum;
        if (dualArr[b] < 0 && dualArr[b] > -EPS) dualArr[b] = 0;
      }
      const int topLevelTCount = static_cast<int>(topLevelTBlsms.size());
      #pragma omp parallel for if(topLevelTCount > 256) schedule(static)
      for (int tb = 0; tb < topLevelTCount; ++tb)
      {
        Size b = topLevelTBlsms[static_cast<Size>(tb)];
        if (rprsArr[b] == cNullItm || blsmParentArr[b] != cNullItm || label[b] != 2)
          continue;
        dualArr[b] -= 2.0 * deltaSum;
        if (dualArr[b] < 0 && dualArr[b] > -EPS) dualArr[b] = 0;
      }
    }
  } // end stage loop

  // Compute cardinality and weight
  *card = 0;
  for (Size u = 0; u < numVtxs; ++u)
  {
    if (mateArr[u] == cNullItm) continue;
    Size v = mateArr[u];
    if (u < v)
    {
      if (mateArr[v] != u)
      {
        mateArr[u] = cNullItm;
        mateArr[v] = cNullItm;
      }
      else
      {
        bool found = false;
        const std::vector<Size>& uNbr = vtxVecArr[u];
        for (Size i = 0; i < uNbr.size(); ++i)
        {
          if (uNbr[i] == v) { found = true; break; }
        }
        if (!found)
        {
          mateArr[u] = cNullItm;
          mateArr[v] = cNullItm;
        }
        else
        {
          (*card)++;
        }
      }
    }
  }

  Val totalW = 0;
  for (Size u = 0; u < numVtxs; ++u)
  {
    const std::vector<Size>& uNbrVec = vtxVecArr[u];
    const std::vector<Val>& uWVec = edgWghtVecArr[u];
    for (Size i = 0; i < uNbrVec.size(); ++i)
    {
      Size v = uNbrVec[i];
      if ((u < v) && (mateArr[u] == v))
        totalW += uWVec[i];
    }
  }
  *edgWght = totalW;
}

}

#endif // MATCHING_ENGINE_H
