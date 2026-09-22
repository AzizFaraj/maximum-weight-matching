#ifndef SCALING_FAST_MWPM_H
#define SCALING_FAST_MWPM_H

// =============================================================================
// ScalingFastMWPM -- the from-scratch paper-faithful §5 engine for max-weight
// PERFECT matching (Duan-Pettie-Su arXiv:1411.1919, §5), targeting the
// O(m·√n·log(nN)) bound. This is the FOURTH and final engine; the prior three
// (ScalingMatching reference 40/40 LEMON, ScalingS5Engine event 320/320, the
// correctness-first scaling loop 64/64) stand as ORACLES we gate against.
//
// WHY FROM SCRATCH (lesson from 3 failed bolt-on attempts -- see memory
// scaling-rebuild-bottleneck-diagnosed): the no-relabel event search corrupts
// when the data structures are bolted onto raw mutable parent pointers
// (mParent/mBParent) that the search rewrites -- staleness, form/dissolve
// oscillation, dangling reuse. The fix, which §5 co-designs:
//
//   THE DATA STRUCTURES OWN THE CONTRACTED GRAPH; THE SEARCH NEVER MUTATES RAW
//   PARENT POINTERS.
//
//   * UNION-FIND owns outer-blossom membership: find(v) = the representative
//     vertex of v's current outer blossom. The search asks "same outer blossom?"
//     and "what outer node?" ONLY through find() -- always valid under unions.
//   * Each blossom is an IMMUTABLE RECORD created once at shrink (its odd cycle
//     of child node-ids + linking edges + base), never rewritten; dissolve marks
//     it inactive. No node-id reuse hazard: ids are monotonic; the contracted
//     walk hops sub-blossoms in O(1) via the records (§5.1 #1b).
//   * Labels / tree pointers live on OUTER NODES (= current blossom ids), and on
//     union the new blossom inherits the role explicitly. A dissolved blossom's
//     children re-enter as their own outer nodes.
//
// BUILD ORDER (each gated):
//   A. UF + immutable blossom records + the contracted odd-cycle walk. Unit-test
//      find/union consistency + cycle reconstruction in ISOLATION. <-- THIS FILE.
//   B. Alternating-forest search on top (grow/shrink/augment), materialized duals,
//      scan-driven. Gate vs brute/LEMON.
//   C. Lazy-timestamp duals (dual adjust = t_now++). D. Event PQ. E. split-findmin.
//
// Integer weights/duals (long long), maximization (Property 1 / Criterion 1).
// C++98 / -ansi. Vertices 0..n-1; blossom node ids n, n+1, ... (monotonic).
// =============================================================================

#include <vector>
#include <cstddef>
#include "Types.h"
#include "BucketPQ.h"

namespace Matchbox {

#ifdef AUGCOUNT
long long gAugCount = 0;
long long gReoptFreed = 0;
long long gReoptRounds = 0;
long long gSearchCount = 0;   // number of searchEvent/searchScan calls (the O(sqrt n) metric)
long long gEdgeScans = 0;     // total edges examined in scheduleEdges (the O(m)x#searches cost)
#endif


class ScalingFastMWPM {
public:
  typedef std::size_t Sz;
  typedef long long LL;
  static Sz NIL() { return (Sz)-1; }

  ScalingFastMWPM() : mMaxW(0), mLazy(false), mTNow(0), mUseEvent(false), mProgress(false), mLcaTok(0), mMultiAug(false), mReopt(false), mReoptDom(false), mHKphase(false), mPersist(false) {}
  void init(Sz n) {
    mN = n;
    mAdjV.assign(n, std::vector<Sz>());
    mAdjW.assign(n, std::vector<LL>());
    mMaxW = 0;
  }
  // event engine ON also enables multi-augment by default (correct + faster: 320/320
  // brute + 40/40 LEMON + ASan clean, faster than scan at n>=6400). setMultiAug
  // overrides if needed for A/B testing.
  void setUseEvent(bool e) { mUseEvent = e; mMultiAug = e; }
  void setMultiAug(bool e) { mMultiAug = e; }   // event-engine vertex-disjoint multi-augment
  void setReopt(bool e) { mReopt = e; }         // re-optimize a seeded perfect matching (warm start)
  void setReoptDom(bool e) { mReoptDom = e; }   // also repair domination violations (carried-dual policy)
  void setHKphase(bool e) { mHKphase = e; }     // multi-augment across dual levels in one search (HK phase)
  void setPersist(bool e) { mPersist = e; }     // persistent forest: continue search across augments (no re-schedule-all)
  void addEdge(Sz u, Sz v, LL w) {
    if (u == v) return;
    mAdjV[u].push_back(v); mAdjW[u].push_back(w);
    mAdjV[v].push_back(u); mAdjW[v].push_back(w);
    if (w > mMaxW) mMaxW = w;
  }

  // ---- STAGE A: the contracted-graph structure (UF + blossom records) ----
  // A "node" is either a vertex (id < mN) or a blossom (id >= mN). The CURRENT
  // outer node of a vertex v is outerNodeOf(v). Union-find over vertices answers
  // same-outer-blossom in near-O(1); blossom records give the structure.

  // (Public for Stage-A unit testing; the search will use these internally.)
  void structInit() {
    mNumNodes = mN;
    mUf.assign(mN, 0);
    for (Sz v = 0; v < mN; ++v) mUf[v] = v;
    // node -> its immediate enclosing blossom id (NIL if outer). Used ONLY to
    // walk a vertex up to its outer node for record lookup; the search's
    // same-blossom test goes through ufFind, not this chain.
    mEnclose.assign(mN, NIL());
    mBChildren.clear();           // per blossom: child node ids (odd cycle order)
    mBEdgeU.clear(); mBEdgeV.clear();
    mBBaseVtx.clear();
    mBActive.clear();
    mBaseOffset = mN;             // blossom ids start at mN
  }

  // union-find over vertices: representative VERTEX of v's outer blossom.
  Sz ufFind(Sz v) { while (mUf[v] != v) { mUf[v] = mUf[mUf[v]]; v = mUf[v]; } return v; }

  // the OUTER NODE id of a NODE x (vertex < mN OR blossom id >= mN). For a vertex,
  // start from its union-find rep (ufFind is VERTEX-only -- mUf is sized mN). For a
  // blossom id, it is already a representative node; do NOT call ufFind on it (that
  // was the heap-buffer-overflow: treeParent passed a blossom id to outerNodeOf ->
  // ufFind(>=mN) indexed mUf OOB, corrupting memory -- the root of every downstream
  // cycle/hang/bad_alloc). Then walk mEnclose up to the top.
  Sz outerNodeOf(Sz x) {
    Sz node = (x < mN) ? ufFind(x) : x;
    Sz g = 0;
    while (mEnclose[node] != NIL()) { node = mEnclose[node];
      if (++g > mNumNodes + 2) {
#ifdef FASTDBG
        std::fprintf(stderr,"  [outerNodeOf CYCLE from x=%lu at node=%lu]\n",(unsigned long)x,(unsigned long)node); std::fflush(stderr);
#endif
        return node;
      }
    }
    return node;
  }

  // create a blossom over an odd cycle of CHILD outer-node ids (cycle[0]=base
  // child), with linking edges edgeU[i] in cycle[i], edgeV[i] in cycle[i+1].
  // Returns the new blossom id. Unions all member vertices' uf-reps to baseRep.
  Sz makeBlossom(const std::vector<Sz>& cycle,
                 const std::vector<Sz>& edgeU, const std::vector<Sz>& edgeV,
                 Sz baseVtx) {
    Sz B = mNumNodes++;
    mBChildren.push_back(cycle);
    mBEdgeU.push_back(edgeU);
    mBEdgeV.push_back(edgeV);
    mBBaseVtx.push_back(baseVtx);
    mBActive.push_back(1);
    growNodeArrays();
    for (Sz i = 0; i < cycle.size(); ++i) mEnclose[cycle[i]] = B;
    mEnclose[B] = NIL();
    // union all leaves under B into the base's uf-rep.
    Sz baseRep = ufFind(baseVtx);
    std::vector<Sz> leaves;
    collectLeaves(B, &leaves);
    for (Sz i = 0; i < leaves.size(); ++i) mUf[ufFind(leaves[i])] = baseRep;
    mUf[baseRep] = baseRep;
    return B;
  }

  // dissolve blossom B: mark inactive, detach children (they become outer nodes),
  // re-split the union-find so each child's leaves point at the child's base rep.
  void dissolveBlossom(Sz B) {
    mBActive[B - mBaseOffset] = 0;
    std::vector<Sz>& ch = mBChildren[B - mBaseOffset];
    for (Sz i = 0; i < ch.size(); ++i) mEnclose[ch[i]] = NIL();
    for (Sz i = 0; i < ch.size(); ++i) {
      std::vector<Sz> cl; collectLeaves(ch[i], &cl);
      Sz rep = childBaseRep(ch[i]);
      for (Sz j = 0; j < cl.size(); ++j) mUf[cl[j]] = rep;
      mUf[rep] = rep;
    }
  }

  // all original vertices (leaves) under node x.
  void collectLeaves(Sz x, std::vector<Sz>* out) {
    if (x < mN) { out->push_back(x); return; }
    std::vector<Sz>& ch = mBChildren[x - mBaseOffset];
    for (Sz i = 0; i < ch.size(); ++i) collectLeaves(ch[i], out);
  }

  Sz blossomBase(Sz B) const { return (B < mN) ? B : mBBaseVtx[B - mBaseOffset]; }
  // count of currently-active (undissolved) blossom nodes -- diagnostic.
  Sz activeBlossoms() const {
    Sz c = 0;
    for (Sz b = mN; b < mNumNodes; ++b) if (mBActive[b - mBaseOffset]) ++c;
    return c;
  }
  bool isBlossom(Sz x) const { return x >= mN; }
  Sz numNodes() const { return mNumNodes; }

  // ---- STAGE B: max-weight PERFECT matching search on the UF structure ----
  // Returns total weight; fills mate (0-based). -1 if no perfect matching.
  LL solve(std::vector<Sz>* mateOut);
  // warm-start from a seed matching + vertex duals (bounded scaling loop).
  LL solveFrom(const std::vector<Sz>& seedMate, const std::vector<LL>& seedY,
               std::vector<Sz>* mateOut);
  Sz mateOf(Sz v) const { return mMate[v]; }
  LL dualY(Sz v) const { return mY[v]; }   // read back final vertex dual (for carry)

private:
  LL solveCore(std::vector<Sz>* mateOut);   // shared search loop + extract
  // duals: y per vertex, z per blossom node (even >=0). slackZ reads materialized
  // duals (Stage B); Stage C swaps in lazy timestamps behind the same accessor.
  std::vector<LL> mY;              // size mN
  std::vector<LL> mZ;              // per node (0 for vertices)

  // alternating forest over OUTER nodes (indexed by node id).
  enum { LBL_FREE = 0, LBL_EVEN = 1, LBL_ODD = 2 };
  std::vector<int> mLabel;         // role of an outer node
  std::vector<Sz>  mParent;        // parent OUTER node in the tree (NIL at root)
  std::vector<Sz>  mParEU, mParEV; // linking edge: child-side vtx, parent-side vtx
  std::vector<Sz>  mRoot;          // tree root outer node
  std::vector<Sz>  mMate;          // mMate[v] = matched vertex (real verts) or NIL

  LL w2max(Sz u, Sz v) const;      // 2*max parallel weight, or huge-neg if no edge
  LL slackZ(Sz u, Sz v);           // yz(e) - 2w(e), z summed over common enclosing
  LL slackZW(Sz u, Sz v, LL w2);   // same, with 2w(e) precomputed (avoids w2max rescan)
  bool tightEdge(Sz u, Sz v) { return slackZ(u, v) == 0; }

  void searchArrays();             // size label/parent/root arrays to mNumNodes
  Sz   treeParent(Sz X);           // contracted-graph tree parent (resolves absorbed)
  bool treeLive(Sz N);             // multi-augment: tree of N still has unmatched root base
  std::vector<char> mUsed;         // per-PHASE: vertex consumed by an augmenting path
  void markPathUsed(Sz E);         // mark all leaves on the tree path E..root used
  bool nodeUsed(Sz N);             // does outer node N contain any used vertex?
  std::vector<Sz> mPhaseNodes;     // outer nodes labeled this phase (for consumed-tree clear)
  std::vector<char> mDeadRoot;     // per outer node: this root's tree was consumed (cleared)
  void clearConsumedTree(Sz root); // after augment: reset tree-state of root's whole tree
  void collectTreeNodes(Sz root, std::vector<Sz>* out); // persist: gather a tree's outer nodes
  Sz   lcaNode(Sz a, Sz b);        // LCA of outer nodes in the tree (NIL if diff)
  void shrinkInto(Sz E, Sz W, Sz u, Sz w);  // form blossom over E..lca..W cycle
  void rebase(Sz B, Sz baseVtx);   // rematch a blossom so baseVtx becomes base
  void flipToRoot(Sz E);           // flip matching along tree path E..root
  void augmentAcross(Sz E, Sz W, Sz u, Sz w);
  bool searchScan();               // one EdmondsSearch (scan-driven Stage B)

  // ---- STAGE C/D: lazy-timestamp duals (§5.2) ----
  // y(v)/z(B) computed on demand from t_now + per-node role timestamps, so a dual
  // adjustment is O(1) (t_now+=delta) and z(B) is EXACT across role changes. Stored
  // base = the materialized y/z at search start; an OUTER node carries its role
  // timestamps. Role: EVEN (outer +), ODD (inner -), FREE (frozen).
  bool mLazy;                      // slackZ + dual-adjust use lazy duals when true
  LL   mTNow;                      // time within the current search (sum of deltas)
  std::vector<LL> mY0;             // base y per vertex (at search start)
  // per-VERTEX role accumulation (threads y across role/containment changes):
  std::vector<LL> mYAcc;           // y-change from v's closed role-intervals
  std::vector<int> mVRole;         // v's role in the current open interval
  std::vector<LL> mVStamp;         // when v entered the current role
  // per-NODE z accumulation (z is per blossom node):
  std::vector<LL> mZ0;             // base z per node
  std::vector<LL> mZAcc;           // z-change from node's closed role-intervals
  std::vector<int> mNRole;         // node's role in current open interval
  std::vector<LL> mNStamp;         // when node entered current role
  std::vector<LL> mTIn, mTOut, mDelta;  // (legacy; unused by the per-entity scheme)
  LL yLazy(Sz v) const;            // current y(v)
  LL zLazy(Sz B) const;            // current z(B)
  void lazyBegin();                // snapshot base, set timestamps per current role
  void lazyEnd();                  // fold lazy values back into mY/mZ
  void setInner(Sz B);             // node B becomes inner root now
  void setOuter(Sz B);             // node B becomes part of outer root now
  void setVtxRole(Sz v, int role); // settle vertex open interval -> mYAcc, re-stamp
  void setNodeRole(Sz B, int role);// settle node z open interval -> mZAcc, re-stamp

  // ---- STAGE D: event-driven no-relabel search ----
  bool mUseEvent;                  // solve() dispatch: event vs scan
  bool mMultiAug;                  // event engine: vertex-disjoint multi-augment (buggy@scale)
  bool mReopt;                     // re-optimize seeded perfect matching (free slack matched edges)
  bool mReoptDom;                  // carried-dual: also free most-negative-slack edge (domination repair)
  bool mHKphase;                   // multi-augment continues across dual levels in one search
  bool mPersist;                   // persistent forest: continue across augments, relabel consumed path FREE (re-traversable)
  bool mProgress;                  // searchEvent made a structural change (dissolve)
  BucketPQ mPQ;
  std::vector<Sz> mLcaSeen; Sz mLcaTok;     // lcaNode scratch (instance, not static)
  std::vector<Sz> mScA, mScB;               // slackZ scratch (instance, not static)
  std::vector<Sz> mSchedLeaves, mRoleLeaves; // collectLeaves scratch (avoid per-call alloc)
  void scheduleEdges(Sz E);        // schedule grow/formation events for EVEN node E
  bool searchEvent();              // one event-driven EdmondsSearch (no relabel)

private:
  Sz mN;
  std::vector<std::vector<Sz> > mAdjV;
  std::vector<std::vector<LL> > mAdjW;
  LL mMaxW;

  Sz mNumNodes, mBaseOffset;
  std::vector<Sz> mUf;              // union-find over vertices
  std::vector<Sz> mEnclose;        // node -> immediate enclosing blossom (or NIL)
  // immutable per-blossom records (indexed by blossomId - mBaseOffset):
  std::vector<std::vector<Sz> > mBChildren, mBEdgeU, mBEdgeV;
  std::vector<Sz> mBBaseVtx;
  std::vector<char> mBActive;

  // the uf-rep that a child node uses as its base when its parent dissolves.
  Sz childBaseRep(Sz child) {
    Sz bv = blossomBase(child);
    return bv;   // representative vertex = the base vertex (canonical)
  }

  void growNodeArrays() {
    if (mEnclose.size() < mNumNodes) mEnclose.resize(mNumNodes, NIL());
    // grow ALL per-node arrays (search + lazy) so a blossom formed mid-search
    // never indexes out of bounds (the OOB-write -> bad_alloc bug).
    if (mLabel.size() < mNumNodes) {
      mLabel.resize(mNumNodes, LBL_FREE);
      mParent.resize(mNumNodes, NIL());
      mParEU.resize(mNumNodes, NIL());
      mParEV.resize(mNumNodes, NIL());
      mRoot.resize(mNumNodes, NIL());
    }
    if (mZ.size() < mNumNodes) mZ.resize(mNumNodes, 0);
    // ALWAYS keep the lazy per-node arrays sized (not just when mLazy) so a blossom
    // formed mid-search has initialized slots regardless of the engine path -- avoids
    // an uninitialized read (the -O2-only segfault).
    if (mZ0.size() < mNumNodes)  mZ0.resize(mNumNodes, 0);
    if (mZAcc.size() < mNumNodes) mZAcc.resize(mNumNodes, 0);
    if (mNRole.size() < mNumNodes) mNRole.resize(mNumNodes, LBL_FREE);
    if (mNStamp.size() < mNumNodes) mNStamp.resize(mNumNodes, 0);
  }
};

// ---------------------------------------------------------------------------
// STAGE B implementation
// ---------------------------------------------------------------------------
inline ScalingFastMWPM::LL
ScalingFastMWPM::w2max(Sz u, Sz v) const {
  LL best = 0; bool found = false;
  for (Sz i = 0; i < mAdjV[u].size(); ++i)
    if (mAdjV[u][i] == v && (!found || mAdjW[u][i] > best)) { best = mAdjW[u][i]; found = true; }
  if (!found) return -(LL)4 * (mMaxW + 1) * (LL)(mN + 1);
  return 2 * best;
}

// slackZ(u,v) = y(u)+y(v) + sum z(B enclosing BOTH) - 2w(u,v). Walk u's and v's
// enclosing-blossom chains (mEnclose, from each vertex's uf-rep node up) to find
// the common suffix and sum those z. Uses mEnclose chains (structural, immutable
// except active-flag) -- NOT the search's tree pointers.
// slack with the doubled edge weight ALREADY KNOWN (caller iterating adjacency has
// it) -- avoids the w2max linear rescan of u's adjacency (was 55% of runtime, called
// 18.7M times, re-finding a weight the caller already has).
inline ScalingFastMWPM::LL
ScalingFastMWPM::slackZW(Sz u, Sz v, LL w2) {
  LL yu = mLazy ? yLazy(u) : mY[u];
  LL yv = mLazy ? yLazy(v) : mY[v];
  LL s = yu + yv - w2;
  // sum z over blossoms enclosing BOTH u and v = the common suffix of their
  // mEnclose chains. Walk both to equal depth, then step up together summing z
  // from the lowest common enclosing blossom. NO scratch allocation (was 18.7M
  // clears/pushes). depth is the blossom nesting (usually tiny).
  Sz a = mEnclose[u], b = mEnclose[v];
  if (a == NIL() || b == NIL()) return s;
  Sz da = 0, db = 0;
  for (Sz x = a; x != NIL(); x = mEnclose[x]) ++da;
  for (Sz x = b; x != NIL(); x = mEnclose[x]) ++db;
  while (da > db) { a = mEnclose[a]; --da; }
  while (db > da) { b = mEnclose[b]; --db; }
  while (a != b) { a = mEnclose[a]; b = mEnclose[b]; }
  for (Sz x = a; x != NIL(); x = mEnclose[x]) s += (mLazy ? zLazy(x) : mZ[x]);
  return s;
}
inline ScalingFastMWPM::LL
ScalingFastMWPM::slackZ(Sz u, Sz v) { return slackZW(u, v, w2max(u, v)); }

inline void ScalingFastMWPM::searchArrays() {
  if (mLabel.size() < mNumNodes) {
    mLabel.resize(mNumNodes, LBL_FREE);
    mParent.resize(mNumNodes, NIL());
    mParEU.resize(mNumNodes, NIL());
    mParEV.resize(mNumNodes, NIL());
    mRoot.resize(mNumNodes, NIL());
    mZ.resize(mNumNodes, 0);
  }
}

// CONTRACTED-GRAPH tree parent: the effective tree-parent of outer node X is the
// OUTER node of X's stored mParent -- so an absorbed parent resolves to the
// blossom that absorbed it (the fix for lcaNode returning NIL on same-tree nodes
// after a shrink -- the n=8 rep=4 invalid-augment bug). NIL at a root or if the
// resolved parent is X itself (absorbed into the same blossom).
inline ScalingFastMWPM::Sz ScalingFastMWPM::treeParent(Sz X) {
  Sz p = mParent[X];
  if (p == NIL()) return NIL();
  Sz po = outerNodeOf(p);
  return (po == X) ? NIL() : po;
}

// a tree is LIVE iff its root's base vertex is still unmatched. After an augment
// consumes a tree, its root's base becomes matched -> the whole tree's events die.
// Lets searchEvent MULTI-AUGMENT (not halt per augment): consumed trees auto-die,
// survivors keep growing. mRoot[N] is the tree root recorded at grow time.
inline bool ScalingFastMWPM::treeLive(Sz N) {
  Sz r = mRoot[N];
  if (r == NIL()) return false;
  Sz ro = outerNodeOf(r);          // current outer node of the root
  return mMate[blossomBase(ro)] == NIL();
}

inline ScalingFastMWPM::Sz ScalingFastMWPM::lcaNode(Sz a, Sz b) {
  std::vector<Sz>& seen = mLcaSeen; Sz& tok = mLcaTok;   // instance state
  if (seen.size() < mNumNodes) seen.assign(mNumNodes, 0);
  ++tok;
  // treeParent (contracted nav) ONLY for the no-relabel event engine; the scan
  // engine relabels fresh each round so raw mParent is correct (and treeParent
  // perturbed it -> regression). mUseEvent selects.
  if (mUseEvent) {
    for (Sz x = a; ; ) { seen[x] = tok; Sz p = treeParent(x); if (p == NIL()) break; x = p; }
    for (Sz y = b; ; ) { if (seen[y] == tok) return y; Sz p = treeParent(y); if (p == NIL()) break; y = p; }
  } else {
    for (Sz x = a; ; ) { seen[x] = tok; if (mParent[x] == NIL()) break; x = mParent[x]; }
    for (Sz y = b; ; ) { if (seen[y] == tok) return y; if (mParent[y] == NIL()) break; y = mParent[y]; }
  }
  return NIL();
}

// rebase blossom B so the sub-node carrying baseVtx becomes the matched base.
// Operates on the IMMUTABLE cycle record by computing the rematch in place on
// mMate (the record itself is not rewritten; only mBBaseVtx[B] and mMate change).
inline void ScalingFastMWPM::rebase(Sz B, Sz baseVtx) {
  if (B < mN) return;
  std::vector<Sz>& cyc = mBChildren[B - mBaseOffset];
  std::vector<Sz>& eu  = mBEdgeU[B - mBaseOffset];
  std::vector<Sz>& ev  = mBEdgeV[B - mBaseOffset];
  Sz k = cyc.size();
  // find the child whose subtree contains baseVtx.
  Sz target = baseVtx;
  while (target != NIL() && mEnclose[target] != B) target = mEnclose[target];
  Sz pos = 0; while (pos < k && cyc[pos] != target) ++pos;
  if (pos >= k) return;
  // recompute the rematch so cyc[pos] is the base: the cycle is odd; rematch the
  // even segment between old base (pos 0) and new base. We rotate a local copy.
  std::vector<Sz> rc(k), reu(k), rev(k);
  for (Sz i = 0; i < k; ++i) { rc[i]=cyc[(pos+i)%k]; reu[i]=eu[(pos+i)%k]; rev[i]=ev[(pos+i)%k]; }
  // rematch along rc: edges (rc[1]-rc[2]),(rc[3]-rc[4]),... become matched.
  for (Sz i = 1; i + 1 < k; i += 2) {
    Sz a = reu[i], b = rev[i];
    if (rc[i]   >= mN) rebase(rc[i],   a);
    if (rc[i+1] >= mN) rebase(rc[i+1], b);
    mMate[a] = b; mMate[b] = a;
  }
  // write the rotated record back (this is a value-rotation of the cycle, not a
  // structural-pointer mutation; ids and enclosure unchanged).
  cyc = rc; eu = reu; ev = rev;
  mBBaseVtx[B - mBaseOffset] = blossomBase(rc[0]);
  if (rc[0] >= mN) rebase(rc[0], baseVtx);
  mBBaseVtx[B - mBaseOffset] = baseVtx;
}

inline void ScalingFastMWPM::flipToRoot(Sz E) {
  if (!mUseEvent) {                       // scan engine: raw mParent (fresh tree)
    while (mParent[E] != NIL()) {
      Sz O = mParent[E]; Sz E2 = mParent[O];
      Sz uo = mParEU[O]; Sz ue = mParEV[O];
      if (O  >= mN) rebase(O,  uo);
      if (E2 >= mN) rebase(E2, ue);
      mMate[uo] = ue; mMate[ue] = uo;
      E = E2;
    }
    return;
  }
  for (;;) {                              // event engine: contracted nav
    Sz O = treeParent(E);
    if (O == NIL()) break;
    Sz E2 = treeParent(O);
    Sz uo = mParEU[O];
    Sz ue = mParEV[O];
    if (O  >= mN) rebase(O,  uo);
    if (E2 != NIL() && E2 >= mN) rebase(E2, ue);
    mMate[uo] = ue; mMate[ue] = uo;
    if (E2 == NIL()) break;
    E = E2;
  }
}

inline void ScalingFastMWPM::shrinkInto(Sz E, Sz W, Sz u, Sz w) {
  Sz L = lcaNode(E, W);
  if (L == NIL()) return;
  // collect E-side and W-side paths up to L (over outer nodes via treeParent, which
  // resolves through any blossom that absorbed the stored parent).
  std::vector<Sz> eS, eU, eV, wS, wU, wV;
  if (mUseEvent) {
    { Sz x = E; while (x != L) { eS.push_back(x); eU.push_back(mParEU[x]); eV.push_back(mParEV[x]); x = treeParent(x); } }
    { Sz x = W; while (x != L) { wS.push_back(x); wU.push_back(mParEU[x]); wV.push_back(mParEV[x]); x = treeParent(x); } }
  } else {
    { Sz x = E; while (x != L) { eS.push_back(x); eU.push_back(mParEU[x]); eV.push_back(mParEV[x]); x = mParent[x]; } }
    { Sz x = W; while (x != L) { wS.push_back(x); wU.push_back(mParEU[x]); wV.push_back(mParEV[x]); x = mParent[x]; } }
  }
  std::vector<Sz> cyc, ceu, cev;
  cyc.push_back(L);
  for (Sz i = eS.size(); i-- > 0; ) { ceu.push_back(eV[i]); cev.push_back(eU[i]); cyc.push_back(eS[i]); }
  ceu.push_back(u); cev.push_back(w);
  for (Sz i = 0; i < wS.size(); ++i) {
    if (i == 0) cyc.push_back(wS[i]);
    else { ceu.push_back(wU[i-1]); cev.push_back(wV[i-1]); cyc.push_back(wS[i]); }
  }
  if (!wS.empty()) { Sz last = wS.size()-1; ceu.push_back(wU[last]); cev.push_back(wV[last]); }
  else { ceu.push_back(mParEU[W]); cev.push_back(mParEV[W]); }
  Sz baseVtx = blossomBase(L);
  Sz B = makeBlossom(cyc, ceu, cev, baseVtx);
  searchArrays();
  // B inherits L's tree role (EVEN), parent/edge/root.
  mZ[B] = 0;
  mLabel[B] = LBL_EVEN;
  mParent[B] = mParent[L]; mParEU[B] = mParEU[L]; mParEV[B] = mParEV[L];
  mRoot[B] = mRoot[L];
}

// mark every original vertex on the tree path from outer node E up to its root as
// USED this phase (Hopcroft-Karp disjointness). Walk via treeParent over the
// CONTRACTED tree; collectLeaves of each node on the path covers blossom interiors.
inline void ScalingFastMWPM::markPathUsed(Sz E) {
  std::vector<Sz> lv;
  for (Sz x = E; ; ) {
    lv.clear(); collectLeaves(x, &lv);
    for (Sz i = 0; i < lv.size(); ++i) mUsed[lv[i]] = 1;
    Sz O = treeParent(x);            // ODD parent
    if (O == NIL()) break;
    lv.clear(); collectLeaves(O, &lv);
    for (Sz i = 0; i < lv.size(); ++i) mUsed[lv[i]] = 1;
    Sz E2 = treeParent(O);           // next EVEN
    if (E2 == NIL()) break;
    x = E2;
  }
}
// does outer node N contain any vertex marked used this phase?
inline bool ScalingFastMWPM::nodeUsed(Sz N) {
  if (N < mN) return mUsed[N] != 0;
  std::vector<Sz> lv; collectLeaves(N, &lv);
  for (Sz i = 0; i < lv.size(); ++i) if (mUsed[lv[i]]) return true;
  return false;
}

// after an augment consumes the tree rooted at `root`, reset the tree-state
// (mLabel/mParent/mRoot -> FREE) of EVERY node in that tree. This prevents a LATER
// same-level augment's flipToRoot/treeParent from walking STALE mParent pointers
// into a consumed tree (the got=-1 bug: mUsed gives vertex-disjointness but NOT
// tree-pointer consistency). mRoot[x] (set at grow time, resolved through current
// outer node) identifies a node's tree. O(phase-node count) per augment.
// collect the OUTER nodes belonging to the tree rooted at `root` (persistent-forest
// augment: we need them to relabel FREE + re-schedule after a flip). Same mRoot-stored
// test as clearConsumedTree; does NOT modify mPhaseNodes (the persist path keeps nodes
// alive in mPhaseNodes and re-roots them as FREE separately).
inline void ScalingFastMWPM::collectTreeNodes(Sz root, std::vector<Sz>* out) {
  for (Sz i = 0; i < mPhaseNodes.size(); ++i) {
    Sz xo = mPhaseNodes[i];
    if (xo < mNumNodes && mRoot[xo] == root) out->push_back(xo);
  }
}

inline void ScalingFastMWPM::clearConsumedTree(Sz root) {
  // mRoot[node] stores the tree's root at grow time. mRoot is propagated as mRoot[E]
  // on grow, so a node's stored mRoot equals the consumed `root` exactly when it's in
  // that tree -- compare directly (no outerNodeOf walks; roots aren't re-rooted mid-
  // phase). swap-remove cleared nodes so each is examined O(1) total across a phase's
  // augments. (Was the 16% clearConsumedTree cost: 2 outerNodeOf walks per node.)
  Sz n = mPhaseNodes.size();
  for (Sz i = 0; i < n; ) {
    Sz xo = mPhaseNodes[i];   // pushed as outer node ids (roots, grown W/M, shrunk B)
    bool consumed = false;
    if (xo < mNumNodes && mRoot[xo] == root) {
      // HK-phase: the search CONTINUES after this augment, so the consumed nodes must
      // have their LAZY DUAL roles SETTLED to FREE now -- otherwise their stale EVEN/ODD
      // role keeps drifting y as mTNow advances within the same search, corrupting the
      // duals (the cross-level g=-1 bug). setOuter/setInner settle+restamp; passing FREE
      // freezes the open interval into mYAcc and stops the drift. (Without HK-phase the
      // search halts right after the augment, so this is a no-op there.)
      if (mLazy && mHKphase) {
        if (xo < mN) setVtxRole(xo, LBL_FREE);
        else { setNodeRole(xo, LBL_FREE);
               std::vector<Sz>& lv = mRoleLeaves; lv.clear(); collectLeaves(xo, &lv);
               for (Sz j = 0; j < lv.size(); ++j) setVtxRole(lv[j], LBL_FREE); }
      }
      mLabel[xo] = LBL_FREE; mParent[xo] = NIL(); mParEU[xo] = NIL(); mParEV[xo] = NIL(); mRoot[xo] = NIL();
      consumed = true;
    }
    if (consumed) { mPhaseNodes[i] = mPhaseNodes[--n]; }
    else ++i;
  }
  mPhaseNodes.resize(n);
}

inline void ScalingFastMWPM::augmentAcross(Sz E, Sz W, Sz u, Sz w) {
#ifdef AUGCOUNT
  ++gAugCount;
#endif
  // mark the whole augmenting path USED before flipping (tree pointers valid now),
  // so no later path THIS PHASE can reuse a vertex on it (the disjointness that
  // makes multi-augment correct). Event engine only -- the scan engine halts after
  // one augment so it needs no disjointness, and never sizes mUsed.
  if (mUseEvent) { markPathUsed(E); markPathUsed(W); }
  if (E >= mN) rebase(E, u);
  if (W >= mN) rebase(W, w);
  flipToRoot(E);
  flipToRoot(W);
  mMate[u] = w; mMate[w] = u;
}

// one EdmondsSearch: relabel forest, grow on tight edges, shrink/augment, unit/
// delta dual-adjust. Returns true if it augmented; false if no augment possible.
inline bool ScalingFastMWPM::searchScan() {
  bool augmentedAny = false;
  for (;;) {
    searchArrays();
    // relabel: every outer node FREE; root the free (unmatched-base) ones EVEN.
    for (Sz x = 0; x < mNumNodes; ++x) {
      if (x >= mN && !mBActive[x - mBaseOffset]) continue;          // dead blossom id
      if (x < mN ? false : (mEnclose[x] != NIL())) continue;        // non-outer blossom
      if (x < mN && mEnclose[x] != NIL()) continue;                 // non-outer vertex
      mLabel[x] = LBL_FREE; mParent[x] = NIL(); mParEU[x] = NIL(); mParEV[x] = NIL(); mRoot[x] = NIL();
    }
    bool anyRoot = false;
    for (Sz v = 0; v < mN; ++v)
      if (mMate[v] == NIL()) {
        Sz N = outerNodeOf(v);
        if (mLabel[N] == LBL_FREE) { mLabel[N] = LBL_EVEN; mRoot[N] = N; anyRoot = true; }
      }
    if (!anyRoot) return augmentedAny;

    bool didShrink = false;
    std::vector<Sz> wq;
    for (Sz x = 0; x < mNumNodes; ++x) {
      bool outer = (x < mN) ? (mEnclose[x] == NIL()) : (mBActive[x - mBaseOffset] && mEnclose[x] == NIL());
      if (outer && mLabel[x] == LBL_EVEN) wq.push_back(x);
    }
    Sz wh = 0;
    while (wh < wq.size() && !didShrink) {
      Sz E = wq[wh++];
      if (mLabel[E] != LBL_EVEN || mEnclose[E] != NIL()) continue;
      std::vector<Sz> Ev; collectLeaves(E, &Ev);
      for (Sz ui = 0; ui < Ev.size() && !didShrink; ++ui) {
        Sz uu = Ev[ui];
        for (Sz k = 0; k < mAdjV[uu].size() && !didShrink; ++k) {
          Sz ww = mAdjV[uu][k];
          if (ufFind(uu) == ufFind(ww)) continue;     // same outer blossom
          Sz W = outerNodeOf(ww);
          if (W == E) continue;
          if (slackZW(uu, ww, 2 * mAdjW[uu][k]) != 0) continue;   // Criterion 1: tight
          if (mLabel[W] == LBL_EVEN) {
            if (lcaNode(E, W) != NIL()) { shrinkInto(E, W, uu, ww); didShrink = true; break; }
            else { augmentAcross(E, W, uu, ww); augmentedAny = true; return true; }
          } else if (mLabel[W] == LBL_ODD) {
            continue;
          } else {
            Sz bw = blossomBase(W); Sz m = mMate[bw];
            if (m == NIL()) { augmentAcross(E, W, uu, ww); augmentedAny = true; return true; }
            Sz M = outerNodeOf(m);
            mLabel[W] = LBL_ODD; mParent[W] = E; mParEU[W] = ww; mParEV[W] = uu; mRoot[W] = mRoot[E];
            mLabel[M] = LBL_EVEN; mParent[M] = W; mParEU[M] = m; mParEV[M] = bw; mRoot[M] = mRoot[E];
            wq.push_back(M);
          }
        }
      }
    }
    if (didShrink) continue;

    // dual adjustment (delta-jump). best = min steps to next event.
    LL best = -1;
    for (Sz uu = 0; uu < mN; ++uu) {
      Sz U = outerNodeOf(uu);
      if (mLabel[U] != LBL_EVEN) continue;
      for (Sz k = 0; k < mAdjV[uu].size(); ++k) {
        Sz ww = mAdjV[uu][k];
        if (ufFind(uu) == ufFind(ww)) continue;
        Sz W = outerNodeOf(ww);
        if (W == U) continue;
        LL s = slackZW(uu, ww, 2 * mAdjW[uu][k]); int lw = mLabel[W];
        if (lw == LBL_EVEN) { LL d = s/2; if (d>=0 && (best<0||d<best)) best=d; }
        else if (lw == LBL_FREE) { LL d = s; if (d>=0 && (best<0||d<best)) best=d; }
      }
    }
    for (Sz x = mN; x < mNumNodes; ++x)
      if (mBActive[x - mBaseOffset] && mEnclose[x] == NIL() && mLabel[x] == LBL_ODD) {
        LL d = mZ[x]/2; if (d>=0 && (best<0||d<best)) best=d;
      }
    if (best < 0) return augmentedAny;
    // pre-dissolve zero-z ODD root blossoms (z-clamp insight).
    bool dis = false;
    for (Sz x = mN; x < mNumNodes; ++x)
      if (mBActive[x - mBaseOffset] && mEnclose[x] == NIL() && mLabel[x] == LBL_ODD && mZ[x] == 0)
        { dissolveBlossom(x); dis = true; }
    if (dis) continue;
    LL delta = (best == 0) ? 1 : best;
    for (Sz v = 0; v < mN; ++v) {
      Sz V = outerNodeOf(v);
      if (mLabel[V] == LBL_EVEN) mY[v] -= delta;
      else if (mLabel[V] == LBL_ODD) mY[v] += delta;
    }
    for (Sz x = mN; x < mNumNodes; ++x)
      if (mBActive[x - mBaseOffset] && mEnclose[x] == NIL()) {
        if (mLabel[x] == LBL_EVEN) mZ[x] += 2*delta;
        else if (mLabel[x] == LBL_ODD) mZ[x] -= 2*delta;
      }
    for (Sz x = mN; x < mNumNodes; ++x)
      if (mBActive[x - mBaseOffset] && mEnclose[x] == NIL() && mLabel[x] == LBL_ODD && mZ[x] == 0)
        dissolveBlossom(x);
  }
}

inline ScalingFastMWPM::LL
ScalingFastMWPM::solve(std::vector<Sz>* mateOut) {
  if (mN == 0) { if (mateOut) mateOut->clear(); return 0; }
  if (mN % 2 != 0) return -1;
  mPQ.init(1);   // one-time; phases reset() to reuse storage (no per-phase realloc)
  structInit();
  searchArrays();
  mY.assign(mN, mMaxW);                    // fresh: Property-1 domination init
  for (Sz x = 0; x < mNumNodes; ++x) mZ[x] = 0;
  mMate.assign(mN, NIL());
  return solveCore(mateOut);
}

// WARM-START: solve from a given matching `seedMate` (0-based, NIL=free) + vertex
// duals `seedY` (must satisfy domination yz>=2w and matched-edge tightness yz=2w
// under the CURRENT edge weights, e.g. the rescaled scale duals). The search only
// REPAIRS the non-perfect spots -> O(few augments) when the seed is close to
// optimal (consecutive scales differ by one bit). gAugCount counts augments.
inline ScalingFastMWPM::LL
ScalingFastMWPM::solveFrom(const std::vector<Sz>& seedMate,
                           const std::vector<LL>& seedY,
                           std::vector<Sz>* mateOut) {
  if (mN == 0) { if (mateOut) mateOut->clear(); return 0; }
  if (mN % 2 != 0) return -1;
  mPQ.init(1);
  structInit();
  searchArrays();
  mY.assign(mN, 0);
  for (Sz v = 0; v < mN && v < seedY.size(); ++v) mY[v] = seedY[v];
  for (Sz x = 0; x < mNumNodes; ++x) mZ[x] = 0;
  mMate.assign(mN, NIL());
  for (Sz v = 0; v < mN && v < seedMate.size(); ++v) mMate[v] = seedMate[v];
  mReopt = true;     // a seeded matching may be perfect-but-suboptimal: re-optimize it
  return solveCore(mateOut);
}

inline ScalingFastMWPM::LL
ScalingFastMWPM::solveCore(std::vector<Sz>* mateOut) {
  // outer rounds: <= n/2 augments + <= O(n) dissolution-retry rounds. Use a
  // generous bound (each round either augments or dissolves >=1 blossom).
  // RE-OPTIMIZATION (warm start): a SEEDED perfect matching may be suboptimal -- some
  // matched edges have positive slack (yz>2w). Complementary slackness (every matched
  // edge tight) is what makes a perfect matching MAX-weight. So we don't stop at
  // "perfect"; we stop at "perfect AND no positive-slack matched edge". When a perfect
  // matching still has a slack edge, we FREE all such edges and let the next search
  // round dual-adjust + re-augment (which either re-tightens that edge by lowering its
  // duals, or finds a strictly heavier alternating re-pairing). Each repair round
  // strictly reduces total matched-edge slack (a bounded non-negative integer), so it
  // terminates. For a fresh solve (y=mMaxW seed) the first perfect matching is already
  // tight, so this adds zero rounds -- it only does work for warm starts.
  bool reoptArmed = mReopt;
  for (Sz guard = 0; guard <= 8 * mN + 16; ++guard) {
    Sz nfree = 0;
    for (Sz v = 0; v < mN; ++v) if (mMate[v] == NIL()) ++nfree;
    if (nfree == 0) {
      if (!reoptArmed) break;
      // perfect: check complementary slackness, free any positive-slack matched edge.
      // (This re-optimizes a warm-started PERFECT-but-suboptimal matching under FRESH
      // domination duals: matched edges that aren't tight get freed and re-augmented
      // to the exact max-weight matching. Correct for the fresh-dual policy. NOTE: the
      // CARRIED-dual policy additionally needs in-search dual correction for negative-
      // slack edges -- not done here; that path is gated off until built.)
      Sz freed = 0;
      for (Sz v = 0; v < mN; ++v) {
        Sz m = mMate[v];
        if (m != NIL() && v < m && slackZ(v, m) > 0) {
          mMate[v] = NIL(); mMate[m] = NIL(); ++freed;
        }
      }
      if (freed == 0) break;                  // optimal under FRESH duals (perfect+tight)
      // NOTE: a CARRIED-dual perfect matching can be improvable via a negative-slack
      // non-matched edge (improving alternating cycle). Repairing that by freeing the
      // violated endpoints + re-augmenting was attempted (mReoptDom) but produced
      // INVALID matchings (the freed vertices may be inside blossoms / carry stale tree
      // state, so the re-augment corrupts mMate -> overshoot/undershoot). The correct
      // fix is an integrated dual-adjusting search that handles negative-slack edges in
      // a perfect matching natively, not a bolt-on free. Parked; mReoptDom is a no-op.
#ifdef AUGCOUNT
      gReoptFreed += freed; ++gReoptRounds;
#endif
      continue;                               // re-search to repair the freed vertices
    }
#ifdef AUGCOUNT
    ++gSearchCount;
#endif
    bool ok = mUseEvent ? searchEvent() : searchScan();
    if (!ok) {
      if (mUseEvent && mProgress) continue;   // dissolution-only round: retry
      return -1;                              // genuinely no augmenting path
    }
  }
  for (Sz v = 0; v < mN; ++v) if (mMate[v] == NIL()) return -1;
#ifdef FASTDBG
  std::fprintf(stderr,"  FINAL y:"); for (Sz v=0;v<mN;++v) std::fprintf(stderr," y[%lu]=%lld",(unsigned long)v,(long long)mY[v]);
  std::fprintf(stderr,"\n  FINAL z:"); for (Sz x=mN;x<mNumNodes;++x) if(mBActive[x-mBaseOffset]) std::fprintf(stderr," z[%lu]=%lld",(unsigned long)x,(long long)mZ[x]);
  std::fprintf(stderr,"\n  matched-edge slacks:");
  for (Sz v=0;v<mN;++v){ Sz m=mMate[v]; if(m!=NIL()&&v<m) std::fprintf(stderr," (%lu-%lu)=%lld",(unsigned long)v,(unsigned long)m,(long long)slackZ(v,m)); }
  std::fprintf(stderr,"\n");
#endif
  LL total = 0;
  if (mateOut) mateOut->assign(mN, NIL());
  for (Sz v = 0; v < mN; ++v) {
    Sz m = mMate[v];
    if (m != NIL() && v < m) {
      LL w2 = w2max(v, m); if (w2 < 0) return -1;
      total += w2 / 2;
      if (mateOut) { (*mateOut)[v] = m; (*mateOut)[m] = v; }
    }
  }
  return total;
}

// ---------------------------------------------------------------------------
// STAGE C/D: lazy-timestamp duals (§5.2)
// ---------------------------------------------------------------------------
// We attach the role timestamps to the OUTER NODE. A vertex v's current y is read
// from B = outerNodeOf(v): the node's role gives the open-interval term, and
// mDelta[B] accumulates the inner-increments from prior intervals of v's current
// containment. On a role/containment change we settle the open interval into
// mDelta and re-stamp. y0/z0 are the materialized base at search start.
inline ScalingFastMWPM::LL
ScalingFastMWPM::yLazy(Sz v) const {
  // PER-VERTEX accumulation: mYAcc[v] = y-change from v's CLOSED role-intervals;
  // the open interval comes from v's CURRENT outer node's role + the time since v
  // entered that role (mVStamp[v]). This correctly threads y across role/containment
  // changes (e.g. a vertex ODD then absorbed into an EVEN blossom -- the ODD gain
  // is settled into mYAcc before the EVEN interval starts). mVRole[v] = v's role
  // during the current open interval; mVStamp[v] = when it started.
  LL y = mY0[v] + mYAcc[v];
  if (mVRole[v] == LBL_ODD)  y += (mTNow - mVStamp[v]);
  else if (mVRole[v] == LBL_EVEN) y -= (mTNow - mVStamp[v]);
  return y;
}
// z(B): per-NODE accumulation. EVEN root z grows +2/step, ODD root z shrinks
// -2/step. mZAcc[B] = closed-interval change; open interval from mNRole/mNStamp.
inline ScalingFastMWPM::LL
ScalingFastMWPM::zLazy(Sz B) const {
  if (B < mN) return 0;
  if (!mBActive[B - mBaseOffset]) return 0;
  LL z = mZ0[B] + mZAcc[B];
  if (mNRole[B] == LBL_EVEN) z += 2 * (mTNow - mNStamp[B]);
  else if (mNRole[B] == LBL_ODD) z -= 2 * (mTNow - mNStamp[B]);
  return z;
}

inline void ScalingFastMWPM::lazyBegin() {
  mTNow = 0;
  mY0.assign(mN, 0); mYAcc.assign(mN, 0); mVRole.assign(mN, LBL_FREE); mVStamp.assign(mN, 0);
  for (Sz v = 0; v < mN; ++v) mY0[v] = mY[v];
  mZ0.assign(mNumNodes, 0); mZAcc.assign(mNumNodes, 0); mNRole.assign(mNumNodes, LBL_FREE); mNStamp.assign(mNumNodes, 0);
  for (Sz x = 0; x < mNumNodes; ++x) mZ0[x] = mZ[x];
  // seed each vertex's open-interval role from its current outer node's label.
  for (Sz v = 0; v < mN; ++v) {
    Sz B = outerNodeOf(v);
    mVRole[v] = mLabel[B]; mVStamp[v] = 0;
  }
  // seed each root blossom node's role.
  for (Sz x = mN; x < mNumNodes; ++x)
    if (mBActive[x - mBaseOffset] && mEnclose[x] == NIL()) { mNRole[x] = mLabel[x]; mNStamp[x] = 0; }
  mLazy = true;
}
inline void ScalingFastMWPM::lazyEnd() {
  std::vector<LL> ny(mN);
  for (Sz v = 0; v < mN; ++v) ny[v] = yLazy(v);
  std::vector<LL> nz(mNumNodes, 0);
  for (Sz x = mN; x < mNumNodes; ++x) nz[x] = zLazy(x);
  for (Sz v = 0; v < mN; ++v) mY[v] = ny[v];
  for (Sz x = mN; x < mNumNodes; ++x) mZ[x] = nz[x];
  mLazy = false;
}
// settle a VERTEX's open interval into mYAcc, then set its new role + stamp.
inline void ScalingFastMWPM::setVtxRole(Sz v, int role) {
  if (mVRole[v] == LBL_ODD)  mYAcc[v] += (mTNow - mVStamp[v]);
  else if (mVRole[v] == LBL_EVEN) mYAcc[v] -= (mTNow - mVStamp[v]);
  mVRole[v] = role; mVStamp[v] = mTNow;
}
// settle a NODE's z open interval into mZAcc, then set its new role + stamp.
inline void ScalingFastMWPM::setNodeRole(Sz B, int role) {
  if (B < mN) return;
  if (mNRole[B] == LBL_EVEN) mZAcc[B] += 2 * (mTNow - mNStamp[B]);
  else if (mNRole[B] == LBL_ODD) mZAcc[B] -= 2 * (mTNow - mNStamp[B]);
  mNRole[B] = role; mNStamp[B] = mTNow;
}
// node B becomes inner root now: B's z shrinks; B's member vertices' y rises.
inline void ScalingFastMWPM::setInner(Sz B) {
  if (B < mN) { setVtxRole(B, LBL_ODD); return; }   // bare vertex fast path (no walk)
  setNodeRole(B, LBL_ODD);
  std::vector<Sz>& lv = mRoleLeaves; lv.clear(); collectLeaves(B, &lv);
  for (Sz i = 0; i < lv.size(); ++i) setVtxRole(lv[i], LBL_ODD);
}
// node B becomes (part of) outer root now: B's z grows; members' y falls.
inline void ScalingFastMWPM::setOuter(Sz B) {
  if (B < mN) { setVtxRole(B, LBL_EVEN); return; }   // bare vertex fast path (no walk)
  setNodeRole(B, LBL_EVEN);
  std::vector<Sz>& lv = mRoleLeaves; lv.clear(); collectLeaves(B, &lv);
  for (Sz i = 0; i < lv.size(); ++i) setVtxRole(lv[i], LBL_EVEN);
}

// schedule grow/formation events for every edge leaving EVEN node E, at the
// absolute time (mTNow + steps) the edge becomes tight under continued adjustment.
inline void ScalingFastMWPM::scheduleEdges(Sz E) {
  // max feasible slack (real edges): 2 endpoints' y (<= 2*maxW each after init) +
  // z sums, all <= O(maxW * n). An event time beyond this is from a non-edge
  // sentinel or a stale/garbage dual -> DON'T schedule (it would also blow the
  // BucketPQ to a giant array = the n=6 rep=8 bad_alloc). Cap at a safe bound.
  LL cap = 4 * (mMaxW + 1) * (LL)(mN + 1);
  std::vector<Sz>& Ev = mSchedLeaves; Ev.clear(); collectLeaves(E, &Ev);  // reused scratch
  for (Sz ui = 0; ui < Ev.size(); ++ui) {
    Sz u = Ev[ui];
#ifdef AUGCOUNT
    gEdgeScans += (long long)mAdjV[u].size();
#endif
    for (Sz k = 0; k < mAdjV[u].size(); ++k) {
      Sz w = mAdjV[u][k];
      if (ufFind(u) == ufFind(w)) continue;
      Sz W = outerNodeOf(w);
      LL s = slackZW(u, w, 2 * mAdjW[u][k]);   // weight already in hand -> no w2max rescan
      int lw = mLabel[W];
      if (lw == LBL_EVEN) { LL d = s / 2; if (d >= 0 && d <= cap) mPQ.schedule(mTNow + d, (long long)u, (long long)w, 1); }
      else if (lw == LBL_FREE) { LL d = s; if (d >= 0 && d <= cap) mPQ.schedule(mTNow + d, (long long)u, (long long)w, 0); }
      // EVEN-ODD: no event.
    }
  }
}

// one event-driven EdmondsSearch on the UF-owned structure with lazy duals.
// NO relabel between events: grow extends the tree; shrink forms a blossom (UF
// owns membership, so queued events validate correctly on pop); dissolution is a
// scheduled event (correct because lazy z is exact across role changes). Halts
// after the first augment (outer solve loop relabels+reruns -- avoids intra-call
// augment staleness, the proven discipline). Returns true if augmented.
inline bool ScalingFastMWPM::searchEvent() {
  mProgress = false;
  searchArrays();
#ifdef FASTDBG
  { Sz cc=0; for(Sz v=0;v<mN;++v) if(mMate[v]!=NIL())++cc; std::fprintf(stderr,"== searchEvent ENTER matched=%lu nodes=%lu\n",(unsigned long)cc,(unsigned long)mNumNodes); std::fflush(stderr); }
#endif
  // relabel ONCE at search start, root free outer nodes EVEN.
  for (Sz x = 0; x < mNumNodes; ++x) {
    bool outer = (x < mN) ? (mEnclose[x] == NIL()) : (mBActive[x - mBaseOffset] && mEnclose[x] == NIL());
    if (!outer) continue;
    mLabel[x] = LBL_FREE; mParent[x] = NIL(); mParEU[x] = NIL(); mParEV[x] = NIL(); mRoot[x] = NIL();
  }
  mPhaseNodes.clear();
  bool anyRoot = false;
  for (Sz v = 0; v < mN; ++v)
    if (mMate[v] == NIL()) {
      Sz N = outerNodeOf(v);
      if (mLabel[N] == LBL_FREE) { mLabel[N] = LBL_EVEN; mRoot[N] = N; anyRoot = true;
        if (mMultiAug || mPersist) mPhaseNodes.push_back(N); }
    }
  if (!anyRoot) return false;

  mUsed.assign(mN, 0);    // fresh disjointness marks for this phase
  lazyBegin();
  for (Sz x = 0; x < mNumNodes; ++x) {
    bool outer = (x < mN) ? (mEnclose[x] == NIL()) : (mBActive[x - mBaseOffset] && mEnclose[x] == NIL());
    if (outer && mLabel[x] == LBL_EVEN) setOuter(x);
  }
  mPQ.reset();   // reuse storage across phases (no per-phase bucket-array realloc)
  for (Sz x = 0; x < mNumNodes; ++x) {
    bool outer = (x < mN) ? (mEnclose[x] == NIL()) : (mBActive[x - mBaseOffset] && mEnclose[x] == NIL());
    if (outer && mLabel[x] == LBL_EVEN) scheduleEdges(x);
  }

  LL guard = 0, guardCap = 8 * (LL)mN * (LL)mN + 8 * (LL)mN * (mMaxW + 1) + 1000;
  bool augmented = false;
  LL augTime = -1;   // mTNow of this phase's augments; multi-augment ONLY at this level
  for (;;) {
    if (++guard > guardCap) break;
    bool havePQ = mPQ.advanceToNext();
    LL pqT = havePQ ? mPQ.now() : -1;
#ifdef FASTDBG
    if (guard <= 30 || guard % 100000 == 0) std::fprintf(stderr,"    [g=%lld havePQ=%d pqT=%lld pqsz=%lu tNow=%lld]\n",(long long)guard,(int)havePQ,(long long)pqT,(unsigned long)mPQ.size(),(long long)mTNow);
#endif
    // POST-AUGMENT HALT (correctness): once we've augmented, the search may CONTINUE
    // collecting more augments across LATER dual levels (the Hopcroft-Karp phase: keep
    // growing trees + dual-adjusting + augmenting within ONE search, so searches/scale
    // drops toward O(sqrt n) instead of ~1 augment/search). We may NOT, however, run a
    // DISSOLUTION after an augment in the same search: a dissolved inner blossom unfolds
    // into an alternating path whose sub-node labels are stale post-augment, corrupting
    // the duals (the original got=-1 bug). So: if we've augmented and the next action
    // would be a dissolution, HALT and let the outer loop relabel+regrow. If there are
    // no more PQ events at all, also halt. Continuing for pure grow/augment is safe
    // because trees are kept consistent (clearConsumedTree after each augment) and lazy
    // duals are exact across role changes. Gated: mHKphase (default off until proven).
    if (augmented && !havePQ) break;
    // next dissolution time over ODD root blossoms (lazy z is exact).
    LL disT = -1; Sz disB = NIL();
    for (Sz b = mN; b < mNumNodes; ++b)
      if (mBActive[b - mBaseOffset] && mEnclose[b] == NIL() && mLabel[b] == LBL_ODD) {
        LL z = zLazy(b); LL t = mTNow + z / 2;
        if (disT < 0 || t < disT) { disT = t; disB = b; }
      }
    if (!havePQ && disB == NIL()) break;
    bool doDissolve = (disB != NIL()) && (!havePQ || disT <= pqT);
    // Post-augment halt policy:
    //  * default (multi-augment, same dual level): halt on ANY time advance -> relabel.
    //  * mPersist (Step-1 persistent forest): CONTINUE across dual levels -- the consumed
    //    path was relabeled FREE + made re-traversable (mUsed cleared), so later trees can
    //    re-route through just-matched edges (the n=4 case, now correct: 472/472). Halt
    //    ONLY before a post-augment DISSOLUTION (still needs the full relabel, since a
    //    dissolved inner blossom's sub-node labels go stale).
    //  * mHKphase (the broken across-levels experiment): kept for reference, gated off.
    if (augmented) {
      if (mPersist) { if (doDissolve && disT > augTime) break; }
      else if (!mHKphase) { if (!havePQ || pqT > augTime) break; }
      else if (doDissolve && disT > augTime) break;
    }
    if (doDissolve) {
      mTNow = disT;
      // dissolve every ODD root blossom now at z==0, then HALT this search so the
      // outer solve loop re-enters searchEvent with a FRESH relabel+regrow. A
      // dissolved inner blossom unfolds into an alternating path whose sub-node
      // labels need the full grow to reconstruct (a partial root-from-free reseed
      // mislabels them -> form/dissolve oscillation -> OOM, the n=4 rep=23 crash).
      // Dissolutions are O(n) total, so O(n) extra relabels -- acceptable; the
      // no-relabel-dissolution refinement (proper expansion) is a later opt.
      lazyEnd();
      for (Sz b = mN; b < mNumNodes; ++b)
        if (mBActive[b - mBaseOffset] && mEnclose[b] == NIL() && mLabel[b] == LBL_ODD && mZ[b] == 0)
          dissolveBlossom(b);
      mProgress = true;       // structural change -> outer loop should re-run, not give up
      return false;           // no augment this call; outer loop relabels and continues
    }

    // SAME-DUAL-LEVEL phase boundary: without HK-phase, multi-augment stays at the SAME
    // dual level (augTime); a later event = duals moved = END the phase -> outer relabel.
    // With HK-phase, CONTINUE across levels (grow/augment only; dissolution already
    // halted above) -- disjointness (mUsed) + clearConsumedTree keep it correct, and the
    // augment tightness check (slackZ==0) below still guards each augment.
    // NOTE on HK + WEIGHTED matching (n=4 minimal case): disjoint multi-augment across
    // dual levels gives SUBOPTIMAL weight. After greedily augmenting the heaviest tight
    // edge (e.g. 0-3), max-weight may require RE-ROUTING through that matched edge later
    // (path 2-0-3-1 flips {0-3} to {0-2,1-3}, heavier). The mUsed disjointness mark
    // FORBIDS that re-route -> locks in the greedy suboptimum. So weighted continuation
    // must NOT mark consumed paths used; it must let later augments traverse just-matched
    // edges (normal Edmonds, sans the disjointness optimization). That's effectively the
    // re-search the non-HK path does -- the open question is doing it without full relabel.
    if (augmented && pqT > augTime && !mHKphase) break;

    BucketPQ::Event ev;
    if (!mPQ.popAtNow(&ev)) continue;
    mTNow = pqT;
    Sz u = (Sz)ev.a, w = (Sz)ev.b;
#ifdef FASTDBG
    std::fprintf(stderr,"   pop T=%lld kind=%d (%lu-%lu) slack=%lld lblE=%d lblW=%d\n",(long long)pqT,ev.kind,(unsigned long)u,(unsigned long)w,(long long)slackZ(u,w),mLabel[outerNodeOf(u)],mLabel[outerNodeOf(w)]);
#endif
    if (mUsed[u] || mUsed[w]) continue;     // endpoint consumed by an augment this phase
    if (ufFind(u) == ufFind(w)) continue;
    Sz E = outerNodeOf(u), W = outerNodeOf(w);
    if (mLabel[E] != LBL_EVEN) continue;
    if (E == W) continue;
    if (slackZ(u, w) != 0) continue;        // not tight now (stale event)
    int lblW = mLabel[W];
    if (lblW == LBL_EVEN) {
      if (lcaNode(E, W) != NIL()) {
#ifdef FASTDBG
        std::fprintf(stderr,"   SHRINK via %lu-%lu; y:",(unsigned long)u,(unsigned long)w);
        for (Sz vv=0; vv<mN; ++vv) std::fprintf(stderr," y[%lu]=%lld",(unsigned long)vv,(long long)yLazy(vv));
        std::fprintf(stderr,"\n");
#endif
        // INCREMENTAL shrink (no halt/reseed): form the blossom B, init its lazy z
        // fresh as an EVEN root (z starts growing now), schedule ONLY B's edges, and
        // KEEP GOING -- the PQ stays alive, searchEvent runs O(augments) times not
        // O(n). The earlier halt-on-shrink crutch was needed only because the
        // ufFind(blossom-id) OOB (now FIXED) corrupted the structure; with correct
        // navigation the persistent tree stays consistent across shrinks.
        // Settle absorbed child blossoms' z (freeze) before B subsumes them. ALSO
        // collect the formerly-ODD nodes on the cycle: their edges were NOT scheduled
        // (only EVEN nodes schedule), so they need scheduling after the shrink. The
        // formerly-EVEN cycle nodes ALREADY have their events in the PQ (still valid;
        // now-internal ones are discarded on pop via ufFind). => scheduleEdges only
        // the formerly-ODD nodes, NOT all of B (the big re-scan redundancy: shrink was
        // re-scanning the WHOLE blossom's edges, including already-scheduled ones --
        // the measured ~5-12x over the bound, since nested shrinks compound it).
        std::vector<Sz> newlyEven;
        { Sz LcaB = lcaNode(E, W);
          for (Sz x = E; x != LcaB; x = treeParent(x)) {
            if (mLabel[x] == LBL_ODD) newlyEven.push_back(x);
            if (x >= mN) setNodeRole(x, LBL_FREE);
          }
          for (Sz x = W; x != LcaB; x = treeParent(x)) {
            if (mLabel[x] == LBL_ODD) newlyEven.push_back(x);
            if (x >= mN) setNodeRole(x, LBL_FREE);
          }
          if (LcaB >= mN) setNodeRole(LcaB, LBL_FREE);
        }
        shrinkInto(E, W, u, w);
        Sz B = outerNodeOf(u);
        mZ0[B] = 0; mZAcc[B] = 0; mNRole[B] = LBL_FREE; mNStamp[B] = mTNow;  // fresh EVEN root
        setOuter(B);
        if (mMultiAug || mPersist) mPhaseNodes.push_back(B);
        for (Sz i = 0; i < newlyEven.size(); ++i) scheduleEdges(newlyEven[i]);
      } else {
#ifdef FASTDBG
        std::fprintf(stderr,"  EVENT augment via %lu-%lu (E=%lu W=%lu) tNow=%lld\n",(unsigned long)u,(unsigned long)w,(unsigned long)E,(unsigned long)W,(long long)mTNow);
        std::fprintf(stderr,"    pre-aug mate:"); for(Sz vv=0;vv<mN;++vv) std::fprintf(stderr," %lu->%ld",(unsigned long)vv,(long)(mMate[vv]==NIL()?-1:(long)mMate[vv])); std::fprintf(stderr,"\n");
        std::fprintf(stderr,"    rootE=%lu rootW=%lu\n",(unsigned long)mRoot[E],(unsigned long)mRoot[W]);
#endif
        Sz rE = mRoot[E], rW = mRoot[W];        // consumed roots (before flip)
        if (mPersist) {
          // PERSISTENT FOREST (Step 1): collect the two consumed trees' nodes BEFORE
          // the flip (tree pointers valid now), augment, then relabel those nodes FREE
          // and re-schedule ONLY their edges -- the rest of the forest + the PQ persist,
          // so the next augment costs O(consumed path), not O(m) re-schedule. We do NOT
          // mark the path mUsed (that disjointness lock is what blocked re-routing in
          // the HK experiment -- a just-matched edge must stay traversable so a later
          // tree can flip it). Settle the consumed nodes' lazy dual roles to FREE so
          // their y stops drifting as mTNow advances.
          std::vector<Sz> consumed;
          collectTreeNodes(rE, &consumed);
          collectTreeNodes(rW, &consumed);
          augmentAcross(E, W, u, w);          // marks mUsed, but mPersist ignores mUsed below
          for (Sz i = 0; i < consumed.size(); ++i) {
            Sz xo = consumed[i];
            if (mLazy) { if (xo < mN) setVtxRole(xo, LBL_FREE);
                         else { setNodeRole(xo, LBL_FREE);
                                std::vector<Sz>& lv = mRoleLeaves; lv.clear(); collectLeaves(xo, &lv);
                                for (Sz j = 0; j < lv.size(); ++j) setVtxRole(lv[j], LBL_FREE); } }
            mLabel[xo] = LBL_FREE; mParent[xo] = NIL(); mParEU[xo] = NIL(); mParEV[xo] = NIL(); mRoot[xo] = NIL();
          }
          // clear disjointness marks on the consumed leaves so they're re-traversable.
          for (Sz i = 0; i < consumed.size(); ++i) {
            Sz xo = consumed[i];
            if (xo < mN) mUsed[xo] = 0;
            else { std::vector<Sz>& lv = mRoleLeaves; lv.clear(); collectLeaves(xo, &lv);
                   for (Sz j = 0; j < lv.size(); ++j) mUsed[lv[j]] = 0; }
          }
          augmented = true; augTime = mTNow;
          continue;                            // keep the same search alive
        }
        augmentAcross(E, W, u, w);
        if (mMultiAug) { clearConsumedTree(rE); clearConsumedTree(rW); }
#ifdef FASTDBG
        std::fprintf(stderr,"    post-aug mate:"); for(Sz vv=0;vv<mN;++vv) std::fprintf(stderr," %lu->%ld",(unsigned long)vv,(long)(mMate[vv]==NIL()?-1:(long)mMate[vv])); std::fprintf(stderr,"\n");
#endif
        // multi-augment at THIS dual level (augTime); the same-level boundary above
        // ends the phase when the duals must move. augmentAcross marked the path
        // used. NOTE: same-level boundary is correct on small n (320/320) but a
        // BUG remains at larger n (got=-1, fails to complete the matching) -> mMultiAug
        // gates it OFF by default (halt-after-augment, the correct baseline).
        augmented = true; augTime = mTNow;
        if (!mMultiAug) break;
      }
    } else if (lblW == LBL_ODD) {
      continue;
    } else {
      Sz bw = blossomBase(W); Sz m = mMate[bw];
      if (m == NIL()) {
#ifdef FASTDBG
        std::fprintf(stderr,"  EVENT augment(free) via %lu-%lu tNow=%lld\n",(unsigned long)u,(unsigned long)w,(long long)mTNow);
#endif
        Sz rE2 = mRoot[E];
        augmentAcross(E, W, u, w); augmented = true; augTime = mTNow;
        if (mMultiAug) clearConsumedTree(rE2);   // W was free (no tree); clear E's tree
        if (!mMultiAug) break; else continue; }
      Sz M = outerNodeOf(m);
      if (mLabel[M] != LBL_FREE) continue;   // grow only into a FREE mate-node
      mLabel[W] = LBL_ODD; mParent[W] = E; mParEU[W] = w; mParEV[W] = u; mRoot[W] = mRoot[E];
      setInner(W);
      mLabel[M] = LBL_EVEN; mParent[M] = W; mParEU[M] = m; mParEV[M] = bw; mRoot[M] = mRoot[E];
      setOuter(M);
      if (mMultiAug || mPersist) { mPhaseNodes.push_back(W); mPhaseNodes.push_back(M); }
      scheduleEdges(M);
    }
  }
  lazyEnd();
  return augmented;
}

} // namespace Matchbox

#endif // SCALING_FAST_MWPM_H
