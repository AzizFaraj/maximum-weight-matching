#ifndef SCALING_SPLIT_FINDMIN_H
#define SCALING_SPLIT_FINDMIN_H

// ============================================================================
// SplitFindmin.h -- the paper's §5.1 #2 "split-findmin for INNER blossoms",
// implemented on the existing ConcatenableQueue (AVL, O(log n)/op). Replaces the
// O(m) per-dissolution edge rescan in ScalingS5Engine: instead of re-scanning all
// m edges to find newly-eligible edges after a blossom dissolves, the min-slack
// edge into each (sub-)blossom is queried in O(log n) via findmin, and a blossom
// dissolves into its odd-cycle components in O(k log n) via split.
//
// Paper ops (Scaling_S5_DataStructures.md L16-18):
//   init(perm of vertices, contiguous per initial blossom)
//   list(u)        = the sublist (root handle) currently containing vertex u
//   split(u)       = split u's sublist right after u into two sublists
//   decreasekey(u, x) = lower key(u) to x (key(u) = min slack of an incident
//                       edge (v,u), v outer; u stores the edge giving key(u))
//   findmin(L)     = (minKey, the vertex achieving it) over sublist L
//
// NOTE on the ConcatenableQueue mapping: the CQ stores LEAVES in left-to-right
// order; each leaf's payload = a vertex id, priority = key(vertex). A "sublist" is
// a CQ queue root. We keep, per vertex, its leaf handle (mLeaf[v]) and its current
// sublist root (mList[v] resolved lazily via findQueue). decreasekey = setPriority;
// findmin = minPriority/minElem; split(u) = CQ split at marker leaf(u).
//
// This is the AVL O(log n) version (per the plan: "start with a simple O(log)
// structure; the fast O(m+n) Gabow-Tarjan version is a later isolated swap").
// Self-contained; depends only on ConcatenableQueue.h + Types.h.
// ============================================================================

#include <vector>
#include "Types.h"
#include "ConcatenableQueue.h"

namespace Matchbox {

class SplitFindmin {
public:
  typedef Size Sz;
#ifdef EVENT_S5_SHARED_SFM_TOPOLOGY
  typedef DualConcatenableQueue Queue;
#else
  typedef ConcatenableQueue Queue;
#endif
  typedef Queue::Handle Handle;
  static const Val INF;   // "no key" sentinel (defined out-of-line)

  SplitFindmin() {}

  // init: the vertex sequence `order` (a permutation) partitioned into contiguous
  // initial sublists by `listEnd[i]==true` marking the LAST vertex of a sublist.
  // All keys start at INF.
  void init(const std::vector<Sz>& order, const std::vector<char>& listEnd) {
    mCQ.clear();
    Sz n = order.size();
    mLeaf.assign(maxId(order) + 1, cNullItm);
    mRoots.clear();
    Handle cur = cNullItm;
    for (Sz i = 0; i < n; ++i) {
      Sz v = order[i];
      Handle h =
#ifdef EVENT_S5_SHARED_SFM_TOPOLOGY
        mCQ.insertRightDual(cur, INF, INF, (Size)v);
#else
        mCQ.insertRight(cur, INF, (Size)v);
#endif
      mLeaf[v] = h;
      if (i < listEnd.size() && listEnd[i]) {           // end of this sublist
        mRoots.push_back(cur);
        cur = cNullItm;
      }
    }
    if (cur != cNullItm) mRoots.push_back(cur);
  }

  // ---- STRUCTURE-AGNOSTIC interface (no backend Handle exposed) ----
  // A "sublist" is identified by ANY vertex currently in it. This is also exactly the
  // Gabow-Tarjan split-findmin interface, so the CQ backend below is a true drop-in:
  // a Gabow-Tarjan SplitFindmin would expose the SAME 4 methods and the engine wiring
  // would not change. (decreasekey/split/findmin/sameList all take vertex ids.)

  // sameList(u, v): are u and v currently in the same sublist?
  bool sameList(Sz u, Sz v) const { return mCQ.findQueue(mLeaf[u]) == mCQ.findQueue(mLeaf[v]); }

  // Opaque identity of u's current list. It is valid until a split
  // touching that list and lets an event client attach a generation to the
  // list rather than to whichever leaf currently attains findmin.
  Sz listId(Sz u) const { return (Sz)mCQ.findQueue(mLeaf[u]); }

  // decreasekey(u, x): lower key(u) to x if x < current.
  Sz decreasekey(Sz u, Val x, bool* changed = 0) {
    Handle h = mLeaf[u];
    return mCQ.decreasePriority(h, x, changed);
  }

  // Begin a new ownership interval for u. Unlike decreasekey, this permits an
  // increase back to INF when a vertex leaves and later re-enters an INNER
  // root; minima from the previous interval must not survive that transition.
  Sz setkey(Sz u, Val x) { return mCQ.setPriority(mLeaf[u], x); }

#ifdef EVENT_S5_SHARED_SFM_TOPOLOGY
  // An independent key channel over the identical laminar-list topology.
  // This avoids maintaining a second AVL tree when two Section-5 roles need
  // different leaf minima but undergo the same split/merge operations.
  Sz decreasekey2(Sz u, Val x, bool* changed = 0) {
    Handle h = mLeaf[u];
    return mCQ.decreasePriority2(h, x, changed);
  }
  Sz setkey2(Sz u, Val x) { return mCQ.setPriority2(mLeaf[u], x); }
  Sz resetkeys2(Val x) { return mCQ.resetPriorities2(x); }
#endif

  // findmin(anyVtxInList): {minKey, vertex achieving it} over the sublist containing it.
  void findmin(Sz anyVtxInList, Val* outKey, Sz* outVtx) const {
    Handle L = mCQ.findQueue(mLeaf[anyVtxInList]);
    if (L == cNullItm) { *outKey = INF; *outVtx = (Sz)cNullItm; return; }
    *outKey = mCQ.minPriority(L);
    *outVtx = (Sz)mCQ.payload(mCQ.minElem(L));
  }
  void findminWithListId(Sz anyVtxInList, Val* outKey, Sz* outVtx,
                         Sz* outListId) const {
    Handle L = mCQ.findQueue(mLeaf[anyVtxInList]);
    *outListId = (Sz)L;
    if (L == cNullItm) { *outKey = INF; *outVtx = (Sz)cNullItm; return; }
    *outKey = mCQ.minPriority(L);
    *outVtx = (Sz)mCQ.payload(mCQ.minElem(L));
  }
#ifdef EVENT_S5_SHARED_SFM_TOPOLOGY
  void findmin2(Sz anyVtxInList, Val* outKey, Sz* outVtx) const {
    Handle L = mCQ.findQueue(mLeaf[anyVtxInList]);
    if (L == cNullItm) { *outKey = INF; *outVtx = (Sz)cNullItm; return; }
    *outKey = mCQ.minPriority2(L);
    *outVtx = (Sz)mCQ.payload(mCQ.minElem2(L));
  }
  void findmin2WithListId(Sz anyVtxInList, Val* outKey, Sz* outVtx,
                          Sz* outListId) const {
    Handle L = mCQ.findQueue(mLeaf[anyVtxInList]);
    *outListId = (Sz)L;
    if (L == cNullItm) { *outKey = INF; *outVtx = (Sz)cNullItm; return; }
    *outKey = mCQ.minPriority2(L);
    *outVtx = (Sz)mCQ.payload(mCQ.minElem2(L));
  }
#endif

  // split(u): split u's sublist immediately AFTER u -> [..,u] becomes one sublist and
  // (u,..] becomes another. Membership is tracked internally (no roots returned).
  void split(Sz u) {
    Handle L = mCQ.findQueue(mLeaf[u]);
    std::vector<Handle> marks; marks.push_back(mLeaf[u]);   // u = last elem of left part
    std::vector<Handle> outs;
    mCQ.split(L, marks, outs);   // outs = [left(..u], right(u..]] -- roots not exposed
  }

  // Legacy/general CQ helper, deliberately outside the paper's split-findmin
  // interface. Canonical ScalingMatching never calls this: Section 5.1 creates
  // one laminar permutation and changes topology only by split as non-outer
  // roots dissolve. Retained for isolated compatibility tests.
  void merge(const std::vector<Sz>& reps) {
    std::vector<Handle> roots; roots.reserve(reps.size());
    for (Sz i = 0; i < reps.size(); ++i) {
      Handle r = mCQ.findQueue(mLeaf[reps[i]]);
      bool dup = false; for (Sz j = 0; j < roots.size(); ++j) if (roots[j] == r) { dup = true; break; }
      if (!dup) roots.push_back(r);
    }
    if (roots.size() > 1) mCQ.merge(roots);
  }

  Val key(Sz u) const { return mCQ.priority(mLeaf[u]); }
#ifdef EVENT_S5_SHARED_SFM_TOPOLOGY
  Val key2(Sz u) const { return mCQ.priority2(mLeaf[u]); }
#endif

  // lastVtx(u): vertex id of the RIGHTMOST leaf in u's current sublist -- the natural pivot to peel off
  // exactly this sublist's run at a later split. O(log n).
  Sz lastVtx(Sz u) const {
    Handle L = mCQ.findQueue(mLeaf[u]);
    Handle r = mCQ.rightmostLeaf(L);
    return (Sz)mCQ.payload(r);
  }

  // sublistVerts(u, out): the vertex ids in u's current sublist, in left-to-right CQ order. Used at a
  // dissolve to find the CURRENT run boundaries (by external membership) and split at them -- robust to
  // any merge-order permutation, since it reads the CQ's actual order. O(sublist size).
  void sublistVerts(Sz u, std::vector<Sz>& out) const {
    out.clear();
    Handle L = mCQ.findQueue(mLeaf[u]);
    if (L == cNullItm) return;
    std::vector<Handle> leaves; mCQ.collectLeaves(L, leaves);
    out.reserve(leaves.size());
    for (Sz i = 0; i < leaves.size(); ++i) out.push_back((Sz)mCQ.payload(leaves[i]));
  }

private:
  static Sz maxId(const std::vector<Sz>& order) {
    Sz m = 0; for (Sz i = 0; i < order.size(); ++i) if (order[i] > m) m = order[i]; return m;
  }
  Queue mCQ;
  std::vector<Handle> mLeaf;    // mLeaf[v] = CQ leaf handle for vertex v
  std::vector<Handle> mRoots;   // initial sublist roots (diagnostics)
};

inline const Val SplitFindmin::INF = (Val)1e18;

} // namespace Matchbox

#endif // SCALING_SPLIT_FINDMIN_H
