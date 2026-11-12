#ifndef BLOSSOM_MERGING_H
#define BLOSSOM_MERGING_H

// ============================================================
// BlossomMerging.h
//
// Implementation of the tree-blossom-merging data structure from
// Section 3.1 of:
//   Harold N. Gabow. "Data Structures for Weighted Matching and
//   Extensions to b-matching and f-factors."
//   ACM Trans. Algorithms 14(3), Article 39 (June 2018).
//
// Solves the blossom-merging problem (Section 2.3) in time
//   O(m + n log n) per phase
// under the restriction that make_edge accepts only back-edges
// (Section 3, simplification 2).
//
// This is a STANDALONE abstract data structure with no matchbox
// engine dependencies. The caller (an Edmonds-style matching
// engine) maintains the supporting forest T externally and
// pushes depth information via set_depth(). Edges are back-edges
// from a vertex to one of its ancestors in T.
//
// API operations (paper Section 2.3, paragraph 6):
//   make_blossom(A)  - add set A to O as a new blossom         O(|A|)
//   merge(X, Y) -> Z - combine X and Y into new blossom Z      O(work)
//   make_edge(vw, t) - add edge vw with cost t                 O(1) amortized
//   min_edge         - return min-cost edge between distinct   O(log n) amortized
//                      current blossoms
//
// Plus standalone helpers needed for the matchbox embedding:
//   reset(n)             - reset for n vertices, no blossoms
//   set_depth(v, d)      - record d as depth of v in T
//   blossom_of(v)        - current blossom containing v
//
// Edge classification (Section 3.1, Table 1) into l/s/d types
// charges work to vertices in stems of rank-maximal blossoms
// (Proposition 3.3) and to credits on loose-list edges (I2),
// giving the O(m + n log n) per-phase bound (Theorem 3.5).
//
// Pre-C++11 compatible. No C++98 extensions. Single-threaded.
// ============================================================

#include <vector>
#include <stdexcept>
#include <cassert>
#include <cmath>
#include <cstdio>
#ifdef MBX_SHARED_D3_TIE_ORDER
#include <stdexcept>
#endif
#ifdef MBX_BM_REPLACEMENT_VERIFY
#include <stdexcept>
#include <unordered_set>
#endif
#include "Types.h"
#include "FibonacciHeap.h"

namespace Matchbox {

// Self-contained diagnostic counters (no external definition needed).
inline unsigned long long& bmPacketScanCounter() {
  static unsigned long long c = 0;
  return c;
}
// FORENSIC watchpoint (R_Scope): GABOW_WATCHA/GABOW_WATCHB name two vertices;
// every BM lifecycle event touching a record with that endpoint pair is logged.
// Compile-gated (-DMBX_BM_WATCH): the static-local reads must not tax clean
// timing builds' hot paths.
#ifdef MBX_BM_WATCH
inline long bmWatchA() {
  static long v = -2;
  if (v == -2) { const char* e = std::getenv("GABOW_WATCHA"); v = e ? atol(e) : -1; }
  return v;
}
inline long bmWatchB() {
  static long v = -2;
  if (v == -2) { const char* e = std::getenv("GABOW_WATCHB"); v = e ? atol(e) : -1; }
  return v;
}
#endif // MBX_BM_WATCH
inline unsigned long long& bmMaxPacketLen() {
  static unsigned long long c = 0;
  return c;
}
// DIAG: refreshBlossomSmallest call count and total packet-scan steps inside it.
// Localizes the hugebubbles slowdown — if refresh is called O(removes) and each
// scans many packets, this is the super-linear cost.
inline unsigned long long& bmRefreshCount() {
  static unsigned long long c = 0;
  return c;
}
inline unsigned long long& bmRefreshPacketSteps() {
  static unsigned long long c = 0;
  return c;
}
inline unsigned long long& bmMaxBlossomPackets() {
  static unsigned long long c = 0;
  return c;
}
// DIAG: peak simultaneously-alive edges in mBM. If this is >> 2m on dense
// graphs, mBM is holding far more than the paper's O(m) edges -> the source of
// the super-linear-with-degree memory blowup. mEdges never shrinks, so its
// final size tracks this peak.
inline unsigned long long& bmPeakAliveEdges() {
  static unsigned long long c = 0;
  return c;
}
// DIAG: split the merge() discard volume by REASON, to decide if the
// "merge-in-place" fix is worth it. The proposed fix only helps edges that
// become INTERNAL (both endpoints in Z = sameBlsm) — those can be dropped
// without re-push. The DEDUP-loser discards (smallest-cost-per-other-blossom
// rule) MUST still be re-pushed under matchbox's lazy slacks, so the fix can't
// touch them. If dedup >> sameBlsm, the fix's ceiling is small.
inline unsigned long long& bmDiscardSameBlsm() { static unsigned long long c=0; return c; }
inline unsigned long long& bmDiscardDedup()    { static unsigned long long c=0; return c; }
// R22 routing census + loose-population diagnostics. Since Phase 3 (2026-07-20)
// merge follows the paper's §3.1 routing: a reclassified R0 survivor whose
// charged endpoint u lies OUTSIDE Z goes to loose(Bu); u in Z gathers into
// packet(Z,r). (The historical deviation -- every survivor to packet(Bu) --
// was A/B-tested on identical traces and retired: equal-or-slower, and its
// unbounded-loose rationale predated the loose pairing heap. Git history:
// commit 422a1f5 has the two-variant experiment.)
inline unsigned long long& bmR22ToPacketZ()  { static unsigned long long c=0; return c; } // u in Z -> packet(Z,r)
inline unsigned long long& bmR22ToLoose()    { static unsigned long long c=0; return c; } // u outside Z -> loose(Bu)
inline unsigned long long& bmLooseToPacket() { static unsigned long long c=0; return c; } // container transitions at reroute
inline unsigned long long& bmPacketToLoose() { static unsigned long long c=0; return c; }
inline unsigned long long& bmReclassCount()  { static unsigned long long c=0; return c; } // Table-1 reclassifications at merge
inline unsigned long long& bmLooseLive()     { static unsigned long long c=0; return c; } // current loose-resident edges
inline unsigned long long& bmLooseLivePeak() { static unsigned long long c=0; return c; }
inline unsigned long long& bmMaxLooseLen()   { static unsigned long long c=0; return c; } // max single loose container
// DIAG (R20): count of internal (same-blossom) edges OBSERVED and discarded at
// the min_edge cleanup point. This is NOT a full I1 assertion: zero means only
// that no internal edge surfaced as a blossom minimum during min_edge -- an
// internal edge buried in a packet would go uncounted. The exhaustive check is
// the compile-gated validator (MBX_BM_VALIDATE, validate() below). In paper-pure
// op sequences Prop 3.2(ii) predicts 0 here; >0 under the engine's continuous
// loop is expected (dissolves break rank monotonicity).
inline unsigned long long& bmMinObservedInternalStale() { static unsigned long long c=0; return c; }
#define gBMPacketScan bmPacketScanCounter()

class BlossomMerging {
public:
  typedef Size BlossomId;
  typedef Size VertexId;
  typedef Size EdgeHandle;

  BlossomMerging() : mN(0), mLogN(0), mEdgeRemovalBatchActive(false) {}

  // Reset for n vertices. All depths zero, no blossoms, no edges.
  // Caller should set_depth() for each vertex before make_blossom().
  void reset(Size n);

  // Record depth d for vertex v in the supporting forest T.
  // Caller must call this for each vertex before that vertex appears
  // in a make_blossom or as an endpoint of a make_edge.
  void set_depth(VertexId v, Size d);

  // Section 2.3 make_blossom: add the set of vertices to O as a new
  // blossom. Vertices must not already be in O (paper precondition
  // A ∩ O = ∅). Returns the new BlossomId.
  // The blossom's "root vertex" is the min-depth vertex in `vertices`.
  BlossomId make_blossom(const std::vector<VertexId>& vertices);

  // Section 2.3 make_edge: add back-edge (v, w) with cost t.
  // PRECONDITION: depth(v) > depth(w) (v is a descendant of w in T,
  // i.e., this is a back-edge).
  // PRECONDITION: v and w are both in some current blossom.
  // If v and w are in the same blossom, the edge is discarded.
  // `payload` is an opaque caller-side identifier (e.g., original edge ID).
  // `occurrenceTag` is an opaque caller-side incarnation identifier. It lets
  // an embedding distinguish a still-live record made for an older supporting
  // forest occurrence from a record for the endpoints' current occurrence.
  // Experimental shared-D3 policy: tieRank is the immutable original-edge
  // ordinal, NOT payload, a reduced endpoint, or an allocated record ID.
  // The caller must explicitly supply it when MBX_SHARED_D3_TIE_ORDER is on;
  // cNullItm is reserved for unavailable ranks and is rejected. Feature-off
  // builds ignore it and retain the historical cost-only layout and policy.
  // Returns an EdgeHandle that can be passed to remove_edge later, or
  // cNullItm if the edge was discarded.
  EdgeHandle make_edge(VertexId v, VertexId w, Val t,
                       Size payload = cNullItm, Size occurrenceTag = 0,
                       Size tieRank = cNullItm);

  // Remove a previously inserted edge by handle. The handle is invalidated.
  // No-op if h == cNullItm or already removed.
  void remove_edge(EdgeHandle h);

  // Exact-removal batch used by whole-tree retirement. Container membership
  // and handle generations are updated immediately for every edge; only the
  // affected owner's cached overall minimum and global-H entry are deferred.
  // end_edge_removal_batch() restores those caches before the next BM query.
  void begin_edge_removal_batch();
  bool remove_edge_if_batched(EdgeHandle h, Size gen);
  void end_edge_removal_batch();

  // Generation of an edge slot. Pair (handle, generation) lets a caller hold
  // a reference and later detect whether the slot was freed/reused.
  Size edge_generation(EdgeHandle h) const
  {
    if (h == cNullItm || h >= mEdges.size()) return cNullItm;
    return mEdges[h].generation;
  }

  // Remove edge h ONLY if it is still in use AND its generation matches gen.
  // Guaranteed-safe no-op for stale handles (slot freed or reused). Returns
  // true iff it actually removed the edge.
  bool remove_edge_if(EdgeHandle h, Size gen)
  {
    if (h == cNullItm || h >= mEdges.size()) return false;
    if (!mEdges[h].inUse || mEdges[h].generation != gen) return false;
    remove_edge(h);
    return true;
  }

  // Read an occurrence tag only while the caller's (handle,generation) pair
  // still names the same live edge slot.
  bool edge_occurrence_tag_if(EdgeHandle h, Size gen, Size* outTag) const
  {
    if (h == cNullItm || h >= mEdges.size() || outTag == NULL) return false;
    if (!mEdges[h].inUse || mEdges[h].generation != gen) return false;
    *outTag = mEdges[h].occurrenceTag;
    return true;
  }

  // Read a stored cost only while the caller's (handle,generation) pair still
  // names the same live edge slot (Stage-3 churn repair: lets the engine keep
  // an occurrence-stale record whose PRICE is still current instead of
  // removing and reconstructing it).
  bool edge_cost_if(EdgeHandle h, Size gen, Val* outCost) const
  {
    if (h == cNullItm || h >= mEdges.size() || outCost == NULL) return false;
    if (!mEdges[h].inUse || mEdges[h].generation != gen) return false;
    *outCost = mEdges[h].cost;
    return true;
  }
#ifdef MBX_GABOW_EDGE_WATCH
  void debug_edge_state(FILE* out, const char* label,
                        EdgeHandle h, Size gen) const
  {
    if (out == NULL) return;
    if (h == cNullItm || h >= mEdges.size() || !mEdges[h].inUse ||
        mEdges[h].generation != gen) {
      std::fprintf(out, "BMEDGE %s h=%lu gen=%lu live=0\n", label,
                   (unsigned long)h, (unsigned long)gen);
      return;
    }
    const Edge& E = mEdges[h];
    const BlossomId owner = E.chargedEnd < mVtxBlossom.size()
      ? mVtxBlossom[E.chargedEnd] : cNullItm;
    const bool ownerLive = owner < mBlossoms.size() && mBlossoms[owner].inUse;
    std::fprintf(out,
      "BMEDGE %s h=%lu gen=%lu live=1 v=%lu w=%lu cost=%.17g "
      "charged=%lu packet=%d rank=%lu owner=%lu ownerLive=%d active=%d "
      "smallest=%lu smallestCost=%.17g hHandle=%lu\n",
      label, (unsigned long)h, (unsigned long)gen,
      (unsigned long)E.v, (unsigned long)E.w, (double)E.cost,
      (unsigned long)E.chargedEnd, E.inPacket ? 1 : 0,
      (unsigned long)E.packetRank, (unsigned long)owner, ownerLive ? 1 : 0,
      ownerLive && mBlossoms[owner].active ? 1 : 0,
      ownerLive ? (unsigned long)mBlossoms[owner].smallestOverall : 0UL,
      ownerLive ? (double)mBlossoms[owner].smallestCost : 0.0,
      ownerLive ? (unsigned long)mBlossoms[owner].hHandle : 0UL);
  }
#endif

  // Section 2.3 merge: combine blossoms X and Y into a new blossom Z.
  // X and Y are destroyed. Returns Z's id.
  // PRECONDITION: X != Y, both currently exist.
  // If `outDiscardedPayloads` is non-NULL, payloads of edges discarded by
  // I[B] dedup (smallest-edge-per-other-blossom rule) get appended. This
  // lets matchbox restore mD3Queued state so refresh can re-push these
  // edges when their stored cost is no longer the cheapest (which happens
  // under matchbox's lazy-delta scheme where slacks change post-merge).
  BlossomId merge(BlossomId X, BlossomId Y,
                  std::vector<Size>* outDiscardedPayloads = NULL);

  // Extension: dispose a blossom (frees its slot, all its edges,
  // all its packets, removes from H, clears mVtxBlossom for its
  // vertices). Used when matchbox's splitBlossom (expand) un-shrinks
  // a previously-merged blossom — the merge cannot be reversed in
  // BlossomMerging's data structure, so we dispose Z and the caller
  // re-make_blossom's the sub-blossoms with their atomic vertex sets.
  // Any EdgeHandle that referred to an edge in Z's containers is
  // invalidated; remove_edge on such handles is a safe no-op.
  // If outDiscardedPayloads is non-NULL, payloads of all destroyed edges
  // are appended so the engine can clear mD3Queued and let the matchbox
  // refresh re-push them when needed.
  void dispose_blossom(BlossomId Z,
                       std::vector<Size>* outDiscardedPayloads = NULL);

  // Section 2.3 min_edge: return the smallest-cost edge currently
  // joining two distinct blossoms. Returns false if no such edge.
  // `outPayload` receives the opaque payload passed to make_edge.
  // `outHandle` (if non-NULL) receives the EdgeHandle for the returned edge —
  // use this to call remove_edge on the exact entry rather than relying on
  // separately-tracked handles that may have been overwritten by subsequent
  // make_edge calls for the same (u, edgeIdx).
  bool min_edge(VertexId& outV, VertexId& outW, Val& outCost,
                Size* outPayload = NULL, EdgeHandle* outHandle = NULL);

  // Find current blossom containing v. Returns cNullItm if v not yet
  // in any blossom.
  BlossomId blossom_of(VertexId v) const
  {
    assert(v < mVtxBlossom.size());
    return mVtxBlossom[v];
  }

  // Size and rank of a current blossom.
  Size blossom_size(BlossomId b) const
  {
    assert(b < mBlossoms.size() && mBlossoms[b].inUse);
    return mBlossoms[b].size;
  }

  // True if blossom id is currently a live blossom.
  bool blossom_alive(BlossomId b) const
  {
    return b < mBlossoms.size() && mBlossoms[b].inUse;
  }
  Size blossom_rank(BlossomId b) const
  {
    assert(b < mBlossoms.size() && mBlossoms[b].inUse);
    return mBlossoms[b].rank;
  }

  // Diagnostics.
  Size total_blossoms_alive() const { return mAliveBlossoms; }
  Size total_edges_alive() const { return mAliveEdges; }

  // A current outer blossom participates in the global min structure. Inner
  // or matched blossoms retain their edges lazily but their header is dormant,
  // so a whole stale packet family is skipped without per-edge eviction.
  void set_vertex_blossom_active(VertexId v, bool active);

#ifdef MBX_BM_VALIDATE
  // ==== EXHAUSTIVE INVARIANT VALIDATOR (compile-gated; Phase 2) ====
  // Absent from production/timing builds: the entire declaration+definition is
  // inside #ifdef MBX_BM_VALIDATE, so clean builds pay zero overhead.
  //
  // Checks, by invariant group (naming used in failure output):
  //  I1i-*   every live edge reachable from EXACTLY one container of one live
  //          blossom; endpoints in bounds; endpoint blossoms live; endpoint
  //          blossoms DISTINCT when paperPure (strict I1(i): internal edges
  //          are impossible in paper-pure op sequences); no freed edge
  //          reachable; no live edge unreachable; free-list/live separation
  //          for edges, blossoms and packets; alive counters exact.
  //  CLASS-* chargedEnd is one of {v,w}; owner (derived mVtxBlossom[chargedEnd])
  //          == the blossom whose container holds the edge (the loose "waiting
  //          area" ownership rule and the packet ownership rule alike);
  //          depth(v) > depth(w); rEdge == logFloor(depth diff) (depths are
  //          immutable, so this must ALWAYS be current). The l/s/d type and
  //          packetRank may be INTENTIONALLY STALE (the paper reclassifies
  //          only when an edge enters R0): recomputing them from CURRENT
  //          blossoms and differing is counted in *outStaleClassifications,
  //          not flagged. What must hold regardless of staleness: packet
  //          edges sit in a packet whose rank equals their stored packetRank,
  //          and that rank is legal for the owner (r(B) < rank <= logN cap).
  //  PH-*    pairing-heap shape: root has null phPrev/phNext; leftmost child's
  //          phPrev is the parent; sibling links doubly consistent; heap order
  //          cost(parent) <= cost(child); no cycles (visit stamping).
  //  PKT-*   packet header ownership (blossom field), one header per nonempty
  //          (B, rank), count == traversed edge count, smallestEdge == the
  //          traversal minimum, empty-header consistency (count==0 <=>
  //          smallestEdge==cNullItm).
  //  RING-*  blossom vertex rings: doubly-linked consistency, size exact,
  //          rank == logFloor(size), rootVertex is a min-depth ring member,
  //          mVtxBlossom agreement, every mapped vertex in exactly one ring.
  //  MIN-*   blossom minimum: smallestOverall/smallestCost equal the true
  //          minimum over the blossom's loose root and packet roots and the
  //          cached handle lives in this blossom's containers.
  //  H-*     the global heap: an ACTIVE blossom with edges has an H entry
  //          whose key EXACTLY equals smallestCost; inactive or edgeless
  //          blossoms have none; H size == number of live handles; the H
  //          minimum equals an independently recomputed global minimum over
  //          active blossoms (dormant blossoms excluded by construction).
  // Complexity: O(slots + m) per call plus O(#packets) bookkeeping; check-mode
  // only. Returns the number of violations (0 = pass); each violation prints
  // one stderr line: VIOLATION[name] ids expected/observed.
  Size validate(bool paperPure, Size* outStaleClassifications = NULL) const;
  // Enumerate all live stored edges, for harness-side oracles (I1(ii)).
  void collect_live_edges(std::vector<VertexId>& outV, std::vector<VertexId>& outW,
                          std::vector<Val>& outC, std::vector<Size>& outPayload) const;
private:
  // Walk one pairing heap (a loose list or one packet), stamping edges and
  // checking PH-, CLASS- and I1i- properties. Returns edges seen; appends
  // costs' minimum via ioMin*/ioHave. `seen` doubles as the cycle guard.
  Size validateHeapWalk(EdgeHandle root, BlossomId owner, bool expectInPacket,
                        Size expectRank, bool paperPure,
                        std::vector<unsigned char>& seen,
                        std::vector<BlossomId>& edgeOwner,
                        Val& ioMin, EdgeHandle& ioMinE, bool& ioHave,
                        Size& viol, Size* staleCls) const;
public:
#endif // MBX_BM_VALIDATE

private:
  // ---- Edge classification (Section 3.1, Table 1) ----
  enum EdgeType {
    TYPE_L = 0,  // long:  r(e) > r(A_v)
    TYPE_S = 1,  // short: r(e) <= r(A_v), r(A_v) <= r(A_w)
    TYPE_D = 2   // down:  A_v != A_w, r(A_v) > r(A_w)
  };

  // Edge record. Edges that have been make_edge'd live in either a
  // packet ring or a loose list of exactly one current blossom (the
  // owner = blossom containing the classification's u endpoint).
  struct Edge {
    VertexId v;          // deeper endpoint (depth(v) > depth(w))
    VertexId w;          // shallower endpoint
    Val cost;            // the cost t passed to make_edge
    Size payload;        // opaque caller-supplied identifier
#ifdef MBX_SHARED_D3_TIE_ORDER
    Size tieRank;        // immutable original-edge ordinal, never slot/endpoint ID
#endif
    Size rEdge;          // floor(log(depth(v) - depth(w))), fixed at make_edge
    EdgeType type;       // current classification
    VertexId chargedEnd; // v or w per Table 1
    Size packetRank;     // current rank classification r (per Table 1)
                         // (owner blossom + packet are DERIVED, never stored: owner =
                         //  mVtxBlossom[chargedEnd]; packet = owner's packet of rank
                         //  packetRank. Storing them forced O(k) rebrand walks at merge,
                         //  which destroyed the paper's O(1) packet concatenation.)
    bool inPacket;       // true if in a packet pairing heap, false if in loose list
    Size occurrenceTag;  // opaque caller-side supporting-forest incarnation
    // Intrusive PAIRING HEAP links (packet membership). Pairing heaps give the three
    // properties the paper's rings + our eviction jointly need: O(1) meld (the paper's
    // O(1) packet concatenation), O(log) amortized arbitrary delete (eviction), O(1)
    // find-min (the cached smallest). phPrev = left sibling, or PARENT if leftmost.
    EdgeHandle phChild;
    EdgeHandle phNext;
    EdgeHandle phPrev;
    Size generation;     // bumped each time this slot is freed; lets callers
                         // hold (handle, gen) and detect stale handles safely
    bool inUse;
  };

  // A packet of a blossom holds edges of a single rank class.
  // Implementation: doubly-linked circular ring of edges (Section 3.1
  // paragraph 14: "a ring, i.e., a circularly linked list").
  // smallestEdge cached for O(1) read by F-heap maintenance.
  struct Packet {
    Size rank;            // packet rank r (in [r(B)+1, log n])
    BlossomId blossom;    // owning blossom (updated O(1) at wholesale packet moves)
    EdgeHandle smallestEdge;  // the PAIRING-HEAP ROOT (== packet minimum, O(1) read)
    Size nextInBlossom;   // next Packet in the blossom's packet list
    Size count;           // #edges in this packet (diagnostics + logical-freight metric)
    bool inUse;
  };

  // A current blossom record.
  struct Blossom {
    Size size;            // number of vertices in this blossom
    Size rank;            // floor(log(size))
    VertexId rootVertex;  // min-depth vertex in V(B)
    VertexId vtxListHead; // head of doubly-linked vertex list (cNullItm if none)
    EdgeHandle looseHead; // head of the loose list (cNullItm if empty)
    EdgeHandle looseSmallest;    // cached smallest edge in the loose list
    Val looseSmallestCost;       // its cost (valid iff looseSmallest != cNullItm)
    Size looseCount;      // #edges in the loose container (diag, like Packet::count)
    EdgeHandle smallestOverall;
    Val smallestCost;
    Size firstPacket;     // head of packet list in mPackets (cNullItm if none)
    Size hHandle;         // handle in F-heap H (cNullItm if not present)
    bool active;          // false => retained lazily but absent from global H
    bool inUse;
  };

  // Internal Design-B contract, used only by the Gabow split mirror below.
  // Caller guarantees that subsequent child lists are an exact disjoint
  // partition of Z's current atoms. No BM ownership/ring consumer, callback,
  // validation, reset, merge or edge removal is allowed until all child makes
  // finish. There is no rollback; a diagnostic failure aborts the computation.
  // Numeric old-owner equality alone is NOT a duplicate check: the first
  // child deliberately reuses Z's slot under the original LIFO allocation.
  friend class FastMatchingEngineV2;

  // Experimental largest-child lifecycle. No BM consumer may run between
  // drain and finish. Parent remains allocated; only smaller atoms change owner.
  void begin_largest_child_partition(BlossomId Z, std::vector<Size>* discarded);
  template<class Members> BlossomId make_smaller_for_largest_partition(
      BlossomId Z, const Members& members);
  void finish_largest_child_partition(BlossomId Z, VertexId canonicalHead,
      Size remaining, VertexId minimumRoot);
#ifdef MBX_BM_REPLACEMENT_TEST_ACCESS
  friend struct BmReplacementTest;
#endif
  BlossomId make_blossom_impl(const std::vector<VertexId>& vertices,
                              BlossomId supersededOwner);
  void dispose_blossom_impl(BlossomId Z,
                            std::vector<Size>* outDiscardedPayloads,
                            bool clearVertexRing);
  BlossomId make_for_frontier_partition(const std::vector<VertexId>& vertices,
                                        BlossomId supersededOwner)
  {
    assert(supersededOwner != cNullItm);
    return make_blossom_impl(vertices, supersededOwner);
  }
  void dispose_for_frontier_partition(BlossomId Z,
                                      std::vector<Size>* outDiscardedPayloads)
  {
    assert(blossom_alive(Z) && !mEdgeRemovalBatchActive);
    dispose_blossom_impl(Z, outDiscardedPayloads, false);
  }
#ifdef MBX_BM_REPLACEMENT_VERIFY
  void verify_frontier_partition(BlossomId Z,
      const std::vector<const std::vector<VertexId>*>& children) const;
  void verify_frontier_result(
      const std::vector<const std::vector<VertexId>*>& children,
      const std::vector<BlossomId>& ids) const;
#endif


  Size mN;          // total vertex count (capacity)
  Size mLogN;       // floor(log mN) + 1, max rank value

  std::vector<Size> mDepth;            // depth in supporting forest T
  std::vector<BlossomId> mVtxBlossom;  // current blossom of each vertex
  std::vector<VertexId> mVtxNext;      // doubly-linked vertex list per blossom
  std::vector<VertexId> mVtxPrev;

  std::vector<Blossom> mBlossoms;
  std::vector<Edge> mEdges;
  std::vector<Packet> mPackets;

  std::vector<BlossomId> mFreeBlossoms;
  std::vector<EdgeHandle> mFreeEdges;
  std::vector<Size> mFreePackets;
  std::vector<BlossomId> mDeferredRefreshBlossoms;
  bool mEdgeRemovalBatchActive;

  Size mAliveBlossoms;
  Size mAliveEdges;

  // The global F-heap H over blossoms.
  // Key = blossom's smallestCost; Payload = BlossomId.
  // Lazy deletion: when a blossom is merged, its H-handle is removed
  // and the new blossom's H-handle is inserted.
#ifdef MBX_SHARED_D3_TIE_ORDER
  // Finite exact lexicographic order. No EPS equivalence appears in a heap
  // comparator. Physical twins may be equivalent if cost and rank coincide.
  struct MinKey {
    Val cost;
    Size tieRank;
    MinKey(Val c = 0, Size r = cNullItm) : cost(c), tieRank(r) {}
    bool operator<(const MinKey& rhs) const {
      if (cost < rhs.cost) return true;
      if (rhs.cost < cost) return false;
      return tieRank < rhs.tieRank;
    }
    bool operator!=(const MinKey& rhs) const {
      return *this < rhs || rhs < *this;
    }
  };
  MinKey edgeKey(EdgeHandle e) const {
    return MinKey(mEdges[e].cost, mEdges[e].tieRank);
  }
  MinKey blossomKey(const Blossom& b) const {
    return MinKey(b.smallestCost, mEdges[b.smallestOverall].tieRank);
  }
  static Val keyCost(const MinKey& k) { return k.cost; }
#else
  typedef Val MinKey;
  MinKey edgeKey(EdgeHandle e) const { return mEdges[e].cost; }
  MinKey blossomKey(const Blossom& b) const { return b.smallestCost; }
  static Val keyCost(MinKey k) { return k; }
#endif
  bool edgeLess(EdgeHandle a, EdgeHandle b) const {
    return edgeKey(a) < edgeKey(b);
  }
  FibonacciHeap<MinKey, BlossomId> mH;

  // Auxiliary I[1..log n] array for gather operations during merge.
  // mGatherHead[r] is the head of a singly-linked list of edges with
  // bucket index r (used to assemble new packet of rank r).
  // mGatherTouched lists the indices that were touched, for O(touched)
  // cleanup.
  std::vector<EdgeHandle> mGatherHead;
  std::vector<Size> mGatherTouched;
  std::vector<EdgeHandle> mGatherEdgeNext;  // singly-linked: next in gather list
  // Scratch buffer for merge dedup. Reused across merges to avoid
  // O(numBlossoms) allocation per merge. Touched entries listed in
  // mDedupTouched for O(|touched|) cleanup at end of merge.
  std::vector<EdgeHandle> mDedupBestPerB;
  // Scratch used by the two-pass pairing-heap child combine.
  std::vector<EdgeHandle> mSweepScratch;

  // ---- internal helpers (defined below) ----

  static Size logFloor(Size x);
  BlossomId allocBlossom();
  void freeBlossom(BlossomId id);
  EdgeHandle allocEdge();
  void freeEdge(EdgeHandle e);
  Size allocPacket();
  void freePacket(Size p);
  void removeEdgeImpl(EdgeHandle h, bool deferRefresh);

  // Compute classification (type, chargedEnd, packetRank) for an edge
  // given its current endpoint blossoms. See Section 3.1, Table 1.
  void classifyEdgeFromBlossoms(EdgeHandle e, BlossomId Bv, BlossomId Bw,
                                EdgeType& outType, VertexId& outU, Size& outR);

  // Add edge e to loose list of blossom b. Updates blossom's smallest
  // if needed and refreshes H entry.
  void addToLoose(BlossomId b, EdgeHandle e);

  // Add edge e to packet of blossom b at rank r. Creates the packet if
  // needed. Updates packet's smallest and blossom's smallest.
  void addToPacket(BlossomId b, Size rank, EdgeHandle e);

  // Per-packet intrusive pairing-heap helpers (O(1) concatenation, Section 3.1).
  EdgeHandle phMeld(EdgeHandle a, EdgeHandle b);        // O(1)
  void phDetachNonRoot(EdgeHandle e);                    // O(1) unlink from parent/sibs
  EdgeHandle phCombineChildren(EdgeHandle first);        // two-pass pairing, amortized O(log)
  void packetInsert(Packet& P, EdgeHandle e);            // meld: O(1)
  void packetRemove(Packet& P, EdgeHandle e);            // arbitrary delete: amortized O(log)
  void packetCollect(EdgeHandle root, std::vector<EdgeHandle>& out);  // O(k) traversal
  Size findPacketOf(const Edge& E);                      // derive: owner's packet of E.packetRank

  // Recompute blossom b's smallestOverall and smallestCost by scanning
  // its packet smallest values and loose list. Then refresh H entry.
  void refreshBlossomSmallest(BlossomId b);
  void recomputeLooseSmallest(BlossomId b);
  // Re-evaluate H entry for blossom b. If b is not in H, insert it.
  // If b is in H and its smallestCost is smaller than the heap key,
  // decreaseKey. If b is in H and the cost is larger, remove and re-insert.
  void refreshHEntry(BlossomId b);

  // Set vertex v's current blossom to b.
  void setVtxBlossom(VertexId v, BlossomId b) { mVtxBlossom[v] = b; }
};

// ---- inline implementations ----

inline Size BlossomMerging::logFloor(Size x)
{
  if (x == 0) return 0;
  Size r = 0;
  while ((x >>= 1) != 0) ++r;
  return r;
}

inline void BlossomMerging::reset(Size n)
{
  mN = n;
  mLogN = (n == 0) ? 0 : (logFloor(n) + 1);
  mDepth.assign(n, 0);
  mVtxBlossom.assign(n, cNullItm);
  mVtxNext.assign(n, cNullItm);
  mVtxPrev.assign(n, cNullItm);
  mBlossoms.clear();
  mEdges.clear();
  mPackets.clear();
  mFreeBlossoms.clear();
  mFreeEdges.clear();
  mFreePackets.clear();
  mDeferredRefreshBlossoms.clear();
  mEdgeRemovalBatchActive = false;
  mAliveBlossoms = 0;
  mAliveEdges = 0;
  mH.clear();
  mGatherHead.assign(mLogN + 2, cNullItm);
  mGatherTouched.clear();
  mGatherEdgeNext.clear();
}

inline void BlossomMerging::set_depth(VertexId v, Size d)
{
  assert(v < mDepth.size());
  mDepth[v] = d;
}

inline BlossomMerging::BlossomId BlossomMerging::allocBlossom()
{
  BlossomId id;
  if (!mFreeBlossoms.empty()) {
    id = mFreeBlossoms.back();
    mFreeBlossoms.pop_back();
  } else {
    Blossom b;
    b.inUse = false;
    mBlossoms.push_back(b);
    id = mBlossoms.size() - 1;
  }
  Blossom& B = mBlossoms[id];
  B.size = 0;
  B.rank = 0;
  B.rootVertex = cNullItm;
  B.vtxListHead = cNullItm;
  B.looseHead = cNullItm;
  B.looseSmallest = cNullItm;
  B.looseSmallestCost = 0;
  B.looseCount = 0;
  B.smallestOverall = cNullItm;
  B.smallestCost = 0;
  B.firstPacket = cNullItm;
  B.hHandle = cNullItm;
  B.active = true;
  B.inUse = true;
  ++mAliveBlossoms;
  return id;
}

inline void BlossomMerging::freeBlossom(BlossomId id)
{
  assert(id < mBlossoms.size() && mBlossoms[id].inUse);
  mBlossoms[id].inUse = false;
  mFreeBlossoms.push_back(id);
  --mAliveBlossoms;
}

inline BlossomMerging::EdgeHandle BlossomMerging::allocEdge()
{
  EdgeHandle e;
  if (!mFreeEdges.empty()) {
    e = mFreeEdges.back();
    mFreeEdges.pop_back();
  } else {
    Edge ed;
    ed.inUse = false;
    ed.generation = 0;
    mEdges.push_back(ed);
    e = mEdges.size() - 1;
  }
  mEdges[e].inUse = true;
#ifdef MBX_SHARED_D3_TIE_ORDER
  // A reused slot must never retain the previous original edge's witness.
  mEdges[e].tieRank = cNullItm;
#endif
  mEdges[e].phChild = cNullItm; mEdges[e].phNext = cNullItm; mEdges[e].phPrev = cNullItm;
  ++mAliveEdges;
  if (mAliveEdges > bmPeakAliveEdges()) bmPeakAliveEdges() = mAliveEdges;
  return e;
}

inline void BlossomMerging::freeEdge(EdgeHandle e)
{
  assert(e < mEdges.size() && mEdges[e].inUse);
#ifdef MBX_BM_WATCH
  if (bmWatchA() >= 0) {
    const Edge& E = mEdges[e];
    if (((long)E.v == bmWatchA() && (long)E.w == bmWatchB()) ||
        ((long)E.v == bmWatchB() && (long)E.w == bmWatchA()))
      std::fprintf(stderr, "W free h=%lu cost=%.0f inPkt=%d\n",
                   (unsigned long)e, (double)E.cost, (int)E.inPacket);
  }
#endif
  mEdges[e].inUse = false;
  ++mEdges[e].generation; // invalidate any outstanding (handle, gen) refs
  mFreeEdges.push_back(e);
  --mAliveEdges;
}

inline Size BlossomMerging::allocPacket()
{
  Size p;
  if (!mFreePackets.empty()) {
    p = mFreePackets.back();
    mFreePackets.pop_back();
  } else {
    Packet pk;
    pk.inUse = false;
    mPackets.push_back(pk);
    p = mPackets.size() - 1;
  }
  Packet& P = mPackets[p];
  P.rank = 0;
  P.blossom = cNullItm;
  P.smallestEdge = cNullItm;
  P.nextInBlossom = cNullItm;
  P.count = 0;
  P.inUse = true;
  return p;
}

inline void BlossomMerging::freePacket(Size p)
{
  assert(p < mPackets.size() && mPackets[p].inUse);
  mPackets[p].inUse = false;
  mPackets[p].smallestEdge = cNullItm;
  mPackets[p].count = 0;
  mFreePackets.push_back(p);
}

inline BlossomMerging::BlossomId
BlossomMerging::make_blossom(const std::vector<VertexId>& vertices)
{
  return make_blossom_impl(vertices, cNullItm);
}

inline BlossomMerging::BlossomId
BlossomMerging::make_blossom_impl(const std::vector<VertexId>& vertices,
                                   BlossomId supersededOwner)
{
  assert(!vertices.empty());
  BlossomId id = allocBlossom();
  Blossom& B = mBlossoms[id];
  B.size = vertices.size();
  B.rank = logFloor(B.size);
  VertexId minV = vertices[0];
  Size minD = mDepth[minV];
  // Build vertex list as doubly-linked circular ring.
  VertexId prev = vertices.back();
  for (Size i = 0; i < vertices.size(); ++i) {
    VertexId v = vertices[i];
    assert(v < mN);
    assert(mVtxBlossom[v] == supersededOwner); // null for ordinary make; old owner for trusted replacement
    mVtxBlossom[v] = id;
    mVtxPrev[v] = prev;
    mVtxNext[prev] = v;
    prev = v;
    if (mDepth[v] < minD) {
      minD = mDepth[v];
      minV = v;
    }
  }
  B.rootVertex = minV;
  B.vtxListHead = vertices[0];
  // Initially no edges. No H insert until first edge.
  return id;
}

// Section 3.1, Table 1 classification.
// Inputs: edge e (with v, w, rEdge already set), Bv = blossom containing v,
//         Bw = blossom containing w. Assumes Bv != Bw.
// Outputs: type, chargedEnd u, packet rank r.
inline void
BlossomMerging::classifyEdgeFromBlossoms(EdgeHandle eh, BlossomId Bv,
                                         BlossomId Bw, EdgeType& outType,
                                         VertexId& outU, Size& outR)
{
  const Edge& e = mEdges[eh];
  Size rEdge = e.rEdge;
  Size rBv = mBlossoms[Bv].rank;
  Size rBw = mBlossoms[Bw].rank;
  // Table 1:
  //   l: r(e) > r(Bv)            -> u = v, r = r(e)
  //   s: r(e) <= r(Bv), r(Bv) <= r(Bw) -> u = v, r = max(r(Bv)+1, r(Bw))
  //   d: Bv != Bw, r(Bv) > r(Bw) -> u = w, r = r(Bv)
  if (rEdge > rBv) {
    outType = TYPE_L;
    outU = e.v;
    outR = rEdge;
  } else if (rBv <= rBw) {
    outType = TYPE_S;
    Size rs = rBv + 1;
    if (rBw > rs) rs = rBw;
    outR = rs;
    outU = e.v;
  } else {
    // rBv > rBw
    outType = TYPE_D;
    outU = e.w;
    outR = rBv;
  }
}

inline void BlossomMerging::refreshHEntry(BlossomId b)
{
  Blossom& B = mBlossoms[b];
  if (!B.active) {
    if (B.hHandle != cNullItm) {
      mH.remove(B.hHandle);
      B.hHandle = cNullItm;
    }
    return;
  }
  if (B.smallestOverall == cNullItm) {
    // No edges. Remove from H if present.
    if (B.hHandle != cNullItm) {
      mH.remove(B.hHandle);
      B.hHandle = cNullItm;
    }
    return;
  }
  MinKey newCost = blossomKey(B);
  if (B.hHandle == cNullItm) {
    B.hHandle = mH.insert(newCost, b);
  } else {
    MinKey curKey = mH.key(B.hHandle);
    if (newCost < curKey) {
      mH.decreaseKey(B.hHandle, newCost);
    } else if (curKey < newCost) {
      // F-heap doesn't support increase-key directly; remove and re-insert.
      mH.remove(B.hHandle);
      B.hHandle = mH.insert(newCost, b);
    }
    // equal: no change
  }
}

inline void BlossomMerging::set_vertex_blossom_active(VertexId v, bool active)
{
  if (v >= mVtxBlossom.size()) return;
  BlossomId b = mVtxBlossom[v];
  if (b == cNullItm || b >= mBlossoms.size() || !mBlossoms[b].inUse) return;
  Blossom& B = mBlossoms[b];
  if (B.active == active) return;
  B.active = active;
  refreshHEntry(b);
}

inline void BlossomMerging::addToLoose(BlossomId b, EdgeHandle e)
{
  Blossom& B = mBlossoms[b];
  Edge& E = mEdges[e];
  E.inPacket = false;
  E.phChild = cNullItm; E.phNext = cNullItm; E.phPrev = cNullItm;
  // Section 3.1 inserts make_edge records into loose(B) in O(1). Use the
  // packet pairing-heap primitive here too, retaining O(1) minimum and safe
  // arbitrary lazy deletion without scanning a long loose ring.
  B.looseHead = phMeld(B.looseHead, e);
  B.looseSmallest = B.looseHead;
  B.looseSmallestCost = mEdges[B.looseHead].cost;
  ++B.looseCount;
  if (B.looseCount > bmMaxLooseLen()) bmMaxLooseLen() = B.looseCount;
  ++bmLooseLive();
  if (bmLooseLive() > bmLooseLivePeak()) bmLooseLivePeak() = bmLooseLive();
  if (B.smallestOverall == cNullItm || edgeLess(e, B.smallestOverall)) {
    B.smallestOverall = e;
    B.smallestCost = E.cost;
    refreshHEntry(b);
  }
}

inline void BlossomMerging::addToPacket(BlossomId b, Size rank, EdgeHandle e)
{
  Blossom& B = mBlossoms[b];
  // Find or create packet of given rank.
  Size pkIdx = B.firstPacket;
  Size found = cNullItm;
  while (pkIdx != cNullItm) {
    if (mPackets[pkIdx].rank == rank) { found = pkIdx; break; }
    pkIdx = mPackets[pkIdx].nextInBlossom;
  }
  if (found == cNullItm) {
    found = allocPacket();
    Packet& P = mPackets[found];
    P.rank = rank;
    P.blossom = b;
    P.nextInBlossom = B.firstPacket;
    B.firstPacket = found;
  }
  Packet& P = mPackets[found];
  if (P.count + 1 > bmMaxPacketLen()) bmMaxPacketLen() = P.count + 1;  // diag wire
  Edge& E = mEdges[e];
  E.inPacket = true;
  packetInsert(P, e);                     // pairing meld: O(1)
  if (B.smallestOverall == cNullItm || edgeLess(e, B.smallestOverall)) {
    B.smallestOverall = e;
    B.smallestCost = E.cost;
    refreshHEntry(b);
  }
}

// ---- Per-packet intrusive PAIRING HEAP (by edge cost) ----
// Restores the paper's O(1) packet concatenation (Section 3.1: packets concatenate by
// header at merge, never touching edges) while supporting the eviction extension the
// continuous-loop integration requires: meld O(1), arbitrary delete amortized O(log),
// find-min O(1) (root == Packet::smallestEdge).

inline BlossomMerging::EdgeHandle BlossomMerging::phMeld(EdgeHandle a, EdgeHandle b)
{
  if (a == cNullItm) return b;
  if (b == cNullItm) return a;
  if (edgeLess(b, a)) { EdgeHandle t = a; a = b; b = t; }
  // attach b as a's first child
  Edge& A = mEdges[a]; Edge& B = mEdges[b];
  B.phPrev = a;                       // leftmost child: phPrev = parent
  B.phNext = A.phChild;
  if (A.phChild != cNullItm) mEdges[A.phChild].phPrev = b;
  A.phChild = b;
  return a;
}

inline void BlossomMerging::phDetachNonRoot(EdgeHandle e)
{
  Edge& E = mEdges[e];
  EdgeHandle p = E.phPrev;
  if (mEdges[p].phChild == e) {       // e is the leftmost child of p
    mEdges[p].phChild = E.phNext;
  } else {                            // p is e's left sibling
    mEdges[p].phNext = E.phNext;
  }
  if (E.phNext != cNullItm) mEdges[E.phNext].phPrev = p;
  E.phPrev = cNullItm; E.phNext = cNullItm;
}

inline BlossomMerging::EdgeHandle BlossomMerging::phCombineChildren(EdgeHandle first)
{
  if (first == cNullItm) return cNullItm;
  // pass 1: meld adjacent pairs left-to-right (collect into scratch)
  mSweepScratch.clear();
  EdgeHandle cur = first;
  while (cur != cNullItm) {
    EdgeHandle a = cur, b = mEdges[a].phNext;
    EdgeHandle nxt = (b == cNullItm) ? cNullItm : mEdges[b].phNext;
    mEdges[a].phPrev = cNullItm; mEdges[a].phNext = cNullItm;
    if (b != cNullItm) { mEdges[b].phPrev = cNullItm; mEdges[b].phNext = cNullItm; }
    mSweepScratch.push_back(phMeld(a, b));
    cur = nxt;
  }
  // pass 2: fold right-to-left
  EdgeHandle root = mSweepScratch.back();
  for (Size i = mSweepScratch.size() - 1; i-- > 0; )
    root = phMeld(mSweepScratch[i], root);
  return root;
}

inline void BlossomMerging::packetInsert(Packet& P, EdgeHandle e)
{
  mEdges[e].phChild = cNullItm; mEdges[e].phNext = cNullItm; mEdges[e].phPrev = cNullItm;
  P.smallestEdge = phMeld(P.smallestEdge, e);
  ++P.count;
}

inline void BlossomMerging::packetRemove(Packet& P, EdgeHandle e)
{
  Edge& E = mEdges[e];
  if (P.smallestEdge == e) {
    P.smallestEdge = phCombineChildren(E.phChild);
  } else {
    phDetachNonRoot(e);
    EdgeHandle sub = phCombineChildren(E.phChild);
    P.smallestEdge = phMeld(P.smallestEdge, sub);
  }
  E.phChild = cNullItm; E.phNext = cNullItm; E.phPrev = cNullItm;
  --P.count;
}

inline void BlossomMerging::packetCollect(EdgeHandle root, std::vector<EdgeHandle>& out)
{
  if (root == cNullItm) return;
  Size stackFrom = out.size();
  out.push_back(root);
  Size head = stackFrom;              // treat out[stackFrom..] as a DFS worklist in place
  while (head < out.size()) {
    EdgeHandle x = out[head++];
    for (EdgeHandle c = mEdges[x].phChild; c != cNullItm; c = mEdges[c].phNext)
      out.push_back(c);
  }
}

inline Size BlossomMerging::findPacketOf(const Edge& E)
{
  BlossomId b = mVtxBlossom[E.chargedEnd];
  Size pk = mBlossoms[b].firstPacket;
  while (pk != cNullItm && mPackets[pk].rank != E.packetRank)
    pk = mPackets[pk].nextInBlossom;
  return pk;                          // O(#packets) = O(log n)
}

// Recompute the cached loose-list smallest by a full scan. Called only when
// the previously-cached loose smallest was removed (lazy invalidation).
inline void BlossomMerging::recomputeLooseSmallest(BlossomId b)
{
  Blossom& B = mBlossoms[b];
  // looseHead is the pairing-heap root and therefore already its minimum.
  B.looseSmallest = B.looseHead;
  B.looseSmallestCost = (B.looseHead == cNullItm) ? 0 : mEdges[B.looseHead].cost;
}

inline void BlossomMerging::refreshBlossomSmallest(BlossomId b)
{
  Blossom& B = mBlossoms[b];
  EdgeHandle bestE = cNullItm;
  Val bestC = 0;
  // Use the cached loose-list smallest. If it was invalidated (cNullItm but
  // the ring is non-empty), recompute it once.
  if (B.looseSmallest == cNullItm && B.looseHead != cNullItm) {
    recomputeLooseSmallest(b);
  }
  if (B.looseSmallest != cNullItm) {
    bestE = B.looseSmallest;
    bestC = B.looseSmallestCost;
  }
  // Scan packet smallest values.
  ++bmRefreshCount();
  unsigned long long pktSteps = 0;
  Size pkIdx = B.firstPacket;
  while (pkIdx != cNullItm) {
    ++pktSteps;
    const Packet& P = mPackets[pkIdx];
    if (P.smallestEdge != cNullItm) {
      Val c = mEdges[P.smallestEdge].cost;
      if (bestE == cNullItm || edgeLess(P.smallestEdge, bestE)) {
        bestE = P.smallestEdge;
        bestC = c;
      }
    }
    pkIdx = P.nextInBlossom;
  }
  bmRefreshPacketSteps() += pktSteps;
  if (pktSteps > bmMaxBlossomPackets()) bmMaxBlossomPackets() = pktSteps;
  B.smallestOverall = bestE;
  B.smallestCost = bestC;
  refreshHEntry(b);
}

inline BlossomMerging::EdgeHandle
BlossomMerging::make_edge(VertexId v, VertexId w, Val t, Size payload,
                          Size occurrenceTag, Size tieRank)
{
#ifdef MBX_SHARED_D3_TIE_ORDER
  // Validate before any edge allocation or container mutation; finite costs
  // are necessary for the strict weak order (in particular reject NaN).
  if (!std::isfinite(t))
    throw std::invalid_argument("shared D3 BM requires a finite stored cost");
  if (tieRank == cNullItm)
    throw std::invalid_argument("shared D3 BM requires an explicit original-edge rank");
#else
  (void)tieRank;
#endif
  assert(v < mN && w < mN);
  // Normalize so v is deeper than w (back-edge orientation).
  if (mDepth[v] < mDepth[w]) {
    VertexId tmp = v; v = w; w = tmp;
  }
  assert(mDepth[v] > mDepth[w]);
  BlossomId Bv = mVtxBlossom[v];
  BlossomId Bw = mVtxBlossom[w];
  assert(Bv != cNullItm && Bw != cNullItm);
  if (Bv == Bw) return cNullItm; // discard per "If Bv = Bw then e is discarded"
  EdgeHandle eh = allocEdge();
  Edge& E = mEdges[eh];
  E.v = v;
  E.w = w;
  E.cost = t;
  E.payload = payload;
#ifdef MBX_SHARED_D3_TIE_ORDER
  E.tieRank = tieRank;
#endif
  E.occurrenceTag = occurrenceTag;
  E.rEdge = logFloor(mDepth[v] - mDepth[w]);
  EdgeType ty; VertexId u; Size r;
  classifyEdgeFromBlossoms(eh, Bv, Bw, ty, u, r);
  E.type = ty;
  E.chargedEnd = u;
  E.packetRank = r;
  BlossomId Bu = (u == v) ? Bv : Bw;
  // Faithful Section 3.1 placement: new edges enter loose(Bu) in O(1) and
  // migrate to rank packets when a merge gathers/reclassifies them.
  addToLoose(Bu, eh);
#ifdef MBX_BM_WATCH
  if (bmWatchA() >= 0 &&
      (((long)v == bmWatchA() && (long)w == bmWatchB()) ||
       ((long)v == bmWatchB() && (long)w == bmWatchA())))
    std::fprintf(stderr,
      "W mk h=%lu chg=%lu cost=%.0f B=%lu act=%d Bsm=%ld BsmCost=%.0f inH=%d\n",
      (unsigned long)eh, (unsigned long)u, (double)t, (unsigned long)Bu,
      (int)mBlossoms[Bu].active, (long)mBlossoms[Bu].smallestOverall,
      (double)mBlossoms[Bu].smallestCost, (int)(mBlossoms[Bu].hHandle != cNullItm));
#endif
  return eh;
}

inline void BlossomMerging::remove_edge(EdgeHandle h)
{
  assert(!mEdgeRemovalBatchActive);
  removeEdgeImpl(h, false);
}

inline void BlossomMerging::removeEdgeImpl(EdgeHandle h, bool deferRefresh)
{
  if (h == cNullItm) return;
  if (h >= mEdges.size() || !mEdges[h].inUse) return;
  Edge& E = mEdges[h];
#ifdef MBX_BM_WATCH
  if (bmWatchA() >= 0 &&
      (((long)E.v == bmWatchA() && (long)E.w == bmWatchB()) ||
       ((long)E.v == bmWatchB() && (long)E.w == bmWatchA())))
    std::fprintf(stderr, "W rm h=%lu\n", (unsigned long)h);
#endif
  BlossomId b = mVtxBlossom[E.chargedEnd];   // owner DERIVED (invariant: chargedEnd in owner)
  bool wasBlossomSmallest = (mBlossoms[b].smallestOverall == h);
  if (E.inPacket) {
    Size pk = findPacketOf(E);               // O(log n) rank scan of b's packet list
    assert(pk != cNullItm);
    Packet& P = mPackets[pk];
    packetRemove(P, h);                      // pairing delete: amortized O(log)
  } else {
    Blossom& B = mBlossoms[b];
    if (B.looseHead == h) {
      B.looseHead = phCombineChildren(E.phChild);
    } else {
      phDetachNonRoot(h);
      EdgeHandle sub = phCombineChildren(E.phChild);
      B.looseHead = phMeld(B.looseHead, sub);
    }
    E.phChild = cNullItm; E.phNext = cNullItm; E.phPrev = cNullItm;
    B.looseSmallest = B.looseHead;
    B.looseSmallestCost = (B.looseHead == cNullItm) ? 0 : mEdges[B.looseHead].cost;
    --B.looseCount;
    --bmLooseLive();
  }
  freeEdge(h);
  // Only recompute the blossom's overall-smallest if the edge we removed WAS
  // that smallest. Removing any other edge cannot change the blossom minimum.
  if (wasBlossomSmallest) {
    if (deferRefresh) mDeferredRefreshBlossoms.push_back(b);
    else refreshBlossomSmallest(b);
  }
}

inline void BlossomMerging::begin_edge_removal_batch()
{
  assert(!mEdgeRemovalBatchActive);
  mDeferredRefreshBlossoms.clear();
  mEdgeRemovalBatchActive = true;
}

inline bool BlossomMerging::remove_edge_if_batched(EdgeHandle h, Size gen)
{
  assert(mEdgeRemovalBatchActive);
  if (h == cNullItm || h >= mEdges.size()) return false;
  if (!mEdges[h].inUse || mEdges[h].generation != gen) return false;
  removeEdgeImpl(h, true);
  return true;
}

inline void BlossomMerging::end_edge_removal_batch()
{
  assert(mEdgeRemovalBatchActive);
  // Preserve the same first-touch order as immediate removal.  An owner can
  // enter this list only once: after its cached minimum is removed, that
  // cached handle deliberately remains unchanged until this loop, so no later
  // removal from the same owner compares equal to smallestOverall.
  for (Size i = 0; i < mDeferredRefreshBlossoms.size(); ++i) {
    const BlossomId b = mDeferredRefreshBlossoms[i];
    if (b < mBlossoms.size() && mBlossoms[b].inUse)
      refreshBlossomSmallest(b);
  }
  mDeferredRefreshBlossoms.clear();
  mEdgeRemovalBatchActive = false;
}

inline void BlossomMerging::dispose_blossom(BlossomId Z,
                                            std::vector<Size>* outDiscardedPayloads)
{
  dispose_blossom_impl(Z, outDiscardedPayloads, true);
}

inline void BlossomMerging::dispose_blossom_impl(BlossomId Z,
    std::vector<Size>* outDiscardedPayloads, bool clearVertexRing)
{
  if (Z == cNullItm) return;
  if (Z >= mBlossoms.size() || !mBlossoms[Z].inUse) return;
  Blossom& B = mBlossoms[Z];

  // Remove Z from H.
  if (B.hHandle != cNullItm) {
    mH.remove(B.hHandle);
    B.hHandle = cNullItm;
  }

  // Free all edges in Z's loose pairing heap (recording payloads for refresh).
  {
    if (B.looseHead != cNullItm) {
      std::vector<EdgeHandle> looseEdges;
      packetCollect(B.looseHead, looseEdges);
      for (Size i = 0; i < looseEdges.size(); ++i) {
        EdgeHandle e = looseEdges[i];
        if (outDiscardedPayloads != NULL)
          outDiscardedPayloads->push_back(mEdges[e].payload);
        freeEdge(e);
      }
    }
    B.looseHead = cNullItm;
    B.looseSmallest = cNullItm;
    bmLooseLive() -= B.looseCount;
    B.looseCount = 0;
  }

  // Free all edges in each packet (collected from its pairing heap), then free each packet.
  {
    std::vector<EdgeHandle> hp;
    Size pk = B.firstPacket;
    while (pk != cNullItm) {
      Size nxt = mPackets[pk].nextInBlossom;
      hp.clear();
      packetCollect(mPackets[pk].smallestEdge, hp);
      for (Size hi = 0; hi < hp.size(); ++hi) {
        if (outDiscardedPayloads != NULL)
          outDiscardedPayloads->push_back(mEdges[hp[hi]].payload);
        freeEdge(hp[hi]);
      }
      freePacket(pk);
      pk = nxt;
    }
    B.firstPacket = cNullItm;
  }

  // Only the private complete-frontier replacement may skip these dead
  // stores. Ordinary disposal still clears every vertex/link exactly as before.
  if (B.vtxListHead != cNullItm) {
    if (clearVertexRing) {
      VertexId head = B.vtxListHead;
      VertexId cur = head;
      do {
        VertexId nxt = mVtxNext[cur];
        mVtxBlossom[cur] = cNullItm;
        mVtxNext[cur] = cNullItm;
        mVtxPrev[cur] = cNullItm;
        cur = nxt;
      } while (cur != head);
    }
    B.vtxListHead = cNullItm;
  }

  B.smallestOverall = cNullItm;
  B.size = 0;
  B.rank = 0;
  freeBlossom(Z);
}


inline void BlossomMerging::begin_largest_child_partition(
    BlossomId Z, std::vector<Size>* outDiscardedPayloads)
{
  if (!blossom_alive(Z) || mEdgeRemovalBatchActive)
    throw std::runtime_error("largest partition invalid parent/batch");
  Blossom& B = mBlossoms[Z];
  // Remove Z from H.
  if (B.hHandle != cNullItm) {
    mH.remove(B.hHandle);
    B.hHandle = cNullItm;
  }

  // Free all edges in Z's loose pairing heap (recording payloads for refresh).
  {
    if (B.looseHead != cNullItm) {
      std::vector<EdgeHandle> looseEdges;
      packetCollect(B.looseHead, looseEdges);
      for (Size i = 0; i < looseEdges.size(); ++i) {
        EdgeHandle e = looseEdges[i];
        if (outDiscardedPayloads != NULL)
          outDiscardedPayloads->push_back(mEdges[e].payload);
        freeEdge(e);
      }
    }
    B.looseHead = cNullItm;
    B.looseSmallest = cNullItm;
    bmLooseLive() -= B.looseCount;
    B.looseCount = 0;
  }

  // Free all edges in each packet (collected from its pairing heap), then free each packet.
  {
    std::vector<EdgeHandle> hp;
    Size pk = B.firstPacket;
    while (pk != cNullItm) {
      Size nxt = mPackets[pk].nextInBlossom;
      hp.clear();
      packetCollect(mPackets[pk].smallestEdge, hp);
      for (Size hi = 0; hi < hp.size(); ++hi) {
        if (outDiscardedPayloads != NULL)
          outDiscardedPayloads->push_back(mEdges[hp[hi]].payload);
        freeEdge(hp[hi]);
      }
      freePacket(pk);
      pk = nxt;
    }
    B.firstPacket = cNullItm;
  }


  // Existing edge destruction order/generations are preserved above. Keep
  // owner/ring/header alive, but expose no stale edge/H/minimum state.
  B.looseSmallestCost = 0;
  B.smallestOverall = cNullItm;
  B.smallestCost = 0;
}

template<class Members>
inline BlossomMerging::BlossomId BlossomMerging::make_smaller_for_largest_partition(
    BlossomId Z, const Members& vertices)
{
  if (!blossom_alive(Z) || vertices.empty() || vertices.size() >= mBlossoms[Z].size)
    throw std::runtime_error("largest partition invalid smaller child");
  // Allocate before modifying vertices. Do not retain Blossom& across this:
  // vector growth may move every header, including the parent.
  const BlossomId id = allocBlossom();
  for (Size i=0; i<vertices.size(); ++i) {
    const VertexId v=vertices[i];
    if (v>=mN || mVtxBlossom[v]!=Z)
      throw std::runtime_error("largest partition foreign/repeated atom");
    const VertexId prev=mVtxPrev[v], next=mVtxNext[v];
    if (prev>=mN || next>=mN || mVtxNext[prev]!=v || mVtxPrev[next]!=v)
      throw std::runtime_error("largest partition old ring links");
    mVtxNext[prev]=next;
    mVtxPrev[next]=prev;
    if (mBlossoms[Z].vtxListHead==v) mBlossoms[Z].vtxListHead=next;
    mVtxBlossom[v]=cNullItm;
    mVtxPrev[v]=mVtxNext[v]=cNullItm;
    --mBlossoms[Z].size;
  }
  // Only after every member is detached may its links be repurposed.
  Blossom& C=mBlossoms[id];
  C.size=vertices.size(); C.rank=logFloor(C.size);
  VertexId root=vertices[0], prev=vertices[vertices.size()-1];
  Size minimum=mDepth[root];
  for (Size i=0; i<vertices.size(); ++i) {
    const VertexId v=vertices[i];
    mVtxBlossom[v]=id; mVtxPrev[v]=prev; mVtxNext[prev]=v; prev=v;
    if (mDepth[v]<minimum) {minimum=mDepth[v];root=v;}
  }
  C.rootVertex=root; C.vtxListHead=vertices[0];
  return id;
}

inline void BlossomMerging::finish_largest_child_partition(
    BlossomId Z, VertexId canonicalHead, Size remaining, VertexId minimumRoot)
{
  if (!blossom_alive(Z) || remaining==0 || mBlossoms[Z].size!=remaining ||
      canonicalHead>=mN || minimumRoot>=mN || mVtxBlossom[canonicalHead]!=Z ||
      mVtxBlossom[minimumRoot]!=Z || mBlossoms[Z].firstPacket!=cNullItm ||
      mBlossoms[Z].looseHead!=cNullItm || mBlossoms[Z].hHandle!=cNullItm)
    throw std::runtime_error("largest partition incomplete coverage/cleanup");
  Blossom& B=mBlossoms[Z];
  B.rank=logFloor(remaining); B.rootVertex=minimumRoot;
  B.vtxListHead=canonicalHead; B.active=true;
}

inline bool BlossomMerging::min_edge(VertexId& outV, VertexId& outW, Val& outCost,
                                     Size* outPayload, EdgeHandle* outHandle)
{
  // F-heap H has lazy deletions: a popped entry is valid only if its
  // payload blossom is still the live current blossom owning the
  // edge. Since we maintain H via insert/decrease/remove on every
  // smallestEdge change (refreshHEntry), the top entry is always live.
  // However, the blossom's smallestEdge might have its endpoints in
  // the same blossom after a merge (because R0 reclassification may
  // not catch every stale edge in non-reclassified packets). We do a
  // validity check on the returned edge.
  while (!mH.empty()) {
    BlossomId b = mH.minPayload();
    assert(mBlossoms[b].inUse);
    EdgeHandle e = mBlossoms[b].smallestOverall;
    if (e == cNullItm) {
      // Stale: blossom has no edges. Remove from H.
      mH.remove(mBlossoms[b].hHandle);
      mBlossoms[b].hHandle = cNullItm;
      continue;
    }
    const Edge& E = mEdges[e];
    BlossomId Bv = mVtxBlossom[E.v];
    BlossomId Bw = mVtxBlossom[E.w];
    if (Bv == Bw) {
      // Stale: edge now joins the same blossom. With the per-packet heap,
      // removing this one edge is O(log packet) and the next-smallest surfaces
      // automatically — no bulk sweep / ring rescan needed. (remove_edge also
      // calls refreshBlossomSmallest if it was the blossom min.)
      ++bmMinObservedInternalStale();   // observed-at-cleanup only; see counter comment
      remove_edge(e);
      continue;
    }
    outV = E.v;
    outW = E.w;
    outCost = E.cost;
    if (outPayload != NULL) *outPayload = E.payload;
    if (outHandle != NULL) *outHandle = e;
    return true;
  }
  return false;
}

inline BlossomMerging::BlossomId BlossomMerging::merge(BlossomId X, BlossomId Y,
    std::vector<Size>* outDiscardedPayloads)
{
  assert(X != Y);
  assert(X < mBlossoms.size() && Y < mBlossoms.size());
  assert(mBlossoms[X].inUse && mBlossoms[Y].inUse);

  // Direct vertex->blossom lookup is deliberately O(1), so merge must rewrite
  // one side's vertex ring. Always rewrite the smaller side: every vertex then
  // changes representative at most O(log n) times instead of the former
  // adversarial Theta(n^2) sequence.
  if (mBlossoms[X].size < mBlossoms[Y].size) {
    BlossomId t = X; X = Y; Y = t;
  }

  // Form Z by absorbing Y into X (reuse X's slot as Z to avoid extra alloc).
  // First compute Z's size and rank.
  Size newSize = mBlossoms[X].size + mBlossoms[Y].size;
  Size newRank = logFloor(newSize);
  bool newActive = mBlossoms[X].active || mBlossoms[Y].active;
  BlossomId Z = X;

  // Pick Z's root as min-depth between X's root and Y's root.
  VertexId rootZ = mBlossoms[X].rootVertex;
  if (mDepth[mBlossoms[Y].rootVertex] < mDepth[rootZ]) {
    rootZ = mBlossoms[Y].rootVertex;
  }

  // Remove X and Y from H now; we'll re-insert Z after building it.
  if (mBlossoms[X].hHandle != cNullItm) {
    mH.remove(mBlossoms[X].hHandle);
    mBlossoms[X].hHandle = cNullItm;
  }
  if (mBlossoms[Y].hHandle != cNullItm) {
    mH.remove(mBlossoms[Y].hHandle);
    mBlossoms[Y].hHandle = cNullItm;
  }

  // Collect R0 = edges that must be reclassified:
  //   - edges in packets of X or Y with rank <= newRank
  //   - all edges in loose(X) and loose(Y)
  // We extract them into a temporary vector, freeing them from their
  // containers as we go.
  std::vector<EdgeHandle> r0;
  r0.reserve(8);

  // Extract the two loose pairing heaps. (No C++11 helper lambda.)
  {
    packetCollect(mBlossoms[X].looseHead, r0);
    mBlossoms[X].looseHead = cNullItm;
    mBlossoms[X].looseSmallest = cNullItm;
    bmLooseLive() -= mBlossoms[X].looseCount;
    mBlossoms[X].looseCount = 0;
  }
  {
    packetCollect(mBlossoms[Y].looseHead, r0);
    mBlossoms[Y].looseHead = cNullItm;
    mBlossoms[Y].looseSmallest = cNullItm;
    bmLooseLive() -= mBlossoms[Y].looseCount;
    mBlossoms[Y].looseCount = 0;
  }

  // Also extract packet contents of rank <= newRank from X and Y.
  // For X:
  {
    Size pkPrev = cNullItm;
    Size pkIdx = mBlossoms[X].firstPacket;
    while (pkIdx != cNullItm) {
      Size nxt = mPackets[pkIdx].nextInBlossom;
      if (mPackets[pkIdx].rank <= newRank) {
        packetCollect(mPackets[pkIdx].smallestEdge, r0);   // extract for reclassification
        if (pkPrev == cNullItm) mBlossoms[X].firstPacket = nxt;
        else mPackets[pkPrev].nextInBlossom = nxt;
        freePacket(pkIdx);
      } else {
        pkPrev = pkIdx;
      }
      pkIdx = nxt;
    }
  }
  // For Y:
  {
    Size pkPrev = cNullItm;
    Size pkIdx = mBlossoms[Y].firstPacket;
    while (pkIdx != cNullItm) {
      Size nxt = mPackets[pkIdx].nextInBlossom;
      if (mPackets[pkIdx].rank <= newRank) {
        packetCollect(mPackets[pkIdx].smallestEdge, r0);
        if (pkPrev == cNullItm) mBlossoms[Y].firstPacket = nxt;
        else mPackets[pkPrev].nextInBlossom = nxt;
        freePacket(pkIdx);
      } else {
        pkPrev = pkIdx;
      }
      pkIdx = nxt;
    }
  }

  // Update vertex -> blossom mapping for all Y vertices to point at Z.
  // Walk Y's vertex ring in O(|Y|). Splice Y's ring into Z's ring.
  if (mBlossoms[Y].vtxListHead != cNullItm) {
    VertexId yHead = mBlossoms[Y].vtxListHead;
    VertexId yCur = yHead;
    do {
      mVtxBlossom[yCur] = Z;
      yCur = mVtxNext[yCur];
    } while (yCur != yHead);
    // Splice Y's circular ring into Z's circular ring.
    if (mBlossoms[Z].vtxListHead == cNullItm) {
      mBlossoms[Z].vtxListHead = yHead;
    } else {
      VertexId zHead = mBlossoms[Z].vtxListHead;
      VertexId zTail = mVtxPrev[zHead];
      VertexId yTail = mVtxPrev[yHead];
      mVtxNext[zTail] = yHead;
      mVtxPrev[yHead] = zTail;
      mVtxNext[yTail] = zHead;
      mVtxPrev[zHead] = yTail;
    }
  }

  // Migrate Y's surviving (high-rank) packets to Z. Walk Y's remaining packet
  // list; for each packet of rank > newRank, either move it to Z (if Z has no
  // packet of that rank) or concatenate its ring into Z's packet (if Z does).
  {
    Size pkY = mBlossoms[Y].firstPacket;
    while (pkY != cNullItm) {
      Size nxt = mPackets[pkY].nextInBlossom;
      Size r = mPackets[pkY].rank;
      // r > newRank guaranteed since rank<=newRank were extracted earlier.
      assert(r > newRank);
      // Find Z's packet of rank r (Z = X's slot, so check Z's surviving packet list).
      Size pkZ = mBlossoms[Z].firstPacket;
      Size foundZ = cNullItm;
      while (pkZ != cNullItm) {
        if (mPackets[pkZ].rank == r) { foundZ = pkZ; break; }
        pkZ = mPackets[pkZ].nextInBlossom;
      }
      // THE PAPER'S O(1) CONCATENATION, RESTORED: no owner rebrand (owner is derived),
      // no per-edge walk in either branch.
      gBMPacketScan += mPackets[pkY].count;  // diag: LOGICAL freight (what a walk would touch)
      if (foundZ == cNullItm) {
        mPackets[pkY].blossom = Z;                          // header move: O(1)
        mPackets[pkY].nextInBlossom = mBlossoms[Z].firstPacket;
        mBlossoms[Z].firstPacket = pkY;
      } else {
        Packet& PZ = mPackets[foundZ];                      // same-rank packets: MELD, O(1)
        PZ.smallestEdge = phMeld(PZ.smallestEdge, mPackets[pkY].smallestEdge);
        PZ.count += mPackets[pkY].count;
        mPackets[pkY].smallestEdge = cNullItm; mPackets[pkY].count = 0;
        freePacket(pkY);
      }
      pkY = nxt;
    }
    mBlossoms[Y].firstPacket = cNullItm;
  }

  // Free Y's slot.
  freeBlossom(Y);

  // Update Z's metadata.
  mBlossoms[Z].size = newSize;
  mBlossoms[Z].rank = newRank;
  mBlossoms[Z].rootVertex = rootZ;
  mBlossoms[Z].active = newActive;

  // Gather: prune R0 to remove edges with both endpoints now in Z, and
  // dedupe (at most one edge per pair (Z, B)) keeping smallest cost.
  // Use mGatherHead indexed by the OTHER blossom id mod something? No —
  // paper uses I[1..n] which is large enough. We use mGatherHead sized to
  // log n for new packets. For the dedupe-by-other-blossom step we need
  // an array sized to max BlossomId, which can be up to n.
  // To keep memory bounded, we'll do a simpler dedupe via a temporary
  // map. For correctness first, optimize later.

  // Prune: drop edges with both endpoints in Z; keep smallest per other-blossom.
  // Group by the OTHER blossom using mDedupBestPerB (reused across calls).
  // Grow on demand to keep memory bounded by max-ever-seen blossom count.
  if (mDedupBestPerB.size() < mBlossoms.size()) {
    mDedupBestPerB.resize(mBlossoms.size(), cNullItm);
  }
  std::vector<BlossomId> touchedB;
  touchedB.reserve(r0.size());

  // Paper Section 3.1's "I[B]" dedup: per other-blossom, keep only the
  // smallest stored-cost edge. For matchbox's lazy-delta scheme (where
  // slacks change over time) the dedup-discarded edges may become valid
  // candidates later — we expose their payloads via outDiscardedPayloads
  // so the engine can clear mD3Queued for those edges, letting matchbox's
  // refreshEvenBlossomEdges re-push them when needed.
  for (Size i = 0; i < r0.size(); ++i) {
    EdgeHandle e = r0[i];
    const Edge& E = mEdges[e];
    BlossomId Bv = mVtxBlossom[E.v];
    BlossomId Bw = mVtxBlossom[E.w];
    if (Bv == Bw) {
      // Both endpoints in same blossom — edge no longer a candidate (INTERNAL).
      ++bmDiscardSameBlsm();
      if (outDiscardedPayloads != NULL)
        outDiscardedPayloads->push_back(E.payload);
      freeEdge(e);
      continue;
    }
    BlossomId other = (Bv == Z) ? Bw : Bv;
    if (other == Z) {
      ++bmDiscardSameBlsm();
      if (outDiscardedPayloads != NULL)
        outDiscardedPayloads->push_back(E.payload);
      freeEdge(e);
      continue;
    }
    EdgeHandle cur = mDedupBestPerB[other];
    if (cur == cNullItm) {
      mDedupBestPerB[other] = e;
      touchedB.push_back(other);
    } else if (edgeLess(e, cur)) {
      ++bmDiscardDedup();
      if (outDiscardedPayloads != NULL)
        outDiscardedPayloads->push_back(mEdges[cur].payload);
      freeEdge(cur);
      mDedupBestPerB[other] = e;
    } else {
      ++bmDiscardDedup();
      if (outDiscardedPayloads != NULL)
        outDiscardedPayloads->push_back(E.payload);
      freeEdge(e);
    }
  }

  // R = the surviving edges. Reclassify each; route per the paper (§3.1):
  //   u in Z      -> gather into packet(Z, r)
  //   u outside Z -> loose(Bu)  (the paper's waiting area; participates in I2 /
  //                  Lemma 3.4 charging. Adopted 2026-07-20 after the identical-
  //                  trace A/B against the historical packet(Bu) deviation.)
  for (Size i = 0; i < touchedB.size(); ++i) {
    BlossomId other = touchedB[i];
    EdgeHandle e = mDedupBestPerB[other];
    if (e == cNullItm) continue;
    Edge& E = mEdges[e];
    BlossomId Bv = mVtxBlossom[E.v];
    BlossomId Bw = mVtxBlossom[E.w];
    EdgeType ty; VertexId u; Size r;
    classifyEdgeFromBlossoms(e, Bv, Bw, ty, u, r);
    ++bmReclassCount();
    E.type = ty;
    E.chargedEnd = u;
    E.packetRank = r;
    BlossomId Bu = (u == E.v) ? Bv : Bw;
    bool wasPacket = E.inPacket;      // container of origin, still valid post-extraction
    if (Bu != Z) {
      ++bmR22ToLoose();
      if (wasPacket) ++bmPacketToLoose();
      addToLoose(Bu, e);
    } else {
      ++bmR22ToPacketZ();
      if (!wasPacket) ++bmLooseToPacket();
      addToPacket(Bu, r, e);
    }
  }
  // Clear touched entries in the scratch buffer for next merge call.
  for (Size i = 0; i < touchedB.size(); ++i) {
    mDedupBestPerB[touchedB[i]] = cNullItm;
  }

  // Refresh Z's smallest (it may have changed from packet additions and
  // from un-reclassified packets of rank > newRank that were kept in place).
  refreshBlossomSmallest(Z);

  return Z;
}

#ifdef MBX_BM_VALIDATE
// ==== EXHAUSTIVE INVARIANT VALIDATOR (compile-gated; see declaration) ====

#define BMV_FAIL(name, fmt, ...) do { \
    if (viol < 50) std::fprintf(stderr, "VIOLATION[" name "] " fmt "\n", __VA_ARGS__); \
    else if (viol == 50) std::fprintf(stderr, "VIOLATION[...] further violations suppressed\n"); \
    ++viol; } while (0)

inline Size BlossomMerging::validateHeapWalk(EdgeHandle root, BlossomId owner,
    bool expectInPacket, Size expectRank, bool paperPure,
    std::vector<unsigned char>& seen, std::vector<BlossomId>& edgeOwner,
    Val& ioMin, EdgeHandle& ioMinE, bool& ioHave,
    Size& viol, Size* staleCls) const
{
  if (root == cNullItm) return 0;
  Size cnt = 0;
  // root shape: a pairing-heap root has no parent and no sibling
  if (root >= mEdges.size()) {
    BMV_FAIL("PH-root-bounds", "owner=%lu root=%lu", (unsigned long)owner, (unsigned long)root);
    return 0;
  }
  if (mEdges[root].phPrev != cNullItm || mEdges[root].phNext != cNullItm)
    BMV_FAIL("PH-root-links", "owner=%lu root=%lu phPrev=%ld phNext=%ld",
             (unsigned long)owner, (unsigned long)root,
             (long)mEdges[root].phPrev, (long)mEdges[root].phNext);
  std::vector<EdgeHandle> stack;
  stack.push_back(root);
  while (!stack.empty()) {
    EdgeHandle e = stack.back(); stack.pop_back();
    if (e >= mEdges.size()) { BMV_FAIL("I1i-bounds", "owner=%lu e=%lu", (unsigned long)owner, (unsigned long)e); continue; }
    const Edge& E = mEdges[e];
    if (!E.inUse) { BMV_FAIL("I1i-freed-reachable", "owner=%lu e=%lu", (unsigned long)owner, (unsigned long)e); continue; }
    if (seen[e]) { BMV_FAIL("I1i-multi-container", "owner=%lu e=%lu (cycle or second container)", (unsigned long)owner, (unsigned long)e); continue; }
    seen[e] = 1; edgeOwner[e] = owner; ++cnt;
    // ---- per-edge CLASS/I1i checks ----
    if (E.inPacket != expectInPacket)
      BMV_FAIL("CLASS-flag", "e=%lu inPacket=%d container=%s", (unsigned long)e, (int)E.inPacket,
               expectInPacket ? "packet" : "loose");
    if (E.chargedEnd != E.v && E.chargedEnd != E.w)
      BMV_FAIL("CLASS-charged", "e=%lu chargedEnd=%lu v=%lu w=%lu", (unsigned long)e,
               (unsigned long)E.chargedEnd, (unsigned long)E.v, (unsigned long)E.w);
    if (E.v >= mN || E.w >= mN) {
      BMV_FAIL("I1i-endpoint-bounds", "e=%lu v=%lu w=%lu n=%lu", (unsigned long)e,
               (unsigned long)E.v, (unsigned long)E.w, (unsigned long)mN);
    } else {
      if (mDepth[E.v] <= mDepth[E.w])
        BMV_FAIL("CLASS-depth", "e=%lu depth(v)=%lu depth(w)=%lu", (unsigned long)e,
                 (unsigned long)mDepth[E.v], (unsigned long)mDepth[E.w]);
      else if (E.rEdge != logFloor(mDepth[E.v] - mDepth[E.w]))
        BMV_FAIL("CLASS-redge", "e=%lu stored=%lu recomputed=%lu", (unsigned long)e,
                 (unsigned long)E.rEdge, (unsigned long)logFloor(mDepth[E.v] - mDepth[E.w]));
      if (E.chargedEnd == E.v || E.chargedEnd == E.w) {
        BlossomId derived = mVtxBlossom[E.chargedEnd];
        if (derived != owner)
          BMV_FAIL("CLASS-owner", "e=%lu derivedOwner=%ld container=%lu", (unsigned long)e,
                   (long)derived, (unsigned long)owner);
      }
      BlossomId Bv = mVtxBlossom[E.v], Bw = mVtxBlossom[E.w];
      bool bvOk = (Bv != cNullItm && Bv < mBlossoms.size() && mBlossoms[Bv].inUse);
      bool bwOk = (Bw != cNullItm && Bw < mBlossoms.size() && mBlossoms[Bw].inUse);
      if (!bvOk || !bwOk)
        BMV_FAIL("I1i-endpoint-blossom", "e=%lu Bv=%ld Bw=%ld", (unsigned long)e, (long)Bv, (long)Bw);
      else {
        if (paperPure && Bv == Bw)
          BMV_FAIL("I1i-internal", "e=%lu both endpoints in blossom %lu cost=%.17g",
                   (unsigned long)e, (unsigned long)Bv, (double)E.cost);
        if (Bv != Bw && staleCls != NULL) {
          // Recompute Table-1 with CURRENT blossoms; a difference is
          // intentionally-stale classification (legal until R0), counted only.
          Size rBv = mBlossoms[Bv].rank, rBw = mBlossoms[Bw].rank;
          EdgeType ty; VertexId u; Size r;
          if (E.rEdge > rBv)      { ty = TYPE_L; u = E.v; r = E.rEdge; }
          else if (rBv <= rBw)    { ty = TYPE_S; u = E.v; r = (rBw > rBv + 1) ? rBw : rBv + 1; }
          else                    { ty = TYPE_D; u = E.w; r = rBv; }
          if (ty != E.type || u != E.chargedEnd || r != E.packetRank) ++*staleCls;
        }
      }
    }
    if (expectInPacket && E.packetRank != expectRank)
      BMV_FAIL("CLASS-rank", "e=%lu packetRank=%lu packet.rank=%lu", (unsigned long)e,
               (unsigned long)E.packetRank, (unsigned long)expectRank);
#ifdef MBX_SHARED_D3_TIE_ORDER
    if (!std::isfinite(E.cost) || E.tieRank == cNullItm)
      BMV_FAIL("KEY-invalid", "e=%lu cost=%.17g rank=%lu", (unsigned long)e,
               (double)E.cost, (unsigned long)E.tieRank);
#endif
    if (!ioHave || edgeLess(e, ioMinE)) { ioMin = E.cost; ioMinE = e; ioHave = true; }
    // ---- child chain: link shape + heap order, then descend ----
    EdgeHandle prevSib = cNullItm;
    for (EdgeHandle c = E.phChild; c != cNullItm; c = mEdges[c].phNext) {
      if (c >= mEdges.size()) { BMV_FAIL("PH-child-bounds", "e=%lu c=%lu", (unsigned long)e, (unsigned long)c); break; }
      EdgeHandle expectPrev = (prevSib == cNullItm) ? e : prevSib;
      if (mEdges[c].phPrev != expectPrev)
        BMV_FAIL("PH-sibling-links", "child=%lu phPrev=%ld expected=%lu", (unsigned long)c,
                 (long)mEdges[c].phPrev, (unsigned long)expectPrev);
      if (mEdges[c].inUse && edgeLess(c, e))
        BMV_FAIL("PH-order", "parent=%lu cost=%.17g child=%lu cost=%.17g", (unsigned long)e,
                 (double)E.cost, (unsigned long)c, (double)mEdges[c].cost);
      stack.push_back(c);
      prevSib = c;
      if (stack.size() > mEdges.size() + 1) { BMV_FAIL("PH-cycle", "owner=%lu", (unsigned long)owner); return cnt; }
    }
  }
  return cnt;
}

inline Size BlossomMerging::validate(bool paperPure, Size* outStaleClassifications) const
{
  Size viol = 0;
  Size staleCnt = 0;
  // ---- RING: blossom vertex rings + vertex coverage ----
  std::vector<unsigned char> vtxSeen(mN, 0);
  Size liveBlossoms = 0;
  for (BlossomId b = 0; b < mBlossoms.size(); ++b) {
    if (!mBlossoms[b].inUse) continue;
    ++liveBlossoms;
    const Blossom& B = mBlossoms[b];
    if (B.vtxListHead == cNullItm) {
      if (B.size != 0)
        BMV_FAIL("RING-size", "b=%lu size=%lu but empty ring", (unsigned long)b, (unsigned long)B.size);
      continue;
    }
    Size cnt = 0, minD = 0; bool haveD = false;
    VertexId cur = B.vtxListHead;
    bool ringOk = true;
    do {
      if (cur >= mN) { BMV_FAIL("RING-bounds", "b=%lu v=%lu", (unsigned long)b, (unsigned long)cur); ringOk = false; break; }
      if (vtxSeen[cur]) { BMV_FAIL("RING-dup", "b=%lu v=%lu seen twice", (unsigned long)b, (unsigned long)cur); ringOk = false; break; }
      vtxSeen[cur] = 1; ++cnt;
      if (mVtxBlossom[cur] != b)
        BMV_FAIL("RING-map", "b=%lu v=%lu mVtxBlossom=%ld", (unsigned long)b, (unsigned long)cur, (long)mVtxBlossom[cur]);
      if (mVtxNext[cur] >= mN || mVtxPrev[mVtxNext[cur]] != cur) {
        BMV_FAIL("RING-links", "b=%lu v=%lu next=%ld", (unsigned long)b, (unsigned long)cur, (long)mVtxNext[cur]);
        ringOk = false; break;
      }
      if (!haveD || mDepth[cur] < minD) { minD = mDepth[cur]; haveD = true; }
      cur = mVtxNext[cur];
      if (cnt > mN) { BMV_FAIL("RING-cycle", "b=%lu ring exceeds n", (unsigned long)b); ringOk = false; break; }
    } while (cur != B.vtxListHead);
    if (ringOk) {
      if (cnt != B.size)
        BMV_FAIL("RING-size", "b=%lu size=%lu ring=%lu", (unsigned long)b, (unsigned long)B.size, (unsigned long)cnt);
      if (B.rank != logFloor(B.size))
        BMV_FAIL("RING-rank", "b=%lu rank=%lu logFloor(size)=%lu", (unsigned long)b,
                 (unsigned long)B.rank, (unsigned long)logFloor(B.size));
      if (B.rootVertex >= mN || mVtxBlossom[B.rootVertex] != b || mDepth[B.rootVertex] != minD)
        BMV_FAIL("RING-root", "b=%lu rootVertex=%ld", (unsigned long)b, (long)B.rootVertex);
    }
  }
  for (VertexId v = 0; v < mN; ++v) {
    BlossomId b = mVtxBlossom[v];
    if (b == cNullItm) {
      if (vtxSeen[v]) BMV_FAIL("RING-orphan", "v=%lu in a ring but unmapped", (unsigned long)v);
      continue;
    }
    if (b >= mBlossoms.size() || !mBlossoms[b].inUse)
      BMV_FAIL("RING-deadmap", "v=%lu mapped to dead blossom %ld", (unsigned long)v, (long)b);
    else if (!vtxSeen[v])
      BMV_FAIL("RING-unreached", "v=%lu mapped to %lu but not in its ring", (unsigned long)v, (unsigned long)b);
  }
  // ---- edges: container traversal (I1i / CLASS / PH / PKT / MIN / per-blossom H) ----
  std::vector<unsigned char> seen(mEdges.size(), 0);
  std::vector<BlossomId> edgeOwner(mEdges.size(), cNullItm);
  std::vector<unsigned char> rankSeen(mLogN + 2, 0);
  Size hHandles = 0;
  bool haveGlobal = false; MinKey globalMin = 0;
  for (BlossomId b = 0; b < mBlossoms.size(); ++b) {
    if (!mBlossoms[b].inUse) continue;
    const Blossom& B = mBlossoms[b];
    Val minC = 0; EdgeHandle minE = cNullItm; bool have = false;
    Size looseWalked = validateHeapWalk(B.looseHead, b, false, 0, paperPure, seen, edgeOwner,
                                        minC, minE, have, viol, &staleCnt);
    if (looseWalked != B.looseCount)
      BMV_FAIL("PKT-loose-count", "b=%lu looseCount=%lu walked=%lu", (unsigned long)b,
               (unsigned long)B.looseCount, (unsigned long)looseWalked);
    // loose cache: looseSmallest is either the root or lazily invalidated
    if (B.looseSmallest != cNullItm) {
      if (B.looseSmallest != B.looseHead)
        BMV_FAIL("MIN-loose", "b=%lu looseSmallest=%ld looseHead=%ld", (unsigned long)b,
                 (long)B.looseSmallest, (long)B.looseHead);
      else if (B.looseSmallestCost != mEdges[B.looseHead].cost)
        BMV_FAIL("MIN-loose-cost", "b=%lu cached=%.17g actual=%.17g", (unsigned long)b,
                 (double)B.looseSmallestCost, (double)mEdges[B.looseHead].cost);
    }
    // packets
    std::vector<Size> touchedRanks;
    Size guard = 0;
    for (Size pk = B.firstPacket; pk != cNullItm; ) {
      if (pk >= mPackets.size()) { BMV_FAIL("PKT-bounds", "b=%lu pk=%lu", (unsigned long)b, (unsigned long)pk); break; }
      const Packet& P = mPackets[pk];
      if (!P.inUse) { BMV_FAIL("PKT-freed-linked", "b=%lu pk=%lu", (unsigned long)b, (unsigned long)pk); break; }
      if (P.blossom != b)
        BMV_FAIL("PKT-owner", "pk=%lu header blossom=%ld list owner=%lu", (unsigned long)pk,
                 (long)P.blossom, (unsigned long)b);
      if (P.rank <= B.rank || P.rank > mLogN)
        BMV_FAIL("PKT-rank-range", "b=%lu rank(B)=%lu pk=%lu rank=%lu logN=%lu", (unsigned long)b,
                 (unsigned long)B.rank, (unsigned long)pk, (unsigned long)P.rank, (unsigned long)mLogN);
      if (P.rank < rankSeen.size()) {
        if (P.count > 0 && rankSeen[P.rank])
          BMV_FAIL("PKT-dup-rank", "b=%lu rank=%lu has two nonempty headers", (unsigned long)b, (unsigned long)P.rank);
        if (P.count > 0) { rankSeen[P.rank] = 1; touchedRanks.push_back(P.rank); }
      }
      if (P.count == 0) {
        if (P.smallestEdge != cNullItm)
          BMV_FAIL("PKT-empty", "pk=%lu count=0 smallestEdge=%ld", (unsigned long)pk, (long)P.smallestEdge);
      } else if (P.smallestEdge == cNullItm) {
        BMV_FAIL("PKT-empty", "pk=%lu count=%lu smallestEdge=null", (unsigned long)pk, (unsigned long)P.count);
      } else {
        Val pMin = 0; EdgeHandle pMinE = cNullItm; bool pHave = false;
        Size walked = validateHeapWalk(P.smallestEdge, b, true, P.rank, paperPure,
                                       seen, edgeOwner, pMin, pMinE, pHave, viol, &staleCnt);
        if (walked != P.count)
          BMV_FAIL("PKT-count", "pk=%lu count=%lu walked=%lu", (unsigned long)pk,
                   (unsigned long)P.count, (unsigned long)walked);
        if (pHave && edgeKey(P.smallestEdge) != edgeKey(pMinE))
          BMV_FAIL("PKT-min", "pk=%lu rootCost=%.17g walkMin=%.17g", (unsigned long)pk,
                   (double)mEdges[P.smallestEdge].cost, (double)pMin);
        if (pHave && (!have || edgeLess(pMinE, minE))) { minC = pMin; minE = pMinE; have = true; }
      }
      pk = P.nextInBlossom;
      if (++guard > mPackets.size() + 1) { BMV_FAIL("PKT-list-cycle", "b=%lu", (unsigned long)b); break; }
    }
    for (Size i = 0; i < touchedRanks.size(); ++i) rankSeen[touchedRanks[i]] = 0;
    // blossom minimum
    if (!have) {
      if (B.smallestOverall != cNullItm)
        BMV_FAIL("MIN-set", "b=%lu smallestOverall=%ld but no edges", (unsigned long)b, (long)B.smallestOverall);
    } else {
      if (B.smallestOverall == cNullItm)
        BMV_FAIL("MIN-set", "b=%lu smallestOverall=null but edges exist (min=%.17g)", (unsigned long)b, (double)minC);
      else {
        if (B.smallestOverall >= mEdges.size() || edgeOwner[B.smallestOverall] != b)
          BMV_FAIL("MIN-member", "b=%lu smallestOverall=%lu not in this blossom's containers",
                   (unsigned long)b, (unsigned long)B.smallestOverall);
        if (B.smallestCost != minC)
          BMV_FAIL("MIN-cost", "b=%lu cached=%.17g recomputed=%.17g", (unsigned long)b,
                   (double)B.smallestCost, (double)minC);
        else if (B.smallestOverall < mEdges.size() && mEdges[B.smallestOverall].cost != minC)
          BMV_FAIL("MIN-handle-cost", "b=%lu handleCost=%.17g min=%.17g", (unsigned long)b,
                   (double)mEdges[B.smallestOverall].cost, (double)minC);
#ifdef MBX_SHARED_D3_TIE_ORDER
        if (B.smallestOverall < mEdges.size() && edgeKey(B.smallestOverall) != edgeKey(minE))
          BMV_FAIL("MIN-rank", "b=%lu cachedRank=%lu recomputedRank=%lu", (unsigned long)b,
                   (unsigned long)mEdges[B.smallestOverall].tieRank,
                   (unsigned long)mEdges[minE].tieRank);
#endif
      }
    }
    // per-blossom H entry
    if (B.hHandle != cNullItm) ++hHandles;
    if (B.active && have) {
      if (B.hHandle == cNullItm)
        BMV_FAIL("H-missing", "b=%lu active with min=%.17g but no H entry", (unsigned long)b, (double)minC);
      else if (!mH.isActive(B.hHandle))
        BMV_FAIL("H-dead-handle", "b=%lu hHandle=%lu", (unsigned long)b, (unsigned long)B.hHandle);
      else if (mH.key(B.hHandle) != blossomKey(B))
        BMV_FAIL("H-key", "b=%lu key=%.17g smallestCost=%.17g", (unsigned long)b,
                 (double)keyCost(mH.key(B.hHandle)), (double)B.smallestCost);
      if (!haveGlobal || blossomKey(B) < globalMin) { globalMin = blossomKey(B); haveGlobal = true; }
    } else {
      if (B.hHandle != cNullItm)
        BMV_FAIL("H-spurious", "b=%lu %s has H entry", (unsigned long)b,
                 B.active ? "edgeless" : "inactive");
    }
  }
  // ---- coverage + free-list separation ----
  Size inUseEdges = 0;
  for (EdgeHandle e = 0; e < mEdges.size(); ++e) {
    if (mEdges[e].inUse) { ++inUseEdges;
      if (!seen[e])
        BMV_FAIL("I1i-unreachable-live", "e=%lu v=%lu w=%lu cost=%.17g", (unsigned long)e,
                 (unsigned long)mEdges[e].v, (unsigned long)mEdges[e].w, (double)mEdges[e].cost);
    }
  }
  if (inUseEdges != mAliveEdges)
    BMV_FAIL("I1i-alive-count", "mAliveEdges=%lu actual=%lu", (unsigned long)mAliveEdges, (unsigned long)inUseEdges);
  if (liveBlossoms != mAliveBlossoms)
    BMV_FAIL("I1i-alive-blossoms", "mAliveBlossoms=%lu actual=%lu", (unsigned long)mAliveBlossoms, (unsigned long)liveBlossoms);
  {
    std::vector<unsigned char> freeSeen(mEdges.size(), 0);
    for (Size i = 0; i < mFreeEdges.size(); ++i) {
      EdgeHandle e = mFreeEdges[i];
      if (e >= mEdges.size()) { BMV_FAIL("FREE-bounds", "edge freelist entry %lu", (unsigned long)e); continue; }
      if (mEdges[e].inUse) BMV_FAIL("FREE-live", "edge %lu on freelist but inUse", (unsigned long)e);
      if (freeSeen[e]) BMV_FAIL("FREE-dup", "edge %lu twice on freelist", (unsigned long)e);
      freeSeen[e] = 1;
    }
    if (mFreeEdges.size() + inUseEdges != mEdges.size())
      BMV_FAIL("FREE-cover", "edges: free=%lu inUse=%lu slots=%lu", (unsigned long)mFreeEdges.size(),
               (unsigned long)inUseEdges, (unsigned long)mEdges.size());
  }
  // ---- global H ----
  if (hHandles != mH.size())
    BMV_FAIL("H-size", "live handles=%lu mH.size=%lu", (unsigned long)hHandles, (unsigned long)mH.size());
  if (haveGlobal) {
    if (mH.empty())
      BMV_FAIL("H-empty", "global min should be %.17g but H empty", (double)keyCost(globalMin));
    else {
      if (mH.minKey() != globalMin)
        BMV_FAIL("H-globalmin", "H.minKey=%.17g recomputed=%.17g", (double)keyCost(mH.minKey()), (double)keyCost(globalMin));
      BlossomId mb = mH.minPayload();
      if (mb >= mBlossoms.size() || !mBlossoms[mb].inUse || !mBlossoms[mb].active
          || mBlossoms[mb].smallestOverall == cNullItm
          || blossomKey(mBlossoms[mb]) != mH.minKey())
        BMV_FAIL("H-minpayload", "minPayload=%ld inconsistent", (long)mb);
    }
  } else if (!mH.empty()) {
    BMV_FAIL("H-nonempty", "no active blossom has edges but H.size=%lu", (unsigned long)mH.size());
  }
  if (outStaleClassifications != NULL) *outStaleClassifications = staleCnt;
  return viol;
}

#undef BMV_FAIL

inline void BlossomMerging::collect_live_edges(std::vector<VertexId>& outV,
    std::vector<VertexId>& outW, std::vector<Val>& outC, std::vector<Size>& outPayload) const
{
  for (EdgeHandle e = 0; e < mEdges.size(); ++e) {
    if (!mEdges[e].inUse) continue;
    outV.push_back(mEdges[e].v);
    outW.push_back(mEdges[e].w);
    outC.push_back(mEdges[e].cost);
    outPayload.push_back(mEdges[e].payload);
  }
}
#endif // MBX_BM_VALIDATE

#ifdef MBX_BM_REPLACEMENT_VERIFY
// Local diagnostic only: expected O(old-ring + child leaves), with scratch
// proportional to old size, not mN. No global BM scan on each expansion.
inline void BlossomMerging::verify_frontier_partition(BlossomId Z,
    const std::vector<const std::vector<VertexId>*>& children) const
{
  if (!blossom_alive(Z) || mEdgeRemovalBatchActive || children.empty())
    throw std::runtime_error("BM replacement: invalid live-super/header/batch");
  const Blossom& B = mBlossoms[Z];
  if (B.size == 0 || B.size > mN || B.vtxListHead >= mN ||
      B.rootVertex >= mN || B.rank != logFloor(B.size))
    throw std::runtime_error("BM replacement: invalid old ring header");
  // Independently enumerate the ACTUAL old ring. Owner equality by itself
  // does not prove ring membership. No n-sized zero array per small expansion.
  std::unordered_set<VertexId> oldRing;
  oldRing.reserve(B.size);
  VertexId cur = B.vtxListHead, minV = cur;
  for (Size j = 0; j < B.size; ++j) {
    if (cur >= mN || mVtxBlossom[cur] != Z || !oldRing.insert(cur).second)
      throw std::runtime_error("BM replacement: old ring bounds/owner/duplicate");
    const VertexId next = mVtxNext[cur], prev = mVtxPrev[cur];
    if (next >= mN || prev >= mN || mVtxPrev[next] != cur || mVtxNext[prev] != cur)
      throw std::runtime_error("BM replacement: old ring broken links");
    if (mDepth[cur] < mDepth[minV]) minV = cur;
    cur = next;
  }
  if (cur != B.vtxListHead || oldRing.find(B.rootVertex) == oldRing.end())
    throw std::runtime_error("BM replacement: old ring size/root mismatch");
  // Consume the separately observed ring identity. erase()==0 catches both
  // foreign atoms and duplicates (including first-child numeric-ID ABA).
  Size total = 0;
  for (Size i = 0; i < children.size(); ++i) {
    if (children[i] == NULL || children[i]->empty())
      throw std::runtime_error("BM replacement: empty/null child");
    const std::vector<VertexId>& members = *children[i];
    for (Size j = 0; j < members.size(); ++j) {
      const VertexId v = members[j];
      if (v >= mN || mVtxBlossom[v] != Z || oldRing.erase(v) != 1)
        throw std::runtime_error("BM replacement: foreign/duplicate/bounds atom");
      ++total;
    }
  }
  if (total != B.size || !oldRing.empty())
    throw std::runtime_error("BM replacement: omitted old-ring atoms");
}

inline void BlossomMerging::verify_frontier_result(
    const std::vector<const std::vector<VertexId>*>& children,
    const std::vector<BlossomId>& ids) const
{
  if (children.size() != ids.size() || children.empty() || mEdgeRemovalBatchActive)
    throw std::runtime_error("BM replacement: incomplete result inventory/batch");
  std::unordered_set<BlossomId> idSeen;
  idSeen.reserve(ids.size());
  for (Size i = 0; i < children.size(); ++i) {
    if (!blossom_alive(ids[i]) || !idSeen.insert(ids[i]).second ||
        children[i] == NULL || children[i]->empty())
      throw std::runtime_error("BM replacement: invalid/repeated child identity");
    const Blossom& B = mBlossoms[ids[i]];
    const std::vector<VertexId>& members = *children[i];
    if (B.size != members.size() || B.vtxListHead != members[0] ||
        B.rank != logFloor(B.size))
      throw std::runtime_error("BM replacement: child size/head/rank mismatch");
    VertexId cur = B.vtxListHead, root = members[0];
    for (Size j = 0; j < members.size(); ++j) {
      if (cur >= mN || cur != members[j] || mVtxBlossom[cur] != ids[i])
        throw std::runtime_error("BM replacement: ordered child ring/owner mismatch");
      const VertexId prev = members[j == 0 ? members.size() - 1 : j - 1];
      const VertexId next = members[j + 1 == members.size() ? 0 : j + 1];
      if (mVtxPrev[cur] != prev || mVtxNext[cur] != next)
        throw std::runtime_error("BM replacement: child prev/next mismatch");
      if (mDepth[cur] < mDepth[root]) root = cur;
      cur = mVtxNext[cur];
    }
    if (cur != B.vtxListHead || B.rootVertex != root)
      throw std::runtime_error("BM replacement: ring end/depth-tie root mismatch");
  }
}
#endif // MBX_BM_REPLACEMENT_VERIFY

} // namespace Matchbox

#endif // BLOSSOM_MERGING_H
