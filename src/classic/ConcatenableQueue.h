#ifndef CONCATENABLE_QUEUE_H
#define CONCATENABLE_QUEUE_H

// ============================================================
// ConcatenableQueue.h
//
// An augmented AVL tree supporting:
//   - insert(queue, priority, payload)        -> O(log n)
//   - find(elem)        -> queue identifier   -> O(log n)
//   - priority(elem)                          -> O(1)
//   - setPriority(elem, p)                    -> O(log n)
//   - minPriority(queue), minElem(queue)      -> O(1)
//   - merge(subQueues)  -> mergedQueue        -> O(k log n)
//   - split(mergedQueue) -> subQueues         -> O(k log n)
//
// Used by jorisvr-style fast weighted matching to maintain
// per-blossom queues of incident S-edges. When blossoms shrink,
// their queues are merged in O(log n); when expand, the merged
// queue splits back to the sub-blossoms in O(log n).
//
// Nodes are pool-allocated by Size index so handles are stable
// across vector growth. The sentinel "no handle / empty queue"
// is the global Matchbox::cNullItm.
//
// Pre-C++11 compatible. Single-threaded.
// ============================================================

#include <vector>
#include <cassert>
#include "Types.h"

namespace Matchbox {

template <bool Dual> struct ConcatenableQueueSecondary {
  void setSecondary(Val, Size) {}
  void setSecondaryPriority(Val) {}
  void setSecondaryLeaf(Size) {}
  Val secondaryPriority() const { return 0; }
  Size secondaryLeaf() const { return cNullItm; }
};

template <> struct ConcatenableQueueSecondary<true> {
  Val minPrio2;
  Size minLeaf2;
  void setSecondary(Val p, Size leaf) { minPrio2 = p; minLeaf2 = leaf; }
  void setSecondaryPriority(Val p) { minPrio2 = p; }
  void setSecondaryLeaf(Size leaf) { minLeaf2 = leaf; }
  Val secondaryPriority() const { return minPrio2; }
  Size secondaryLeaf() const { return minLeaf2; }
};

template <bool Dual> class ConcatenableQueueT {
public:
  typedef Size Handle;

  // Sentinel for "no queue / empty queue" is the global Matchbox::cNullItm.

  ConcatenableQueueT() : mFreeHead(cNullItm) {}

  explicit ConcatenableQueueT(Size initCapacity) : mFreeHead(cNullItm)
  {
    if (initCapacity > 0) mNodes.reserve(initCapacity);
  }

  // Reset everything. After this, all handles are invalid.
  void clear()
  {
    mNodes.clear();
    mFreeHead = cNullItm;
  }

  // ---- Node accessors (read-only on the leaf side) ----
  bool isLeaf(Handle h) const { return mNodes[h].leaf; }
  Val  priority(Handle h) const { assert(mNodes[h].leaf); return mNodes[h].minPrio; }
  Val  priority2(Handle h) const { assert(mNodes[h].leaf); return mNodes[h].secondaryPriority(); }
  Size payload(Handle h)  const { assert(mNodes[h].leaf); return mNodes[h].payload; }
  bool empty(Handle queueRoot) const { return queueRoot == cNullItm; }

  // ---- Min-tracking on a queue ----
  Val minPriority(Handle queueRoot) const
  {
    assert(queueRoot != cNullItm);
    return mNodes[queueRoot].minPrio;
  }
  Handle minElem(Handle queueRoot) const
  {
    assert(queueRoot != cNullItm);
    return mNodes[queueRoot].minLeaf;
  }
  Val minPriority2(Handle queueRoot) const
  {
    assert(queueRoot != cNullItm);
    return mNodes[queueRoot].secondaryPriority();
  }
  Handle minElem2(Handle queueRoot) const
  {
    assert(queueRoot != cNullItm);
    return mNodes[queueRoot].secondaryLeaf();
  }

  // ---- Queue root from a leaf ----
  // Walks parent pointers; O(log n).
  Handle findQueue(Handle leaf) const
  {
    Handle cur = leaf;
    while (mNodes[cur].parent != cNullItm)
      cur = mNodes[cur].parent;
    return cur;
  }

  // ---- Insert a leaf into the rightmost position of a queue ----
  // queueRoot is updated by reference. Returns the leaf handle.
  Handle insertRight(Handle& queueRoot, Val prio, Size pay)
  {
    Handle leaf = rAllocLeaf(prio, prio, pay);
    if (queueRoot == cNullItm) {
      queueRoot = leaf;
      return leaf;
    }
    queueRoot = rInsertAtRightmost(queueRoot, leaf);
    return leaf;
  }

  // Two independent minimum channels over the same concatenable topology.
  // Existing clients use insertRight(); SplitFindmin uses this form so its
  // INNER and Criteria 1--3 grow payloads can share all split/merge work.
  Handle insertRightDual(Handle& queueRoot, Val prio, Val prio2, Size pay)
  {
    Handle leaf = rAllocLeaf(prio, prio2, pay);
    if (queueRoot == cNullItm) {
      queueRoot = leaf;
      return leaf;
    }
    queueRoot = rInsertAtRightmost(queueRoot, leaf);
    return leaf;
  }

  // Symmetric: insert at leftmost.
  Handle insertLeft(Handle& queueRoot, Val prio, Size pay)
  {
    Handle leaf = rAllocLeaf(prio, prio, pay);
    if (queueRoot == cNullItm) {
      queueRoot = leaf;
      return leaf;
    }
    queueRoot = rInsertAtLeftmost(queueRoot, leaf);
    return leaf;
  }

  // ---- Update the priority of a leaf ----
  // Walks up to the root re-computing min augmentation. O(log n).
  Size setPriority(Handle leaf, Val newPrio)
  {
    assert(mNodes[leaf].leaf);
    mNodes[leaf].minPrio = newPrio;
    return rRefreshUp(mNodes[leaf].parent);
  }

  Size setPriority2(Handle leaf, Val newPrio)
  {
    assert(mNodes[leaf].leaf);
    mNodes[leaf].setSecondaryPriority(newPrio);
    return rRefreshUp(mNodes[leaf].parent);
  }

  // Monotone decrease-key fast paths. Refresh ancestors only until the
  // aggregate (minimum,witness) pair stops changing; no higher ancestor can
  // then change. Calling rRefreshNode preserves the queue's original
  // deterministic left-before-right tie witness exactly. setPriority and
  // setPriority2 remain the general paths for ownership resets that may
  // increase a key.
  Size decreasePriority(Handle leaf, Val newPrio, bool* changed = 0)
  {
    assert(mNodes[leaf].leaf);
    if (changed) *changed = false;
    if (newPrio >= mNodes[leaf].minPrio) return 0;
    if (changed) *changed = true;
    mNodes[leaf].minPrio = newPrio;
    Handle h = mNodes[leaf].parent;
    Size visits = 0;
    while (h != cNullItm) {
      ++visits;
      Val oldPrio = mNodes[h].minPrio;
      Handle oldLeaf = mNodes[h].minLeaf;
      rRefreshNode(h);
      if (mNodes[h].minPrio == oldPrio && mNodes[h].minLeaf == oldLeaf) break;
      h = mNodes[h].parent;
    }
    return visits;
  }

  Size decreasePriority2(Handle leaf, Val newPrio, bool* changed = 0)
  {
    assert(mNodes[leaf].leaf);
    if (changed) *changed = false;
    if (newPrio >= mNodes[leaf].secondaryPriority()) return 0;
    if (changed) *changed = true;
    mNodes[leaf].setSecondaryPriority(newPrio);
    Handle h = mNodes[leaf].parent;
    Size visits = 0;
    while (h != cNullItm) {
      ++visits;
      Val oldPrio = mNodes[h].secondaryPriority();
      Handle oldLeaf = mNodes[h].secondaryLeaf();
      rRefreshNode(h);
      if (mNodes[h].secondaryPriority() == oldPrio &&
          mNodes[h].secondaryLeaf() == oldLeaf) break;
      h = mNodes[h].parent;
    }
    return visits;
  }

  // Start a new epoch for the independent second channel in linear time.
  // Every live subtree has the same minimum afterwards, so its existing
  // minLeaf2 remains a valid witness and no per-leaf upward refresh is needed.
  Size resetPriorities2(Val newPrio)
  {
    Size visits = 0;
    for (Size i = 0; i < mNodes.size(); ++i)
      if (mNodes[i].inUse) {
        mNodes[i].setSecondaryPriority(newPrio);
        ++visits;
      }
    return visits;
  }

  // ---- Remove a leaf, returning the new root for that queue ----
  // O(log n). This is not strictly needed for the matching algorithm but
  // helpful for tests and for one-off deletions.
  Handle remove(Handle leaf)
  {
    assert(mNodes[leaf].leaf);
    Handle parent = mNodes[leaf].parent;
    if (parent == cNullItm) {
      // Single-leaf queue. Free the leaf, queue becomes empty.
      rFreeNode(leaf);
      return cNullItm;
    }
    Handle gp = mNodes[parent].parent;
    Handle sibling = (mNodes[parent].left == leaf)
                     ? mNodes[parent].right
                     : mNodes[parent].left;
    // Replace parent with sibling in the tree.
    mNodes[sibling].parent = gp;
    if (gp == cNullItm) {
      // Sibling becomes the new root.
      rFreeNode(leaf);
      rFreeNode(parent);
      return sibling;
    }
    if (mNodes[gp].left == parent) mNodes[gp].left = sibling;
    else mNodes[gp].right = sibling;
    rFreeNode(leaf);
    rFreeNode(parent);
    // Rebalance up from gp.
    return rRebalanceUp(gp);
  }

  // ---- Merge a list of sub-queues into one ----
  // subRoots may include cNullItm entries (empty queues), which are skipped.
  // Returns the merged queue root. After merge, subRoots are invalid as
  // independent queues; they are subtrees of the merged root.
  // O(k log n) where k = #non-empty subQueues, n = total elements.
  Handle merge(const std::vector<Handle>& subRoots)
  {
    Handle merged = cNullItm;
    for (Size i = 0; i < subRoots.size(); ++i) {
      if (subRoots[i] == cNullItm) continue;
      if (merged == cNullItm) merged = subRoots[i];
      else merged = rConcat(merged, subRoots[i]);
    }
    return merged;
  }

  // ---- Split a merged queue back into k sub-queues ----
  // splitLeaves[i] is a leaf handle that marks the LAST element of sub-queue i.
  // i.e., sub-queue i contains all elements up-to-and-including splitLeaves[i],
  // starting after splitLeaves[i-1]. The last sub-queue takes the remainder.
  //
  // Returns the resulting sub-queue roots into outRoots (resized to k+1).
  // After split, mergedRoot is invalid.
  //
  // O(k log n).
  void split(Handle mergedRoot,
             const std::vector<Handle>& splitLeaves,
             std::vector<Handle>& outRoots)
  {
    outRoots.clear();
    Handle cur = mergedRoot;
    for (Size i = 0; i < splitLeaves.size(); ++i) {
      Handle pivot = splitLeaves[i];
      // After this call, leftPart contains elements up to & incl. pivot,
      // and `cur` holds the rest.
      Handle leftPart;
      rSplitAfter(cur, pivot, leftPart, cur);
      outRoots.push_back(leftPart);
    }
    outRoots.push_back(cur);
  }

  // ---- Convenience: collect all leaves of a queue (for diagnostics/tests) ----
  void collectLeaves(Handle queueRoot, std::vector<Handle>& out) const
  {
    out.clear();
    if (queueRoot == cNullItm) return;
    rCollectLeaves(queueRoot, out);
  }

  // Number of elements (leaves) in a queue.
  Size size(Handle queueRoot) const
  {
    if (queueRoot == cNullItm) return 0;
    return mNodes[queueRoot].subtreeSize;
  }

  // Total nodes in pool (for diagnostics).
  Size poolSize() const { return mNodes.size(); }
  static Size nodeFootprint() { return sizeof(Node); }

  // Rightmost (last) leaf of a queue -- walk the right spine. O(log n).
  Handle rightmostLeaf(Handle queueRoot) const {
    Handle cur = queueRoot;
    if (cur == cNullItm) return cNullItm;
    while (!mNodes[cur].leaf) cur = mNodes[cur].right;
    return cur;
  }
/* CQ_BALANCE_DIAGNOSTIC_BEGIN */
#ifdef MBX_CQ_BALANCE_DIAGNOSTIC
  // Read-only, iterative diagnostic. No queue fields or operations are changed.
  struct DiagnosticShape {
    Size visited, leaves, actualHeight, storedRootHeight, avlBad;
    Size badHandles, badParents, badUse, badStructure, badSize, badHeight;
    Size badPrimary, badSecondary, cycles, duplicates;
    std::vector<Handle> leafHandles;
    std::vector<Size> payloads;
    DiagnosticShape() : visited(0), leaves(0), actualHeight(0), storedRootHeight(0),
      avlBad(0), badHandles(0), badParents(0), badUse(0), badStructure(0),
      badSize(0), badHeight(0), badPrimary(0), badSecondary(0), cycles(0), duplicates(0) {}
    Size errorsExceptBalance() const {
      return badHandles + badParents + badUse + badStructure + badSize + badHeight +
             badPrimary + badSecondary + cycles + duplicates;
    }
  };

  DiagnosticShape diagnosticInspect(Handle root) const {
    DiagnosticShape out;
    if (root == cNullItm) return out;
    if (root >= mNodes.size()) { ++out.badHandles; return out; }
    out.storedRootHeight = mNodes[root].height;
    const Size count = mNodes.size();
    std::vector<unsigned char> color(count, 0);
    std::vector<Size> sizes(count, 0), heights(count, 0);
    std::vector<Size> first(count, 0), last(count, 0), position(count, cNullItm);
    std::vector<Val> minimum(count, 0), minimum2(count, 0);
    struct Frame {
      Handle h, parent;
      bool exit;
      Frame(Handle a, Handle b, bool c) : h(a), parent(b), exit(c) {}
    };
    std::vector<Frame> pending;
    pending.push_back(Frame(root, cNullItm, false));
    while (!pending.empty()) {
      const Frame frame = pending.back();
      pending.pop_back();
      const Handle h = frame.h;
      if (h == cNullItm || h >= count) { ++out.badHandles; continue; }
      const Node& node = mNodes[h];
      if (!frame.exit) {
        if (node.parent != frame.parent) ++out.badParents;
        if (color[h] == 1) { ++out.cycles; continue; }
        if (color[h] == 2) { ++out.duplicates; continue; }
        color[h] = 1;
        ++out.visited;
        if (!node.inUse) ++out.badUse;
        if (node.leaf) {
          if (node.left != cNullItm || node.right != cNullItm) ++out.badStructure;
          sizes[h] = heights[h] = 1;
          minimum[h] = node.minPrio;
          minimum2[h] = node.secondaryPriority();
          first[h] = position[h] = out.leafHandles.size();
          out.leafHandles.push_back(h);
          out.payloads.push_back(node.payload);
          last[h] = out.leafHandles.size();
          ++out.leaves;
          if (node.subtreeSize != 1) ++out.badSize;
          if (node.height != 1) ++out.badHeight;
          if (node.minLeaf != h) ++out.badPrimary;
          if (Dual && node.secondaryLeaf() != h) ++out.badSecondary;
          color[h] = 2;
        } else {
          if (node.left == node.right) ++out.badStructure;
          pending.push_back(Frame(h, frame.parent, true));
          // LIFO right-then-left preserves the complete in-order leaf sequence.
          pending.push_back(Frame(node.right, h, false));
          pending.push_back(Frame(node.left, h, false));
        }
      } else {
        const Handle l = node.left, r = node.right;
        if (l == cNullItm || r == cNullItm || l >= count || r >= count ||
            color[l] != 2 || color[r] != 2) {
          ++out.badStructure;
          color[h] = 2;
          continue;
        }
        sizes[h] = sizes[l] + sizes[r];
        heights[h] = 1 + (heights[l] > heights[r] ? heights[l] : heights[r]);
        first[h] = first[l];
        last[h] = last[r];
        minimum[h] = minimum[l] <= minimum[r] ? minimum[l] : minimum[r];
        minimum2[h] = minimum2[l] <= minimum2[r] ? minimum2[l] : minimum2[r];
        if ((heights[l] > heights[r] && heights[l] - heights[r] > 1) ||
            (heights[r] > heights[l] && heights[r] - heights[l] > 1)) ++out.avlBad;
        if (node.subtreeSize != sizes[h]) ++out.badSize;
        if (node.height != heights[h]) ++out.badHeight;
        const Handle witness = node.minLeaf;
        if (node.minPrio != minimum[h] || witness == cNullItm || witness >= count ||
            !mNodes[witness].leaf || !mNodes[witness].inUse ||
            position[witness] == cNullItm || position[witness] < first[h] ||
            position[witness] >= last[h] || mNodes[witness].minPrio != minimum[h])
          ++out.badPrimary;
        if (Dual) {
          const Handle witness2 = node.secondaryLeaf();
          // resetPriorities2 may retain any valid tied witness, not necessarily
          // the leftmost one. Validate value and subtree membership instead.
          if (node.secondaryPriority() != minimum2[h] || witness2 == cNullItm || witness2 >= count ||
              !mNodes[witness2].leaf || !mNodes[witness2].inUse ||
              position[witness2] == cNullItm || position[witness2] < first[h] ||
              position[witness2] >= last[h] || mNodes[witness2].secondaryPriority() != minimum2[h])
            ++out.badSecondary;
        }
        color[h] = 2;
      }
    }
    out.actualHeight = heights[root];
    return out;
  }
#endif
/* CQ_BALANCE_DIAGNOSTIC_END */

private:
  struct Node : public ConcatenableQueueSecondary<Dual> {
    // Common fields:
    Handle parent;
    Val    minPrio;       // for leaves: own priority. for internal: min over subtree
    Handle minLeaf;       // leaf achieving minPrio in this subtree
    Size   subtreeSize;   // # leaves in subtree
    Size   height;        // AVL height (1 for leaf, 1 + max(child heights) for internal)
    bool   leaf;
    bool   inUse;

    // Leaf-only:
    Size payload;

    // Internal-only:
    Handle left, right;
  };

  std::vector<Node> mNodes;
  Handle mFreeHead;       // head of free-list, linked via Node::parent

  // ---- Pool management ----
  Handle rAllocNode()
  {
    if (mFreeHead != cNullItm) {
      Handle h = mFreeHead;
      mFreeHead = mNodes[h].parent;
      mNodes[h].inUse = true;
      return h;
    }
    Node n;
    n.parent = cNullItm; n.left = cNullItm; n.right = cNullItm;
    n.minPrio = 0; n.minLeaf = cNullItm;
    n.setSecondary(0, cNullItm);
    n.subtreeSize = 0; n.height = 0;
    n.leaf = false; n.inUse = true;
    n.payload = cNullItm;
    mNodes.push_back(n);
    return mNodes.size() - 1;
  }

  Handle rAllocLeaf(Val prio, Val prio2, Size pay)
  {
    Handle h = rAllocNode();
    Node& n = mNodes[h];
    n.parent = cNullItm;
    n.left = cNullItm; n.right = cNullItm;
    n.minPrio = prio; n.minLeaf = h;
    n.setSecondary(prio2, h);
    n.subtreeSize = 1; n.height = 1;
    n.leaf = true;
    n.payload = pay;
    return h;
  }

  Handle rAllocInternal(Handle l, Handle r)
  {
    Handle h = rAllocNode();
    Node& n = mNodes[h];
    n.parent = cNullItm;
    n.left = l; n.right = r;
    n.leaf = false;
    mNodes[l].parent = h;
    mNodes[r].parent = h;
    rRefreshNode(h);
    return h;
  }

  void rFreeNode(Handle h)
  {
    mNodes[h].inUse = false;
    mNodes[h].parent = mFreeHead;
    mFreeHead = h;
  }

  // ---- Augmentation / balancing helpers ----
  Size rHeight(Handle h) const { return (h == cNullItm) ? 0 : mNodes[h].height; }

  void rRefreshNode(Handle h)
  {
    Node& n = mNodes[h];
    if (n.leaf) {
      n.subtreeSize = 1;
      n.height = 1;
      n.minLeaf = h;
      n.setSecondaryLeaf(h);
      // minPrio already set by caller
      return;
    }
    const Node& L = mNodes[n.left];
    const Node& R = mNodes[n.right];
    n.subtreeSize = L.subtreeSize + R.subtreeSize;
    Size hl = L.height, hr = R.height;
    n.height = 1 + (hl > hr ? hl : hr);
    if (L.minPrio <= R.minPrio) {
      n.minPrio = L.minPrio;
      n.minLeaf = L.minLeaf;
    } else {
      n.minPrio = R.minPrio;
      n.minLeaf = R.minLeaf;
    }
    if (L.secondaryPriority() <= R.secondaryPriority()) {
      n.setSecondary(L.secondaryPriority(), L.secondaryLeaf());
    } else {
      n.setSecondary(R.secondaryPriority(), R.secondaryLeaf());
    }
  }

  Size rRefreshUp(Handle h)
  {
    Size visits = 0;
    while (h != cNullItm) {
      ++visits;
      rRefreshNode(h);
      h = mNodes[h].parent;
    }
    return visits;
  }

  // Standard AVL rotations. Each rotation returns the new root of the
  // rotated subtree (for the caller to attach to its parent).
  Handle rRotateLeft(Handle x)
  {
    Handle y = mNodes[x].right;
    Handle yL = mNodes[y].left;
    mNodes[y].left = x;
    mNodes[x].right = yL;
    if (yL != cNullItm) mNodes[yL].parent = x;
    mNodes[y].parent = mNodes[x].parent;
    mNodes[x].parent = y;
    rRefreshNode(x);
    rRefreshNode(y);
    return y;
  }

  Handle rRotateRight(Handle x)
  {
    Handle y = mNodes[x].left;
    Handle yR = mNodes[y].right;
    mNodes[y].right = x;
    mNodes[x].left = yR;
    if (yR != cNullItm) mNodes[yR].parent = x;
    mNodes[y].parent = mNodes[x].parent;
    mNodes[x].parent = y;
    rRefreshNode(x);
    rRefreshNode(y);
    return y;
  }

  // Rebalance node h. Returns new root of subtree at h.
  Handle rRebalance(Handle h)
  {
    if (mNodes[h].leaf) return h;
    Size hl = rHeight(mNodes[h].left);
    Size hr = rHeight(mNodes[h].right);
    if (hl > hr + 1) {
      Handle c = mNodes[h].left;
      Size cl = rHeight(mNodes[c].left);
      Size cr = rHeight(mNodes[c].right);
      if (cr > cl) {
        // left-right: rotate left on c first.
        Handle newLeft = rRotateLeft(c);
        mNodes[h].left = newLeft;
        mNodes[newLeft].parent = h;
        rRefreshNode(h);
      }
      return rRotateRight(h);
    }
    if (hr > hl + 1) {
      Handle c = mNodes[h].right;
      Size cl = rHeight(mNodes[c].left);
      Size cr = rHeight(mNodes[c].right);
      if (cl > cr) {
        Handle newRight = rRotateRight(c);
        mNodes[h].right = newRight;
        mNodes[newRight].parent = h;
        rRefreshNode(h);
      }
      return rRotateLeft(h);
    }
    rRefreshNode(h);
    return h;
  }

  // Walk up from h, rebalancing each ancestor. Returns the (possibly new) root.
  Handle rRebalanceUp(Handle h)
  {
    Handle cur = h;
    Handle root = cNullItm;
    while (cur != cNullItm) {
      Handle parent = mNodes[cur].parent;
      Handle newSub = rRebalance(cur);
      if (parent == cNullItm) {
        root = newSub;
        mNodes[root].parent = cNullItm;
        break;
      }
      // Re-parent newSub if rotation changed the subtree root.
      if (mNodes[parent].left == cur) mNodes[parent].left = newSub;
      else mNodes[parent].right = newSub;
      mNodes[newSub].parent = parent;
      cur = parent;
    }
    return root;
  }

  // ---- Insert at extremal position (helpers used by insertLeft/Right) ----
  Handle rInsertAtRightmost(Handle root, Handle leaf)
  {
    // Walk to rightmost path, attach.
    Handle cur = root;
    while (!mNodes[cur].leaf) cur = mNodes[cur].right;
    return rInsertBeside(cur, leaf, /*onRight=*/true);
  }

  Handle rInsertAtLeftmost(Handle root, Handle leaf)
  {
    Handle cur = root;
    while (!mNodes[cur].leaf) cur = mNodes[cur].left;
    return rInsertBeside(cur, leaf, /*onRight=*/false);
  }

  // Insert leaf as new sibling of an existing leaf, on the specified side.
  Handle rInsertBeside(Handle existing, Handle leaf, bool onRight)
  {
    Handle parent = mNodes[existing].parent;
    Handle newInternal;
    if (onRight) newInternal = rAllocInternal(existing, leaf);
    else         newInternal = rAllocInternal(leaf, existing);
    if (parent == cNullItm) {
      mNodes[newInternal].parent = cNullItm;
      return newInternal;
    }
    if (mNodes[parent].left == existing) mNodes[parent].left = newInternal;
    else mNodes[parent].right = newInternal;
    mNodes[newInternal].parent = parent;
    return rRebalanceUp(parent);
  }

  // ---- Concatenate two trees: A's leaves come before B's leaves ----
  // Both A and B are non-NULL. Returns new root.
  // O(|height(A) - height(B)|) plus rebalancing.
  Handle rConcat(Handle A, Handle B)
  {
    Size hA = mNodes[A].height;
    Size hB = mNodes[B].height;
    if (hA <= hB + 1 && hB <= hA + 1) {
      // Heights are close enough — just create a new internal parent.
      Handle root = rAllocInternal(A, B);
      mNodes[root].parent = cNullItm;
      return root;
    }
    if (hA > hB) {
      // A is taller. Walk down A's right spine until height matches hB+1, splice B there.
      Handle cur = A;
      while (!mNodes[cur].leaf && rHeight(cur) > hB + 1)
        cur = mNodes[cur].right;
      // cur is internal node whose right child has height ≤ hB+1, OR cur is a leaf (height 1)
      Handle parentOfCur = mNodes[cur].parent;
      Handle newInternal = rAllocInternal(cur, B);
      if (parentOfCur == cNullItm) {
        // cur was the root.
        mNodes[newInternal].parent = cNullItm;
        return newInternal;
      }
      // cur was a right-child of parentOfCur (since we walked right).
      // Replace cur with newInternal.
      mNodes[parentOfCur].right = newInternal;
      mNodes[newInternal].parent = parentOfCur;
      return rRebalanceUp(parentOfCur);
    }
    // hB > hA: walk down B's left spine.
    Handle cur = B;
    while (!mNodes[cur].leaf && rHeight(cur) > hA + 1)
      cur = mNodes[cur].left;
    Handle parentOfCur = mNodes[cur].parent;
    Handle newInternal = rAllocInternal(A, cur);
    if (parentOfCur == cNullItm) {
      mNodes[newInternal].parent = cNullItm;
      return newInternal;
    }
    mNodes[parentOfCur].left = newInternal;
    mNodes[newInternal].parent = parentOfCur;
    return rRebalanceUp(parentOfCur);
  }

  // ---- Split tree at a leaf ----
  // After call: leftOut contains all leaves up-to-and-including pivot (in-order);
  //             rightOut contains the rest.
  // Splits in O(log n).
  void rSplitAfter(Handle root, Handle pivot, Handle& leftOut, Handle& rightOut)
  {
    // Walk from pivot up to root, accumulating left and right pieces.
    // Standard AVL split via path tracing.
    Handle leftAcc = cNullItm;
    Handle rightAcc = cNullItm;

    // Start at pivot leaf — leftAcc contains pivot.
    leftAcc = pivot;
    // Detach pivot from its parent context as we walk up.
    Handle prev = pivot;
    Handle cur = mNodes[pivot].parent;
    while (cur != cNullItm) {
      Handle p = mNodes[cur].parent;
      if (mNodes[cur].left == prev) {
        // We came up from the left. The right subtree of cur joins rightAcc.
        Handle rightSubtree = mNodes[cur].right;
        // Detach
        mNodes[rightSubtree].parent = cNullItm;
        if (rightAcc == cNullItm) rightAcc = rightSubtree;
        else rightAcc = rConcatTopLevel(rightAcc, rightSubtree);
      } else {
        // We came up from the right. The left subtree of cur joins leftAcc on the left.
        Handle leftSubtree = mNodes[cur].left;
        mNodes[leftSubtree].parent = cNullItm;
        if (leftAcc == cNullItm) leftAcc = leftSubtree;
        else leftAcc = rConcatTopLevel(leftSubtree, leftAcc);
      }
      // The internal node `cur` is no longer needed; free it.
      // Caller should not have any handle to it beyond what we've already processed.
      rFreeNode(cur);
      prev = cur;
      cur = p;
    }
    // Detach leftAcc, rightAcc from anything.
    if (leftAcc != cNullItm) mNodes[leftAcc].parent = cNullItm;
    if (rightAcc != cNullItm) mNodes[rightAcc].parent = cNullItm;
    leftOut = leftAcc;
    rightOut = rightAcc;
    (void)root; // root is now decomposed; we drove the walk from pivot.
  }

  // Wrapper around rConcat that handles cNullItm inputs and ensures result.parent=NULL.
  Handle rConcatTopLevel(Handle A, Handle B)
  {
    if (A == cNullItm) return B;
    if (B == cNullItm) return A;
    Handle r = rConcat(A, B);
    if (r != cNullItm) mNodes[r].parent = cNullItm;
    return r;
  }

  // ---- Diagnostics ----
  void rCollectLeaves(Handle h, std::vector<Handle>& out) const
  {
    if (mNodes[h].leaf) { out.push_back(h); return; }
    rCollectLeaves(mNodes[h].left, out);
    rCollectLeaves(mNodes[h].right, out);
  }
};

typedef ConcatenableQueueT<false> ConcatenableQueue;
typedef ConcatenableQueueT<true> DualConcatenableQueue;

} // namespace Matchbox

#endif // CONCATENABLE_QUEUE_H
