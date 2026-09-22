#ifndef FIBONACCI_HEAP_H
#define FIBONACCI_HEAP_H

#include <vector>
#include <cassert>
#include <cmath>
#include "Types.h"

namespace Matchbox {

// ============================================================
// FibonacciHeap: A Fibonacci heap with O(1) amortized insert
// and decrease-key, O(log n) amortized delete-min.
//
// Uses index-based node handles into a pre-allocated pool.
// Handle = Size index; cNullItm = null sentinel.
//
// Required by Gabow's 2018 weighted matching algorithm for
// the global blossom-merging heap, per-blossom grow/expand
// heaps, and the main D1/D2/D4 event heaps.
//
// Pre-C++98 compatible. No C++11 features used.
// ============================================================
template<typename Key, typename Payload>
class FibonacciHeap {
public:
  typedef Size Handle;

  struct Node {
    Key key;
    Payload payload;
    Size parent;   // cNullItm if root
    Size child;    // one child (cNullItm if leaf)
    Size left;     // circular doubly-linked sibling
    Size right;    // circular doubly-linked sibling
    Size degree;   // number of children
    bool mark;     // lost a child since becoming non-root
    bool inUse;    // false if on free list
  };

  // Construct with initial pool capacity (grows if needed).
  explicit FibonacciHeap(Size initCapacity = 0)
    : mMin(cNullItm), mSize(0), mFreeHead(cNullItm)
  {
    if (initCapacity > 0) {
      mNodes.reserve(initCapacity);
    }
  }

  // Reset heap to empty state, recycling all nodes.
  void clear()
  {
    for (Size i = 0; i < mNodes.size(); ++i) {
      mNodes[i].inUse = false;
    }
    mMin = cNullItm;
    mSize = 0;
    // Rebuild free list.
    mFreeHead = cNullItm;
    for (Size i = mNodes.size(); i > 0; --i) {
      mNodes[i - 1].right = mFreeHead;
      mFreeHead = i - 1;
    }
  }

  bool empty() const { return mSize == 0; }
  Size size() const { return mSize; }

  // Access the minimum without removal.
  // Precondition: !empty().
  Key minKey() const
  {
    assert(mMin != cNullItm);
    return mNodes[mMin].key;
  }

  Payload minPayload() const
  {
    assert(mMin != cNullItm);
    return mNodes[mMin].payload;
  }

  Handle minHandle() const
  {
    assert(mMin != cNullItm);
    return mMin;
  }

  // Read key/payload of any live node.
  Key key(Handle h) const { return mNodes[h].key; }
  Payload payload(Handle h) const { return mNodes[h].payload; }

  // Insert a new element. Returns a handle for future decrease-key.
  // O(1) amortized.
  Handle insert(Key k, Payload p)
  {
    Size idx = rAllocNode();
    Node& nd = mNodes[idx];
    nd.key = k;
    nd.payload = p;
    nd.parent = cNullItm;
    nd.child = cNullItm;
    nd.left = idx;
    nd.right = idx;
    nd.degree = 0;
    nd.mark = false;
    nd.inUse = true;

    if (mMin == cNullItm) {
      mMin = idx;
    } else {
      rInsertIntoList(mMin, idx);
      if (nd.key < mNodes[mMin].key) {
        mMin = idx;
      }
    }
    ++mSize;
    return idx;
  }

  // Remove and return the minimum element.
  // O(log n) amortized.
  // Precondition: !empty().
  void deleteMin(Key* outKey, Payload* outPayload)
  {
    assert(mMin != cNullItm);
    Size z = mMin;
    const Node& zn = mNodes[z];
    if (outKey != NULL) *outKey = zn.key;
    if (outPayload != NULL) *outPayload = zn.payload;

    // Add all children of z to root list.
    if (zn.child != cNullItm) {
      // Walk child ring, set parent = cNullItm.
      Size c = zn.child;
      do {
        mNodes[c].parent = cNullItm;
        c = mNodes[c].right;
      } while (c != zn.child);
      // Splice child ring into root list.
      rSpliceLists(mMin, zn.child);
    }

    // Remove z from root list.
    if (z == mNodes[z].right) {
      // z was the only root (and has no children after splice --
      // but children were spliced in, so check again).
      // Actually: if z.right == z after removing children, heap is now empty.
      // But we spliced children first, so z.right may have changed.
      // Need to re-read since splice modifies links.
      if (mNodes[z].right == z) {
        mMin = cNullItm;
      } else {
        mMin = mNodes[z].right;
        rRemoveFromList(z);
        rConsolidate();
      }
    } else {
      mMin = mNodes[z].right;
      rRemoveFromList(z);
      rConsolidate();
    }

    rFreeNode(z);
    --mSize;
  }

  // Decrease the key of a node.
  // O(1) amortized.
  // Precondition: newKey <= current key.
  void decreaseKey(Handle h, Key newKey)
  {
    assert(mNodes[h].inUse);
    Node& x = mNodes[h];
    assert(!(x.key < newKey)); // newKey must be <= current
    x.key = newKey;

    Size p = x.parent;
    if (p != cNullItm && x.key < mNodes[p].key) {
      rCut(h, p);
      rCascadingCut(p);
    }
    if (x.key < mNodes[mMin].key) {
      mMin = h;
    }
  }

  // Remove an arbitrary node.
  // O(log n) amortized.
  void remove(Handle h)
  {
    // Cut to root with artificially minimal key, then delete-min.
    // We implement this without needing a "negative infinity" key:
    // just cut to root and manually set as min, then deleteMin.
    Node& x = mNodes[h];
    Size p = x.parent;
    if (p != cNullItm) {
      rCut(h, p);
      rCascadingCut(p);
    }
    mMin = h;
    Key dummyK;
    Payload dummyP;
    deleteMin(&dummyK, &dummyP);
  }

  // Merge another heap into this one. The other heap becomes empty.
  // O(1).
  void merge(FibonacciHeap& other)
  {
    if (other.mMin == cNullItm) return;

    // We need to absorb other's node pool into ours.
    // Remap indices: other's node i becomes i + offset in our pool.
    Size offset = mNodes.size();
    Size otherPoolSize = other.mNodes.size();
    mNodes.resize(mNodes.size() + otherPoolSize);
    for (Size i = 0; i < otherPoolSize; ++i) {
      mNodes[offset + i] = other.mNodes[i];
      Node& nd = mNodes[offset + i];
      if (nd.inUse) {
        if (nd.parent != cNullItm) nd.parent += offset;
        if (nd.child != cNullItm) nd.child += offset;
        nd.left += offset;
        nd.right += offset;
      }
    }

    Size otherMin = other.mMin + offset;

    // Remap other's free list too.
    if (other.mFreeHead != cNullItm) {
      Size f = other.mFreeHead;
      while (f != cNullItm) {
        Size next = other.mNodes[f].right;
        Size remapped = f + offset;
        mNodes[remapped].right = (next != cNullItm) ? next + offset : mFreeHead;
        f = next;
      }
      mFreeHead = other.mFreeHead + offset;
    }

    if (mMin == cNullItm) {
      mMin = otherMin;
    } else {
      rSpliceLists(mMin, otherMin);
      if (mNodes[otherMin].key < mNodes[mMin].key) {
        mMin = otherMin;
      }
    }
    mSize += other.mSize;

    other.mMin = cNullItm;
    other.mSize = 0;
    other.mNodes.clear();
    other.mFreeHead = cNullItm;
  }

  // Check if a handle refers to a live node.
  bool isActive(Handle h) const
  {
    return h < mNodes.size() && mNodes[h].inUse;
  }

  // Update the payload of a node (key unchanged).
  void setPayload(Handle h, Payload p)
  {
    mNodes[h].payload = p;
  }

private:
  std::vector<Node> mNodes;
  Size mMin;       // index of min root (cNullItm if empty)
  Size mSize;      // number of active elements
  Size mFreeHead;  // head of free list (linked via right pointers)

  // Allocate a node from the pool.
  Size rAllocNode()
  {
    if (mFreeHead != cNullItm) {
      Size idx = mFreeHead;
      mFreeHead = mNodes[idx].right;
      return idx;
    }
    Size idx = mNodes.size();
    mNodes.resize(mNodes.size() + 1);
    return idx;
  }

  // Return a node to the free list.
  void rFreeNode(Size idx)
  {
    mNodes[idx].inUse = false;
    mNodes[idx].right = mFreeHead;
    mFreeHead = idx;
  }

  // Insert node b into the circular list containing node a.
  // b must be a singleton (b.left == b.right == b).
  void rInsertIntoList(Size a, Size b)
  {
    Size aRight = mNodes[a].right;
    mNodes[a].right = b;
    mNodes[b].left = a;
    mNodes[b].right = aRight;
    mNodes[aRight].left = b;
  }

  // Splice two circular lists together.
  // a and b are members of different circular lists.
  void rSpliceLists(Size a, Size b)
  {
    if (a == cNullItm || b == cNullItm) return;
    // a -- aRight  becomes  a -- b -- ... -- bLeft -- aRight
    // b -- bRight           (and bRight ... remains connected through b)
    Size aRight = mNodes[a].right;
    Size bLeft = mNodes[b].left;
    mNodes[a].right = b;
    mNodes[b].left = a;
    mNodes[bLeft].right = aRight;
    mNodes[aRight].left = bLeft;
  }

  // Remove node x from its circular list.
  // Precondition: x is not the only element in its list,
  // or caller handles the empty case.
  void rRemoveFromList(Size x)
  {
    Size xl = mNodes[x].left;
    Size xr = mNodes[x].right;
    mNodes[xl].right = xr;
    mNodes[xr].left = xl;
    mNodes[x].left = x;
    mNodes[x].right = x;
  }

  // Make y a child of x (both roots; y.key >= x.key).
  void rLink(Size y, Size x)
  {
    // Remove y from root list.
    rRemoveFromList(y);
    // Make y a child of x.
    mNodes[y].parent = x;
    mNodes[y].mark = false;
    if (mNodes[x].child == cNullItm) {
      mNodes[x].child = y;
      mNodes[y].left = y;
      mNodes[y].right = y;
    } else {
      rInsertIntoList(mNodes[x].child, y);
    }
    ++(mNodes[x].degree);
  }

  // Consolidate: combine roots of the same degree.
  void rConsolidate()
  {
    if (mMin == cNullItm) return;

    // Max degree is O(log n). For n elements, max degree <= 1.4405 * log2(n) + 1.
    // Use 64 as safe upper bound (handles up to 2^44 elements).
    static const Size MAX_DEGREE = 64;
    Size degTable[MAX_DEGREE];
    for (Size i = 0; i < MAX_DEGREE; ++i) {
      degTable[i] = cNullItm;
    }

    // Collect all current roots into a flat list first,
    // because consolidation modifies the root list.
    std::vector<Size> roots;
    {
      Size r = mMin;
      do {
        roots.push_back(r);
        r = mNodes[r].right;
      } while (r != mMin);
    }

    for (Size ri = 0; ri < roots.size(); ++ri) {
      Size x = roots[ri];
      Size d = mNodes[x].degree;
      while (d < MAX_DEGREE && degTable[d] != cNullItm) {
        Size y = degTable[d];
        if (mNodes[y].key < mNodes[x].key) {
          // swap x and y so x has smaller key
          Size tmp = x; x = y; y = tmp;
        }
        rLink(y, x);
        degTable[d] = cNullItm;
        ++d;
      }
      assert(d < MAX_DEGREE);
      degTable[d] = x;
    }

    // Rebuild root list and find new min.
    mMin = cNullItm;
    for (Size i = 0; i < MAX_DEGREE; ++i) {
      if (degTable[i] != cNullItm) {
        Size nd = degTable[i];
        mNodes[nd].left = nd;
        mNodes[nd].right = nd;
        mNodes[nd].parent = cNullItm;
        if (mMin == cNullItm) {
          mMin = nd;
        } else {
          rInsertIntoList(mMin, nd);
          if (mNodes[nd].key < mNodes[mMin].key) {
            mMin = nd;
          }
        }
      }
    }
  }

  // Cut x from its parent p, adding x to the root list.
  void rCut(Size x, Size p)
  {
    // Remove x from p's child list.
    if (mNodes[x].right == x) {
      // x is p's only child.
      mNodes[p].child = cNullItm;
    } else {
      if (mNodes[p].child == x) {
        mNodes[p].child = mNodes[x].right;
      }
      rRemoveFromList(x);
    }
    --(mNodes[p].degree);

    // Add x to root list.
    mNodes[x].parent = cNullItm;
    mNodes[x].mark = false;
    mNodes[x].left = x;
    mNodes[x].right = x;
    rInsertIntoList(mMin, x);
  }

  // Cascading cut from node y upward.
  void rCascadingCut(Size y)
  {
    Size p = mNodes[y].parent;
    if (p != cNullItm) {
      if (!mNodes[y].mark) {
        mNodes[y].mark = true;
      } else {
        rCut(y, p);
        rCascadingCut(p);
      }
    }
  }
};

} // namespace Matchbox

#endif // FIBONACCI_HEAP_H
