#ifndef DYNAMIC_NCA_H
#define DYNAMIC_NCA_H

#include <vector>
#include <cassert>
#include "Types.h"

namespace Matchbox {

// ============================================================
// DynamicNCA: NCA on an incrementally-growing forest.
//
// Implementation: binary-lifting / jump-pointer tables, maintained
// incrementally in addLeaf. Each addLeaf is O(log depth) time AND
// space (it materializes the full jump table for the new leaf, one
// vector per node -> ~24 B vector header + 8 B/level, cache-hostile
// vector-of-vectors layout); each nca query is O(log depth). There
// is NO global rebuild step, so interleaved addLeaf/nca sequences
// (matchbox's access pattern) stay cheap in practice.
//
// COMPLEXITY HONESTY (Phase 4A): this does NOT meet the bound Gabow's
// search analysis assumes. The paper charges total O(n) for all
// add_leaf calls and O(m) for all NCA queries per search (a
// specialized incremental-tree NCA, e.g. Gabow's own static-tree /
// union-find style structures). Binary lifting gives O(n log n) +
// O(m log n) per search — a log-factor deviation on BOTH terms.
// Whether that log factor matters in practice is an empirical
// question (see the MBX_GABOW_SEARCH_AUDIT counters below and
// docs/Gabow_Code_Structure_Audit.md); do not cite this structure as
// satisfying the paper's stated NCA bound.
//
// (A Bender-Farach-Colton sparse-table RMQ gives O(1) queries but
// needs an O(n log n) rebuild after any addLeaf. Because matchbox
// interleaves addLeaf and nca within a search, that rebuild fired on
// nearly every query and dominated runtime on dense graphs. The
// jump-pointer structure is the better of the two tried options for
// interleaved use; it is still not the paper's structure.)
//
// Pre-C++98 compatible.
// ============================================================

#ifdef MBX_GABOW_SEARCH_AUDIT
// Phase-4A audit counters (compile-gated; zero cost in normal builds).
inline unsigned long long& dncaNcaCalls()     { static unsigned long long c=0; return c; }
inline unsigned long long& dncaNcaSteps()     { static unsigned long long c=0; return c; } // lifting iterations
inline unsigned long long& dncaAddLeafCalls() { static unsigned long long c=0; return c; }
inline unsigned long long& dncaAddLeafSteps() { static unsigned long long c=0; return c; } // jump-table entries built
inline unsigned long long& dncaMaxLift()      { static unsigned long long c=0; return c; } // max lifting steps in one query
#endif
class DynamicNCA {
public:
  explicit DynamicNCA(Size maxNodes = 0) : mMaxDepth(0)
  {
    if (maxNodes > 0) {
      mParent.resize(maxNodes, cNullItm);
      mDepth.resize(maxNodes, 0);
      mAnc.resize(maxNodes);
    }
  }

  void reset(Size maxNodes)
  {
    mParent.assign(maxNodes, cNullItm);
    mDepth.assign(maxNodes, 0);
    mAnc.assign(maxNodes, std::vector<Size>());
    mMaxDepth = 0;
  }

  void clearNode(Size x)
  {
    if (x >= mParent.size()) return;
    mParent[x] = cNullItm;
    mDepth[x] = 0;
    mAnc[x].clear();
  }

  void addLeaf(Size x, Size par)
  {
    if (x >= mParent.size()) {
      Size newSize = x + 1;
      mParent.resize(newSize, cNullItm);
      mDepth.resize(newSize, 0);
      mAnc.resize(newSize);
    }
    mParent[x] = par;
#ifdef MBX_GABOW_SEARCH_AUDIT
    ++dncaAddLeafCalls();
#endif
    if (par == cNullItm) {
      mDepth[x] = 0;
      mAnc[x].clear();
    } else {
      mDepth[x] = mDepth[par] + 1;
      if (mDepth[x] > mMaxDepth) mMaxDepth = mDepth[x];
      std::vector<Size>& tbl = mAnc[x];
      tbl.clear();
      tbl.push_back(par);
      Size cur = par;
      while (true) {
        const std::vector<Size>& curTbl = mAnc[cur];
        Size lvl = tbl.size() - 1;
        if (lvl >= curTbl.size()) break;
        Size jump = curTbl[lvl];
        if (jump == cNullItm) break;
        tbl.push_back(jump);
        cur = jump;
      }
#ifdef MBX_GABOW_SEARCH_AUDIT
      dncaAddLeafSteps() += tbl.size();
#endif
    }
  }

  // Nearest common ancestor of u and v. O(log depth).
  Size nca(Size u, Size v)
  {
    if (u >= mParent.size() || v >= mParent.size()) return cNullItm;
#ifdef MBX_GABOW_SEARCH_AUDIT
    ++dncaNcaCalls();
    unsigned long long lifts = 0;
#endif
    Size a = u, b = v;
    if (mDepth[a] < mDepth[b]) { Size t = a; a = b; b = t; }
    // Lift a up to b's depth.
    Size diff = mDepth[a] - mDepth[b];
    for (Size k = 0; diff > 0; ++k, diff >>= 1) {
      if (diff & 1) {
        if (k >= mAnc[a].size()) return cNullItm;
        a = mAnc[a][k];
#ifdef MBX_GABOW_SEARCH_AUDIT
        ++lifts;
#endif
        if (a == cNullItm) return cNullItm;
      }
    }
    if (a == b) {
#ifdef MBX_GABOW_SEARCH_AUDIT
      dncaNcaSteps() += lifts;
      if (lifts > dncaMaxLift()) dncaMaxLift() = lifts;
#endif
      return a;
    }
    // Lift both until parents converge.
    Size maxK = mAnc[a].size();
    if (mAnc[b].size() < maxK) maxK = mAnc[b].size();
    while (maxK > 0) {
      --maxK;
      if (maxK < mAnc[a].size() && maxK < mAnc[b].size() &&
          mAnc[a][maxK] != mAnc[b][maxK]) {
        a = mAnc[a][maxK];
        b = mAnc[b][maxK];
#ifdef MBX_GABOW_SEARCH_AUDIT
        lifts += 2;
#endif
      }
    }
#ifdef MBX_GABOW_SEARCH_AUDIT
    dncaNcaSteps() += lifts;
    if (lifts > dncaMaxLift()) dncaMaxLift() = lifts;
#endif
    return mParent[a];
  }

  Size depth(Size x) const { return mDepth[x]; }
  Size parent(Size x) const { return mParent[x]; }
  Size maxDepth() const { return mMaxDepth; }  // DIAG: deepest node ever added

  bool exists(Size x) const
  {
    return x < mParent.size() &&
           (mParent[x] != cNullItm || mDepth[x] == 0);
  }

private:
  std::vector<Size> mParent;
  std::vector<Size> mDepth;
  std::vector<std::vector<Size> > mAnc;  // mAnc[v][k] = 2^k-th ancestor of v
  Size mMaxDepth;  // DIAG: max depth ever assigned (tracks T growth w/o reset)
};

} // namespace Matchbox

#endif // DYNAMIC_NCA_H
