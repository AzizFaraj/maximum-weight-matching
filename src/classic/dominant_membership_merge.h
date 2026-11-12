#ifndef MBX_DOMINANT_MEMBERSHIP_MERGE_20261026_H
#define MBX_DOMINANT_MEMBERSHIP_MERGE_20261026_H

#include <algorithm>
#include <cstddef>
#include <utility>
#include <vector>

// Sort a newly concatenated (vertex, immediate-child-index) map by reusing
// the largest child's already sorted vertex keys. This helper is independent
// of MatchingEngine and is compatible with C++98.
//
// Required engine invariant: childMap's first components are exactly the
// vertex keys in whole[blockStart, blockStart + blockCount). A child's old
// second components are deliberately ignored: they index a different cycle.
// The shape, sortedness, and uniqueness checks below do not prove set equality.
// The caller must not use a stale/incomplete child cache.
//
// The result is identical to std::sort(whole.begin(), whole.end()) under that
// invariant. Rejected preconditions leave whole unchanged and return false,
// so the caller must perform its ordinary sort. True means the sorted-child
// shortcut was used. No source child map
// is modified, and whole's size and capacity are unchanged.
template <class Size>
bool mbxDominantMembershipMerge(
  std::vector<std::pair<Size, Size> >& whole,
  const std::vector<std::pair<Size, Size> >& childMap,
  Size childIndex, std::size_t blockStart, std::size_t blockCount)
{
  typedef std::pair<Size, Size> Pair;
  const std::size_t n = whole.size();
  bool eligible = n > 1 && blockCount > 0 &&
                  &whole != &childMap &&
                  blockStart <= n && blockCount <= n - blockStart &&
                  blockCount >= n - blockCount &&
                  childMap.size() == blockCount;
  if (eligible) {
    for (std::size_t j = 0; j < blockCount; ++j) {
      if (whole[blockStart + j].second != childIndex ||
          (j > 0 && !(childMap[j - 1].first < childMap[j].first))) {
        eligible = false;
        break;
      }
    }
  }
  if (!eligible) {
    return false;
  }

  // Preserve all other child blocks before overwriting any part of whole.
  std::vector<Pair> remainder;
  remainder.reserve(n - blockCount);
  remainder.insert(remainder.end(), whole.begin(), whole.begin() + blockStart);
  remainder.insert(remainder.end(), whole.begin() + blockStart + blockCount,
                   whole.end());
  std::sort(remainder.begin(), remainder.end());

  for (std::size_t j = 0; j < blockCount; ++j)
    whole[j] = Pair(childMap[j].first, childIndex);
  std::copy(remainder.begin(), remainder.end(), whole.begin() + blockCount);
  std::inplace_merge(whole.begin(), whole.begin() + blockCount, whole.end());
  return true;
}

#endif
