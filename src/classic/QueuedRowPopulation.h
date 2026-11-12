#ifndef MATCHBOX_EXACT_QUEUED_ROW_POPULATION_H
#define MATCHBOX_EXACT_QUEUED_ROW_POPULATION_H
#include <vector>
#include <stdexcept>

namespace Matchbox {
// Exact representation invariant: population[u] equals the number of nonzero
// flags in atomic-vertex row u. The caller initializes those flags to zero.
template<class Index> class QueuedRowPopulation {
 public:
  void resetZeroRows(Index rows) { population_.assign(rows, 0); }
  Index population(Index row) const { return population_[row]; }
  template<class Rows> void set(Rows& flags, Index row, Index col, char value) {
    const bool before = flags[row][col] != 0;
    const bool after = value != 0;
    if (before != after) {
      if (after) ++population_[row];
      else {
        if (population_[row] == 0) throw std::logic_error("queued-row-population-underflow");
        --population_[row];
      }
    }
    flags[row][col] = value;
#ifdef MBX_GABOW_EMPTY_ROW_VERIFY
    verify(flags, row);
#endif
  }
  template<class Rows> bool clearRow(Rows& flags, Index row) {
#ifdef MBX_GABOW_EMPTY_ROW_VERIFY
    verify(flags, row);
#endif
    if (population_[row] == 0) return false;
    for (Index col = 0; col < flags[row].size(); ++col) flags[row][col] = 0;
    population_[row] = 0;
#ifdef MBX_GABOW_EMPTY_ROW_VERIFY
    verify(flags, row);
#endif
    return true;
  }
  template<class Rows> void verify(const Rows& flags, Index row) const {
    Index actual = 0;
    for (Index col = 0; col < flags[row].size(); ++col) actual += flags[row][col] != 0;
    if (actual != population_[row]) throw std::logic_error("queued-row-population-mismatch");
  }
 private:
  std::vector<Index> population_;
};
}
#endif
