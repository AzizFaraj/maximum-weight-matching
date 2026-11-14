#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

#include "Graph.h"
#include "MatchingAlgorithms.h"

using namespace Matchbox;

namespace {

void Validate(const Graph& graph, const MatchingResult& result, Size n,
              Size expectedCardinality, Val expectedWeight) {
  assert(result.status == SolveStatus::kSuccess);
  assert(result.mate.size() == n);
  assert(result.cardinality == expectedCardinality);
  assert(std::fabs(result.weight - expectedWeight) < 0.5);
  Size pairs = 0;
  Val weight = 0;
  for (Size u = 0; u < n; ++u) {
    Size degree = 0;
    std::vector<Size> neighbors;
    std::vector<Val> weights;
    assert(graph.GetAdj(u, &degree, &neighbors, &weights) == eErrNone);
    const Size v = result.mate[u];
    if (degree == 0) assert(v == cNullItm);
    if (v == cNullItm) continue;
    assert(v < n && v != u && result.mate[v] == u);
    bool found = false;
    for (Size i = 0; i < degree; ++i) {
      if (neighbors[i] == v) {
        found = true;
        if (u < v) weight += weights[i];
        break;
      }
    }
    assert(found);
    if (u < v) ++pairs;
  }
  assert(pairs == expectedCardinality);
  assert(std::fabs(weight - expectedWeight) < 0.5);
}

void Check(Size n, const std::vector<Size>& source, const std::vector<Size>& target,
           const std::vector<Val>& weight, Size cardinality, Val optimum) {
  const Graph graph(n, source.size(), source, target, weight);
  Err error = eErrNone;
  graph.GetInitErr(&error);
  assert(error == eErrNone);
  const MatchingResult edmonds = Edmonds(graph);
  const MatchingResult gabow = Gabow(graph);
  Validate(graph, edmonds, n, cardinality, optimum);
  Validate(graph, gabow, n, cardinality, optimum);
  assert(edmonds.mate == gabow.mate);
}

}  // namespace

int main() {
  // Original CI failure: the last CSR row starts at the backing vector's end.
  Check(3, {0}, {1}, {35}, 1, 35);
  Check(5, {0}, {1}, {35}, 1, 35);  // Several trailing empty rows.
  Check(3, {1}, {2}, {35}, 1, 35);  // Leading empty row.
  Check(3, {0}, {2}, {35}, 1, 35);  // Middle empty row.
  Check(6, {0, 2}, {1, 3}, {35, 20}, 2, 55);  // Disconnected + trailing.
  Check(8, {1, 4}, {2, 5}, {35, 20}, 2, 55);  // Leading/middle/trailing.
  for (Size n : {Size(0), Size(1), Size(4), Size(12)}) {
    Check(n, {}, {}, {}, 0, 0);  // Completely empty CSR backing storage.
  }
  std::cout << "empty-row regression checks passed (both classical backends)\n";
  return 0;
}
