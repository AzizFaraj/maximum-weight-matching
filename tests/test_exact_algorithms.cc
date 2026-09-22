#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

#include "Graph.h"
#include "MatchingAlgorithms.h"

using namespace Matchbox;

namespace {

Graph PerfectFixture() {
  const std::vector<Size> source{0, 2, 4, 0, 0, 1, 1, 2, 3};
  const std::vector<Size> target{1, 3, 5, 2, 4, 3, 5, 4, 5};
  const std::vector<Val> weight{10, 9, 8, 1, 1, 1, 1, 1, 1};
  return Graph(6, source.size(), source, target, weight);
}

void CheckResult(const Graph& graph, const MatchingResult& result) {
  assert(result.status == SolveStatus::kSuccess);
  assert(result.cardinality == 3);
  assert(std::fabs(result.weight - 27) < 0.5);
  for (Size v = 0; v < result.mate.size(); ++v) {
    const Size mate = result.mate[v];
    if (mate != cNullItm) {
      assert(mate < result.mate.size());
      assert(result.mate[mate] == v);
    }
  }
}

}  // namespace

int main() {
  const Graph graph = PerfectFixture();
  CheckResult(graph, Edmonds(graph));
  CheckResult(graph, Gabow(graph));
  CheckResult(graph, Scaling(graph));

  const std::vector<Size> source{0, 1, 0};
  const std::vector<Size> target{1, 2, 2};
  const std::vector<Val> weight{3, 2, 1};
  const Graph noPerfect(3, source.size(), source, target, weight);
  assert(Scaling(noPerfect).status == SolveStatus::kNoPerfectMatching);

  std::cout << "exact algorithm checks passed\n";
  return 0;
}
