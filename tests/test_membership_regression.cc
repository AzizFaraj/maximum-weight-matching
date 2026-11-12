#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

#include "Graph.h"
#include "MatchingAlgorithms.h"

using namespace Matchbox;

namespace {

void Check(const Graph& graph, const MatchingResult& result) {
  Size n = 0;
  graph.GetNumVtxs(&n);
  assert(result.status == SolveStatus::kSuccess);
  assert(result.mate.size() == n);
  assert(result.cardinality == 332);
  assert(std::fabs(result.weight - 273865) < 0.5);
  Size pairs = 0;
  for (Size u = 0; u < n; ++u) {
    const Size v = result.mate[u];
    if (v == cNullItm) continue;
    assert(v < n && v != u && result.mate[v] == u);
    if (u < v) ++pairs;
  }
  assert(pairs == result.cardinality);

  std::vector<Size> source, target;
  std::vector<Val> weights;
  graph.GetEdgs(&source, &target, &weights);
  Val total = 0;
  Size selected = 0;
  for (Size e = 0; e < source.size(); ++e) {
    if (result.mate[source[e]] == target[e]) {
      total += weights[e];
      ++selected;
    }
  }
  assert(selected == result.cardinality);
  assert(std::fabs(total - result.weight) < 0.5);
}

}  // namespace

int main(int argc, char** argv) {
  assert(argc == 2);
  Graph graph;
  assert(graph.Read(argv[1]) == eErrNone);
  const MatchingResult edmonds = Edmonds(graph);
  const MatchingResult gabow = Gabow(graph);
  Check(graph, edmonds);
  Check(graph, gabow);
  assert(edmonds.mate == gabow.mate);
  std::cout << "676-vertex membership regression passed; identical classic mates\n";
  return 0;
}
