#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

#include "Graph.h"
#include "MatchingAlgorithms.h"

using namespace Matchbox;

namespace {

unsigned int NextRandom(unsigned int* state) {
  *state = 1664525u * (*state) + 1013904223u;
  return *state;
}

Val BruteForce(const std::vector<std::vector<Val>>& weight,
               std::vector<bool>* used) {
  Size first = 0;
  while (first < used->size() && (*used)[first]) ++first;
  if (first == used->size()) return 0;

  (*used)[first] = true;
  Val best = BruteForce(weight, used);
  for (Size other = first + 1; other < used->size(); ++other) {
    if (!(*used)[other] && weight[first][other] > 0) {
      (*used)[other] = true;
      best = std::max(best, weight[first][other] + BruteForce(weight, used));
      (*used)[other] = false;
    }
  }
  (*used)[first] = false;
  return best;
}

Val BruteForce(const std::vector<std::vector<Val>>& weight) {
  std::vector<bool> used(weight.size(), false);
  return BruteForce(weight, &used);
}

void CheckRandomSparseGraphs() {
  for (Size n = 2; n <= 10; ++n) {
    for (unsigned int trial = 0; trial < 12; ++trial) {
      unsigned int state = 1009u + 97u * static_cast<unsigned int>(n) + trial;
      std::vector<Size> source;
      std::vector<Size> target;
      std::vector<Val> edgeWeight;
      std::vector<std::vector<Val>> weight(n, std::vector<Val>(n, 0));
      for (Size u = 0; u < n; ++u) {
        for (Size v = u + 1; v < n; ++v) {
          if ((NextRandom(&state) % 100u) < 55u) {
            const Val w = static_cast<Val>(1u + NextRandom(&state) % 50u);
            source.push_back(u);
            target.push_back(v);
            edgeWeight.push_back(w);
            weight[u][v] = weight[v][u] = w;
          }
        }
      }
      const Graph graph(n, source.size(), source, target, edgeWeight);
      const Val optimum = BruteForce(weight);
      assert(std::fabs(Edmonds(graph).weight - optimum) < 0.5);
      assert(std::fabs(Gabow(graph).weight - optimum) < 0.5);
    }
  }
}

void CheckAllThreeOnCompleteGraphs() {
  for (Size n = 2; n <= 10; n += 2) {
    for (unsigned int trial = 0; trial < 8; ++trial) {
      unsigned int state = 7919u + 193u * static_cast<unsigned int>(n) + trial;
      std::vector<Size> source;
      std::vector<Size> target;
      std::vector<Val> edgeWeight;
      std::vector<std::vector<Val>> weight(n, std::vector<Val>(n, 0));
      for (Size u = 0; u < n; ++u) {
        for (Size v = u + 1; v < n; ++v) {
          const Val w = static_cast<Val>(1u + NextRandom(&state) % 100u);
          source.push_back(u);
          target.push_back(v);
          edgeWeight.push_back(w);
          weight[u][v] = weight[v][u] = w;
        }
      }
      const Graph graph(n, source.size(), source, target, edgeWeight);
      const Val optimum = BruteForce(weight);
      for (const MatchingResult& result :
           {Edmonds(graph), Gabow(graph), Scaling(graph)}) {
        assert(result.status == SolveStatus::kSuccess);
        assert(result.cardinality == n / 2);
        assert(std::fabs(result.weight - optimum) < 0.5);
      }
    }
  }
}

}  // namespace

int main() {
  CheckRandomSparseGraphs();
  CheckAllThreeOnCompleteGraphs();
  std::cout << "randomized brute-force checks passed\n";
  return 0;
}
