#include <chrono>
#include <iomanip>
#include <iostream>
#include <string>

#include "Graph.h"
#include "MatchingAlgorithms.h"

using Matchbox::Graph;
using Matchbox::MatchingResult;
using Matchbox::Size;
using Matchbox::SolveStatus;

namespace {

void Usage(const char* program) {
  std::cerr << "usage: " << program
            << " <edmonds|gabow|scaling> <graph.mtx>\n";
}

}  // namespace

int main(int argc, char** argv) {
  if (argc != 3) {
    Usage(argv[0]);
    return 2;
  }

  const std::string algorithm = argv[1];
  Graph graph;
  if (graph.Read(argv[2]) != Matchbox::eErrNone) {
    std::cerr << "input-error: could not read " << argv[2] << '\n';
    return 2;
  }

  const auto start = std::chrono::steady_clock::now();
  MatchingResult result;
  if (algorithm == "edmonds") {
    result = Matchbox::Edmonds(graph);
  } else if (algorithm == "gabow") {
    result = Matchbox::Gabow(graph);
  } else if (algorithm == "scaling") {
    result = Matchbox::Scaling(graph);
  } else {
    Usage(argv[0]);
    return 2;
  }
  const double seconds = std::chrono::duration<double>(
      std::chrono::steady_clock::now() - start).count();

  Size vertices = 0;
  graph.GetNumVtxs(&vertices);
  std::cout << std::setprecision(17)
            << "algorithm=" << algorithm
            << " status=" << Matchbox::SolveStatusName(result.status)
            << " wall=" << seconds
            << " matching_size=" << result.cardinality
            << " matching_weight=" << result.weight
            << " perfect=" << (2 * result.cardinality == vertices ? 1 : 0)
            << '\n';

  return result.status == SolveStatus::kSuccess ? 0 : 3;
}
