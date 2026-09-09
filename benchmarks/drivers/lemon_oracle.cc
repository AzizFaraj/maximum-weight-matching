// Independent Suite B objective-eligibility oracle.
// Runs either LEMON unrestricted MWM or LEMON perfect MWPM on the same
// normalized stored-weight Matrix Market graph used by the measured solvers.
#include <algorithm>
#include <cerrno>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <string>
#include <vector>

#include <lemon/list_graph.h>
#include <lemon/matching.h>

typedef std::chrono::steady_clock OracleClock;

#ifndef HPC_SOURCE_COMMIT
#define HPC_SOURCE_COMMIT "unknown"
#endif
#ifndef HPC_SOURCE_STATE
#define HPC_SOURCE_STATE "unknown"
#endif
#ifndef HPC_BUILD_CONFIG
#define HPC_BUILD_CONFIG "unknown"
#endif

struct InputEdge {
  std::size_t u, v;
  long long weight;
};

static bool edgeLess(const InputEdge& a, const InputEdge& b) {
  return a.u < b.u || (a.u == b.u && a.v < b.v);
}

static void hashWord(std::uint64_t* hash, std::uint64_t word) {
  for (int i = 0; i < 8; ++i) {
    *hash ^= (word >> (8 * i)) & 0xffU;
    *hash *= 1099511628211ULL;
  }
}

static bool readGeneralMatrixMarket(const char* path, std::size_t* n,
                                    std::vector<InputEdge>* edges,
                                    std::string* error) {
  FILE* file = std::fopen(path, "r");
  if (!file) {
    *error = std::string("cannot open graph: ") + std::strerror(errno);
    return false;
  }
  char line[1 << 16];
  unsigned long long rows = 0, cols = 0, announced = 0;
  bool haveSize = false;
  while (std::fgets(line, sizeof(line), file)) {
    if (line[0] == '%') continue;
    if (std::sscanf(line, "%llu %llu %llu", &rows, &cols, &announced) == 3) {
      haveSize = true;
      break;
    }
  }
  if (!haveSize || rows != cols || rows > std::numeric_limits<std::size_t>::max()) {
    std::fclose(file);
    *error = "Suite B oracle requires valid square Matrix Market input";
    return false;
  }
  *n = (std::size_t)rows;
  edges->clear();
  if (announced <= std::numeric_limits<std::size_t>::max())
    edges->reserve((std::size_t)announced);

  unsigned long long a = 0, b = 0;
  long double value = 1.0L;
  while (std::fgets(line, sizeof(line), file)) {
    if (line[0] == '%') continue;
    value = 1.0L;
    const int fields = std::sscanf(line, "%llu %llu %Lf", &a, &b, &value);
    if (fields < 2) continue;
    if (a == 0 || b == 0 || a > rows || b > cols) {
      std::fclose(file); *error = "edge endpoint outside declared dimensions";
      return false;
    }
    --a; --b;
    if (a == b) continue;
    if (a > b) std::swap(a, b);
    if (!std::isfinite((double)value) || value < 0 ||
        value > (long double)std::numeric_limits<long long>::max()) {
      std::fclose(file); *error = "weight is negative, non-finite, or too large";
      return false;
    }
    const long long rounded = (long long)std::llround(value);
    if (std::fabs(value - (long double)rounded) > 1e-9L) {
      std::fclose(file); *error = "stored-weight mode requires integral weights";
      return false;
    }
    InputEdge edge = {(std::size_t)a, (std::size_t)b, rounded};
    edges->push_back(edge);
  }
  std::fclose(file);

  std::sort(edges->begin(), edges->end(), edgeLess);
  std::size_t out = 0;
  for (std::size_t i = 0; i < edges->size(); ++i) {
    if (out && (*edges)[out - 1].u == (*edges)[i].u &&
        (*edges)[out - 1].v == (*edges)[i].v) {
      if ((*edges)[i].weight > (*edges)[out - 1].weight)
        (*edges)[out - 1].weight = (*edges)[i].weight;
    } else {
      (*edges)[out++] = (*edges)[i];
    }
  }
  edges->resize(out);
  return true;
}

template <class Solver>
static bool summarizeMatching(const lemon::ListGraph& graph,
                              const Solver& solver, std::size_t n,
                              const lemon::ListGraph::EdgeMap<long long>& weights,
                              std::size_t* cardinality,
                              long long* exactWeight) {
  std::vector<unsigned char> used(n, 0);
  std::size_t card = 0;
  __int128 sum = 0;
  for (lemon::ListGraph::EdgeIt edge(graph); edge != lemon::INVALID; ++edge) {
    if (!solver.matching(edge)) continue;
    const int a = graph.id(graph.u(edge));
    const int b = graph.id(graph.v(edge));
    if (a < 0 || b < 0 || (std::size_t)a >= n || (std::size_t)b >= n ||
        used[(std::size_t)a] || used[(std::size_t)b]) return false;
    used[(std::size_t)a] = used[(std::size_t)b] = 1;
    ++card;
    sum += (__int128)weights[edge];
  }
  if (sum < std::numeric_limits<long long>::min() ||
      sum > std::numeric_limits<long long>::max()) return false;
  *cardinality = card;
  *exactWeight = (long long)sum;
  return true;
}

int main(int argc, char** argv) {
  if (argc != 3 || (std::strcmp(argv[1], "mwm") != 0 &&
                    std::strcmp(argv[1], "mwpm") != 0)) {
    std::fprintf(stderr, "usage: %s <mwm|mwpm> <general-graph.mtx>\n", argv[0]);
    return 2;
  }
  const bool perfectFormulation = std::strcmp(argv[1], "mwpm") == 0;
  std::size_t n = 0;
  std::vector<InputEdge> edges;
  std::string error;
  if (!readGeneralMatrixMarket(argv[2], &n, &edges, &error)) {
    std::fprintf(stderr, "INPUT_ERROR %s\n", error.c_str());
    return 2;
  }

  std::uint64_t digest = 1469598103934665603ULL;
  hashWord(&digest, n); hashWord(&digest, edges.size());
  lemon::ListGraph graph;
  std::vector<lemon::ListGraph::Node> nodes(n);
  for (std::size_t i = 0; i < n; ++i) nodes[i] = graph.addNode();
  lemon::ListGraph::EdgeMap<long long> weights(graph);
  for (std::size_t i = 0; i < edges.size(); ++i) {
    lemon::ListGraph::Edge edge = graph.addEdge(nodes[edges[i].u], nodes[edges[i].v]);
    weights[edge] = edges[i].weight;
    hashWord(&digest, edges[i].u); hashWord(&digest, edges[i].v);
    hashWord(&digest, (std::uint64_t)edges[i].weight);
  }

  bool feasible = true;
  bool valid = false;
  std::size_t card = 0;
  long long exactWeight = 0, reportedWeight = 0;
  const OracleClock::time_point begin = OracleClock::now();
  if (perfectFormulation) {
    lemon::MaxWeightedPerfectMatching<lemon::ListGraph,
      lemon::ListGraph::EdgeMap<long long> > solver(graph, weights);
    feasible = solver.run();
    const double seconds = std::chrono::duration<double>(OracleClock::now() - begin).count();
    if (feasible) {
      valid = summarizeMatching(graph, solver, n, weights, &card, &exactWeight);
      reportedWeight = solver.matchingWeight();
      valid = valid && card * 2 == n && exactWeight == reportedWeight;
    }
    std::printf("ORACLE schema=suite-b-lemon-v1 commit=%s source-state=%s build=%s "
                "formulation=mwpm graph=%s n=%zu m=%zu digest=%016llx\n",
                HPC_SOURCE_COMMIT, HPC_SOURCE_STATE, HPC_BUILD_CONFIG,
                argv[2], n, edges.size(), (unsigned long long)digest);
    std::printf("ORACLE_RESULT status=%s feasible=%d wall=%.6f card=%zu perfect=%d "
                "exact-weight=%lld reported-weight=%lld valid=%d\n",
                !feasible ? "infeasible" : (valid ? "success" : "invalid-result"),
                feasible ? 1 : 0, seconds, card, feasible && card * 2 == n ? 1 : 0,
                exactWeight, reportedWeight, valid ? 1 : 0);
  } else {
    lemon::MaxWeightedMatching<lemon::ListGraph,
      lemon::ListGraph::EdgeMap<long long> > solver(graph, weights);
    solver.run();
    const double seconds = std::chrono::duration<double>(OracleClock::now() - begin).count();
    valid = summarizeMatching(graph, solver, n, weights, &card, &exactWeight);
    reportedWeight = solver.matchingWeight();
    valid = valid && exactWeight == reportedWeight;
    std::printf("ORACLE schema=suite-b-lemon-v1 commit=%s source-state=%s build=%s "
                "formulation=mwm graph=%s n=%zu m=%zu digest=%016llx\n",
                HPC_SOURCE_COMMIT, HPC_SOURCE_STATE, HPC_BUILD_CONFIG,
                argv[2], n, edges.size(), (unsigned long long)digest);
    std::printf("ORACLE_RESULT status=%s feasible=1 wall=%.6f card=%zu perfect=%d "
                "exact-weight=%lld reported-weight=%lld valid=%d\n",
                valid ? "success" : "invalid-result", seconds, card,
                card * 2 == n ? 1 : 0, exactWeight, reportedWeight, valid ? 1 : 0);
  }
  std::fflush(stdout);
  return (!feasible || valid) ? 0 : 3;
}
