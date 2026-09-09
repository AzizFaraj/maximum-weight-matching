// One-process, one-engine benchmark for the paper's Suite B campaign.
// Edmonds, Gabow, and Scaling consume the same normalized stored-weight graph.
// No weight offset, padding, vertex deletion, or two-copy reduction is applied.
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

#include "Types.h"
#include "Graph.h"
#include "MatchingEngine.h"
#include "ScalingMatching.h"

using namespace Matchbox;
typedef std::chrono::steady_clock BenchClock;

#ifndef HPC_SOURCE_COMMIT
#define HPC_SOURCE_COMMIT "unknown"
#endif
#ifndef HPC_BUILD_CONFIG
#define HPC_BUILD_CONFIG "unknown"
#endif
#ifndef HPC_SOURCE_STATE
#define HPC_SOURCE_STATE "unknown"
#endif

struct InputEdge {
  Size u, v;
  long long fileWeight;
};

static bool edgeLess(const InputEdge& a, const InputEdge& b) {
  return a.u < b.u || (a.u == b.u && a.v < b.v);
}

static std::uint64_t splitmix64(std::uint64_t& state) {
  std::uint64_t z = (state += 0x9e3779b97f4a7c15ULL);
  z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
  z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
  return z ^ (z >> 31);
}

static void hashWord(std::uint64_t* hash, std::uint64_t word) {
  for (int i = 0; i < 8; ++i) {
    *hash ^= (word >> (8 * i)) & 0xffU;
    *hash *= 1099511628211ULL;
  }
}

static bool readMatrixMarket(const char* path, const std::string& graphKind,
                             Size* n,
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
  if (!haveSize || rows > std::numeric_limits<Size>::max() ||
      cols > std::numeric_limits<Size>::max()) {
    std::fclose(file); *error = "invalid MatrixMarket dimensions"; return false;
  }
  if (graphKind == "general") {
    if (rows != cols) {
      std::fclose(file);
      *error = "general graph MatrixMarket input must be square";
      return false;
    }
    *n = (Size)rows;
  } else {
    if (rows > (unsigned long long)std::numeric_limits<Size>::max() - cols) {
      std::fclose(file);
      *error = "bipartite row+column vertex count overflows Size";
      return false;
    }
    *n = (Size)(rows + cols);
  }
  edges->clear();
  if (announced <= std::numeric_limits<size_t>::max())
    edges->reserve((size_t)announced);
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
    if (graphKind == "bipartite") {
      b += rows;
    } else if (a == b) {
      continue;
    }
    if (a > b) std::swap(a, b);
    if (!std::isfinite((double)value) || value < 0 ||
        value > (long double)std::numeric_limits<long long>::max()) {
      std::fclose(file); *error = "file weight is negative, non-finite, or too large";
      return false;
    }
    const long long rounded = (long long)std::llround(value);
    if (std::fabs(value - (long double)rounded) > 1e-9L) {
      std::fclose(file); *error = "file-weight mode requires integral weights";
      return false;
    }
    InputEdge e; e.u = (Size)a; e.v = (Size)b; e.fileWeight = rounded;
    edges->push_back(e);
  }
  std::fclose(file);
  std::sort(edges->begin(), edges->end(), edgeLess);
  size_t out = 0;
  for (size_t i = 0; i < edges->size(); ++i) {
    if (out && (*edges)[out - 1].u == (*edges)[i].u &&
        (*edges)[out - 1].v == (*edges)[i].v) {
      if ((*edges)[i].fileWeight > (*edges)[out - 1].fileWeight)
        (*edges)[out - 1].fileWeight = (*edges)[i].fileWeight;
    } else {
      (*edges)[out++] = (*edges)[i];
    }
  }
  edges->resize(out);
  return true;
}

static void usage(const char* argv0) {
  std::fprintf(stderr,
    "usage: %s <edmonds|gabow|scaling> <graph.mtx> <general|bipartite> "
    "<file|uniform> <W> <seed> <require-perfect:0|1> <trial>\n", argv0);
}

int main(int argc, char** argv) {
  if (argc != 9) { usage(argv[0]); return 2; }
  const std::string engine = argv[1];
  const std::string graphPath = argv[2];
  const std::string graphKind = argv[3];
  const std::string weightMode = argv[4];
  const long long W = std::strtoll(argv[5], 0, 10);
  const std::uint64_t seed = std::strtoull(argv[6], 0, 10);
  const bool requirePerfect = std::atoi(argv[7]) != 0;
  const unsigned long long trial = std::strtoull(argv[8], 0, 10);
  if ((engine != "edmonds" && engine != "gabow" && engine != "scaling") ||
      (graphKind != "general" && graphKind != "bipartite") ||
      (weightMode != "file" && weightMode != "uniform") ||
      (weightMode == "uniform" && W <= 0) || trial == 0) {
    usage(argv[0]); return 2;
  }
  if (engine == "scaling" && graphKind != "general") {
    std::fprintf(stderr, "INPUT_ERROR Scaling accepts Suite B general-graph inputs only\n");
    return 2;
  }

  Size n = 0; std::vector<InputEdge> edges; std::string error;
  if (!readMatrixMarket(graphPath.c_str(), graphKind, &n, &edges, &error)) {
    std::fprintf(stderr, "INPUT_ERROR %s\n", error.c_str()); return 2;
  }
  std::vector<Size> source(edges.size()), target(edges.size());
  std::vector<long long> exactWeights(edges.size());
  std::vector<Val> weights(edges.size());
  std::uint64_t randomState = seed ^ 0xd1b54a32d192ed03ULL;
  std::uint64_t digest = 1469598103934665603ULL;
  hashWord(&digest, n); hashWord(&digest, edges.size());
  for (size_t i = 0; i < edges.size(); ++i) {
    source[i] = edges[i].u; target[i] = edges[i].v;
    exactWeights[i] = weightMode == "file" ? edges[i].fileWeight :
      1 + (long long)(splitmix64(randomState) % (std::uint64_t)W);
    weights[i] = (Val)exactWeights[i];
    hashWord(&digest, source[i]); hashWord(&digest, target[i]);
    hashWord(&digest, (std::uint64_t)exactWeights[i]);
  }

  Graph graph(n, (Size)edges.size(), source, target, weights);
  const char* requestedAlgorithm = engine == "edmonds" ? "edmonds" :
    (engine == "gabow" ? "gabow" : "hybrid-figure6");
  const char* requestedBackend = engine == "edmonds" ? "gmg1986" :
    (engine == "gabow" ? "gabow2018" : "integrated-event-materialized-shell");
#ifdef SCALING_S5_FULL_KEY_DRIVER_DISABLE
  const char* requestedFullKeyDriver = engine == "scaling" ? "disabled" : "n/a";
#else
  const char* requestedFullKeyDriver = engine == "scaling" ? "enabled" : "n/a";
#endif
  const char* requestedConfiguration = engine == "scaling" ? "practical" : "production";
  // Emit and flush immutable input/build identity before entering the solver.
  // A timeout, SIGKILL, or crash must still leave enough evidence to identify
  // the exact graph and implementation that was attempted.
  std::printf("RUN_START schema=suite-b-hpc-v2 commit=%s source-state=%s build=%s "
              "engine=%s algorithm=%s backend=%s full-key-driver=%s configuration=%s "
              "graph=%s graph-kind=%s n=%lu m=%lu digest=%016llx "
              "weight-mode=%s W=%lld seed=%llu trial=%llu require-perfect=%d\n",
              HPC_SOURCE_COMMIT, HPC_SOURCE_STATE, HPC_BUILD_CONFIG,
              engine.c_str(), requestedAlgorithm, requestedBackend,
              requestedFullKeyDriver, requestedConfiguration,
              graphPath.c_str(), graphKind.c_str(),
              (unsigned long)n, (unsigned long)edges.size(),
              (unsigned long long)digest, weightMode.c_str(), W,
              (unsigned long long)seed, trial, requirePerfect ? 1 : 0);
  std::fflush(stdout);

  MatchingEngine matching;
  std::vector<Size> mate; Size engineCard = 0; Val engineWeight = 0;
  bool solverOk = true;
  const char* solverStatus = "success";
  const char* scalingAlgorithm = requestedAlgorithm;
  const char* scalingBackend = requestedBackend;
  const char* scalingFullKeyDriver = requestedFullKeyDriver;
  const char* scalingLeafRepresentation = "n/a";
  Size scalingSqrtN = 0, scalingTau = 0;
  int scalingFvrStage1Criterion = 0, scalingFvrStage2Criterion = 0;
  int scalingDismantleStage2Criterion = 0, scalingFinalizationCriterion = 0;
  bool scalingIncrementalShrink = false, scalingIncrementalDissolve = false;
  bool scalingPersistentSearch = false, scalingEventSearch = false;
  bool scalingLazyDualClock = false, scalingEventBatchSameTime = false;
  bool scalingEventCompactInitialization = false, scalingEventBoundaryReuse = false;
  bool scalingCriterion2FreeFrontier = false;
  bool scalingSplitFindminBoundary = false, scalingSplitFindminScanOracle = false;
  bool scalingSplitFindminGrow = false, scalingSplitFindminGrowScanOracle = false;
  bool scalingSharedSplitFindminTopology = false, scalingSplitFindminListFrontier = false;
  bool scalingDetailedScheduleTiming = false, scalingDynamicShell = false;
  bool scalingSensitivityParameters = false;
  const BenchClock::time_point begin = BenchClock::now();
  if (engine == "edmonds")
    matching.ComputeMaxEdgWghtMatchingEdmondsFastV2(
      graph, &mate, &engineCard, &engineWeight);
  else if (engine == "gabow")
    matching.ComputeMaxEdgWghtMatchingGabow(
      graph, &mate, &engineCard, &engineWeight);
  else {
    ScalingMatching scaling;
    solverOk = scaling.computeHybrid(graph, &mate, &engineCard, &engineWeight);
    solverStatus = ScalingMatching::statusName(scaling.lastStatus());
    const ScalingMatching::RunConfiguration cfg = scaling.runConfiguration();
    scalingAlgorithm = ScalingMatching::algorithmName(cfg.algorithm);
    scalingBackend = cfg.backend;
    scalingLeafRepresentation = ScalingMatching::leafRepresentationName();
#ifdef SCALING_S5_FULL_KEY_DRIVER_DISABLE
    scalingFullKeyDriver = "disabled";
#else
    scalingFullKeyDriver = "enabled";
#endif
    scalingSqrtN = cfg.sqrtN;
    scalingTau = cfg.tau;
    scalingFvrStage1Criterion = cfg.fvrStage1Criterion;
    scalingFvrStage2Criterion = cfg.fvrStage2Criterion;
    scalingDismantleStage2Criterion = cfg.dismantleStage2Criterion;
    scalingFinalizationCriterion = cfg.finalizationCriterion;
    scalingIncrementalShrink = cfg.incrementalShrink;
    scalingIncrementalDissolve = cfg.incrementalDissolve;
    scalingPersistentSearch = cfg.persistentSearch;
    scalingEventSearch = cfg.eventSearch;
    scalingLazyDualClock = cfg.lazyDualClock;
    scalingEventBatchSameTime = cfg.eventBatchSameTime;
    scalingEventCompactInitialization = cfg.eventCompactInitialization;
    scalingEventBoundaryReuse = cfg.eventBoundaryReuse;
    scalingCriterion2FreeFrontier = cfg.eventCriterion2FreeFrontier;
    scalingSplitFindminBoundary = cfg.eventSplitFindminBoundary;
    scalingSplitFindminScanOracle = cfg.eventSplitFindminScanOracle;
    scalingSplitFindminGrow = cfg.eventSplitFindminGrow;
    scalingSplitFindminGrowScanOracle = cfg.eventSplitFindminGrowScanOracle;
    scalingSharedSplitFindminTopology = cfg.eventSharedSplitFindminTopology;
    scalingSplitFindminListFrontier = cfg.eventSplitFindminListFrontier;
    scalingDetailedScheduleTiming = cfg.detailedScheduleTiming;
    scalingDynamicShell = cfg.dynamicShell;
    scalingSensitivityParameters = cfg.sensitivityParameters;
  }
  const double seconds = std::chrono::duration<double>(
    BenchClock::now() - begin).count();

  Size mutual = 0, asymmetric = 0;
  for (Size v = 0; v < n; ++v) {
    const Size m = v < mate.size() ? mate[v] : cNullItm;
    if (m >= n) continue;
    if (m < mate.size() && mate[m] == v) {
      if (v < m) ++mutual;
    } else {
      ++asymmetric;
    }
  }
  long long exactWeight = 0;
  Size matchedInputEdges = 0;
  for (size_t i = 0; i < edges.size(); ++i)
    if (source[i] < mate.size() && target[i] < mate.size() &&
        mate[source[i]] == target[i] && mate[target[i]] == source[i]) {
      exactWeight += exactWeights[i];
      ++matchedInputEdges;
    }
  const long double engineRounded = std::llround((long double)engineWeight);
  const bool integralEngineWeight =
    std::fabs((long double)engineWeight - engineRounded) <= 1e-7L;
  const bool weightAgrees = integralEngineWeight &&
    engineRounded == (long double)exactWeight;
  const bool valid = solverOk && asymmetric == 0 && mutual == engineCard &&
    matchedInputEdges == mutual && weightAgrees;
  const bool perfect = valid && mutual * 2 == n;
  const bool accepted = valid && (!requirePerfect || perfect);

  std::printf("RUN schema=suite-b-hpc-v2 commit=%s source-state=%s build=%s engine=%s graph=%s "
              "graph-kind=%s n=%lu m=%lu digest=%016llx weight-mode=%s "
              "W=%lld seed=%llu trial=%llu require-perfect=%d "
              "algorithm=%s backend=%s full-key-driver=%s configuration=%s "
              "sqrtN=%lu tau=%lu leaf-representation=%s "
              "fvr-stage1-criterion=%d fvr-stage2-criterion=%d "
              "dismantle-stage2-criterion=%d finalization-criterion=%d "
              "incremental-shrink=%d incremental-dissolve=%d "
              "persistent-search=%d event-search=%d lazy-dual-clock=%d "
              "event-batch-same-time=%d event-compact-initialization=%d "
              "event-boundary-reuse=%d criterion2-free-frontier=%d "
              "split-findmin-boundary=%d split-findmin-scan-oracle=%d "
              "split-findmin-grow=%d split-findmin-grow-scan-oracle=%d "
              "shared-split-findmin-topology=%d split-findmin-list-frontier=%d "
              "detailed-schedule-timing=%d dynamic-shell=%d "
              "sensitivity-parameters=%d\n",
              HPC_SOURCE_COMMIT, HPC_SOURCE_STATE, HPC_BUILD_CONFIG,
              engine.c_str(), graphPath.c_str(),
              graphKind.c_str(),
              (unsigned long)n, (unsigned long)edges.size(),
              (unsigned long long)digest, weightMode.c_str(), W,
              (unsigned long long)seed, trial, requirePerfect ? 1 : 0,
              scalingAlgorithm, scalingBackend, scalingFullKeyDriver,
              requestedConfiguration,
              (unsigned long)scalingSqrtN, (unsigned long)scalingTau,
              scalingLeafRepresentation,
              scalingFvrStage1Criterion, scalingFvrStage2Criterion,
              scalingDismantleStage2Criterion, scalingFinalizationCriterion,
              scalingIncrementalShrink ? 1 : 0,
              scalingIncrementalDissolve ? 1 : 0,
              scalingPersistentSearch ? 1 : 0, scalingEventSearch ? 1 : 0,
              scalingLazyDualClock ? 1 : 0,
              scalingEventBatchSameTime ? 1 : 0,
              scalingEventCompactInitialization ? 1 : 0,
              scalingEventBoundaryReuse ? 1 : 0,
              scalingCriterion2FreeFrontier ? 1 : 0,
              scalingSplitFindminBoundary ? 1 : 0,
              scalingSplitFindminScanOracle ? 1 : 0,
              scalingSplitFindminGrow ? 1 : 0,
              scalingSplitFindminGrowScanOracle ? 1 : 0,
              scalingSharedSplitFindminTopology ? 1 : 0,
              scalingSplitFindminListFrontier ? 1 : 0,
              scalingDetailedScheduleTiming ? 1 : 0,
              scalingDynamicShell ? 1 : 0,
              scalingSensitivityParameters ? 1 : 0);
  const char* resultStatus = accepted ? "success" :
    (!solverOk ? solverStatus : (valid ? "nonperfect" : "invalid-result"));
  std::printf("RESULT status=%s solver-status=%s wall=%.6f card=%lu mutual=%lu "
              "matched-input-edges=%lu asymmetric=%lu perfect=%d "
              "exact-weight=%lld engine-weight=%.17g weight-agrees=%d\n",
              resultStatus, solverStatus,
              seconds, (unsigned long)engineCard, (unsigned long)mutual,
              (unsigned long)matchedInputEdges, (unsigned long)asymmetric,
              perfect ? 1 : 0, exactWeight, (double)engineWeight,
              weightAgrees ? 1 : 0);
  std::fflush(stdout);
  return accepted ? 0 : 3;
}
