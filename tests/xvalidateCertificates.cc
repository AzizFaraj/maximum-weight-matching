// Deterministically search for small graphs where FastV2's returned matching or
// retained primal/dual certificate disagrees between the GMG and Gabow backends.
#include <cstdio>
#include <cstdlib>
#include <vector>
#include "Types.h"
#include "Graph.h"
#include "MatchingEngine.h"

using namespace Matchbox;

static unsigned long gSeed = 12345;
static unsigned long nextRandom() {
  gSeed = gSeed * 1103515245UL + 12345UL;
  return (gSeed >> 16) & 0x7fffffffUL;
}

static bool writeGraph(const char* path, Size n, const std::vector<Size>& u,
                       const std::vector<Size>& v, const std::vector<Val>& w) {
  FILE* out = std::fopen(path, "w");
  if (!out) return false;
  std::fprintf(out, "%%%%MatrixMarket matrix coordinate integer symmetric\n");
  std::fprintf(out, "%lu %lu %lu\n", (unsigned long)n, (unsigned long)n,
               (unsigned long)u.size());
  for (Size i = 0; i < u.size(); ++i)
    std::fprintf(out, "%lu %lu %.0f\n", (unsigned long)(u[i] + 1),
                 (unsigned long)(v[i] + 1), (double)w[i]);
  return std::fclose(out) == 0;
}

int main(int argc, char** argv) {
  const Size numGraphs = argc > 1 ? (Size)std::strtoul(argv[1], 0, 10) : 10000;
  const Size maxN = argc > 2 ? (Size)std::strtoul(argv[2], 0, 10) : 80;
  const char* failurePath = argc > 3 ? argv[3] : "certificate_failure.mtx";
  const Size expectedDegree = argc > 4 ? (Size)std::strtoul(argv[4], 0, 10) : 4;
  const Size weightRange = argc > 5 ? (Size)std::strtoul(argv[5], 0, 10) : 1000;
  if (maxN < 4) return 2;
  if (expectedDegree == 0 || weightRange == 0) return 2;

  for (Size g = 0; g < numGraphs; ++g) {
    const Size n = 4 + nextRandom() % (maxN - 3);
    std::vector<Size> u, v;
    std::vector<Val> w;
    for (Size a = 0; a < n; ++a) {
      for (Size b = a + 1; b < n; ++b) {
        if (nextRandom() % n < expectedDegree) {
          u.push_back(a); v.push_back(b);
          w.push_back((Val)(1 + nextRandom() % weightRange));
        }
      }
    }
    if (u.empty()) continue;
    Graph graph(n, u.size(), u, v, w);
    std::vector<Size> mateGmg, mateGab;
    Size cardGmg = 0, cardGab = 0, certCard = 0;
    Val weightGmg = 0, weightGab = 0, certWeight = 0, dual = 0, slack = 0;

    FastMatchingEngineV2 gmg(FastMatchingEngineV2::kGMG1986);
    gmg.compute(graph, &mateGmg, &cardGmg, &weightGmg);
    const bool certGmg = gmg.verifyCertificate(
      mateGmg, &certCard, &certWeight, &dual, &slack);

    FastMatchingEngineV2 gab(FastMatchingEngineV2::kGabow2018);
    gab.compute(graph, &mateGab, &cardGab, &weightGab);
    const bool certGab = gab.verifyCertificate(
      mateGab, &certCard, &certWeight, &dual, &slack);

    if (!certGmg || !certGab || weightGmg != weightGab) {
      if (!writeGraph(failurePath, n, u, v, w)) return 2;
      std::printf("FAIL graph=%lu n=%lu m=%lu gmg-weight=%.17g "
                  "gabow-weight=%.17g gmg-cert=%s gabow-cert=%s "
                  "gabow-dual=%.17g gabow-min-slack=%.17g fixture=%s\n",
                  (unsigned long)g, (unsigned long)n, (unsigned long)u.size(),
                  (double)weightGmg, (double)weightGab,
                  certGmg ? "pass" : "FAIL", certGab ? "pass" : "FAIL",
                  (double)dual, (double)slack, failurePath);
      return 1;
    }
  }
  std::printf("PASS graphs=%lu max-n=%lu expected-degree=%lu weight-range=%lu\n",
              (unsigned long)numGraphs, (unsigned long)maxN,
              (unsigned long)expectedDegree, (unsigned long)weightRange);
  return 0;
}
