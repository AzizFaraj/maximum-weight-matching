// Regression for the scale-count endpoint when max((n/2+1) * w) is a power
// of two.  On this six-cycle, the two alternating perfect matchings have
// weights 6 and 3.  With n=6 and max input weight 2, the extended maximum is
// (6/2+1)*2 = 8.  Omitting bit 3 reverses the working-weight preference.
#include <cstdio>
#include <vector>

#include "Types.h"
#include "Graph.h"
#include "ScalingMatching.h"

using namespace Matchbox;

typedef bool (ScalingMatching::*ScalingEntry)(
  const Graph&, std::vector<Size>*, Size*, Val*);

static Graph makeAlternatingCycle() {
  const Size sourceData[] = {0, 1, 2, 3, 4, 5};
  const Size targetData[] = {1, 2, 3, 4, 5, 0};
  const Val weightData[] = {2, 1, 2, 1, 2, 1};
  std::vector<Size> source(sourceData, sourceData + 6);
  std::vector<Size> target(targetData, targetData + 6);
  std::vector<Val> weight(weightData, weightData + 6);
  return Graph(6, 6, source, target, weight);
}

static bool checkEntry(const char* name, ScalingEntry entry, Size expectedTau) {
  Graph graph = makeAlternatingCycle();
  ScalingMatching scaling;
  std::vector<Size> mate;
  Size card = 0;
  Val weight = 0;
  const bool solved = (scaling.*entry)(graph, &mate, &card, &weight);
  const bool preferredHeavyEdges = mate.size() == 6 &&
    mate[0] == 1 && mate[1] == 0 &&
    mate[2] == 3 && mate[3] == 2 &&
    mate[4] == 5 && mate[5] == 4;
  const ScalingMatching::RunConfiguration config = scaling.runConfiguration();
  const bool parametersPreserved = config.sqrtN == 3 &&
                                   config.tau == expectedTau;
  const bool ok = solved &&
    scaling.lastStatus() == ScalingMatching::STATUS_SUCCESS &&
    card == 3 && (long long)(weight + (Val)0.5) == 6 &&
    preferredHeavyEdges && parametersPreserved;
  std::printf(
    "%s solved=%d status=%s card=%lu weight=%.0f heavy-perfect=%d "
    "sqrtN=%lu tau=%lu %s\n",
    name, solved ? 1 : 0,
    ScalingMatching::statusName(scaling.lastStatus()),
    (unsigned long)card, weight, preferredHeavyEdges ? 1 : 0,
    (unsigned long)config.sqrtN, (unsigned long)config.tau,
    ok ? "PASS" : "FAIL");
  return ok;
}

int main() {
  bool ok = true;
  ok = checkEntry("hybrid", &ScalingMatching::computeHybrid, 9) && ok;
  ok = checkEntry("hybrid-reference",
                  &ScalingMatching::computeHybridReference, 9) && ok;
  ok = checkEntry("liquidationist",
                  &ScalingMatching::computeLiquidationist, 3) && ok;
  std::printf("SCALING EXACT-POWER ENDPOINT GATE: %s\n",
              ok ? "PASS" : "FAIL");
  return ok ? 0 : 1;
}
