#include <cstdio>
#include <vector>
#include "Types.h"
#include "Graph.h"
#include "ScalingMatching.h"

using namespace Matchbox;

typedef bool (ScalingMatching::*ScalingEntry)(
  const Graph&, std::vector<Size>*, Size*, Val*);

static Graph makeGraph(Val firstWeight) {
  std::vector<Size> source, target;
  std::vector<Val> weight;
  const Size s[] = {0, 1, 2, 3, 0, 1};
  const Size t[] = {1, 2, 3, 0, 2, 3};
  const Val w[] = {firstWeight, 7, 6, 5, 4, 3};
  for (Size i = 0; i < 6; ++i) {
    source.push_back(s[i]);
    target.push_back(t[i]);
    weight.push_back(w[i]);
  }
  return Graph(4, 6, source, target, weight);
}

static bool checkEntry(const char* name, ScalingEntry entry) {
  Graph integral = makeGraph(8);
  Graph fractional = makeGraph((Val)7.5);

  ScalingMatching accepted;
  std::vector<Size> mate;
  Size card = 99;
  Val total = (Val)-1;
  bool integralOk = (accepted.*entry)(integral, &mate, &card, &total);
  bool acceptedOk = integralOk && card == 2 && mate.size() == 4 &&
    accepted.lastStatus() == ScalingMatching::STATUS_SUCCESS;

  ScalingMatching rejected;
  mate.assign(9, 0);
  card = 99;
  total = (Val)123;
  bool fractionalOk = (rejected.*entry)(fractional, &mate, &card, &total);
  bool rejectedOk = !fractionalOk &&
    rejected.lastStatus() ==
      ScalingMatching::STATUS_UNSUPPORTED_OR_OVERFLOWING_INPUT &&
    mate.size() == 4 && card == 0 && total == 0;

  std::printf("%s integral=%s fractional-rejected=%s\n", name,
    acceptedOk ? "PASS" : "FAIL", rejectedOk ? "PASS" : "FAIL");
  return acceptedOk && rejectedOk;
}

int main() {
  bool ok = true;
  ok = checkEntry("hybrid", &ScalingMatching::computeHybrid) && ok;
  ok = checkEntry("hybrid-reference",
    &ScalingMatching::computeHybridReference) && ok;
  ok = checkEntry("liquidationist",
    &ScalingMatching::computeLiquidationist) && ok;
  std::printf("SCALING INPUT-DOMAIN GATE: %s\n", ok ? "PASS" : "FAIL");
  return ok ? 0 : 1;
}
