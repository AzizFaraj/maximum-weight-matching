// Canonical Figure-6 identity and scan-reference differential gate.
#define MATCHED_GROW_AUDIT
#include <cstdio>
#include <cstdlib>
#include <vector>
#include "Types.h"
#include "Graph.h"
#include "ScalingMatching.h"

using namespace Matchbox;


static unsigned long gSeed;
static unsigned long nextRand() {
  gSeed = gSeed * 1103515245UL + 12345UL;
  return (gSeed >> 16) & 0x7fffffffUL;
}

static void makeGraph(int n, int rep, std::vector<Size>* u,
                      std::vector<Size>* v, std::vector<Val>* w) {
  gSeed = 97UL + (unsigned long)n * 31UL + (unsigned long)rep * 101UL;
  u->clear(); v->clear(); w->clear();
  std::vector<std::vector<char> > used(n, std::vector<char>(n, 0));
  for (int i = 0; i < n; ++i) {
    int j = (i + 1) % n;
    if (i < j && !used[i][j]) {
      used[i][j] = used[j][i] = 1;
      u->push_back(i); v->push_back(j); w->push_back((Val)(1 + nextRand() % 200));
    }
  }
  for (int e = 0; e < 4 * n; ++e) {
    int a = (int)(nextRand() % n), b = (int)(nextRand() % n);
    if (a == b || used[a][b]) continue;
    used[a][b] = used[b][a] = 1;
    u->push_back(a); v->push_back(b); w->push_back((Val)(1 + nextRand() % 1000));
  }
}

int main(int argc, char** argv) {
  static char envReweight[] = "MBX_REWEIGHT=1";
  static char envEvent[] = "MBX_EVENTSEARCH=1";
  static char envAdjust[] = "MBX_MAXADJ=999";
  putenv(envReweight); putenv(envEvent); putenv(envAdjust);
  int ok = 0, fail = 0;
  unsigned long long referenceIneligibleMatchedGrows = 0;
  unsigned long long parkedVertices = 0, nonTauCorrections = 0;
  unsigned long long tightnessFailures = 0;
  for (int n = 4; n <= 30; n += 2) for (int rep = 0; rep < 30; ++rep) {
    if ((argc > 1 && std::atoi(argv[1]) != n) ||
        (argc > 2 && std::atoi(argv[2]) != rep)) continue;
    std::vector<Size> u, v; std::vector<Val> w;
    makeGraph(n, rep, &u, &v, &w);
    Graph g((Size)n, (Size)u.size(), u, v, w);
    ScalingMatching fast, ref;
    std::vector<Size> mf, mr;
    Size cf = 0, cr = 0; Val wf = 0, wr = 0;
    bool af = fast.computeHybrid(g, &mf, &cf, &wf);
    bool ar = ref.computeHybridReference(g, &mr, &cr, &wr);
    referenceIneligibleMatchedGrows +=
      ref.operationCounters().matchedGrowIneligible;
    parkedVertices += fast.operationCounters().perfectionParkedVertices;
    nonTauCorrections += fast.operationCounters().perfectionNonTauCorrections;
    tightnessFailures += fast.operationCounters().perfectionTightnessFailures;
    ScalingMatching::RunConfiguration fc = fast.runConfiguration();
    ScalingMatching::RunConfiguration rc = ref.runConfiguration();
    bool same = af == ar && cf == cr && (long long)(wf + .5) == (long long)(wr + .5);
    bool identity = fc.algorithm == ScalingMatching::ALGORITHM_HYBRID_FIGURE6 &&
                    rc.algorithm == ScalingMatching::ALGORITHM_HYBRID_REFERENCE &&
                    fc.eventSearch && !rc.eventSearch &&
                    std::string(fc.backend) == "integrated-event-materialized-shell" &&
                    std::string(rc.backend) == "scan-reference" &&
                    fc.eventSplitFindminBoundary &&
                    !fc.eventSplitFindminScanOracle &&
                    fc.eventSplitFindminGrow &&
                    !fc.eventSplitFindminGrowScanOracle &&
                    fc.eventSplitFindminListFrontier &&
                    !rc.eventSplitFindminBoundary &&
                    !rc.eventSplitFindminGrow &&
                    fc.fvrStage1Criterion == 2 && fc.fvrStage2Criterion == 3 &&
                    fc.dismantleStage2Criterion == 1 &&
                    fc.finalizationCriterion == 3 &&
                    rc.fvrStage1Criterion == 2 && rc.fvrStage2Criterion == 3 &&
                    rc.dismantleStage2Criterion == 1 &&
                    rc.finalizationCriterion == 3 &&
                    fc.tau == rc.tau && fc.sqrtN == rc.sqrtN &&
                    fast.lastStatus() == ref.lastStatus();
    if (same && identity) ++ok;
    else {
      ++fail;
      if (fail <= 8)
        std::printf("FAIL n=%d rep=%d af/ar=%d/%d card=%lu/%lu weight=%.0f/%.0f status=%s/%s\n",
          n, rep, af ? 1 : 0, ar ? 1 : 0, (unsigned long)cf,
          (unsigned long)cr, wf, wr,
          ScalingMatching::statusName(fast.lastStatus()),
          ScalingMatching::statusName(ref.lastStatus()));
    }
  }
  bool activated = referenceIneligibleMatchedGrows > 0;
  bool perfectionActivated = parkedVertices > 0 && nonTauCorrections > 0 &&
                             tightnessFailures == 0;
  std::printf("HYBRID IDENTITY/REFERENCE GATE: oks=%d fails=%d deferred-activations=%llu parked=%llu non-tau-corrections=%llu tightness-failures=%llu\n",
              ok, fail, referenceIneligibleMatchedGrows, parkedVertices,
              nonTauCorrections, tightnessFailures);
  return fail || !activated || !perfectionActivated ? 1 : 0;
}
