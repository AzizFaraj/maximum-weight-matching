// Small exact probe of the standard MWM-to-MWPM two-copy reduction.
//
// Original graph: a four-vertex weighted star.  It has no perfect matching,
// and its unrestricted maximum-weight matching is the single weight-10 edge.
// The reduced graph has eight vertices, ten edges, and MWPM weight 20.
#include <cstdio>
#include <vector>

#include "Types.h"
#include "Graph.h"
#include "MatchingEngine.h"
#include "ScalingMatching.h"

using namespace Matchbox;

struct Result {
  const char* algorithm;
  bool ok;
  Size card;
  Val weight;
};

static bool validMatching(const std::vector<Size>& mate, Size n, Size card) {
  Size counted = 0;
  for (Size v = 0; v < n; ++v) {
    if (mate[v] == cNullItm) continue;
    if (mate[v] >= n || mate[mate[v]] != v) return false;
    if (v < mate[v]) ++counted;
  }
  return counted == card;
}

int main() {
  const Size originalN = 4;
  const Size reducedN = 2 * originalN;
  const Size originalU[] = {0, 0, 0};
  const Size originalV[] = {1, 2, 3};
  const long long originalW[] = {10, 9, 8};

  std::vector<Size> u, v;
  std::vector<Val> w;
  for (Size e = 0; e < 3; ++e) {
    u.push_back(originalU[e]);
    v.push_back(originalV[e]);
    w.push_back((Val)originalW[e]);
    u.push_back(originalU[e] + originalN);
    v.push_back(originalV[e] + originalN);
    w.push_back((Val)originalW[e]);
  }
  for (Size x = 0; x < originalN; ++x) {
    u.push_back(x);
    v.push_back(x + originalN);
    w.push_back((Val)0);
  }

  Graph reduced(reducedN, (Size)u.size(), u, v, w);
  std::printf("REDUCTION original_n=%lu original_m=3 reduced_n=%lu reduced_m=%lu "
              "expected_mwpm_card=4 expected_mwpm_weight=20\n",
              (unsigned long)originalN, (unsigned long)reducedN,
              (unsigned long)u.size());

  std::vector<Size> mate;
  std::vector<Result> results;
  MatchingEngine classics;
  Size card = 0;
  Val weight = 0;

  classics.ComputeMaxEdgWghtMatchingEdmondsFastV2(
      reduced, &mate, &card, &weight);
  results.push_back((Result){"Edmonds-native-MWM",
      validMatching(mate, reducedN, card), card, weight});

  classics.ComputeMaxEdgWghtMatchingGabow(reduced, &mate, &card, &weight);
  results.push_back((Result){"Gabow-native-MWM",
      validMatching(mate, reducedN, card), card, weight});

  ScalingMatching scaling;
  bool scalingOk = scaling.computeHybrid(reduced, &mate, &card, &weight);
  results.push_back((Result){"Scaling-MWPM",
      scalingOk && validMatching(mate, reducedN, card), card, weight});

  for (Size i = 0; i < results.size(); ++i) {
    const bool perfect = results[i].card * 2 == reducedN;
    std::printf("RESULT algorithm=%s valid=%d perfect=%d card=%lu weight=%.0f\n",
                results[i].algorithm, results[i].ok ? 1 : 0,
                perfect ? 1 : 0, (unsigned long)results[i].card,
                results[i].weight);
  }

  const Result& scalingResult = results.back();
  bool zeroWeightReductionOk = scalingResult.ok && scalingResult.card == 4 &&
                               scalingResult.weight == (Val)20;

  // Our current Edmonds and Gabow entry points optimize unrestricted weight,
  // not perfect-matching weight.  For a directly runnable common-input probe,
  // add the same cardinality-dominating constant to every reduced edge.  This
  // does not alter the ordering of perfect matchings.  Here C=55 is one more
  // than the sum of all positive edge weights in the reduced graph.
  const Val cardinalityBonus = (Val)55;
  std::vector<Val> commonWeights(w.size());
  for (Size e = 0; e < w.size(); ++e)
    commonWeights[e] = w[e] + cardinalityBonus;
  Graph common(reducedN, (Size)u.size(), u, v, commonWeights);
  std::vector<Result> commonResults;

  classics.ComputeMaxEdgWghtMatchingEdmondsFastV2(
      common, &mate, &card, &weight);
  commonResults.push_back((Result){"Edmonds-common-input",
      validMatching(mate, reducedN, card), card, weight});

  classics.ComputeMaxEdgWghtMatchingGabow(common, &mate, &card, &weight);
  commonResults.push_back((Result){"Gabow-common-input",
      validMatching(mate, reducedN, card), card, weight});

  ScalingMatching commonScaling;
  scalingOk = commonScaling.computeHybrid(common, &mate, &card, &weight);
  commonResults.push_back((Result){"Scaling-common-input",
      scalingOk && validMatching(mate, reducedN, card), card, weight});

  bool commonInputOk = true;
  for (Size i = 0; i < commonResults.size(); ++i) {
    const bool perfect = commonResults[i].card * 2 == reducedN;
    const Val decodedOriginalWeight =
        (commonResults[i].weight - cardinalityBonus * commonResults[i].card) / 2;
    std::printf("COMMON algorithm=%s valid=%d perfect=%d card=%lu "
                "encoded_weight=%.0f decoded_original_mwm_weight=%.0f\n",
                commonResults[i].algorithm, commonResults[i].ok ? 1 : 0,
                perfect ? 1 : 0, (unsigned long)commonResults[i].card,
                commonResults[i].weight, decodedOriginalWeight);
    commonInputOk = commonInputOk && commonResults[i].ok && perfect &&
                    commonResults[i].weight == (Val)240 &&
                    decodedOriginalWeight == (Val)10;
  }

  return zeroWeightReductionOk && commonInputOk ? 0 : 1;
}
