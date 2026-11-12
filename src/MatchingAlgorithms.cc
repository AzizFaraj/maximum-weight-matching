#include "MatchingAlgorithms.h"

// Freeze the validated shared-order/ownership policy for the three-function
// adapter even when it is compiled without the repository's CMake target.
#ifndef MBX_SHARED_D3_TIE_ORDER
#define MBX_SHARED_D3_TIE_ORDER 1
#endif
#ifndef MBX_GMG_D3_OWNERSHIP_EXPERIMENT
#define MBX_GMG_D3_OWNERSHIP_EXPERIMENT 1
#endif
#include "classic/MatchingEngine.h"

// This translation unit permanently selects the practical Scaling
// configuration used for the paper's headline experiments.
#ifndef SCALING_S5_FULL_KEY_DRIVER_DISABLE
#define SCALING_S5_FULL_KEY_DRIVER_DISABLE 1
#endif
#ifndef SCALING_EVENT_INITIAL_DISSOLUTION_FRONTIER
#define SCALING_EVENT_INITIAL_DISSOLUTION_FRONTIER 1
#endif
#ifndef SCALING_ZERO_ROOT_CLEANUP_FRONTIER
#define SCALING_ZERO_ROOT_CLEANUP_FRONTIER 1
#endif
#include "scaling/ScalingMatching.h"

namespace Matchbox {

MatchingResult Edmonds(const Graph& graph) {
  MatchingResult result;
  MatchingEngine engine;
  engine.ComputeMaxEdgWghtMatchingEdmondsFastV2(
      graph, &result.mate, &result.cardinality, &result.weight);
  return result;
}

MatchingResult Gabow(const Graph& graph) {
  MatchingResult result;
  MatchingEngine engine;
  engine.ComputeMaxEdgWghtMatchingGabow(
      graph, &result.mate, &result.cardinality, &result.weight);
  return result;
}

MatchingResult Scaling(const Graph& graph) {
  MatchingResult result;
  ScalingMatching scaling;
  scaling.computeHybrid(
      graph, &result.mate, &result.cardinality, &result.weight);

  switch (scaling.lastStatus()) {
    case ScalingMatching::STATUS_SUCCESS:
      result.status = SolveStatus::kSuccess;
      break;
    case ScalingMatching::STATUS_SUCCESS_EXPERIMENTAL_PARTIAL:
      result.status = SolveStatus::kExperimentalPartial;
      break;
    case ScalingMatching::STATUS_NO_PERFECT_MATCHING:
      result.status = SolveStatus::kNoPerfectMatching;
      break;
    case ScalingMatching::STATUS_UNSUPPORTED_OR_OVERFLOWING_INPUT:
      result.status = SolveStatus::kUnsupportedOrOverflowingInput;
      break;
    case ScalingMatching::STATUS_INVARIANT_VIOLATION:
      result.status = SolveStatus::kInvariantViolation;
      break;
    case ScalingMatching::STATUS_UNEXPECTED_NO_PROGRESS:
      result.status = SolveStatus::kUnexpectedNoProgress;
      break;
  }
  return result;
}

const char* SolveStatusName(SolveStatus status) {
  switch (status) {
    case SolveStatus::kSuccess:
      return "success";
    case SolveStatus::kExperimentalPartial:
      return "experimental-partial";
    case SolveStatus::kNoPerfectMatching:
      return "no-perfect-matching";
    case SolveStatus::kUnsupportedOrOverflowingInput:
      return "unsupported-or-overflowing-input";
    case SolveStatus::kInvariantViolation:
      return "invariant-violation";
    case SolveStatus::kUnexpectedNoProgress:
      return "unexpected-no-progress";
  }
  return "unknown";
}

}  // namespace Matchbox
