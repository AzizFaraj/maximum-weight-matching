#ifndef MATCHING_ALGORITHMS_H
#define MATCHING_ALGORITHMS_H

#include <vector>

#include "Types.h"

namespace Matchbox {

class Graph;

enum class SolveStatus {
  kSuccess,
  kExperimentalPartial,
  kNoPerfectMatching,
  kUnsupportedOrOverflowingInput,
  kInvariantViolation,
  kUnexpectedNoProgress
};

struct MatchingResult {
  std::vector<Size> mate;
  Size cardinality = 0;
  Val weight = 0;
  SolveStatus status = SolveStatus::kSuccess;
};

/// Exact unrestricted maximum-weight matching using the Edmonds backend.
MatchingResult Edmonds(const Graph& graph);

/// Exact unrestricted maximum-weight matching using the Gabow backend.
MatchingResult Gabow(const Graph& graph);

/// Exact maximum-weight perfect matching using practical Figure-6 Scaling.
MatchingResult Scaling(const Graph& graph);

const char* SolveStatusName(SolveStatus status);

}  // namespace Matchbox

#endif  // MATCHING_ALGORITHMS_H
