#ifndef SCALING_MATCHING_H
#define SCALING_MATCHING_H

// The linked leaf-order representation is the production default.  The old
// per-blossom materialized leaf vectors have quadratic retained storage on
// deeply nested Hybrid instances (observed as a 20+ GiB working set at
// n=500,000).  Keep that representation only for explicit compatibility and
// differential experiments.
#if defined(SCALING_LEAF_NOVEC) && defined(SCALING_LEAF_VECTOR_LEGACY)
#error "Choose linked or legacy-vector ScalingMatching leaves, not both"
#elif !defined(SCALING_LEAF_NOVEC) && !defined(SCALING_LEAF_VECTOR_LEGACY)
#define SCALING_LEAF_NOVEC 1
#endif

// Select the canonical secondary split-findmin channel before including its
// implementation. SplitFindmin chooses the dual-key concrete queue from this
// flag; classical ConcatenableQueue clients retain their single-key layout.
#if !defined(SCALING_S5_FULL_KEY_DRIVER_DISABLE) && \
    !defined(EVENT_S5_FULL_KEY_DRIVER_AUDIT) && \
    !defined(EVENT_S5_FULL_KEY_DRIVER_NO_ORACLE)
#define SCALING_S5_FULL_KEY_DRIVER_CANONICAL
#define EVENT_S5_FULL_KEY_DRIVER_NO_ORACLE
#endif
#if defined(EVENT_S5_FULL_KEY_DRIVER_AUDIT) || \
    defined(EVENT_S5_FULL_KEY_DRIVER_NO_ORACLE)
#ifndef EVENT_S5_TARGET_KEY_PERSIST_FRESH
#define EVENT_S5_TARGET_KEY_PERSIST_FRESH
#endif
#ifndef EVENT_S5_SHARED_SFM_TOPOLOGY
#define EVENT_S5_SHARED_SFM_TOPOLOGY
#endif
#endif
#if defined(SCALING_S5_FULL_KEY_DRIVER_CANONICAL) && \
    !defined(EVENT_S5_FULL_KEY_DRIVER_VERIFY)
#define EVENT_S5_TARGET_KEY_SKIP_EXACT_AUDIT
#endif
#if defined(EVENT_S5_SHARED_SFM_TOPOLOGY) && \
    !defined(SCALING_S5_SINGLETON_FAST_DISABLE)
// Most Section-5 non-outer lists are singleton roots.  Keep their secondary
// key in flat storage and enter the shared AVL only while an inherited/current
// blossom genuinely owns a multi-vertex interval.  The topology only splits,
// so a list's transition to this representation is permanent for the search.
#define EVENT_S5_SINGLETON_TARGET_FAST
#endif

// ============================================================================
// ScalingMatching — Duan-Pettie-Su (SODA 2017) LIQUIDATIONIST scaling algorithm
// for maximum-weight PERFECT matching on general graphs with non-negative
// integer edge weights.
//
// STANDALONE engine (Path A): self-contained, does NOT touch matchbox's GMG /
// Gabow engines. Those are used only as a correctness ORACLE in tests.
//
// Status: compute() is the verified single-scale exact blossom solver.
// computeHybrid() is the canonical Figure-6/Section-4.2 scaling entry point;
// computeHybridReference() is its scan-backed semantic oracle. The integrated
// Section-5 backend is practical/adapted, so Theorem 4.3's exact bound is not
// claimed. See docs/Hybrid_Current_Conformance.md.
//
// Design notes:
//  * Integer duals/weights. Internally we use `long long` (LL): the algorithm
//    forms an extended weight wbar(e) = (n/2 + 1) * what(e). For n ~ 1-2M and
//    what <= ~1000 this reaches ~5e8 per edge and dual sums reach ~1e15 — well
//    past 32-bit. (The existing LEMON driver also used long long.)
//  * C++98 / -ansi compatible: no auto, no range-for, no move semantics.
//  * Maintains Relaxed Complementary Slackness (paper Property 2): even integer
//    weights w, integer y, even z >= 0, near-domination yz(e) >= w(e)-2,
//    near-tightness yz(e) <= w(e) on matched/blossom edges.
// ============================================================================

#include <vector>
#include <cassert>
#include <cstddef>
#include <algorithm>   // std::reverse, std::swap
#include <utility>     // std::move
#include <limits>
#include <cmath>
#if defined(SCDBG) || defined(EVENTTRACE)
#include <cstdio>
#endif
#ifdef HYBPROF
#include <cmath>       // pow (Lemma 4.10 budget m(B) n(B)^{3/4})
#include <cstdio>
#endif

#include <time.h>
#include "Types.h"   // Matchbox::Size, Matchbox::Val, cNullItm
#include "Graph.h"   // Matchbox::Graph
#include "BucketPQ.h" // Matchbox::BucketPQ (S5 event scheduling)
#include "BinaryEventHeap.h" // unbounded practical PQSearch event queue
#include "SplitFindmin.h" // AVL-backed inner-blossom event ownership

// Canonical Section-5 grow scheduling.  The maintained split-findmin list
// minima drive Criteria 1--3; the former per-edge OUTER--FREE queue can be
// restored for controlled A/B tests with SCALING_S5_FULL_KEY_DRIVER_DISABLE.
// Define EVENT_S5_FULL_KEY_DRIVER_VERIFY to retain the O(m) end-of-search scan
// audit in a canonical build.
#if defined(EVENT_S5_FULL_KEY_DRIVER_NO_ORACLE) && \
    !defined(EVENT_S5_FULL_KEY_DRIVER_AUDIT)
#define EVENT_S5_FULL_KEY_DRIVER_AUDIT
#endif
#if defined(EVENT_S5_FULL_KEY_DRIVER_AUDIT) && \
    !defined(EVENT_S5_TARGET_KEY_DRIVER_AUDIT)
#define EVENT_S5_TARGET_KEY_DRIVER_AUDIT
#endif
#if defined(EVENT_S5_FULL_KEY_DRIVER_NO_ORACLE) && \
    !defined(EVENT_S5_TARGET_KEY_DRIVER_NO_ORACLE)
#define EVENT_S5_TARGET_KEY_DRIVER_NO_ORACLE
#endif
#if defined(EVENT_S5_TARGET_KEY_DRIVER_NO_ORACLE) && \
    !defined(EVENT_S5_TARGET_KEY_DRIVER_AUDIT)
#define EVENT_S5_TARGET_KEY_DRIVER_AUDIT
#endif
#if defined(EVENT_S5_TARGET_KEY_DRIVER_AUDIT) && \
    !defined(EVENT_S5_TARGET_KEY_SHADOW)
#define EVENT_S5_TARGET_KEY_SHADOW
#endif
#if defined(EVENT_S5_SHARED_SFM_TOPOLOGY) && \
    !defined(EVENT_S5_TARGET_KEY_SHADOW)
#error EVENT_S5_SHARED_SFM_TOPOLOGY requires the Section-5 target-key state
#endif

// Integrated Section-5 boundary-frontier candidate.  The scan-free consumer
// is canonical only after every current INNER root has been seeded at fresh
// search initialization; later role changes maintain the same keys online.
#if !defined(SCALING_EVENT_SFM_DISABLE) && \
    !defined(EVENT_SFM_BOUNDARY_FRONTIER) && \
    !defined(EVENT_SFM_BOUNDARY_FRONTIER_FAST)
// Defining SCALING_EVENT_SFM_DISABLE restores the counted dissolution-boundary
// scan; defining EVENT_SFM_BOUNDARY_FRONTIER explicitly retains that scan as an
// in-process oracle for development builds.
#define EVENT_SFM_BOUNDARY_FRONTIER_FAST
#endif
#ifdef EVENT_SFM_BOUNDARY_FRONTIER_FAST
#ifndef EVENT_SFM_BOUNDARY_FRONTIER
#define EVENT_SFM_BOUNDARY_FRONTIER
#endif
#ifndef EVENT_SFM_REUSE
#define EVENT_SFM_REUSE
#endif
#ifndef EVENT_SFM_LAZY_ACTIVE
#define EVENT_SFM_LAZY_ACTIVE
#endif
#ifndef EVENT_SFM_NO_ORACLE
#define EVENT_SFM_NO_ORACLE
#endif
#endif
#ifdef EVENT_SFM_BOUNDARY_FRONTIER
#ifndef EVENT_SFM_SHADOW
#define EVENT_SFM_SHADOW
#endif
#ifndef EVENT_SFM_LAZY_REPAIR
#define EVENT_SFM_LAZY_REPAIR
#endif
#ifndef EVENT_SFM_CONSUME
#define EVENT_SFM_CONSUME
#endif
#ifndef EVENT_SFM_BOUNDARY_FRONTIER_FAST
#ifndef EVENT_SFM_TOPOLOGY_AUDIT
#define EVENT_SFM_TOPOLOGY_AUDIT
#endif
#endif
#endif
#if defined(EVENT_SFM_CONSUME) && !defined(EVENT_SFM_SHADOW)
#error EVENT_SFM_CONSUME requires EVENT_SFM_SHADOW validation state
#endif
#include "ScalingFastMWPM.h" // the fast UF+lazy-dual+event-PQ engine (Layer 2 subroutine)
#include "ScalingS5Engine.h" // the no-relabel S5 engine: O(m sqrt n)/search (Layer 2, fast path)

namespace Matchbox {

static inline unsigned long long scalingMonotonicNanoseconds() {
  struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t);
  return (unsigned long long)t.tv_sec * 1000000000ULL +
         (unsigned long long)t.tv_nsec;
}

struct ScalingNanosecondGuard {
  explicit ScalingNanosecondGuard(unsigned long long* counter)
    : mCounter(counter), mStart(scalingMonotonicNanoseconds()) {}
  ~ScalingNanosecondGuard() {
    if (mCounter) *mCounter += scalingMonotonicNanoseconds() - mStart;
  }
  unsigned long long* mCounter;
  unsigned long long mStart;
};

#ifdef OPCNT
// Operation counters to PROVE where the cold-scale time goes (not infer). Per-scale dump.
long long gOpEdgeScans = 0;   // edges examined in the grow inner loop
long long gOpOuterSteps = 0;  // mBParent-walk steps in outerOf (the O(depth) cost)
long long gOpSlackSteps = 0;  // mBParent-walk steps in slackZ
long long gOpUfSteps = 0;     // ufFind steps
long long gOpGrowIters = 0;   // grow outer-loop iterations
long long gOpUfMiss = 0;      // edges where UF said "different" but outerOf says "same" (false-neg)
long long gOpDissolves = 0;   // total blossom dissolutions (distinct events)
long long gDisIn42 = 0, gDisFvr = 0;   // dissolve split: §4.2 small-blossom vs FVR search
long long gOpShrinks = 0;     // total shrinks (blossom formations)
long long gOpUfSlowPath = 0;  // exact O(depth) laminar resolution when UF map is unavailable
long long gOpEligible = 0;    // PROBE: edges that PASSED the eligibility test in the grow (vs gOpEdgeScans = all scanned)
#endif
#ifdef LOOPTIME
#include <time.h>
// accurate wall-time attribution of the per-call O(n) loops in edmondsSearchEx.
// (defined in the test TU to avoid header multiple-definition; declared extern here.)
extern double gLtFlatten, gLtReset, gLtReroot, gLtWorklist, gLtEventScan, gLtDualApply, gLtGrow;
static inline double ltNow(){ struct timespec t; clock_gettime(CLOCK_MONOTONIC,&t); return t.tv_sec + t.tv_nsec/1e9; }
#define LT_START double _lt0 = ltNow()
#define LT_END(acc) do { (acc) += ltNow() - _lt0; } while(0)
#else
#define LT_START
#define LT_END(acc)
#endif

#ifdef HYBPROF
// HYBRID §4.2 instrumentation: confirm the Gabow dissolution path actually runs (vs falling
// through to liquidate), and that shell searches/translations execute. Proves the gate is
// EXERCISING the new code, not a false pass.
long long gHybDissolveCalls = 0;   // gabowDissolveSmallBlossom invocations
long long gHybShellSearches = 0;   // shellSearch invocations
long long gHybStage2Calls   = 0;   // dismantleStage2 invocations
long long gHybStage2Empty   = 0;   // attempt after path already fully dissolved
long long gHybStage2Free0   = 0;   // illegal stage-2 attempt with no free vertex in R
long long gHybStage2Free1   = 0;   // legal unique-omega attempt
long long gHybStage2FreeMany= 0;   // illegal attempt with >=2 free vertices in R
long long gHybStage1Stalls  = 0;   // stage 1 failed its legal-progress condition
long long gHybDPathEntries  = 0;
long long gHybDPathTopEmpty = 0;
long long gHybShellAugCalls = 0;
long long gHybShellNoAugCalls = 0;
long long gHybPostShellEmpty = 0;
long long gHybPostShellFree1 = 0;
long long gHybPostShellOther = 0;
long long gHybShellCUpdates = 0;
long long gHybShellDUpdates = 0;
long long gHybShellAbsorbs = 0;
long long gHybFlattenDuringDissolve = 0;
long long gHybFlattenInherited = 0;
long long gHybFlattenCurrent = 0;
long long gHybFlattenInheritedPos = 0;
long long gHybFlattenCurrentPos = 0;
long long gHybFlattenInheritedZCleared = 0;
long long gHybFlattenInheritedZBefore = 0;
long long gHybFreeParityViol = 0;
long long gHybCurrentBaseRotations = 0;
long long gHybCurrentMultiFree = 0;
long long gHybMaxFreeInSmall = 0;  // max #free interior vtxs in a small blossom at dissolution
// Complexity (Lemma 4.10: DismantlePath(R) = O(m(R) n(R)^{3/4}); whole dissolution over a small
// blossom B = O(m(B) n(B)^{3/4})). gHybDissolveAdjusts = total dual adjustments inside §4.2;
// gHybDissolveScans = total edge scans; gHybBudget = Sum over dissolved B of m(B)*n(B)^{0.75}.
// Ratio (adjusts+scans)/budget should stay BOUNDED (flat) as n grows => bound met.
long long gHybDissolveAdjusts = 0;
long long gHybDissolveScans   = 0;
double     gHybBudget         = 0.0;
bool       gHybInDissolve     = false;  // gate: count scans/adjusts only while inside §4.2
#endif

#ifdef SCPROF
static long long gProfIters = 0, gProfAug = 0, gProfAdj = 0, gProfShrink = 0;
static long long gProfCalls = 0, gProfRelabels = 0;  // edmondsSearchEx calls + O(n) relabel passes
static long long gRelShrink = 0, gRelDissolve = 0, gRelFirst = 0;  // relabel CAUSE: shrink / dissolve / first-iter(per-call entry)
long long gRelinkOK = 0, gRelinkBail = 0;            // relinkOddCycle success / rejected states
long long gRelabelShell = 0, gRelabelMain = 0;       // relabels inside ShellSearch vs main FVR
long long gRelinkB[8] = {0,0,0,0,0,0,0,0};           // bail reason histogram
// Split adjustment accounting to test the sqrt(n)-per-scale structure honestly:
//   gProfScaleAdj   = adjustments during the per-scale step-7 SearchOne loop (all scales)
//   gProfFinalAdj   = adjustments during Finalization (the post-last-scale Crit3 to perfection)
//   gProfFinalSrch  = number of Finalization searches (should be ~O(n/tau)=O(sqrt n) per the paper)
// gProfInFinal flips on when Finalization starts so the shared edmondsSearch routes its counts.
static long long gProfScaleAdj = 0, gProfFinalAdj = 0, gProfFinalSrch = 0;
static bool gProfInFinal = false;
#endif

// CRUMB: a no-op breadcrumb hook (kept as a marker at the entry of blossom sub-calls
// shrink/augmentFrom/rotateToBase/flipUpToRoot/flattenNode/liquidate). Currently expands to
// nothing; a future debug session can redefine it to record/abort on a runaway sub-call cycle
// (see edm-cost-and-lemon-comparison memory for the technique). Zero cost as-is.
#define CRUMB(op,a,b) do {} while(0)

// Local null sentinel for ScalingMatching indices. A free function avoids the
// C++98 ODR-definition requirement of an in-class `static const` member whose
// address is taken.
inline Size scNull() { return (Size)-1; }

class ScalingMatching {
public:
  typedef long long LL;                 // integer weights and duals

  enum ResultStatus {
    STATUS_SUCCESS = 0,
    STATUS_SUCCESS_EXPERIMENTAL_PARTIAL,
    STATUS_NO_PERFECT_MATCHING,
    STATUS_UNSUPPORTED_OR_OVERFLOWING_INPUT,
    STATUS_INVARIANT_VIOLATION,
    STATUS_UNEXPECTED_NO_PROGRESS
  };

  enum AlgorithmIdentity {
    ALGORITHM_UNSPECIFIED = 0,
    ALGORITHM_HYBRID_FIGURE6,
    ALGORITHM_HYBRID_REFERENCE,
    ALGORITHM_HYBRID_EVENT_EXPERIMENTAL,
    ALGORITHM_HYBRID_DROP_PERFECT_EXPERIMENTAL,
    ALGORITHM_HYBRID_SKIP_FINALIZATION_EXPERIMENTAL,
    ALGORITHM_HYBRID_LEGACY_EXPERIMENTAL,
    ALGORITHM_S5_EXPERIMENTAL,
    ALGORITHM_LIQUIDATIONIST_FIGURE5
  };

  struct RunConfiguration {
    AlgorithmIdentity algorithm;
    const char* backend;
    Size sqrtN, tau;
    // Name the paper sites explicitly: Hybrid has both a two-stage FVR and a
    // separate stage 2 inside DismantlePath.
    int fvrStage1Criterion, fvrStage2Criterion;
    int dismantleStage2Criterion, finalizationCriterion;
    bool incrementalShrink, incrementalDissolve, persistentSearch, eventSearch;
    bool lazyDualClock, eventBatchSameTime, eventCompactInitialization;
    bool eventBoundaryReuse, eventCriterion2FreeFrontier;
    bool eventSplitFindminBoundary, eventSplitFindminScanOracle;
    bool eventSplitFindminGrow, eventSplitFindminGrowScanOracle;
    bool eventSharedSplitFindminTopology;
    bool eventSplitFindminListFrontier;
    bool detailedScheduleTiming;
    bool dynamicShell;
    bool sensitivityParameters;
  };

  enum SearchRegion { REGION_ALL_REAL = 0, REGION_SHELL, REGION_VERTEX_SET };
  enum SearchStoppingRule {
    STOP_AFTER_AUGMENTATION = 0,
    STOP_AT_DUAL_BUDGET,
    STOP_AT_SHELL_EVENT,
    STOP_AT_PERFECT_OR_NO_PATH
  };
  enum SearchBackend { BACKEND_SCAN = 0, BACKEND_EVENT };

  struct SearchRequest {
    SearchRequest() : roots(NULL), regionVertices(NULL), criterion(1), region(REGION_ALL_REAL),
      dualBudget(-1), maxAdjustments(-1), maximalAugment(false),
      batchEventFreeTime(false), batchSameTimeAugment(false), persistAcrossAdjust(false),
      persistAcrossAugment(false), exhaustBudget(false),
      stoppingRule(STOP_AFTER_AUGMENTATION), shellTranslation(false) {}
    const std::vector<Size>* roots;
    const std::vector<Size>* regionVertices;
    int criterion;
    SearchRegion region;
    LL dualBudget;
    LL maxAdjustments;
    bool maximalAugment;
    bool batchEventFreeTime;
    // Event backend only: after its first augmentation, retain the queue and
    // consume further vertex-disjoint augmenting events at the SAME absolute
    // dual time.  It never crosses a dual adjustment with the changed matching.
    bool batchSameTimeAugment;
    bool persistAcrossAdjust;
    bool persistAcrossAugment;
    bool exhaustBudget;
    SearchStoppingRule stoppingRule;
    bool shellTranslation;
  };

  struct SearchResult {
    SearchResult() : status(STATUS_SUCCESS), augmented(false),
      dualUnits(0), reachedBudget(false), noProgress(false) {}
    ResultStatus status;
    bool augmented;
    LL dualUnits;
    bool reachedBudget;
    bool noProgress;
  };

  struct OperationCounters {
    OperationCounters() : edgeScans(0), queueInsertions(0), queuePops(0),
      bucketQueueInsertions(0), heapQueueInsertions(0),
      bucketQueuePops(0), heapQueuePops(0),
      splitFindminWork(0), augmentations(0), dualUnits(0),
      blossomFormations(0), blossomDissolutions(0), boundaryChanges(0),
      absorptionHalts(0), shellAugmentationHalts(0),
      shellOutermostDissolutionHalts(0),
      dismantleStage1AllDissolvedExits(0),
      dismantleStage1SingletonExits(0), dismantleStage2Calls(0),
      shellMcmNanoseconds(0), shellSearchNanoseconds(0),
      searchInitializations(0),
      shellSearchCalls(0), shellBoundaryConfigurations(0),
      shellMaskVertices(0), shellMaskNanoseconds(0),
      shellMcmPhaseNanoseconds(0), shellMcmCalls(0), shellMcmPasses(0),
      shellMcmAugmentingPhases(0), shellMcmCertificationPhases(0),
      shellMcmPhaseAugmentations(0), shellMcmMaxPhaseAugmentations(0),
      shellMcmInitialScanNanoseconds(0), shellMcmGreedyNanoseconds(0),
      shellMcmAugmentSearchNanoseconds(0),
      shellMcmCertificationNanoseconds(0),
      shellMcmRootFilterNanoseconds(0), shellMcmAuditNanoseconds(0),
      shellDualAdjustmentNanoseconds(0), shellDualVertexVisits(0),
      zeroRootCleanupNanoseconds(0), zeroRootCleanupCalls(0),
      zeroRootCleanupNodeVisits(0),
      dismantleNanoseconds(0), dismantleChainNanoseconds(0),
      dismantleFreeCountNanoseconds(0), dismantleShellNanoseconds(0),
      dismantleProgressNanoseconds(0), dismantleCalls(0),
      dismantleIterations(0),
      eventCriterion1Searches(0), eventCriterion2Searches(0),
      eventCriterion3Searches(0),
      eventSearchNanoseconds(0), eventSchedulingNanoseconds(0),
      eventForestInitializationNanoseconds(0), dissolutionNanoseconds(0),
      eventFreeDiscoveryNanoseconds(0), eventLabelResetNanoseconds(0),
      eventRootInitializationNanoseconds(0), eventSfmInitializationNanoseconds(0),
      eventTargetKeyInitializationNanoseconds(0), eventLazyInitializationNanoseconds(0),
      s5TargetInitFrontierNanoseconds(0), s5TargetInitTargetsNanoseconds(0),
      s5TargetInitOwnersNanoseconds(0), s5TargetInitTokensNanoseconds(0),
      s5TargetInitResizeNanoseconds(0), s5TargetInitFreshNanoseconds(0),
      s5TargetInitRootsNanoseconds(0), s5TargetInitTargetVisits(0),
      s5TargetInitOwnerVisits(0), s5TargetInitTokenVisits(0),
      s5TargetInitRootVisits(0),
      eventInitializationNodeVisits(0), eventInitializationNodesReset(0),
      eventInitializationDormantDeferrals(0),
      eventDissolutionResetCandidates(0),
      eventScaleResetCandidates(0),
      eventInitializationFreeCandidates(0), eventInitializationRoots(0),
      eventFreeCandidateFullScans(0), eventFreeCandidateCacheHits(0),
      eventFreeCandidateVerticesExamined(0),
      eventCriterion2AdjustmentOnlyReturns(0),
      eventCriterion2ZeroDualAugmentReturns(0),
      eventCriterion2AdjustmentToZeroAugmentPairs(0),
      eventBoundaryReuseHits(0), eventBoundaryReuseMisses(0),
      eventBoundaryEventsSaved(0), eventBoundaryMaxEventsSaved(0),
      dissolutionRelinkNanoseconds(0), blossomRotationNanoseconds(0),
      blossomRotationCalls(0), blossomRotationParentSteps(0),
      blossomRotationFrames(0), blossomRotationMaxDepth(0),
      blossomRotationNoopCalls(0),
      eventUsefulPops(0), eventStaleOwnerPops(0),
      eventStaleRolePops(0), eventInternalPops(0),
      eventIneligiblePops(0), eventControlPops(0),
      eventEdgeOffers(0), eventEdgeFreshOffers(0),
      eventEdgeTransitionOffers(0), eventEdgeShellRejects(0),
      eventEdgeInternalRejects(0), eventEdgeOddRejects(0),
      eventEdgeFreshDedupRejects(0), eventEdgeNoEventRejects(0),
      eventEdgeCandidates(0), eventEdgeImmediateCandidates(0),
      eventEdgeBudgetDrops(0), eventEdgeBoundaryCoalesced(0),
      eventEdgeSeedMarkersAvoided(0),
      eventOwnerMinOwners(0), eventOwnerMinPayloads(0),
      eventOwnerMinReplacements(0), eventOwnerMinLaterRejected(0),
      eventOwnerMinMaxPayloadsPerSearch(0),
      eventOwnerFutureMatchedGrow(0), eventOwnerFutureMatchedBlossom(0),
      eventOwnerFutureDissolve(0), eventOwnerFutureOtherControl(0),
      eventOwnerPayloadChecks(0), eventOwnerPayloadFailures(0),
      eventOwnerBoundaryPayloads(0), eventOwnerBoundaryExpired(0),
      eventOwnerBoundaryStaleOwner(0), eventOwnerBoundaryStaleRole(0),
      eventOwnerFutureMaxMatched(0), eventOwnerFutureMaxControl(0),
      eventOwnerBoundaryMatched(0), eventOwnerBoundaryMatchedExpired(0),
      eventOwnerBoundaryControl(0), eventOwnerBoundaryControlExpired(0),
      matchedGrowTransitions(0), matchedGrowIneligible(0),
      matchedGrowDelayOne(0), matchedGrowDelayTwo(0),
      matchedGrowDelayOther(0),
      algorithmFallbackAttempts(0), dissolutionRelinkFailures(0),
      malformedFreeRootAttempts(0),
      postInitializationFullReseeds(0), dissolutionBoundaryScans(0),
      dynamicDissolutionReconstructions(0),
      dynamicDissolutionReconstructionChildren(0),
      dynamicDissolutionReconstructionVertices(0),
      dynamicDissolutionReconstructionScans(0),
      dynamicDissolutionReconstructionNanoseconds(0),
      dissolutionTimestampTransfers(0), dissolutionRetainedBaseLeaves(0),
      dissolutionContainmentSteps(0), dualMaterializations(0),
      lazyDualVertexActivations(0), lazyDualNodeActivations(0),
      lazyDualVerticesMaterialized(0), lazyDualNodesMaterialized(0),
      eventBatchedAugmentations(0), sfmBoundaryQueries(0),
      eventInvariantChecks(0), eventInvariantFailures(0),
      sfmBoundaryExact(0), sfmBoundaryStale(0),
      sfmBoundaryTooLow(0), sfmBoundaryTooHigh(0),
      sfmBoundaryInvalidOwner(0), sfmBoundaryKeyDrift(0),
      sfmBoundaryInternalKey(0), sfmBoundaryUnexplainedLow(0),
      sfmBoundaryTopologyLeak(0), sfmBoundaryMissingSource(0),
      sfmBoundaryRepairs(0), sfmBoundaryRepairScans(0),
      sfmBoundaryRepairFailures(0),
      sfmInitializations(0),
      sfmMergeTopologyFailures(0), sfmSplitTopologyFailures(0),
      sfmInitialLists(0), sfmInitialPivots(0),
      sfmFormationTopologyNoops(0), sfmActualSplits(0),
      sfmGrowSetCalls(0), sfmGrowSetNoops(0),
      sfmGrowSetAncestorVisits(0),
      sfmGrowDecreaseCalls(0), sfmGrowDecreaseAncestorVisits(0),
      sfmGrowDecreaseImmediateStops(0), sfmGrowResetNodeVisits(0),
      s5TargetKeyInitializations(0), s5TargetKeyOffers(0),
      s5TargetKeyDecreases(0), s5TargetKeySeedOffers(0),
      s5TargetKeySeedDecreases(0), s5TargetKeyVertexChecks(0),
      s5TargetKeyVertexExact(0), s5TargetKeyVertexTooLow(0),
      s5TargetKeyVertexTooHigh(0), s5TargetKeyListChecks(0),
      s5TargetKeyListExact(0), s5TargetKeyListTooLow(0),
      s5TargetKeyListTooHigh(0), s5TargetKeyInvalidSource(0),
      s5TargetKeyRepairs(0), s5TargetKeyRepairIncidences(0),
      s5TargetKeyMissingRepairs(0),
      s5TargetKeyOwnerTransfers(0), s5TargetKeyDissolveResets(0),
      s5TargetKeyRetiredTrees(0), s5TargetKeyRetiredVertices(0),
      s5TargetKeyTransitionReveals(0),
      s5TargetKeyTransitionRevealIncidences(0),
      s5TargetKeyInvalidRange(0), s5TargetKeyInvalidOwnerChange(0),
      s5TargetKeyInvalidInternal(0), s5TargetKeyInvalidTargetInterval(0),
      s5TargetKeyInvalidMatchingRole(0), s5TargetKeyInvalidSourceRole(0),
      s5TargetMinEventsScheduled(0), s5TargetMinEventsPopped(0),
      s5TargetMinEventsExact(0), s5TargetMinEventsStale(0),
      s5TargetMinEventsDriven(0), s5TargetFrontierCoalesced(0),
      s5TargetFrontierSuperseded(0), s5TargetFrontierInvalidated(0),
      freshRowAuditRows(0), freshRowAuditReusable(0),
      freshRowAuditReusableIncidences(0), freshRowAuditScaleMisses(0),
      freshRowAuditContextMisses(0), freshRowAuditSourceDualMisses(0),
      freshRowAuditSemanticMisses(0), perfectionParkedVertices(0),
      perfectionNonTauCorrections(0), perfectionTightnessFailures(0),
      liquidationistStep6Regions(0),
      liquidationistStep6Searches(0), liquidationistStep6Augmentations(0),
      liquidationistStep6DualUnits(0) {}
    unsigned long long edgeScans, queueInsertions, queuePops;
    unsigned long long bucketQueueInsertions, heapQueueInsertions;
    unsigned long long bucketQueuePops, heapQueuePops;
    unsigned long long splitFindminWork;
    unsigned long long augmentations, dualUnits, blossomFormations;
    unsigned long long blossomDissolutions, boundaryChanges, absorptionHalts;
    unsigned long long shellAugmentationHalts;
    unsigned long long shellOutermostDissolutionHalts;
    unsigned long long dismantleStage1AllDissolvedExits;
    unsigned long long dismantleStage1SingletonExits;
    unsigned long long dismantleStage2Calls;
    // shellMcmNanoseconds is retained as a compatibility alias for the old
    // inclusive counter. New reporting must use shellSearchNanoseconds for
    // the whole Figure-7 call and shellMcmPhaseNanoseconds for exact MCM work.
    unsigned long long shellMcmNanoseconds, shellSearchNanoseconds;
    unsigned long long searchInitializations;
    unsigned long long shellSearchCalls, shellBoundaryConfigurations;
    unsigned long long shellMaskVertices, shellMaskNanoseconds;
    unsigned long long shellMcmPhaseNanoseconds, shellMcmCalls, shellMcmPasses;
    unsigned long long shellMcmAugmentingPhases;
    unsigned long long shellMcmCertificationPhases;
    unsigned long long shellMcmPhaseAugmentations;
    unsigned long long shellMcmMaxPhaseAugmentations;
    unsigned long long shellMcmInitialScanNanoseconds;
    unsigned long long shellMcmGreedyNanoseconds;
    unsigned long long shellMcmAugmentSearchNanoseconds;
    unsigned long long shellMcmCertificationNanoseconds;
    unsigned long long shellMcmRootFilterNanoseconds;
    unsigned long long shellMcmAuditNanoseconds;
    unsigned long long shellDualAdjustmentNanoseconds, shellDualVertexVisits;
    unsigned long long zeroRootCleanupNanoseconds, zeroRootCleanupCalls;
    unsigned long long zeroRootCleanupNodeVisits;
    unsigned long long dismantleNanoseconds, dismantleChainNanoseconds;
    unsigned long long dismantleFreeCountNanoseconds;
    unsigned long long dismantleShellNanoseconds, dismantleProgressNanoseconds;
    unsigned long long dismantleCalls, dismantleIterations;
    unsigned long long eventCriterion1Searches, eventCriterion2Searches;
    unsigned long long eventCriterion3Searches;
    unsigned long long eventSearchNanoseconds, eventSchedulingNanoseconds;
    unsigned long long eventForestInitializationNanoseconds;
    unsigned long long eventFreeDiscoveryNanoseconds;
    unsigned long long eventLabelResetNanoseconds;
    unsigned long long eventRootInitializationNanoseconds;
    unsigned long long eventSfmInitializationNanoseconds;
    unsigned long long eventTargetKeyInitializationNanoseconds;
    unsigned long long eventLazyInitializationNanoseconds;
    unsigned long long s5TargetInitFrontierNanoseconds;
    unsigned long long s5TargetInitTargetsNanoseconds;
    unsigned long long s5TargetInitOwnersNanoseconds;
    unsigned long long s5TargetInitTokensNanoseconds;
    unsigned long long s5TargetInitResizeNanoseconds;
    unsigned long long s5TargetInitFreshNanoseconds;
    unsigned long long s5TargetInitRootsNanoseconds;
    unsigned long long s5TargetInitTargetVisits;
    unsigned long long s5TargetInitOwnerVisits;
    unsigned long long s5TargetInitTokenVisits;
    unsigned long long s5TargetInitRootVisits;
    unsigned long long eventInitializationNodeVisits;
    unsigned long long eventInitializationNodesReset;
    unsigned long long eventInitializationDormantDeferrals;
    unsigned long long eventDissolutionResetCandidates;
    unsigned long long eventScaleResetCandidates;
    unsigned long long eventInitializationFreeCandidates;
    unsigned long long eventInitializationRoots;
    unsigned long long eventFreeCandidateFullScans;
    unsigned long long eventFreeCandidateCacheHits;
    unsigned long long eventFreeCandidateVerticesExamined;
    unsigned long long eventCriterion2AdjustmentOnlyReturns;
    unsigned long long eventCriterion2ZeroDualAugmentReturns;
    unsigned long long eventCriterion2AdjustmentToZeroAugmentPairs;
    unsigned long long eventBoundaryReuseHits, eventBoundaryReuseMisses;
    unsigned long long eventBoundaryEventsSaved, eventBoundaryMaxEventsSaved;
    unsigned long long dissolutionNanoseconds;
    unsigned long long dissolutionRelinkNanoseconds;
    unsigned long long blossomRotationNanoseconds;
    unsigned long long blossomRotationCalls;
    unsigned long long blossomRotationParentSteps;
    unsigned long long blossomRotationFrames;
    unsigned long long blossomRotationMaxDepth;
    unsigned long long blossomRotationNoopCalls;
    unsigned long long eventUsefulPops, eventStaleOwnerPops;
    unsigned long long eventStaleRolePops, eventInternalPops;
    unsigned long long eventIneligiblePops, eventControlPops;
    unsigned long long eventEdgeOffers, eventEdgeFreshOffers;
    unsigned long long eventEdgeTransitionOffers, eventEdgeShellRejects;
    unsigned long long eventEdgeInternalRejects, eventEdgeOddRejects;
    unsigned long long eventEdgeFreshDedupRejects, eventEdgeNoEventRejects;
    unsigned long long eventEdgeCandidates, eventEdgeImmediateCandidates;
    unsigned long long eventEdgeBudgetDrops, eventEdgeBoundaryCoalesced;
    unsigned long long eventEdgeSeedMarkersAvoided;
    unsigned long long eventOwnerMinOwners, eventOwnerMinPayloads;
    unsigned long long eventOwnerMinReplacements;
    unsigned long long eventOwnerMinLaterRejected;
    unsigned long long eventOwnerMinMaxPayloadsPerSearch;
    unsigned long long eventOwnerFutureMatchedGrow;
    unsigned long long eventOwnerFutureMatchedBlossom;
    unsigned long long eventOwnerFutureDissolve;
    unsigned long long eventOwnerFutureOtherControl;
    unsigned long long eventOwnerPayloadChecks, eventOwnerPayloadFailures;
    unsigned long long eventOwnerBoundaryPayloads, eventOwnerBoundaryExpired;
    unsigned long long eventOwnerBoundaryStaleOwner;
    unsigned long long eventOwnerBoundaryStaleRole;
    unsigned long long eventOwnerFutureMaxMatched, eventOwnerFutureMaxControl;
    unsigned long long eventOwnerBoundaryMatched;
    unsigned long long eventOwnerBoundaryMatchedExpired;
    unsigned long long eventOwnerBoundaryControl;
    unsigned long long eventOwnerBoundaryControlExpired;
    unsigned long long matchedGrowTransitions, matchedGrowIneligible;
    unsigned long long matchedGrowDelayOne, matchedGrowDelayTwo;
    unsigned long long matchedGrowDelayOther;
    unsigned long long algorithmFallbackAttempts;
    unsigned long long dissolutionRelinkFailures;
    unsigned long long malformedFreeRootAttempts;
    unsigned long long postInitializationFullReseeds;
    unsigned long long dissolutionBoundaryScans;
    unsigned long long dynamicDissolutionReconstructions;
    unsigned long long dynamicDissolutionReconstructionChildren;
    unsigned long long dynamicDissolutionReconstructionVertices;
    unsigned long long dynamicDissolutionReconstructionScans;
    unsigned long long dynamicDissolutionReconstructionNanoseconds;
    unsigned long long dissolutionTimestampTransfers;
    unsigned long long dissolutionRetainedBaseLeaves;
    unsigned long long dissolutionContainmentSteps;
    unsigned long long dualMaterializations;
    unsigned long long lazyDualVertexActivations, lazyDualNodeActivations;
    unsigned long long lazyDualVerticesMaterialized, lazyDualNodesMaterialized;
    unsigned long long eventBatchedAugmentations;
    unsigned long long eventInvariantChecks, eventInvariantFailures;
    unsigned long long sfmBoundaryQueries, sfmBoundaryExact;
    unsigned long long sfmBoundaryStale;
    unsigned long long sfmBoundaryTooLow, sfmBoundaryTooHigh;
    unsigned long long sfmBoundaryInvalidOwner;
    unsigned long long sfmBoundaryKeyDrift, sfmBoundaryInternalKey;
    unsigned long long sfmBoundaryUnexplainedLow;
    unsigned long long sfmBoundaryTopologyLeak, sfmBoundaryMissingSource;
    unsigned long long sfmBoundaryRepairs, sfmBoundaryRepairScans;
    unsigned long long sfmBoundaryRepairFailures;
    unsigned long long sfmInitializations;
    unsigned long long sfmMergeTopologyFailures, sfmSplitTopologyFailures;
    unsigned long long sfmInitialLists, sfmInitialPivots;
    unsigned long long sfmFormationTopologyNoops, sfmActualSplits;
    unsigned long long sfmGrowSetCalls, sfmGrowSetNoops;
    unsigned long long sfmGrowSetAncestorVisits;
    unsigned long long sfmGrowDecreaseCalls, sfmGrowDecreaseAncestorVisits;
    unsigned long long sfmGrowDecreaseImmediateStops;
    unsigned long long sfmGrowResetNodeVisits;
    unsigned long long s5TargetKeyInitializations, s5TargetKeyOffers;
    unsigned long long s5TargetKeyDecreases, s5TargetKeySeedOffers;
    unsigned long long s5TargetKeySeedDecreases, s5TargetKeyVertexChecks;
    unsigned long long s5TargetKeyVertexExact, s5TargetKeyVertexTooLow;
    unsigned long long s5TargetKeyVertexTooHigh, s5TargetKeyListChecks;
    unsigned long long s5TargetKeyListExact, s5TargetKeyListTooLow;
    unsigned long long s5TargetKeyListTooHigh, s5TargetKeyInvalidSource;
    unsigned long long s5TargetKeyRepairs, s5TargetKeyRepairIncidences;
    unsigned long long s5TargetKeyMissingRepairs;
    unsigned long long s5TargetKeyOwnerTransfers, s5TargetKeyDissolveResets;
    unsigned long long s5TargetKeyRetiredTrees, s5TargetKeyRetiredVertices;
    unsigned long long s5TargetKeyTransitionReveals;
    unsigned long long s5TargetKeyTransitionRevealIncidences;
    unsigned long long s5TargetKeyInvalidRange;
    unsigned long long s5TargetKeyInvalidOwnerChange;
    unsigned long long s5TargetKeyInvalidInternal;
    unsigned long long s5TargetKeyInvalidTargetInterval;
    unsigned long long s5TargetKeyInvalidMatchingRole;
    unsigned long long s5TargetKeyInvalidSourceRole;
    unsigned long long s5TargetMinEventsScheduled, s5TargetMinEventsPopped;
    unsigned long long s5TargetMinEventsExact, s5TargetMinEventsStale;
    unsigned long long s5TargetMinEventsDriven;
    unsigned long long s5TargetFrontierCoalesced;
    unsigned long long s5TargetFrontierSuperseded;
    unsigned long long s5TargetFrontierInvalidated;
    unsigned long long freshRowAuditRows, freshRowAuditReusable;
    unsigned long long freshRowAuditReusableIncidences;
    unsigned long long freshRowAuditScaleMisses, freshRowAuditContextMisses;
    unsigned long long freshRowAuditSourceDualMisses;
    unsigned long long freshRowAuditSemanticMisses;
    unsigned long long perfectionParkedVertices;
    unsigned long long perfectionNonTauCorrections;
    unsigned long long perfectionTightnessFailures;
    unsigned long long liquidationistStep6Regions;
    unsigned long long liquidationistStep6Searches;
    unsigned long long liquidationistStep6Augmentations;
    unsigned long long liquidationistStep6DualUnits;
  };

  ScalingMatching() : mN(0), mM(0), mMaxW(0), mNumScales(0), mTau(1),
                      mUseScaleW(false), mUseFastEngine(false), mCarryDuals(false), mCarryReduced(false),
                      mFaithfulReweight(false), mFaithfulLoop(false), mBoundedSearch(false),
                      mPersistSearch(true),   // persistent-forest step-7: now SAFE (default ON), 3.7-8x faster
                      mEventSearch(false),    // opt-in: event-PQ step-7 search (no O(m) rescan)
                      mCarryM(false),         // EXPERIMENT: carry real-real M across scales (cliff probe)
                      mContinueAug(false),    // CLIFF FIX: continue across dual levels (cliff probe)
                      mIncrRelabel(false),    // CLIFF FIX: incremental entry relabel (cliff probe)
                      mNoRescanGrow(false), mNoRescanDual(false), mGrowScanClock(0),
                      mLastSearchDualUnits(0), mSearchRootClock(0),
                      mHybrid(false), mHybridCanonical(false), mHybridReference(false),
                      mDropPerfectExperimental(false), mSkipFinalizationExperimental(false),
                      mHybNoTranslate(false), mIncrShrink(false), mIncrDissolve(false), mMultiLevel(false), mBridgeFvr(false),
                      mFvrUnitReference(false), mDynamicShell(false),
                      mShellActive(false), mShellCstar(scNull()), mShellDstar(scNull()),
                      mShellLocalResetReady(false),
#ifdef EVENT_SFM_BOUNDARY_FRONTIER
                      mEventSfmFeedFreeBoundary(false),
#endif
#ifdef EVENT_S5_TARGET_KEY_SHADOW
                      mS5TargetKeyRevealClock(0), mS5TargetKeyClockBase(0),
                      mS5TargetKeyTouchClock(1), mS5TargetFrontierClock(1),
                      mS5TargetKeyTopologyScale(scNull()),
                      mS5TargetKeyCriterion(CRIT2),
                      mS5TargetKeyTopologyReady(false), mS5TargetKeyActive(false),
                      mS5TargetKeySeeding(false),
#endif
                      mEventStructuresActive(false), mEventSfmInitialized(false),
                      mEventSfmScale(scNull()), mEventUseBinary(false),
                      mEventLazyClock(false), mEventBatchSameTime(false),
                      mEventCompactInitialization(false),
                      mEventBoundaryReuse(false), mEventBoundaryReuseReady(false),
                      mEventOneAdjustment(false),
                      mEventCriterion2FreeScale(scNull()),
                      mEventCriterion2FreeUniverse(scNull()),
                      mEventDeferringSeedBoundary(false),
                      mEventEarliestDeferredEdge(std::numeric_limits<LL>::max()),
#ifdef EVENT_OWNER_MIN_AUDIT
                      mEventOwnerMinClock(1), mEventOwnerMinCurrentPayloads(0),
#endif
                      mTNow(0), mLazyEpochClock(0), mLcaClock(0),
                      mEventTouchedClock(1),
                      mShellDissolveClock(1),
                      mPreviousCriterion2AdjustmentOnly(false), mLazyDuals(false),
                      mEventDualBudget(-1), mEventEpochDualBase(0),
                      mDismEpochClock(0), mAuditScale(0), mAuditOldBegin(0),
                      mAuditOldEnd(0), mAuditFlattenTraceCount(0),
                      mMatchedWCacheScale(scNull()), mMatchedWCacheMode(false),
                      mOldForestSeparated(false), mLastStatus(STATUS_SUCCESS),
                      mAlgorithmIdentity(ALGORITHM_UNSPECIFIED) {}  // HYBRID §4.2 ShellSearch
  void setUseFastEngine(bool e) { mUseFastEngine = e; }  // Layer 2: per-scale fast engine
  // FAITHFUL step-4: materialize the reweight (w-=y(u)+y(v); y<-0) between scales so
  // per-scale magnitudes stay bounded (weight-independence). Default OFF (gated).
  void setFaithfulReweight(bool e) { mFaithfulReweight = e; }
  // computeBoundedScaling: carry duals (2y'+3, O(sqrt n)/scale) vs fresh (correct, slow).
  void setCarryDuals(bool e) { mCarryDuals = e; }
  // EXPERIMENT: reduced-cost step-4 (feed w-y(u)-y(v) so the search sees small numbers).
  void setCarryReduced(bool e) { mCarryReduced = e; }
  // FAITHFUL multi-scale: per-scale = real-edge search + park-on-dummies (the proven
  // single-scale mechanism), carried duals, Finalization at end. NOT unmatch-all+resolve.
  void setFaithfulLoop(bool e) { mFaithfulLoop = e; }
  // HYBRID: in the per-scale inherited-blossom step, dissolve SMALL blossoms via Gabow's
  // §4.2 O(m n^{3/4}) dissolution instead of liquidation (large ones still liquidate). Off
  // by default; when off the loop is byte-identical Liquidationist.
  void setHybrid(bool e) { mHybrid = e; }
  void setHybNoTranslate(bool e) { mHybNoTranslate = e; }  // DEBUG isolation
  // WEIGHT-INDEPENDENT per-scale (spec steps 6-8): dummies present during the search +
  // tau-bounded SearchOne(Crit2) (maximal augment + ONE dual adjust, tau times) instead of
  // unbounded Crit1-grow-to-perfect. Keeps z = O(tau) so the carried reduced cost w' stays
  // O(sqrt n)-bounded => weight-independent. Requires setFaithfulLoop + setFaithfulReweight.
  void setBoundedSearch(bool e) { mBoundedSearch = e; }
  // EDM optimization (DEFAULT ON): persistent-forest step-7 search (no per-search relabel+regrow),
  // 3.7-8x faster. The former instance-specific hang (n=1600) was unconditional cycle guards
  // missing on the lcaOuter/rotateToBase/flipUpToRoot chain walks (they were SCDBG/CYCDBG-only);
  // now unconditional + bail-safe, so persist terminates + is exact everywhere (PERSIST-SEARCH
  // gate 180/180 vs brute; full green baseline intact). setPersistSearch(false) reverts to the
  // slower from-scratch relabel per search if ever needed.
  void setPersistSearch(bool e) { mPersistSearch = e; }
  // EDM (opt-in): route the bounded step-7 search through the Stage-3 event-PQ BucketSearch
  // (no O(m) per-round rescan). Validated 150/150 as a whole-solve; this gates it into the
  // weight-independent scaling loop's per-scale step.
  void setEventSearch(bool e) { mEventSearch = e; }
  // INCREMENTAL SHRINK (the cold-scale cliff fix): after a blossom shrink, do NOT relabel the
  // whole forest from scratch (shrink already leaves the new blossom EVEN with valid tree
  // pointers -- the forest is consistent). Instead keep growing from the new blossom. Kills the
  // ~940 full-forest relabels/scale at cold scales (3.8s->~0.07s at n=3200, proven blossom-less).
  void setIncrShrink(bool e) { mIncrShrink = e; }
  // Remember exhausted adjacency rows across incremental shrinks at one dual level. Cursors are
  // invalidated conservatively after an augmentation, a dissolution, or a dual adjustment.
  void setNoRescanGrow(bool e) { mNoRescanGrow = e; }
  void setNoRescanDual(bool e) { mNoRescanDual = e; }
  LL lastSearchDualUnits() const { return mLastSearchDualUnits; }
  const OperationCounters& operationCounters() const { return mCounters; }
  ResultStatus lastStatus() const { return mLastStatus; }
  static const char* statusName(ResultStatus s) {
    switch (s) {
      case STATUS_SUCCESS: return "success";
      case STATUS_SUCCESS_EXPERIMENTAL_PARTIAL: return "success-experimental-partial";
      case STATUS_NO_PERFECT_MATCHING: return "no-perfect-matching";
      case STATUS_UNSUPPORTED_OR_OVERFLOWING_INPUT: return "unsupported-or-overflowing-input";
      case STATUS_INVARIANT_VIOLATION: return "invariant-violation";
      case STATUS_UNEXPECTED_NO_PROGRESS: return "unexpected-no-progress";
    }
    return "unknown";
  }
  static const char* algorithmName(AlgorithmIdentity a) {
    switch (a) {
      case ALGORITHM_HYBRID_FIGURE6: return "hybrid-figure6";
      case ALGORITHM_HYBRID_REFERENCE: return "hybrid-reference";
      case ALGORITHM_HYBRID_EVENT_EXPERIMENTAL: return "hybrid-event-experimental";
      case ALGORITHM_HYBRID_DROP_PERFECT_EXPERIMENTAL: return "hybrid-drop-perfect-experimental";
      case ALGORITHM_HYBRID_SKIP_FINALIZATION_EXPERIMENTAL: return "hybrid-skip-finalization-experimental";
      case ALGORITHM_HYBRID_LEGACY_EXPERIMENTAL: return "hybrid-legacy-experimental";
      case ALGORITHM_S5_EXPERIMENTAL: return "s5-experimental";
      case ALGORITHM_LIQUIDATIONIST_FIGURE5: return "liquidationist-figure5";
      default: return "unspecified";
    }
  }
  static const char* leafRepresentationName() {
#ifdef SCALING_LEAF_NOVEC
    return "linked";
#else
    return "legacy-vector";
#endif
  }
  RunConfiguration runConfiguration() const {
    RunConfiguration c;
    c.algorithm = mAlgorithmIdentity;
    c.backend = mAlgorithmIdentity == ALGORITHM_LIQUIDATIONIST_FIGURE5 ?
                  "integrated-scan" :
                (mHybridReference ? "scan-reference" :
                (mHybridCanonical && mEventSearch ?
                   "integrated-event-materialized-shell" :
                (mHybridCanonical ? "integrated-scan" :
                 (mAlgorithmIdentity == ALGORITHM_S5_EXPERIMENTAL ?
                    "detached-s5-experimental" : "legacy"))));
    c.sqrtN = mSqrtN; c.tau = mTau;
    c.fvrStage1Criterion = CRIT2;
    c.fvrStage2Criterion = mAlgorithmIdentity == ALGORITHM_LIQUIDATIONIST_FIGURE5 ?
                            0 : CRIT3;
    c.dismantleStage2Criterion =
      (mAlgorithmIdentity == ALGORITHM_HYBRID_FIGURE6 ||
       mAlgorithmIdentity == ALGORITHM_HYBRID_REFERENCE ||
       mAlgorithmIdentity == ALGORITHM_HYBRID_DROP_PERFECT_EXPERIMENTAL ||
       mAlgorithmIdentity == ALGORITHM_HYBRID_SKIP_FINALIZATION_EXPERIMENTAL) ? CRIT1 : 0;
    c.finalizationCriterion = CRIT3;
    c.incrementalShrink = mIncrShrink;
    c.incrementalDissolve = mIncrDissolve;
    c.persistentSearch = mPersistSearch;
    c.eventSearch = mEventSearch;
    c.lazyDualClock = mEventLazyClock;
    c.eventBatchSameTime = mEventBatchSameTime;
    c.eventCompactInitialization = mEventCompactInitialization;
    c.eventBoundaryReuse = mEventBoundaryReuse;
#ifdef EVENT_CRIT2_FREE_FRONTIER_DISABLE
    c.eventCriterion2FreeFrontier = false;
#else
    c.eventCriterion2FreeFrontier = true;
#endif
#ifdef EVENT_SFM_BOUNDARY_FRONTIER
    c.eventSplitFindminBoundary = mEventSearch;
#ifdef EVENT_SFM_NO_ORACLE
    c.eventSplitFindminScanOracle = false;
#else
    c.eventSplitFindminScanOracle = mEventSearch;
#endif
#else
    c.eventSplitFindminBoundary = false;
    c.eventSplitFindminScanOracle = false;
#endif
#ifdef EVENT_S5_SHARED_SFM_TOPOLOGY
    c.eventSharedSplitFindminTopology = mEventSearch;
#else
    c.eventSharedSplitFindminTopology = false;
#endif
#ifdef EVENT_S5_FULL_KEY_DRIVER_NO_ORACLE
    c.eventSplitFindminGrow = mEventSearch;
    c.eventSplitFindminGrowScanOracle = false;
#elif defined(EVENT_S5_FULL_KEY_DRIVER_AUDIT)
    c.eventSplitFindminGrow = mEventSearch;
    c.eventSplitFindminGrowScanOracle = mEventSearch;
#else
    c.eventSplitFindminGrow = false;
    c.eventSplitFindminGrowScanOracle = false;
#endif
#if defined(EVENT_S5_FULL_KEY_DRIVER_AUDIT) && \
    !defined(EVENT_S5_LIST_FRONTIER_DISABLE)
    c.eventSplitFindminListFrontier = mEventSearch;
#else
    c.eventSplitFindminListFrontier = false;
#endif
#ifdef SCALING_DETAILED_SCHEDULE_TIMING
    c.detailedScheduleTiming = true;
#else
    c.detailedScheduleTiming = false;
#endif
    c.dynamicShell = mDynamicShell;
    c.sensitivityParameters = false;
    return c;
  }
  void setFvrUnitReference(bool e) { mFvrUnitReference = e; }
  void setDynamicShell(bool e) { mDynamicShell = e; } // faithful C*/D* updates + absorption halt
  void setIncrDissolve(bool e) { mIncrDissolve = e; }  // EXPERIMENT: no-relabel live-ODD dissolution
  void setEventLazyClock(bool e) { mEventLazyClock = e; } // differential/performance referee
  // Experimental safe first step toward persistent event phases: after one
  // augmentation, consume additional events at the SAME dual time, but never
  // cross a dual level. Consumed trees are excluded, so augmented paths remain
  // vertex-disjoint within the batch. Default off until differential gates pass.
  void setEventBatchSameTime(bool e) { mEventBatchSameTime = e; }
  void setEventCompactInitialization(bool e) { mEventCompactInitialization = e; }
  void setEventBoundaryReuse(bool e) { mEventBoundaryReuse = e; }
  void setMultiLevel(bool e) { mMultiLevel = e; }      // D step3: multi-dual-level batching per call
  void setBridgeFvr(bool e) { mBridgeFvr = e; }        // D path(a): run the FVR via the S5 engine bridge

  // D PATH (a): run the per-scale FVR via the S5 engine (lazy duals + sqrt(n) per-level batching).
  // Builds a fresh S5 from the loop's CURRENT mAdj (mAdj grows with dummies per scale), loads the
  // carried forest+duals+matching, runs S5's event search (mode-3, Crit2/3) for the tau-bounded FVR,
  // then reads matching+duals back into the loop. Returns #free-real reduction or just runs in place.
  // crit = CRIT2 (stage 1) or CRIT3 (stage 2). stageBudget caps adjustments like the tau-loop.
  void runFvrBridged(int crit) {
    // CALLED at FVR time when mN == mN0 (the loop's rooting switch). The loop's vtx/blossom boundary
    // is blossomBase = 2*mN0 = 2*mN; blossom ids live at >= 2*mN. S5 ties its boundary to ITS mN, so
    // init S5 with 2*mN (= blossomBase) so blossom ids ALIGN, and tell S5 only the first mN (=mN0)
    // vertices are REAL (must be rooted+matched) -- the dummy slots [mN, 2*mN) have no FVR edges.
    Size reals = mN;                                // = mN0 at FVR time
    Size vBoundary = 2 * mN;                         // = blossomBase: S5's vtx/blossom boundary
    ScalingS5Engine s5;
    s5.setUseEvent(true); s5.setLazyEvent(true);
    s5.setNoRelabelShrink(true); s5.setNoRelabelAug(true);
    s5.setCrit(crit == CRIT3 ? 3 : 2);
    s5.init(vBoundary);
    s5.setReals(reals);
    // adjacency over the FULL [0, 2*mN) vertex range (reals + dummy slots). Loop mAdj is symmetric
    // (both dirs stored); add each undirected edge ONCE (u<w). Dummy slots may have 0-weight parked
    // edges (from a prior scale's Perfection) -- include them so parked reals stay matched.
    for (Size u = 0; u < vBoundary; ++u)
      for (Size k = 0; k < mAdj[u].size(); ++k) {
        Size w = mAdj[u][k];
        if (u < w && w < vBoundary) s5.addEdge(u, w, (mUseScaleW ? mW[u][k] : 2 * mWhat[u][k]) / 2);
      }
    // load carried state. (mNumNodes >= 2*mN; blossom ids >= 2*mN now align with S5's boundary.)
#ifdef SCMAG
    { Size nb=0,nbz=0; for (Size b=vBoundary; b<mNumNodes; ++b) if (mBParent[b]==scNull()) { ++nb; if (mZ[b]!=0) ++nbz; }
      Size fr=0; for (Size v=0; v<reals; ++v) if (mMate[v]==scNull()) ++fr;
      std::fprintf(stderr,"[BRIDGE-IN] reals=%lu freeIn=%lu liveBlossoms=%lu (z>0:%lu) mNumNodes=%lu vBoundary=%lu\n",
        (unsigned long)reals,(unsigned long)fr,(unsigned long)nb,(unsigned long)nbz,(unsigned long)mNumNodes,(unsigned long)vBoundary); }
#endif
    s5.loadState(mNumNodes, mMate, mY, mZ, mBParent, mBBase, mBCycle, mBCycEdgeU, mBCycEdgeV, mNodeLeaves, mUf);
#ifdef PS5_LEAF_SPLICE
    s5.loadLeafLists(mLfNextX, mLfPrevX, mLfHeadX, mLfTailX, mLfCntX);
#endif
    // tau-BOUNDED (NOT to-perfect): the FVR must leave reals parked on dummies for Finalization to
    // re-optimize. searchBounded(reals checked) leaves dummy slots alone via setReals.
    // NOTE: one searchEvent multi-augments to a FULL PM (freeAfter=0 measured), so bounding CALLS
    // (searchBounded(mSqrtN)) is INERT -- it still drives to perfect = the decorative mistake. Use
    // the PHASE BUDGET (mSqrtN MV-phases) so still-free reals remain for Finalization. (MBX_BRIDGE_PB)
    if (std::getenv("MBX_BRIDGE_PB")) {
      s5.setPhaseBudget((ScalingS5Engine::Sz)mSqrtN); s5.searchBounded(1); s5.setPhaseBudget(0);
    } else s5.searchBounded(mSqrtN);
#ifdef SCMAG
    { std::vector<Size> mm; s5.readMate(&mm); Size fr=0; for (Size v=0; v<reals; ++v) if (mm[v]==scNull()) ++fr;
      std::fprintf(stderr, "[BRIDGE] reals=%lu freeAfter=%lu\n", (unsigned long)reals, (unsigned long)fr); }
#endif
    // (dissolveZeroRoots cleanup REVERTED -- dissolving z=0 EVEN-roots broke the carry -> hang.
    //  The 17 z=0 blossoms are mostly harmless junk; the real bug is the matching trajectory.)
    s5.readMate(&mMate); std::vector<LL> ny, nz; s5.readDuals(&ny, &nz);
    // Read back the BLOSSOM FOREST S5 built (blossoms formed in the FVR carry to the next scale --
    // their z>0 is carried dual state; without this the loop loses them = the sumZ-divergence bug).
    Size newNumNodes = s5.readNumNodes();
    { std::vector<Size> bp, bb, uf; std::vector<std::vector<Size> > bc, bcu, bcv, nl;
      s5.readForest(&bp, &bb, &bc, &bcu, &bcv, &nl, &uf);
#ifdef PS5_LEAF_SPLICE
      s5.readLeafLists(&mLfNextX, &mLfPrevX, &mLfHeadX, &mLfTailX, &mLfCntX);
#endif
      Size cap = newNumNodes + 1; if (cap < mBParent.size()) cap = mBParent.size();
      // grow loop arrays to cap if S5 created more nodes than the loop currently sizes
      if (newNumNodes > mNumNodes) {
        mBParent.resize(cap, scNull()); mBBase.resize(cap, scNull()); mZ.resize(cap, 0);
        mBCycle.resize(cap); mBCycEdgeU.resize(cap); mBCycEdgeV.resize(cap); mNodeLeaves.resize(cap);
        mNLabel.resize(cap, LBL_FREE_OF_TREE); mNParent.resize(cap, scNull());
        mNParEdgeU.resize(cap, scNull()); mNParEdgeV.resize(cap, scNull()); mNRoot.resize(cap, scNull());
      }
      for (Size b = 0; b < newNumNodes; ++b) { mBParent[b] = bp[b]; mBBase[b] = bb[b];
        mBCycle[b] = bc[b]; mBCycEdgeU[b] = bcu[b]; mBCycEdgeV[b] = bcv[b]; mNodeLeaves[b] = nl[b]; }
      for (Size v = 0; v < mN && v < uf.size(); ++v) mUf[v] = uf[v];
      mNumNodes = newNumNodes;
    }
    for (Size v = 0; v < vBoundary && v < ny.size(); ++v) mY[v] = ny[v];
    for (Size b = 0; b < mNumNodes && b < nz.size(); ++b) mZ[b] = nz[b];
#ifdef SCMAG
    { LL ys=0; for (Size v=0; v<reals; ++v) ys+=mY[v]; std::fprintf(stderr,"[BRIDGE-OUT] sum y(reals)=%lld y[0..3]=%lld,%lld,%lld,%lld\n",(long long)ys,(long long)mY[0],(long long)mY[1],(long long)mY[2],(long long)mY[3]); }
#endif
  }
private:
  // EXPERIMENT: re-thread the alternating tree through an ODD blossom B's cycle in place of a full
  // relabel, so a live-ODD dissolution costs O(|B|) not O(m). Returns true if it re-threaded (caller
  // then dissolves B no-relabel); false is a structured invariant failure.
  bool relinkOddCycle(Size B);
public:

  // Build the integer working graph from a matchbox Graph. Paper-facing
  // scaling entries accept exactly the declared non-negative integral domain;
  // unsupported numeric values are rejected rather than rounded. Sets up adjacency,
  // per-vertex blossom-forest leaves, the dual arrays, wbar, and scale count.
  void init(const Graph& graph);

  // Run the verified single-scale exact blossom solver and fill outputs.
  // Use computeLiquidationist() or computeHybrid() for the scaling algorithms.
  void compute(const Graph& graph,
               std::vector<Size>* mateVec,
               Size* card,
               Val* edgWght);

  // ---- read-only accessors (used by drivers and tests) ----
  Size numVtxs() const { return mN; }
  Size numEdgs() const { return mM; }
  LL   maxWeight() const { return mMaxW; }
  Size numScales() const { return mNumScales; }
  Size tau() const { return mTau; }
  // ith bit (i counted from the most significant downwards, scale 1 = MSB) of
  // the extended weight wbar(e) for directed edge slot at (u, idx).
  // Exposed for scaling-weight sanity checks.
  LL wbarOf(Size u, Size idx) const { return mWbar[u][idx]; }

private:
  // ---- problem size / scaling parameters ----
  Size mN;          // number of (original) vertices
  Size mM;          // number of undirected edges
  LL   mMaxW;       // N = max original integer weight
  Size mNumScales;  // bit length of the inclusive maximum extended weight, >= 1
  Size mTau;        // FVR threshold: Liquidationist sqrt(n); Hybrid sqrt(n) log n (>= small/large split)
  Size mSqrtN;      // ceil(sqrt(n)); Hybrid stage-1 SearchOne count + the small/large blossom split

  // ---- integer working graph (CSR-ish jagged adjacency, mirrors matchbox) ----
  // mAdj[u] = neighbor list; mWhat[u][i] = original integer weight of edge
  // (u, mAdj[u][i]); mWbar[u][i] = extended weight (n/2+1)*what.
  std::vector<std::vector<Size> > mAdj;
  std::vector<std::vector<LL> >   mWhat;
  std::vector<std::vector<LL> >   mWbar;

  // ---- matching: mMate[v] = partner of v, or kNull if free ----
  std::vector<Size> mMate;

  // ---- duals (per scale; integer) ----
  // y[v]: vertex dual (integer). z[B]: blossom dual (even, >= 0).
  std::vector<LL> mY;
  std::vector<LL> mZ;            // indexed by blossom id

  // ---- S5 LAZY-DUAL timestamp scheme (paper S5.2; see docs/Scaling_S5_DataStructures.md)
  // Makes a dual adjustment O(1) (just ++mTNow) instead of an O(n) relabel, so the
  // alternating forest can PERSIST across dual adjustments (no per-iteration rebuild
  // = the fair-vs-GMG/Gabow goal). Current y(v)/z(B) are computed ON DEMAND from
  // these timestamps. Validated against the materialized mY/mZ before being wired in.
  BucketPQ mPQ;                  // S5 event scheduler (grow/formation/dissolution)
  BinaryEventHeap mUnboundedPQ;  // Criterion-1/3 PQSearch without a bounded horizon
  SplitFindmin mInnerSfm;         // inner-blossom min incident event ownership
  std::vector<std::vector<Size> > mSfmSplitPivots;
#ifdef EVENT_SFM_SHADOW
  std::vector<Size> mSfmKeySource;
  std::vector<Size> mSfmKeyOwner;
#endif
#ifdef EVENT_S5_TARGET_KEY_SHADOW
  // Non-driving Section 5 key(u)=min slack*(v,u), v outer. The default keeps
  // a separate AVL oracle; EVENT_S5_SHARED_SFM_TOPOLOGY uses an independent
  // key channel in mInnerSfm while sharing its validated laminar topology.
#ifndef EVENT_S5_SHARED_SFM_TOPOLOGY
  SplitFindmin mS5TargetKeySfm;
#endif
  std::vector<std::vector<Size> > mS5TargetKeySplitPivots;
  std::vector<Size> mS5TargetKeySource;
  std::vector<Size> mS5TargetKeyOwner;
  std::vector<Size> mS5TargetKeyTargetOwner;
  std::vector<Size> mS5TargetKeySourceToken;
  std::vector<Size> mS5TargetKeyEdgeIndex;
  std::vector<char> mS5TargetKeyMatched;
  std::vector<std::vector<Size> > mS5TargetKeyOwnerTargets;
  std::vector<std::vector<Size> > mS5TargetKeyTargetOwnerTargets;
  std::vector<std::vector<Size> > mS5TargetKeyTreeNodes;
  std::vector<std::vector<Size> > mS5TargetKeyTreeVertices;
  std::vector<std::vector<Size> > mS5TargetKeyTreeAllVertices;
  std::vector<std::vector<Size> > mS5TargetKeyTreeOwnedTargets;
  std::vector<Size> mS5TargetKeyTrackedRoot;
  std::vector<Size> mS5TargetKeyTrackedVertexToken;
  std::vector<Size> mS5TargetKeyAllVertexToken;
  std::vector<unsigned int> mS5TargetKeyTrackedRootStamp;
  std::vector<unsigned int> mS5TargetKeyTrackedVertexStamp;
  std::vector<unsigned int> mS5TargetKeyAllVertexStamp;
  std::vector<unsigned int> mS5TargetKeyRevealStamp;
  unsigned int mS5TargetKeyRevealClock;
  LL mS5TargetKeyClockBase;
  std::vector<Size> mS5TargetKeyTouchedTargets;
  std::vector<Size> mS5TargetKeyTouchedOwners;
  std::vector<Size> mS5TargetKeyTouchedTokens;
  std::vector<unsigned int> mS5TargetKeyTargetTouchStamp;
  std::vector<unsigned int> mS5TargetKeyOwnerTouchStamp;
  std::vector<unsigned int> mS5TargetKeyTokenTouchStamp;
  unsigned int mS5TargetKeyTouchClock;
#ifdef EVENT_S5_SINGLETON_TARGET_FAST
  std::vector<char> mS5TargetKeySingleton;
  std::vector<Val> mS5TargetKeySingletonValue;
  std::vector<Size> mS5TargetKeySingletonListId;
#endif
  struct S5TargetFrontierEntry {
    Val fire;
    Size winner, source;
    unsigned int version;
    bool touched;
    S5TargetFrontierEntry()
      : fire(SplitFindmin::INF), winner(cNullItm), source(cNullItm),
        version(0), touched(false) {}
  };
  std::vector<S5TargetFrontierEntry> mS5TargetFrontier;
  std::vector<Size> mS5TargetFrontierTouched;
  unsigned int mS5TargetFrontierClock;
  Size mS5TargetKeyTopologyScale;
  int mS5TargetKeyCriterion;
  bool mS5TargetKeyTopologyReady;
  bool mS5TargetKeyActive;
  bool mS5TargetKeySeeding;
#endif
#ifdef EVENT_SFM_BOUNDARY_FRONTIER
  // True only while a dissolution exposes EVEN siblings whose formerly
  // internal edges must seed the freshly split FREE-child intervals.
  bool mEventSfmFeedFreeBoundary;
#endif
  bool mEventStructuresActive;
  bool mEventSfmInitialized;
  Size mEventSfmScale;
  bool mEventUseBinary;
  bool mEventLazyClock;
  bool mEventBatchSameTime;
  bool mEventCompactInitialization;
  bool mEventBoundaryReuse;
  bool mEventBoundaryReuseReady;
  bool mEventOneAdjustment;
  // Initial SearchOne seeding can discover a decreasing sequence of boundary
  // minima. Delay its one control marker until the complete root pass ends.
  bool mEventDeferringSeedBoundary;
  LL mEventEarliestDeferredEdge;
#ifdef EVENT_OWNER_MIN_AUDIT
  // Shadow-only compact frontier: one minimum future ordinary edge event per
  // current owner, with all ties. Epoch stamps avoid an O(nodes) clear at each
  // SearchOne. This does not drive the queue until its lifecycle is verified.
  std::vector<unsigned int> mEventOwnerMinEpoch;
  std::vector<LL> mEventOwnerMinTime;
  std::vector<unsigned int> mEventOwnerMinTies;
  std::vector<unsigned int> mEventOwnerMinVersion;
  std::vector<BucketPQ::Event> mEventOwnerMinPrimary;
  struct EventOwnerMinTie {
    Size owner;
    unsigned int version;
    BucketPQ::Event event;
  };
  std::vector<Size> mEventOwnerMinTouched;
  std::vector<EventOwnerMinTie> mEventOwnerMinTieLog;
  std::vector<BucketPQ::Event> mEventOwnerMatchedFuture;
  std::vector<BucketPQ::Event> mEventOwnerControlFuture;
  unsigned int mEventOwnerMinClock;
  unsigned long long mEventOwnerMinCurrentPayloads;
#endif
#ifdef EVENT_FRESH_ROW_FRONTIER_AUDIT
  // Oracle-only cache census. A digest describes every value that determines
  // the ordinary fresh-row event offers at the moment the exact scan runs.
  // It never drives scheduling; it measures which rows a future versioned
  // frontier could reuse without retaining an alternating forest.
  struct FreshRowAuditEntry {
    FreshRowAuditEntry() : valid(false), scale(0), criterion(0), owner(scNull()),
                           sourceY(0), time(0), degree(0), digest(0) {}
    bool valid;
    Size scale;
    int criterion;
    Size owner;
    LL sourceY;
    LL time;
    Size degree;
    unsigned long long digest;
  };
  std::vector<FreshRowAuditEntry> mFreshRowAudit;
#endif
  LL mTNow;                      // current time = # dual adjustments so far
  std::vector<LL>   mY0;         // y0(v): initial y of each vertex
  std::vector<LL>   mZ0;         // z0(B): initial z of each blossom (B may be non-root)
  std::vector<LL>   mTRoot;      // t_root(B): time B became a root blossom (scNull=undef)
  std::vector<LL>   mTIn;        // t_in(B): time B became an INNER root blossom
  std::vector<LL>   mTOut;       // t_out(B): time B became (part of) an OUTER root blossom
  std::vector<LL>   mDelta;      // Delta(B): #dual adj B's verts saw as INNER in [0,max(troot,tout)]
  std::vector<LL>   mYAcc;       // closed vertex-role intervals within the current search
  std::vector<LL>   mZAcc;       // closed blossom-role intervals within the current search
  std::vector<LL>   mVStamp, mNStamp;
  std::vector<int>  mVRole, mNRole;
  // Each timestamp search remains a fresh logical epoch. Epoch stamps avoid
  // copying and clearing the whole graph for singleton DismantlePath searches:
  // an untouched entry inherits its already-materialized mY/mZ value.
  std::vector<Size> mLazyVEpoch, mLazyNEpoch;
  std::vector<Size> mLazyActiveVertices, mLazyActiveNodes;
  Size mLazyEpochClock;
  // Alternating-forest LCA marks. Grow geometrically: mNumNodes increases one
  // blossom at a time, so exact-size growth would refill an ever-larger array
  // once per formation.
  std::vector<Size> mLcaSeen;
  Size mLcaClock;
  // Nodes labeled during the current event search.  The next SearchOne resets
  // exactly this set instead of sweeping every allocated/absorbed blossom.
  // Epoch stamps deduplicate transitions (a node can change role repeatedly).
  std::vector<Size> mEventTouchedStamp;
  std::vector<Size> mEventTouchedNodes;
  std::vector<char> mEventDormantFlag;
  std::vector<Size> mEventDormantNodes;
  std::vector<Size> mEventFreeCandidates;
  // During Figure-6 stage-1 FVR, Criterion-2 SearchOne calls only augment:
  // the set of exposed real roots monotonically shrinks until the scale ends.
  // Retain that exact set and compact matched endpoints instead of rediscovering
  // it with an O(n) vertex sweep at every one-adjustment boundary.
  std::vector<Size> mEventCriterion2Free;
  Size mEventCriterion2FreeScale;
  Size mEventCriterion2FreeUniverse;
  std::vector<Size> mEventInitialRoots;
  std::vector<BucketPQ::Event> mEventDeferredEdgeEvents;
  std::vector<BucketPQ::Event> mEventBoundaryReuseEvents;
  Size mEventTouchedClock;
  // Deduplicate live INNER root blossoms when initializing a shell search.
  // Shell initialization must be proportional to the current shell, not to
  // every historical blossom id allocated by earlier scales.
  std::vector<Size> mShellDissolveStamp;
  Size mShellDissolveClock;
  bool mPreviousCriterion2AdjustmentOnly;
  bool mLazyDuals;               // when true, y()/z() use the timestamp scheme
  // Event times are local to the current queue epoch.  These fields convert a
  // local fire time to the request's cumulative bounded-search clock.
  LL mEventDualBudget, mEventEpochDualBase;

  // ---- blossom forest scaffolding ----
  // Blossoms 0..mN-1 are the trivial (singleton) blossoms = vertices.
  // Non-trivial blossoms get ids >= mN. mBlossomParent[b] = enclosing blossom
  // or kNull if b is a root. (Full structure filled in Phase 1+.)
  std::vector<Size> mBlossomParent;

  // current scale's even integer weight w(e), filled per scale in Phase 3.
  // mW[u][idx] mirrors mAdj[u][idx]. When mUseScaleW is true the search uses
  // mW as the edge weight w(e); otherwise it uses 2*mWhat (the Gate-1a/1b path).
  std::vector<std::vector<LL> > mW;
  bool mUseScaleW;
  bool mUseFastEngine;   // Layer 2: use ScalingFastMWPM per scale instead of edmondsSearch
  bool mCarryDuals;      // computeBoundedScaling: carry duals 2y'+3 (vs fresh per scale)
  bool mCarryReduced;    // EXPERIMENT: reduced-cost step-4 (feed w-y(u)-y(v), small numbers)
  bool mFaithfulReweight; // FAITHFUL step-4: materialize reweight between scales (bounded magnitude)
  bool mFaithfulLoop;    // FAITHFUL multi-scale per-scale step (real-edge search + park, not resolve)
  bool mHybrid;          // HYBRID: small inherited blossoms dissolved via Gabow §4.2 (vs liquidate)
  bool mHybridCanonical; // Figure 6 path: separated old/current forests and dynamic shells
  bool mHybridReference; // same semantics using the deliberately simple scan backend
  bool mDropPerfectExperimental; // accept a terminal partial matching after no augmenting path
  bool mSkipFinalizationExperimental; // publish the post-scale matching without final searches
  bool mHybNoTranslate;  // DEBUG: disable the C*/D* unit translation (to isolate its effect)
  bool mIncrShrink;      // incremental shrink: skip relabel-on-shrink, keep growing (cliff fix)
  bool mIncrDissolve;    // EXPERIMENT (default OFF): no-relabel live-ODD dissolution via cycle re-thread
  bool mMultiLevel;      // D step3 (default OFF): edmondsSearchEx processes MULTIPLE dual levels/call
  bool mBridgeFvr;       // D path(a) (default OFF): run the FVR via the S5 engine bridge
  bool mFvrUnitReference; // validation: disable event-free batching, apply one dual unit/call
  bool mDynamicShell;     // experimental faithful ShellSearch boundary/absorption policy
                         // (batch like S5: multi-augment per level so trees get consumed before
                         // dissolutions go live -> fewer live-odd relabels = the n>=51200 cliff fix)
  // S5-style per-search disjointness marks for the no-relabel AUGMENT: a vertex consumed by an
  // augmenting path this search is marked, and any later edge touching a marked vertex is skipped
  // (S5 :927). Robust replacement for the fragile treeLive consumed-tree retirement. Sized per
  // search when mIncrShrink (the no-relabel-augment path).
  std::vector<char> mUsedV;
  bool mBoundedSearch;   // WEIGHT-INDEPENDENT per-scale: dummies-present + tau-bounded SearchOne(Crit2)
  bool mPersistSearch;   // EDM opt-in: persistent-forest step-7 search (fast, has a latent hang)
  bool mEventSearch;     // EDM opt-in: event-PQ step-7 search (no O(m) rescan)
  bool mCarryM;          // EXPERIMENT: carry real-real M across scales (don't empty) -- cliff probe
  bool mContinueAug;     // CLIFF FIX: continue one edmondsSearchEx call across dual levels (retire dead trees) instead of returning after each augment -- removes the per-call O(n) entry relabel
  bool mIncrRelabel;     // CLIFF FIX (safe): entry relabel resets only PREVIOUSLY-LABELED nodes (mLabeledList) not all O(mNumNodes) -- identical result, O(labeled) cost
  bool mNoRescanGrow;    // per-vertex adjacency cursor survives shrink-only worklist re-entry
  bool mNoRescanDual;    // retain exhausted rows across pure dual changes; wake min-event rows only
  Size mGrowScanClock;   // generation for mGrowScanEpoch; zero means "not scanned"
  std::vector<Size> mGrowScanEpoch;
  std::vector<Size> mGrowScanPos;
  std::vector<Size> mDualWakeU; // rows containing an edge that attains the next dual event
  std::vector<Size> mDualWakeK; // matching adjacency index for mDualWakeU
  LL mLastSearchDualUnits;      // actual y-unit distance consumed by the most recent search call
  // Epoch-marked explicit Edmonds root set.  A zero local epoch means the traditional
  // policy (all unmatched vertices); a nonzero epoch means root only marked vertices.
  // This avoids an O(n) mask clear for singleton/confined searches.
  std::vector<Size> mSearchRootEpoch;
  Size mSearchRootClock;
  std::vector<Size> mLabeledList;  // nodes given a non-FREE label since last reset (for incremental entry relabel)
  std::vector<std::vector<LL> > mWoff;  // accumulated reweight offset per edge slot (to recover true weight)
  std::vector<Size> mGrowWq;            // EDM: reused grow-worklist buffer (avoid per-pass vector alloc)
  // Reused scratch for shrink's cycle collection (was 6 local std::vectors malloc'd EVERY shrink ->
  // ~17K allocs/scale at the cold-scale spike, the ~16% stl_vector cost). clear() keeps capacity.
  std::vector<Size> mShrEside, mShrEu, mShrEv, mShrWside, mShrWu, mShrWv;

  // Helpers are deliberately distinct. ceilLog2 is used by parameter choices
  // such as tau; bitLength is the number of binary weight positions and must
  // retain the leading 1 when the inclusive maximum is an exact power of two.
  static Size ceilLog2(LL x);
  static Size bitLength(LL x);
  static bool checkedAddLL(LL a, LL b, LL* out);
  static bool checkedMulLL(LL a, LL b, LL* out);

  // ========================================================================
  // PHASE 1a — MAX-weight PERFECT matching, primal-dual, NO blossoms
  // (correctness-first; paper-faithful MAXIMIZATION direction).
  //
  // Goal: MAXIMUM-weight PERFECT matching (every vertex matched), the same
  // problem the Duan-Pettie-Su scaling paper solves. Implemented in the
  // maximization direction the paper uses (no min-reduction).
  //
  // Even-weight discipline: work with w2(e) = 2*what(e) so all slacks are even.
  // Duals y[v] (integer). Domination (max): y(u)+y(v) >= w2(e) for all edges;
  //   slack(e) = y(u)+y(v) - w2(e) >= 0; tight = 0 (Criterion 1).
  // Search: grow alternating trees from ALL free vertices using tight edges;
  //   NO vertex-drop (perfect matching => everyone must match); augment when two
  //   trees meet via a tight edge; when stuck, dual-adjust (EVEN y -= delta,
  //   ODD y += delta) to expose a new tight edge. Repeat until no free vertices.
  // If an odd cycle would be needed -> mBlossomNeeded=true (1a bail; Phase 1b
  // adds shrink/expand). If free vertices remain but no dual step can help and
  // no blossom is needed -> NO perfect matching exists (mNoPerfect=true).
  // ========================================================================

  enum Label { LBL_FREE_OF_TREE = 0, LBL_EVEN, LBL_ODD };
  std::vector<int>  mLabel;     // LBL_* per vertex
  std::vector<Size> mTreeRoot;  // which free-root's tree v belongs to (or scNull)
  std::vector<Size> mParent;    // tree parent vertex (the vertex we grew from)
  bool mBlossomNeeded;          // set if a blossom would be required (1a bail)
  bool mNoPerfect;              // set if no perfect matching exists

  // slack under DOUBLED weights (maximization domination): >= 0, tight == 0.
  LL slack1a(Size u, Size idx) const {
    return mY[u] + mY[mAdj[u][idx]] - 2 * mWhat[u][idx];
  }

  // Run the no-blossom max-weight PERFECT search. Fills mMate. Returns true on
  // success; sets mBlossomNeeded / mNoPerfect on the two failure modes.
  bool runPerfectNoBlossom();

  // flip one side of the augmenting path from EVEN vertex `endpt` up to its root.
  void augmentPath1a(Size endpt);

  // ========================================================================
  // PAPER'S integer Edmonds search (DPS arXiv:1411.1919, Fig 2/3). Built ONCE;
  // the scaling loop (Phase 3) reuses it. See docs/Scaling_Paper_Exact_Spec.md
  // (transcribed verbatim from the paper) -- THAT file is the contract.
  //
  // Key paper facts (not my reconstruction):
  //  * yz(e) = y(u)+y(v) + sum_{B contains e} z(B).  y integer, z even >=0.
  //  * Eligibility by CRITERION (1: tight yz=w; 2: E_B or (e notM & yz=w-2) or
  //    (e in M & yz=w); 3: yz=w or yz=w-2). G_elig = eligible edges; contract Ω.
  //  * Dual adjustment is UNIT STEPS, looped: y-=1 on V_out (EVEN-reachable),
  //    y+=1 on V_in (ODD-reachable), z+=2 on root blossom in V_out, z-=2 on root
  //    blossom in V_in. NOT a single delta jump. NO d1 vertex-drop.
  //  * EdmondsSearch loops {Augment, Shrink, DualAdjust, Dissolve}, halting after
  //    the first Augment that finds >=1 augmenting path. SearchOne = one maximal
  //    augment-set + one dual adjustment (Criterion 2).
  //  * A single search does NOT reach the optimum -- it preserves Property 1/2
  //    and grows |M|. Optimality emerges from the SCALING loop. Validate the
  //    search by INVARIANT preservation, not single-search oracle weight.
  //
  // Blossom forest (ids 0..n-1 = vertices; >= n = blossoms):
  //   mBParent[i], mBBase[i], mBCycle[B], mBCycEdgeU/V[B][k], mZ[B] (even >=0).
  // Weights here are DOUBLED (w2=2*what) so all slacks are even. The criterion
  // compares yz(e) against w(e)=w2 and w2-2.
  // ========================================================================

  enum Crit { CRIT1 = 1, CRIT2 = 2, CRIT3 = 3 };

  std::vector<Size> mBParent;
  std::vector<Size> mBBase;
  // UNION-FIND over original vertices -> the BASE vertex of their current outer node (the S5
  // mechanism for O(alpha) outer-node resolution that stays correct under no-relabel shrink).
  // mUf[v] = uf-parent; ufFind(v) -> rep = base vertex of v's outer blossom. Maintained on
  // shrink (union B's leaves into base rep) + dissolve (re-split). Foundation for the cliff fix.
  mutable std::vector<Size> mUf;   // mutable: ufFind does path-halving, callable from const slackZ
  // Per-blossom: has this blossom EVER had z>0 (been EVEN-adjusted since formation)? A z=0 ODD
  // blossom that was NEVER positive was just born at z=0 and flipped ODD without earning
  // dissolution -- dissolving it undoes the shrink and re-forms next regrow (the cold-scale spin).
  // Only dissolve z=0 ODD blossoms that WERE positive (z adjusted down to 0). Sized like mZ.
  std::vector<char> mBEverPos;
  // base-vertex -> current outer-node id (the inverse of mBBase for outer nodes). Lets
  // outerOf(w) = mBaseToNode[ufFind(w)] in O(alpha) instead of the O(depth) mBParent walk
  // (the 1.58B-step cold-scale cost). Maintained wherever mBBase/mBParent of an outer node change.
  std::vector<Size> mBaseToNode;
  std::vector<std::vector<Size> > mBCycle;
  std::vector<std::vector<Size> > mBCycEdgeU;
  std::vector<std::vector<Size> > mBCycEdgeV;
  // mNodeLeaves[B] = all original vertices (leaves) under node B. For a vertex
  // v it is {v}; for a blossom B it is the union of its children's leaves.
  // Lets the grow scan iterate a node's members in O(|B|) instead of scanning
  // all n vertices (kills the O(n^2) grow). Built on shrink; children keep theirs.
  std::vector<std::vector<Size> > mNodeLeaves;
  std::vector<char> mBActive; // explicit liveness; cycle storage survives generic dissolve
#ifdef PS5_LEAF_SPLICE
  // Stage-3 mailbox: the S5 leaf-list arrays carried between scales (O(n); replaces shipping
  // the leaf vectors once PS5_LEAF_NOVEC lands). Empty until the first readLeafLists.
  std::vector<Size> mLfNextX, mLfPrevX, mLfHeadX, mLfTailX, mLfCntX;
  // During Figure 6 dissolution, inherited Omega' and mutable current Omega
  // coexist. They cannot share one mutable intrusive order: splicing a current
  // blossom may split an inherited segment. Snapshot Omega' once per scale as
  // one DFS leaf order plus half-open intervals, still O(n + #nodes).
  std::vector<Size> mOldLfOrderX, mOldLfBeginX, mOldLfEndX;
  void scalingLfInit(Size vertexCount, Size nodeCapacity);
  void scalingLfSnapshotOld(Size oldBegin, Size oldEnd);
  void scalingLfEnsureFromVectors();
  void scalingLfOnShrink(Size B);
  void scalingLfOnDissolve(Size B);
  void scalingLfCheckNode(Size B, const char* where) const;
  Size scalingLfCount(Size B) const;
  Size directChildContaining(Size B, Size v,
                             unsigned long long* steps = NULL) const;
  void scalingLfCollect(Size B, std::vector<Size>* out) const;
#endif
  Size nodeLeafCount(Size B) const;
  void collectNodeLeaves(Size B, std::vector<Size>* out) const;
  Size mNumNodes;                       // total allocated nodes (verts+blossoms)

  std::vector<int>  mNLabel;            // LBL_* per outer node (EVEN=V_out, ODD=V_in)
  std::vector<Size> mNParent;           // parent OUTER node in the alternating forest
  std::vector<Size> mNParEdgeU;         // edge (orig u in this node,
  std::vector<Size> mNParEdgeV;         //       orig v in parent node) that grew it
  std::vector<Size> mNRoot;             // tree root (a free outer node)

  // HYBRID §4.2 ShellSearch confinement (OFF by default; existing edmondsSearchEx callers
  // see mShellActive=false => byte-identical behavior). When ON, the search roots ONLY at
  // free vertices in the current atomic shell and treats STRADDLING edges (leaving the
  // shell) as ineligible. mInShell[v]=1 iff vertex v is in the current G(C*,D*).
  bool mShellActive;
  std::vector<char> mInShell;           // per-vertex shell membership (set on shell entry)
  Size mShellCstar, mShellDstar;        // current shell's C* (outer undissolved) / D* (inner; scNull if none)
  // mInDstar[v]=1 iff v in D* (gets y+=2*delta); a shell vertex not in D* is in C*\D* (y+=1*delta).
  std::vector<char> mInDstar;
  std::vector<Size> mShellMarked;        // vertices whose shell masks must be cleared
  std::vector<LL> mShellNextY;           // checked translated y values before atomic commit
  std::vector<Size> mShellNextZNodes;    // labeled current roots updated by shell dual step
  std::vector<LL> mShellNextZ;           // checked z values paired with mShellNextZNodes
  bool mShellLocalResetReady;            // gated: first shell reset full, later resets touched outers
  std::vector<Size> mDismSeenEpoch;      // same-iteration absorption detection
  std::vector<Size> mDismSeenOwner;
  Size mDismEpochClock;
  // HYBPROF scale-entry boundary: inherited Omega_0 ids are exactly
  // [mAuditOldBegin,mAuditOldEnd); later ids belong to current-scale Omega.
  Size mAuditScale;
  Size mAuditOldBegin;
  Size mAuditOldEnd;
  Size mAuditFlattenTraceCount;
  // A matched edge is queried repeatedly while its endpoints move through the
  // event forest.  Its scale weight is immutable, so cache the heaviest
  // parallel-edge value by current mate pair for the lifetime of one scale.
  std::vector<Size> mMatchedWCacheMate;
  std::vector<LL> mMatchedWCacheValue;
  Size mMatchedWCacheScale;
  bool mMatchedWCacheMode;
  bool mOldForestSeparated;
  ResultStatus mLastStatus;
  AlgorithmIdentity mAlgorithmIdentity;
  OperationCounters mCounters;
  std::vector<Size> mOldParent;
  std::vector<std::vector<Size> > mOldContaining;

  Size outerOf(Size v) const;           // outermost node containing vertex v
  // UNION-FIND helpers (S5 mechanism). ufFind(v) returns the BASE vertex of v's current outer
  // node in O(alpha). ufUnion merges. ufOuter(v) -> the outer NODE id of v (= the node whose
  // base is ufFind(v)); equals outerOf(v) once UF is maintained, in O(alpha) not O(depth).
  Size ufFind(Size v) const;
  Size outerOfUF(Size v) const;         // O(alpha) outer-node via UF + base->node map
  void ufUnionInto(Size member, Size baseRep);   // mUf[find(member)] = baseRep
  Size ufOuter(Size v);                 // outer node of v via UF (O(alpha))
  // CONTRACTED-forest tree parent (S5 treeParent :168): outer node of mNParent[X], resolving
  // through any blossom that absorbed the stored parent; scNull at a root or if it resolves to X
  // (absorbed into the same blossom). Used by the no-relabel augment so the tree stays consistent
  // under incremental shrink without a relabel.
  Size treeParentR(Size X) {
    Size p = mNParent[X];
    if (p == scNull()) return scNull();
    Size po = outerOf(p);   // KEEP O(depth): outerOfUF here broke correctness (mid-search map staleness)
    return (po == X) ? scNull() : po;
  }
#ifdef UFCHECK
  void ufVerify(const char* where);     // assert ufOuter(v)==outerOf(v) for all v (debug)
#endif
  LL   slackZ(Size u, Size v) const;    // yz(e) - w2(e)  (>= -2 under Property 2)
  LL   oldCommonZ(Size u, Size v) const;
  LL   slackZW(Size u, Size v, LL w2) const; // slackZ core given the edge weight (no w2max re-scan)
  LL   w2max(Size u, Size v) const;     // doubled weight of heaviest parallel (u,v)
  LL   matchedEdgeW2(Size u, Size v);   // scale-local cache for the current mate edge
  bool eligible(Size u, Size v, int crit) const;   // criterion test on edge (u,v)
  bool eligibleAt(Size u, Size idx, int crit) const; // eligible reading weight from mAdj[u][idx] (O(1))
  bool eligibleAtCross(Size u, Size idx, int crit) const; // eligibleAt for a PROVEN cross-blossom edge (no ufFind/z-sum)
  Size outerOfUFRep(Size rep, Size v) const;        // outerOfUF with ufFind(v) precomputed (rep)
  bool treeLive(Size N) const;          // 2b: tree of N still has an unmatched root base
  void nextGrowScanEpoch();             // invalidate per-row grow cursors after a semantic change

  // ---- S5 lazy-dual accessors (paper S5.2 formulas) ----
  // listOf(v): the root blossom B used in the y(v) formula = last root blossom
  //   containing v just before v became outer, or the current root if v non-outer.
  //   We use the search-label of that node to pick the formula branch.
  LL yOf(Size v) const;                 // current y(v) from timestamps
  LL yOfOuter(Size v, Size outer) const;// current y(v), outer node already resolved
  LL zOf(Size B) const;                 // current z(B) from timestamps
  // maintenance hooks (called by the persistent search at events):
  void tsInit();                        // init timestamps from materialized mY/mZ at t=0
  void tsEnsureVertex(Size v);          // activate one vertex in this epoch
  void tsEnsureNode(Size B);            // activate a blossom and its leaves
  void tsOnDissolveInner(Size Bprime);  // inner blossom Bprime dissolves
  void tsOnFormOuter(Size Bprime);      // outer blossom Bprime created (cycle = its subs)
  void materializeLazyDuals();          // close the current timestamp epoch into mY/mZ
  void lazySetVertexRole(Size v, int role);
  void lazySetNodeRole(Size B, int role);

  // EdmondsSearch(F): F = the free OUTER nodes to search from. Loops the four
  // steps until >=1 augmenting path is found and augmented (returns true), or no
  // progress possible (returns false). maxAdjust caps applied adjustment events
  // (<0 = unbounded); maxDualUnits independently caps their total delta. A literal
  // SearchOne uses both caps equal to 1. pqDelta=true batches an event-free interval
  // into one adjustment whose full delta is charged to maxDualUnits.
  bool edmondsSearch(int crit, bool maximalAugment, LL maxAdjust, bool pqDelta = false,
                     LL maxDualUnits = -1,
                     const std::vector<Size>* explicitRoots = 0,
                     bool exhaustDualBudget = false);
  // edmondsSearch + persistAcrossAdjust: skip the repair+relabel after a pure
  // dual-adjust round (forest persists; no structural change). Step 2a of the
  // persistent-forest engine. edmondsSearch() = edmondsSearchEx(...,false).
  bool edmondsSearchEx(int crit, bool maximalAugment, LL maxAdjust, bool pqDelta,
                       bool persistAcrossAdjust, bool persistAcrossAugment = false,
                       LL maxDualUnits = -1,
                       const std::vector<Size>* explicitRoots = 0,
                       bool exhaustDualBudget = false);
  SearchResult executeSearch(const SearchRequest& request, SearchBackend backend);
  bool validateSearchEntry(const SearchRequest& request, bool full) const;
  bool validateEventInvariantState(int criterion);

  // building blocks
  void shrink(Size uOuter, Size vOuter, Size uOrig, Size vOrig);
  void augmentFrom(Size uOuter, Size vOuter, Size uOrig, Size vOrig);
  void flipUpToRoot(Size startOuter);
  Size lcaOuter(Size a, Size b);
  void rotateToBase(Size B, Size newBaseVtx);
  void dissolve(Size B);
  void flattenNode(Size B);             // fully dissolve B + descendants to repair

  // Strip dummy (id>=mN0) edges from every real vertex's adjacency, IN PLACE (compact +
  // resize, no per-vertex temp-vector allocation), AND clear the dummy vertices' adjacency.
  // Replaces the old `na/nw/nh/nwb` push_back-into-fresh-vectors pattern that was ~40% of the
  // n>=3200 runtime in malloc/memcpy/realloc (callgrind). Keeps mAdj/mW/mWhat/mWbar in sync.
  void stripDummyAdj(Size mN0);
  bool runPerfectBlossom();             // false if no perfect matching
  bool runPerfectBlossomPQ();           // same, but delta-jump (PQSearch) dual-adjust
  bool runPerfectBlossomPersist();      // delta-jump + persist forest across dual-adjusts
  bool runPerfectBlossomEvent();        // STAGE 3: event-PQ search, no O(m) edge re-scan
  // one event-driven search: seed events for EVEN nodes' edges, pop min-time event,
  // advance duals, process grow/formation/dissolution. Returns true if augmented.
  bool edmondsSearchEvent(int crit, const std::vector<Size>* explicitRoots = 0,
                          LL maxDualUnits = -1, bool exhaustDualBudget = false,
                          bool batchSameTimeAugment = false,
                          LL maxAdjustments = -1);
  bool eventShellMcmAtCurrentDual();
  void scheduleNodeEdges(Size E, int crit, LL tNow,
                         bool dedupFreshEven = false); // schedule events for EVEN node E
  bool scheduleMatchedEdgeFromInner(Size W, int crit, LL tNow);
  void scheduleEdgeFromEven(Size u, Size w, Size E, int crit, LL tNow, LL edgeW2,
                             bool dedupFreshEven = false,
                             Size edgeIndex = scNull(),
                             LL sourceY = std::numeric_limits<LL>::max());
  void initializeEventStructures();                    // one per-search O(m) setup
  void eventSfmOnShrink(Size B);                        // outer formation: topology-neutral
  void eventSfmOnDissolve(Size B);                      // split into odd-cycle members
#ifdef EVENT_S5_TARGET_KEY_SHADOW
  void initializeS5TargetKeyShadow(int crit);
  void s5TargetKeyOnShrink(Size B);
  void s5TargetKeyOnDissolve(Size B);
  void s5TargetKeyAfterDissolve(Size B, LL tNow);
  void s5TargetKeyTrackTreeNode(Size node);
  void s5TargetKeyRevealFree(Size v, LL tNow);
  void s5TargetKeyRetireTree(Size root, LL tNow);
  void auditS5TargetKeys(int crit, LL tNow);
  void s5TargetKeyRememberTarget(Size target);
  void s5TargetKeyRememberOwner(Size owner);
  void s5TargetKeyRememberToken(Size token);
  Size s5TargetFrontierListId(Size vertex) const;
  void s5TargetFrontierReset();
  void s5TargetFrontierInvalidateList(Size listId);
  bool s5TargetKeyGap(LL slack, bool matched, int crit, LL* gap) const;
  Val s5TargetSfmKey(Size v) const;
  void s5TargetSfmSetKey(Size v, Val key);
  void s5TargetSfmDecreaseKey(Size v, Val key);
  void s5TargetSfmFindMin(Size v, Val* key, Size* at) const;
  void s5TargetSfmFindMinWithList(Size v, Val* key, Size* at,
                                  Size* listId) const;
  void s5TargetSfmSublistVerts(Size v, std::vector<Size>& out) const;
#ifdef EVENT_S5_SINGLETON_TARGET_FAST
  bool s5TargetSfmIsSingleton(Size v) const;
  void s5TargetSfmPromoteSingleton(Size v);
#endif
  void s5TargetKeyReindexPayload(Size v, Size owner, Size targetOwner,
                                 Size sourceToken, bool oldLive, bool newLive);
#ifdef EVENT_S5_TARGET_KEY_DRIVER_AUDIT
  void scheduleS5TargetMinEvent(Size target);
  void scheduleS5TargetPersistentFrontier(Size target);
#endif
#endif
#ifdef EVENT_SFM_SHADOW
  void sfmShadowBeginInner(Size B);                     // reset/reseed one INNER interval
#endif
#ifdef EVENT_SFM_TOPOLOGY_AUDIT
  bool sfmShadowListMatches(Size B) const;
#endif
  void resetPracticalEventQueue();
  void schedulePracticalEvent(LL time, long long a, long long b, int kind,
                              long long owner = -1);
  void scheduleShellBoundaryEvent(LL tNow);
  bool popPracticalEvent(BucketPQ::Event* event);
  void applyDualTo(LL newTime, LL* tNow);             // advance duals to absolute time
  void beginCompactEventLabels();
  void rememberCompactEventLabel(Size node);

  // ---- S5 event-driven PERSISTENT engine (paper S5; docs/Scaling_S5_DataStructures.md)
  // Same MWPM result as runPerfectBlossom but NO per-iteration rebuild: the
  // alternating forest persists; the search advances by processing scheduled
  // events (grow/formation/dissolution) in time order via a bucket PQ, and dual
  // adjustments are applied incrementally to labeled nodes only. Reuses the
  // blossom helpers (shrink/augmentFrom/rotateToBase/lcaOuter/eligible/slackZ).
  bool runPerfectS5();                  // S5 persistent-forest engine (no rebuild)
  // one persistent search: relabel ONCE, then grow-extend / delta-jump dual-adjust
  // / dissolve, relabel ONLY after augment or dissolve. Returns true if augmented.
  bool edmondsSearchS5(int crit);

  // ---- Phase 2/3: the Liquidationist scaling loop (paper Fig 5) ----
  // All blossom ids >= mN; vertices 0..mN-1 (+ dummy vertices appended). The
  // scaling loop works on the doubled-even scale weight mW[u][idx]. See
  // docs/Scaling_Paper_Exact_Spec.md.
  void liquidate(Size B);               // y(u)+=z(B)/2 forall u in B; z(B)=0; dissolve
  void separateInheritedForest(Size oldBegin, Size oldEnd);
  void finishInheritedForest();
  bool checkCanonicalPostDissolution(Size realN) const;
  bool checkCanonicalSeparatedForests(Size realN) const;
  bool canonicalDualObjective(Size realN, LL* objective) const;
  void liquidateInherited(Size B);
  void collectBlossomVerts(Size B, std::vector<Size>* out) const;  // leaves of B
  Size blossomSize(Size B) const;       // # original vertices in B

  // ---- HYBRID §4.2 (Gabow's O(m(B) n(B)^{3/4}) small-blossom dissolution) ----
  // Reuses the existing blossom forest (mBParent/mBCycle/mNodeLeaves) + shrink/dissolve/
  // liquidate + edmondsSearch. Built C-first (this dissolution), gated vs Lemma 4.1.
  // docs/Scaling_Hybrid_Exact_Spec.md is ground truth.
  // majorChild(B): the sub-blossom B1 of B with |B1| > |B|/2, or scNull if none (a node
  //   has at most one major child since two would each exceed half).
  Size majorChild(Size B) const;
  // isMajorPathRoot(R): R is NOT a major child of its parent (so it heads a major path).
  bool isMajorPathRoot(Size R) const;
  // gabowDissolveSmallBlossom(B): dissolve all old sub-blossoms in the subtree rooted at B
  //   (a maximal small inherited blossom) per §4.2. Returns false on guard failure.
  bool gabowDissolveSmallBlossom(Size B);
  // collectMajorPathRoots(B, out): postorder list of major-path roots in B's subtree.
  void collectMajorPathRoots(Size B, std::vector<Size>* out) const;
  // undissolvedPathChain(R, out): the undissolved blossoms on major path P(R), innermost
  //   first (B1 subset B2 subset ... = deepest-to-R). "Undissolved" = still in mBCycle
  //   structure with z(B)!=0 OR alive node (we treat a node B>=mN with mBParent==scNull or
  //   parent on the path as undissolved until dissolve() is called on it).
  void undissolvedPathChain(Size R, std::vector<Size>* out) const;
  // dismantleStage2(R): §4.2 DismantlePath STAGE 2 (the |F|=1 branch). omega = the unique
  //   free vertex; liquidate B1..Bl on P(R); PQSearch({omega}) Crit1 halting after exactly
  //   T=sum z(Bi)/2 adjustments (restores yz, finds no augmenting path). Reuses liquidate +
  //   edmondsSearch(CRIT1, maxAdjust=T). Returns false on guard failure.
  bool dismantleStage2(Size R);
  bool resolveOldShellBounds(Size C, Size D, Size pathRoot,
                             Size* Cstar, Size* Dstar) const;
  Size nextLiveOldParent(Size B, Size pathRoot) const;
  Size nextLiveOldMajorChild(Size B) const;
  void dissolveZeroCurrentRoots();
  // shellSearch(Cstar, Dstar): §4.2 ShellSearch (Fig 7). Sets mInShell for the atomic shell
  //   G(Cstar,Dstar) (= leaves(Cstar) minus leaves(Dstar)), runs a confined Criterion-1
  //   maximal-augmentation search in the eligible shell, then clears shell state.
  //   Dstar may be scNull (odd shell). Returns true if it augmented >=1 path.
  bool shellSearch(Size Cstar, Size Dstar, Size pathRoot = scNull(),
                   Size iterEpoch = 0, Size callOwner = 0);
  // dismantlePath(R): §4.2 DismantlePath(R) = STAGE 1 (atomic-shell ShellSearch loop while
  //   |F|>=2) then STAGE 2 (|F|=1 liquidate cleanup). Dissolves all old blossoms on P(R).
  bool dismantlePath(Size R);
  // #free vertices in shell G(Cstar,Dstar) (leaves(Cstar)\leaves(Dstar) that are unmatched).
  Size shellFreeCount(Size Cstar, Size Dstar) const;
  bool runScaling();                    // full Fig-5 loop; fills mMate (real verts)
  bool runLiquidationistLeveling(
      const std::vector<std::vector<Size> >& maximalOldSmall, Size realN);
  bool runFaithfulSingleScale();        // ISOLATION TEST: one scale @ full weight w=2*what,
                                        // dummies+SearchOne(Crit2)+Finalization(Crit3) to
                                        // optimum. Validates the re-opt mechanism vs LEMON.

public:
#ifdef HYBTEST
  void hybTestInitSearchState(const Graph& graph) {
    init(graph);
    Size cap = 2 * mN + 1;
    mBParent.assign(cap, scNull()); mBBase.assign(cap, scNull());
    mBCycle.assign(cap, std::vector<Size>());
    mBCycEdgeU.assign(cap, std::vector<Size>());
    mBCycEdgeV.assign(cap, std::vector<Size>());
    mNodeLeaves.assign(cap, std::vector<Size>());
    mZ.assign(cap, 0); mY.assign(mN, mMaxW);
    for (Size v = 0; v < mN; ++v) {
      mBBase[v] = v; mNodeLeaves[v].assign(1, v);
    }
#ifdef PS5_LEAF_SPLICE
    scalingLfInit(mN, cap);
#endif
    mMate.assign(mN, scNull()); mNumNodes = mN;
    mNLabel.assign(cap, LBL_FREE_OF_TREE);
    mNParent.assign(cap, scNull()); mNParEdgeU.assign(cap, scNull());
    mNParEdgeV.assign(cap, scNull()); mNRoot.assign(cap, scNull());
  }

  // Test-only probe for EdmondsSearch(F): initialize a graph, perform only the
  // initially eligible maximal augmentation round, and return its cardinality.
  // A null root list means all free vertices; otherwise exactly the listed free
  // vertices seed alternating trees.
  Size hybTestExplicitRootCard(const Graph& graph,
                               const std::vector<Size>* roots) {
    hybTestInitSearchState(graph);
    edmondsSearch(CRIT1, /*maximalAugment=*/true,
                  /*maxAdjust=*/0, /*pqDelta=*/false,
                  /*maxDualUnits=*/0, roots);
    Size card = 0;
    for (Size v = 0; v < mN; ++v)
      if (mMate[v] != scNull() && v < mMate[v]) ++card;
    return card;
  }

  // Test-only shell-MCM oracle hook. Tests use uniform-weight graphs, so every
  // shell edge is Criterion-1 eligible at initialization. Vertices outside the
  // supplied region must remain unmatched.
  Size hybTestShellMcmCard(const Graph& graph,
                           const std::vector<Size>& shellVertices) {
    hybTestInitSearchState(graph);
    mInShell.assign(mN, 0); mInDstar.assign(mN, 0); mShellMarked.clear();
    for (Size i = 0; i < shellVertices.size(); ++i) {
      Size v = shellVertices[i];
      if (v < mN && !mInShell[v]) {
        mInShell[v] = 1; mShellMarked.push_back(v);
      }
    }
    mShellActive = true; mShellCstar = scNull(); mShellDstar = scNull();
    SearchRequest request;
    request.criterion = CRIT1; request.region = REGION_SHELL;
    request.maximalAugment = true; request.maxAdjustments = 0;
    request.dualBudget = 0; request.shellTranslation = true;
    request.stoppingRule = STOP_AT_SHELL_EVENT;
    executeSearch(request, BACKEND_SCAN);
    Size card = 0;
    for (Size v = 0; v < mN; ++v)
      if (mMate[v] != scNull() && v < mMate[v]) ++card;
    mShellActive = false; mShellMarked.clear();
    std::fill(mInShell.begin(), mInShell.end(), 0);
    return card;
  }

  // Same Figure-7 shell-MCM probe, but begin from a caller-supplied matching.
  // This is the state that matters inside DismantlePath: M need not be empty
  // when an Augmentation step starts.  Reject malformed/cross-shell seed
  // matchings by returning scNull().
  Size hybTestShellMcmCardFromMate(
      const Graph& graph, const std::vector<Size>& shellVertices,
      const std::vector<Size>& initialMate) {
    hybTestInitSearchState(graph);
    if (initialMate.size() != mN) return scNull();
    mInShell.assign(mN, 0); mInDstar.assign(mN, 0); mShellMarked.clear();
    for (Size i = 0; i < shellVertices.size(); ++i) {
      Size v = shellVertices[i];
      if (v < mN && !mInShell[v]) {
        mInShell[v] = 1; mShellMarked.push_back(v);
      }
    }
    for (Size v = 0; v < mN; ++v) {
      Size w = initialMate[v];
      if (w == scNull()) continue;
      if (w >= mN || initialMate[w] != v || !mInShell[v] || !mInShell[w])
        return scNull();
      bool edgeFound = false;
      for (Size i = 0; i < mAdj[v].size(); ++i)
        if (mAdj[v][i] == w) { edgeFound = true; break; }
      if (!edgeFound) return scNull();
    }
    mMate = initialMate;
    mShellActive = true; mShellCstar = scNull(); mShellDstar = scNull();
    SearchRequest request;
    request.criterion = CRIT1; request.region = REGION_SHELL;
    request.maximalAugment = true; request.maxAdjustments = 0;
    request.dualBudget = 0; request.shellTranslation = true;
    request.stoppingRule = STOP_AT_SHELL_EVENT;
    executeSearch(request, BACKEND_SCAN);
    Size card = 0;
    for (Size v = 0; v < mN; ++v)
      if (mMate[v] != scNull() && v < mMate[v]) ++card;
    mShellActive = false; mShellMarked.clear();
    std::fill(mInShell.begin(), mInShell.end(), 0);
    return card;
  }

  // Test-only oracle hook for the event-backed Figure 7 MCM. This synthetic
  // shell has no inherited C*/D* boundaries, so it exercises only the exact
  // zero-adjustment eligible-shell augmentation phase.
  Size hybTestEventShellMcmCard(const Graph& graph,
                                const std::vector<Size>& shellVertices,
                                bool* symmetric) {
    hybTestInitSearchState(graph);
    mIncrShrink = true;  // production event searches use contracted parents
    mInShell.assign(mN, 0); mInDstar.assign(mN, 0); mShellMarked.clear();
    for (Size i = 0; i < shellVertices.size(); ++i) {
      Size v = shellVertices[i];
      if (v < mN && !mInShell[v]) {
        mInShell[v] = 1; mShellMarked.push_back(v);
      }
    }
    mShellActive = true; mShellCstar = scNull(); mShellDstar = scNull();
    eventShellMcmAtCurrentDual();
    bool mateOk = true;
    for (Size v = 0; v < mN; ++v) {
      Size w = mMate[v];
      if (w != scNull() && (w >= mN || mMate[w] != v)) {
        mateOk = false; break;
      }
    }
    if (symmetric) *symmetric = mateOk;
    Size card = 0;
    if (mateOk) for (Size v = 0; v < mN; ++v)
      if (mMate[v] != scNull() && v < mMate[v]) ++card;
    mShellActive = false; mShellCstar = scNull(); mShellDstar = scNull();
    mShellMarked.clear(); std::fill(mInShell.begin(), mInShell.end(), 0);
    return card;
  }

  Size hybTestSearchContractCard(const Graph& graph,
                                 const std::vector<Size>* roots,
                                 SearchBackend backend, LL dualBudget,
                                 bool exhaustBudget, LL* dualUnits,
                                 ResultStatus* status,
                                 bool batchSameTimeAugment = false,
                                 bool maximalAugment = false,
                                 std::vector<Size>* mateOut = NULL,
                                 int criterion = CRIT1,
                                 LL maxAdjustments = -1) {
    hybTestInitSearchState(graph);
    SearchRequest request;
    request.roots = roots; request.criterion = criterion;
    request.dualBudget = dualBudget; request.maxAdjustments = maxAdjustments;
    request.batchEventFreeTime = true; request.exhaustBudget = exhaustBudget;
    request.batchSameTimeAugment = batchSameTimeAugment;
    request.maximalAugment = maximalAugment;
    request.stoppingRule = dualBudget >= 0 ? STOP_AT_DUAL_BUDGET
                                           : STOP_AFTER_AUGMENTATION;
    SearchResult result = executeSearch(request, backend);
    if (dualUnits) *dualUnits = result.dualUnits;
    if (status) *status = result.status;
    if (mateOut) *mateOut = mMate;
    Size card = 0;
    for (Size v = 0; v < mN; ++v)
      if (mMate[v] != scNull() && v < mMate[v]) ++card;
    return card;
  }

  // Test-only Figure-3 probe.  The fixture supplies a Property-2 state in
  // which y==0, every initially matched edge has weight 0 (slack 0), and
  // every other edge has weight 1 (slack -2).  Consequently every graph edge
  // is Criterion-2 eligible before the call, while every edge flipped by an
  // augmenting path becomes ineligible afterward, exactly as in Lemma 2.3.
  // Perform only SearchOne's augmentation step; the caller independently
  // checks that the symmetric difference is a vertex-disjoint augmenting-path
  // set and that no eligible augmenting path remains.
  Size hybTestSearchOneCardFromMate(const Graph& graph,
                                    const std::vector<Size>& initialMate,
                                    SearchBackend backend,
                                    std::vector<Size>* mateOut,
                                    ResultStatus* status,
                                    std::vector<Size>* outerOut = NULL) {
    hybTestInitSearchState(graph);
    if (initialMate.size() != mN) {
      if (status) *status = STATUS_INVARIANT_VIOLATION;
      return scNull();
    }
    for (Size v = 0; v < mN; ++v) {
      Size w = initialMate[v];
      if (w == scNull()) continue;
      if (w >= mN || w == v || initialMate[w] != v) {
        if (status) *status = STATUS_INVARIANT_VIOLATION;
        return scNull();
      }
      bool edgeFound = false;
      for (Size i = 0; i < mAdj[v].size(); ++i)
        if (mAdj[v][i] == w) { edgeFound = true; break; }
      if (!edgeFound) {
        if (status) *status = STATUS_INVARIANT_VIOLATION;
        return scNull();
      }
    }
    mMate = initialMate;
    mY.assign(mN, 0);
    mHybridCanonical = true;
    mHybridReference = backend == BACKEND_SCAN;
    mAlgorithmIdentity = backend == BACKEND_SCAN
      ? ALGORITHM_HYBRID_REFERENCE : ALGORITHM_HYBRID_EVENT_EXPERIMENTAL;
    mIncrShrink = backend == BACKEND_EVENT;
    SearchRequest request;
    request.criterion = CRIT2;
    request.maximalAugment = true;
    request.maxAdjustments = 0;
    request.stoppingRule = STOP_AFTER_AUGMENTATION;
    SearchResult result = executeSearch(request, backend);
    if (status) *status = result.status;
    if (mateOut) *mateOut = mMate;
    if (outerOut) {
      outerOut->resize(mN);
      for (Size v = 0; v < mN; ++v) (*outerOut)[v] = outerOf(v);
    }
    Size card = 0;
    for (Size v = 0; v < mN; ++v)
      if (mMate[v] != scNull() && v < mMate[v]) ++card;
    return card;
  }

  // Deterministic Section-5 matched-grow timing probe.  The unmatched edge
  // 0--1 is Criterion-2 eligible immediately, while 1--2 is matched with
  // slack -2.  Thus vertex 1 becomes inner at time zero, but its mate 2 must
  // remain unreached for two dual units; only then may the matched edge grow.
  // Run those units as two distinct SearchOne calls, as Figure 6 does.
  bool hybTestMatchedGrowTiming(const Graph& graph, SearchBackend backend,
                                std::vector<LL>* dualsOut) {
    hybTestInitSearchState(graph);
    if (mN != 4 || mAdj[0].size() != 1 || mAdj[1].size() != 2 ||
        mAdj[2].size() != 1 || mAdj[3].size() != 0)
      return false;
    mY.assign(mN, 1);
    mMate[1] = 2; mMate[2] = 1;
    mHybridCanonical = true;
    mHybridReference = backend == BACKEND_SCAN;
    mAlgorithmIdentity = backend == BACKEND_SCAN
      ? ALGORITHM_HYBRID_REFERENCE : ALGORITHM_HYBRID_EVENT_EXPERIMENTAL;
    mIncrShrink = backend == BACKEND_EVENT;
    std::vector<Size> roots(1, 0);
    for (int pass = 0; pass < 2; ++pass) {
      SearchRequest request;
      request.roots = &roots; request.criterion = CRIT2;
      request.dualBudget = 1; request.maxAdjustments = 1;
      request.batchEventFreeTime = true; request.exhaustBudget = true;
      request.stoppingRule = STOP_AT_DUAL_BUDGET;
      SearchResult result = executeSearch(request, backend);
      if (result.status != STATUS_SUCCESS || result.augmented ||
          result.dualUnits != 1 || !result.reachedBudget)
        return false;
    }
    if (dualsOut) *dualsOut = mY;
    return true;
  }

  // Synthetic DismantlePath stage-2 fixture. R={0,1,2} has the unique free
  // vertex omega=0 and z(R)=2, while vertices 3 and 4 form an unrelated eligible
  // free edge. Correct singleton scoping must leave that outside edge unmatched.
  bool hybTestStage2Singleton(const Graph& graph) {
    hybTestInitSearchState(graph);
    if (mN < 5) return false;
    Size R = mN;
    mNumNodes = R + 1;
    mBBase[R] = 0; mZ[R] = 2;
    mBCycle[R].push_back(0); mBCycle[R].push_back(1); mBCycle[R].push_back(2);
    mNodeLeaves[R].push_back(0); mNodeLeaves[R].push_back(1); mNodeLeaves[R].push_back(2);
    mBParent[0] = mBParent[1] = mBParent[2] = R;
    mMate[1] = 2; mMate[2] = 1;
    LL dualBefore = 0;
    for (Size v = 0; v < mN; ++v) dualBefore += mY[v];
    dualBefore += mZ[R];                 // floor(|R|/2)=1
    // Enter through the real DismantlePath transition: it must observe |F|=1
    // and dispatch stage 2 itself (the former fixture called stage 2 directly).
#ifdef HYBPROF
    long long before = gHybStage2Calls;
#endif
    bool ok = dismantlePath(R);
    LL dualAfter = 0;
    for (Size v = 0; v < mN; ++v) dualAfter += mY[v];
    return ok && mLastSearchDualUnits == 1 && mY[0] == mMaxW &&
           dualAfter == dualBefore &&
           mZ[R] == 0 && mBParent[0] == scNull() &&
           mMate[1] == 2 && mMate[2] == 1 &&
           mMate[3] == scNull() && mMate[4] == scNull()
#ifdef HYBPROF
           && gHybStage2Calls == before + 1
#endif
           ;
  }

  // Nested major-path transition fixture: R={B,3,4}, B={0,1,2}, with B the
  // major child of R and omega=0 the unique free vertex in the whole path.
  bool hybTestStage2Nested(const Graph& graph) {
    hybTestInitSearchState(graph);
    if (mN < 8) return false;
    Size B = mN, R = mN + 1;
    mNumNodes = R + 1;
    mBBase[B] = mBBase[R] = 0;
    mZ[B] = 2; mZ[R] = 4;
    mBCycle[B].push_back(0); mBCycle[B].push_back(1); mBCycle[B].push_back(2);
    mNodeLeaves[B].push_back(0); mNodeLeaves[B].push_back(1); mNodeLeaves[B].push_back(2);
    mBParent[0] = mBParent[1] = mBParent[2] = B;
    mBCycle[R].push_back(B); mBCycle[R].push_back(3); mBCycle[R].push_back(4);
    for (Size v = 0; v <= 4; ++v) mNodeLeaves[R].push_back(v);
    mBParent[B] = mBParent[3] = mBParent[4] = R;
    mMate[1] = 2; mMate[2] = 1; mMate[3] = 4; mMate[4] = 3;
    LL dualBefore = 0;
    for (Size v = 0; v < mN; ++v) dualBefore += mY[v];
    dualBefore += mZ[B] + 2 * mZ[R];     // floor(3/2)=1, floor(5/2)=2
#ifdef HYBPROF
    long long before = gHybStage2Calls;
#endif
    bool ok = dismantlePath(R);
    LL dualAfter = 0;
    for (Size v = 0; v < mN; ++v) dualAfter += mY[v];
    bool yMonotone = true;
    for (Size v = 0; v <= 4; ++v) if (mY[v] < mMaxW) yMonotone = false;
    return ok && mLastSearchDualUnits == 3 && mY[0] == mMaxW &&
           yMonotone && dualAfter == dualBefore &&
           mZ[B] == 0 && mZ[R] == 0 &&
           mBParent[0] == scNull() && mBParent[3] == scNull() &&
           mMate[1] == 2 && mMate[2] == 1 &&
           mMate[3] == 4 && mMate[4] == 3 &&
           mMate[5] == scNull() && mMate[6] == scNull()
#ifdef HYBPROF
           && gHybStage2Calls == before + 1
#endif
           ;
  }

  // Unit-test hook for the HYBRID §4.2 major-path decomposition. Builds a synthetic blossom
  // forest directly into the members (no solve needed), runs collectMajorPathRoots, and
  // self-checks the structural invariants. Returns true iff all invariants hold.
  // Forest passed as: parent[i] (scNull=root), cyc[i] (sub-blossoms in order), leafCount[i].
  bool hybTestDecomp(Size nVerts, Size nNodes,
                     const std::vector<Size>& parent,
                     const std::vector<std::vector<Size> >& cyc,
                     const std::vector<Size>& leafCount,
                     Size topB) {
    mN = nVerts; mNumNodes = nNodes;
    mBParent = parent; mBCycle = cyc;
    mNodeLeaves.assign(nNodes, std::vector<Size>());
    for (Size i = 0; i < nNodes; ++i) mNodeLeaves[i].assign(leafCount[i], 0); // size only
    // (1) majorChild uniqueness + correctness: at most one child > |B|/2.
    for (Size b = nVerts; b < nNodes; ++b) {
      Size half = leafCount[b] / 2, found = scNull(), cnt = 0;
      for (Size i = 0; i < cyc[b].size(); ++i)
        if (leafCount[cyc[b][i]] > half) { ++cnt; found = cyc[b][i]; }
      if (cnt > 1) return false;                         // can't have 2 major children
      if (majorChild(b) != found) return false;          // matches our computation
    }
    // (2) collectMajorPathRoots: postorder + every non-trivial node lies on exactly ONE
    //     major path (i.e. each node maps to exactly one major-path root by walking up
    //     through major-child links until a major-path root).
    std::vector<Size> roots; collectMajorPathRoots(topB, &roots);
    // postorder check: a root's index in `roots` comes after all its sub-tree roots.
    // coverage: walk each non-trivial node up via "am I my parent's major child?" to its
    // path root, assert that root is in `roots` exactly once.
    std::vector<char> inRoots(nNodes, 0);
    for (Size i = 0; i < roots.size(); ++i) { if (inRoots[roots[i]]) return false; inRoots[roots[i]] = 1; }
    for (Size b = nVerts; b < nNodes; ++b) {
      // only consider nodes inside topB's subtree (reachable down from topB)
      Size pr = b;
      while (!isMajorPathRoot(pr)) pr = mBParent[pr];
      // pr is b's major-path root; it must be a collected root (if b is under topB)
      // (we accept nodes outside topB's subtree having pr not in roots; skip those by
      //  checking topB-ancestry cheaply: walk up to a forest root, must reach topB.)
      Size anc = b; bool underTop = false;
      while (anc != scNull() && anc < nNodes) { if (anc == topB) { underTop = true; break; } anc = mBParent[anc]; }
      if (underTop && !inRoots[pr]) return false;
    }
    return true;
  }

  // Deterministic Section-5.1 topology gate. The initial current forest is
  // B11=(B10,5,6), B10=(B9,3,4), B9=(0,1,2), plus roots 7 and 8.
  bool hybTestSplitOnlySfmLifecycle() {
    const Size B0 = 9, B1 = 10, B2 = 11, cap = 12;
    mCounters = OperationCounters();
    mN = 9; mNumNodes = cap;
    mBParent.assign(cap, scNull());
    mBBase.assign(cap, scNull());
    mBCycle.assign(cap, std::vector<Size>());
    mNodeLeaves.assign(cap, std::vector<Size>());
    mBActive.assign(cap, 0);
    mZ.assign(cap, 0);
    for (Size v = 0; v < mN; ++v) {
      mBBase[v] = v;
      mNodeLeaves[v].assign(1, v);
    }
    mBCycle[B0].push_back(0); mBCycle[B0].push_back(1);
    mBCycle[B0].push_back(2);
    mBCycle[B1].push_back(B0); mBCycle[B1].push_back(3);
    mBCycle[B1].push_back(4);
    mBCycle[B2].push_back(B1); mBCycle[B2].push_back(5);
    mBCycle[B2].push_back(6);
    mBBase[B0] = mBBase[B1] = mBBase[B2] = 0;
    mBActive[B0] = mBActive[B1] = mBActive[B2] = 1;
    for (Size i = 0; i < mBCycle[B0].size(); ++i)
      mBParent[mBCycle[B0][i]] = B0;
    for (Size i = 0; i < mBCycle[B1].size(); ++i)
      mBParent[mBCycle[B1][i]] = B1;
    for (Size i = 0; i < mBCycle[B2].size(); ++i)
      mBParent[mBCycle[B2][i]] = B2;
#ifdef PS5_LEAF_SPLICE
    scalingLfInit(mN, cap);
    scalingLfOnShrink(B0); scalingLfOnShrink(B1); scalingLfOnShrink(B2);
#endif
    mOldForestSeparated = false;
    mEventStructuresActive = false;
    mEventSfmInitialized = false;
    mS5TargetKeyActive = false;
    mS5TargetKeyTopologyReady = false;
    mLastStatus = STATUS_SUCCESS;
    initializeEventStructures();
    if (mLastStatus != STATUS_SUCCESS ||
        mCounters.sfmInitialLists != 3 ||
        mCounters.sfmInitialPivots != 6 ||
        mSfmSplitPivots.size() <= B2 ||
        mSfmSplitPivots[B0] != std::vector<Size>({0, 1}) ||
        mSfmSplitPivots[B1] != std::vector<Size>({2, 3}) ||
        mSfmSplitPivots[B2] != std::vector<Size>({4, 5})) return false;
    std::vector<Size> got;
    mInnerSfm.sublistVerts(0, got);
    if (got != std::vector<Size>({0, 1, 2, 3, 4, 5, 6}) ||
        mInnerSfm.sameList(0, 7) || mInnerSfm.sameList(7, 8)) return false;

    eventSfmOnShrink(B2);
    if (!mInnerSfm.sameList(0, 6) ||
        mCounters.sfmFormationTopologyNoops != 1) return false;

    eventSfmOnDissolve(B2);
    if (!mInnerSfm.sameList(0, 4) || mInnerSfm.sameList(4, 5) ||
        mInnerSfm.sameList(5, 6)) return false;
#ifdef PS5_LEAF_SPLICE
    scalingLfOnDissolve(B2);
#endif
    eventSfmOnDissolve(B1);
    if (!mInnerSfm.sameList(0, 2) || mInnerSfm.sameList(2, 3) ||
        mInnerSfm.sameList(3, 4)) return false;
#ifdef PS5_LEAF_SPLICE
    scalingLfOnDissolve(B1);
#endif
    eventSfmOnDissolve(B0);
    if (mInnerSfm.sameList(0, 1) || mInnerSfm.sameList(1, 2)) return false;
#ifdef PS5_LEAF_SPLICE
    scalingLfOnDissolve(B0);
#endif
    return mCounters.sfmActualSplits == 6 &&
           mCounters.splitFindminWork == 6 &&
           mCounters.sfmSplitTopologyFailures == 0;
  }

  // Exact Section-5.2 timestamp formulas for inner dissolution and outer
  // formation. These quantities are audited independently of the equivalent
  // role-accumulator channel that currently drives y/z materialization.
  bool hybTestTimestampTransitions() {
    const Size B = 3, cap = 4;
    mN = 3; mNumNodes = cap; mTNow = 7; mLazyDuals = true;
    mBCycle.assign(cap, std::vector<Size>());
    mBCycle[B].push_back(0); mBCycle[B].push_back(1);
    mBCycle[B].push_back(2);
    mBBase.assign(cap, scNull());
    for (Size v = 0; v < mN; ++v) mBBase[v] = v;
    mBBase[B] = 0;
    mNodeLeaves.assign(cap, std::vector<Size>());
    for (Size v = 0; v < mN; ++v) mNodeLeaves[v].assign(1, v);
#ifdef PS5_LEAF_SPLICE
    scalingLfInit(mN, cap); scalingLfOnShrink(B);
#endif
    mY.assign(mN, 0); mZ.assign(cap, 0);
    mY0 = mY; mZ0 = mZ;
    mYAcc.assign(mN, 0); mZAcc.assign(cap, 0);
    mVStamp.assign(mN, 0); mNStamp.assign(cap, 0);
    mVRole.assign(mN, LBL_FREE_OF_TREE);
    mNRole.assign(cap, LBL_FREE_OF_TREE);
    mLazyEpochClock = 1;
    mLazyVEpoch.assign(mN, mLazyEpochClock);
    mLazyNEpoch.assign(cap, 0); mLazyNEpoch[B] = mLazyEpochClock;
    mLazyActiveVertices.clear();
    for (Size v = 0; v < mN; ++v) mLazyActiveVertices.push_back(v);
    mLazyActiveNodes.assign(1, B);
    mTRoot.assign(cap, -1); mTIn.assign(cap, -1);
    mTOut.assign(cap, -1); mDelta.assign(cap, 0);
    mNRole[B] = LBL_ODD; mNStamp[B] = 2; mTIn[B] = 2;
    mDelta[B] = 11;
    // Deliberately poison the child values: dissolution must overwrite them
    // with Delta(B)+(t_now-t_in(B)), not increment their previous values.
    mDelta[0] = 101; mDelta[1] = 202; mDelta[2] = 303;
    mNLabel.assign(cap, LBL_FREE_OF_TREE);
    mNLabel[0] = LBL_ODD; mNLabel[1] = LBL_EVEN;
    tsOnDissolveInner(B);
    if (mLastStatus != STATUS_SUCCESS) return false;
    for (Size c = 0; c < 3; ++c)
      if (mDelta[c] != 16 || mTRoot[c] != 7) return false;
    if (mTIn[0] != 7 || mTOut[1] != 7) return false;

    // Formation updates only formerly INNER constituents, then starts the
    // new blossom's OUTER interval at the current time.
    mTNow = 13;
    mVRole[0] = LBL_ODD; mVStamp[0] = 9; mTIn[0] = 9;
    mVRole[1] = LBL_EVEN; mVStamp[1] = 9;
    mVRole[2] = LBL_FREE_OF_TREE; mVStamp[2] = 9;
    mDelta[0] = 4; mDelta[1] = 5; mDelta[2] = 6;
    mNRole[B] = LBL_FREE_OF_TREE; mNStamp[B] = 13;
    mTRoot[B] = -1; mTOut[B] = -1;
    tsOnFormOuter(B);
    return mLastStatus == STATUS_SUCCESS &&
           mDelta[0] == 8 && mTOut[0] == 13 &&
           mDelta[1] == 5 && mDelta[2] == 6 &&
           mTRoot[B] == 13 && mTOut[B] == 13 &&
           mNRole[B] == LBL_EVEN;
  }
#endif
  // S5 engine entry (for gating vs runPerfectBlossom/LEMON + timing).
  bool computeS5(const Graph& graph, std::vector<Size>* mateVec,
                 Size* card, Val* edgWght) {
    init(graph);
    Size mN0 = mN;
    bool ok = runPerfectS5();
    mateVec->assign(mN0, cNullItm);
    Size c = 0; Val wsum = 0;
    if (ok) for (Size v = 0; v < mN0; ++v)
      if (mMate[v] != scNull() && mMate[v] < mN0 && mMate[v] > v) {
        Size mate = mMate[v]; (*mateVec)[v] = mate; (*mateVec)[mate] = v; ++c;
        LL best = 0; bool f = false;
        for (Size i = 0; i < mAdj[v].size(); ++i)
          if (mAdj[v][i] == mate && (!f || mWhat[v][i] > best)) { best = mWhat[v][i]; f = true; }
        if (f) wsum += (Val)best;
      }
    if (card) *card = c; if (edgWght) *edgWght = wsum;
    return ok;
  }

  // PQSearch (delta-jump) engine entry: same MWPM result as compute(), but the
  // dual adjustment jumps by the largest safe delta, so the relabel runs per
  // EVENT not per UNIT step. Gated vs compute()/LEMON for correctness + timed.
  bool computeEvent(const Graph& graph, std::vector<Size>* mateVec,
                    Size* card, Val* edgWght) {
    init(graph);
    setIncrShrink(true); // event forest persists across a shrink transition
    Size mN0 = mN;
    bool ok = runPerfectBlossomEvent();
    // Never publish a malformed experimental result. In particular, an
    // incremental blossom transition must not be allowed to hide an
    // asymmetric mate array behind the output loop's `mate > v` filter.
    if (ok) for (Size v = 0; v < mN0; ++v) {
      Size mate = mMate[v];
      if (mate == scNull() || mate >= mN0 || mMate[mate] != v) {
        mLastStatus = STATUS_INVARIANT_VIOLATION;
        ok = false;
        break;
      }
    }
    mateVec->assign(mN0, cNullItm);
    Size c = 0; Val wsum = 0;
    if (ok) for (Size v = 0; v < mN0; ++v)
      if (mMate[v] != scNull() && mMate[v] < mN0 && mMate[v] > v) {
        Size mate = mMate[v]; (*mateVec)[v] = mate; (*mateVec)[mate] = v; ++c;
        LL best = 0; bool f = false;
        for (Size i = 0; i < mAdj[v].size(); ++i)
          if (mAdj[v][i] == mate && (!f || mWhat[v][i] > best)) { best = mWhat[v][i]; f = true; }
        if (f) wsum += (Val)best;
      }
    if (card) *card = c; if (edgWght) *edgWght = wsum;
    return ok;
  }

  bool computePersist(const Graph& graph, std::vector<Size>* mateVec,
                      Size* card, Val* edgWght) {
    init(graph);
    Size mN0 = mN;
    bool ok = runPerfectBlossomPersist();
    mateVec->assign(mN0, cNullItm);
    Size c = 0; Val wsum = 0;
    if (ok) for (Size v = 0; v < mN0; ++v)
      if (mMate[v] != scNull() && mMate[v] < mN0 && mMate[v] > v) {
        Size mate = mMate[v]; (*mateVec)[v] = mate; (*mateVec)[mate] = v; ++c;
        LL best = 0; bool f = false;
        for (Size i = 0; i < mAdj[v].size(); ++i)
          if (mAdj[v][i] == mate && (!f || mWhat[v][i] > best)) { best = mWhat[v][i]; f = true; }
        if (f) wsum += (Val)best;
      }
    if (card) *card = c; if (edgWght) *edgWght = wsum;
    return ok;
  }

  bool computePQ(const Graph& graph, std::vector<Size>* mateVec,
                 Size* card, Val* edgWght) {
    init(graph);
    Size mN0 = mN;
    bool ok = runPerfectBlossomPQ();
    mateVec->assign(mN0, cNullItm);
    Size c = 0; Val wsum = 0;
    if (ok) for (Size v = 0; v < mN0; ++v)
      if (mMate[v] != scNull() && mMate[v] < mN0 && mMate[v] > v) {
        Size mate = mMate[v]; (*mateVec)[v] = mate; (*mateVec)[mate] = v; ++c;
        LL best = 0; bool f = false;
        for (Size i = 0; i < mAdj[v].size(); ++i)
          if (mAdj[v][i] == mate && (!f || mWhat[v][i] > best)) { best = mWhat[v][i]; f = true; }
        if (f) wsum += (Val)best;
      }
    if (card) *card = c; if (edgWght) *edgWght = wsum;
    return ok;
  }

  // expose for Gate-1 testing: run the no-blossom max-weight PERFECT search.
  // Returns true on success (perfect matching found). On failure, exactly one
  // of *blossomNeeded / *noPerfect is set.
  bool computeNoBlossom(const Graph& graph,
                        std::vector<Size>* mateVec,
                        Size* card, Val* edgWght,
                        bool* blossomNeeded, bool* noPerfect);

  // ISOLATION TEST entry: build from graph, run the single-scale faithful
  // re-optimization (Crit2 + dummies + Finalization Crit3) and report card/weight.
  bool computeFaithful(const Graph& graph, std::vector<Size>* mateVec,
                       Size* card, Val* edgWght) {
    init(graph);
    Size mN0 = mN;
    bool ok = runFaithfulSingleScale();
    mateVec->assign(mN0, cNullItm);
    Size c = 0; Val wsum = 0;
    if (ok) for (Size v = 0; v < mN0; ++v)
      if (mMate[v] != scNull() && mMate[v] < mN0 && mMate[v] > v) {
        Size mate = mMate[v]; (*mateVec)[v] = mate; (*mateVec)[mate] = v; ++c;
        LL best = 0; bool f = false;
        for (Size i = 0; i < mAdj[v].size(); ++i)
          if (mAdj[v][i] == mate && (!f || mWhat[v][i] > best)) { best = mWhat[v][i]; f = true; }
        if (f) wsum += (Val)best;
      }
    if (card) *card = c; if (edgWght) *edgWght = wsum;
    return ok;
  }

  // LAYER 2 (TIME-BOUND): warm-started bounded scaling loop on the FAST engine.
  // For the paper's O(m sqrt n log nN): log(nN) scales, each WARM-STARTED from the
  // previous scale's matching + rescaled duals (so each scale repairs in O(few)
  // augments, not a full re-solve). The scale weight is the top-i bits of wbar; the
  // dual carries as 2*y'. Returns the exact MWPM. Gated vs computeScaling (the
  // correct re-solve oracle) + brute + LEMON. mUseBoundedLoop default OFF.
  bool computeBoundedScaling(const Graph& graph, std::vector<Size>* mateVec,
                             Size* card, Val* edgWght) {
    init(graph);
    Size n = mN;
    if (n == 0) { if (mateVec) mateVec->clear(); if (card)*card=0; if(edgWght)*edgWght=0; return true; }
    if (n % 2 != 0) { return false; }
    Size nbits = mNumScales;
    // build a static edge list (u<v) with the wbar value, indexed for fast addEdge.
    std::vector<Size> eu, ev; std::vector<LL> ewbar;
    for (Size v = 0; v < n; ++v)
      for (Size i = 0; i < mAdj[v].size(); ++i) {
        Size w = mAdj[v][i];
        if (w > v) { eu.push_back(v); ev.push_back(w); ewbar.push_back(mWbar[v][i]); }
      }
    // NOTE: the fast engine's free sentinel is NIL()==(Size)-1, which DIFFERS from
    // matchbox cNullItm==SIZE_MAX-1. Keep `mate` in the fast engine's NIL form for
    // the whole loop (the seed is passed straight through to solveFrom); translate
    // to cNullItm only at final extraction.
    const Size FNIL = ScalingFastMWPM::NIL();
    std::vector<Size> mate(n, FNIL);
    std::vector<LL> y(n, 0);             // carried duals (used iff mCarryDuals)
    for (Size scale = 1; scale <= nbits; ++scale) {
      Size bitpos = nbits - scale;       // MSB-first
      // scale weight for the fast engine = top-`scale` bits of wbar = wbar >> bitpos.
      ScalingFastMWPM fe; fe.setUseEvent(true); fe.init(n);
      if (mCarryDuals) fe.setReoptDom(true);   // carried duals: repair domination violations
      LL maxw = 0;
      // EXPERIMENT (mCarryReduced): reduced-cost step-4. Feed the engine the REDUCED
      // weight wb - y(u) - y(v) (carried duals from the PREVIOUS scale subtracted out),
      // shifted by a per-edge constant to stay non-negative. For a PERFECT matching this
      // shifts every perfect matching's weight by the same constant Sum_v y(v), so the
      // argmax is UNCHANGED -- but the numbers the search sees are SMALL (slacks, not
      // full-magnitude weights). Tests whether reduced costs (a) kill the large-N blowup
      // and (b) unfreeze the carried-dual loop. yShift accumulates the subtracted potential.
      std::vector<LL> wfed(eu.size());
      LL yshift = 0;
      if (mCarryReduced) {
        for (Size v = 0; v < n; ++v) yshift += y[v];  // y = carried duals (prev scale)
        for (Size e = 0; e < eu.size(); ++e) {
          LL wb = ewbar[e] >> bitpos;
          // reduced cost = wb - y(u) - y(v); can be negative -> shift by +maxAbs later.
          wfed[e] = wb - y[eu[e]] - y[ev[e]];
        }
        // shift all reduced weights so the minimum is 0 (keeps them non-negative; a
        // uniform additive shift on a PERFECT matching is constant -> argmax preserved).
        LL mn = 0; for (Size e = 0; e < eu.size(); ++e) if (wfed[e] < mn) mn = wfed[e];
        for (Size e = 0; e < eu.size(); ++e) {
          wfed[e] -= mn;                              // now >= 0
          fe.addEdge(eu[e], ev[e], wfed[e]);
          if (wfed[e] > maxw) maxw = wfed[e];
        }
      } else {
        for (Size e = 0; e < eu.size(); ++e) {
          LL wb = ewbar[e] >> bitpos;
          fe.addEdge(eu[e], ev[e], wb);
          if (wb > maxw) maxw = wb;
        }
      }
      // TWO warm-start dual policies (both warm-start the MATCHING M):
      //  * FRESH (default): y = maxw (full domination). solveFrom's reopt frees the
      //    untight matched edges and re-augments to the EXACT max-weight matching.
      //    Provably correct (264/264) but ~n/2 augments/scale -- the reopt frees ~all
      //    edges because fresh domination makes every matched edge untight, so it
      //    degenerates to a fresh re-solve (Theta(n)/scale, no speedup).
      //  * CARRIED (mCarryDuals): y = 2y'+3 (paper Fig-5 step 2). Matched edges stay
      //    near-tight (slack 4-6) so FEW are freed/scale (~n/2 -> a handful). The
      //    newly-exposed bit drives the improving edge's slack NEGATIVE, which reopt's
      //    domination-violation rule (b) detects and frees -- this is what un-freezes
      //    the matching while keeping augments/scale small (the O(sqrt n) regime).
      std::vector<LL> sy(n);
      if (mCarryReduced)    for (Size v = 0; v < n; ++v) sy[v] = maxw;  // small (reduced wts)
      else if (mCarryDuals) for (Size v = 0; v < n; ++v) sy[v] = 2 * y[v] + 3;
      else                  for (Size v = 0; v < n; ++v) sy[v] = maxw;
      std::vector<Size> nm;
      LL g = fe.solveFrom(mate, sy, &nm);
      if (g < 0) { return false; }
      mate = nm;
      if (mCarryReduced) {
        // carry the ORIGINAL-space duals for the next scale's reduced-cost subtraction:
        // y_new(v) = y_old(v) + dualY_reduced(v). (dualY is in reduced space; adding the
        // prior potential reconstructs an original-space dual that the next scale subtracts.)
        for (Size v = 0; v < n; ++v) y[v] = y[v] + fe.dualY(v);
      } else if (mCarryDuals) for (Size v = 0; v < n; ++v) y[v] = fe.dualY(v);
    }
    // extract
    if (mateVec) mateVec->assign(n, cNullItm);
    Size c = 0; Val wsum = 0;
    for (Size v = 0; v < n; ++v)
      if (mate[v] != FNIL && mate[v] < n && mate[v] > v) {
        if (mateVec) { (*mateVec)[v] = mate[v]; (*mateVec)[mate[v]] = v; }
        ++c;
        LL best = 0; bool f = false;
        for (Size i = 0; i < mAdj[v].size(); ++i)
          if (mAdj[v][i] == mate[v] && (!f || mWhat[v][i] > best)) { best = mWhat[v][i]; f = true; }
        if (f) wsum += (Val)best;
      }
    if (card) *card = c; if (edgWght) *edgWght = wsum;
    return true;
  }

  // LAYER 2 (TIME-BOUND, S5): the scaling loop on the NO-RELABEL S5 engine. Each scale is
  // a fresh O(m sqrt n) max-weight perfect matching over the top-i bits of wbar; log(nN)
  // scales => O(m sqrt n log(nN)) = THE PAPER'S BOUND. (The S5 engine's per-level multi-
  // augment gives O(sqrt n) searches/scale, each O(m); the lazy duals make a dual
  // adjustment O(1); no-relabel shrink/augment + guarded role-sync keep per-edge work flat.)
  // Fresh per scale (no warm-start needed for the asymptotic bound; warm-start would only
  // shave constants). Returns exact MWPM. Gated vs computeScaling + brute + LEMON.
  bool computeS5Bounded(const Graph& graph, std::vector<Size>* mateVec,
                        Size* card, Val* edgWght) {
    init(graph);
    Size n = mN;
    if (n == 0) { if (mateVec) mateVec->clear(); if (card)*card=0; if(edgWght)*edgWght=0; return true; }
    if (n % 2 != 0) { return false; }
    Size nbits = mNumScales;
    std::vector<Size> eu, ev; std::vector<LL> ewbar;
    for (Size v = 0; v < n; ++v)
      for (Size i = 0; i < mAdj[v].size(); ++i) {
        Size w = mAdj[v][i];
        if (w > v) { eu.push_back(v); ev.push_back(w); ewbar.push_back(mWbar[v][i]); }
      }
    const Size FNIL = ScalingS5Engine::NIL();
    std::vector<Size> mate(n, FNIL);
    for (Size scale = 1; scale <= nbits; ++scale) {
      Size bitpos = nbits - scale;       // top-`scale` bits of wbar
      ScalingS5Engine fe;
      fe.setUseEvent(true); fe.setLazyEvent(true);
      fe.setNoRelabelShrink(true); fe.setNoRelabelAug(true);   // full fast path
      fe.init(n);
      for (Size e = 0; e < eu.size(); ++e)
        fe.addEdge(eu[e], ev[e], ewbar[e] >> bitpos);
      std::vector<Size> nm;
      LL g = fe.solve(&nm);
      if (g < 0) { return false; }
      mate = nm;
    }
    if (mateVec) mateVec->assign(n, cNullItm);
    Size c = 0; Val wsum = 0;
    for (Size v = 0; v < n; ++v)
      if (mate[v] != FNIL && mate[v] < n && mate[v] > v) {
        if (mateVec) { (*mateVec)[v] = mate[v]; (*mateVec)[mate[v]] = v; }
        ++c;
        LL best = 0; bool f = false;
        for (Size i = 0; i < mAdj[v].size(); ++i)
          if (mAdj[v][i] == mate[v] && (!f || mWhat[v][i] > best)) { best = mWhat[v][i]; f = true; }
        if (f) wsum += (Val)best;
      }
    if (card) *card = c; if (edgWght) *edgWght = wsum;
    return true;
  }

  // Carried multi-scale loop shared by the named scaling entries. It
  // re-optimizes via bounded search, dummies, and finalization. The Hybrid
  // branch implements Figure 6 semantics; the named Liquidationist branch
  // includes Figure 5's confined maximum-y step-6 leveling searches.
  bool computeScaling(const Graph& graph, std::vector<Size>* mateVec,
                      Size* card, Val* edgWght, bool useFast = false) {
    mLastStatus = STATUS_SUCCESS;
    mCounters = OperationCounters();
    if (mateVec == NULL || card == NULL || edgWght == NULL) {
      mLastStatus = STATUS_UNSUPPORTED_OR_OVERFLOWING_INPUT;
      return false;
    }
    Size inputN = 0, inputM = 0;
    graph.GetSize(&inputN, &inputM);
    (void)inputM;
    mateVec->assign(inputN, cNullItm);
    *card = 0;
    *edgWght = 0;
    init(graph);
    if (mLastStatus == STATUS_UNSUPPORTED_OR_OVERFLOWING_INPUT) return false;
    setUseFastEngine(useFast);
    Size mN0 = mN;
    bool ok = runScaling();
    mateVec->assign(mN0, cNullItm);
    Size c = 0; Val wsum = 0;
    if (ok) for (Size v = 0; v < mN0; ++v)
      if (mMate[v] != scNull() && mMate[v] < mN0 && mMate[v] > v) {
        Size mate = mMate[v]; (*mateVec)[v] = mate; (*mateVec)[mate] = v; ++c;
        LL best = 0; bool f = false;
        for (Size i = 0; i < mAdj[v].size(); ++i)
          if (mAdj[v][i] == mate && (!f || mWhat[v][i] > best)) { best = mWhat[v][i]; f = true; }
        if (f) wsum += (Val)best;
      }
    if (card) *card = c; if (edgWght) *edgWght = wsum;
    return ok;
  }

  bool computeLiquidationist(const Graph& graph, std::vector<Size>* mateVec,
                             Size* card, Val* edgWght) {
    mAlgorithmIdentity = ALGORITHM_LIQUIDATIONIST_FIGURE5;
    mDropPerfectExperimental = false; mSkipFinalizationExperimental = false;
    mHybridCanonical = false; mHybridReference = false;
    setFaithfulLoop(true); setFaithfulReweight(true);
    setBoundedSearch(true); setHybrid(false);
    setIncrShrink(false); setPersistSearch(false);
    setNoRescanGrow(false); setNoRescanDual(false); setEventSearch(false);
    mCarryM = false; mContinueAug = false; mIncrRelabel = false;
    setIncrDissolve(false); setMultiLevel(false); setBridgeFvr(false);
    setHybNoTranslate(false); setDynamicShell(false);
    return computeScaling(graph, mateVec, card, edgWght);
  }

  // HYBRID (Duan-Pettie-Su Figure 6 semantic entry point): the Liquidationist
  // scaling loop with the three Hybrid changes -- (1) no reweight (already a simulated no-op),
  // (2) SMALL inherited blossoms dissolved via Gabow's §4.2 O(m n^{3/4}) routine (vs liquidate),
  // (3) two-stage Free Vertex Reduction (sqrt n SearchOne(Crit2) then BucketSearch(Crit3) to
  // tau = sqrt n log n). Computes the exact max-weight PERFECT matching. Gated: testHybridLoop
  // The semantic and sanitizer promotion gates pass. The canonical path uses
  // the integrated practical Section-5 event backend; the exact Theorem-4.3
  // bound is still not claimed because standard UF/heap and repeated shell-MCM
  // remain implementation adaptations.
  bool computeHybrid(const Graph& graph, std::vector<Size>* mateVec, Size* card, Val* edgWght) {
    mAlgorithmIdentity = ALGORITHM_HYBRID_FIGURE6;
    mDropPerfectExperimental = false; mSkipFinalizationExperimental = false;
    mHybridCanonical = true; mHybridReference = false;
    setFaithfulLoop(true);     // the faithful multi-scale loop (Fig 5/6 structure)
    setBoundedSearch(true);    // tau-bounded FVR step (the two-stage path for Hybrid)
    setHybrid(true);           // §4.2 dissolution + larger tau + two-stage FVR
    setIncrShrink(true);
    setPersistSearch(true);
    setNoRescanGrow(true);
    setNoRescanDual(false);
    setFaithfulReweight(false); setEventSearch(true);
    mCarryM = false; mContinueAug = false; mIncrRelabel = false;
    setIncrDissolve(true); setMultiLevel(false); setBridgeFvr(false);
    setHybNoTranslate(false); setDynamicShell(true);
    setEventLazyClock(true); setEventBatchSameTime(false);
#ifdef SCALING_EVENT_FULL_INITIALIZATION
    setEventCompactInitialization(false);
#else
    setEventCompactInitialization(true);
#endif
#ifdef SCALING_EVENT_NO_BOUNDARY_REUSE
    setEventBoundaryReuse(false);
#else
    // Reuse only processes events tied at the current dual time. Any positive
    // time boundary forces the ordinary fresh initialization path below.
    setEventBoundaryReuse(true);
#endif
    return computeScaling(graph, mateVec, card, edgWght);
  }

  // Figure 6 / Section 4.2 semantic oracle using conservative scan searches.
  bool computeHybridReference(const Graph& graph, std::vector<Size>* mateVec,
                              Size* card, Val* edgWght) {
    mAlgorithmIdentity = ALGORITHM_HYBRID_REFERENCE;
    mDropPerfectExperimental = false; mSkipFinalizationExperimental = false;
    mHybridCanonical = true; mHybridReference = true;
    setFaithfulLoop(true); setBoundedSearch(true); setHybrid(true);
    setFaithfulReweight(false); setEventSearch(false);
    setIncrShrink(false); setPersistSearch(false);
    setNoRescanGrow(false); setNoRescanDual(false);
    mCarryM = false; mContinueAug = false; mIncrRelabel = false;
    setIncrDissolve(false); setMultiLevel(false); setBridgeFvr(false);
    setHybNoTranslate(false); setDynamicShell(true);
    return computeScaling(graph, mateVec, card, edgWght);
  }

  // Deprecated compatibility entry point for pre-promotion experiments. It
  // executes the same Figure-6/event policy as computeHybrid(), while retaining
  // its historical algorithm identity so archived benchmark drivers remain
  // distinguishable.
  bool computeHybridEventExperimental(const Graph& graph,
                                      std::vector<Size>* mateVec,
                                      Size* card, Val* edgWght) {
    mAlgorithmIdentity = ALGORITHM_HYBRID_EVENT_EXPERIMENTAL;
    mDropPerfectExperimental = false; mSkipFinalizationExperimental = false;
    mHybridCanonical = true; mHybridReference = false;
    setFaithfulLoop(true); setBoundedSearch(true); setHybrid(true);
    setFaithfulReweight(false); setEventSearch(true);
    setIncrShrink(true); setPersistSearch(true);
    setNoRescanGrow(true); setNoRescanDual(false);
    mCarryM = false; mContinueAug = false; mIncrRelabel = false;
    setIncrDissolve(true); setMultiLevel(false); setBridgeFvr(false);
    setHybNoTranslate(false); setDynamicShell(true);
    // Section 5.2 lazy timestamps apply to ordinary BucketSearch/PQSearch.
    // edmondsSearchEvent deliberately keeps Figure 7 shell requests on their
    // exact materialized C*/D* translation clock.
    setEventLazyClock(true); setEventBatchSameTime(false);
#ifdef SCALING_EVENT_FULL_INITIALIZATION
    // Reproducible performance-control build.  This changes only the named
    // experimental event backend and is reported by RunConfiguration.
    setEventCompactInitialization(false);
#else
    setEventCompactInitialization(true);
#endif
#ifdef SCALING_EVENT_NO_BOUNDARY_REUSE
    setEventBoundaryReuse(false);
#else
    setEventBoundaryReuse(true);
#endif
    return computeScaling(graph, mateVec, card, edgWght);
  }

  // Deliberately minimal non-perfect ablation requested for experimental
  // measurement. It runs the canonical Figure-6 machinery unchanged, but:
  //   (1) accepts odd-order input, and
  //   (2) after Finalization cannot find another augmenting path, publishes
  //       the remaining real-vertex matching instead of reporting failure.
  // This is NOT advertised as an exact unrestricted-MWM algorithm: it adds no
  // D1/free-dual-zero event and supplies no MWM complementary-slackness
  // certificate. Callers must compare its weight/cardinality to an exact MWM
  // oracle and report it as a heuristic/ablation.
  bool computeHybridDropPerfectExperimental(const Graph& graph,
                                            std::vector<Size>* mateVec,
                                            Size* card, Val* edgWght) {
    mAlgorithmIdentity = ALGORITHM_HYBRID_DROP_PERFECT_EXPERIMENTAL;
    mDropPerfectExperimental = true;
    mSkipFinalizationExperimental = false;
    mHybridCanonical = true; mHybridReference = false;
    setFaithfulLoop(true); setBoundedSearch(true); setHybrid(true);
    setFaithfulReweight(false); setEventSearch(true);
    setIncrShrink(true); setPersistSearch(true);
    setNoRescanGrow(true); setNoRescanDual(false);
    mCarryM = false; mContinueAug = false; mIncrRelabel = false;
    setIncrDissolve(true); setMultiLevel(false); setBridgeFvr(false);
    setHybNoTranslate(false); setDynamicShell(true);
    setEventLazyClock(true); setEventBatchSameTime(false);
#ifdef SCALING_EVENT_FULL_INITIALIZATION
    setEventCompactInitialization(false);
#else
    setEventCompactInitialization(true);
#endif
#ifdef SCALING_EVENT_NO_BOUNDARY_REUSE
    setEventBoundaryReuse(false);
#else
    setEventBoundaryReuse(true);
#endif
    return computeScaling(graph, mateVec, card, edgWght);
  }

  // More aggressive companion ablation: run all Figure-6 scales through
  // Perfection, then delete the final dummy edges and publish immediately,
  // without Criterion-3 Finalization. This is the variant that can actually
  // save the final searches; it is intentionally a heuristic and may be below
  // both unrestricted-MWM and maximum-cardinality weighted quality.
  bool computeHybridSkipFinalizationExperimental(const Graph& graph,
                                                  std::vector<Size>* mateVec,
                                                  Size* card, Val* edgWght) {
    mAlgorithmIdentity = ALGORITHM_HYBRID_SKIP_FINALIZATION_EXPERIMENTAL;
    mDropPerfectExperimental = true;
    mSkipFinalizationExperimental = true;
    mHybridCanonical = true; mHybridReference = false;
    setFaithfulLoop(true); setBoundedSearch(true); setHybrid(true);
    setFaithfulReweight(false); setEventSearch(true);
    setIncrShrink(true); setPersistSearch(true);
    setNoRescanGrow(true); setNoRescanDual(false);
    mCarryM = false; mContinueAug = false; mIncrRelabel = false;
    setIncrDissolve(true); setMultiLevel(false); setBridgeFvr(false);
    setHybNoTranslate(false); setDynamicShell(true);
    setEventLazyClock(true); setEventBatchSameTime(false);
#ifdef SCALING_EVENT_FULL_INITIALIZATION
    setEventCompactInitialization(false);
#else
    setEventCompactInitialization(true);
#endif
#ifdef SCALING_EVENT_NO_BOUNDARY_REUSE
    setEventBoundaryReuse(false);
#else
    setEventBoundaryReuse(true);
#endif
    return computeScaling(graph, mateVec, card, edgWght);
  }

  // Preserved pre-paper implementation. Behavior-changing environment flags
  // are intentionally confined to this explicitly experimental entry point.
  bool computeHybridLegacyExperimental(const Graph& graph,
                                       std::vector<Size>* mateVec,
                                       Size* card, Val* edgWght) {
    mAlgorithmIdentity = ALGORITHM_HYBRID_LEGACY_EXPERIMENTAL;
    mDropPerfectExperimental = false; mSkipFinalizationExperimental = false;
    mHybridCanonical = false; mHybridReference = false;
    setFaithfulLoop(true); setBoundedSearch(true); setHybrid(true);
    setIncrShrink(true); setPersistSearch(true);
    setNoRescanGrow(!std::getenv("MBX_RESCAN_GROW"));
    setNoRescanDual(mNoRescanGrow && std::getenv("MBX_NO_RESCAN_DUAL") &&
                    !std::getenv("MBX_RESCAN_DUAL"));
    setFaithfulReweight(std::getenv("MBX_REWEIGHT") != 0);
    setEventSearch(std::getenv("MBX_EVENTSEARCH") != 0);
    mCarryM = std::getenv("MBX_CARRYM") != 0;
    mContinueAug = std::getenv("MBX_CONTINUEAUG") != 0;
    mIncrRelabel = std::getenv("MBX_INCRRELABEL") != 0;
    setIncrDissolve(std::getenv("MBX_INCRDISSOLVE") != 0);
    setMultiLevel(std::getenv("MBX_MULTILEVEL") != 0);
    setBridgeFvr(std::getenv("MBX_BRIDGEFVR") != 0);
    setHybNoTranslate(false);
    setDynamicShell(std::getenv("MBX_DYNAMIC_SHELL") != 0);
    return computeScaling(graph, mateVec, card, edgWght);
  }

  // D FAST PATH (the cliff fix, MEASURED): for weights that fit (any realistic N), the S5 no-relabel
  // sqrt(n) engine solves the full max-weight PERFECT matching in ONE solve -- exact AND far faster
  // than the scaling loop (n=51200: 17.4s->4.2s; n=204800: intractable->55s; smooth ~2.5x/doubling).
  // The scaling loop's per-scale re-search is the cliff; a single S5 solve sidesteps it entirely.
  // (Scaling is only needed for HUGE weights where one solve's dual magnitude overflows -- not
  // reached at realistic sizes: dual magnitude ~ maxW*n fits long long for maxW*n < 9e18.)
  // Failure is reported directly; this named algorithm never switches to a
  // different solver behind the caller's back.
  bool computeS5Experimental(const Graph& graph, std::vector<Size>* mateVec,
                             Size* card, Val* edgWght) {
    mAlgorithmIdentity = ALGORITHM_S5_EXPERIMENTAL;
    init(graph);                                   // mN, mAdj, mWhat (rounded int weights)
    ScalingS5Engine s5;
    s5.setUseEvent(true); s5.setLazyEvent(true);
    s5.setNoRelabelShrink(true); s5.setNoRelabelAug(true);
    s5.init(mN);
    for (Size u = 0; u < mN; ++u)
      for (Size k = 0; k < mAdj[u].size(); ++k) { Size w = mAdj[u][k]; if (u < w) s5.addEdge(u, w, mWhat[u][k]); }
    std::vector<ScalingS5Engine::Sz> sm;
    LL g = s5.solve(&sm);
    if (g < 0) {
      if (mateVec) mateVec->assign(mN, scNull());
      if (card) *card = 0;
      if (edgWght) *edgWght = 0;
      mNoPerfect = true;
      mLastStatus = STATUS_NO_PERFECT_MATCHING;
      return false;
    }
    if (mateVec) { mateVec->assign(mN, scNull()); for (Size v = 0; v < mN; ++v) (*mateVec)[v] = (sm[v]==(ScalingS5Engine::Sz)-1)?scNull():(Size)sm[v]; }
    Size c = 0; for (Size v = 0; v < mN; ++v) if (sm[v] != (ScalingS5Engine::Sz)-1) ++c; c /= 2;
    if (card) *card = c;
    if (edgWght) *edgWght = (Val)g;
    return c * 2 == mN;
  }

  // Deprecated source-compatible alias: this is a one-shot S5 solver, not
  // Figure 6 Hybrid.
  bool computeHybridFast(const Graph& graph, std::vector<Size>* mateVec,
                         Size* card, Val* edgWght) {
    return computeS5Experimental(graph, mateVec, card, edgWght);
  }

  // ===========================================================================================
  // D PATH (a) OPTION 1: PERSISTENT-S5. One ScalingS5Engine owns (y,z,M,Omega) across all scales;
  // this loop is a thin driver. D2 SKELETON: build ONE S5 from the input graph + solve (single
  // full-weight solve = the scale-1/computePerfect case), prove the driver shell + S5-owns-state.
  // D3+ will wrap the scale loop (scaling/liquidate/§4.2/FVR/perfection) around the SAME S5 state.
  // Gated: only runs when explicitly called (computeHybrid UNCHANGED = the reference).
  // ===========================================================================================
  bool computePersistS5(const Graph& graph, std::vector<Size>* mateVec, Size* card, Val* edgWght) {
    typedef ScalingS5Engine::Sz S5z; const S5z S5NIL = (S5z)-1;
    init(graph);                                  // populates mN, mAdj, mWhat, mWbar, mNumScales
    Size mN0 = mN;
    if (mN0 == 0) { if (mateVec) mateVec->clear(); if (card) *card = 0; if (edgWght) *edgWght = 0; return true; }
    // D3-real: persistent-S5 Liquidationist loop WITH dummies/perfection/finalization (the cross-
    // scale re-opt engine). ONE S5 holds y/z/M/Omega across scales. Vertex space = [0,2mN0): reals
    // [0,mN0) + dummy slots [mN0,2mN0); blossoms >= 2mN0. setReals(mN0) so only reals are rooted.
    // Per scale: weights=top-i bits, y<-2y+3 (reals AND dummies), empty M, FVR (bounded -> leave
    // free reals), Perfection (park each free real on its dummy, y(d)=tau). After all scales:
    // Finalization (un-park, strip dummy edges, re-search reals to perfection on TRUE weights).
    Size V = 2 * mN0;                             // S5 vertex universe (reals + dummy slots)
    ScalingS5Engine s5;
    s5.setUseEvent(true); s5.setLazyEvent(true);
    s5.setNoRelabelShrink(true); s5.setNoRelabelAug(true);
#ifndef PS5_NO_NRD_DISSOLVE
    // THE MISSING SETTER (2026-07-07): the loop enabled no-relabel SHRINK + AUGMENT but NEVER
    // no-relabel DISSOLUTION -> every dissolution fell to the O(n) relabel+full-reseed path (line
    // ~2921), firing 920x @204800 = ~9.2M of the 11.3M scans = the giant-blossom cliff. The whole
    // in-place dissolution machinery (relinkOddCycleS5 + nested-intact, default-flag-enabled) was
    // built but unused because THIS call was absent. Opt out: -DPS5_NO_NRD_DISSOLVE.
    s5.setNoRelabelDissolve(true);
#endif
    s5.init(V); s5.setReals(mN0); s5.setIgnoreDummies(true);
    s5.setAllowOvertight(true);   // carried duals are Property-2 (slack can be <0) -- over-tight eligible
    // FAITHFUL: dummies NEVER searchable.
    // The loop's mFaithfulLoop branch searches REAL edges only, parks free reals on dummies AFTER
    // (lazily), strips them next scale. Making dummies searchable = the "premature-park" bug. So we
    // add NO dummy edges; dummies are pure y-carry slots [mN0,2mN0). (setIgnoreDummies = belt+braces.)
    for (Size u = 0; u < mN0; ++u)
      for (Size k = 0; k < mAdj[u].size(); ++k) { Size w = mAdj[u][k]; if (u < w) s5.addEdge(u, w, 0); }
    std::vector<std::vector<ScalingS5Engine::Sz> >& AV = s5.adjV();
    std::vector<std::vector<LL> > s5wbar(mN0);                  // only reals have edges/wbar
    for (Size u = 0; u < mN0; ++u) { s5wbar[u].assign(AV[u].size(), 0);
      for (Size k = 0; k < AV[u].size(); ++k) { Size w = AV[u][k];
        for (Size j = 0; j < mAdj[u].size(); ++j) if (mAdj[u][j] == w) { s5wbar[u][k] = mWbar[u][j]; break; } } }
    Size nbits = mNumScales;
    std::vector<std::vector<LL> > wPrev(mN0);
    for (Size u = 0; u < mN0; ++u) wPrev[u].assign(AV[u].size(), 0);
    std::vector<LL>& Y = s5.y();
    Y.assign(V, 0);
#ifdef PS5MEAS
    // Instrumentation: count REAL-REAL matched pairs at phase boundaries to measure
    // augments-per-scale (Phase A) and Finalization augments. A "Phase A augment" = a
    // newly-formed real-real matched edge during searchBounded/searchToPerfect.
    auto countRR = [&](void)->Size { Size c=0; std::vector<ScalingS5Engine::Sz>& MM=s5.mate();
      for (Size v=0; v<mN0; ++v) if (MM[v]!=S5NIL && (Size)MM[v]<mN0 && (Size)MM[v]>v) ++c; return c; };
    std::fprintf(stderr,"[PS5MEAS] n=%lu scales=%lu tau=%lu\n",
      (unsigned long)mN0,(unsigned long)nbits,(unsigned long)mTau);
#endif
    for (Size scale = 1; scale <= nbits; ++scale) {
      Size bitpos = nbits - scale;
      std::vector<std::vector<LL> >& EW = s5.edgeW();
      for (Size u = 0; u < mN0; ++u)
        for (Size k = 0; k < AV[u].size(); ++k) {
          LL bit = (s5wbar[u][k] >> bitpos) & 1;
          LL w2 = 2 * (wPrev[u][k] + bit);     // w2 = doubled scale weight (= loop's mW)
          EW[u][k] = w2 / 2;                   // S5 stores w (slack uses 2*w)
          wPrev[u][k] = w2;                    // BUGFIX: w'(this scale) = FULL w2 (= loop's mW), NOT w2/2
        }
      for (Size v = 0; v < V; ++v) Y[v] = 2 * Y[v] + 3;   // rescale ALL y (reals + dummy slots)
#ifndef PS5_NO_LIQ_FOLD
      // DEFAULT-ON (2026-07-05, 594-gated). Spec step 2 also doubles z (z <- 2z'). With z doubled
      // here, the liquidation fold inside resetForestKeepDuals is uniformly y += z/2 -- correct at
      // scale boundaries (= z_s, the spec's Liquidate on doubled z). Opt out: -DPS5_NO_LIQ_FOLD
      // (must match the engine-side fold flag).
      { std::vector<LL>& Z5 = s5.z(); for (Size b = 0; b < Z5.size(); ++b) Z5[b] *= 2; }
#endif
#ifdef PS5_DOMINATE_S1
      // EXPERIMENT: on scale 1 (fresh duals, NO carry) the y<-2*0+3=3 seed leaves EVERY edge massively
      // OVER-TIGHT (weight ~2^MSB >> y=3) -> overtight mode schedules ALL m edges to fire at tNow=0 ->
      // O(m^2)-ish single-level blowup (measured: scale1 = 125s @ n=25600 vs solve()=2.85s). Seed scale
      // 1 with proper DOMINATION (y = max scale weight) so slacks are >=0 and edges fire at staggered
      // times = the clean solve() regime. Scales 2+ keep the carried 2y'+3.
      if (scale == 1) { LL mw=0; for (Size u=0;u<mN0;++u) for (Size k=0;k<EW[u].size();++k) if (EW[u][k]>mw) mw=EW[u][k];
        for (Size v=0; v<mN0; ++v) Y[v]=mw; }
#endif
#ifdef PS5_CARRY_M
      // LEGACY M-carry (pre-2026-07-05 default). Fig 5 step 1 is M<-EMPTY per scale; under a phase
      // budget the carry FREEZES the loop (scale 1 completes M, later scales have no free vertices,
      // duals never differentiate -> 578/594 wrong). Kept opt-in for the drive-to-PM A/B only.
      std::vector<ScalingS5Engine::Sz> savedM = s5.mate();
#endif
#ifdef PS5MEAS
      std::fprintf(stderr,"[PS5MEAS] pre-reset scale %lu: numNodes=%lu (mN=%lu)\n",
        (unsigned long)scale,(unsigned long)s5.numNodes(),(unsigned long)V);
#endif
#ifdef PS5_STEP6
      // Fig 5 step 6 (LEVELING) needs the OLD scale's maximal small blossoms' member lists --
      // capture them BEFORE the reset wipes the blossom structure. Spec: small = |B'| < tau.
      std::vector<std::vector<ScalingS5Engine::Sz> > s6Blossoms;
      s5.collectMaximalSmall((ScalingS5Engine::Sz)(mTau > 0 ? mTau : 1), &s6Blossoms);
#endif
#ifdef PS5TIME
      { struct timespec _t0,_t1; clock_gettime(CLOCK_MONOTONIC,&_t0);
        s5.resetForestKeepDuals();
        clock_gettime(CLOCK_MONOTONIC,&_t1);
        std::fprintf(stderr,"[PS5TIME] scale %lu reset=%.3fs\n",(unsigned long)scale,
          (double)((_t1.tv_sec-_t0.tv_sec)+(_t1.tv_nsec-_t0.tv_nsec)/1e9)); }
#else
      s5.resetForestKeepDuals();                          // fresh forest, KEEP duals; restore M' next
#endif
#ifdef PS5_CARRY_M
      // (legacy restore -- see the PS5_CARRY_M note above; default is the paper's M<-EMPTY)
      { std::vector<ScalingS5Engine::Sz>& MM=s5.mate(); for(Size v=0;v<mN0;++v) if(v<savedM.size()&&savedM[v]!=S5NIL&&(Size)savedM[v]<mN0){ MM[v]=savedM[v]; } }
#endif
#ifdef PS5_STEP6
      // ---- Fig 5 STEP 6 (LEVELING), spec L75-81 -- verbatim semantics: -------------------------
      //   for each MAXIMAL old small blossom B':
      //     while max{ y(u) | u in B' free } > 0:
      //       Y  = that max;  F = { free u in B' | y(u)=Y };  Y' = max{0, max y over B'\F free}
      //       PQSearch(F) (Criterion 1) until an augmenting LEVEL completes or (Y-Y') dual
      //       adjustments happen.  (One multi-augment level approximates the spec's "an
      //       augmenting path is found+augmented" -- MV batching within one dual level.)
      // The z-fold at liquidation raised members' y ABOVE the surroundings; leveling burns the
      // excess down (matching members internally as paths appear) so the step-7 tau-bounded FVR
      // starts from leveled duals -- the paper's precondition for the tau-tick accounting.
      { std::vector<LL>& Y6 = s5.y();
        std::vector<ScalingS5Engine::Sz>& M6 = s5.mate();
        s5.setCrit(ScalingS5Engine::CRIT1);
        for (Size bi = 0; bi < s6Blossoms.size(); ++bi) {
          std::vector<ScalingS5Engine::Sz>& L = s6Blossoms[bi];
          Size s6Guard = 4 * L.size() + 64;              // each round augments or levels Y down
          while (s6Guard--) {
            LL Ymax = 0, Y2 = 0;                          // Y and Y' (Y' floored at 0 per spec)
            for (Size j = 0; j < L.size(); ++j) { Size v = L[j];
              if (v < mN0 && M6[v] == S5NIL && Y6[v] > Ymax) Ymax = Y6[v]; }
            if (Ymax <= 0) break;
            std::vector<ScalingS5Engine::Sz> F6;
            for (Size j = 0; j < L.size(); ++j) { Size v = L[j];
              if (v >= mN0 || M6[v] != S5NIL) continue;
              if (Y6[v] == Ymax) F6.push_back(v);
              else if (Y6[v] > Y2) Y2 = Y6[v]; }
            if (F6.empty()) break;
            s5.setRootSet(F6);
            s5.setTickBudget(Ymax - Y2 > 0 ? Ymax - Y2 : 1);
            s5.setPhaseBudget(1);                        // stop after the first augmenting level
            s5.searchBounded(1);
            s5.setPhaseBudget(0);
            s5.setTickBudget(0);
            s5.clearRootSet();
          }
        }
      }
#endif
      // Phase A: BOUNDED real-edge search under the carried duals (Crit1). tau=mTau (~sqrt n,
      // Liquidationist). PAPER-FAITHFUL (Fig 5 step 7): do NOT drive each scale to a real PM --
      // run at most tau SearchOne's and leave still-free reals to be PARKED on dummies by Phase B,
      // so the next scale has free reals to reoptimize incrementally. This spreads the n/2 reopt
      // augments across the log(nN) scales (each O(sqrt n)) instead of dumping them all on
      // Finalization (the "mistake to avoid"). PS5_DRIVE_TO_PM restores the old decorative
      // drive-to-PM for A/B measurement.
#ifdef PS5_CRIT2A
      s5.setCrit(ScalingS5Engine::CRIT2);   // spec Fig 5 step 7: SearchOne uses Criterion 2
#else
      s5.setCrit(ScalingS5Engine::CRIT1);
#endif
#ifdef PS5MEAS
      Size rrBefore = countRR();
#endif
// DEFAULT (2026-07-05): constant phase budget 32. The full sweep on the fixed stack (EMPTY_M +
// LIQ_FOLD + FIN_KEEPB) is exact at every budget (B=4/8/16/32 all 594/594) and B=32 drains the
// free vertices entirely IN-scale (Finalization = 0 augs), making the loop the FASTEST engine
// measured: 1.38/3.19/10.03s at n=25600/51200/102400 vs Hybrid-loop 4.14/8.82/174.96s (the
// Hybrid 102400 cliff!) and single-shot 1.82/4.39/11.34s. Budget=tau (~sqrt n, the spec letter)
// is ALSO exact but re-solves each scale to completion = 6-15x slower; use -DPS5_BUDGET_TAU for
// the spec-letter configuration, -DPS5_DRIVE_TO_PM for the legacy unbudgeted mode.
#if !defined(PS5_PHASE_BUDGET) && !defined(PS5_TICK_BUDGET) && !defined(PS5_BUDGET_TAU) && !defined(PS5_DRIVE_TO_PM)
#define PS5_PHASE_BUDGET 32
#endif
#if defined(PS5_TICK_BUDGET)
      // PAPER-FAITHFUL tau bound (Fig 5 step 7): bound the scale's TOTAL dual adjustment to
      // PS5_TICK_BUDGET * tau clock ticks (multiplier 1 = the spec's tau SearchOne calls = tau unit
      // adjustments). Unlike the phase budget, this bounds exactly what the paper bounds; the search
      // stops mid-phase when the budget is exhausted, leaving still-free reals for Phase B parking.
      s5.setTickBudget((LL)mTau * (LL)(PS5_TICK_BUDGET));
      s5.searchToPerfect();               // stops early on budget exhaustion
      s5.setTickBudget(0);
#elif defined(PS5_BUDGET_TAU)
      // PAPER-FAITHFUL per-scale allowance: budget = tau (~sqrt n, the Fig 5 step-7 allowance),
      // as a PHASE budget. Scales with n where a compile-time constant cannot; at n=25600 the
      // measured sweep (B=2..32, all 594-exact post-fix) shows fin cost ~ freeLeft * O(m) --
      // one aug/phase with a full reseed each -- so draining free vertices in-scale is the win.
      s5.setPhaseBudget((ScalingS5Engine::Sz)(mTau > 0 ? mTau : 1));
#ifdef PS5TIME
      { struct timespec _t0,_t1; clock_gettime(CLOCK_MONOTONIC,&_t0);
        s5.searchBounded(1);
        clock_gettime(CLOCK_MONOTONIC,&_t1);
        double _dt=(_t1.tv_sec-_t0.tv_sec)+(_t1.tv_nsec-_t0.tv_nsec)/1e9;
        Size _fr=0; { std::vector<ScalingS5Engine::Sz>& MM=s5.mate(); for(Size v=0;v<mN0;++v) if(MM[v]==S5NIL)++_fr; }
        std::fprintf(stderr,"[PS5TIME] scale %lu phaseA=%.3fs freeAfter=%lu\n",(unsigned long)scale,_dt,(unsigned long)_fr); }
#else
      s5.searchBounded(1);
#endif
      s5.setPhaseBudget(0);
#elif defined(PS5_PHASE_BUDGET)
      // THE DEFAULT PATH (PS5_PHASE_BUDGET=32 unless overridden; see the default block above).
      // Bound Phase A to PS5_PHASE_BUDGET MV-phases per scale. Bounding searchEvent CALLS (the old
      // searchBounded(tau)) was inert -- one searchEvent drains the whole matching; the PHASE is
      // the real lever. The old finding "small budgets break weight-optimality at Finalization"
      // was an artifact of the pre-fix loop (M-carry + missing fold + fin-reset); on the fixed
      // stack every budget tested (2..32) is 594-exact.
      s5.setPhaseBudget((ScalingS5Engine::Sz)(PS5_PHASE_BUDGET));
#ifdef PS5TIME
      { struct timespec _t0,_t1; clock_gettime(CLOCK_MONOTONIC,&_t0);
        s5.searchBounded(1);
        clock_gettime(CLOCK_MONOTONIC,&_t1);
        double _dt=(_t1.tv_sec-_t0.tv_sec)+(_t1.tv_nsec-_t0.tv_nsec)/1e9;
        Size _fr=0; { std::vector<ScalingS5Engine::Sz>& MM=s5.mate(); for(Size v=0;v<mN0;++v) if(MM[v]==S5NIL)++_fr; }
        std::fprintf(stderr,"[PS5TIME] scale %lu phaseA=%.3fs freeAfter=%lu\n",(unsigned long)scale,_dt,(unsigned long)_fr); }
#else
      s5.searchBounded(1);                 // ONE phase-budgeted searchEvent = up to budget MV-phases
#endif
      s5.setPhaseBudget(0);                // restore unbounded for any later/Finalization search
#else
#ifdef PS5TIME
      { struct timespec _t0,_t1; clock_gettime(CLOCK_MONOTONIC,&_t0);
        s5.searchToPerfect();
        clock_gettime(CLOCK_MONOTONIC,&_t1);
        double _dt=(_t1.tv_sec-_t0.tv_sec)+(_t1.tv_nsec-_t0.tv_nsec)/1e9;
        Size _fr=0; { std::vector<ScalingS5Engine::Sz>& MM=s5.mate(); for(Size v=0;v<mN0;++v) if(MM[v]==S5NIL)++_fr; }
        std::fprintf(stderr,"[PS5TIME] scale %lu searchToPerfect=%.3fs freeAfter=%lu\n",(unsigned long)scale,_dt,(unsigned long)_fr); }
#else
      s5.searchToPerfect();                // DEFAULT: drive each scale to a real PM (proven, exact 594)
#endif
#endif
#ifdef PS5MEAS
      { Size rrAfter = countRR(); Size freeReals=0; std::vector<ScalingS5Engine::Sz>& MM=s5.mate();
        for (Size v=0; v<mN0; ++v) if (MM[v]==S5NIL) ++freeReals;
        std::fprintf(stderr,"[PS5MEAS] scale %lu: phaseA_augs=%ld rr=%lu freeReals=%lu phases=%lu\n",
          (unsigned long)scale,(long)((long)rrAfter-(long)rrBefore),
          (unsigned long)rrAfter,(unsigned long)freeReals,
          (unsigned long)s5.lastPhaseCount()); }
#endif
#ifdef SCMAG
      if (scale <= 4) { LL ys=0; for (Size v=0; v<mN0; ++v) ys+=Y[v];
        std::fprintf(stderr,"[PS5 scale %lu post-search] sumY=%lld y:",(unsigned long)scale,(long long)ys);
        for (Size v=0;v<mN0;++v) std::fprintf(stderr," %lld",(long long)Y[v]);
        std::fprintf(stderr,"\n  mate:"); { std::vector<ScalingS5Engine::Sz>& MM=s5.mate(); for(Size v=0;v<mN0;++v) std::fprintf(stderr," %ld",(long)(MM[v]==S5NIL?-1:(long)MM[v])); } std::fprintf(stderr,"\n"); }
#endif
      // Phase B (Perfection): park each STILL-free real on its dummy: M[u]=d, y(d)=tau. No edge.
      std::vector<ScalingS5Engine::Sz>& M = s5.mate();
      for (Size u = 0; u < mN0; ++u) if (M[u] == S5NIL) { Size d = mN0 + u; Y[d] = (LL)mTau; M[u] = d; M[d] = u; }
#ifdef PS5_SCALECHK
      // Property-2 invariant check after each scale's Phase A+B: near-domination (slack >= -2 on
      // every real edge) and near-tightness (matched slack in [-2,0]). First scale to violate = the
      // step that breaks the invariant. Slack convention: y(u)+y(v) - 2*mAdjW (blossom z omitted;
      // live-blossom count printed for context).
      { std::vector<std::vector<LL> >& EWc = s5.edgeW();
        LL minS = 0, maxS = 0; long violDom = 0, violM = 0; bool first = true;
        for (Size u = 0; u < mN0; ++u) for (Size k = 0; k < AV[u].size(); ++k) { Size v2 = AV[u][k]; if (v2 <= u || v2 >= mN0) continue;
          LL sl = Y[u] + Y[v2] - 2 * EWc[u][k];
          if (first) { minS = maxS = sl; first = false; }
          if (sl < minS) minS = sl; if (sl > maxS) maxS = sl;
          if (sl < -2) ++violDom;
          if ((Size)M[u] == v2 && (sl < -2 || sl > 0)) ++violM;
        }
        std::fprintf(stderr,"[SCALECHK] scale %lu: slack[%lld,%lld] domViol=%ld matchedViol=%ld blossoms=%ld\n",
          (unsigned long)scale,(long long)minS,(long long)maxS,violDom,violM,
          (long)((long)s5.numNodes()-(long)V));
      }
#endif
    }
    // Finalization: free every real parked on a dummy, then Crit3 re-search reals to perfection.
    { std::vector<ScalingS5Engine::Sz>& M = s5.mate();
#if defined(PS5_DUMMY_SEARCH)
      // PAPER-FAITHFUL Finalization (spec L85-86,109-112): DELETE dummies -> parked reals become FREE ->
      // Crit3 search re-matches them, augmenting paths routing real->dummy->real ROTATE M to optimum.
      // STEP A: add SEARCHABLE 0-weight dummy edges NOW (only at Finalization, so no premature-park
      // during the scale loop), drop ignore-dummies, KEEP the real-real matches, free parked reals,
      // and let Crit3 searchToPerfect rotate through the dummies. Tests the dummy-search plumbing.
      for (Size u = 0; u < mN0; ++u) { Size d = mN0 + u; s5.addEdge(u, d, 0); }   // 0-weight, searchable
      s5.setIgnoreDummies(false);
      for (Size u = 0; u < mN0; ++u) if (M[u] != S5NIL && (Size)M[u] >= mN0) { M[(Size)M[u]] = S5NIL; M[u] = S5NIL; }
      { std::vector<ScalingS5Engine::Sz> finSaved = s5.mate();
        s5.resetForestKeepDuals();
        std::vector<ScalingS5Engine::Sz>& MM=s5.mate();
        for(Size v=0;v<mN0;++v) if(v<finSaved.size()&&finSaved[v]!=S5NIL&&(Size)finSaved[v]<mN0){ MM[v]=finSaved[v]; } }
      s5.setCrit(ScalingS5Engine::CRIT3);
#elif !defined(PS5_FIN_CARRY) && !defined(PS5_FIN_RESEED) && !defined(PS5_FIN_RESET)
      // DEFAULT (2026-07-05, was PS5_FIN_KEEPB; opt out with -DPS5_FIN_RESET for the old
      // reset-at-finalization behavior).
      // PAPER-FAITHFUL Finalization (spec L85-86: delete dummies, PQSearch Crit3 -- Omega is KEPT).
      // searchEvent() begins with relabel(), which rebuilds the search forest from the current
      // mMate + blossom structure, so NO reset is needed here at all: unpark the parked reals,
      // keep live blossoms (their z intact -- liquidating them here made matched cross-boundary
      // edges over-tight and Crit3-ineligible, the n=8 rep=1 off-by-14), keep the real-real
      // matching, and let the Crit3 search finish the freed vertices.
      for (Size u = 0; u < mN0; ++u) if (M[u] != S5NIL && (Size)M[u] >= mN0) { M[(Size)M[u]] = S5NIL; M[u] = S5NIL; }
#ifdef PS5_PARITY_ALIGN
      // Paper precondition (§2/ar5iv verbatim): "the y-values of all F vertices have the same
      // parity". Scale starts guarantee it (2y+3 = odd), but by Finalization the freed reals carry
      // mixed parities (each spent different tick-counts in trees), producing ODD slacks -> the
      // event queue's slack/2 truncation misfires by one -> the off-by-2 weight class. Align: drop
      // the off-parity free duals by 1 (lowering y keeps near-domination >= -2 slack-neutral-or-
      // -looser on incident edges and cannot break matched near-tightness <= 0).
      { Size odd = 0, even = 0;
        for (Size u = 0; u < mN0; ++u) if (M[u] == S5NIL) { if (Y[u] & 1) ++odd; else ++even; }
        LL wantOdd = (odd >= even) ? 1 : 0;
        for (Size u = 0; u < mN0; ++u) if (M[u] == S5NIL && ((Y[u] & 1) != wantOdd)) Y[u] -= 1;
      }
#endif
      s5.setCrit(ScalingS5Engine::CRIT3);
#elif defined(PS5_FIN_CARRY)
      // CARRY the real-real matches into Finalization (don't re-solve from empty): resetForestKeepDuals
      // wipes M, so snapshot + restore the real-real part, leaving only the freed parked reals unmatched.
      for (Size u = 0; u < mN0; ++u) if (M[u] != S5NIL && (Size)M[u] >= mN0) { M[(Size)M[u]] = S5NIL; M[u] = S5NIL; }
      { std::vector<ScalingS5Engine::Sz> finSaved = s5.mate();
        s5.resetForestKeepDuals();
        std::vector<ScalingS5Engine::Sz>& MM=s5.mate();
        for(Size v=0;v<mN0;++v) if(v<finSaved.size()&&finSaved[v]!=S5NIL&&(Size)finSaved[v]<mN0){ MM[v]=finSaved[v]; } }
      s5.setCrit(ScalingS5Engine::CRIT3);
#elif defined(PS5_FIN_RESEED)
      // Re-seed FULL DOMINATION (y=maxw on reals) before the Crit3 re-solve from empty -- the proven
      // 264/264 "fresh dual" policy. Tests whether the reset-resolve off-by-3 is carried-dual staleness.
      for (Size u = 0; u < mN0; ++u) if (M[u] != S5NIL && (Size)M[u] >= mN0) { M[(Size)M[u]] = S5NIL; M[u] = S5NIL; }
      { s5.resetForestKeepDuals();
        std::vector<LL>& YY = s5.y(); std::vector<std::vector<LL> >& EWf = s5.edgeW();
        LL mw = 0; for (Size u=0;u<mN0;++u) for (Size k=0;k<EWf[u].size();++k) if (EWf[u][k]>mw) mw=EWf[u][k];
        for (Size v=0; v<mN0; ++v) YY[v] = mw;   // reals dominated; dummy slots [mN0,2mN0) stay 0
        s5.setCrit(ScalingS5Engine::CRIT3); }
#else
      for (Size u = 0; u < mN0; ++u) if (M[u] != S5NIL && (Size)M[u] >= mN0) { M[(Size)M[u]] = S5NIL; M[u] = S5NIL; }
      s5.resetForestKeepDuals(); s5.setCrit(ScalingS5Engine::CRIT3);
#endif
#ifdef PS5MEAS
      Size rrBeforeFin = countRR();
#endif
#ifdef PS5TIME
      { struct timespec _t0,_t1; clock_gettime(CLOCK_MONOTONIC,&_t0);
        s5.searchToPerfect();
        clock_gettime(CLOCK_MONOTONIC,&_t1);
        std::fprintf(stderr,"[PS5TIME] FINALIZATION searchToPerfect=%.3fs phases=%lu\n",
          (double)((_t1.tv_sec-_t0.tv_sec)+(_t1.tv_nsec-_t0.tv_nsec)/1e9),
          (unsigned long)s5.lastPhaseCount()); }
#else
      s5.searchToPerfect();
#endif
#ifdef PS5MEAS
      std::fprintf(stderr,"[PS5MEAS] FINALIZATION: rrBefore=%lu rrAfter=%lu fin_augs=%ld (n/2=%lu)\n",
        (unsigned long)rrBeforeFin,(unsigned long)countRR(),
        (long)((long)countRR()-(long)rrBeforeFin),(unsigned long)(mN0/2));
#endif
#ifdef PS5_QSTAT
      s5.qstatPrint();   // whole-loop event-queue totals (inserts/pops/acted + distance histograms)
#endif
#ifdef S5DBG
      { unsigned long mm=0; for (Size u=0;u<mN0;++u) mm+=AV[u].size(); s5.dbgDumpCounters(mm/2); }
#endif
#ifdef PS5_STAGE0
      s5.s0Print();      // persistent-clock Stage 0: redundant-regrowth fraction
#endif
#ifdef PS5_FINDUMP
      // Post-finalization dual/slack dump (small n): slack(e) = y(u)+y(v) - 2*mAdjW (S5 convention;
      // blossom z terms omitted -- flagged if any live blossom exists). Crit3 can only use edges
      // with slack in {0,-2}; a matched edge outside that set blocks rerouting (parity suspect).
      { std::vector<LL>& YY = s5.y(); std::vector<std::vector<LL> >& EWd = s5.edgeW();
        std::vector<ScalingS5Engine::Sz>& MM = s5.mate();
        std::fprintf(stderr,"[FINDUMP] numNodes=%lu (mN=%lu; blossoms live=%ld)\n",
          (unsigned long)s5.numNodes(),(unsigned long)V,(long)((long)s5.numNodes()-(long)V));
        std::fprintf(stderr,"[FINDUMP] y:"); for (Size v=0; v<mN0; ++v) std::fprintf(stderr," %lld",(long long)YY[v]);
        std::fprintf(stderr,"\n[FINDUMP] edges (u,v,2w,slack,M?):\n");
        for (Size u=0; u<mN0; ++u) for (Size k=0; k<AV[u].size(); ++k) { Size w2v=AV[u][k]; if (w2v<=u||w2v>=mN0) continue;
          LL sl = YY[u]+YY[w2v]-2*EWd[u][k];
          std::fprintf(stderr,"  (%lu,%lu) 2w=%lld slack=%lld %s%s\n",(unsigned long)u,(unsigned long)w2v,
            (long long)(2*EWd[u][k]),(long long)sl,((Size)MM[u]==w2v)?"M":"-",
            (sl!=0&&sl!=-2)?"  <-- Crit3-INELIGIBLE":"");
        }
      }
#endif
    }
    std::vector<ScalingS5Engine::Sz>& M = s5.mate();
    Size c = 0; LL tw = 0;
    for (Size v = 0; v < mN0; ++v) if (M[v] != S5NIL && (Size)M[v] > v && (Size)M[v] < mN0) {
      ++c; Size w = (Size)M[v];
      for (Size j = 0; j < mAdj[v].size(); ++j) if (mAdj[v][j] == w) { tw += mWhat[v][j]; break; }
    }
    if (mateVec) { mateVec->assign(mN0, scNull()); for (Size v = 0; v < mN0; ++v) (*mateVec)[v] = (M[v]==S5NIL||(Size)M[v]>=mN0)?scNull():(Size)M[v]; }
    if (card) *card = c;
    if (edgWght) *edgWght = (Val)tw;
    return c * 2 == mN0;
  }

  // expose for Gate-1b: run the WITH-blossom max-weight PERFECT search.
  // Returns true on success; *noPerfect set if no perfect matching exists.
  bool computePerfect(const Graph& graph,
                      std::vector<Size>* mateVec,
                      Size* card, Val* edgWght,
                      bool* noPerfect);
private:
};

// ---------------------------------------------------------------------------
// Implementation (header-only)
// ---------------------------------------------------------------------------

inline Size ScalingMatching::ceilLog2(LL x) {
  if (x <= 1) return 1;            // at least one scale
  Size bits = 0;
  // ceil(log2(x)) = bit_length(x-1). This form never left-shifts into
  // LL's sign bit when x is close to the supported 64-bit limit.
  unsigned long long v = (unsigned long long)(x - 1);
  while (v != 0) { v >>= 1; ++bits; }
  return bits;
}

inline Size ScalingMatching::bitLength(LL x) {
  if (x <= 0) return 1;            // retain one zero scale for all-zero input
  Size bits = 0;
  unsigned long long v = (unsigned long long)x;
  while (v != 0) { v >>= 1; ++bits; }
  return bits;
}

inline bool ScalingMatching::checkedAddLL(LL a, LL b, LL* out) {
  const LL hi = std::numeric_limits<LL>::max();
  const LL lo = std::numeric_limits<LL>::min();
  if ((b > 0 && a > hi - b) || (b < 0 && a < lo - b)) return false;
  *out = a + b; return true;
}

inline bool ScalingMatching::checkedMulLL(LL a, LL b, LL* out) {
  if (a == 0 || b == 0) { *out = 0; return true; }
  const LL hi = std::numeric_limits<LL>::max();
  const LL lo = std::numeric_limits<LL>::min();
  if (a == -1) { if (b == lo) return false; *out = -b; return true; }
  if (b == -1) { if (a == lo) return false; *out = -a; return true; }
  if (a > 0) {
    if (b > 0) { if (a > hi / b) return false; }
    else { if (b < lo / a) return false; }
  } else {
    if (b > 0) { if (a < lo / b) return false; }
    else { if (a < hi / b) return false; }
  }
  *out = a * b; return true;
}

inline void ScalingMatching::init(const Graph& graph) {
  mCounters = OperationCounters();
  mEventStructuresActive = false;
  mEventSfmInitialized = false;
  mEventSfmScale = scNull();
  mSfmSplitPivots.clear();
#ifdef EVENT_S5_TARGET_KEY_SHADOW
  mS5TargetKeyActive = false;
  mS5TargetKeySplitPivots.clear();
  mS5TargetKeySource.clear();
  mS5TargetKeyOwner.clear();
  mS5TargetKeyTargetOwner.clear();
  mS5TargetKeySourceToken.clear();
  mS5TargetKeyEdgeIndex.clear();
  mS5TargetKeyMatched.clear();
  mS5TargetKeyOwnerTargets.clear();
  mS5TargetKeyTargetOwnerTargets.clear();
  mS5TargetKeyTreeNodes.clear();
  mS5TargetKeyTreeVertices.clear();
  mS5TargetKeyTreeAllVertices.clear();
  mS5TargetKeyTreeOwnedTargets.clear();
  mS5TargetKeyTrackedRoot.clear();
  mS5TargetKeyTrackedVertexToken.clear();
  mS5TargetKeyAllVertexToken.clear();
  mS5TargetKeyTrackedRootStamp.clear();
  mS5TargetKeyTrackedVertexStamp.clear();
  mS5TargetKeyAllVertexStamp.clear();
  mS5TargetKeyRevealStamp.clear();
  mS5TargetKeyRevealClock = 0;
  mS5TargetKeyClockBase = 0;
  mS5TargetKeyTouchedTargets.clear();
  mS5TargetKeyTouchedOwners.clear();
  mS5TargetKeyTouchedTokens.clear();
  mS5TargetKeyTargetTouchStamp.clear();
  mS5TargetKeyOwnerTouchStamp.clear();
  mS5TargetKeyTokenTouchStamp.clear();
  mS5TargetKeyTouchClock = 1;
#ifdef EVENT_S5_SINGLETON_TARGET_FAST
  mS5TargetKeySingleton.clear();
  mS5TargetKeySingletonValue.clear();
  mS5TargetKeySingletonListId.clear();
#endif
  mS5TargetFrontier.clear();
  mS5TargetFrontierTouched.clear();
  mS5TargetFrontierClock = 1;
  mS5TargetKeyTopologyScale = scNull();
  mS5TargetKeyTopologyReady = false;
  mS5TargetKeySeeding = false;
#endif
  mEventTouchedStamp.clear();
  mEventTouchedNodes.clear();
  mEventDormantFlag.clear();
  mEventDormantNodes.clear();
  mEventFreeCandidates.clear();
  mEventCriterion2Free.clear();
  mEventCriterion2FreeScale = scNull();
  mEventCriterion2FreeUniverse = scNull();
  mEventInitialRoots.clear();
  mEventDeferredEdgeEvents.clear();
  mEventBoundaryReuseEvents.clear();
  mEventBoundaryReuseReady = false;
  mEventTouchedClock = 1;
  mShellDissolveStamp.clear();
  mShellDissolveClock = 1;
  mPreviousCriterion2AdjustmentOnly = false;
  Size numEdgs = 0;
  graph.GetSize(&mN, &numEdgs);
  mLastStatus = STATUS_SUCCESS;
  mBActive.clear();
  if (mN > (std::numeric_limits<Size>::max() - 4) / 4) {
    mLastStatus = STATUS_UNSUPPORTED_OR_OVERFLOWING_INPUT;
    return;
  }

  mAdj.assign(mN, std::vector<Size>());
  mWhat.assign(mN, std::vector<LL>());
  mMaxW = 0;
  mM = 0;

  // Read adjacency + weights. The DPS scaling algorithm is specified for
  // non-negative integral weights. Silently rounding here changes the problem
  // instance and makes the benchmark identity ambiguous, so reject values
  // outside that domain.
  for (Size v = 0; v < mN; ++v) {
    Size deg = 0;
    std::vector<Size> rawAdj;
    std::vector<Val> rawWgt;
    graph.GetAdj(v, &deg, &rawAdj, &rawWgt);
    mAdj[v].reserve(rawAdj.size());
    mWhat[v].reserve(rawWgt.size());
    for (Size i = 0; i < rawAdj.size(); ++i) {
      if (rawAdj[i] == v) continue;            // skip self-loops
      Size w = rawAdj[i];
      if (w >= mN || i >= rawWgt.size()) {
        mLastStatus = STATUS_UNSUPPORTED_OR_OVERFLOWING_INPUT;
        return;
      }
      Val wd = rawWgt[i];
      if (!std::isfinite((double)wd) || wd < 0 ||
          std::floor((double)wd) != (double)wd ||
          wd >= (Val)std::numeric_limits<LL>::max()) {
        mLastStatus = STATUS_UNSUPPORTED_OR_OVERFLOWING_INPUT;
        return;
      }
      LL wi = (LL)wd;
      mAdj[v].push_back(w);
      mWhat[v].push_back(wi);
      if (wi > mMaxW) mMaxW = wi;
    }
  }
  // each undirected edge appears in both endpoints' lists -> count once
  {
    Size dir = 0;
    for (Size v = 0; v < mN; ++v) {
      if (dir > std::numeric_limits<Size>::max() - mAdj[v].size()) {
        mLastStatus = STATUS_UNSUPPORTED_OR_OVERFLOWING_INPUT;
        return;
      }
      dir += mAdj[v].size();
    }
    mM = dir / 2;
  }

  // Extended weight wbar(e) = (n/2 + 1) * what(e).
  LL factor = (LL)(mN / 2) + 1;
  mWbar.assign(mN, std::vector<LL>());
  for (Size v = 0; v < mN; ++v) {
    mWbar[v].assign(mAdj[v].size(), 0);
    for (Size i = 0; i < mAdj[v].size(); ++i)
      if (!checkedMulLL(factor, mWhat[v][i], &mWbar[v][i])) {
        mLastStatus = STATUS_UNSUPPORTED_OR_OVERFLOWING_INPUT;
        return;
      }
  }

  // The observed maximum is inclusive, so the scale count is its bit length.
  // ceilLog2(maxWbar) is one too small when maxWbar is an exact power of two.
  LL maxWbar = 0;
  if (!checkedMulLL(factor, mMaxW, &maxWbar)) {
    mLastStatus = STATUS_UNSUPPORTED_OR_OVERFLOWING_INPUT;
    return;
  }
  mNumScales = bitLength(maxWbar);

  // tau ~ sqrt(n) (>= 1).  Liquidationist uses tau=sqrt(n) (Fig 5 step 7: SearchOne tau times).
  // HYBRID (Thm 4.3) uses a LARGER tau in [sqrt(n) log n, n^{2/3}] -- the two-stage FVR (sqrt n
  // SearchOne, then BucketSearch Crit3 up to tau) keeps it within O(m sqrt n). We pick
  // tau = sqrt(n) * ceil(log2 n) (the low end of the range; >= sqrt n so never smaller than
  // Liquidationist). mSqrtN = floor-ish sqrt(n) is the stage-1 count for Hybrid.
  {
    Size t = mN == 0 ? 1 : (Size)std::sqrt((double)mN);
    if (t < 1) t = 1;
    while (t < mN / t || (t == mN / t && mN % t != 0)) ++t;
    while (t > 1) {
      Size q = t - 1;
      if (q > mN / q || (q == mN / q && mN % q == 0)) --t;
      else break;
    }
    mSqrtN = (t < 1) ? 1 : t;
    if (mHybrid) {
      Size lg = ceilLog2((LL)(mN > 1 ? mN : 2));
      Size mult = lg < 1 ? 1 : lg;
      if (mSqrtN > std::numeric_limits<Size>::max() / mult) {
        mLastStatus = STATUS_UNSUPPORTED_OR_OVERFLOWING_INPUT;
        return;
      }
      mTau = mSqrtN * mult;                        // ceil(sqrt(n))*ceil(log2(n))
    } else {
      mTau = mSqrtN;                               // Liquidationist tau = sqrt(n)
    }
  }

  // Matching + dual + forest scaffolding (empty / trivial for Phase 0).
  mMate.assign(mN, scNull());
  mY.assign(mN, 0);
  mZ.clear();
  mBlossomParent.assign(mN, scNull());   // each vertex is its own root blossom
  mW.assign(mN, std::vector<LL>());
  for (Size v = 0; v < mN; ++v) mW[v].assign(mAdj[v].size(), 0);
}

inline void ScalingMatching::compute(const Graph& graph,
                                     std::vector<Size>* mateVec,
                                     Size* card,
                                     Val* edgWght) {
  if (mateVec == NULL || card == NULL || edgWght == NULL) return;
  init(graph);
  Size mN0 = mN;
  // Use runPerfectBlossom (the verified Gate-1b engine: Criterion 1, exact
  // tightness => Property 1 => exact MWPM). With the z-clamp fix it is exactly
  // optimal vs LEMON on all tested graphs up to n=400, weights to 1000. (The
  // dummy-based runScaling/computeFaithful path had a maximal-augment bug and is
  // superseded for the single-scale exact engine; the multi-scale speed loop is
  // future work.)
  bool ok = runPerfectBlossom();
  mateVec->assign(mN0, cNullItm);
  Size c = 0; Val wsum = 0;
  if (ok) {
    for (Size v = 0; v < mN0; ++v) {
      if (mMate[v] != scNull() && mMate[v] < mN0 && mMate[v] > v) {
        Size mate = mMate[v];
        (*mateVec)[v] = mate; (*mateVec)[mate] = v;
        ++c;
        LL best = 0; bool found = false;
        for (Size i = 0; i < mAdj[v].size(); ++i)
          if (mAdj[v][i] == mate && (!found || mWhat[v][i] > best)) { best = mWhat[v][i]; found = true; }
        if (found) wsum += (Val)best;
      }
    }
  }
  *card = c;
  *edgWght = wsum;
}

// ---- Phase 1a: weighted primal-dual search, NO blossoms ----

inline void ScalingMatching::augmentPath1a(Size endpt) {
  // Flip one side of the augmenting path, from EVEN vertex `endpt` up to its
  // free EVEN root. The caller matches the connecting edge (endpt, other)
  // separately. Tree invariant:
  //   parent[ODD o]  = EVEN e   via an UNMATCHED (grow) edge (e,o)
  //   parent[EVEN m] = ODD o    via a MATCHED edge (o,m)  [m is o's old mate]
  // Path from endpt to root:  endpt(E) - o1(O) - e1(E) - o2(O) - e2(E) - ... - root(E)
  // After augment, the new matched edges along this side are each ODD vertex to
  // its EVEN parent: (o1,e1), (o2,e2), ...  (the connecting edge + these).
  // We walk in (EVEN even-vertex, its ODD parent) steps: at EVEN `e` whose
  // parent is ODD `o`, match (o, parent[o]) and continue from parent[o].
  Size e = endpt;
  while (true) {
    Size o = mParent[e];            // ODD parent of e, or scNull if e is root
    if (o == scNull()) break;       // e is the free EVEN root: done this side
    Size pe = mParent[o];           // EVEN parent of o (the grow edge target)
    mMate[o] = pe;                  // new matched edge (o, pe)
    mMate[pe] = o;
    e = pe;                         // continue up from pe (now matched to o)
  }
}

// ===========================================================================
// PAPER DUAL FORMULAS (Duan-Pettie-Su, SODA 2017) — kept here as the contract
// this code (and Phase 1b/2/3) must obey. We are in the MAXIMIZATION direction.
//
//   Duals:  y(v) integer (one per vertex); z(B) even >= 0 (one per blossom).
//   yz(e)  = y(u) + y(v) + sum_{B: e in B} z(B).      [no blossoms in 1a: yz=y(u)+y(v)]
//   Domination (max):  yz(e) >= w(e)   for every edge e.
//   Tightness:         yz(e) == w(e)   (matched / blossom edges).
//   slack(e)           = yz(e) - w(e) >= 0;  eligible (Criterion 1) iff slack==0.
//
//   Even-weight discipline: we run on w2(e) = 2*what(e) so every slack is even,
//   and the EVEN-EVEN augment step (slack/2) is exact. y-values live in this
//   doubled scale. (The paper enforces even w throughout for the same reason.)
//
//   Relaxed Comp. Slackness for the full scaling loop (Property 2, NOT needed in
//   1a but recorded for Phase 2/3): even w, integer y, even z>=0,
//     near-domination yz(e) >= w(e) - 2,  near-tightness yz(e) <= w(e).
//   Per-scale rescale (Phase 3): w <- 2(w' + bit), y <- 2y'+3, z <- 2z'.
//
//   PERFECT matching => NO vertex-drop: every free vertex is always an EVEN tree
//   root and stays in the search; duals only shift, never release a vertex.
// ===========================================================================

inline bool ScalingMatching::runPerfectNoBlossom() {
  mBlossomNeeded = false;
  mNoPerfect = false;
  mMate.assign(mN, scNull());

  if (mN == 0) return true;
  if (mN % 2 != 0) { mNoPerfect = true; return false; }  // odd n: no perfect matching

  // Dual init (DOUBLED scale): y[v] = maxW gives y[u]+y[v] = 2*maxW >= 2*what,
  // so domination holds with even slacks from the start.
  mY.assign(mN, mMaxW);

  mLabel.assign(mN, LBL_FREE_OF_TREE);
  mTreeRoot.assign(mN, scNull());
  mParent.assign(mN, scNull());

  // Each augmenting phase matches at least one previously-free vertex pair, so
  // at most n/2 phases. We stop early as soon as the matching is perfect.
  for (Size iter = 0; iter <= mN; ++iter) {
    // (re)label: EVERY free vertex is an EVEN root (perfect => no drop).
    mLabel.assign(mN, LBL_FREE_OF_TREE);
    mTreeRoot.assign(mN, scNull());
    mParent.assign(mN, scNull());
    bool anyRoot = false;
    for (Size v = 0; v < mN; ++v) {
      if (mMate[v] == scNull()) {
        mLabel[v] = LBL_EVEN;
        mTreeRoot[v] = v;
        anyRoot = true;
      }
    }
    if (!anyRoot) break;   // no free vertex left -> matching is perfect -> done

    bool augmented = false;
    // grow until augment, or stuck (then dual-adjust), repeated.
    while (!augmented && !mBlossomNeeded) {
      // 1) scan for a tight edge (slack==0, Criterion 1) from an EVEN vertex.
      bool grewOrAugmented = false;
      for (Size u = 0; u < mN && !augmented; ++u) {
        if (mLabel[u] != LBL_EVEN) continue;
        for (Size i = 0; i < mAdj[u].size(); ++i) {
          if (slack1a(u, i) != 0) continue;     // only tight edges (Criterion 1)
          Size w = mAdj[u][i];
          if (mLabel[w] == LBL_EVEN) {
            if (mTreeRoot[w] != mTreeRoot[u]) {
              // EVEN-EVEN, different trees -> augmenting path. Flip both sides,
              // then match the connecting edge (u,w).
              augmentPath1a(u);
              augmentPath1a(w);
              mMate[u] = w; mMate[w] = u;
              augmented = true;
              grewOrAugmented = true;
              break;
            } else {
              // EVEN-EVEN, same tree -> would need a blossom. 1a bails.
              mBlossomNeeded = true;
              grewOrAugmented = true;
              break;
            }
          } else if (mLabel[w] == LBL_ODD) {
            continue;                            // already ODD, nothing to do
          } else {
            // w is outside the trees.
            if (mMate[w] == scNull()) {
              // tight edge EVEN -> free vertex => direct augmenting path.
              // Flip u's side up to its root, then match (u,w).
              augmentPath1a(u);
              mMate[u] = w; mMate[w] = u;
              augmented = true;
              grewOrAugmented = true;
              break;
            } else {
              // GROW: w becomes ODD (parent u), its mate becomes EVEN.
              Size m = mMate[w];
              mLabel[w] = LBL_ODD;  mParent[w] = u;  mTreeRoot[w] = mTreeRoot[u];
              mLabel[m] = LBL_EVEN; mParent[m] = w;  mTreeRoot[m] = mTreeRoot[u];
              grewOrAugmented = true;
            }
          }
        }
      }
      if (augmented || mBlossomNeeded) break;
      if (grewOrAugmented) continue;   // made progress; rescan for new tight edges

      // 2) stuck: dual adjustment. For PERFECT matching there is NO vertex-drop
      //    term -- delta is purely the min slack that exposes a new tight edge:
      //      (a) EVEN-EVEN edges across DIFFERENT trees: slack drops by 2*delta
      //      (b) EVEN-> outside (matched, non-tree) edges: slack drops by delta
      //    EVEN y -= delta, ODD y += delta.
      LL delta = -1;
      for (Size u = 0; u < mN; ++u) {
        if (mLabel[u] != LBL_EVEN) continue;
        for (Size i = 0; i < mAdj[u].size(); ++i) {
          Size w = mAdj[u][i];
          if (mLabel[w] == LBL_EVEN) {
            if (mTreeRoot[w] != mTreeRoot[u]) {
              LL half = slack1a(u, i) / 2;       // both EVEN: slack drops by 2 delta
              if (delta < 0 || half < delta) delta = half;
            }
            // same tree EVEN-EVEN: blossom on tightness; handled in scan above
          } else if (mLabel[w] != LBL_ODD) {
            LL s = slack1a(u, i);                 // EVEN->outside: slack drops by delta
            if (delta < 0 || s < delta) delta = s;
          }
        }
      }
      if (delta < 0) {
        // No EVEN-incident edge can ever become tight, yet free vertices remain:
        // the free EVEN roots cannot reach a partner -> NO perfect matching.
        mNoPerfect = true;
        return false;
      }
      if (delta == 0) {  // safety: a tight edge slipped through the scan; rescan
        continue;
      }
      // apply: EVEN y -= delta, ODD y += delta
      for (Size v = 0; v < mN; ++v) {
        if (mLabel[v] == LBL_EVEN) mY[v] -= delta;
        else if (mLabel[v] == LBL_ODD) mY[v] += delta;
      }
    }
    if (mBlossomNeeded) return false;
    if (!augmented) {
      // a phase that couldn't augment but free vertices remain -> no perfect M
      mNoPerfect = true;
      return false;
    }
  }
  return !mBlossomNeeded && !mNoPerfect;
}

inline bool ScalingMatching::computeNoBlossom(const Graph& graph,
                                              std::vector<Size>* mateVec,
                                              Size* card, Val* edgWght,
                                              bool* blossomNeeded, bool* noPerfect) {
  init(graph);
  bool ok = runPerfectNoBlossom();
  if (blossomNeeded) *blossomNeeded = mBlossomNeeded;
  if (noPerfect)     *noPerfect     = mNoPerfect;
  // extract matching
  mateVec->assign(mN, cNullItm);
  Size c = 0;
  Val wsum = 0;
  for (Size v = 0; v < mN; ++v) {
    if (mMate[v] != scNull() && mMate[v] > v) {
      (*mateVec)[v] = mMate[v];
      (*mateVec)[mMate[v]] = v;
      ++c;
      // find original weight of edge (v, mate)
      // The matched pair (v, mate) may have PARALLEL edges with different
      // weights; the search effectively uses the heaviest one (its slack goes
      // tight first), so the matched edge's contribution is the MAX parallel.
      Size mate = mMate[v];
      LL best = 0; bool found = false;
      for (Size i = 0; i < mAdj[v].size(); ++i)
        if (mAdj[v][i] == mate && (!found || mWhat[v][i] > best))
          { best = mWhat[v][i]; found = true; }
      if (found) wsum += (Val)best;
    }
  }
  if (card) *card = c;
  if (edgWght) *edgWght = wsum;
  return ok;
}

inline bool ScalingMatching::computePerfect(const Graph& graph,
                                            std::vector<Size>* mateVec,
                                            Size* card, Val* edgWght,
                                            bool* noPerfect) {
  init(graph);
  bool ok = runPerfectBlossom();
  if (noPerfect) *noPerfect = mNoPerfect;
  mateVec->assign(mN, cNullItm);
  Size c = 0; Val wsum = 0;
  if (ok) {
    for (Size v = 0; v < mN; ++v) {
      if (mMate[v] != scNull() && mMate[v] > v) {
        (*mateVec)[v] = mMate[v];
        (*mateVec)[mMate[v]] = v;
        ++c;
        Size mate = mMate[v];
        LL best = 0; bool found = false;
        for (Size i = 0; i < mAdj[v].size(); ++i)
          if (mAdj[v][i] == mate && (!found || mWhat[v][i] > best))
            { best = mWhat[v][i]; found = true; }
        if (found) wsum += (Val)best;
      }
    }
  }
  if (card) *card = c;
  if (edgWght) *edgWght = wsum;
  return ok;
}


// ===========================================================================
// PAPER-FAITHFUL integer Edmonds search (DPS arXiv:1411.1919, Fig 2/3).
// Contract: docs/Scaling_Paper_Exact_Spec.md. Weights DOUBLED (w2=2*what).
// ===========================================================================

inline Size ScalingMatching::outerOf(Size v) const {
  Size c = v;
#ifdef CYCDBG
  Size steps = 0;
  while (mBParent[c] != scNull()) {
    c = mBParent[c];
    if (++steps > mNumNodes + 8) {
      std::fprintf(stderr, "[CYCDBG] outerOf(%lu) parent-chain CYCLE: ",
                   (unsigned long)v);
      Size d = v; for (Size k = 0; k < 12 && d != scNull(); ++k) { std::fprintf(stderr, "%lu->", (unsigned long)d); d = mBParent[d]; }
      std::fprintf(stderr, "...\n"); std::fflush(stderr);
      return c;   // bail to avoid the hang while debugging
    }
  }
  return c;
#else
  while (mBParent[c] != scNull()) { c = mBParent[c];
#ifdef OPCNT
    ++gOpOuterSteps;
#endif
  }
  return c;
#endif
}

// UNION-FIND (S5 mechanism). ufFind(v) = representative VERTEX of v's current outer blossom (the
// base vertex), with path-halving, O(alpha). Used for the O(alpha) same-outer-blossom test on the
// edge-scan hot path (ufFind(u)==ufFind(w)) and as the foundation for no-relabel shrink. Maintained
// in shrink (union B's leaves into base rep) + dissolve (re-split). When UF is empty (size 0),
// callers use the exact outerOf slow path, so this is inert until initialized.
inline Size ScalingMatching::ufFind(Size v) const {
  if (v >= mUf.size()) return v;
  while (mUf[v] != v) { mUf[v] = mUf[mUf[v]]; v = mUf[v]; }  // path-halving (mUf is mutable)
  return v;
}
// O(alpha) outer-node resolution: ufFind(v) -> base-vertex rep -> mBaseToNode -> node id.
// Replaces the O(depth) outerOf mBParent-walk in the hot grow loop. Falls back to outerOf if the
// map is empty or stale (rep not mapped). Correct iff mBaseToNode is maintained (shrink/dissolve);
// UFCHECK-verified equal to outerOf.
inline Size ScalingMatching::outerOfUF(Size v) const {
  if (mUf.empty() || mBaseToNode.empty()) return outerOf(v);
  return outerOfUFRep(ufFind(v), v);
}
// Same as outerOfUF but with ufFind(v) ALREADY computed by the caller (the grow loop computes
// ufFind(w) for the same-blossom test, then needs the outer node -- this avoids re-running ufFind(w)
// a SECOND time per edge, a measured redundancy). Precondition: rep == ufFind(v), UF/map non-empty.
inline Size ScalingMatching::outerOfUFRep(Size rep, Size v) const {
  Size node = (rep < mBaseToNode.size()) ? mBaseToNode[rep] : scNull();
  if (node != scNull() && node < mNumNodes &&
      mBParent[node] == scNull() &&
      (node < mN || (node < mBActive.size() && mBActive[node])))
    return node;
#ifdef OPCNT
  ++gOpUfSlowPath;
#endif
  return outerOf(v);   // exact representation slow path
}
inline void ScalingMatching::ufUnionInto(Size member, Size baseRep) {
  if (member >= mUf.size() || baseRep >= mUf.size()) return;
  Size r = ufFind(member); mUf[r] = baseRep; mUf[baseRep] = baseRep;
}
inline Size ScalingMatching::ufOuter(Size v) {
  // outer node of v: from v's uf-rep (base vertex), walk mBParent to the top. (Increment 1 keeps
  // this for parity-checking against outerOf; the hot path uses ufFind for the same-blossom test.)
  Size base = ufFind(v);
  return outerOf(base);
}

#ifdef UFCHECK
// Parity check (debug): for every real vertex, ufOuter(v) must equal outerOf(v). Run after each
// search to validate the UF foundation before the augment relies on it. Prints + aborts on mismatch.
inline void ScalingMatching::ufVerify(const char* where) {
  for (Size v = 0; v < mN; ++v) {
    Size b = outerOf(v);
    Size a = ufOuter(v);
    if (a != b) { std::fprintf(stderr, "[UFCHECK %s] v=%lu ufOuter=%lu outerOf=%lu MISMATCH\n",
                  where,(unsigned long)v,(unsigned long)a,(unsigned long)b); std::fflush(stderr); std::abort(); }
    Size c = outerOfUF(v);
    if (c != b) { std::fprintf(stderr, "[UFCHECK %s] v=%lu outerOfUF=%lu outerOf=%lu MISMATCH\n",
                  where,(unsigned long)v,(unsigned long)c,(unsigned long)b); std::fflush(stderr); std::abort(); }
  }
}
#endif

// Current weight w(e) of the heaviest parallel edge (u,v); 0 if no real edge.
// Source = mW[u][idx] when mUseScaleW (scaling loop), else 2*mWhat (Gate-1a/1b).
inline ScalingMatching::LL
ScalingMatching::w2max(Size u, Size v) const {
  LL w2 = 0; bool found = false;
  for (Size i = 0; i < mAdj[u].size(); ++i) {
#ifdef EDMPROF
    extern long long gW2Calls, gW2Scan; ++gW2Scan;
#endif
    if (mAdj[u][i] == v) {
      LL wuv = mUseScaleW ? mW[u][i] : 2 * mWhat[u][i];
      if (!found || wuv > w2) { w2 = wuv; found = true; }
    }
  }
#ifdef EDMPROF
  { extern long long gW2Calls; ++gW2Calls; }
#endif
  return found ? w2 : 0;
}

inline ScalingMatching::LL
ScalingMatching::matchedEdgeW2(Size u, Size v) {
  // The scale weight arrays are immutable during one scale. A matching may
  // leave and later return to the same endpoint pair without invalidating the
  // cached edge weight.
  if (mMatchedWCacheScale != mAuditScale ||
      mMatchedWCacheMode != mUseScaleW ||
      mMatchedWCacheMate.size() < mMate.size()) {
    mMatchedWCacheMate.assign(mMate.size(), scNull());
    mMatchedWCacheValue.resize(mMate.size(), 0);
    mMatchedWCacheScale = mAuditScale;
    mMatchedWCacheMode = mUseScaleW;
  }
  if (u < mMatchedWCacheMate.size() && mMatchedWCacheMate[u] == v)
    return mMatchedWCacheValue[u];
  LL value = w2max(u, v);
  if (u < mMatchedWCacheMate.size()) {
    mMatchedWCacheMate[u] = v;
    mMatchedWCacheValue[u] = value;
  }
  if (v < mMatchedWCacheMate.size()) {
    mMatchedWCacheMate[v] = u;
    mMatchedWCacheValue[v] = value;
  }
  return value;
}

// Section 5.2 timestamp state on the canonical forest. Closed role intervals
// accumulate in mYAcc/mZAcc; the current interval is represented by role+stamp.
// This keeps a dual jump O(1), and lets grow/shrink settle only affected nodes.
inline void ScalingMatching::tsInit() {
  Size cap = mZ.size();
  if (mY0.size() < mN) mY0.resize(mN);
  if (mYAcc.size() < mN) mYAcc.resize(mN);
  if (mVStamp.size() < mN) mVStamp.resize(mN);
  if (mVRole.size() < mN) mVRole.resize(mN, LBL_FREE_OF_TREE);
  if (mLazyVEpoch.size() < mN) mLazyVEpoch.resize(mN, 0);
  if (mZ0.size() < cap) mZ0.resize(cap);
  if (mZAcc.size() < cap) mZAcc.resize(cap);
  if (mNStamp.size() < cap) mNStamp.resize(cap);
  if (mNRole.size() < cap) mNRole.resize(cap, LBL_FREE_OF_TREE);
  if (mTRoot.size() < cap) mTRoot.resize(cap);
  if (mTIn.size() < cap) mTIn.resize(cap);
  if (mTOut.size() < cap) mTOut.resize(cap);
  if (mDelta.size() < cap) mDelta.resize(cap);
  if (mLazyNEpoch.size() < cap) mLazyNEpoch.resize(cap, 0);
  if (++mLazyEpochClock == 0) {
    std::fill(mLazyVEpoch.begin(), mLazyVEpoch.end(), 0);
    std::fill(mLazyNEpoch.begin(), mLazyNEpoch.end(), 0);
    ++mLazyEpochClock;
  }
  mLazyActiveVertices.clear();
  mLazyActiveNodes.clear();
  mLazyDuals = true;

  // beginCompactEventLabels has already built precisely the roots of this
  // fresh search. Activate those roots only; later grow/shrink operations
  // activate additional records on first use.
  for (Size i = 0; i < mEventInitialRoots.size(); ++i) {
    Size root = mEventInitialRoots[i];
    if (root < mN) lazySetVertexRole(root, LBL_EVEN);
    else if (root < cap) lazySetNodeRole(root, LBL_EVEN);
  }
}

inline void ScalingMatching::tsEnsureVertex(Size v) {
  if (!mLazyDuals || v >= mN ||
      (v < mLazyVEpoch.size() && mLazyVEpoch[v] == mLazyEpochClock)) return;
  mLazyVEpoch[v] = mLazyEpochClock;
  mY0[v] = mY[v]; mYAcc[v] = 0;
  mVStamp[v] = mTNow; mVRole[v] = LBL_FREE_OF_TREE;
  mTRoot[v] = -1; mTIn[v] = -1; mTOut[v] = -1; mDelta[v] = 0;
  mLazyActiveVertices.push_back(v);
  ++mCounters.lazyDualVertexActivations;
}

inline void ScalingMatching::tsEnsureNode(Size B) {
  if (!mLazyDuals || B < mN) return;
  if (B >= mNRole.size()) {
    Size cap = mZ.size();
    mZ0.resize(cap); mZAcc.resize(cap); mNStamp.resize(cap);
    mNRole.resize(cap, LBL_FREE_OF_TREE);
    mTRoot.resize(cap); mTIn.resize(cap); mTOut.resize(cap); mDelta.resize(cap);
    mLazyNEpoch.resize(cap, 0);
  }
  if (mLazyNEpoch[B] == mLazyEpochClock) return;
  mLazyNEpoch[B] = mLazyEpochClock;
  mZ0[B] = mZ[B]; mZAcc[B] = 0;
  mNStamp[B] = mTNow; mNRole[B] = LBL_FREE_OF_TREE;
  mTRoot[B] = -1; mTIn[B] = -1; mTOut[B] = -1; mDelta[B] = 0;
  mLazyActiveNodes.push_back(B);
  ++mCounters.lazyDualNodeActivations;
  // Node-role intervals change every leaf dual. Activate exactly those leaves
  // now so final materialization remains proportional to participating state.
#ifdef SCALING_LEAF_NOVEC
  Size end = B < mLfTailX.size() ? mLfTailX[B] : scNull();
  for (Size v = B < mLfHeadX.size() ? mLfHeadX[B] : scNull(), next = scNull();
       v != scNull(); v = next) {
    next = (v == end) ? scNull() : mLfNextX[v];
    tsEnsureVertex(v);
  }
#else
  if (B < mNodeLeaves.size()) {
    const std::vector<Size>& leaves = mNodeLeaves[B];
    for (Size i = 0; i < leaves.size(); ++i) tsEnsureVertex(leaves[i]);
  }
#endif
}

inline ScalingMatching::LL ScalingMatching::yOf(Size v) const {
  if (!mLazyDuals || v >= mY0.size()) return mY[v];
  Size B = outerOf(v);
  return yOfOuter(v, B);
}

// Materialize one vertex dual when its current outer node is already known.
// Event scheduling and pop validation resolve that node for role/ownership
// checks first; repeating the potentially deep laminar-parent walk in yOf was
// pure representation overhead.
inline ScalingMatching::LL
ScalingMatching::yOfOuter(Size v, Size B) const {
  if (!mLazyDuals || v >= mN) return mY[v];
  bool vertexActive = v < mLazyVEpoch.size() &&
                      mLazyVEpoch[v] == mLazyEpochClock;
  LL y = vertexActive ? mY0[v] + mYAcc[v] : mY[v];
  int role; LL stamp;
  if (B == v) {
    role = vertexActive ? mVRole[v] : LBL_FREE_OF_TREE;
    stamp = vertexActive ? mVStamp[v] : mTNow;
  } else {
    bool nodeActive = B < mLazyNEpoch.size() &&
                      mLazyNEpoch[B] == mLazyEpochClock;
    role = nodeActive ? mNRole[B] : LBL_FREE_OF_TREE;
    stamp = nodeActive ? mNStamp[B] : mTNow;
  }
  if (role == LBL_ODD) y += mTNow - stamp;
  else if (role == LBL_EVEN) y -= mTNow - stamp;
  return y;
}

inline ScalingMatching::LL ScalingMatching::zOf(Size B) const {
  if (!mLazyDuals || B >= mZ.size()) return mZ[B];
  if (B >= mLazyNEpoch.size() || mLazyNEpoch[B] != mLazyEpochClock)
    return mZ[B];
  LL z = mZ0[B] + mZAcc[B];
  int role = mNRole[B]; LL stamp = mNStamp[B];
  if (role == LBL_ODD) z -= 2 * (mTNow - stamp);
  else if (role == LBL_EVEN) z += 2 * (mTNow - stamp);
  return z;
}

inline void ScalingMatching::lazySetVertexRole(Size v, int role) {
  if (!mLazyDuals || v >= mN) return;
  tsEnsureVertex(v);
  LL span = mTNow - mVStamp[v];
  if (mVRole[v] == LBL_ODD) mYAcc[v] += span;
  else if (mVRole[v] == LBL_EVEN) mYAcc[v] -= span;
  mVRole[v] = role; mVStamp[v] = mTNow;
}

inline void ScalingMatching::lazySetNodeRole(Size B, int role) {
  if (!mLazyDuals || B < mN) return;
  tsEnsureNode(B);
  LL span = mTNow - mNStamp[B];
  int old = mNRole[B];
  if (old == LBL_ODD || old == LBL_EVEN) {
    LL rate = old == LBL_ODD ? 1 : -1;
#ifdef SCALING_LEAF_NOVEC
    Size end = B < mLfTailX.size() ? mLfTailX[B] : scNull();
    for (Size v = B < mLfHeadX.size() ? mLfHeadX[B] : scNull(), next = scNull();
         v != scNull(); v = next) {
      next = (v == end) ? scNull() : mLfNextX[v];
      mYAcc[v] += rate * span;
    }
#else
    const std::vector<Size>& leaves = mNodeLeaves[B];
    for (Size i = 0; i < leaves.size(); ++i) mYAcc[leaves[i]] += rate * span;
#endif
  }
  if (old == LBL_ODD) mZAcc[B] -= 2 * span;
  else if (old == LBL_EVEN) mZAcc[B] += 2 * span;
  mNRole[B] = role; mNStamp[B] = mTNow;
  if (role == LBL_ODD || role == LBL_EVEN) {
    if (mTRoot[B] < 0) mTRoot[B] = mTNow;
    if (role == LBL_ODD) mTIn[B] = mTNow;
    else mTOut[B] = mTNow;
  }
}

inline void ScalingMatching::materializeLazyDuals() {
  if (!mLazyDuals) return;
  ++mCounters.dualMaterializations;
  mCounters.lazyDualVerticesMaterialized += mLazyActiveVertices.size();
  mCounters.lazyDualNodesMaterialized += mLazyActiveNodes.size();
  for (Size i = 0; i < mLazyActiveVertices.size(); ++i) {
    Size v = mLazyActiveVertices[i];
    mY[v] = yOf(v);
  }
  for (Size i = 0; i < mLazyActiveNodes.size(); ++i) {
    Size B = mLazyActiveNodes[i];
    mZ[B] = zOf(B);
  }
  mLazyDuals = false;
}

inline void ScalingMatching::tsOnDissolveInner(Size B) {
  if (!mLazyDuals || B >= mBCycle.size()) return;
  tsEnsureNode(B);
  ++mCounters.dissolutionTimestampTransfers;
  LL innerSpan = (B < mTIn.size() && mTIn[B] >= 0)
               ? mTNow - mTIn[B] : 0;
  lazySetNodeRole(B, LBL_FREE_OF_TREE);
  const std::vector<Size>& children = mBCycle[B];
  for (Size i = 0; i < children.size(); ++i) {
    Size c = children[i]; int role = mNLabel[c];
    if (c < mN) tsEnsureVertex(c); else tsEnsureNode(c);
    if (c < mDelta.size()) {
      LL nextDelta = 0;
      // Section 5.2: Delta(child) <- Delta(B) +
      // (t_now - t_in(B)).  The child's previous value describes an older
      // root interval and must be replaced, not incremented.
      if (!checkedAddLL(mDelta[B], innerSpan, &nextDelta)) {
        mLastStatus = STATUS_UNSUPPORTED_OR_OVERFLOWING_INPUT;
        return;
      }
      mDelta[c] = nextDelta;
      mTRoot[c] = mTNow;
    }
    if (c < mN) lazySetVertexRole(c, role);
    else lazySetNodeRole(c, role);
    if (c < mTIn.size()) {
      if (role == LBL_ODD) mTIn[c] = mTNow;
      else if (role == LBL_EVEN) mTOut[c] = mTNow;
    }
  }
}
inline void ScalingMatching::tsOnFormOuter(Size B) {
  if (!mLazyDuals || B >= mBCycle.size()) return;
  const std::vector<Size>& children = mBCycle[B];
  for (Size i = 0; i < children.size(); ++i) {
    Size c = children[i];
    if (c < mN) tsEnsureVertex(c); else tsEnsureNode(c);
    int oldRole = c < mN ? mVRole[c] : mNRole[c];
    if (oldRole == LBL_ODD &&
        c < mTIn.size() && mTIn[c] >= 0) {
      LL nextDelta = 0;
      if (!checkedAddLL(mDelta[c], mTNow - mTIn[c], &nextDelta)) {
        mLastStatus = STATUS_UNSUPPORTED_OR_OVERFLOWING_INPUT;
        return;
      }
      mDelta[c] = nextDelta;
      mTOut[c] = mTNow;
    }
    if (c < mN) lazySetVertexRole(c, LBL_FREE_OF_TREE);
    else lazySetNodeRole(c, LBL_FREE_OF_TREE);
  }
  tsEnsureNode(B);
  mTRoot[B] = mTNow;
  mNRole[B] = LBL_FREE_OF_TREE; mNStamp[B] = mTNow;
  lazySetNodeRole(B, LBL_EVEN);
}

// slackZ(u,v) = yz(e) - w2(e) = y(u)+y(v)+sum z(B enclosing both) - w2(u,v).
inline ScalingMatching::LL
ScalingMatching::slackZ(Size u, Size v) const {
  return slackZW(u, v, w2max(u, v));
}

inline ScalingMatching::LL
ScalingMatching::oldCommonZ(Size u, Size v) const {
  if (!mOldForestSeparated || u >= mOldContaining.size() ||
      v >= mOldContaining.size()) return 0;
  const std::vector<Size>& a = mOldContaining[u];
  const std::vector<Size>& b = mOldContaining[v];
  Size i = a.size(), j = b.size();
  LL zsum = 0;
  // Ancestors are stored inner-to-outer, so common enclosing old blossoms
  // are exactly the common suffix of the two lists.
  while (i > 0 && j > 0 && a[i - 1] == b[j - 1]) {
    Size B = a[i - 1];
    if (B < mZ.size()) zsum += mZ[B];
    --i; --j;
  }
  return zsum;
}

// slackZ core given the edge weight w2 already in hand (EDM WIN: callers in the grow loop and
// dual-adjust scan iterate mAdj[u] and KNOW the weight is mW[u][k] (scaling) / 2*mWhat[u][k];
// passing it here avoids w2max's O(deg) RE-SCAN of mAdj[u] to re-find that same edge -- ~400M
// redundant scan-iters at n=6400, a quadratic O(n*deg^2) hot path). The blossom-z-sum part is
// O(blossom-depth) and unchanged.
inline ScalingMatching::LL
ScalingMatching::slackZW(Size u, Size v, LL w2) const {
  LL s = yOf(u) + yOf(v) + oldCommonZ(u, v) - w2;
  // O(alpha) SHORT-CIRCUIT (the cold-scale fix): blossoms enclosing BOTH u and v are common
  // ancestors in the blossom forest -> they exist ONLY if u and v are in the SAME outer node
  // (outer nodes are roots: mBParent==scNull, so distinct outer nodes share NO ancestor). If
  // ufFind(u)!=ufFind(v) (different outer node, O(alpha)), the z-sum is 0 -> return immediately,
  // skipping the O(depth) chain walk. (UF maintained == outerOf; UFCHECK-verified. This killed
  // the 1.69B slackSteps at n=12800 scale 2 -- cross-blossom edges no longer walk to root.)
  if (!mUf.empty()) {
    Size ru = (u < mUf.size()) ? ufFind(u) : u;
    Size rv = (v < mUf.size()) ? ufFind(v) : v;
    if (ru != rv) return s;   // different outer nodes -> no common enclosing blossom -> z-sum 0
  }
  // Sum z(B) over blossoms B that enclose BOTH u and v. Such B form a common
  // suffix of u's and v's parent chains (the chain from the lowest common
  // enclosing blossom up to the outermost). Walk both chains to equal depth, then
  // step up together until they meet; every node from the meet upward encloses
  // both. This is O(depth), not O(depth^2) (the old per-b containment rescan).
  Size a = mBParent[u], b = mBParent[v];
  if (a == scNull() || b == scNull()) return s;   // a bare vertex shares no blossom
  Size da = 0, db = 0;
#ifdef CYCDBG
  for (Size x = a; x != scNull(); x = mBParent[x]) { if (++da > mNumNodes + 8) { std::fprintf(stderr, "[CYCDBG] slackZ a-chain CYCLE from u=%lu (a=%lu)\n", (unsigned long)u, (unsigned long)a); std::fflush(stderr); break; } }
  db = 0;
  for (Size x = b; x != scNull(); x = mBParent[x]) { if (++db > mNumNodes + 8) { std::fprintf(stderr, "[CYCDBG] slackZ b-chain CYCLE from v=%lu (b=%lu)\n", (unsigned long)v, (unsigned long)b); std::fflush(stderr); break; } }
#else
  for (Size x = a; x != scNull(); x = mBParent[x]) ++da;
  for (Size x = b; x != scNull(); x = mBParent[x]) ++db;
#ifdef OPCNT
  gOpSlackSteps += (long long)da + (long long)db;
#endif
#endif
  while (da > db) { a = mBParent[a]; --da; }
  while (db > da) { b = mBParent[b]; --db; }
#ifdef CYCDBG
  { Size mg = 0; while (a != b) { if (++mg > mNumNodes + 8) { std::fprintf(stderr, "[CYCDBG] slackZ MERGE never converges u=%lu v=%lu\n", (unsigned long)u, (unsigned long)v); std::fflush(stderr); break; } a = mBParent[a]; b = mBParent[b]; } }
#else
  while (a != b) { a = mBParent[a]; b = mBParent[b]; }
#endif
  // a == b == lowest common enclosing blossom (or scNull if none). Sum from here up.
  for (Size x = a; x != scNull(); x = mBParent[x]) s += zOf(x);
  return s;
}

// Eligibility on a REAL edge (u,v) under the given criterion. yz=w2 means
// slackZ==0; yz=w2-2 means slackZ==-2. (e in E_B handled by caller via outerOf
// equality -- an edge inside a blossom is implicitly eligible.)
inline bool ScalingMatching::eligible(Size u, Size v, int crit) const {
  LL s = slackZ(u, v);
  if (crit == CRIT1) return s == 0;
  if (crit == CRIT3) return s == 0 || s == -2;
  // CRIT2: (e notin M and yz=w-2) or (e in M and yz=w).  [E_B case: caller skips]
  bool matched = (mMate[u] == v);
  if (matched) return s == 0;
  return s == -2;
}

// eligibleAt: same as eligible(u, mAdj[u][idx], crit) but reads the edge WEIGHT directly from
// mW[u][idx] (scaling) / 2*mWhat[u][idx] -- O(1), no w2max re-scan. NOTE: uses the weight at
// THIS slot, not w2max's "heaviest parallel edge"; matchbox Graph dedups parallels on input
// (the scaling engine builds one slot per (u,v)), so this is the same value. The EDM-hot
// grow/dual-adjust loops use this.
inline bool ScalingMatching::eligibleAt(Size u, Size idx, int crit) const {
  Size v = mAdj[u][idx];
  LL w2 = mUseScaleW ? mW[u][idx] : 2 * mWhat[u][idx];
  LL s = slackZW(u, v, w2);
  if (crit == CRIT1) return s == 0;
  if (crit == CRIT3) return s == 0 || s == -2;
  bool matched = (mMate[u] == v);
  if (matched) return s == 0;
  return s == -2;
}
// CROSS-BLOSSOM eligibility: caller (grow loop) already PROVED ufFind(v)!=ufFind(u) (the edge
// passed the UF same-blossom `continue`), so slackZW's z-sum is 0 and slack == mY[u]+mY[v]-w2
// exactly. Skips slackZW's 2 ufFind calls + early-return chain walk per edge (ufFind = the #1
// hotspot, ~11-14% of runtime). MUST only be called when that guard has run with UF active.
inline bool ScalingMatching::eligibleAtCross(Size u, Size idx, int crit) const {
  Size v = mAdj[u][idx];
  LL w2 = mUseScaleW ? mW[u][idx] : 2 * mWhat[u][idx];
  // Distinct current Omega nodes can still share an inherited Omega_0 node.
  LL s = mY[u] + mY[v] + oldCommonZ(u, v) - w2;
  if (crit == CRIT1) return s == 0;
  if (crit == CRIT3) return s == 0 || s == -2;
  bool matched = (mMate[u] == v);
  if (matched) return s == 0;
  return s == -2;
}

// LCA of outer nodes a,b in the alternating forest (via mNParent). scNull if
// different trees.
inline Size ScalingMatching::lcaOuter(Size a, Size b) {
  if (mLcaSeen.size() < mNumNodes) {
    Size target = mLcaSeen.empty() ? (Size)1024 : mLcaSeen.size();
    const Size limit = mBParent.size();
    while (target < mNumNodes) {
      Size doubled = target <= std::numeric_limits<Size>::max() / 2
                   ? target * 2 : mNumNodes;
      target = doubled < mNumNodes ? mNumNodes : doubled;
      if (target > limit) { target = limit; break; }
    }
    if (target < mNumNodes) target = mNumNodes;
    mLcaSeen.resize(target, 0);
  }
  if (++mLcaClock == 0) {
    std::fill(mLcaSeen.begin(), mLcaSeen.end(), 0);
    ++mLcaClock;
  }
  // Walk a's mNParent chain to its root, marking. The cycle guards below are NOT debug-only:
  // under the PERSISTENT forest (persistAcrossAugment), mNParent can transiently contain a
  // cycle (a tree pointer left pointing into a structure that was rewired by an in-grow
  // augment/dissolve). Walking such a chain in plain for(;;) NEVER TERMINATES and reads only
  // VALID in-bounds indices, so it is invisible to ASan -- THIS was the persist n=1600 hang
  // (gdb backtrace: lcaOuter:1081 b-chain). A cycle => a,b are not joined by a simple
  // mNParent path => no clean LCA => return scNull (caller treats as different-tree: augment
  // or skip, never shrink), which is the safe, correctness-preserving answer.
  // Under incremental shrink (mIncrShrink) navigate the CONTRACTED tree via treeParent =
  // outerOf(mNParent[x]) (S5 :451), so the LCA is an OUTER node reachable by the shrink cycle
  // walk's identical stepping. Non-incr path keeps the raw mNParent walk => byte-identical.
  Size x = a; LL g1 = 0;
  for (;;) { mLcaSeen[x] = mLcaClock;
    Size p = mNParent[x]; if (p == scNull()) break;
    Size nx = mIncrShrink ? outerOf(p) : p;
    if (nx == x) break;                                 // absorbed into self -> root
    x = nx;
    if (++g1 > (LL)mNumNodes + 5) return scNull(); }   // a-chain cycle: no clean LCA
  Size y = b; LL g2 = 0;
  for (;;) { if (mLcaSeen[y] == mLcaClock) return y;
    Size p = mNParent[y]; if (p == scNull()) break;
    Size ny = mIncrShrink ? outerOf(p) : p;
    if (ny == y) break;
    y = ny;
    if (++g2 > (LL)mNumNodes + 5) return scNull(); }   // b-chain cycle: no clean LCA
  return scNull();
}

// Rematch around blossom B so the sub-node carrying vertex newBaseVtx becomes
// the base (cycle index 0), flipping the even-length cycle segment.  This is an
// explicit post-order walk: valid inputs can contain thousands of nested
// blossoms, so using the process call stack here is not safe at HPC sizes.
inline void ScalingMatching::rotateToBase(Size B, Size newBaseVtx) {
  ScalingNanosecondGuard rotationTimer(&mCounters.blossomRotationNanoseconds);
  if (B < mN) return;
  ++mCounters.blossomRotationCalls;
  // rotateToBase is idempotent.  The base invariant is recursive
  // (mBBase[B] == mBBase[mBCycle[B][0]]), so equality here means the complete
  // nested internal matching is already oriented to newBaseVtx.  Relinking a
  // deep dissolution chain asks for this same orientation frequently; walking
  // thousands of tail frames in that case performs no matching change.
  if (B < mBBase.size() && mBBase[B] == newBaseVtx) {
    ++mCounters.blossomRotationNoopCalls;
    return;
  }
  struct RotateFrame {
    Size B, newBase, pos, pair;
    unsigned char stage;
    // Path from newBase upward through (and including) the direct child of B.
    // Tail recursion moves the remaining prefix into the child frame, avoiding
    // a fresh root-to-ancestor climb at every level of a nested chain.
    std::vector<Size> basePath;
    RotateFrame(Size b, Size nb) : B(b), newBase(nb), pos(0), pair(0), stage(0) {}
  };
  std::vector<RotateFrame> todo;
  todo.reserve(64);
  todo.push_back(RotateFrame(B, newBaseVtx));
  ++mCounters.blossomRotationFrames;

  while (!todo.empty()) {
    RotateFrame& f = todo.back();
    std::vector<Size>& cyc = mBCycle[f.B];
    std::vector<Size>& eu  = mBCycEdgeU[f.B];
    std::vector<Size>& ev  = mBCycEdgeV[f.B];

    if (f.stage == 0) {
      CRUMB("rotateToBase", f.B, f.newBase);
      Size k = cyc.size();
      if (k == 0) { todo.pop_back(); continue; }
      Size target = f.basePath.empty() ? f.newBase : f.basePath.back();
      // Under the persistent forest mBParent can transiently contain a cycle.
      // Retain the old unconditional guard, but apply it to each explicit frame.
      LL gg = 0;
      if (f.basePath.empty()) {
        while (target != scNull() && mBParent[target] != f.B) {
          f.basePath.push_back(target);
          target = mBParent[target];
          if (++gg > (LL)mNumNodes + 8) break;
        }
        mCounters.blossomRotationParentSteps += (unsigned long long)gg;
        if (target != scNull() && gg <= (LL)mNumNodes + 8)
          f.basePath.push_back(target);
      }
      if (target == scNull() || gg > (LL)mNumNodes + 8) {
#ifdef SCDBG
        std::printf("    !!rotateToBase B=%lu newBaseVtx=%lu NOT under B\n",
                    (unsigned long)f.B, (unsigned long)f.newBase);
#endif
        todo.pop_back();
        continue;
      }
      // The discovered base path is also an upper bound on the additional
      // tail frames this rotation can require.  Grow the explicit stack once
      // to that bound instead of repeatedly moving every live RotateFrame as
      // a deeply nested path crosses capacities 64,128,256,... .  Restarting
      // stage 0 is safe: basePath is now populated, so no ancestry is rescanned.
      Size frameNeed = todo.size() + f.basePath.size();
      if (frameNeed > todo.capacity() && frameNeed <= mNumNodes + 16) {
        todo.reserve(frameNeed);
        continue;
      }
      if (todo.size() > mCounters.blossomRotationMaxDepth)
        mCounters.blossomRotationMaxDepth = todo.size();
      Size pos = 0; while (pos < k && cyc[pos] != target) ++pos;
      if (pos == k) { todo.pop_back(); continue; }
      // The base must move to an EVEN position.  If it is odd, traverse the
      // symmetric direction by reversing the cycle and its oriented edge rows.
      if (pos % 2 == 1) {
        std::vector<Size> rc(k), reu(k), rev2(k);
        rc[0] = cyc[0];
        for (Size j = 1; j < k; ++j) rc[j] = cyc[k - j];
        for (Size j = 0; j < k; ++j) {
          Size oe = (k - 1 - j) % k;
          reu[j] = ev[oe];
          rev2[j] = eu[oe];
        }
        cyc.swap(rc); eu.swap(reu); ev.swap(rev2);
        pos = 0; while (pos < k && cyc[pos] != target) ++pos;
        if (pos == k) { todo.pop_back(); continue; }
      }
      f.pos = pos;
      f.pair = 0;
      f.stage = (pos == 0) ? 4 : 1;
      continue;
    }

    if (f.stage == 1) {                  // descend into the pair's left sub-node
      Size child = cyc[f.pair], base = eu[f.pair];
      f.stage = 2;
      if (child >= mN && todo.size() <= mNumNodes + 16)
        { todo.push_back(RotateFrame(child, base));
          ++mCounters.blossomRotationFrames; }
      continue;
    }
    if (f.stage == 2) {                  // then its right sub-node
      Size child = cyc[f.pair + 1], base = ev[f.pair];
      f.stage = 3;
      if (child >= mN && todo.size() <= mNumNodes + 16)
        { todo.push_back(RotateFrame(child, base));
          ++mCounters.blossomRotationFrames; }
      continue;
    }
    if (f.stage == 3) {                  // recursive children are complete: rematch pair
      Size a = eu[f.pair], b = ev[f.pair];
      mMate[a] = b; mMate[b] = a;
      f.pair += 2;
      if (f.pair + 1 <= f.pos) { f.stage = 1; continue; }

      Size k = cyc.size();
      std::vector<Size> nc(k), neu(k), nev(k);
      for (Size i = 0; i < k; ++i) {
        nc[i] = cyc[(f.pos + i) % k];
        neu[i] = eu[(f.pos + i) % k];
        nev[i] = ev[(f.pos + i) % k];
      }
      cyc.swap(nc); eu.swap(neu); ev.swap(nev);
      f.stage = 4;
      continue;
    }

    if (f.stage == 4) {                  // tail descent through the new base sub-node
      Size child = cyc[0];
      f.stage = 5;
      if (child >= mN && todo.size() <= mNumNodes + 16) {
        RotateFrame childFrame(child, f.newBase);
        childFrame.basePath.swap(f.basePath);
        if (!childFrame.basePath.empty() && childFrame.basePath.back() == child)
          childFrame.basePath.pop_back();
        todo.push_back(std::move(childFrame));
        ++mCounters.blossomRotationFrames;
      }
      continue;
    }

    mBBase[f.B] = mBBase[cyc[0]];        // all child rotations are now complete
    todo.pop_back();
  }
}

// Flip matched<->unmatched along the tree path from EVEN outer node startEven up
// to its free root. Invariant: climbing from EVEN E, parent is ODD O via a
// MATCHED edge; O's parent is EVEN E2 via the UNMATCHED grow edge stored on O.
inline void ScalingMatching::flipUpToRoot(Size startEven) {
  CRUMB("flipUpToRoot", startEven, 0);
  Size flipSteps = 0;
  if (!mIncrShrink) {
    // ---- ORIGINAL (relabel-based) walk: EVEN E -> ODD O -> EVEN E2, read O's parent edge.
    // Untouched so the proven baseline (relabel keeps E outer + the EVEN-ODD-EVEN invariant) is
    // BYTE-IDENTICAL. ----
    Size E = startEven;
    while (mNParent[E] != scNull()) {
      if (++flipSteps > mNumNodes + 8) return;
      Size O  = mNParent[E];
      Size E2 = mNParent[O];
      Size uo = mNParEdgeU[O];      // endpoint in O
      Size ue = mNParEdgeV[O];      // endpoint in E2
      if (O  >= mN) rotateToBase(O,  uo);
      if (E2 >= mN) rotateToBase(E2, ue);
      mMate[uo] = ue; mMate[ue] = uo;
      E = E2;
    }
    return;
  }
  // ---- S5's EXACT relabel-independent walk (ScalingS5Engine::flipUpToRoot :518), now that
  // lcaOuter + shrink use the SAME treeParent resolution so the forest is S5-consistent. ----
  // treeParent(X) = outerOf(mNParent[X]) resolved through absorbed nodes, NIL if it resolves to X
  // itself (absorbed into the same blossom) -- the contracted-forest parent. Flip O's grow edge.
  Size E = outerOf(startEven);
  // S5 disjointness marking, done HERE (O(path), not O(n)): mark each PATH node's leaves used as
  // we walk. (The old O(n)-per-augment whole-node-scan in augmentFrom was 74% of runtime = the
  // super-linear blowup.) Non-path vertices of the consumed tree are skipped via treeLive (dead
  // tree). markUsed of a node's leaves is O(|node|); summed over the path = O(path size).
  bool mark = mIncrShrink && !mUsedV.empty();
#ifdef SCALING_LEAF_NOVEC
  std::vector<Size> markLeaves;
  if (mark) { scalingLfCollect(E, &markLeaves); for (Size j=0;j<markLeaves.size();++j) if (markLeaves[j]<mUsedV.size()) mUsedV[markLeaves[j]]=1; }
#else
  if (mark) { std::vector<Size>& lv = mNodeLeaves[E]; for (Size j=0;j<lv.size();++j) if (lv[j]<mUsedV.size()) mUsedV[lv[j]]=1; }
#endif
  for (;;) {
    if (++flipSteps > 2 * mNumNodes + 16) return;
    Size O = treeParentR(E);
    if (O == scNull()) break;
    Size E2 = treeParentR(O);
    Size uo = mNParEdgeU[O];
    Size ue = mNParEdgeV[O];
    if (uo == scNull() || ue == scNull()) break;
    if (O  >= mN) rotateToBase(O,  uo);
    if (E2 != scNull() && E2 >= mN) rotateToBase(E2, ue);
    mMate[uo] = ue; mMate[ue] = uo;
    if (mark) {   // mark O's and E2's leaves (the path continues through them)
#ifdef SCALING_LEAF_NOVEC
      scalingLfCollect(O, &markLeaves); for (Size j=0;j<markLeaves.size();++j) if (markLeaves[j]<mUsedV.size()) mUsedV[markLeaves[j]]=1;
      if (E2 != scNull()) { scalingLfCollect(E2, &markLeaves); for (Size j=0;j<markLeaves.size();++j) if (markLeaves[j]<mUsedV.size()) mUsedV[markLeaves[j]]=1; }
#else
      std::vector<Size>& lo = mNodeLeaves[O]; for (Size j=0;j<lo.size();++j) if (lo[j]<mUsedV.size()) mUsedV[lo[j]]=1;
      if (E2 != scNull()) { std::vector<Size>& le = mNodeLeaves[E2]; for (Size j=0;j<le.size();++j) if (le[j]<mUsedV.size()) mUsedV[le[j]]=1; }
#endif
    }
    if (E2 == scNull()) break;
    E = E2;
  }
}

#ifdef PS5_LEAF_SPLICE
inline void ScalingMatching::scalingLfInit(Size vertexCount, Size nodeCapacity) {
  mLfNextX.assign(vertexCount, scNull());
  mLfPrevX.assign(vertexCount, scNull());
  mLfHeadX.assign(nodeCapacity, scNull());
  mLfTailX.assign(nodeCapacity, scNull());
  mLfCntX.assign(nodeCapacity, 0);
  mOldLfOrderX.clear();
  mOldLfBeginX.assign(nodeCapacity, scNull());
  mOldLfEndX.assign(nodeCapacity, scNull());
  for (Size v = 0; v < vertexCount && v < nodeCapacity; ++v) {
    mLfHeadX[v] = v; mLfTailX[v] = v; mLfCntX[v] = 1;
  }
}

inline void ScalingMatching::scalingLfSnapshotOld(Size oldBegin, Size oldEnd) {
  Size cap = mNodeLeaves.size();
  mOldLfOrderX.clear(); mOldLfOrderX.reserve(mN);
  mOldLfBeginX.assign(cap, scNull());
  mOldLfEndX.assign(cap, scNull());
  std::vector<char> seenNode(cap, 0), seenVertex(mN, 0);
  struct Frame { Size node, next; bool entered; };
  std::vector<Size> roots;
  // Omega' consists of positive-z records. A live blossom may retain zero-z
  // historical ancestors, so mBParent==null is not the right root test: choose
  // exactly the positive blossoms with no positive inherited ancestor.
  for (Size b = oldBegin; b < oldEnd && b < cap; ++b) {
    if (b >= mZ.size() || mZ[b] <= 0) continue;
    bool positiveAncestor = false;
    Size p = b < mBParent.size() ? mBParent[b] : scNull(), guard = 0;
    while (p != scNull() && p < oldEnd && ++guard <= oldEnd + 2) {
      if (p >= oldBegin && p < mZ.size() && mZ[p] > 0) {
        positiveAncestor = true; break;
      }
      p = p < mBParent.size() ? mBParent[p] : scNull();
    }
    if (!positiveAncestor) roots.push_back(b);
  }
  for (Size ri = 0; ri < roots.size(); ++ri) {
    Size root = roots[ri];
    if (root >= cap || seenNode[root]) continue;
    std::vector<Frame> stack;
    stack.push_back(Frame{root, 0, false});
    while (!stack.empty()) {
      Frame& f = stack.back();
      if (!f.entered) {
        if (f.node >= cap || seenNode[f.node]) { stack.pop_back(); continue; }
        seenNode[f.node] = 1; f.entered = true;
        mOldLfBeginX[f.node] = mOldLfOrderX.size();
        if (f.node < mN) {
          if (!seenVertex[f.node]) {
            seenVertex[f.node] = 1; mOldLfOrderX.push_back(f.node);
          }
          mOldLfEndX[f.node] = mOldLfOrderX.size();
          stack.pop_back();
          continue;
        }
      }
      if (f.node < mBCycle.size() && f.next < mBCycle[f.node].size()) {
        Size child = mBCycle[f.node][f.next++];
        if (child < cap && !seenNode[child])
          stack.push_back(Frame{child, 0, false});
      } else {
        mOldLfEndX[f.node] = mOldLfOrderX.size();
        stack.pop_back();
      }
    }
  }
#ifdef SCALING_LEAF_SHADOW
  for (Size b = oldBegin; b < oldEnd && b < cap; ++b) {
    if (b >= mOldLfBeginX.size() || mOldLfBeginX[b] == scNull()) continue;
    std::vector<Size> snapshot(
      mOldLfOrderX.begin() + mOldLfBeginX[b],
      mOldLfOrderX.begin() + mOldLfEndX[b]);
    std::vector<Size> expected = mNodeLeaves[b];
    std::sort(snapshot.begin(), snapshot.end());
    std::sort(expected.begin(), expected.end());
    expected.erase(std::unique(expected.begin(), expected.end()), expected.end());
    if (snapshot != expected) {
      std::fprintf(stderr,
        "[SCALING-LEAF-SHADOW old-snapshot B=%lu snapshot=%lu vector=%lu]\n",
        (unsigned long)b, (unsigned long)snapshot.size(),
        (unsigned long)expected.size());
      std::abort();
    }
  }
#endif
}

inline void ScalingMatching::scalingLfEnsureFromVectors() {
  if (!mLfHeadX.empty()) return;
  Size vertexCount = !mUf.empty() ? mUf.size() : mN;
  Size cap = mNodeLeaves.size();
  scalingLfInit(vertexCount, cap);
  std::vector<char> seen(vertexCount, 0);
  for (Size b = 0; b < mNumNodes; ++b) {
    if (b < mBParent.size() && mBParent[b] != scNull()) continue;
    const std::vector<Size>& leaves = mNodeLeaves[b];
    if (leaves.empty()) continue;
    Size head = scNull(), tail = scNull(), count = 0;
    for (Size i = 0; i < leaves.size(); ++i) {
      Size v = leaves[i];
      if (v >= vertexCount || seen[v]) continue;
      seen[v] = 1;
      if (head == scNull()) head = tail = v;
      else { mLfNextX[tail] = v; mLfPrevX[v] = tail; tail = v; }
      ++count;
    }
    mLfHeadX[b] = head; mLfTailX[b] = tail; mLfCntX[b] = count;
  }
  // Recover the segment endpoints of interior blossoms from their top-level
  // owner's list. This path is only for legacy/test entry points; runScaling
  // initializes singleton lists before any blossom forms.
  for (Size b = mN; b < mNumNodes; ++b) {
    if (b >= mBParent.size() || mBParent[b] == scNull() ||
        mNodeLeaves[b].empty()) continue;
    std::vector<char> mine(vertexCount, 0);
    for (Size i = 0; i < mNodeLeaves[b].size(); ++i)
      if (mNodeLeaves[b][i] < vertexCount) mine[mNodeLeaves[b][i]] = 1;
    Size top = b, guard = 0;
    while (top < mBParent.size() && mBParent[top] != scNull() &&
           ++guard <= mNumNodes + 2) top = mBParent[top];
    Size head = scNull(), tail = scNull(), count = 0;
    if (top < mLfHeadX.size()) {
      Size end = mLfTailX[top], steps = 0;
      for (Size v = mLfHeadX[top]; v != scNull() && steps++ <= vertexCount;
           v = (v == end ? scNull() : mLfNextX[v]))
        if (mine[v]) { if (head == scNull()) head = v; tail = v; ++count; }
    }
    mLfHeadX[b] = head; mLfTailX[b] = tail; mLfCntX[b] = count;
  }
}

inline void ScalingMatching::scalingLfCheckNode(Size B,
                                                 const char* where) const {
#ifdef SCALING_LEAF_SHADOW
  if (B >= mLfHeadX.size() || B >= mNodeLeaves.size()) return;
  std::vector<Size> listLeaves;
  Size end = mLfTailX[B], guard = 0;
  for (Size v = mLfHeadX[B]; v != scNull() && guard++ <= mLfNextX.size();
       v = (v == end ? scNull() : mLfNextX[v])) listLeaves.push_back(v);
  std::vector<Size> vectorLeaves = mNodeLeaves[B];
  std::sort(listLeaves.begin(), listLeaves.end());
  std::sort(vectorLeaves.begin(), vectorLeaves.end());
  vectorLeaves.erase(std::unique(vectorLeaves.begin(), vectorLeaves.end()),
                     vectorLeaves.end());
  if (listLeaves != vectorLeaves || listLeaves.size() != mLfCntX[B]) {
    std::fprintf(stderr,
      "[SCALING-LEAF-SHADOW %s B=%lu list=%lu vector=%lu count=%lu active=%d parent=%ld z=%lld oldSep=%d old=[%lu,%lu)]\n",
      where, (unsigned long)B, (unsigned long)listLeaves.size(),
      (unsigned long)vectorLeaves.size(), (unsigned long)mLfCntX[B],
      (B < mBActive.size() && mBActive[B]) ? 1 : 0,
      (long)(B < mBParent.size() && mBParent[B] != scNull() ? mBParent[B] : -1),
      (long long)(B < mZ.size() ? mZ[B] : 0), mOldForestSeparated ? 1 : 0,
      (unsigned long)mAuditOldBegin, (unsigned long)mAuditOldEnd);
    std::abort();
  }
#else
  (void)B; (void)where;
#endif
}

inline Size ScalingMatching::scalingLfCount(Size B) const {
  if (B < mLfCntX.size() && mLfHeadX[B] != scNull()) return mLfCntX[B];
  return B < mNodeLeaves.size() ? mNodeLeaves[B].size() : 0;
}

// Return B's direct laminar child containing original vertex v.  While B is
// active, mBParent is precisely this ancestry relation, so walking ancestors
// is proportional to nesting depth rather than to the sizes of all preceding
// sibling leaf ranges.  relinkOddCycle calls this before dissolve(B), while
// the complete parent chain is still intact.
inline Size ScalingMatching::directChildContaining(
    Size B, Size v, unsigned long long* steps) const {
  if (B < mN || B >= mNumNodes || v >= mN) return scNull();
  Size x = v;
  Size guard = 0;
  while (x != B && guard++ <= mNumNodes) {
    if (steps) ++*steps;
    if (x >= mBParent.size()) return scNull();
    Size p = mBParent[x];
    if (p == B) return x;
    if (p == scNull() || p == x) return scNull();
    x = p;
  }
  return scNull();
}

inline void ScalingMatching::scalingLfCollect(
    Size B, std::vector<Size>* out) const {
  if (!out) return;
  out->clear();
  if (B >= mLfHeadX.size() || mLfHeadX[B] == scNull()) {
    if (B < mNodeLeaves.size()) *out = mNodeLeaves[B];
    return;
  }
  out->reserve(mLfCntX[B]);
  Size end = mLfTailX[B], guard = 0;
  for (Size v = mLfHeadX[B]; v != scNull() && guard++ <= mLfNextX.size();
       v = (v == end ? scNull() : mLfNextX[v])) out->push_back(v);
}

inline void ScalingMatching::scalingLfOnShrink(Size B) {
  scalingLfEnsureFromVectors();
  if (B >= mLfHeadX.size()) {
    Size cap = mZ.size();
    mLfHeadX.resize(cap, scNull()); mLfTailX.resize(cap, scNull());
    mLfCntX.resize(cap, 0);
  }
  Size head = scNull(), tail = scNull(), count = 0;
  for (Size i = 0; i < mBCycle[B].size(); ++i) {
    Size child = mBCycle[B][i];
    if (child >= mLfHeadX.size() || mLfHeadX[child] == scNull()) continue;
    if (head == scNull()) head = mLfHeadX[child];
    else {
      mLfNextX[tail] = mLfHeadX[child];
      mLfPrevX[mLfHeadX[child]] = tail;
    }
    tail = mLfTailX[child]; count += mLfCntX[child];
  }
  mLfHeadX[B] = head; mLfTailX[B] = tail; mLfCntX[B] = count;
  scalingLfCheckNode(B, "shrink");
#ifdef SCALING_LEAF_SHADOW
  for (Size b = 0; b < mNumNodes; ++b)
    if (b != B && b < mNodeLeaves.size() && !mNodeLeaves[b].empty() &&
        (b < mN || b >= mBActive.size() || mBActive[b]))
      if (!(mOldForestSeparated && b >= mAuditOldBegin && b < mAuditOldEnd))
      if (!(b >= mAuditOldBegin && b < mAuditOldEnd &&
            b < mZ.size() && mZ[b] == 0))
      scalingLfCheckNode(b, "shrink-sweep");
#endif
}

inline void ScalingMatching::scalingLfOnDissolve(Size B) {
  if (mLfHeadX.empty() || B >= mLfHeadX.size()) return;
  for (Size i = 0; i < mBCycle[B].size(); ++i) {
    Size child = mBCycle[B][i];
    if (child >= mLfHeadX.size() || mLfHeadX[child] == scNull()) continue;
    mLfPrevX[mLfHeadX[child]] = scNull();
    mLfNextX[mLfTailX[child]] = scNull();
    scalingLfCheckNode(child, "dissolve-child");
  }
  mLfHeadX[B] = scNull(); mLfTailX[B] = scNull(); mLfCntX[B] = 0;
#ifdef SCALING_LEAF_SHADOW
  for (Size b = 0; b < mNumNodes; ++b)
    if (b != B && b < mNodeLeaves.size() && !mNodeLeaves[b].empty() &&
        (b < mN || b >= mBActive.size() || mBActive[b]))
      if (!(mOldForestSeparated && b >= mAuditOldBegin && b < mAuditOldEnd))
      if (!(b >= mAuditOldBegin && b < mAuditOldEnd &&
            b < mZ.size() && mZ[b] == 0))
      scalingLfCheckNode(b, "dissolve-sweep");
#endif
}
#endif

inline Size ScalingMatching::nodeLeafCount(Size B) const {
#ifdef SCALING_LEAF_NOVEC
  if (mOldForestSeparated && B >= mAuditOldBegin && B < mAuditOldEnd &&
      B < mOldLfBeginX.size() && mOldLfBeginX[B] != scNull())
    return mOldLfEndX[B] - mOldLfBeginX[B];
  return scalingLfCount(B);
#else
  return B < mNodeLeaves.size() ? mNodeLeaves[B].size() : 0;
#endif
}

inline void ScalingMatching::collectNodeLeaves(
    Size B, std::vector<Size>* out) const {
#ifdef SCALING_LEAF_NOVEC
  if (out && mOldForestSeparated && B >= mAuditOldBegin && B < mAuditOldEnd &&
      B < mOldLfBeginX.size() && mOldLfBeginX[B] != scNull()) {
    out->assign(mOldLfOrderX.begin() + mOldLfBeginX[B],
                mOldLfOrderX.begin() + mOldLfEndX[B]);
    return;
  }
  scalingLfCollect(B, out);
#else
  if (!out) return;
  out->clear();
  if (B < mNodeLeaves.size()) *out = mNodeLeaves[B];
#endif
}

// SHRINK: EVEN outer nodes E,W (same tree) joined by eligible edge (uOrig in E,
// vOrig in W). Build a new blossom B over the odd cycle E..lca..W; z(B)=0.
inline void ScalingMatching::shrink(Size E, Size W, Size uOrig, Size vOrig) {
  CRUMB("shrink", E, W);
  ++mCounters.blossomFormations;
  Size L = lcaOuter(E, W);
#ifdef SCDBG
  std::fprintf(stderr,"    SHRINK E=%lu W=%lu lca=%ld via %lu-%lu\n",(unsigned long)E,(unsigned long)W,(long)(L==scNull()?-1:(long)L),(unsigned long)uOrig,(unsigned long)vOrig); std::fflush(stderr);
#endif
  if (L == scNull()) return;   // not same tree: caller bug; bail (no shrink)
#ifdef OPCNT
  ++gOpShrinks;
#endif
  Size B = mNumNodes++;
  // Blossom ids are allocated monotonically (not reclaimed on dissolve), so a
  // search with many shrink/dissolve cycles can exceed the initial capacity.
  // Grow all per-node arrays on demand. (Phase 5 will reclaim ids via free-list.)
  if (B >= mZ.size()) {
    Size ncap = 2 * B + 16;
    mZ.resize(ncap, 0);
    mBParent.resize(ncap, scNull());
    mBBase.resize(ncap, scNull());
    mBCycle.resize(ncap);
    mBCycEdgeU.resize(ncap);
    mBCycEdgeV.resize(ncap);
    mNodeLeaves.resize(ncap);
    mBActive.resize(ncap, 0);
    mNLabel.resize(ncap, LBL_FREE_OF_TREE);
    mNParent.resize(ncap, scNull());
    mNParEdgeU.resize(ncap, scNull());
    mNParEdgeV.resize(ncap, scNull());
    mNRoot.resize(ncap, scNull());
  }
  if (B >= mBActive.size()) mBActive.resize(B + 1, 0);
  mBActive[B] = 1;
  mZ[B] = 0; mBParent[B] = scNull();
  mBCycle[B].clear(); mBCycEdgeU[B].clear(); mBCycEdgeV[B].clear();
  // reused member buffers (no per-shrink malloc); shrink is non-reentrant so sharing is safe.
  std::vector<Size>& eSide = mShrEside; std::vector<Size>& eU = mShrEu; std::vector<Size>& eV = mShrEv;
  std::vector<Size>& wSide = mShrWside; std::vector<Size>& wU = mShrWu; std::vector<Size>& wV = mShrWv;
  eSide.clear(); eU.clear(); eV.clear(); wSide.clear(); wU.clear(); wV.clear();
#ifdef CYCDBG
  { Size x = E; Size sc = 0; while (x != L) { if (++sc > mNumNodes + 8) { std::fprintf(stderr, "[CYCDBG] shrink E-side mNParent never reaches L=%lu from E=%lu (chain corrupt)\n", (unsigned long)L, (unsigned long)E); std::fflush(stderr); break; } eSide.push_back(x); eU.push_back(mNParEdgeU[x]); eV.push_back(mNParEdgeV[x]); x = mNParent[x]; } }
  { Size x = W; Size sc = 0; while (x != L) { if (++sc > mNumNodes + 8) { std::fprintf(stderr, "[CYCDBG] shrink W-side mNParent never reaches L=%lu from W=%lu (chain corrupt)\n", (unsigned long)L, (unsigned long)W); std::fflush(stderr); break; } wSide.push_back(x); wU.push_back(mNParEdgeU[x]); wV.push_back(mNParEdgeV[x]); x = mNParent[x]; } }
#else
  // CYCLE WALK to L. Steps must use the SAME resolution as lcaOuter (above) or the walk misses L.
  // Under mIncrShrink: treeParent step = outerOf(mNParent[x]) (S5 shrink :548); else raw mNParent
  // (byte-identical non-incr path). Edges pushed are mNParEdge*[x] of the node x being left -- in
  // sync with x, exactly as S5 (push x's edge, THEN step).
  if (mIncrShrink) {
    LL cg = 0;
    { Size x = E; while (x != L) { if (++cg > (LL)mNumNodes + 8) break; eSide.push_back(x); eU.push_back(mNParEdgeU[x]); eV.push_back(mNParEdgeV[x]); Size p = mNParent[x]; if (p == scNull()) break; Size nx = outerOf(p); if (nx == x) break; x = nx; } }
    cg = 0;
    { Size x = W; while (x != L) { if (++cg > (LL)mNumNodes + 8) break; wSide.push_back(x); wU.push_back(mNParEdgeU[x]); wV.push_back(mNParEdgeV[x]); Size p = mNParent[x]; if (p == scNull()) break; Size nx = outerOf(p); if (nx == x) break; x = nx; } }
  } else {
  { Size x = E; while (x != L) { eSide.push_back(x); eU.push_back(mNParEdgeU[x]); eV.push_back(mNParEdgeV[x]); x = mNParent[x]; } }
  { Size x = W; while (x != L) { wSide.push_back(x); wU.push_back(mNParEdgeU[x]); wV.push_back(mNParEdgeV[x]); x = mNParent[x]; } }
  }
#endif
  mBCycle[B].push_back(L);
  for (Size i = eSide.size(); i-- > 0; ) {
    mBCycEdgeU[B].push_back(eV[i]); mBCycEdgeV[B].push_back(eU[i]);
    mBCycle[B].push_back(eSide[i]);
  }
  mBCycEdgeU[B].push_back(uOrig); mBCycEdgeV[B].push_back(vOrig);
  for (Size i = 0; i < wSide.size(); ++i) {
    if (i == 0) { mBCycle[B].push_back(wSide[i]); }
    else { mBCycEdgeU[B].push_back(wU[i-1]); mBCycEdgeV[B].push_back(wV[i-1]); mBCycle[B].push_back(wSide[i]); }
  }
  if (!wSide.empty()) {
    Size last = wSide.size()-1;
    mBCycEdgeU[B].push_back(wU[last]);
    mBCycEdgeV[B].push_back(wV[last]);
  }
  mBBase[B] = mBBase[L];
  for (Size i = 0; i < mBCycle[B].size(); ++i) mBParent[mBCycle[B][i]] = B;
  // B's leaf set = union of its children's leaf sets (children keep their own,
  // so dissolve restores them for free).
  mNodeLeaves[B].clear();
#ifndef SCALING_LEAF_NOVEC
  Size blossomLeafCount = 0;
  for (Size i = 0; i < mBCycle[B].size(); ++i)
    blossomLeafCount += mNodeLeaves[mBCycle[B][i]].size();
  mNodeLeaves[B].reserve(blossomLeafCount);
  for (Size i = 0; i < mBCycle[B].size(); ++i) {
    std::vector<Size>& cl = mNodeLeaves[mBCycle[B][i]];
    for (Size j = 0; j < cl.size(); ++j) mNodeLeaves[B].push_back(cl[j]);
  }
#endif
#ifdef PS5_LEAF_SPLICE
  scalingLfOnShrink(B);
#endif
  // UNION-FIND: every child already has one representative for all of its
  // leaves. Merge those child representatives into B's base representative;
  // revisiting every leaf here made a nested chain quadratic even though the
  // union-find invariant is identical. This is the standard outer-blossom
  // union operation from Section 5.1.
  if (!mUf.empty() && mBBase[B] < mUf.size()) {
    Size baseRep = ufFind(mBBase[B]);
    for (Size i = 0; i < mBCycle[B].size(); ++i) {
      Size child = mBCycle[B][i];
      if (child >= mBBase.size() || mBBase[child] >= mUf.size()) continue;
      Size childRep = ufFind(mBBase[child]);
      if (childRep != baseRep) mUf[childRep] = baseRep;
    }
    mUf[baseRep] = baseRep;
    // base-vertex -> outer-node map: B is now the outer node whose uf-rep is baseRep.
    if (!mBaseToNode.empty() && baseRep < mBaseToNode.size()) mBaseToNode[baseRep] = B;
  }
  // B inherits L's tree role (EVEN, same parent/parent-edge/root as L). Resolve
  // L's parent to its current OUTER node: in the persistent forest L's stored
  // mNParent may be a node that just got absorbed into B (or B itself), which
  // would make B its own ancestor (the n=20 self-loop mNParent[20]=20). If the
  // resolved parent is B (or scNull), B is a ROOT blossom.
  mNLabel[B] = LBL_EVEN;
  rememberCompactEventLabel(B);
  if (mIncrRelabel ||
      (mShellActive && !mHybridCanonical &&
       std::getenv("MBX_SHELL_LOCAL_RESET")))
    mLabeledList.push_back(B);   // record for incremental/shell-local entry relabel
  Size lp = mNParent[L];
  Size lpOuter = (lp == scNull()) ? scNull() : outerOf(lp);
  if (lpOuter == B) lpOuter = scNull();      // parent absorbed into B -> B is root
  mNParent[B] = lpOuter;
  mNParEdgeU[B] = mNParEdgeU[L]; mNParEdgeV[B] = mNParEdgeV[L];
  mNRoot[B] = (lpOuter == scNull()) ? B : mNRoot[L];
  if (mEventStructuresActive
#ifdef EVENT_SFM_REUSE
      || mEventSfmInitialized
#endif
      ) eventSfmOnShrink(B);
#ifdef EVENT_S5_TARGET_KEY_SHADOW
#ifdef EVENT_S5_TARGET_KEY_PERSIST_FRESH
  if (mS5TargetKeyTopologyReady && mS5TargetKeyTopologyScale == mAuditScale)
    s5TargetKeyOnShrink(B);
#else
  if (mS5TargetKeyActive) s5TargetKeyOnShrink(B);
#endif
#endif
}

// EXPERIMENT (gated, default OFF): re-thread the alternating tree through ODD blossom B's cycle so a
// live-ODD dissolution does NOT need a full relabel. Edmonds odd-blossom expansion: B entered the
// tree at its base sub-node (matched up to B's parent). One child sub-tree hangs off some cycle
// sub-node via an unmatched edge. Expanding B: the even-length arc from base to the child-attach
// sub-node alternates ODD(base),EVEN,ODD,...,EVEN(attach), each linked by the cycle's matched edges;
// the other arc's sub-nodes pair up matched-internally and leave the tree (become FREE). We rewire
// mNLabel/mNParent/mNParEdge/mNRoot for the on-arc sub-nodes and FREE the off-arc ones. Returns
// false (caller falls back to full relabel) if the structure isn't the clean expandable case.
inline bool ScalingMatching::relinkOddCycle(Size B) {
  ScalingNanosecondGuard relinkTimer(&mCounters.dissolutionRelinkNanoseconds);
#ifdef SCPROF
  extern long long gRelinkB[8];
#endif
  const std::vector<Size>& cyc = mBCycle[B];
  Size k = cyc.size();
#ifdef SCPROF
  if (k < 3 || (k % 2) == 0) { ++gRelinkB[0]; return false; }
#endif
  if (k < 3 || (k % 2) == 0) return false;   // a real blossom cycle is odd-length >=3
  // Nested members remain contracted while B itself is expanded one laminar
  // level. A nested zero-z member exposed by dissolve(B) receives its own
  // dissolution event below; positive-z members remain ordinary tree nodes.
  Size root = mNRoot[B];
  Size rawPar = mNParent[B];
  Size par = rawPar == scNull() ? scNull() : outerOf(rawPar);
  if (par == B) par = scNull();
  Size peU = mNParEdgeU[B], peV = mNParEdgeV[B];   // B's tree edge: peU in B, peV in parent
  if (par == scNull()) return false; // an INNER blossom must have an entry arc
  // base sub-node = cyc[0] (shrink stores L=base first). The base carries B's upward match/parent.
  // Find the cycle index whose sub-node contains peU (B's parent-edge endpoint inside B) = entry.
  Size entryIdx = scNull();
  Size entryChild = directChildContaining(
      B, peU, &mCounters.dissolutionContainmentSteps);
  for (Size i = 0; i < k; ++i) {
    if (cyc[i] >= mNumNodes) return false;
    if (cyc[i] == entryChild) { entryIdx = i; break; }
  }
#ifdef SCPROF
  if (par != scNull() && entryIdx == scNull()) { ++gRelinkB[1]; return false; }
#endif
  if (par != scNull() && entryIdx == scNull()) return false;   // can't locate entry -> bail safely
  // Find the child sub-tree: a node C with parent resolving to B, attached at vertex aU in B.
  // Scan is O(live nodes) -- acceptable vs the O(m) relabel it replaces. Collect ALL such children.
  // The child-attach cycle index = the sub-node containing aU.
  Size childAttachIdx = entryIdx;   // default: child hangs at entry (degenerate, even arc len 0)
  bool haveChild = false;
  Size hangingChild = scNull();
  Size baseVtx = mBBase[B];
  Size childMate = baseVtx < mMate.size() ? mMate[baseVtx] : scNull();
  if (childMate != scNull() && outerOf(childMate) != B) {
    Size child = outerOf(childMate);
    if (child != scNull() && child != B && mNParent[child] != scNull() &&
        outerOf(mNParent[child]) == B) {
      Size aU = mNParEdgeV[child], attach = scNull();
      Size attachChild = directChildContaining(
          B, aU, &mCounters.dissolutionContainmentSteps);
      for (Size i = 0; i < k; ++i)
        if (cyc[i] == attachChild) { attach = i; break; }
      if (attach == scNull()) return false;
      childAttachIdx = attach; hangingChild = child; haveChild = true;
    }
  }
  (void)haveChild;
  // Walk the cycle both directions from entryIdx to childAttachIdx; the EVEN-edge-count arc is the
  // tree path. (k is odd, so exactly one arc has even length.) dist forward = (childAttachIdx-entry) mod k.
  Size fwd = (childAttachIdx + k - entryIdx) % k;
  Size dir; Size arcLen;
  if ((fwd % 2) == 0) { dir = 1; arcLen = fwd; }
  else { dir = k - 1; arcLen = k - fwd; }   // backward step = +(k-1) mod k; backward arc len = k-fwd
#ifdef RELINKDBG
  { static int rc=0; if(++rc<=6) {
      // is cyc[0] the base sub-node? does peU live in base?
      bool peuInBase=false; { const std::vector<Size>& lv=mNodeLeaves[cyc[0]]; for(Size j=0;j<lv.size();++j) if(lv[j]==peU){peuInBase=true;break;} }
      bool baseSnIs0=false; { const std::vector<Size>& lv=mNodeLeaves[cyc[0]]; for(Size j=0;j<lv.size();++j) if(lv[j]==mBBase[B]){baseSnIs0=true;break;} }
      std::fprintf(stderr,"[RL B=%lu k=%lu entry=%lu attach=%lu fwd=%lu dir=%lu arcLen=%lu par=%lu peU=%lu base=%lu peuInBase=%d baseAt0=%d mate[peU]=%lu mate[base]=%lu]\n",
        (unsigned long)B,(unsigned long)k,(unsigned long)entryIdx,(unsigned long)childAttachIdx,(unsigned long)fwd,(unsigned long)dir,(unsigned long)arcLen,(unsigned long)par,(unsigned long)peU,(unsigned long)mBBase[B],(int)peuInBase,(int)baseSnIs0,(unsigned long)mMate[peU],(unsigned long)mMate[mBBase[B]]); }
  }
#endif
  // Relabel on-arc sub-nodes alternating ODD(entry),EVEN,ODD,...; link each to the previous via the
  // cycle matched edge. Off-arc sub-nodes become FREE (their internal match is restored by dissolve).
  // First mark which indices are on the arc.
  std::vector<char> onArc(k, 0);
  { Size i = entryIdx; for (Size s = 0; s <= arcLen; ++s) { onArc[i] = 1; i = (i + dir) % k; } }
  // entry sub-node: ODD, parent = B's parent (the tree above), parent-edge = B's parent edge.
  Size prev = scNull();
  { Size i = entryIdx;
    for (Size s = 0; s <= arcLen; ++s) {
      Size sn = cyc[i];
      mNLabel[sn] = (s % 2 == 0) ? LBL_ODD : LBL_EVEN;
      rememberCompactEventLabel(sn);
      mNRoot[sn] = root;
      if (s == 0) { mNParent[sn] = par; mNParEdgeU[sn] = peU; mNParEdgeV[sn] = peV; }
      else {
        // linked to prev via the cycle edge between (i-dir) and i. mBCycEdge stores edge between
        // consecutive cycle members; index convention: edge e[j] connects cyc[j] and cyc[j+1].
        Size eIdx = (dir == 1) ? ((i + k - 1) % k) : i;   // edge between prev and sn
        Size eu = mBCycEdgeU[B][eIdx], ev = mBCycEdgeV[B][eIdx];
        // orient: parent endpoint in prev, child endpoint in sn.
        mNParent[sn] = prev;
        if (dir == 1) { mNParEdgeU[sn] = ev; mNParEdgeV[sn] = eu; }
        else          { mNParEdgeU[sn] = eu; mNParEdgeV[sn] = ev; }
#ifdef RELINKDBG
        { static int ec=0; if(++ec<=24) { bool m = (mMate[eu]==ev); std::fprintf(stderr,"   edge s=%lu sn=%lu(lbl%d) prev=%lu eu=%lu ev=%lu matched=%d  (EVEN-needs-unmatched-to-ODDparent)\n",(unsigned long)s,(unsigned long)sn,(int)mNLabel[sn],(unsigned long)prev,(unsigned long)eu,(unsigned long)ev,(int)m); } }
#endif
      }
      if (sn >= mN) {
        Size requiredBase = mBBase[sn];
        if (mNLabel[sn] == LBL_EVEN) {
          // The edge to the preceding ODD member is matched.
          requiredBase = mNParEdgeU[sn];
        } else if (s < arcLen) {
          // This ODD member is matched to the following EVEN member.
          Size next = (i + dir) % k;
          Size eIdx = (dir == 1) ? i : next;
          requiredBase = (dir == 1) ? mBCycEdgeU[B][eIdx]
                                    : mBCycEdgeV[B][eIdx];
        } else if (haveChild) {
          // The final ODD member is matched to B's hanging EVEN subtree.
          requiredBase = mNParEdgeV[hangingChild];
        }
        rotateToBase(sn, requiredBase);
      }
      prev = sn; i = (i + dir) % k;
    }
  }
  // Off-arc sub-nodes: FREE (not in tree). dissolve() restores their internal matching.
  for (Size i = 0; i < k; ++i) if (!onArc[i]) {
    Size sn = cyc[i];
    mNLabel[sn] = LBL_FREE_OF_TREE; mNParent[sn] = scNull();
    mNParEdgeU[sn] = scNull(); mNParEdgeV[sn] = scNull(); mNRoot[sn] = scNull();
  }
  // The matched EVEN subtree below B used to point at the contracted node B.
  // Expansion removes B, so splice that subtree onto the ODD cycle member
  // containing B's base/matched endpoint. Its stored matched parent edge stays
  // valid; only the contracted-tree parent and root ownership change.
  if (haveChild) {
    Size attachNode = cyc[childAttachIdx];
    if (mNLabel[attachNode] != LBL_ODD ||
        mNLabel[hangingChild] != LBL_EVEN)
      return false;
    mNParent[hangingChild] = attachNode;
    mNRoot[hangingChild] = root;
  }
#ifdef RELINKDBG
  // Check the MATCHING of every cycle-base vertex: each sub-node's base must be matched (to its
  // tree parent if EVEN/in-arc-odd-base, or internally if off-arc). A vertex with mMate==null after
  // relink = the collapse bug. Print any free base among cycle members.
  { static int mc=0; if (++mc<=8) for (Size i = 0; i < k; ++i) { Size bs = mBBase[cyc[i]];
      std::fprintf(stderr,"   [match] cyc[%lu]=%lu base=%lu mate[base]=%lu onArc=%d lbl=%d\n",
        (unsigned long)i,(unsigned long)cyc[i],(unsigned long)bs,(unsigned long)mMate[bs],(int)onArc[i],(int)mNLabel[cyc[i]]); } }
#endif
  return true;
}

// DISSOLVE inner root blossom B with z==0: free its sub-nodes (they become outer).
inline void ScalingMatching::dissolve(Size B) {
  ScalingNanosecondGuard dissolutionTimer(&mCounters.dissolutionNanoseconds);
  ++mCounters.blossomDissolutions;
  if (mEventStructuresActive
#ifdef EVENT_SFM_REUSE
      || mEventSfmInitialized
#endif
      ) eventSfmOnDissolve(B);
#ifdef EVENT_S5_TARGET_KEY_SHADOW
#ifdef EVENT_S5_TARGET_KEY_PERSIST_FRESH
  if (mS5TargetKeyTopologyReady && mS5TargetKeyTopologyScale == mAuditScale)
    s5TargetKeyOnDissolve(B);
#else
  if (mS5TargetKeyActive) s5TargetKeyOnDissolve(B);
#endif
#endif
#ifdef OPCNT
  ++gOpDissolves;
#ifdef DISSPLIT
  { extern long long gDisIn42, gDisFvr; extern bool gHybInDissolve; if (gHybInDissolve) ++gDisIn42; else ++gDisFvr; }
#endif
#endif
  if (B < mBActive.size()) mBActive[B] = 0;
#ifdef PS5_LEAF_SPLICE
  scalingLfOnDissolve(B);
#endif
  // Resolve UF ownership from the representative's authoritative laminar
  // ancestry immediately before breaking B apart. Cached ownership can become
  // stale across nested rotations/expansions: retaining the wrong child leaves
  // a bare representative mapped to an unrelated active blossom.
  Size retainedRep = (!mUf.empty() && mBBase[B] < mUf.size())
                   ? ufFind(mBBase[B]) : scNull();
  Size retainedChild = retainedRep < mN
                     ? directChildContaining(B, retainedRep) : scNull();
  for (Size i = 0; i < mBCycle[B].size(); ++i) {
    Size child = mBCycle[B][i];
    mBParent[child] = scNull();
    bool wasDormant = child < mEventDormantFlag.size() &&
                      mEventDormantFlag[child] != 0;
    rememberCompactEventLabel(child);
    if (wasDormant) ++mCounters.eventDissolutionResetCandidates;
  }
  // UNION-FIND re-split (S5 dissolve :585-595): each freed child's leaves re-point to that
  // child's base rep, so ufFind resolves to the now-smaller outer node. (Guarded by mUf sized.)
  if (!mUf.empty()) {
    for (Size i = 0; i < mBCycle[B].size(); ++i) {
      Size child = mBCycle[B][i];
      Size rep = mBBase[child];
      if (rep >= mUf.size()) continue;
      // shrink() unions every cycle child into B's base representative.  On
      // expansion, the child that contains that same base can retain the
      // representative: once the other children are split away, every leaf
      // left on baseRep belongs to this child already.  Rewriting the retained
      // child's complete leaf set is redundant and is quadratic on a nested
      // chain (large base child plus two small siblings at every level).
      if (child == retainedChild && retainedRep != scNull()) {
        mUf[retainedRep] = retainedRep;
        if (!mBaseToNode.empty() && retainedRep < mBaseToNode.size())
          mBaseToNode[retainedRep] = child;
        mCounters.dissolutionRetainedBaseLeaves += nodeLeafCount(child);
        continue;
      }
#ifdef SCALING_LEAF_NOVEC
      Size end = child < mLfTailX.size() ? mLfTailX[child] : scNull();
      for (Size v = child < mLfHeadX.size() ? mLfHeadX[child] : scNull(), next = scNull();
           v != scNull(); v = next) {
        next = (v == end) ? scNull() : mLfNextX[v];
        if (v < mUf.size()) mUf[v] = rep;
      }
#else
      std::vector<Size>& cl = mNodeLeaves[child];
      for (Size j = 0; j < cl.size(); ++j) if (cl[j] < mUf.size()) mUf[cl[j]] = rep;
#endif
      mUf[rep] = rep;
      // base-vertex -> outer-node map: child is now an outer node (rep = its base).
      if (!mBaseToNode.empty() && rep < mBaseToNode.size()) mBaseToNode[rep] = child;
    }
  }
#ifdef EVENT_S5_TARGET_KEY_SHADOW
  if (mS5TargetKeyActive) s5TargetKeyAfterDissolve(B, mTNow);
#endif
}

// FULLY FLATTEN blossom B down to bare vertices: recursively dissolve B and all
// its descendant blossoms, set z=0 throughout. Used to REPAIR a malformed blossom
// that holds a FREE non-base member (its cycle invariant is broken). After this
// every original vertex under B is a bare outer node; the search rebuilds clean
// blossoms. Distributes B's z into member y's (like liquidate) to preserve
// domination -- yz on edges only rises, never falls.
inline void ScalingMatching::flattenNode(Size B) {
  CRUMB("flattenNode", B, 0);
  if (B < mN) return;
  if (mZ[B] != 0) { liquidate(B); }   // y(u)+=z/2; z=0; dissolve (1 level)
  else dissolve(B);
  // recurse into former children (now bare-rooted nodes) that are themselves blossoms
  std::vector<Size> kids = mBCycle[B];
  for (Size i = 0; i < kids.size(); ++i)
    if (kids[i] >= mN) flattenNode(kids[i]);
}

inline void ScalingMatching::augmentFrom(Size E, Size W, Size uOrig, Size vOrig) {
  ++mCounters.augmentations;
  CRUMB("augmentFrom", E, W);
#ifdef SCDBG
  Size c0=0; for (Size v=0; v<mN; ++v) if (mMate[v]!=scNull()) ++c0;
#endif
  if (E >= mN) rotateToBase(E, uOrig);
  if (W >= mN) rotateToBase(W, vOrig);
  flipUpToRoot(E);    // marks E-side path leaves used (O(path), when mIncrShrink)
  // The event loop may locally coerce a dead-tree target to FREE while its
  // stored label and parent chain remain stale. Following that chain can flip
  // into an unrelated live tree and leave asymmetric mates. A real cross-tree
  // target is EVEN; a genuinely free target has no path to flip.
  if (mNLabel[W] == LBL_EVEN) flipUpToRoot(W);
  mMate[uOrig] = vOrig; mMate[vOrig] = uOrig;
  // mark the connecting edge endpoints (the augmenting path's two ends) -- O(1).
  if (mIncrShrink && !mUsedV.empty()) {
    if (uOrig < mUsedV.size()) mUsedV[uOrig] = 1;
    if (vOrig < mUsedV.size()) mUsedV[vOrig] = 1;
  }
#ifdef SCDBG
  Size c1=0; for (Size v=0; v<mN; ++v) if (mMate[v]!=scNull()) ++c1;
  if (c1 <= c0) std::fprintf(stderr,"  !!AUGMENT DID NOT GROW card %lu->%lu (E=%lu W=%lu via %lu-%lu)\n",(unsigned long)c0,(unsigned long)c1,(unsigned long)E,(unsigned long)W,(unsigned long)uOrig,(unsigned long)vOrig);
  for (Size v = 0; v < mN; ++v) if (mMate[v] != scNull()) {
    if (mMate[v] >= mN) std::fprintf(stderr,"    !!RANGE mate[%lu]=%lu\n",(unsigned long)v,(unsigned long)mMate[v]);
    else if (mMate[mMate[v]] != v) std::fprintf(stderr,"    !!ASYMM mate[%lu]=%lu\n",(unsigned long)v,(unsigned long)mMate[v]);
  }
#endif
}

// Validate the paper's EdmondsSearch/SearchOne entry contract. Explicit root
// sets are always checked because this costs only O(|F|). A
// SCALING_INVARIANTS build additionally checks inferred roots, granularity,
// domination/near-domination, matched-edge tightness, and blossom-dual signs.
// The full pass is intentionally absent from headline builds because it scans
// all edges at every search entry.
inline bool ScalingMatching::validateSearchEntry(
    const SearchRequest& request, bool full) const {
  int rootParity = -1;
  if (request.roots != NULL) {
    for (Size i = 0; i < request.roots->size(); ++i) {
      Size v = (*request.roots)[i];
      if (v >= mN || v >= mMate.size() || mMate[v] != scNull() ||
          (mShellActive && (v >= mInShell.size() || !mInShell[v])))
        return false;
      int p = (int)(mY[v] & 1);
      if (rootParity < 0) rootParity = p;
      else if (p != rootParity) {
#ifdef SCALING_INVARIANT_TRACE
        std::fprintf(stderr,
          "[SCINV] inferred-root-parity v=%lu y=%lld parity=%d expected=%d crit=%d\n",
          (unsigned long)v, (long long)mY[v], p, rootParity,
          request.criterion);
#endif
        return false;
      }
    }
  } else if (full) {
    for (Size v = 0; v < mN; ++v) {
      if (v >= mMate.size() || mMate[v] != scNull() ||
          (mShellActive && (v >= mInShell.size() || !mInShell[v])))
        continue;
      int p = (int)(mY[v] & 1);
      if (rootParity < 0) rootParity = p;
      else if (p != rootParity) return false;
    }
  }
  if (!full) return true;

  const LL dominationFloor = request.criterion == CRIT1 ? 0 : -2;
  for (Size u = 0; u < mN; ++u) {
    for (Size i = 0; i < mAdj[u].size(); ++i) {
      Size v = mAdj[u][i];
      if (v >= mN || v <= u) continue;
      if ((mW[u][i] & 1) != 0) {
#ifdef SCALING_INVARIANT_TRACE
        std::fprintf(stderr, "[SCINV] odd-weight edge=%lu-%lu w=%lld crit=%d\n",
          (unsigned long)u, (unsigned long)v, (long long)mW[u][i],
          request.criterion);
#endif
        return false;
      }
      LL slack = slackZ(u, v);
      if (slack < dominationFloor) {
#ifdef SCALING_INVARIANT_TRACE
        std::fprintf(stderr,
          "[SCINV] domination edge=%lu-%lu slack=%lld floor=%lld crit=%d\n",
          (unsigned long)u, (unsigned long)v, (long long)slack,
          (long long)dominationFloor, request.criterion);
#endif
        return false;
      }
      if (u < mMate.size() && mMate[u] == v) {
        if (request.criterion == CRIT1) {
          if (slack != 0) return false;
        } else if (slack > 0) {
#ifdef SCALING_INVARIANT_TRACE
          std::fprintf(stderr,
            "[SCINV] matched-near-tight edge=%lu-%lu slack=%lld crit=%d\n",
            (unsigned long)u, (unsigned long)v, (long long)slack,
            request.criterion);
#endif
          return false;
        }
      }
    }
  }
  for (Size b = mN; b < mNumNodes && b < mZ.size(); ++b) {
    if (b < mBActive.size() && !mBActive[b] &&
        !(mOldForestSeparated && b >= mAuditOldBegin && b < mAuditOldEnd))
      continue;
    if (mZ[b] < 0 || (mZ[b] & 1) != 0) {
#ifdef SCALING_INVARIANT_TRACE
      std::fprintf(stderr, "[SCINV] blossom-dual B=%lu z=%lld crit=%d\n",
        (unsigned long)b, (long long)mZ[b], request.criterion);
#endif
      return false;
    }

    // Omega' records are dual containers during Section 4.2, not active
    // current-scale blossoms. Their old cycle edges obey the preceding scale's
    // near-tightness statement and are deliberately excluded here.
    if (mOldForestSeparated && b >= mAuditOldBegin && b < mAuditOldEnd)
      continue;
    if (b >= mBActive.size() || !mBActive[b]) continue;
    if (b >= mBCycle.size() || b >= mBCycEdgeU.size() ||
        b >= mBCycEdgeV.size() ||
        mBCycEdgeU[b].size() != mBCycle[b].size() ||
        mBCycEdgeV[b].size() != mBCycle[b].size()) {
#ifdef SCALING_INVARIANT_TRACE
      std::fprintf(stderr,
        "[SCINV] cycle-shape B=%lu active=%d z=%lld cycle=%lu eu=%lu ev=%lu crit=%d\n",
        (unsigned long)b, b < mBActive.size() ? (int)mBActive[b] : -1,
        (long long)mZ[b],
        (unsigned long)(b < mBCycle.size() ? mBCycle[b].size() : 0),
        (unsigned long)(b < mBCycEdgeU.size() ? mBCycEdgeU[b].size() : 0),
        (unsigned long)(b < mBCycEdgeV.size() ? mBCycEdgeV[b].size() : 0),
        request.criterion);
#endif
      return false;
    }
    for (Size k = 0; k < mBCycEdgeU[b].size(); ++k) {
      Size u = mBCycEdgeU[b][k], v = mBCycEdgeV[b][k];
      if (u >= mN || v >= mN) {
#ifdef SCALING_INVARIANT_TRACE
        std::fprintf(stderr,
          "[SCINV] cycle-endpoint B=%lu k=%lu edge=%lu-%lu n=%lu crit=%d\n",
          (unsigned long)b, (unsigned long)k, (unsigned long)u,
          (unsigned long)v, (unsigned long)mN, request.criterion);
#endif
        return false;
      }
      bool found = false;
      for (Size ei = 0; ei < mAdj[u].size(); ++ei)
        if (mAdj[u][ei] == v) { found = true; break; }
      if (!found) {
#ifdef SCALING_INVARIANT_TRACE
        std::fprintf(stderr,
          "[SCINV] cycle-missing-edge B=%lu k=%lu edge=%lu-%lu crit=%d\n",
          (unsigned long)b, (unsigned long)k, (unsigned long)u,
          (unsigned long)v, request.criterion);
#endif
        return false;
      }
      LL slack = slackZ(u, v);
      if (request.criterion == CRIT1) {
        if (slack != 0) {
#ifdef SCALING_INVARIANT_TRACE
          std::fprintf(stderr,
            "[SCINV] cycle-slack B=%lu k=%lu edge=%lu-%lu slack=%lld crit=1 z=%lld parent=%ld\n",
            (unsigned long)b, (unsigned long)k, (unsigned long)u,
            (unsigned long)v, (long long)slack, (long long)mZ[b],
            (long)(b < mBParent.size() && mBParent[b] != scNull()
                     ? (long)mBParent[b] : -1L));
#endif
          return false;
        }
      } else if (slack < -2 || slack > 0) {
#ifdef SCALING_INVARIANT_TRACE
        std::fprintf(stderr,
          "[SCINV] cycle-slack B=%lu k=%lu edge=%lu-%lu slack=%lld crit=%d z=%lld parent=%ld\n",
          (unsigned long)b, (unsigned long)k, (unsigned long)u,
          (unsigned long)v, (long long)slack, request.criterion,
          (long long)mZ[b],
          (long)(b < mBParent.size() && mBParent[b] != scNull()
                   ? (long)mBParent[b] : -1L));
#endif
        return false;
      }
    }
  }
  return true;
}

inline bool ScalingMatching::validateEventInvariantState(int criterion) {
  ++mCounters.eventInvariantChecks;
  for (Size v = 0; v < mMate.size(); ++v) {
    Size w = mMate[v];
    if (w == scNull()) continue;
    if (w >= mMate.size() || w == v || mMate[w] != v) {
      ++mCounters.eventInvariantFailures;
      mLastStatus = STATUS_INVARIANT_VIOLATION;
      return false;
    }
  }
  SearchRequest request;
  std::vector<Size> noExplicitRoots;
  request.criterion = criterion;
  request.roots = &noExplicitRoots;
  if (validateSearchEntry(request, true)) return true;
  ++mCounters.eventInvariantFailures;
  mLastStatus = STATUS_INVARIANT_VIOLATION;
  return false;
}

inline ScalingMatching::SearchResult
ScalingMatching::executeSearch(const SearchRequest& request, SearchBackend backend) {
  SearchResult result;
  if (request.criterion < CRIT1 || request.criterion > CRIT3 ||
      (request.region == REGION_SHELL) != request.shellTranslation ||
      (request.region == REGION_SHELL && !mShellActive) ||
      (request.region == REGION_VERTEX_SET &&
       (request.regionVertices == NULL || mShellActive || backend != BACKEND_SCAN)) ||
      (request.region != REGION_VERTEX_SET && request.regionVertices != NULL) ||
      (request.stoppingRule == STOP_AT_DUAL_BUDGET && request.dualBudget < 0) ||
      (request.stoppingRule == STOP_AT_SHELL_EVENT &&
       request.region != REGION_SHELL)) {
    result.status = STATUS_INVARIANT_VIOLATION;
    mLastStatus = result.status;
    return result;
  }

  bool vertexSetRegion = request.region == REGION_VERTEX_SET;
  bool savedNoTranslate = mHybNoTranslate;
  if (vertexSetRegion) {
    mInShell.assign(mN, 0);
    mInDstar.assign(mN, 0);
    mShellMarked.clear();
    for (Size i = 0; i < request.regionVertices->size(); ++i) {
      Size v = (*request.regionVertices)[i];
      if (v < mN && !mInShell[v]) {
        mInShell[v] = 1;
        mShellMarked.push_back(v);
      }
    }
    if (mShellMarked.empty()) {
      result.status = STATUS_INVARIANT_VIOLATION;
      mLastStatus = result.status;
      return result;
    }
    mShellActive = true;
    mShellCstar = scNull();
    mShellDstar = scNull();
    mShellLocalResetReady = false;
    mHybNoTranslate = true;
  }

  bool fullEntryCheck = false;
#ifdef SCALING_INVARIANTS
  fullEntryCheck = true;
#endif
  if (!validateSearchEntry(request, fullEntryCheck)) {
    if (vertexSetRegion) {
      mShellActive = false;
      mShellCstar = scNull();
      mShellDstar = scNull();
      mShellLocalResetReady = false;
      mShellMarked.clear();
      std::fill(mInShell.begin(), mInShell.end(), 0);
      mHybNoTranslate = savedNoTranslate;
    }
    result.status = STATUS_INVARIANT_VIOLATION;
    mLastStatus = result.status;
    return result;
  }

  // Both implementations operate directly on this object's canonical
  // M/y/z/Omega state. The event driver is exposed through the same contract;
  // shell requests are accepted only with Figure 7's exact criterion, stopping
  // rule, translation flag, and MCM identity.
  if (backend == BACKEND_EVENT) {
    if (request.region == REGION_SHELL) {
      if (request.criterion != CRIT1 || !request.maximalAugment ||
          request.stoppingRule != STOP_AT_SHELL_EVENT ||
          !request.shellTranslation || request.roots != NULL) {
        result.status = STATUS_INVARIANT_VIOLATION;
        mLastStatus = result.status;
        return result;
      }
      // Figure 7's Augmentation step is an exact MCM operation in the
      // currently eligible shell. Saturate at the current dual first. If no
      // path exists, advance with the translated event clock to the next
      // augmentation or inherited-boundary dissolution, then saturate again.
      result.augmented = eventShellMcmAtCurrentDual();
      if (!result.augmented && mLastStatus == STATUS_SUCCESS) {
        bool first = edmondsSearchEvent(
          CRIT1, NULL, -1, false, false, request.maxAdjustments);
        LL spent = mLastSearchDualUnits;
        if (first && mLastStatus == STATUS_SUCCESS) {
          result.augmented = true;
          eventShellMcmAtCurrentDual();
          mLastSearchDualUnits = spent;
        }
      }
    } else {
      // SearchOne's maximal augmentation is a maximal set of mutually
      // vertex-disjoint paths at one dual time, not Figure 7's MCM operation.
      result.augmented = edmondsSearchEvent(
        request.criterion, request.roots, request.dualBudget,
        request.exhaustBudget,
        request.batchSameTimeAugment || request.maximalAugment,
        request.maxAdjustments);
    }
  } else {
    result.augmented = edmondsSearchEx(
      request.criterion, request.maximalAugment, request.maxAdjustments,
      request.batchEventFreeTime, request.persistAcrossAdjust,
      request.persistAcrossAugment, request.dualBudget, request.roots,
      request.exhaustBudget);
  }
  result.dualUnits = mLastSearchDualUnits;
  result.reachedBudget = request.dualBudget >= 0 &&
                         result.dualUnits == request.dualBudget;
  result.noProgress = !result.augmented && result.dualUnits == 0;
  result.status = mLastStatus;
  if (vertexSetRegion) {
    mShellActive = false;
    mShellCstar = scNull();
    mShellDstar = scNull();
    mShellLocalResetReady = false;
    mShellMarked.clear();
    std::fill(mInShell.begin(), mInShell.end(), 0);
    mHybNoTranslate = savedNoTranslate;
  }
  return result;
}

// EdmondsSearch(F) per Fig 2 (and SearchOne via maximalAugment + maxAdjust=1).
// Relabels the forest once, roots all FREE outer nodes EVEN, then loops:
//   grow eligible forest (BFS) -> on EVEN-EVEN same-tree eligible edge SHRINK,
//   on EVEN-EVEN cross-tree (or EVEN-free) eligible edge AUGMENT. When growth
//   reaches a fixpoint with no augment: do ONE unit DUAL ADJUSTMENT (y-=1 V_out,
//   y+=1 V_in, z+=2/-=2 on root blossoms), dissolve zero-z inner blossoms, and
//   continue. Halt (return true) after the first round that augments >=1 path;
//   if maximalAugment, keep augmenting vertex-disjoint paths within that round.
//   Return false if no augment is possible (no eligible event, no dual move).
// maxAdjust caps the number of applied adjustment events (<0 = unbounded).
// maxDualUnits caps the sum of their deltas, which is the paper's FVR clock.
inline bool ScalingMatching::edmondsSearch(int crit, bool maximalAugment,
                                           LL maxAdjust, bool pqDelta,
                                           LL maxDualUnits,
                                           const std::vector<Size>* explicitRoots,
                                           bool exhaustDualBudget) {
  return edmondsSearchEx(crit, maximalAugment, maxAdjust, pqDelta,
                         /*persistAcrossAdjust=*/false,
                         /*persistAcrossAugment=*/false, maxDualUnits, explicitRoots,
                         exhaustDualBudget);
}

// edmondsSearchEx: edmondsSearch + the persistence flag. When
// persistAcrossAdjust is true, a PURE DUAL-ADJUST round does NOT trigger the
// repair+relabel+regrow at the loop top on the next iteration -- the forest
// (labels, tree pointers) persists, because a dual adjustment changes y/z values
// but NO structure. This removes the per-dual-adjust relabel (the safe half of
// the persistent-forest win; see scaling-rebuild-bottleneck-diagnosed). Augment
// and shrink/dissolve still relabel (Step 2b will persist across those too).
// Gated identical-weight vs the persistAcrossAdjust=false reference.
// 2b liveness: the alternating tree containing outer node N is LIVE iff its root
// still has an UNMATCHED base vertex. After an augment consumes a tree, its root's
// base becomes matched, so every node in that tree returns false here -- the
// persistent grow then ignores those stale labels without a relabel pass. Resolve
// the stored root through outerOf (it may have been absorbed into a blossom).
inline bool ScalingMatching::treeLive(Size N) const {
  Size r = mNRoot[N];
  if (r == scNull()) return false;
  Size ro = outerOf(r);                 // current outer node of the root
  Size base = mBBase[ro];               // root tree's exposed (free) vertex
  return mMate[base] == scNull();
}

inline void ScalingMatching::nextGrowScanEpoch() {
  if (!mNoRescanGrow) return;
  if (mGrowScanEpoch.size() != mN) {
    mGrowScanEpoch.assign(mN, 0);
    mGrowScanPos.assign(mN, 0);
    mGrowScanClock = 0;
  }
  ++mGrowScanClock;
  if (mGrowScanClock == 0) {             // unsigned wrap: reserve zero as the invalid generation
    std::fill(mGrowScanEpoch.begin(), mGrowScanEpoch.end(), 0);
    ++mGrowScanClock;
  }
}

inline bool ScalingMatching::edmondsSearchEx(int crit, bool maximalAugment,
                                             LL maxAdjust, bool pqDelta,
                                             bool persistAcrossAdjust,
                                             bool persistAcrossAugment,
                                             LL maxDualUnits,
                                             const std::vector<Size>* explicitRoots,
                                             bool exhaustDualBudget) {
  ++mCounters.searchInitializations;
  LL adjustments = 0;
  const bool shellLocalReset = mShellActive && !mHybridCanonical &&
                               std::getenv("MBX_SHELL_LOCAL_RESET");
  // A non-shell search labels nodes without the shell-local touched list.
  // Force the next shell to establish a fresh full-reset baseline.
  if (!mShellActive && mOldForestSeparated) mShellLocalResetReady = false;
  mLastSearchDualUnits = 0;
  // Build an epoch-marked root set once per search call.  The common null case has
  // no marking cost and preserves the all-free-root policy byte for byte.
  Size rootEpoch = 0;
  if (explicitRoots && !mShellActive) {
    if (mSearchRootEpoch.size() < mN) mSearchRootEpoch.resize(mN, 0);
    if (++mSearchRootClock == 0) {
      std::fill(mSearchRootEpoch.begin(), mSearchRootEpoch.end(), 0);
      ++mSearchRootClock;
    }
    rootEpoch = mSearchRootClock;
    for (Size i = 0; i < explicitRoots->size(); ++i) {
      Size v = (*explicitRoots)[i];
      if (v < mN) mSearchRootEpoch[v] = rootEpoch;
    }
  }
  bool augmentedAny = false;
  nextGrowScanEpoch();
  bool needRelabel = true;        // first iteration always relabels (fresh roots)
  // S5-style disjointness marks for no-relabel multi-augment (mIncrShrink + persistAcrossAugment):
  // cleared on each relabel (fresh phase); a vertex on an augmenting path is marked; later edges
  // touching a marked vertex are skipped so the multi-augment stays VERTEX-DISJOINT (S5 :831/927).
  bool useUsed = mIncrShrink && persistAcrossAugment;
  if (useUsed && mUsedV.size() < mN) mUsedV.assign(mN, 0);
  // BOUNDED-CALL no-progress guard: a BOUNDED call (maxAdjust>=0, i.e. step-7 SearchOne) does
  // <=maxAdjust dual adjustments, so a healthy call needs very few loop iters. If it ever cycles
  // far past that it is in a grow/shrink/relabel spin -- we abandon the call cleanly (return
  // augmentedAny; the matching/duals are valid at loop top, so the tau-loop just sees no progress
  // this call -- never a wrong answer). Cheap insurance for the bounded path. (NOTE: this does
  // NOT cover the known setPersistSearch n=1600 hang -- that spin is downstream in the UNBOUNDED
  // Finalization Crit3 search, which this guard deliberately does not touch; persist stays gated
  // OFF for that reason. The guard only protects bounded SearchOne calls.)
  LL cycleGuard = 0;
  LL preDissolveCont = 0;   // count of pre-adjust-dissolve `continue`s this call (see guard below)
  // TIGHT bound for BOUNDED calls (step-7 SearchOne, maxAdjust small): a healthy maximal-augment
  // + <=maxAdjust-adjust search visits the forest O(maxAdjust+1) times, so a few * (n + maxAdjust)
  // outer iters is a generous ceiling. The persistent-forest variant can otherwise do a LARGE
  // (finite) number of no-progress relabel/grow cycles on rare instances (the n=1600 "hang" was
  // pathological SLOWNESS, not a true infinite loop -- no single sub-loop spins; the outer loop
  // just iterates ~tens of thousands of times making microscopic progress). Bailing the call
  // cleanly (return augmentedAny -- valid matching/duals at loop top) lets the tau-loop proceed;
  // correctness is preserved (a missed augment this call is fine, the next call / Finalization
  // re-tries). Bound scales with maxAdjust so unbounded callers (maxAdjust<0) are NOT limited here.
  const LL cycleMax = 16 * (LL)mN + 64 + 8 * (maxAdjust < 0 ? 0 : maxAdjust);
#ifdef SCPROF
  LL iters = 0; ++gProfCalls;
#endif
  for (;;) {
#ifdef SCPROF
    ++iters; ++gProfIters;
#endif
    // persistent-path no-progress guard (see cycleMax above): ONLY for BOUNDED calls
    // (maxAdjust>=0, i.e. step-7 SearchOne) -- those do <=maxAdjust dual adjustments so a
    // healthy call needs very few loop iters; vastly more = a stale-label spin. We abandon
    // the call cleanly (return augmentedAny -- the matching/duals are valid at loop top, the
    // tau-loop just sees no progress this call: never a wrong answer, at worst a missed
    // speedup). NOT applied to unbounded calls (maxAdjust<0, e.g. computePerfect oracle),
    // which legitimately loop many times to drive one augmenting path -- their own
    // adjustments-based safety bound (below) governs them, so the proven paths are untouched.
    if (maxAdjust >= 0 && ++cycleGuard > cycleMax) {
      if (mHybridCanonical) mLastStatus = STATUS_UNEXPECTED_NO_PROGRESS;
      return augmentedAny;
    }
    bool anyRoot = false;
    if (needRelabel) {
#ifdef SCPROF
      ++gProfRelabels;
      { extern long long gRelabelShell, gRelabelMain; if (mShellActive) ++gRelabelShell; else ++gRelabelMain; }
#endif
      // REPAIR malformed blossoms FIRST: a FREE vertex must be a blossom's BASE
      // (or a bare node). Flatten any blossom holding a free NON-BASE member, then
      // RELABEL the alternating forest from the current matching + persistent
      // blossom structure.
      { LT_START;
      for (Size v = 0; v < mN; ++v)
        if (mMate[v] == scNull()) {
          if (rootEpoch != 0 && mSearchRootEpoch[v] != rootEpoch) continue;
          Size N = outerOf(v);
#ifdef CYCDBG
          Size fc = 0;
          while (N >= mN && mBBase[N] != v) {
            if (++fc > mNumNodes + 8) {
              std::fprintf(stderr, "[CYCDBG] FLATTEN-SPIN free v=%lu stuck under N=%lu (mBBase[N]=%lu) -- flattenNode not dissolving\n",
                           (unsigned long)v, (unsigned long)N, (unsigned long)mBBase[N]); std::fflush(stderr);
              break;
            }
            flattenNode(N); N = outerOf(v);
          }
#else
          while (N >= mN && mBBase[N] != v) {
#ifdef HYBPROF
            const Size auditN = N;
            const bool auditInherited = auditN >= mAuditOldBegin && auditN < mAuditOldEnd;
            const LL auditZBefore = auditN < mZ.size() ? mZ[auditN] : 0;
            if (gHybInDissolve) {
              ++gHybFlattenDuringDissolve;
              if (auditInherited) {
                ++gHybFlattenInherited;
                if (auditZBefore > 0) {
                  ++gHybFlattenInheritedPos;
                  gHybFlattenInheritedZBefore += auditZBefore;
                }
              } else {
                ++gHybFlattenCurrent;
                if (auditZBefore > 0) ++gHybFlattenCurrentPos;
              }
              if (std::getenv("MBX_TRACE_OLD_FLATTEN") && auditInherited &&
                  mAuditFlattenTraceCount < 24)
                std::fprintf(stderr,
                  "[OLD-FLATTEN before] scale=%lu freeV=%lu node=%lu old=[%lu,%lu) base=%lu z=%lld\n",
                  (unsigned long)mAuditScale, (unsigned long)v, (unsigned long)auditN,
                  (unsigned long)mAuditOldBegin, (unsigned long)mAuditOldEnd,
                  (unsigned long)mBBase[auditN], (long long)auditZBefore);
            }
#endif
            if (mOldForestSeparated && N >= mAuditOldEnd) {
              // A current Omega root with one exposed non-base vertex is still
              // a valid critical blossom: rotate its internal matching/base to
              // that exposed vertex. Liquidating it changes dual parity and
              // discards a positive-z active blossom.
              Size freeInside = 0;
              std::vector<Size> nodeLeaves; collectNodeLeaves(N, &nodeLeaves);
              for (Size qi = 0; qi < nodeLeaves.size(); ++qi) {
                Size qv = nodeLeaves[qi];
                if (qv < mN && mMate[qv] == scNull()) ++freeInside;
              }
              if (freeInside == 1) {
#ifdef HYBPROF
                ++gHybCurrentBaseRotations;
#endif
                rotateToBase(N, v);
                N = outerOf(v);
              } else {
#ifdef HYBPROF
                ++gHybCurrentMultiFree;
                if (std::getenv("MBX_TRACE_CURRENT_MULTI")) {
                  std::fprintf(stderr,
                    "[CURRENT-MULTI] scale=%lu shell=%d freeV=%lu node=%lu base=%lu z=%lld leaves=",
                    (unsigned long)mAuditScale, mShellActive ? 1 : 0,
                    (unsigned long)v, (unsigned long)N,
                    (unsigned long)mBBase[N], (long long)mZ[N]);
                  for (Size qi = 0; qi < nodeLeaves.size(); ++qi) {
                    Size qv = nodeLeaves[qi];
                    std::fprintf(stderr, " %lu(m=%ld,y=%lld)",
                      (unsigned long)qv,
                      (long)(qv < mN && mMate[qv] != scNull() ? mMate[qv] : -1),
                      (long long)(qv < mN ? mY[qv] : 0));
                  }
                  std::fprintf(stderr, "\n");
                }
#endif
                flattenNode(N); N = outerOf(v);
              }
            } else {
              flattenNode(N); N = outerOf(v);
            }
#ifdef HYBPROF
            if (gHybInDissolve && auditInherited && auditZBefore > 0 &&
                auditN < mZ.size() && mZ[auditN] == 0)
              ++gHybFlattenInheritedZCleared;
            if (gHybInDissolve && std::getenv("MBX_TRACE_OLD_FLATTEN") &&
                auditInherited && mAuditFlattenTraceCount < 24) {
              std::fprintf(stderr,
                "[OLD-FLATTEN after ] scale=%lu node=%lu z:%lld->%lld parent=%ld\n",
                (unsigned long)mAuditScale, (unsigned long)auditN,
                (long long)auditZBefore,
                (long long)(auditN < mZ.size() ? mZ[auditN] : 0),
                (long)(auditN < mBParent.size() && mBParent[auditN] != scNull()
                       ? mBParent[auditN] : -1));
              ++mAuditFlattenTraceCount;
            }
#endif
          }
#endif
        }
      LT_END(gLtFlatten); }
      { LT_START;
      if (mIncrRelabel || (shellLocalReset && mShellLocalResetReady)) {
        // INCREMENTAL reset: only PREVIOUSLY-LABELED nodes were non-FREE; resetting just those is
        // identical to the full O(mNumNodes) sweep (every other node is already FREE) but O(labeled).
        // This is the per-call entry-relabel cliff fix (measured: ~970 relabels x O(n)/scale).
        // Reset EVERY previously-labeled node (NOT just current outer nodes): a node can be ABSORBED
        // into a blossom while labeled (mBParent set), then later DISSOLVED back to outer carrying a
        // stale EVEN/ODD label -> the card=0/hang bug. Clearing its label unconditionally is safe
        // (an absorbed node's label is not read while inside the blossom).
        Size retainedAbsorbed = 0;
        for (Size j = 0; j < mLabeledList.size(); ++j) { Size i = mLabeledList[j];
          if (shellLocalReset && !mIncrRelabel && i < mBParent.size() &&
              mBParent[i] != scNull()) {
            // Match the reference full reset (which skips absorbed nodes), but
            // retain the id so a later dissolution exposing it can reset it.
            mLabeledList[retainedAbsorbed++] = i;
            continue;
          }
          if (i < mNLabel.size()) {
            mNLabel[i] = LBL_FREE_OF_TREE;
            mNParent[i] = scNull(); mNParEdgeU[i] = scNull(); mNParEdgeV[i] = scNull();
            mNRoot[i] = scNull();
          }
        }
        if (shellLocalReset && !mIncrRelabel)
          mLabeledList.resize(retainedAbsorbed);
        else
          mLabeledList.clear();
      } else {
        for (Size i = 0; i < mNumNodes; ++i) {
          if (mBParent[i] != scNull()) continue;
          mNLabel[i] = LBL_FREE_OF_TREE;
          mNParent[i] = scNull(); mNParEdgeU[i] = scNull(); mNParEdgeV[i] = scNull();
          mNRoot[i] = scNull();
        }
        if (shellLocalReset) {
          mLabeledList.clear();
          mShellLocalResetReady = true;
        }
      }
      LT_END(gLtReset); }
      { LT_START;
      Size rootScanN = mShellActive ? mShellMarked.size() : mN;
      for (Size rootScanI = 0; rootScanI < rootScanN; ++rootScanI) {
        Size v = mShellActive ? mShellMarked[rootScanI] : rootScanI;
        if (v >= mN || (mShellActive && !mInShell[v])) continue;
        if (mMate[v] == scNull()) {
          if (rootEpoch != 0 && mSearchRootEpoch[v] != rootEpoch) continue;
          Size N = outerOf(v);
          if (mNLabel[N] == LBL_FREE_OF_TREE) { mNLabel[N]=LBL_EVEN; mNRoot[N]=N; anyRoot=true;
            if (mIncrRelabel || shellLocalReset) mLabeledList.push_back(N); }
        }
      }
      LT_END(gLtReroot); }
#ifdef UFCHECK
      ufVerify("post-relabel");
#endif
      if (useUsed) std::fill(mUsedV.begin(), mUsedV.begin() + mN, (char)0);  // fresh phase
    } else {
      // persistent path: forest kept from the previous (dual-adjust) iteration.
      // A root exists iff some free vertex's outer node is still EVEN-rooted.
      for (Size i = 0; i < mNumNodes; ++i)
        if (mBParent[i] == scNull() && mNLabel[i] == LBL_EVEN && mNRoot[i] == i)
          { anyRoot = true; break; }
    }
    if (!anyRoot) return augmentedAny;

    // ---- GROW to fixpoint via worklist of EVEN outer nodes ----
    // 2b persistence: a tree is LIVE iff its root's base vertex is unmatched.
    // After an augment, the two consumed trees' roots' bases become matched, so
    // their nodes auto-die (treeLive==false) -- no relabel needed. We re-derive
    // liveness on the fly during grow (cheap; only over labeled nodes).
    bool didAugment = false, didShrink = false;
    // EDM: reuse a member buffer for the grow worklist instead of allocating a fresh
    // std::vector every grow pass (callgrind: vector malloc/free/realloc/memcpy was ~13-18%
    // of the run -- this wq, rebuilt millions of times, was the bulk). clear() keeps capacity.
    std::vector<Size>& wq = mGrowWq; wq.clear();
    { LT_START;
    for (Size i = 0; i < mNumNodes; ++i)
      if (mBParent[i] == scNull() && mNLabel[i] == LBL_EVEN) {
        wq.push_back(i); ++mCounters.queueInsertions;
      }
    LT_END(gLtWorklist); }
    double _ltgrow0; (void)_ltgrow0;
#ifdef LOOPTIME
    _ltgrow0 = ltNow();
#endif
    Size wh = 0;
    while (wh < wq.size() && !didShrink) {
      Size E = wq[wh++];
      ++mCounters.queuePops;
      if (mBParent[E] != scNull() || mNLabel[E] != LBL_EVEN) continue;
      if (persistAcrossAugment && !treeLive(E)) continue;
      // iterate the vertices INSIDE node E directly (O(|E|)) instead of scanning
      // all n vertices for membership (the old O(n^2) grow). mNodeLeaves[E] holds
      // exactly E's original vertices.
      // NO COPY (the copy was O(|E|) per EVEN node -> O(n^2)/grow at cold scales where blossoms
      // are huge; callgrind: 80% of the n=12800 spike was here). Index mNodeLeaves[E] FRESH each
      // iteration: a shrink/blossom-formation that resizes the OUTER mNodeLeaves vector would
      // dangle a held reference, but re-subscripting mNodeLeaves[E][ui] re-resolves through the
      // (possibly reallocated) outer vector -- safe. And a shrink sets didShrink/incrBroke, which
      // breaks the loop before the next access, so E's own leaf list never changes under us.
      bool incrBroke = false;   // incremental shrink absorbed E -> stop scanning its stale leaves
      std::vector<Size> growLeaves; collectNodeLeaves(E, &growLeaves);
      for (Size ui = 0; ui < growLeaves.size() && !didShrink && !incrBroke; ++ui) {
        Size u = growLeaves[ui];
        // ufFind(u) is INVARIANT across the inner edge loop (u's outer node only changes on a
        // shrink, which sets didShrink/incrBroke and exits). Hoist it ONCE per u -- it was being
        // recomputed per-edge (x mAdj[u].size()) in the UF same-blossom test, the single hottest
        // line in the program (ufFind path-halving = 14% of total runtime, callgrind n=6400).
        Size ruU = (!mUf.empty() && u < mUf.size()) ? ufFind(u) : u;
        // Hoist the adjacency-row pointers (mAdj is the static input graph -- never mutated during
        // search; only blossom/forest arrays change), so the inner loop indexes a raw pointer rather
        // than re-resolving nested std::vector operator[]/.size() per edge (~13% of runtime was
        // std::vector indexing internals, callgrind). mW/mWhat rows are parallel to mAdj[u].
        const std::vector<Size>& adjU = mAdj[u];
        const Size* adjUp = adjU.data(); Size degU = adjU.size();
        const LL* wRow = mUseScaleW ? mW[u].data() : mWhat[u].data();
        Size scanBegin = 0;
        if (mNoRescanGrow) {
          if (mGrowScanEpoch[u] == mGrowScanClock) scanBegin = mGrowScanPos[u];
          else { mGrowScanEpoch[u] = mGrowScanClock; mGrowScanPos[u] = 0; }
        }
        for (Size k = scanBegin; k < degU && !didShrink && !incrBroke; ++k) {
          // Mark before processing: if this edge forms a blossom, the replacement blossom resumes
          // at k+1 rather than rediscovering the same shrink edge. Edges before k keep unchanged
          // slack and endpoint roles across a shrink-only structural update.
          if (mNoRescanGrow) mGrowScanPos[u] = k + 1;
          Size w = adjUp[k];
          ++mCounters.edgeScans;
#ifdef OPCNT
          ++gOpEdgeScans;
#endif
#ifdef HYBPROF
          if (gHybInDissolve) ++gHybDissolveScans;
#endif
          // UF same-blossom fast-path (S5 mechanism): if w and u share an outer node, the edge
          // is internal -> skip in O(alpha), avoiding the O(depth) outerOf(w) walk. Equivalent to
          // the W==E check below (UFCHECK proved ufFind agrees with outerOf) but cheaper. u is
          // inside the current EVEN node E (u in Everts), so ufFind(u)==base(E).
          // Compute ufFind(w) ONCE: reused by the same-blossom test AND outerOfUFRep below (was
          // computed twice/edge -- once here, once inside outerOfUF, a measured redundancy).
          bool ufOnG = !mUf.empty();
          Size rw = ufOnG ? ufFind(w) : w;
#ifdef OPCNT
          if (ufOnG && rw == ruU) { ++gOpUfSteps; continue; }   // UF skip fired
#else
          if (ufOnG && rw == ruU) continue;
#endif
          // S5 disjointness: skip an edge whose endpoint was consumed by an augmenting path this
          // phase (keeps the no-relabel multi-augment vertex-disjoint). (S5 :927.) u is fresh
          // (we'd not be scanning E if E were used), w may be on a prior path this phase.
          if (useUsed && (mUsedV[w] || mUsedV[u])) continue;
          Size W = ufOnG ? outerOfUFRep(rw, w) : outerOf(w);   // reuse rw=ufFind(w), no recompute
#ifdef OPCNT
          if (W == E) ++gOpUfMiss;   // UF said different but outerOf says same = UF false-negative
#endif
          if (W == E) continue;                 // inside same blossom (UF miss safety net)
          // HYBRID ShellSearch: a STRADDLING edge (w outside the current atomic shell) is
          // ineligible -- the search is confined to G(C*,D*). (No-op when !mShellActive.)
          if (mShellActive && !mInShell[w]) continue;
          // UF on => the same-blossom `continue` above proved this edge is cross-blossom, so use the
          // ufFind-free eligibility (slack = mY[u]+mY[w]-w2). Use slackZW when UF is off.
          if (!(mUf.empty() ? eligibleAt(u, k, crit) : eligibleAtCross(u, k, crit))) continue;
#ifdef OPCNT
          ++gOpEligible;   // PROBE: this edge passed eligibility (became a real grow/shrink/augment candidate)
#endif
          // 2b: a labeled node whose tree was consumed by a prior augment this
          // round is effectively FREE (re-growable). PHYSICALLY RETIRE it -- not just a
          // local relabel: clear its stale tree pointers (mNLabel/mNParent/mNRoot) so that
          // if we augment ONTO it, flipUpToRoot(W) sees a true isolated free root and walks
          // NO stale parent chain. (THE n=1600 HANG: previously only the LOCAL lblW was set
          // to FREE while mNParent[W]/mNRoot[W] still pointed into the dead tree, so an
          // augment onto W flipped the wrong edges, didn't increase cardinality, and left W
          // re-augmentable -> in-grow multi-augment loop. Retiring W here makes the augment
          // valid and the consumed node truly free.)
          int lblW = mNLabel[W];
          if (persistAcrossAugment && lblW != LBL_FREE_OF_TREE && !treeLive(W)) {
            lblW = LBL_FREE_OF_TREE;
            mNLabel[W] = LBL_FREE_OF_TREE;
            mNParent[W] = scNull(); mNParEdgeU[W] = scNull(); mNParEdgeV[W] = scNull();
            mNRoot[W] = scNull();
          }
          if (lblW == LBL_EVEN) {
            // Same tree? Use lcaOuter as the SOURCE OF TRUTH (mNRoot can go stale
            // across blossom form/dissolve). lca != null => same tree => shrink.
            if (lcaOuter(E, W) != scNull()) {
              // INCREMENTAL SHRINK (mIncrShrink): shrink() leaves the new blossom B EVEN with
              // valid tree pointers, so the forest is consistent -- NO full relabel needed.
              // Capture B (= mNumNodes before the ++ inside shrink), push it onto the grow
              // worklist, and keep growing. break the edge scan (E is now absorbed into B,
              // its Everts list is stale); the wq while-loop continues with B (and other
              // still-EVEN nodes). Stale absorbed entries in wq self-skip (mBParent!=null guard).
              if (mIncrShrink) {
                Size B = mNumNodes;            // shrink does B = mNumNodes++ exactly once
                shrink(E, W, u, w);
                wq.push_back(B); ++mCounters.queueInsertions;
                incrBroke = true; break;       // leave both loops; do NOT set didShrink (no relabel)
              }
              shrink(E, W, u, w); didShrink = true; break;
            }
            else {
              // cross-tree augmenting path. If maximalAugment, record and keep
              // going to find a vertex-disjoint set; else augment immediately.
              augmentFrom(E, W, u, w);
              didAugment = true; augmentedAny = true;
              nextGrowScanEpoch();       // matching/tree roles changed; prior row decisions are stale
              if (persistAcrossAugment) {
                // PERSIST: do NOT relabel. The two consumed trees auto-die
                // (their roots' bases are now matched). Keep growing the survivors
                // from the current worklist (other EVEN nodes still queued).
                // KNOWN ISSUE: this in-grow multi-augment persist is the SOURCE of BOTH the
                // big speedup (it is the whole 3.7x-22x Edm win) AND the rare hang -- the 2b
                // dead-tree liveness (treeLive) does not reliably retire consumed EVEN nodes
                // on some instances (n=1600), so they get re-augmented in a cycle. Fixing
                // treeLive to correctly retire consumed trees would make persist default-safe;
                // until then setPersistSearch is gated OFF. (Forcing a relabel here removes the
                // hang but ALSO the entire speedup -- measured -- so they are inseparable as-is.)
                break;     // break neighbor loop; outer while continues with wq
              }
              if (!maximalAugment) { return true; }
              didShrink = true; break;            // maximal: re-loop (relabel+regrow)
            }
          } else if (lblW == LBL_ODD) {
            continue;
          } else {
            // W FREE_OF_TREE (or 2b-dead). matched -> GROW; unmatched -> augment.
            Size bw = mBBase[W];
            Size m = mMate[bw];
            if (m == scNull()) {
              augmentFrom(E, W, u, w);
              didAugment = true; augmentedAny = true;
              nextGrowScanEpoch();       // matching/tree roles changed; prior row decisions are stale
              if (persistAcrossAugment) { break; }  // persist: keep growing survivors
              if (!maximalAugment) return true;
              didShrink = true; break;
            }
            Size M = outerOf(m);
#ifdef MATCHED_GROW_AUDIT
            if (crit == CRIT2 || crit == CRIT3) {
              ++mCounters.matchedGrowTransitions;
              LL ms = slackZ(bw, m);
              bool me = crit == CRIT2 ? ms == 0 : (ms == 0 || ms == -2);
              if (!me) {
                ++mCounters.matchedGrowIneligible;
                LL target = (crit == CRIT3 && ms <= -2) ? -2 : 0;
                LL delay = target - ms;
                if (delay == 1) ++mCounters.matchedGrowDelayOne;
                else if (delay == 2) ++mCounters.matchedGrowDelayTwo;
                else ++mCounters.matchedGrowDelayOther;
              }
            }
#endif
            mNLabel[W] = LBL_ODD; mNParent[W] = E;
            mNParEdgeU[W] = w; mNParEdgeV[W] = u;
            mNRoot[W] = mNRoot[E];
            if (mHybridCanonical && (crit == CRIT2 || crit == CRIT3)) {
              bool matchedEligible = eligible(bw, m, crit);
              if (mNLabel[M] == LBL_FREE_OF_TREE && matchedEligible) {
                mNLabel[M] = LBL_EVEN; mNParent[M] = W;
                mNParEdgeU[M] = m; mNParEdgeV[M] = bw;
                mNRoot[M] = mNRoot[E];
                wq.push_back(M); ++mCounters.queueInsertions;
              } else if (mNLabel[M] == LBL_ODD && matchedEligible) {
                // Both ends were reached independently before their matched
                // edge became eligible. Rebuild the auxiliary forest at this
                // same dual time; the eligible edge is then handled in the
                // ordinary outer/inner ordering without a persistent
                // inner--inner event representation.
                didShrink = true;
              }
              if (mIncrRelabel || shellLocalReset) mLabeledList.push_back(W);
              if (didShrink) break;
              continue;
            }
            mNLabel[M] = LBL_EVEN; mNParent[M] = W;
            mNParEdgeU[M] = m; mNParEdgeV[M] = bw;
            mNRoot[M] = mNRoot[E];
            if (mIncrRelabel || shellLocalReset) {
              mLabeledList.push_back(W); mLabeledList.push_back(M);
            }
            wq.push_back(M); ++mCounters.queueInsertions;
          }
        }
      }
    }
#ifdef LOOPTIME
    gLtGrow += ltNow() - _ltgrow0;
#endif
    // A shrink or (maximal) augment changed the structure: loop -> relabel+regrow.
    if (didShrink) { needRelabel = true;
#ifdef SCPROF
      ++gRelShrink;
#endif
      continue; }
    // growth fixpoint, no augment this round.
    // D step3 path-(b) ATTEMPT (multi-level batching in the LOOP engine) FAILED: crossing a dual
    // level here without retiring consumed-tree labels + without lazy duals = the cross-level
    // staleness bug (broke 0/100; documented in bounded-scaling-loop-status: continuation across a
    // dual level on NON-LAZY duals produces invalid matchings). mMultiLevel kept as a gated no-op
    // pending path-(a) (use the S5 engine, which HAS lazy duals, as the FVR searcher). See roadmap.
    // CLIFF FIX (mContinueAug, opt-in): the paper halts after an augment (return true) -> the tau-loop
    // re-calls edmondsSearchEx, and EACH re-call pays an O(n) ENTRY RELABEL (measured: ~970 calls x O(n)/
    // scale = THE cliff). With persistAcrossAugment the forest SURVIVES the augment, so instead of
    // returning we can RETIRE only the just-consumed (dead) trees -- O(consumed), not O(n) -- and
    // CONTINUE the same call across the next dual level. This collapses ~970 calls/scale to ~1, removing
    // the per-call entry relabel. RISK: cross-level staleness (the old mMultiLevel broke 0/100 on
    // non-lazy duals) -- the retire-dead-trees below is exactly the missing piece that broke it.
    if (augmentedAny) {
      if (!(mContinueAug && persistAcrossAugment)) return true;   // default: paper-faithful halt
      // INCREMENTAL RETIRE: a tree is dead iff its root's base is now matched (treeLive==false). Clear
      // ONLY dead nodes' labels (re-rootable); live trees keep their valid labels. O(labeled nodes),
      // and in practice O(consumed) since only the augmented trees died. Then fall through to the dual
      // adjustment + continue (no return -> no entry relabel next iteration).
      for (Size i = 0; i < mNumNodes; ++i) {
        if (mBParent[i] != scNull()) continue;
        if (mNLabel[i] == LBL_FREE_OF_TREE) continue;
        if (!treeLive(i)) {
          mNLabel[i] = LBL_FREE_OF_TREE; mNParent[i] = scNull();
          mNParEdgeU[i] = scNull(); mNParEdgeV[i] = scNull(); mNRoot[i] = scNull();
        }
      }
      // a still-free vertex whose outer node went FREE must re-root EVEN for the next level's grow.
      for (Size v = 0; v < mN; ++v)
        if (mMate[v] == scNull()) {
          if (rootEpoch != 0 && mSearchRootEpoch[v] != rootEpoch) continue;
          if (mShellActive && !mInShell[v]) continue;
          Size N = outerOf(v);
          if (mNLabel[N] == LBL_FREE_OF_TREE) { mNLabel[N] = LBL_EVEN; mNRoot[N] = N; }
        }
      if (useUsed) std::fill(mUsedV.begin(), mUsedV.begin() + mN, (char)0);  // fresh phase for next level
      // fall through to DUAL ADJUSTMENT, continue the same call (needRelabel stays false)
    }

    // ---- DUAL ADJUSTMENT ----
    if (maxAdjust >= 0 && adjustments >= maxAdjust) return augmentedAny;
    // EVEN outer nodes' vertices move y down; ODD up; EVEN blossom z up, ODD down.
    // PQSEARCH (pqDelta=true): jump by the LARGEST safe delta = min over the next
    // eligibility/dissolution events, instead of one unit. This is the paper's
    // PQSearch (EdmondsSearch with unbounded, priority-queued dual adjustments,
    // O(m+n log n)); it collapses up to maxW unit steps into O(#events) jumps.
    // Unit-step (pqDelta=false) is the slow correctness-reference path.
    // EVENT FEASIBILITY (computed in BOTH modes -- this is the fix for the dense
    // non-convergence: unit-step previously adjusted blindly with delta=1 and
    // NEVER checked whether any event could fire, so when the last free vertices
    // were structurally unable to expose an augmenting path the duals drifted
    // forever. `best` = min steps until SOME eligibility/dissolution event fires;
    // best<0 => no event possible => no augmenting path => return.
    // pqDelta=true jumps by `best` (PQSearch); pqDelta=false steps by 1 (correct
    // reference) but uses `best` only to DETECT the no-event termination.
    double _ltdual0; (void)_ltdual0;
#ifdef LOOPTIME
    _ltdual0 = ltNow();
#endif
    LL best = -1;
    if (mNoRescanDual) { mDualWakeU.clear(); mDualWakeK.clear(); }
    // EDM (the big one): iterate ONLY the EVEN outer nodes and their member vertices, instead of
    // `for u < mN` + outerOf(u) per vertex. The old full sweep called outerOf O(n) times per
    // outer iteration JUST to find the EVEN vertices -- and `while(mBParent[c]!=null)` (outerOf)
    // was 18% of the whole run at n=3200 by sheer call volume (callgrind). EVEN nodes are exactly
    // the roots in mNRoot[E]==E... no: any EVEN outer node. We scan EVEN root-blossom nodes via
    // mNodeLeaves (their original vertices), matching the grow loop's iteration. This makes the
    // event scan O(sum of EVEN-tree degrees) instead of O(n + m) -- the per-iteration O(n) tax
    // that made the wall-clock grow super-linearly while op-counts stayed flat.
    for (Size E = 0; E < mNumNodes; ++E) {
      if (mBParent[E] != scNull() || mNLabel[E] != LBL_EVEN) continue;  // EVEN outer nodes only
      std::vector<Size> Everts; collectNodeLeaves(E, &Everts);
      for (Size ui = 0; ui < Everts.size(); ++ui) {
        Size u = Everts[ui];
        Size ruU = (!mUf.empty() && u < mUf.size()) ? ufFind(u) : u;  // invariant across edge loop (no shrink here)
        for (Size k = 0; k < mAdj[u].size(); ++k) {
          ++mCounters.edgeScans;
          Size w = mAdj[u][k];
          // UF O(alpha) same-blossom skip (mirror the grow loop): internal edge -> not an event.
          bool ufOn = !mUf.empty();
          Size rw = ufOn ? ufFind(w) : w;   // compute ufFind(w) ONCE (reused by outerOfUFRep below)
          if (ufOn && rw == ruU) continue;
          Size W = ufOn ? outerOfUFRep(rw, w) : outerOf(w);   // reuse rw, no recompute
          if (W == E) continue;
          // When UF is on, the continue above PROVED ufFind(w)!=ufFind(u) = different outer nodes =>
          // slackZW's z-sum is 0 (no common enclosing blossom). So slack == mY[u]+mY[w]-w2 EXACTLY,
          // and calling slackZW would only re-run 2 ufFinds + a chain walk that returns early. Inline
          // it (skips 2 ufFind calls/edge -- ufFind was 14% of runtime). Use the exact path when UF is off.
          LL w2e = mUseScaleW ? mW[u][k] : 2 * mWhat[u][k];
          LL s = ufOn ? (mY[u] + mY[w] + oldCommonZ(u, w) - w2e)
                      : slackZW(u, w, w2e);
          int lw = mNLabel[W];
          LL target = 0;
          bool hasFutureTarget = true;
          if (crit == CRIT2) {
            if (mMate[u] == w) {
              target = 0;
              // slack*=-slack for a matched edge.  In the two endpoint roles
              // considered below slack only decreases, so negative slack moves
              // away from eligibility rather than toward it.
              if (s < 0) hasFutureTarget = false;
            } else target = -2;
          } else if (crit == CRIT3 && s < 0) {
            // Paper Section 5: slack*=slack for slack>=0, and slack+2 for
            // slack in {-1,-2}.  Thus positive slack first reaches 0; only a
            // negative near-dominating edge targets -2.
            target = -2;
          }
          LL gap = s - target;
          if (gap < 0) hasFutureTarget = false;
          if (lw == LBL_EVEN && hasFutureTarget) {
            // both EVEN: slack drops by 2/step; an odd gap cannot hit the
            // exact eligibility value under the fixed forest roles.
            if ((gap & 1) != 0) continue;
            LL d = gap / 2;
            if (d >= 0) {
              if (best < 0 || d < best) {
                best = d;
                if (mNoRescanDual) { mDualWakeU.clear(); mDualWakeK.clear();
                  mDualWakeU.push_back(u); mDualWakeK.push_back(k); }
              } else if (mNoRescanDual && d == best) {
                mDualWakeU.push_back(u); mDualWakeK.push_back(k);
              }
            }
          } else if (lw == LBL_FREE_OF_TREE && hasFutureTarget) {
            LL d = gap;                            // EVEN->unlabeled: drops by 1/step
            if (d >= 0) {
              if (best < 0 || d < best) {
                best = d;
                if (mNoRescanDual) { mDualWakeU.clear(); mDualWakeK.clear();
                  mDualWakeU.push_back(u); mDualWakeK.push_back(k); }
              } else if (mNoRescanDual && d == best) {
                mDualWakeU.push_back(u); mDualWakeK.push_back(k);
              }
            }
          }
          // EVEN->ODD: yz unchanged by adjustment (no event from this edge).
        }
      }
    }
    if (mHybridCanonical && (crit == CRIT2 || crit == CRIT3)) {
      // In the paper-reference forest a reached INNER endpoint does not pull
      // its mate into the tree until the matched edge is eligible. If the mate
      // is unreached the slack rises by one per unit; if it was independently
      // reached INNER, it rises by two and becomes a blossom event.
      for (Size W = 0; W < mNumNodes; ++W) {
        if (mBParent[W] != scNull() || mNLabel[W] != LBL_ODD) continue;
        Size bw = mBBase[W];
        if (bw >= mN || mMate[bw] == scNull()) continue;
        Size m = mMate[bw], M = outerOf(m);
        int mateLabel = mNLabel[M];
        if (M == W || (mateLabel != LBL_FREE_OF_TREE && mateLabel != LBL_ODD))
          continue;
        if (mateLabel == LBL_ODD && bw > m) continue; // one event per pair
        LL s = slackZ(bw, m);
        LL target = (crit == CRIT3 && s <= -2) ? -2 : 0;
        LL gap = target - s;
        LL rate = mateLabel == LBL_ODD ? 2 : 1;
        if (gap < 0 || (gap % rate) != 0) continue;
        LL d = gap / rate;
        if (best < 0 || d < best) {
          best = d;
          if (mNoRescanDual) { mDualWakeU.clear(); mDualWakeK.clear(); }
        }
      }
    }
    // d4: ODD root blossom z hits 0 (dissolution) after z/2 steps.
    for (Size b = mN; b < mNumNodes; ++b)
      if (mBParent[b] == scNull() && mNLabel[b] == LBL_ODD) {
        LL d = mZ[b] / 2;
        if (d >= 0 && (best < 0 || d < best)) {
          best = d;
          if (mNoRescanDual) { mDualWakeU.clear(); mDualWakeK.clear(); }
        }
      }
    if (best < 0) {
      // DismantlePath stage 2 must perform exactly T units even if no structural
      // event remains. Batch that final inert interval to its requested boundary.
      // Ordinary searches retain the no-event early exit.
      if (exhaustDualBudget && maxDualUnits >= 0 &&
          mLastSearchDualUnits < maxDualUnits)
        best = maxDualUnits - mLastSearchDualUnits;
      else
        return augmentedAny;
    }
    // Dissolve any ODD blossom that ALREADY has z==0 BEFORE adjusting -- otherwise
    // the adjustment would push its z negative and (clamped at 0) raise yz on its
    // internal edges by 2*delta, breaking tightness of internal MATCHED edges
    // (THE n=8/n=60 BUG: matched edges drift to slackZ=2 -> suboptimal). After
    // dissolving, relabel/re-grow so the freed sub-nodes get proper labels.
    bool dissolvedAny = false;
    bool preDissolvedLiveOdd = false;   // a dissolved z=0 ODD blossom that was on a LIVE tree
    for (Size b = mN; b < mNumNodes; ++b)
      if (mBParent[b] == scNull() && mNLabel[b] == LBL_ODD && mZ[b] == 0)
        {
          // (SPIN-FIX attempt: skipping never-z-positive dissolutions broke correctness 831/81 --
          // some freshly-born z=0 ODD dissolutions ARE genuinely needed expansions, not spin. The
          // everPos heuristic over-corrects. Reverted. The spin is subtler. See hybrid-roadmap.)
          { bool plive = mIncrShrink && treeLive(b);
            bool prelinked = (plive && mIncrDissolve) ? relinkOddCycle(b) : false;
#ifdef SCPROF
            if (plive && mIncrDissolve) { extern long long gRelinkOK, gRelinkBail; if (prelinked) ++gRelinkOK; else ++gRelinkBail; }
#endif
            if (plive && !prelinked) preDissolvedLiveOdd = true;
            // ALWAYS retire b's own label (FREE): b is expanded into its sub-nodes (relink wired the
            // sub-nodes' tree labels; b itself is no longer an outer node). Leaving mNLabel[b]=ODD
            // would re-detect b next iteration (the re-dissolve spin). relink sets sub-node labels,
            // not b's, so FREE'ing b here is safe in both paths.
            dissolve(b); mNLabel[b] = LBL_FREE_OF_TREE; dissolvedAny = true; } }
        // RETIRE the dissolved blossom's OWN label (set FREE): dissolve(b) only clears its
        // CHILDREN's mBParent; b itself keeps mBParent[b]==scNull() && mNLabel[b]==ODD && mZ[b]==0,
        // so without this the re-detect condition stays TRUE and b is re-dissolved every iteration.
        // This is a correct invariant fix (a dissolved blossom is no longer an active ODD node) and
        // is harmless to the non-persist path (which relabels fully anyway). NOTE: it removes ONE
        // layer of the persist re-dissolve spin but NOT the full n=1600 persist hang -- the deeper
        // cause is that the dissolve->relabel->regrow cycle re-creates a dissolvable ODD-z=0 blossom
        // (persist forest-state consistency); persist stays gated OFF. See edm-cost memory.
    // PRE-DISSOLVE CONTINUE GUARD (cheap correct safety net): each genuine dissolve permanently
    // removes a blossom, so legitimate pre-adjust-dissolve `continue`s in ONE call are bounded by
    // the blossom count. Capping at 2*mNumNodes+8 and returning cleanly on exceed prevents a
    // per-call dissolve->relabel->regrow non-termination from spinning forever (returns
    // augmentedAny -- valid matching/duals at this point, never a wrong answer). NOTE: this bounds
    // each CALL but does NOT fully fix the persist n=1600 hang (which recurs ACROSS calls -- see
    // edm-cost memory, BUG C); persist stays gated OFF. Harmless to all callers (a healthy search
    // dissolves far fewer than 2*mNumNodes times per call).
    if (dissolvedAny) {
      nextGrowScanEpoch();               // exposed/relabelled nodes require fresh adjacency decisions
      // NO-RELABEL pre-dissolution (mIncrShrink): a DEAD-tree z=0 ODD blossom dissolving doesn't
      // need a relabel (it is off any live path) -- keep growing in place (UF re-split in dissolve
      // resolves the freed nodes). This BREAKS the dissolve->relabel->regrow->re-form SPIN (the
      // cold-scale 589-relabel spike): no relabel => no regrow => no re-created z=0 blossom. Only a
      // LIVE-tree z=0 ODD dissolution still relabels (re-thread the path). Non-incr path unchanged.
      if (mIncrShrink && !preDissolvedLiveOdd) { continue; }   // dead-tree: keep growing, no relabel
      if (++preDissolveCont <= 2 * (LL)mNumNodes + 8) { needRelabel = true; continue; }
      if (mHybridCanonical) mLastStatus = STATUS_UNEXPECTED_NO_PROGRESS;
      return augmentedAny;
    }

    // Unit mode steps by 1 (safe: all zero-z ODD blossoms were dissolved above,
    // so the min ODD-blossom z is >=2 and z-=2 stays >=0 -- no clamp, no slack
    // drift on internal matched edges). PQSearch jumps by `best`.
    LL delta = pqDelta ? (best == 0 ? 1 : best) : 1;
    // Figure 3/6 budgets count UNITS of dual adjustment, not event jumps.  A delta jump is
    // an equivalent batching of event-free unit steps, but it consumes `delta` units.  Clamp
    // to the caller's remaining budget; when the next event lies beyond that boundary we
    // materialize only the remaining inert units and return without processing the future event.
    if (maxDualUnits >= 0) {
      LL remaining = maxDualUnits - mLastSearchDualUnits;
      if (remaining <= 0) return augmentedAny;
      if (delta > remaining) delta = remaining;
    }
    // SHELL: the unit translation decrements z(C*)/z(D*) by 2*delta. Clamp delta so neither
    // goes BELOW 0 (the paper dissolves a shell blossom exactly when its z hits 0; overshoot
    // to negative z = the domination-break bug). After this adjust we halt if z reached 0.
    bool shellExhausted = false;
    if (mShellActive && !mHybNoTranslate) {
      if (mShellCstar != scNull()) { LL room = mZ[mShellCstar] / 2; if (room < delta) delta = room; }
      if (mShellDstar != scNull()) { LL room = mZ[mShellDstar] / 2; if (room < delta) delta = room; }
      // delta==0 => z(C*) or z(D*) already 0 => the shell blossom is exhausted; do NOT force
      // a +1 step (that overshoots z to -2 = the domination bug). Halt the search instead.
      if (delta <= 0) shellExhausted = true;
    }
    if (shellExhausted) return augmentedAny;   // §4.2 halting condition 3: C*/D* dissolved

    ScalingNanosecondGuard shellDualTimer(
      mShellActive ? &mCounters.shellDualAdjustmentNanoseconds : NULL);
    Size dualVCount = mShellActive ? mShellMarked.size() : mN;
    if (mShellActive) mCounters.shellDualVertexVisits += dualVCount;
    for (Size dualVI = 0; dualVI < dualVCount; ++dualVI) {
      Size v = mShellActive ? mShellMarked[dualVI] : dualVI;
      if (v >= mN) continue;
      Size V = outerOfUF(v);   // O(alpha) (was O(depth) outerOf -- the dual-adjust hot loop)
      if (mNLabel[V] == LBL_EVEN) mY[v] -= delta;
      else if (mNLabel[V] == LBL_ODD) mY[v] += delta;
      // Fold Figure 7's boundary translation into the same shell pass. It is
      // algebraically independent of the alternating-tree role adjustment.
      if (mShellActive && !mHybNoTranslate) {
        if (mInDstar[v]) mY[v] += 2 * delta;
        else if (mInShell[v]) mY[v] += delta;
      }
    }
    for (Size b = mN; b < mNumNodes; ++b)
      if (mBParent[b] == scNull()) {
        // SHELL: C*/D* z is owned by the UNIT TRANSLATION below; if they're also tree-labeled
        // (C* is typically ODD), the standard z-adjust would DOUBLE-decrement -> z overshoots
        // negative (z(C*) 60->-2 bug). Skip C*/D* here; translation handles them.
        if (mShellActive && (b == mShellCstar || b == mShellDstar)) continue;
        if (mNLabel[b] == LBL_EVEN) { mZ[b] += 2 * delta; if (b < mBEverPos.size()) mBEverPos[b] = 1; }
        else if (mNLabel[b] == LBL_ODD) mZ[b] -= 2 * delta;  // bounded by best, stays >=0
      }
    // HYBRID §4.2 ShellSearch UNIT TRANSLATION of C*/D* (per unit of dual adjustment, so
    // scaled by delta): z(C*)-=2; z(D*)-=2 (if D*); y+=2 on D*, y+=1 on C*\D*. Keeps yz of
    // edges straddling C*\D* non-decreasing => preserves domination (Fig 8). No-op when
    // !mShellActive. (mShellCstar/mShellDstar are the shell's enclosing undissolved blossoms.)
    if (mShellActive && !mHybNoTranslate) {
      if (mShellCstar != scNull()) mZ[mShellCstar] -= 2 * delta;
      if (mShellDstar != scNull()) mZ[mShellDstar] -= 2 * delta;
#ifdef HYBPROF
      if (gHybInDissolve) {
        int parity = -1; Size first = scNull();
        for (Size v = 0; v < mN; ++v) if (mMate[v] == scNull()) {
          int p = (int)(mY[v] & 1);
          if (parity < 0) { parity = p; first = v; }
          else if (p != parity) {
            ++gHybFreeParityViol;
            if (std::getenv("MBX_TRACE_HYB_PARITY") && gHybFreeParityViol <= 24)
              std::fprintf(stderr,
                "[HYB-PARITY] scale=%lu after-adjust first=v%lu y=%lld p=%d other=v%lu y=%lld p=%d C=%lu D=%ld delta=%lld outer=%lu uf=%lu label=%d\n",
                (unsigned long)mAuditScale, (unsigned long)first,
                (long long)mY[first], parity, (unsigned long)v, (long long)mY[v], p,
                (unsigned long)mShellCstar,
                (long)(mShellDstar == scNull() ? -1 : mShellDstar),
                (long long)delta, (unsigned long)outerOf(v),
                (unsigned long)outerOfUF(v), mNLabel[outerOfUF(v)]);
            break;
          }
        }
      }
#endif
      // §4.2 HALTING CONDITION 3: C* (the outermost undissolved shell blossom) dissolves when
      // z(C*) hits 0 -> the shell search MUST stop (running on would over-reduce duals and
      // break domination -- the bug). Also stop if D* dissolves (z(D*)=0). Signal halt by
      // returning; caller (shellSearch) treats it as completion.
      if ((mShellCstar != scNull() && mZ[mShellCstar] <= 0) ||
          (mShellDstar != scNull() && mZ[mShellDstar] <= 0)) {
        ++adjustments;
        mLastSearchDualUnits += delta;
        mCounters.dualUnits += (unsigned long long)delta;
        return augmentedAny;     // halt: C*/D* exhausted (Blossom Dissolution handled by caller)
      }
    }
    ++adjustments;
    mLastSearchDualUnits += delta;
    mCounters.dualUnits += (unsigned long long)delta;
#ifdef HYBPROF
    if (gHybInDissolve) ++gHybDissolveAdjusts;
#endif
#ifdef SCPROF
    ++gProfAdj;
    if (gProfInFinal) ++gProfFinalAdj; else ++gProfScaleAdj;
#endif
    // ---- BLOSSOM DISSOLUTION: inner (ODD) root blossoms with z==0 ----
    bool dissolvedPost = false;
    bool dissolvedLiveOdd = false;   // a dissolved ODD blossom that was on a LIVE tree path
    for (Size b = mN; b < mNumNodes; ++b)
      if (mBParent[b] == scNull() && mNLabel[b] == LBL_ODD && mZ[b] == 0)
        { // (SPIN-FIX everPos guard reverted -- broke correctness; some never-positive dissolutions
          //  are genuine expansions. See pre-adjust site note + hybrid-roadmap.)
          // a LIVE-tree ODD blossom's dissolution breaks the alternating path -> needs relabel to
          // re-thread (incremental keep-FREE loses the path -> suboptimal, the n=120 case). A
          // dead-tree (consumed) ODD blossom is safe to dissolve no-relabel.
          bool live = mIncrShrink && treeLive(b);
          // EXPERIMENT: re-thread the tree through b's cycle so a LIVE-ODD dissolution skips the
          // full relabel (the n>=51200 cliff: 3153 relabels x O(m)). If it re-threads cleanly,
          // this b does NOT count as dissolvedLiveOdd (no relabel needed). Gated default OFF.
          bool relinked = (live && mIncrDissolve) ? relinkOddCycle(b) : false;
#ifdef SCPROF
          if (live && mIncrDissolve) { extern long long gRelinkOK, gRelinkBail; if (relinked) ++gRelinkOK; else ++gRelinkBail; }
#endif
          if (live && !relinked) dissolvedLiveOdd = true;
          dissolve(b); mNLabel[b] = LBL_FREE_OF_TREE; dissolvedPost = true; }   // ALWAYS retire b's label (relink sets sub-node labels, not b's) -- prevents the re-dissolve spin
    // PERSISTENCE: if this was a PURE dual adjustment (no dissolve), the forest
    // structure is unchanged -- only y/z shifted -- so skip the repair+relabel on
    // the next iteration and keep growing with the persisted labels. A dissolve
    // changed structure, so relabel.
    // NO-RELABEL DISSOLUTION (mIncrShrink): a dissolve frees its ODD sub-nodes (they become FREE
    // outer nodes via dissolve()'s UF re-split + mBParent->scNull). The REST of the forest is
    // unchanged, and freed nodes resolve correctly via treeParent/UF. So instead of a full O(n)
    // relabel (the n=25600 cold-scale cliff: 600+ relabels x O(n)), KEEP the forest: the freed
    // sub-nodes will be re-grown as children when an EVEN neighbor scans them (a missed grow this
    // round is re-tried next round -- safe, never a wrong matching, just possibly one extra round).
    // Re-seed the grow worklist from EVEN roots so growth continues without the relabel.
    // Skip the relabel ONLY when no LIVE-tree ODD blossom dissolved (those need re-threading).
    needRelabel = mIncrShrink ? dissolvedLiveOdd : dissolvedPost;
    if (mHybridCanonical && (crit == CRIT2 || crit == CRIT3))
      needRelabel = true;  // discard independent-inner scheduling state
    if (!mNoRescanDual || dissolvedPost || mShellActive || best <= 0) {
      // Dissolution/translation changes endpoint roles beyond the minimum-slack edge set.
      nextGrowScanEpoch();
    } else {
      // Before this adjustment every grow row was exhausted.  With unchanged forest roles, only
      // edges attaining `best` can have newly reached the eligibility target.  Rewind precisely
      // those rows to their earliest attaining edge; all other exhausted cursors remain valid.
      for (Size i = 0; i < mDualWakeU.size(); ++i) {
        Size u = mDualWakeU[i], k = mDualWakeK[i];
        if (mGrowScanEpoch[u] != mGrowScanClock) {
          mGrowScanEpoch[u] = mGrowScanClock; mGrowScanPos[u] = k;
        } else if (k < mGrowScanPos[u]) mGrowScanPos[u] = k;
      }
    }
#ifdef SCPROF
    if (needRelabel) ++gRelDissolve;
#endif
    // SearchOne ends at the unit-budget boundary.  In particular, if this adjustment made
    // an edge eligible exactly at the boundary, that event belongs to the next SearchOne call.
    if ((maxDualUnits >= 0 && mLastSearchDualUnits >= maxDualUnits) ||
        (maxAdjust >= 0 && adjustments >= maxAdjust))
      return augmentedAny;
#ifdef LOOPTIME
    gLtDualApply += ltNow() - _ltdual0;
#endif
    // safety bound. UNIT-step adjustment count scales with the WEIGHTS (each
    // augment can need O(maxW) unit steps to tighten an edge), so the bound must
    // include a weight term, not just n^2 -- otherwise large-weight instances
    // bail prematurely and return a partial matching (the n=60 weights-up-to-1000
    // card=0). (PQSearch's delta-batching is what removes the O(maxW) factor.)
    if (maxAdjust < 0 && adjustments > 8 * (LL)mN * (LL)mN + 4 * (LL)mN * (mMaxW + 1) + 1000) {
      if (mHybridCanonical) mLastStatus = STATUS_UNEXPECTED_NO_PROGRESS;
      return augmentedAny;
    }
    // loop back to grow with the new duals (labels persist).
  }
}

// Phase-1b milestone driver: run the paper's EdmondsSearch (Criterion 1,
// Property 1) repeatedly to a PERFECT matching. NOTE per the paper this single
// full-weight search is NOT guaranteed globally max-weight (that needs the
// scaling loop) -- it produces a valid perfect matching that PRESERVES Property
// 1. We validate by invariant checks; the optimality gate comes with Phase 3.
inline bool ScalingMatching::runPerfectBlossom() {
  mNoPerfect = false;
  if (mN == 0) return true;
  if (mN % 2 != 0) { mNoPerfect = true; return false; }

  Size cap = 2 * mN + 1;
  mBParent.assign(cap, scNull());
  mBBase.assign(cap, scNull());
  mBCycle.assign(cap, std::vector<Size>());
  mBCycEdgeU.assign(cap, std::vector<Size>());
  mBCycEdgeV.assign(cap, std::vector<Size>());
  mNodeLeaves.assign(cap, std::vector<Size>());
  mZ.assign(cap, 0);
  mY.assign(mN, mMaxW);                      // Property-1 domination init
  for (Size v = 0; v < mN; ++v) { mBBase[v] = v; mNodeLeaves[v].assign(1, v); }
#ifdef PS5_LEAF_SPLICE
  scalingLfInit(mN, cap);
#endif
  mMate.assign(mN, scNull());
  mNumNodes = mN;
  mNLabel.assign(cap, LBL_FREE_OF_TREE);
  mNParent.assign(cap, scNull());
  mNParEdgeU.assign(cap, scNull());
  mNParEdgeV.assign(cap, scNull());
  mNRoot.assign(cap, scNull());

  for (Size guard = 0; guard <= 2 * mN + 2; ++guard) {
    bool perfect = true;
    for (Size v = 0; v < mN; ++v) if (mMate[v] == scNull()) { perfect = false; break; }
    if (perfect) return true;
    if (!edmondsSearch(CRIT1, /*maximalAugment=*/false, /*maxAdjust=*/-1)) {
      mNoPerfect = true; return false;
    }
  }
  mNoPerfect = true;
  return false;
}

// Identical to runPerfectBlossom but uses the delta-jump (PQSearch) dual
// adjustment. The relabel/regrow at the top of edmondsSearch then runs once per
// EVENT (O(n)/augment) instead of once per UNIT step (O(maxW)/augment) -- which
// is the n=12800 -> 20s blowup. Gated vs runPerfectBlossom/brute-force to confirm
// the delta-jump produces the SAME max-weight perfect matching.
inline bool ScalingMatching::runPerfectBlossomPQ() {
  mNoPerfect = false;
  if (mN == 0) return true;
  if (mN % 2 != 0) { mNoPerfect = true; return false; }

  Size cap = 2 * mN + 1;
  mBParent.assign(cap, scNull());
  mBBase.assign(cap, scNull());
  mBCycle.assign(cap, std::vector<Size>());
  mBCycEdgeU.assign(cap, std::vector<Size>());
  mBCycEdgeV.assign(cap, std::vector<Size>());
  mNodeLeaves.assign(cap, std::vector<Size>());
  mZ.assign(cap, 0);
  mY.assign(mN, mMaxW);
  for (Size v = 0; v < mN; ++v) { mBBase[v] = v; mNodeLeaves[v].assign(1, v); }
#ifdef PS5_LEAF_SPLICE
  scalingLfInit(mN, cap);
#endif
  mMate.assign(mN, scNull());
  mNumNodes = mN;
  mNLabel.assign(cap, LBL_FREE_OF_TREE);
  mNParent.assign(cap, scNull());
  mNParEdgeU.assign(cap, scNull());
  mNParEdgeV.assign(cap, scNull());
  mNRoot.assign(cap, scNull());

  for (Size guard = 0; guard <= 2 * mN + 2; ++guard) {
    bool perfect = true;
    for (Size v = 0; v < mN; ++v) if (mMate[v] == scNull()) { perfect = false; break; }
    if (perfect) return true;
    if (!edmondsSearch(CRIT1, /*maximalAugment=*/false, /*maxAdjust=*/-1, /*pqDelta=*/true)) {
      mNoPerfect = true; return false;
    }
  }
  mNoPerfect = true;
  return false;
}

// Step 2a persistent engine: delta-jump + persist the forest across pure
// dual-adjust rounds (no repair+relabel after an adjustment that changed no
// structure). Same MWPM result as runPerfectBlossomPQ; gated identical.
inline bool ScalingMatching::runPerfectBlossomPersist() {
  mNoPerfect = false;
  if (mN == 0) return true;
  if (mN % 2 != 0) { mNoPerfect = true; return false; }

  Size cap = 2 * mN + 1;
  mBParent.assign(cap, scNull());
  mBBase.assign(cap, scNull());
  mBCycle.assign(cap, std::vector<Size>());
  mBCycEdgeU.assign(cap, std::vector<Size>());
  mBCycEdgeV.assign(cap, std::vector<Size>());
  mNodeLeaves.assign(cap, std::vector<Size>());
  mZ.assign(cap, 0);
  mY.assign(mN, mMaxW);
  for (Size v = 0; v < mN; ++v) { mBBase[v] = v; mNodeLeaves[v].assign(1, v); }
#ifdef PS5_LEAF_SPLICE
  scalingLfInit(mN, cap);
#endif
  mMate.assign(mN, scNull());
  mNumNodes = mN;
  mNLabel.assign(cap, LBL_FREE_OF_TREE);
  mNParent.assign(cap, scNull());
  mNParEdgeU.assign(cap, scNull());
  mNParEdgeV.assign(cap, scNull());
  mNRoot.assign(cap, scNull());

  for (Size guard = 0; guard <= 2 * mN + 2; ++guard) {
    bool perfect = true;
    for (Size v = 0; v < mN; ++v) if (mMate[v] == scNull()) { perfect = false; break; }
    if (perfect) return true;
    if (!edmondsSearchEx(CRIT1, /*maximalAugment=*/false, /*maxAdjust=*/-1,
                         /*pqDelta=*/true, /*persistAcrossAdjust=*/true,
                         /*persistAcrossAugment=*/true)) {
      mNoPerfect = true; return false;
    }
  }
  mNoPerfect = true;
  return false;
}

// STAGE 3 event-PQ engine driver. Same init as the others; repeatedly run the
// event-driven search to a perfect matching.
inline bool ScalingMatching::runPerfectBlossomEvent() {
  mNoPerfect = false;
  if (mN == 0) return true;
  if (mN % 2 != 0) { mNoPerfect = true; return false; }

  Size cap = 2 * mN + 1;
  mBParent.assign(cap, scNull());
  mBBase.assign(cap, scNull());
  mBCycle.assign(cap, std::vector<Size>());
  mBCycEdgeU.assign(cap, std::vector<Size>());
  mBCycEdgeV.assign(cap, std::vector<Size>());
  mNodeLeaves.assign(cap, std::vector<Size>());
  mZ.assign(cap, 0);
  mY.assign(mN, mMaxW);
  for (Size v = 0; v < mN; ++v) { mBBase[v] = v; mNodeLeaves[v].assign(1, v); }
#ifdef PS5_LEAF_SPLICE
  scalingLfInit(mN, cap);
#endif
  mMate.assign(mN, scNull());
  mNumNodes = mN;
  mNLabel.assign(cap, LBL_FREE_OF_TREE);
  mNParent.assign(cap, scNull());
  mNParEdgeU.assign(cap, scNull());
  mNParEdgeV.assign(cap, scNull());
  mNRoot.assign(cap, scNull());

  for (Size guard = 0; guard <= 2 * mN + 2; ++guard) {
    bool perfect = true;
    for (Size v = 0; v < mN; ++v) if (mMate[v] == scNull()) { perfect = false; break; }
    if (perfect) return true;
    if (!edmondsSearchEvent(CRIT1, 0, -1, false, mEventBatchSameTime, -1)) {
      mNoPerfect = true; return false;
    }
  }
  mNoPerfect = true;
  return false;
}

// ===========================================================================
// STAGE 3 — event-driven Edmonds search (paper §5 BucketSearch). Replaces the
// O(m)-per-round edge re-scan (grow-eligibility + delta-computation; the 49.6M
// slackZ calls that gprof showed dominate at n=6400) with a BucketPQ of edge
// events. Optional Section-5.2 timestamps make a dual jump O(1); grow/shrink
// transfer affected roles incrementally; no structural event may rebuild the
// forest or event queue as a recovery path.
// The win is that we POP the next tight-edge event instead of scanning all m edges.
// Built on the proven 2b persistent forest
// (treeLive). Events are scheduled at absolute time tNow + steps-to-tight; stale
// events (endpoint labels changed since scheduling) are filtered on pop by
// re-checking slack. Gated identical-weight vs the 2b/reference engine.
// ===========================================================================

// Schedule events for every edge out of EVEN node E's vertices, at the absolute
// time the edge becomes tight under continued dual adjustment.
inline void ScalingMatching::initializeEventStructures() {
  ++mCounters.sfmInitializations;
  std::vector<Size> order;
  std::vector<char> listEnd;
  std::vector<char> seenVertex(mN, 0);
  std::vector<char> seenBlossom(mNumNodes, 0);
  order.reserve(mN);
  listEnd.reserve(mN);
  mSfmSplitPivots.assign(mZ.size(), std::vector<Size>());

  // Section 5.1 initializes one permutation in which every initial blossom,
  // maximal or not, is contiguous. It then SPLITS that order as inner roots
  // dissolve; forming an outer blossom belongs to union-find and never merges
  // split-findmin lists. Build the laminar DFS order explicitly and record the
  // right boundary of every direct child while all nested intervals are still
  // visible. These leaf pivots remain stable for the complete search.
  struct InitialSfmFrame {
    Size node, next;
    InitialSfmFrame(Size n, Size i) : node(n), next(i) {}
  };
  for (Size root = mN; root < mNumNodes; ++root) {
    if (root >= mBCycle.size() || mBCycle[root].empty() ||
        (root < mBActive.size() && !mBActive[root]) ||
        mBParent[root] != scNull() ||
        (mOldForestSeparated && root >= mAuditOldBegin &&
         root < mAuditOldEnd)) continue;
    std::vector<InitialSfmFrame> stack;
    stack.push_back(InitialSfmFrame(root, 0));
    if (seenBlossom[root]) {
      mLastStatus = STATUS_INVARIANT_VIOLATION;
      return;
    }
    seenBlossom[root] = 1;
    while (!stack.empty()) {
      InitialSfmFrame& frame = stack.back();
      if (frame.node < mN) {
        Size v = frame.node;
        if (seenVertex[v]) {
          mLastStatus = STATUS_INVARIANT_VIOLATION;
          return;
        }
        seenVertex[v] = 1;
        order.push_back(v);
        listEnd.push_back(0);
        stack.pop_back();
      } else if (frame.node >= mBCycle.size() ||
                 mBCycle[frame.node].empty()) {
        mLastStatus = STATUS_INVARIANT_VIOLATION;
        return;
      } else if (frame.next < mBCycle[frame.node].size()) {
        Size child = mBCycle[frame.node][frame.next++];
        if (child >= mNumNodes ||
            (child >= mN && seenBlossom[child])) {
          mLastStatus = STATUS_INVARIANT_VIOLATION;
          return;
        }
        if (child >= mN) seenBlossom[child] = 1;
        stack.push_back(InitialSfmFrame(child, 0));
        continue;
      } else {
        stack.pop_back();
      }
      if (!stack.empty()) {
        InitialSfmFrame& parent = stack.back();
        if (parent.node >= mN &&
            parent.next < mBCycle[parent.node].size()) {
          if (order.empty()) {
            mLastStatus = STATUS_INVARIANT_VIOLATION;
            return;
          }
          if (parent.node >= mSfmSplitPivots.size())
            mSfmSplitPivots.resize(parent.node + 1);
          mSfmSplitPivots[parent.node].push_back(order.back());
        }
      }
    }
    if (order.empty()) {
      mLastStatus = STATUS_INVARIANT_VIOLATION;
      return;
    }
    listEnd.back() = 1;
  }
  // Real vertices outside current Omega are trivial root blossoms. Inherited
  // Omega' is intentionally absent from this Section-5 partition.
  for (Size v = 0; v < mN; ++v) if (!seenVertex[v]) {
    order.push_back(v);
    listEnd.push_back(1);
    seenVertex[v] = 1;
  }
  if (order.size() != mN) {
    mLastStatus = STATUS_INVARIANT_VIOLATION;
    return;
  }
  mInnerSfm.init(order, listEnd);
#ifdef EVENT_S5_SINGLETON_TARGET_FAST
  mS5TargetKeySingleton.assign(mN, 0);
  mS5TargetKeySingletonValue.assign(mN, SplitFindmin::INF);
  mS5TargetKeySingletonListId.assign(mN, scNull());
  Size singletonRunBegin = 0;
  for (Size i = 0; i < order.size(); ++i) {
    if (!listEnd[i]) continue;
    if (i == singletonRunBegin) {
      Size v = order[i];
      mS5TargetKeySingleton[v] = 1;
      mS5TargetKeySingletonListId[v] = mInnerSfm.listId(v);
    }
    singletonRunBegin = i + 1;
  }
#endif
  for (Size i = 0; i < listEnd.size(); ++i)
    if (listEnd[i]) ++mCounters.sfmInitialLists;
  for (Size b = mN; b < mSfmSplitPivots.size(); ++b)
    mCounters.sfmInitialPivots += mSfmSplitPivots[b].size();
#ifdef EVENT_SFM_SHADOW
  mSfmKeySource.assign(mN, scNull());
  mSfmKeyOwner.assign(mN, scNull());
#endif
  mEventSfmInitialized = true;
  mEventSfmScale = mAuditScale;
}

inline void ScalingMatching::eventSfmOnShrink(Size B) {
  (void)B;
  ++mCounters.sfmFormationTopologyNoops;
  // Deliberately topology-neutral. Section 5.1 assigns growing OUTER
  // blossoms to union-find. Split-findmin represents the dissolving NON-OUTER
  // partition and supports init/list/split/decreasekey/findmin--not merge.
}

inline void ScalingMatching::eventSfmOnDissolve(Size B) {
  // The Section-5 split-findmin represents current Omega only. Dissolving an
  // inherited Omega' boundary must not split a current list that crosses that
  // old boundary; Figure 7 handles inherited boundary changes separately.
  if (mOldForestSeparated && B >= mAuditOldBegin && B < mAuditOldEnd) return;
  if (B >= mBCycle.size() || mBCycle[B].size() < 2) return;
  // Initial DFS order makes every member a contiguous run. Peel each direct
  // child at its precomputed stable leaf boundary; nested members stay intact
  // and become independent lists after expansion.
  if (B >= mSfmSplitPivots.size() ||
      mSfmSplitPivots[B].size() + 1 != mBCycle[B].size()) return;
#if defined(EVENT_S5_TARGET_KEY_SHADOW) && \
    defined(EVENT_S5_SHARED_SFM_TOPOLOGY)
  if (mS5TargetKeyActive && B < mBBase.size() && mBBase[B] < mN)
    s5TargetFrontierInvalidateList(
      s5TargetFrontierListId(mBBase[B]));
#endif
  const std::vector<Size>& pivots = mSfmSplitPivots[B];
  for (Size i = 0; i < pivots.size(); ++i) {
    mInnerSfm.split(pivots[i]);
    ++mCounters.splitFindminWork;
    ++mCounters.sfmActualSplits;
  }
#ifdef EVENT_S5_SINGLETON_TARGET_FAST
  // Each direct child is now an independent list.  A one-leaf child can use
  // the flat secondary-key channel for the rest of this topology epoch.
  for (Size i = 0; i < mBCycle[B].size(); ++i) {
    Size child = mBCycle[B][i];
    if (nodeLeafCount(child) == 1 && child < mBBase.size())
      s5TargetSfmPromoteSingleton(mBBase[child]);
  }
#endif
#ifdef EVENT_SFM_TOPOLOGY_AUDIT
  for (Size i = 0; i < mBCycle[B].size(); ++i)
    if (!sfmShadowListMatches(mBCycle[B][i])) {
      ++mCounters.sfmSplitTopologyFailures;
    }
#endif
}

#ifdef EVENT_S5_TARGET_KEY_SHADOW
inline Val ScalingMatching::s5TargetSfmKey(Size v) const {
#ifdef EVENT_S5_SINGLETON_TARGET_FAST
  if (s5TargetSfmIsSingleton(v)) return mS5TargetKeySingletonValue[v];
#endif
#ifdef EVENT_S5_SHARED_SFM_TOPOLOGY
  return mInnerSfm.key2(v);
#else
  return mS5TargetKeySfm.key(v);
#endif
}

inline void ScalingMatching::s5TargetSfmSetKey(Size v, Val key) {
  ++mCounters.sfmGrowSetCalls;
  if (s5TargetSfmKey(v) == key) {
    ++mCounters.sfmGrowSetNoops;
    return;
  }
#ifdef EVENT_S5_SINGLETON_TARGET_FAST
  if (s5TargetSfmIsSingleton(v)) {
    mS5TargetKeySingletonValue[v] = key;
    return;
  }
#endif
  Size visits = 0;
#ifdef EVENT_S5_SHARED_SFM_TOPOLOGY
  visits = mInnerSfm.setkey2(v, key);
#else
  visits = mS5TargetKeySfm.setkey(v, key);
#endif
  mCounters.sfmGrowSetAncestorVisits += visits;
}

inline void ScalingMatching::s5TargetSfmDecreaseKey(Size v, Val key) {
#ifdef EVENT_S5_SINGLETON_TARGET_FAST
  if (s5TargetSfmIsSingleton(v)) {
    if (key < mS5TargetKeySingletonValue[v]) {
      mS5TargetKeySingletonValue[v] = key;
      ++mCounters.sfmGrowDecreaseCalls;
      ++mCounters.sfmGrowDecreaseImmediateStops;
    }
    return;
  }
#endif
  bool changed = false;
  Size visits = 0;
#ifdef EVENT_S5_SHARED_SFM_TOPOLOGY
  visits = mInnerSfm.decreasekey2(v, key, &changed);
#else
  visits = mS5TargetKeySfm.decreasekey(v, key, &changed);
#endif
  if (changed) {
    ++mCounters.sfmGrowDecreaseCalls;
    mCounters.sfmGrowDecreaseAncestorVisits += visits;
    if (visits == 0) ++mCounters.sfmGrowDecreaseImmediateStops;
  }
}

inline void ScalingMatching::s5TargetSfmFindMin(
    Size v, Val* key, Size* at) const {
#ifdef EVENT_S5_SINGLETON_TARGET_FAST
  if (s5TargetSfmIsSingleton(v)) {
    *key = mS5TargetKeySingletonValue[v];
    *at = v;
    return;
  }
#endif
#ifdef EVENT_S5_SHARED_SFM_TOPOLOGY
  mInnerSfm.findmin2(v, key, at);
#else
  mS5TargetKeySfm.findmin(v, key, at);
#endif
}

inline void ScalingMatching::s5TargetSfmFindMinWithList(
    Size v, Val* key, Size* at, Size* listId) const {
#ifdef EVENT_S5_SINGLETON_TARGET_FAST
  if (s5TargetSfmIsSingleton(v)) {
    *key = mS5TargetKeySingletonValue[v];
    *at = v;
    *listId = mS5TargetKeySingletonListId[v];
    return;
  }
#endif
#ifdef EVENT_S5_SHARED_SFM_TOPOLOGY
  mInnerSfm.findmin2WithListId(v, key, at, listId);
#else
  mS5TargetKeySfm.findminWithListId(v, key, at, listId);
#endif
}

inline void ScalingMatching::s5TargetSfmSublistVerts(
    Size v, std::vector<Size>& out) const {
  // Both channels use the same list membership in shared mode, so the primary
  // channel's topology traversal is authoritative for either payload.
#ifdef EVENT_S5_SHARED_SFM_TOPOLOGY
  mInnerSfm.sublistVerts(v, out);
#else
  mS5TargetKeySfm.sublistVerts(v, out);
#endif
}

#ifdef EVENT_S5_SINGLETON_TARGET_FAST
inline bool ScalingMatching::s5TargetSfmIsSingleton(Size v) const {
  return v < mS5TargetKeySingleton.size() && mS5TargetKeySingleton[v] != 0;
}

inline void ScalingMatching::s5TargetSfmPromoteSingleton(Size v) {
  if (v >= mN || v >= mS5TargetKeySingleton.size() ||
      mS5TargetKeySingleton[v]) return;
  // Preserve the leaf payload installed while it belonged to a larger list.
  // After the topology split no later operation can merge this leaf back.
  mS5TargetKeySingletonValue[v] = mInnerSfm.key2(v);
  mS5TargetKeySingletonListId[v] = mInnerSfm.listId(v);
  mS5TargetKeySingleton[v] = 1;
}

#endif

inline void ScalingMatching::s5TargetKeyReindexPayload(
    Size v, Size owner, Size targetOwner, Size sourceToken,
    bool oldLive, bool newLive) {
  if (v >= mN) return;
  const Size oldOwner = mS5TargetKeyOwner[v];
  const Size oldTargetOwner = mS5TargetKeyTargetOwner[v];
  const Size oldSourceToken = mS5TargetKeySourceToken[v];
  mS5TargetKeyOwner[v] = owner;
  mS5TargetKeyTargetOwner[v] = targetOwner;
  mS5TargetKeySourceToken[v] = sourceToken;
  if (!newLive) return;
  if (targetOwner < mNumNodes &&
      (!oldLive || oldTargetOwner != targetOwner)) {
    if (targetOwner >= mS5TargetKeyTargetOwnerTargets.size())
      mS5TargetKeyTargetOwnerTargets.resize(targetOwner + 1);
    s5TargetKeyRememberOwner(targetOwner);
    mS5TargetKeyTargetOwnerTargets[targetOwner].push_back(v);
  }
  if (owner < mNumNodes && (!oldLive || oldOwner != owner)) {
    if (owner >= mS5TargetKeyOwnerTargets.size())
      mS5TargetKeyOwnerTargets.resize(owner + 1);
    s5TargetKeyRememberOwner(owner);
    mS5TargetKeyOwnerTargets[owner].push_back(v);
  }
  if (sourceToken < mN && (!oldLive || oldSourceToken != sourceToken)) {
    s5TargetKeyRememberToken(sourceToken);
    mS5TargetKeyTreeOwnedTargets[sourceToken].push_back(v);
  }
}

inline void ScalingMatching::s5TargetKeyRememberTarget(Size target) {
  if (target >= mN) return;
  if (mS5TargetKeyTargetTouchStamp[target] == mS5TargetKeyTouchClock) return;
  mS5TargetKeyTargetTouchStamp[target] = mS5TargetKeyTouchClock;
  mS5TargetKeyTouchedTargets.push_back(target);
}

inline void ScalingMatching::s5TargetKeyRememberOwner(Size owner) {
  if (owner >= mS5TargetKeyOwnerTouchStamp.size())
    mS5TargetKeyOwnerTouchStamp.resize(owner + 1, 0);
  if (mS5TargetKeyOwnerTouchStamp[owner] == mS5TargetKeyTouchClock) return;
  mS5TargetKeyOwnerTouchStamp[owner] = mS5TargetKeyTouchClock;
  mS5TargetKeyTouchedOwners.push_back(owner);
}

inline void ScalingMatching::s5TargetKeyRememberToken(Size token) {
  if (token >= mN) return;
  if (mS5TargetKeyTokenTouchStamp[token] == mS5TargetKeyTouchClock) return;
  mS5TargetKeyTokenTouchStamp[token] = mS5TargetKeyTouchClock;
  mS5TargetKeyTouchedTokens.push_back(token);
}

inline Size ScalingMatching::s5TargetFrontierListId(Size vertex) const {
  if (vertex >= mN) return scNull();
#ifdef EVENT_S5_SINGLETON_TARGET_FAST
  if (s5TargetSfmIsSingleton(vertex))
    return mS5TargetKeySingletonListId[vertex];
#endif
#ifdef EVENT_S5_SHARED_SFM_TOPOLOGY
  return mInnerSfm.listId(vertex);
#else
  return mS5TargetKeySfm.listId(vertex);
#endif
}

inline void ScalingMatching::s5TargetFrontierReset() {
  for (Size i = 0; i < mS5TargetFrontierTouched.size(); ++i) {
    Size id = mS5TargetFrontierTouched[i];
    if (id < mS5TargetFrontier.size()) {
      mS5TargetFrontier[id].version = 0;
      mS5TargetFrontier[id].touched = false;
    }
  }
  mS5TargetFrontierTouched.clear();
  mS5TargetFrontierClock = 1;
}

inline void ScalingMatching::s5TargetFrontierInvalidateList(Size listId) {
  if (listId >= mS5TargetFrontier.size() ||
      mS5TargetFrontier[listId].version == 0) return;
  mS5TargetFrontier[listId].version = 0;
  ++mCounters.s5TargetFrontierInvalidated;
}

inline bool ScalingMatching::s5TargetKeyGap(
    LL slack, bool matched, int crit, LL* gap) const {
  if (gap == 0) return false;
  LL candidate = 0;
  if (crit == CRIT1 || (crit == CRIT2 && matched) ||
      (crit == CRIT3 && slack >= 0)) {
    candidate = slack;
  } else if (crit == CRIT2 || crit == CRIT3) {
    if (!checkedAddLL(slack, 2, &candidate)) return false;
  } else {
    return false;
  }
  if (candidate < 0) return false;
  *gap = candidate;
  return true;
}

inline void ScalingMatching::s5TargetKeyOnShrink(Size B) {
  if (B >= mBCycle.size()) return;
#ifndef EVENT_S5_SHARED_SFM_TOPOLOGY
  std::vector<Size> reps;
  reps.reserve(mBCycle[B].size());
  for (Size i = 0; i < mBCycle[B].size(); ++i) {
    Size child = mBCycle[B][i];
    if (child < mBBase.size() && mBBase[child] < mN)
      reps.push_back(mBBase[child]);
  }
  if (reps.empty()) return;
  if (mS5TargetKeyActive)
    for (Size i = 0; i < reps.size(); ++i)
      s5TargetFrontierInvalidateList(
        s5TargetFrontierListId(reps[i]));
  if (B >= mS5TargetKeySplitPivots.size())
    mS5TargetKeySplitPivots.resize(B + 1);
  std::vector<Size>& pivots = mS5TargetKeySplitPivots[B];
  pivots.clear();
  pivots.reserve(reps.size() - 1);
  for (Size i = 0; i + 1 < reps.size(); ++i)
    pivots.push_back(mS5TargetKeySfm.lastVtx(reps[i]));
  mS5TargetKeySfm.merge(reps);
#endif

  // Persistent topology is maintained even between Edmonds forests. Payload
  // ownership, however, belongs to one active forest and must not be carried
  // through an inactive matching transition.
  if (!mS5TargetKeyActive) return;

  // Every direct cycle member is now represented by outer owner B. Transfer
  // only live payloads; duplicate/stale reverse entries are discarded here.
  if (B >= mS5TargetKeyOwnerTargets.size())
    mS5TargetKeyOwnerTargets.resize(B + 1);
  s5TargetKeyRememberOwner(B);
  for (Size i = 0; i < mBCycle[B].size(); ++i) {
    Size child = mBCycle[B][i];
    if (child >= mS5TargetKeyOwnerTargets.size()) continue;
    std::vector<Size>& owned = mS5TargetKeyOwnerTargets[child];
    for (Size j = 0; j < owned.size(); ++j) {
      Size target = owned[j];
      if (target >= mN || mS5TargetKeyOwner[target] != child
          ) continue;
      Size source = mS5TargetKeySource[target];
      if (source >= mN || outerOfUF(source) != B) continue;
      s5TargetKeyReindexPayload(target, B,
        mS5TargetKeyTargetOwner[target],
        mS5TargetKeySourceToken[target], true, true);
      ++mCounters.s5TargetKeyOwnerTransfers;
    }
    std::vector<Size>().swap(owned);

    // A newly formed OUTER blossom consumes every former FREE target interval
    // of its cycle members. Invalidate those leaf payloads through the reverse
    // target-owner index; if the tree is later retired, revealFree() starts a
    // new interval with the then-current outer owner.
    if (child < mS5TargetKeyTargetOwnerTargets.size()) {
      std::vector<Size>& targeted = mS5TargetKeyTargetOwnerTargets[child];
      for (Size j = 0; j < targeted.size(); ++j) {
        Size target = targeted[j];
        if (target >= mN || mS5TargetKeyTargetOwner[target] != child
            ) continue;
        const bool oldLive = s5TargetSfmKey(target) != SplitFindmin::INF;
        s5TargetSfmSetKey(target, SplitFindmin::INF);
        mS5TargetKeySource[target] = scNull();
        s5TargetKeyReindexPayload(target, scNull(), B, scNull(),
                                  oldLive, false);
        mS5TargetKeyEdgeIndex[target] = scNull();
        mS5TargetKeyMatched[target] = 0;
      }
      std::vector<Size>().swap(targeted);
    }
  }
  s5TargetKeyTrackTreeNode(B);
}

inline void ScalingMatching::s5TargetKeyOnDissolve(Size B) {
#ifndef EVENT_S5_SHARED_SFM_TOPOLOGY
  const bool inheritedBoundary =
    mOldForestSeparated && B >= mAuditOldBegin && B < mAuditOldEnd;
#endif
  if (B >= mBCycle.size() || mBCycle[B].size() < 2) return;
#ifndef EVENT_S5_SHARED_SFM_TOPOLOGY
  if (!inheritedBoundary) {
    if (B >= mS5TargetKeySplitPivots.size() ||
        mS5TargetKeySplitPivots[B].size() + 1 != mBCycle[B].size()) return;
    if (mS5TargetKeyActive && B < mBBase.size() && mBBase[B] < mN)
      s5TargetFrontierInvalidateList(
        s5TargetFrontierListId(mBBase[B]));
    const std::vector<Size>& pivots = mS5TargetKeySplitPivots[B];
    for (Size i = 0; i < pivots.size(); ++i)
      mS5TargetKeySfm.split(pivots[i]);
  }
#endif

  if (!mS5TargetKeyActive) return;

  // relinkOddCycle assigns the exposed child roles before dissolve() removes
  // B. A FREE child starts a new target interval; no key from its former ODD
  // interval or from B's contracted ownership may survive.
  for (Size i = 0; i < mBCycle[B].size(); ++i) {
    Size child = mBCycle[B][i];
    if (child >= mNLabel.size()) continue;
    if (mNLabel[child] != LBL_FREE_OF_TREE) {
      s5TargetKeyTrackTreeNode(child);
      continue;
    }
#ifdef SCALING_LEAF_NOVEC
    std::vector<Size> leaves; scalingLfCollect(child, &leaves);
#else
    const std::vector<Size>& leaves = mNodeLeaves[child];
#endif
    for (Size j = 0; j < leaves.size(); ++j) {
      Size v = leaves[j];
      if (v >= mN) continue;
      const bool oldLive = s5TargetSfmKey(v) != SplitFindmin::INF;
      s5TargetSfmSetKey(v, SplitFindmin::INF);
      mS5TargetKeySource[v] = scNull();
      s5TargetKeyReindexPayload(v, scNull(), child, scNull(),
                                oldLive, false);
      mS5TargetKeyEdgeIndex[v] = scNull();
      mS5TargetKeyMatched[v] = 0;
      ++mCounters.s5TargetKeyDissolveResets;
    }
  }
}

inline void ScalingMatching::s5TargetKeyTrackTreeNode(Size node) {
  if (!mS5TargetKeyActive || node >= mNumNodes || node >= mNRoot.size()) return;
  Size root = mNRoot[node];
  if (root == scNull() || root >= mNumNodes) return;
  // A freshly installed root names itself and is already an outer node.
  // Avoid a union-find walk for the millions of singleton roots initialized
  // across Figure-6 SearchOne calls; grown nodes retain the general lookup.
  Size rootOuter = root == node && mBParent[node] == scNull()
    ? node : outerOf(root);
  if (rootOuter >= mBBase.size() || mBBase[rootOuter] >= mN) return;
  // The exposed root base is stable even when a root blossom forms and the
  // root node id changes. It is therefore the tree membership token.
  Size token = mBBase[rootOuter];
  if (node >= mS5TargetKeyTrackedRoot.size())
    mS5TargetKeyTrackedRoot.resize(node + 1, scNull());
  if (node >= mS5TargetKeyTrackedRootStamp.size())
    mS5TargetKeyTrackedRootStamp.resize(node + 1, 0);
  if (mS5TargetKeyTrackedRootStamp[node] == mS5TargetKeyTouchClock &&
      mS5TargetKeyTrackedRoot[node] == token) return;
  mS5TargetKeyTrackedRoot[node] = token;
  mS5TargetKeyTrackedRootStamp[node] = mS5TargetKeyTouchClock;
  if (token >= mS5TargetKeyTreeNodes.size())
    mS5TargetKeyTreeNodes.resize(token + 1);
  s5TargetKeyRememberToken(token);
  mS5TargetKeyTreeNodes[token].push_back(node);
  // Retirement needs every real vertex of both EVEN and ODD tree nodes. Keep
  // one token-deduplicated flat list so nested blossom leaves are not traversed
  // repeatedly when an augmentation consumes the tree.
#ifdef SCALING_LEAF_NOVEC
  std::vector<Size> leaves;
  scalingLfCollect(node, &leaves);
#else
  const std::vector<Size>& leaves = mNodeLeaves[node];
#endif
  for (Size i = 0; i < leaves.size(); ++i) {
    Size v = leaves[i];
    if (v >= mN) continue;
    if (mS5TargetKeyAllVertexStamp[v] == mS5TargetKeyTouchClock &&
        mS5TargetKeyAllVertexToken[v] == token) continue;
    mS5TargetKeyAllVertexToken[v] = token;
    mS5TargetKeyAllVertexStamp[v] = mS5TargetKeyTouchClock;
    mS5TargetKeyTreeAllVertices[token].push_back(v);
  }
}

inline void ScalingMatching::s5TargetKeyAfterDissolve(Size B, LL tNow) {
  if (!mS5TargetKeyActive || B >= mBCycle.size()) return;
  if (mS5TargetKeyOwnerTargets.size() < mNumNodes)
    mS5TargetKeyOwnerTargets.resize(mNumNodes);
  // Once B's parent links and UF ownership have been removed, payloads that
  // named B can be transferred to the exact exposed source owner. If that
  // source is no longer live OUTER, reveal the target's next candidate.
  if (B < mS5TargetKeyOwnerTargets.size()) {
    std::vector<Size>& owned = mS5TargetKeyOwnerTargets[B];
    for (Size i = 0; i < owned.size(); ++i) {
      Size v = owned[i];
      if (v >= mN || mS5TargetKeyOwner[v] != B
          ) continue;
      Size source = mS5TargetKeySource[v];
      Size owner = source < mN ? outerOfUF(source) : scNull();
      Size targetOwner = outerOfUF(v);
      if (owner < mNumNodes && owner != targetOwner &&
        mNLabel[owner] == LBL_EVEN && treeLive(owner)) {
        s5TargetKeyReindexPayload(v, owner,
          mS5TargetKeyTargetOwner[v], mS5TargetKeySourceToken[v], true, true);
        ++mCounters.s5TargetKeyOwnerTransfers;
      } else {
        s5TargetKeyRevealFree(v, tNow);
      }
    }
    std::vector<Size>().swap(owned);
  }
  // A FREE blossom can carry a target key before it dissolves. Rebind each
  // such payload to its exposed child interval after UF restoration. FREE
  // children are revealed exactly; non-FREE children cannot own grow keys.
  if (B < mS5TargetKeyTargetOwnerTargets.size()) {
    std::vector<Size>& targeted = mS5TargetKeyTargetOwnerTargets[B];
#ifdef EVENT_S5_TARGET_KEY_TRACE
    if (!targeted.empty())
      std::fprintf(stderr,
        "[S5TK-TARGET-DISSOLVE] scale=%lu B=%lu count=%lu t=%lld\n",
        (unsigned long)mAuditScale, (unsigned long)B,
        (unsigned long)targeted.size(), (long long)tNow);
#endif
    for (Size i = 0; i < targeted.size(); ++i) {
      Size v = targeted[i];
      if (v >= mN || mS5TargetKeyTargetOwner[v] != B
          ) continue;
      Size targetOwner = outerOfUF(v);
      int targetLabel = targetOwner < mNLabel.size()
        ? mNLabel[targetOwner] : LBL_ODD;
      if (targetLabel != LBL_FREE_OF_TREE && !treeLive(targetOwner))
        targetLabel = LBL_FREE_OF_TREE;
      if (targetLabel == LBL_FREE_OF_TREE) {
        s5TargetKeyRevealFree(v, tNow);
      } else {
        const bool oldLive = s5TargetSfmKey(v) != SplitFindmin::INF;
        s5TargetSfmSetKey(v, SplitFindmin::INF);
        mS5TargetKeySource[v] = scNull();
        s5TargetKeyReindexPayload(v, scNull(), targetOwner, scNull(),
                                  oldLive, false);
        mS5TargetKeyEdgeIndex[v] = scNull();
        mS5TargetKeyMatched[v] = 0;
      }
    }
    std::vector<Size>().swap(targeted);
  }
  // The pre-dissolve hook reset new FREE intervals before UF was restored.
  // Reveal them now against the authoritative exposed ownership.
  for (Size i = 0; i < mBCycle[B].size(); ++i) {
    Size child = mBCycle[B][i];
    if (child >= mNLabel.size() || mNLabel[child] != LBL_FREE_OF_TREE) continue;
#ifdef SCALING_LEAF_NOVEC
    std::vector<Size> leaves; scalingLfCollect(child, &leaves);
#else
    const std::vector<Size>& leaves = mNodeLeaves[child];
#endif
    for (Size j = 0; j < leaves.size(); ++j)
      s5TargetKeyRevealFree(leaves[j], tNow);
  }
}

inline void ScalingMatching::s5TargetKeyRevealFree(Size v, LL tNow) {
  if (!mS5TargetKeyActive || v >= mN) return;
  Size V = outerOfUF(v);
  int lblV = mNLabel[V];
  if (lblV != LBL_FREE_OF_TREE && !treeLive(V)) lblV = LBL_FREE_OF_TREE;
  if (lblV != LBL_FREE_OF_TREE) return;
  Val best = SplitFindmin::INF;
  Size bestSource = scNull(), bestOwner = scNull(), bestIndex = scNull();
  for (Size ai = 0; ai < mAdj[v].size(); ++ai) {
    ++mCounters.s5TargetKeyTransitionRevealIncidences;
    Size u = mAdj[v][ai];
    if (u >= mN || (mShellActive && !mInShell[u])) continue;
    Size U = outerOfUF(u);
    if (U == V || mNLabel[U] != LBL_EVEN || !treeLive(U)) continue;
    LL edgeW2 = mUseScaleW ? mW[v][ai] : 2 * mWhat[v][ai];
    LL s = yOfOuter(u, U) + yOfOuter(v, V) + oldCommonZ(u, v) - edgeW2;
    LL gap = 0, absoluteNow = 0, fire = 0;
    if (!s5TargetKeyGap(s, mMate[u] == v, mS5TargetKeyCriterion, &gap) ||
        !checkedAddLL(mS5TargetKeyClockBase, tNow, &absoluteNow) ||
        !checkedAddLL(absoluteNow, gap, &fire)) continue;
    if ((Val)fire < best) {
      best = (Val)fire; bestSource = u; bestOwner = U; bestIndex = ai;
    }
  }
  const Val oldKey = s5TargetSfmKey(v);
  if (best != SplitFindmin::INF && oldKey == SplitFindmin::INF)
    s5TargetKeyRememberTarget(v);
  s5TargetSfmSetKey(v, best);
  mS5TargetKeySource[v] = bestSource;
  mS5TargetKeyEdgeIndex[v] = bestIndex;
  mS5TargetKeyMatched[v] = bestSource < mN && mMate[bestSource] == v;
  Size newSourceToken = scNull();
  if (bestOwner < mNRoot.size() && mNRoot[bestOwner] != scNull()) {
    Size ro = outerOf(mNRoot[bestOwner]);
    if (ro < mBBase.size() && mBBase[ro] < mN) {
      newSourceToken = mBBase[ro];
    }
  }
  s5TargetKeyReindexPayload(v, bestOwner, V, newSourceToken,
                            oldKey != SplitFindmin::INF,
                            best != SplitFindmin::INF);
  ++mCounters.s5TargetKeyTransitionReveals;
#ifdef EVENT_S5_TARGET_KEY_TRACE
  std::fprintf(stderr,
    "[S5TK-REVEAL] v=%lu V=%lu label=%d source=%ld owner=%ld key=%lld t=%lld\n",
    (unsigned long)v, (unsigned long)V, lblV,
    (long)(bestSource == scNull() ? -1 : bestSource),
    (long)(bestOwner == scNull() ? -1 : bestOwner),
    (long long)(best == SplitFindmin::INF ? -1 : best), (long long)tNow);
#endif
#ifdef EVENT_S5_TARGET_KEY_DRIVER_AUDIT
  scheduleS5TargetMinEvent(v);
#endif
}

#ifdef EVENT_S5_TARGET_KEY_DRIVER_AUDIT
inline void ScalingMatching::scheduleS5TargetMinEvent(Size target) {
  if (!mS5TargetKeyActive || target >= mN) return;
  // SearchOne is the overwhelmingly hot path and discards its queue at the
  // first positive adjustment. Keep it on the versionless list-min payload:
  // key/source validation is a complete cancellation boundary for that queue.
  if (mEventOneAdjustment) {
    Val fire = SplitFindmin::INF;
    Size at = scNull();
    s5TargetSfmFindMin(target, &fire, &at);
    target = at;
    if (target >= mN) return;
    Size source = mS5TargetKeySource[target];
    if (fire == SplitFindmin::INF || source >= mN) return;
    LL localFire = 0;
    if (!checkedAddLL((LL)fire, -mS5TargetKeyClockBase, &localFire) ||
        localFire < 0) return;
    schedulePracticalEvent(localFire, (long long)source, (long long)target,
                           10 /*SearchOne list minimum*/, -1);
    if (mLastStatus == STATUS_SUCCESS)
      ++mCounters.s5TargetMinEventsScheduled;
    return;
  }
  scheduleS5TargetPersistentFrontier(target);
}

inline void ScalingMatching::scheduleS5TargetPersistentFrontier(Size target) {
  if (!mS5TargetKeyActive || target >= mN) return;
  // Section 5 exposes key(u) through findmin(list(u)).  Queue the aggregate
  // minimum, not merely the leaf that happened to receive the latest
  // decrease.  A later decrease, split, or merge changes this aggregate and
  // lazily invalidates the older payload at pop.
  Val fire = SplitFindmin::INF;
  Size at = scNull();
  Size listId = scNull();
  s5TargetSfmFindMinWithList(target, &fire, &at, &listId);
  target = at;
  if (target >= mN) return;
  Size source = mS5TargetKeySource[target];
  if (fire == SplitFindmin::INF || source >= mN) return;
  LL localFire = 0;
  if (!checkedAddLL((LL)fire, -mS5TargetKeyClockBase, &localFire) ||
      localFire < 0) return;
#ifdef EVENT_S5_LIST_FRONTIER_DISABLE
  schedulePracticalEvent(localFire, (long long)source, (long long)target,
                         10 /*versionless list minimum*/, -1);
  if (mLastStatus == STATUS_SUCCESS)
    ++mCounters.s5TargetMinEventsScheduled;
  return;
#else
  if (listId == scNull()) return;
  if (listId >= mS5TargetFrontier.size())
    mS5TargetFrontier.resize(listId + 1);
  S5TargetFrontierEntry& frontier = mS5TargetFrontier[listId];
  if (frontier.version != 0 && frontier.fire == fire &&
      frontier.winner == target && frontier.source == source) {
    ++mCounters.s5TargetFrontierCoalesced;
    return;
  }
  if (mS5TargetFrontierClock >=
      std::numeric_limits<unsigned int>::max() - 1U) {
    mLastStatus = STATUS_UNSUPPORTED_OR_OVERFLOWING_INPUT;
    return;
  }
  unsigned int version = mS5TargetFrontierClock++;
  const bool superseded = frontier.version != 0;
  // Kind 10 carries the logical list-frontier generation in owner. Replacing
  // findmin(B) cancels the older generation lazily and leaves exactly one
  // current grow frontier for B.
  schedulePracticalEvent(localFire, (long long)source, (long long)target,
                         10 /*Section 5 non-outer-list minimum*/,
                         (long long)version);
  if (mLastStatus == STATUS_SUCCESS) {
    if (!frontier.touched) {
      frontier.touched = true;
      mS5TargetFrontierTouched.push_back(listId);
    }
    frontier.version = version;
    frontier.fire = fire;
    frontier.winner = target;
    frontier.source = source;
    if (superseded) ++mCounters.s5TargetFrontierSuperseded;
    ++mCounters.s5TargetMinEventsScheduled;
  }
#endif
}
#endif

inline void ScalingMatching::s5TargetKeyRetireTree(Size token, LL tNow) {
  if (!mS5TargetKeyActive || token >= mS5TargetKeyTreeNodes.size()) return;
  // A fresh singleton root's only real vertex is identical to its stable tree
  // token. Initialisation leaves that membership implicit; materialize it only
  // for the small subset of trees that are actually retired. Grown/blossom
  // nodes have already appended their additional leaves through TrackTreeNode.
  if (token < mN &&
      (mS5TargetKeyAllVertexStamp[token] != mS5TargetKeyTouchClock ||
       mS5TargetKeyAllVertexToken[token] != token)) {
    mS5TargetKeyAllVertexStamp[token] = mS5TargetKeyTouchClock;
    mS5TargetKeyAllVertexToken[token] = token;
    s5TargetKeyRememberToken(token);
    mS5TargetKeyTreeAllVertices[token].push_back(token);
  }
  if (++mS5TargetKeyRevealClock == 0) {
    std::fill(mS5TargetKeyRevealStamp.begin(),
              mS5TargetKeyRevealStamp.end(), 0);
    ++mS5TargetKeyRevealClock;
  }
  ++mCounters.s5TargetKeyRetiredTrees;
  const std::vector<Size>& treeVertices = mS5TargetKeyTreeAllVertices[token];
#ifdef EVENT_S5_TARGET_KEY_TRACE
  std::fprintf(stderr,
    "[S5TK-RETIRE] token=%lu vertices=%lu owned=%lu t=%lld\n",
    (unsigned long)token, (unsigned long)treeVertices.size(),
    (unsigned long)mS5TargetKeyTreeOwnedTargets[token].size(),
    (long long)tNow);
#endif
  // The augmentation consumes only its marked contracted path, not every
  // branch of either alternating tree.  Retire the dead tree labels on each
  // wholly-unused outer node before revealing its target keys.  Otherwise an
  // edge that was OUTER->INNER before augmentation remains hidden behind a
  // stale INNER label even though that branch is now available to another
  // surviving root in the same Figure-3 phase.
  if (mS5TargetKeyCriterion == CRIT2 && !mUsedV.empty()) {
    for (Size i = 0; i < treeVertices.size(); ++i) {
      Size v = treeVertices[i];
      if (v >= mN || (v < mUsedV.size() && mUsedV[v])) continue;
      Size V = outerOfUF(v);
      if (V >= mNumNodes || mNLabel[V] == LBL_FREE_OF_TREE || treeLive(V))
        continue;
      // treeVertices is a flat leaf list.  The base is present exactly once,
      // so let it represent the contracted outer node instead of collecting
      // the same blossom leaves once per member.
      if (V >= mBBase.size() || mBBase[V] != v) continue;
#ifdef SCALING_LEAF_NOVEC
      std::vector<Size> leaves; scalingLfCollect(V, &leaves);
#else
      const std::vector<Size>& leaves = mNodeLeaves[V];
#endif
      bool pathNode = false;
      for (Size j = 0; j < leaves.size(); ++j)
        if (leaves[j] < mUsedV.size() && mUsedV[leaves[j]]) {
          pathNode = true; break;
        }
      if (pathNode) continue;
      if (V < mN) lazySetVertexRole(V, LBL_FREE_OF_TREE);
      else lazySetNodeRole(V, LBL_FREE_OF_TREE);
      mNLabel[V] = LBL_FREE_OF_TREE;
      mNParent[V] = scNull(); mNParEdgeU[V] = scNull();
      mNParEdgeV[V] = scNull(); mNRoot[V] = scNull();
    }
  }
  for (Size i = 0; i < treeVertices.size(); ++i) {
    Size v = treeVertices[i];
    if (v >= mN || (v < mUsedV.size() && mUsedV[v]) ||
        mS5TargetKeyRevealStamp[v] == mS5TargetKeyRevealClock)
      continue;
    mS5TargetKeyRevealStamp[v] = mS5TargetKeyRevealClock;
    ++mCounters.s5TargetKeyRetiredVertices;
    s5TargetKeyRevealFree(v, tNow);
  }
  const std::vector<Size>& ownedTargets = mS5TargetKeyTreeOwnedTargets[token];
  const Size ownedCount = ownedTargets.size();
  for (Size i = 0; i < ownedCount; ++i) {
    Size v = ownedTargets[i];
    if (v >= mN || mS5TargetKeyRevealStamp[v] == mS5TargetKeyRevealClock)
      continue;
    if (mS5TargetKeySourceToken[v] != token) continue;
    mS5TargetKeyRevealStamp[v] = mS5TargetKeyRevealClock;
    s5TargetKeyRevealFree(v, tNow);
  }
}

inline void ScalingMatching::initializeS5TargetKeyShadow(int crit) {
#ifdef S5_TARGET_INIT_PROFILE
  const auto startTargetInitPart = []() -> unsigned long long {
    return scalingMonotonicNanoseconds();
  };
  const auto finishTargetInitPart = [](unsigned long long* counter,
                                       unsigned long long start) {
    *counter += scalingMonotonicNanoseconds() - start;
  };
#else
  const auto startTargetInitPart = []() -> unsigned long long { return 0; };
  const auto finishTargetInitPart = [](unsigned long long*,
                                       unsigned long long) {};
#endif
  ++mCounters.s5TargetKeyInitializations;
  mS5TargetKeyCriterion = crit;
  bool reuseTopology = false;
#ifdef EVENT_S5_TARGET_KEY_PERSIST_FRESH
#ifdef EVENT_S5_SHARED_SFM_TOPOLOGY
  reuseTopology = mEventSfmInitialized &&
                  mEventSfmScale == mAuditScale &&
                  mS5TargetKeySource.size() == mN;
#else
  reuseTopology = mS5TargetKeyTopologyReady &&
                  mS5TargetKeyTopologyScale == mAuditScale &&
                  mS5TargetKeySource.size() == mN;
#endif
#endif
  // A fresh Edmonds search owns a fresh event queue. Boundary reuse bypasses
  // this initializer and therefore deliberately preserves live generations.
  unsigned long long initPartStart = startTargetInitPart();
  s5TargetFrontierReset();
  finishTargetInitPart(&mCounters.s5TargetInitFrontierNanoseconds,
                       initPartStart);
  if (reuseTopology) {
    // End the prior forest's ownership intervals without rebuilding the
    // laminar split-findmin topology. Work is proportional to state actually
    // touched by that forest, not to all graph vertices or allocated blossoms.
    initPartStart = startTargetInitPart();
    mCounters.s5TargetInitTargetVisits += mS5TargetKeyTouchedTargets.size();
    for (Size i = 0; i < mS5TargetKeyTouchedTargets.size(); ++i) {
      Size v = mS5TargetKeyTouchedTargets[i];
      if (v >= mN) continue;
      const bool oldLive = s5TargetSfmKey(v) != SplitFindmin::INF;
      s5TargetSfmSetKey(v, SplitFindmin::INF);
      mS5TargetKeySource[v] = scNull();
      s5TargetKeyReindexPayload(v, scNull(), scNull(), scNull(),
                                oldLive, false);
      mS5TargetKeyEdgeIndex[v] = scNull();
      mS5TargetKeyMatched[v] = 0;
    }
    finishTargetInitPart(&mCounters.s5TargetInitTargetsNanoseconds,
                         initPartStart);
    initPartStart = startTargetInitPart();
    mCounters.s5TargetInitOwnerVisits += mS5TargetKeyTouchedOwners.size();
    for (Size i = 0; i < mS5TargetKeyTouchedOwners.size(); ++i) {
      Size owner = mS5TargetKeyTouchedOwners[i];
      if (owner < mS5TargetKeyOwnerTargets.size())
        mS5TargetKeyOwnerTargets[owner].clear();
      if (owner < mS5TargetKeyTargetOwnerTargets.size())
        mS5TargetKeyTargetOwnerTargets[owner].clear();
    }
    finishTargetInitPart(&mCounters.s5TargetInitOwnersNanoseconds,
                         initPartStart);
    initPartStart = startTargetInitPart();
    mCounters.s5TargetInitTokenVisits += mS5TargetKeyTouchedTokens.size();
    for (Size i = 0; i < mS5TargetKeyTouchedTokens.size(); ++i) {
      Size token = mS5TargetKeyTouchedTokens[i];
      if (token >= mN) continue;
      // Deduplication markers are epoch-stamped. Advancing the search epoch
      // below retires every marker at once; only the per-token payload lists
      // themselves need to release their logical contents here.
      mS5TargetKeyTreeNodes[token].clear();
      mS5TargetKeyTreeVertices[token].clear();
      mS5TargetKeyTreeAllVertices[token].clear();
      mS5TargetKeyTreeOwnedTargets[token].clear();
    }
    finishTargetInitPart(&mCounters.s5TargetInitTokensNanoseconds,
                         initPartStart);
    initPartStart = startTargetInitPart();
    mS5TargetKeyTouchedTargets.clear();
    mS5TargetKeyTouchedOwners.clear();
    mS5TargetKeyTouchedTokens.clear();
    if (++mS5TargetKeyTouchClock == 0) {
      std::fill(mS5TargetKeyTargetTouchStamp.begin(),
                mS5TargetKeyTargetTouchStamp.end(), 0);
      std::fill(mS5TargetKeyOwnerTouchStamp.begin(),
                mS5TargetKeyOwnerTouchStamp.end(), 0);
      std::fill(mS5TargetKeyTokenTouchStamp.begin(),
                mS5TargetKeyTokenTouchStamp.end(), 0);
      std::fill(mS5TargetKeyTrackedRootStamp.begin(),
                mS5TargetKeyTrackedRootStamp.end(), 0);
      std::fill(mS5TargetKeyTrackedVertexStamp.begin(),
                mS5TargetKeyTrackedVertexStamp.end(), 0);
      std::fill(mS5TargetKeyAllVertexStamp.begin(),
                mS5TargetKeyAllVertexStamp.end(), 0);
      ++mS5TargetKeyTouchClock;
    }
    if (mS5TargetKeySplitPivots.size() < mZ.size())
      mS5TargetKeySplitPivots.resize(mZ.size());
    if (mS5TargetKeyOwnerTargets.size() < mZ.size())
      mS5TargetKeyOwnerTargets.resize(mZ.size());
    if (mS5TargetKeyTargetOwnerTargets.size() < mZ.size())
      mS5TargetKeyTargetOwnerTargets.resize(mZ.size());
    if (mS5TargetKeyTrackedRoot.size() < mZ.size())
      mS5TargetKeyTrackedRoot.resize(mZ.size(), scNull());
    if (mS5TargetKeyTrackedRootStamp.size() < mZ.size())
      mS5TargetKeyTrackedRootStamp.resize(mZ.size(), 0);
    if (mS5TargetKeyOwnerTouchStamp.size() < mZ.size())
      mS5TargetKeyOwnerTouchStamp.resize(mZ.size(), 0);
    finishTargetInitPart(&mCounters.s5TargetInitResizeNanoseconds,
                         initPartStart);
  } else {
    initPartStart = startTargetInitPart();
#ifndef EVENT_S5_SHARED_SFM_TOPOLOGY
    std::vector<Size> order(mN);
    std::vector<char> listEnd(mN, 1);
    for (Size v = 0; v < mN; ++v) order[v] = v;
    mS5TargetKeySfm.init(order, listEnd);
#else
    // initializeEventStructures() owns the shared topology, but a fresh
    // target-payload epoch must explicitly initialize its independent key2
    // channel.  Primary INNER-key activity and topology reconstruction do not
    // imply anything about these absolute grow-event keys.
    mCounters.sfmGrowResetNodeVisits +=
      mInnerSfm.resetkeys2(SplitFindmin::INF);
#ifdef EVENT_S5_SINGLETON_TARGET_FAST
    std::fill(mS5TargetKeySingletonValue.begin(),
              mS5TargetKeySingletonValue.end(), SplitFindmin::INF);
#endif
#endif
    mS5TargetKeySplitPivots.assign(mZ.size(), std::vector<Size>());
    mS5TargetKeySource.assign(mN, scNull());
    mS5TargetKeyOwner.assign(mN, scNull());
    mS5TargetKeyTargetOwner.assign(mN, scNull());
    mS5TargetKeySourceToken.assign(mN, scNull());
    mS5TargetKeyEdgeIndex.assign(mN, scNull());
    mS5TargetKeyMatched.assign(mN, 0);
    mS5TargetKeyOwnerTargets.assign(mZ.size(), std::vector<Size>());
    mS5TargetKeyTargetOwnerTargets.assign(mZ.size(), std::vector<Size>());
    mS5TargetKeyTreeNodes.assign(mN, std::vector<Size>());
    mS5TargetKeyTreeVertices.assign(mN, std::vector<Size>());
    mS5TargetKeyTreeAllVertices.assign(mN, std::vector<Size>());
    mS5TargetKeyTreeOwnedTargets.assign(mN, std::vector<Size>());
    mS5TargetKeyTrackedRoot.assign(mZ.size(), scNull());
    mS5TargetKeyTrackedVertexToken.assign(mN, scNull());
    mS5TargetKeyAllVertexToken.assign(mN, scNull());
    mS5TargetKeyTrackedRootStamp.assign(mZ.size(), 0);
    mS5TargetKeyTrackedVertexStamp.assign(mN, 0);
    mS5TargetKeyAllVertexStamp.assign(mN, 0);
    mS5TargetKeyRevealStamp.assign(mN, 0);
    mS5TargetKeyRevealClock = 0;
    mS5TargetKeyTouchedTargets.clear();
    mS5TargetKeyTouchedOwners.clear();
    mS5TargetKeyTouchedTokens.clear();
    mS5TargetKeyTargetTouchStamp.assign(mN, 0);
    mS5TargetKeyOwnerTouchStamp.assign(mZ.size(), 0);
    mS5TargetKeyTokenTouchStamp.assign(mN, 0);
    mS5TargetKeyTouchClock = 1;
    mS5TargetKeyTopologyReady = false;
    mS5TargetKeyTopologyScale = mAuditScale;
    // Build current Omega bottom-up. Omega' remains a distinct inherited
    // forest and is never merged into Section 5's current-blossom lists.
#ifndef EVENT_S5_SHARED_SFM_TOPOLOGY
    for (Size b = mN; b < mNumNodes; ++b)
      if (b < mBCycle.size() && !mBCycle[b].empty() &&
          (b >= mBActive.size() || mBActive[b]) &&
          !(mOldForestSeparated && b >= mAuditOldBegin && b < mAuditOldEnd))
        s5TargetKeyOnShrink(b);
#endif
    mS5TargetKeyTopologyReady = true;
    finishTargetInitPart(&mCounters.s5TargetInitFreshNanoseconds,
                         initPartStart);
  }
  // In shared mode initializeEventStructures() owns the topology epoch. The
  // target payload layer may reuse its touched-state storage across a scale
  // change, but its shrink/dissolve hooks must follow the new authoritative
  // topology immediately.
  mS5TargetKeyTopologyScale = mAuditScale;
  mS5TargetKeyTopologyReady = true;
  mS5TargetKeyClockBase = 0;
  mS5TargetKeyActive = true;
  initPartStart = startTargetInitPart();
  mCounters.s5TargetInitRootVisits += mEventInitialRoots.size();
  for (Size i = 0; i < mEventInitialRoots.size(); ++i) {
    Size root = mEventInitialRoots[i];
    // Singleton membership is the identity token itself and is materialized
    // lazily by s5TargetKeyRetireTree. Non-singleton roots still need their
    // complete leaf set registered now.
    if (root < mN && root < mNRoot.size() && mNRoot[root] == root &&
        mBParent[root] == scNull()) continue;
    s5TargetKeyTrackTreeNode(root);
  }
  finishTargetInitPart(&mCounters.s5TargetInitRootsNanoseconds,
                       initPartStart);
#ifdef EVENT_S5_TARGET_KEY_TRACE
  for (Size v = 0; v < mN; ++v)
    if (mS5TargetKeySource[v] == scNull() &&
        s5TargetSfmKey(v) != SplitFindmin::INF)
      std::fprintf(stderr,
        "[S5TK-INIT-STALE] scale=%lu crit=%d v=%lu key=%lld reuse=%d touched=%lu\n",
        (unsigned long)mAuditScale, crit, (unsigned long)v,
        (long long)s5TargetSfmKey(v), reuseTopology ? 1 : 0,
        (unsigned long)mS5TargetKeyTouchedTargets.size());
#endif
}

inline void ScalingMatching::auditS5TargetKeys(int crit, LL tNow) {
  if (!mS5TargetKeyActive || crit != mS5TargetKeyCriterion) return;
  std::vector<Val> exact(mN, SplitFindmin::INF);
  std::vector<Size> exactSource(mN, scNull());
  std::vector<Size> exactOwner(mN, scNull());
  for (Size v = 0; v < mN; ++v) {
    if (mShellActive && !mInShell[v]) continue;
    Size V = outerOfUF(v);
    int lblV = mNLabel[V];
    if (lblV != LBL_FREE_OF_TREE && !treeLive(V)) lblV = LBL_FREE_OF_TREE;
    if (lblV != LBL_FREE_OF_TREE) continue;
    for (Size ai = 0; ai < mAdj[v].size(); ++ai) {
      Size u = mAdj[v][ai];
      if (u >= mN || (mShellActive && !mInShell[u])) continue;
      Size U = outerOfUF(u);
      if (U == V || mNLabel[U] != LBL_EVEN || !treeLive(U)) continue;
      LL edgeW2 = mUseScaleW ? mW[v][ai] : 2 * mWhat[v][ai];
      LL s = yOfOuter(u, U) + yOfOuter(v, V) + oldCommonZ(u, v) - edgeW2;
      LL gap = 0;
      LL absoluteNow = 0, fire = 0;
      if (s5TargetKeyGap(s, mMate[u] == v, crit, &gap) &&
          checkedAddLL(mS5TargetKeyClockBase, tNow, &absoluteNow) &&
          checkedAddLL(absoluteNow, gap, &fire) &&
          (Val)fire < exact[v]) {
        exact[v] = (Val)fire;
        exactSource[v] = u;
        exactOwner[v] = U;
      }
    }
    ++mCounters.s5TargetKeyVertexChecks;
    Val retained = s5TargetSfmKey(v);
    bool invalidSource = false;
    if (retained != SplitFindmin::INF) {
      Size src = mS5TargetKeySource[v];
      Size owner = mS5TargetKeyOwner[v];
      Size targetOwner = mS5TargetKeyTargetOwner[v];
      Size currentOwner = src < mN ? outerOfUF(src) : scNull();
      invalidSource = src >= mN || owner >= mNumNodes ||
          currentOwner != owner || currentOwner == V ||
          targetOwner != V ||
          (bool)mS5TargetKeyMatched[v] != (mMate[src] == v) ||
          mNLabel[currentOwner] != LBL_EVEN || !treeLive(currentOwner);
      // A search may stop with future queue entries belonging to a forest it
      // has just consumed.  If the independent scan also finds no live source,
      // this is an expired payload, not a representation miss.  Clear it in
      // the audit view without an incidence repair; the next search performs
      // the same touched-key reset before activating its new forest.
      if (invalidSource && exact[v] == SplitFindmin::INF) {
        const bool oldLive = retained != SplitFindmin::INF;
        s5TargetSfmSetKey(v, SplitFindmin::INF);
        mS5TargetKeySource[v] = scNull();
        s5TargetKeyReindexPayload(v, scNull(), V, scNull(),
                                  oldLive, false);
        mS5TargetKeyEdgeIndex[v] = scNull();
        mS5TargetKeyMatched[v] = 0;
        retained = SplitFindmin::INF;
        invalidSource = false;
      }
      if (invalidSource) {
        ++mCounters.s5TargetKeyInvalidSource;
#ifdef EVENT_S5_TARGET_KEY_TRACE
        std::fprintf(stderr,
          "[S5TK-INVALID] scale=%lu v=%lu V=%lu src=%lu owner=%lu nowOwner=%lu "
          "targetOwner=%lu ownerLabel=%d ownerLive=%d oldMatched=%d matchedNow=%d "
          "crit=%d retained=%lld exact=%lld exactSource=%lu\n",
          (unsigned long)mAuditScale, (unsigned long)v, (unsigned long)V,
          (unsigned long)src, (unsigned long)owner,
          (unsigned long)currentOwner, (unsigned long)targetOwner,
          currentOwner < mNLabel.size() ? mNLabel[currentOwner] : -9,
          currentOwner < mNumNodes && treeLive(currentOwner) ? 1 : 0,
          (int)mS5TargetKeyMatched[v],
          src < mN && mMate[src] == v ? 1 : 0, crit,
          (long long)retained, (long long)exact[v],
          (unsigned long)exactSource[v]);
#endif
        if (src >= mN || owner >= mNumNodes || currentOwner >= mNumNodes)
          ++mCounters.s5TargetKeyInvalidRange;
        else if (currentOwner != owner)
          ++mCounters.s5TargetKeyInvalidOwnerChange;
        else if (currentOwner == V)
          ++mCounters.s5TargetKeyInvalidInternal;
        else if (targetOwner != V)
          ++mCounters.s5TargetKeyInvalidTargetInterval;
        else if ((bool)mS5TargetKeyMatched[v] != (mMate[src] == v))
          ++mCounters.s5TargetKeyInvalidMatchingRole;
        else
          ++mCounters.s5TargetKeyInvalidSourceRole;
        ++mCounters.s5TargetKeyRepairs;
        mCounters.s5TargetKeyRepairIncidences += mAdj[v].size();
        // Conservative reveal-next: the exact scan above supplies the next
        // live owner. This repair exists only in explicit verification builds;
        // canonical production scheduling never executes the scan audit.
        const bool oldLive = retained != SplitFindmin::INF;
        s5TargetSfmSetKey(v, exact[v]);
        if (exact[v] != SplitFindmin::INF) s5TargetKeyRememberTarget(v);
        mS5TargetKeySource[v] = exactSource[v];
        mS5TargetKeyEdgeIndex[v] = scNull();
        mS5TargetKeyMatched[v] = exactSource[v] < mN &&
                                 mMate[exactSource[v]] == v;
        Size sourceToken = scNull();
        if (exactOwner[v] < mNRoot.size() &&
            mNRoot[exactOwner[v]] != scNull()) {
          Size ro = outerOf(mNRoot[exactOwner[v]]);
          if (ro < mBBase.size()) sourceToken = mBBase[ro];
        }
        s5TargetKeyReindexPayload(v, exactOwner[v], V, sourceToken,
                                  oldLive, exact[v] != SplitFindmin::INF);
        retained = exact[v];
      }
    }
    if (retained > exact[v]) {
#ifdef EVENT_S5_TARGET_KEY_TRACE
        std::fprintf(stderr,
          "[S5TK-MISSING] scale=%lu v=%lu V=%lu retained=%lld exact=%lld "
          "label=%d root=%lu live=%d source=%lu owner=%lu tracked=%lu\n",
          (unsigned long)mAuditScale, (unsigned long)v, (unsigned long)V,
          (long long)retained, (long long)exact[v],
          V < mNLabel.size() ? mNLabel[V] : -9,
          (unsigned long)(V < mNRoot.size() ? mNRoot[V] : scNull()),
          V < mNumNodes && treeLive(V) ? 1 : 0,
          (unsigned long)exactSource[v], (unsigned long)exactOwner[v],
          (unsigned long)(V < mS5TargetKeyTrackedRoot.size() ?
                          mS5TargetKeyTrackedRoot[V] : scNull()));
#endif
      // A same-time augmentation can retire a tree, turning vertices that
      // were not FREE when source rows were offered into FREE targets. Reveal
      // their first key conservatively here; the driving implementation will
      // invoke this repair from the tree-retirement transition itself.
      ++mCounters.s5TargetKeyRepairs;
      ++mCounters.s5TargetKeyMissingRepairs;
      mCounters.s5TargetKeyRepairIncidences += mAdj[v].size();
      const bool oldLive = retained != SplitFindmin::INF;
      s5TargetSfmSetKey(v, exact[v]);
      if (exact[v] != SplitFindmin::INF) s5TargetKeyRememberTarget(v);
      mS5TargetKeySource[v] = exactSource[v];
      mS5TargetKeyEdgeIndex[v] = scNull();
      mS5TargetKeyMatched[v] = exactSource[v] < mN &&
                               mMate[exactSource[v]] == v;
      Size sourceToken = scNull();
      if (exactOwner[v] < mNRoot.size() &&
          mNRoot[exactOwner[v]] != scNull()) {
        Size ro = outerOf(mNRoot[exactOwner[v]]);
        if (ro < mBBase.size()) sourceToken = mBBase[ro];
      }
      s5TargetKeyReindexPayload(v, exactOwner[v], V, sourceToken,
                                oldLive, exact[v] != SplitFindmin::INF);
      retained = exact[v];
    }
    if (retained == exact[v]) ++mCounters.s5TargetKeyVertexExact;
    else if (retained < exact[v]) ++mCounters.s5TargetKeyVertexTooLow;
    else ++mCounters.s5TargetKeyVertexTooHigh;
  }

  // Verify the actual list aggregate, not merely each leaf key. This catches
  // incorrect split/merge boundaries independently of edge ownership errors.
  for (Size B = 0; B < mNumNodes; ++B) {
    if (B >= mBParent.size() || mBParent[B] != scNull() ||
        B >= mBBase.size() || mBBase[B] >= mN) continue;
    if (B >= mN && B < mBActive.size() && !mBActive[B]) continue;
    if (mOldForestSeparated && B >= mAuditOldBegin && B < mAuditOldEnd) continue;
    int lblB = mNLabel[B];
    if (lblB != LBL_FREE_OF_TREE && !treeLive(B)) lblB = LBL_FREE_OF_TREE;
    if (lblB != LBL_FREE_OF_TREE) continue;
    std::vector<Size> leaves;
    s5TargetSfmSublistVerts(mBBase[B], leaves);
    Val want = SplitFindmin::INF;
    for (Size i = 0; i < leaves.size(); ++i)
      if (leaves[i] < exact.size() && exact[leaves[i]] < want)
        want = exact[leaves[i]];
    Val got = SplitFindmin::INF; Size at = scNull();
    s5TargetSfmFindMin(mBBase[B], &got, &at);
    ++mCounters.s5TargetKeyListChecks;
    if (got == want) ++mCounters.s5TargetKeyListExact;
    else if (got < want) ++mCounters.s5TargetKeyListTooLow;
    else ++mCounters.s5TargetKeyListTooHigh;
  }
}
#endif

#ifdef EVENT_SFM_TOPOLOGY_AUDIT
inline bool ScalingMatching::sfmShadowListMatches(Size B) const {
  if (B >= mNumNodes || B >= mBBase.size() || mBBase[B] >= mN) return false;
  std::vector<Size> represented;
  mInnerSfm.sublistVerts(mBBase[B], represented);
#ifdef SCALING_LEAF_NOVEC
  std::vector<Size> expected;
  scalingLfCollect(B, &expected);
#else
  std::vector<Size> expected = mNodeLeaves[B];
#endif
  std::sort(represented.begin(), represented.end());
  std::sort(expected.begin(), expected.end());
  represented.erase(std::unique(represented.begin(), represented.end()),
                    represented.end());
  expected.erase(std::unique(expected.begin(), expected.end()), expected.end());
  return represented == expected;
}
#endif

#ifdef EVENT_SFM_SHADOW
inline void ScalingMatching::sfmShadowBeginInner(Size B) {
#if defined(EVENT_SFM_BOUNDARY_FRONTIER) && defined(EVENT_SFM_REUSE) && \
    defined(EVENT_SFM_LAZY_ACTIVE)
  // The laminar-list topology persists independently of an Edmonds forest.
  // Activate per-search keys only when a real blossom enters an INNER
  // interval; vertex GROWs never need a split-findmin dissolution query.
  if (!mEventStructuresActive && B >= mN && mEventSfmInitialized)
    mEventStructuresActive = true;
#endif
  if (!mEventStructuresActive || B >= mNumNodes) return;
#ifdef EVENT_SFM_PRESEED_FREE
  // Ordinary OUTER-row scheduling has already maintained keys for FREE
  // targets. Carry them into the INNER interval; lazy validation at split
  // removes sources whose role/owner changed in the meantime.
  (void)B;
  return;
#else
#ifdef SCALING_LEAF_NOVEC
  std::vector<Size> leafScratch;
  scalingLfCollect(B, &leafScratch);
  const std::vector<Size>& leaves = leafScratch;
#else
  const std::vector<Size>& leaves = mNodeLeaves[B];
#endif
  // A new INNER interval invalidates every absolute minimum retained from an
  // earlier interval. Seed the current OUTER boundary once; subsequent newly
  // OUTER rows feed ordinary decrease-key updates through scheduleEdgeFromEven.
  for (Size li = 0; li < leaves.size(); ++li) {
    Size v = leaves[li];
    if (v >= mN) continue;
    mInnerSfm.setkey(v, SplitFindmin::INF);
    mSfmKeySource[v] = scNull();
    mSfmKeyOwner[v] = scNull();
    // v remains in the same INNER root throughout this complete row.
    // Materializing its lazy dual per incident edge was pure repeated work.
    const LL targetY = yOfOuter(v, B);
    for (Size ai = 0; ai < mAdj[v].size(); ++ai) {
      Size u = mAdj[v][ai];
      if (u >= mN) continue;
      Size U = outerOfUF(u);
      if (U == B || mNLabel[U] != LBL_EVEN || !treeLive(U)) continue;
      LL edgeW2 = mUseScaleW ? mW[v][ai] : 2 * mWhat[v][ai];
      LL s = yOfOuter(u, U) + targetY + oldCommonZ(u, v) - edgeW2;
      if ((Val)s < mInnerSfm.key(v)) {
        mInnerSfm.decreasekey(v, (Val)s);
        mSfmKeySource[v] = u;
        mSfmKeyOwner[v] = U;
      }
    }
  }
#endif
}
#endif

inline void ScalingMatching::resetPracticalEventQueue() {
  mEventEarliestDeferredEdge = std::numeric_limits<LL>::max();
  mEventDeferredEdgeEvents.clear();
  if (mEventUseBinary) mUnboundedPQ.clear();
  else mPQ.reset();
}

inline void ScalingMatching::schedulePracticalEvent(LL time, long long a,
                                                     long long b, int kind,
                                                     long long owner) {
  const bool edgeEvent = kind == 0 || kind == 1 || kind == 5 || kind == 6 ||
                         kind == 7 || kind == 8 || kind == 10;
  if (mEventDualBudget >= 0) {
    LL cumulative = 0;
    if (!checkedAddLL(mEventEpochDualBase, time, &cumulative)) {
      mLastStatus = STATUS_UNSUPPORTED_OR_OVERFLOWING_INPUT;
      return;
    }
    // Events strictly beyond a bounded request cannot be observed. Events at
    // the boundary remain scheduled so zero-z dissolutions can be processed;
    // edge events there are deferred by the SearchOne boundary rule.
    if (cumulative > mEventDualBudget) {
#ifdef EVENT_OFFER_AUDIT
      if (edgeEvent) ++mCounters.eventEdgeBudgetDrops;
#endif
      return;
    }
  }
#ifdef EVENT_OWNER_MIN_AUDIT
  // Kinds 0/1 and their indexed 5/6 forms are ordinary edge events whose
  // owner is an outer node. Matched kinds 7/8 use a different payload contract
  // and require their own cancellable ownership representation.
  const bool ordinaryEdge = kind == 0 || kind == 1 || kind == 5 || kind == 6;
  if (mEventOneAdjustment && ordinaryEdge && time > mTNow && owner >= 0 &&
      (Size)owner < mEventOwnerMinEpoch.size()) {
    const Size o = (Size)owner;
    BucketPQ::Event offered;
    offered.time = time; offered.a = a; offered.b = b;
    offered.kind = kind; offered.owner = owner;
    if (mEventOwnerMinEpoch[o] != mEventOwnerMinClock) {
      mEventOwnerMinEpoch[o] = mEventOwnerMinClock;
      mEventOwnerMinTime[o] = time;
      mEventOwnerMinTies[o] = 1;
      ++mEventOwnerMinVersion[o];
      if (mEventOwnerMinVersion[o] == 0) ++mEventOwnerMinVersion[o];
      mEventOwnerMinPrimary[o] = offered;
      mEventOwnerMinTouched.push_back(o);
      ++mCounters.eventOwnerMinOwners;
      ++mCounters.eventOwnerMinPayloads;
      ++mEventOwnerMinCurrentPayloads;
    } else if (time < mEventOwnerMinTime[o]) {
      mCounters.eventOwnerMinPayloads -= mEventOwnerMinTies[o];
      mEventOwnerMinCurrentPayloads -= mEventOwnerMinTies[o];
      mEventOwnerMinTime[o] = time;
      mEventOwnerMinTies[o] = 1;
      ++mEventOwnerMinVersion[o];
      if (mEventOwnerMinVersion[o] == 0) ++mEventOwnerMinVersion[o];
      mEventOwnerMinPrimary[o] = offered;
      ++mCounters.eventOwnerMinPayloads;
      ++mEventOwnerMinCurrentPayloads;
      ++mCounters.eventOwnerMinReplacements;
    } else if (time == mEventOwnerMinTime[o]) {
      ++mEventOwnerMinTies[o];
      EventOwnerMinTie tie;
      tie.owner = o; tie.version = mEventOwnerMinVersion[o];
      tie.event = offered;
      mEventOwnerMinTieLog.push_back(tie);
      ++mCounters.eventOwnerMinPayloads;
      ++mEventOwnerMinCurrentPayloads;
    } else {
      ++mCounters.eventOwnerMinLaterRejected;
    }
    if (mEventOwnerMinCurrentPayloads >
        mCounters.eventOwnerMinMaxPayloadsPerSearch)
      mCounters.eventOwnerMinMaxPayloadsPerSearch =
        mEventOwnerMinCurrentPayloads;
  }
  if (mEventOneAdjustment && time > mTNow) {
    if (kind == 7 || kind == 8 || kind == 2) {
      BucketPQ::Event future;
      future.time = time; future.a = a; future.b = b;
      future.kind = kind; future.owner = owner;
      if (kind == 7 || kind == 8) {
        mEventOwnerMatchedFuture.push_back(future);
        if (kind == 7) ++mCounters.eventOwnerFutureMatchedGrow;
        else ++mCounters.eventOwnerFutureMatchedBlossom;
        if (mEventOwnerMatchedFuture.size() >
            mCounters.eventOwnerFutureMaxMatched)
          mCounters.eventOwnerFutureMaxMatched =
            mEventOwnerMatchedFuture.size();
      } else {
        mEventOwnerControlFuture.push_back(future);
        ++mCounters.eventOwnerFutureDissolve;
        if (mEventOwnerControlFuture.size() >
            mCounters.eventOwnerFutureMaxControl)
          mCounters.eventOwnerFutureMaxControl =
            mEventOwnerControlFuture.size();
      }
    }
    else if (!ordinaryEdge && kind != 4)
      ++mCounters.eventOwnerFutureOtherControl;
  }
#endif
  if (mEventUseBinary) {
    mUnboundedPQ.schedule(time, a, b, kind, owner);
    ++mCounters.heapQueueInsertions;
  } else {
    // SearchOne consumes zero-time edge transitions, then stops immediately
    // after its first positive dual adjustment. Only the earliest positive
    // edge time can therefore affect this request; one marker preserves that
    // clock boundary without storing every deferred edge payload.
    if (mEventOneAdjustment && edgeEvent && time > mTNow) {
      if (!mEventBoundaryReuse) {
        if (time >= mEventEarliestDeferredEdge) {
#ifdef EVENT_OFFER_AUDIT
          ++mCounters.eventEdgeBoundaryCoalesced;
#endif
          return;
        }
        mEventEarliestDeferredEdge = time;
        a = b = owner = -1;
        kind = 4;
      } else {
      BucketPQ::Event deferred;
      deferred.time = time; deferred.a = a; deferred.b = b;
      deferred.kind = kind; deferred.owner = owner;
      if (time < mEventEarliestDeferredEdge) {
        if (mEventDeferringSeedBoundary &&
            mEventEarliestDeferredEdge != std::numeric_limits<LL>::max())
#ifdef EVENT_OFFER_AUDIT
          ++mCounters.eventEdgeSeedMarkersAvoided;
#else
          (void)0;
#endif
        mEventEarliestDeferredEdge = time;
        mEventDeferredEdgeEvents.clear();
        mEventDeferredEdgeEvents.push_back(deferred);
        a = b = owner = -1;
        kind = 4;
      } else if (time == mEventEarliestDeferredEdge) {
        mEventDeferredEdgeEvents.push_back(deferred);
#ifdef EVENT_OFFER_AUDIT
        ++mCounters.eventEdgeBoundaryCoalesced;
#endif
        return;
      } else {
#ifdef EVENT_OFFER_AUDIT
        ++mCounters.eventEdgeBoundaryCoalesced;
#endif
        return;
      }
      if (mEventDeferringSeedBoundary) return;
      }
    }
    if (!mPQ.schedule(time, a, b, kind, owner)) {
      mLastStatus = STATUS_UNSUPPORTED_OR_OVERFLOWING_INPUT;
      return;
    }
    ++mCounters.bucketQueueInsertions;
  }
  ++mCounters.queueInsertions;
}

inline bool ScalingMatching::popPracticalEvent(BucketPQ::Event* event) {
  if (!event) return false;
  if (mEventUseBinary) {
    BinaryEventHeap::Event binary;
    if (!mUnboundedPQ.pop(&binary)) return false;
    event->time = binary.time; event->a = binary.a;
    event->b = binary.b; event->kind = binary.kind;
    event->owner = binary.owner;
    ++mCounters.heapQueuePops;
  } else {
    if (!mPQ.advanceToNext() || !mPQ.popAtNow(event)) return false;
    ++mCounters.bucketQueuePops;
  }
  ++mCounters.queuePops;
  return true;
}

inline void ScalingMatching::scheduleNodeEdges(Size E, int crit, LL tNow,
                                                bool dedupFreshEven) {
#ifdef SCALING_DETAILED_SCHEDULE_TIMING
  ScalingNanosecondGuard schedulingTimer(&mCounters.eventSchedulingNanoseconds);
#endif
#ifdef SCALING_LEAF_NOVEC
  Size leafEnd = E < mLfTailX.size() ? mLfTailX[E] : scNull();
  Size ui = 0;
  for (Size u = E < mLfHeadX.size() ? mLfHeadX[E] : scNull(), next = scNull();
       u != scNull(); u = next, ++ui) {
    next = (u == leafEnd) ? scNull() : mLfNextX[u];
#else
  std::vector<Size>& Everts = mNodeLeaves[E];
  for (Size ui = 0; ui < Everts.size(); ++ui) {
    Size u = Everts[ui];
#endif
    // Keep the source representative path compressed once per adjacency row.
    // slackZW used to do this again for every incident edge; omitting it
    // entirely lets deep union chains accumulate and slows later outerOfUF
    // queries, while one find here retains that benefit at row granularity.
    if (!mUf.empty() && u < mUf.size()) (void)ufFind(u);
    // E and the lazy clock are constant across this complete adjacency row.
    // Materializing y(u) per incident edge was pure repeated timestamp work.
    const LL sourceY = yOfOuter(u, E);
    // These rows and the weight mode are immutable for the complete scan.
    // Keep their storage and scale multiplier out of the per-edge path; this
    // loop accounts for hundreds of millions of incidence visits on the
    // large event benchmarks.
    const std::vector<Size>& adjRow = mAdj[u];
    const Size degree = adjRow.size();
    const Size* targets = degree ? adjRow.data() : NULL;
    const std::vector<LL>& weightRow = mUseScaleW ? mW[u] : mWhat[u];
    const LL* weights = degree ? weightRow.data() : NULL;
    const LL weightMultiplier = mUseScaleW ? 1 : 2;
#ifdef EVENT_FRESH_ROW_FRONTIER_AUDIT
    if (dedupFreshEven && !mShellActive && !mEventStructuresActive && u < mN) {
      if (mFreshRowAudit.size() < mN) mFreshRowAudit.resize(mN);
      unsigned long long digest = 1469598103934665603ULL;
      const auto mix = [&digest](unsigned long long x) {
        digest ^= x + 0x9e3779b97f4a7c15ULL + (digest << 6) + (digest >> 2);
        digest *= 1099511628211ULL;
      };
      mix((unsigned long long)mAuditScale);
      mix((unsigned long long)(unsigned int)crit);
      mix((unsigned long long)E);
      mix((unsigned long long)degree);
      mix(mUseScaleW ? 1ULL : 0ULL);
      for (Size ak = 0; ak < degree; ++ak) {
        Size aw = targets[ak];
        Size aW = outerOfUF(aw);
        int alabel = mNLabel[aW];
        bool alive = treeLive(aW);
        if (alabel != LBL_FREE_OF_TREE && !alive)
          alabel = LBL_FREE_OF_TREE;
        LL aw2 = weightMultiplier * weights[ak];
        mix((unsigned long long)ak);
        mix((unsigned long long)aw);
        mix((unsigned long long)aW);
        mix((unsigned long long)(unsigned int)(alabel + 2));
        mix(alive ? 1ULL : 0ULL);
        mix((unsigned long long)yOfOuter(aw, aW));
        mix((unsigned long long)oldCommonZ(u, aw));
        mix((unsigned long long)aw2);
        mix((u < mMate.size() && mMate[u] == aw) ? 1ULL : 0ULL);
        mix((aw < mMate.size()) ? (unsigned long long)mMate[aw]
                               : (unsigned long long)scNull());
      }
      FreshRowAuditEntry& prior = mFreshRowAudit[u];
      ++mCounters.freshRowAuditRows;
      if (prior.valid && prior.scale == mAuditScale &&
          prior.criterion == crit && prior.owner == E &&
          prior.degree == degree && prior.digest == digest &&
          prior.sourceY == sourceY && prior.time == tNow) {
        ++mCounters.freshRowAuditReusable;
        mCounters.freshRowAuditReusableIncidences += degree;
      } else if (prior.valid) {
        if (prior.scale != mAuditScale) ++mCounters.freshRowAuditScaleMisses;
        else if (prior.criterion != crit || prior.owner != E ||
                 prior.degree != degree)
          ++mCounters.freshRowAuditContextMisses;
        else if (prior.digest != digest)
          ++mCounters.freshRowAuditSemanticMisses;
        else
          ++mCounters.freshRowAuditSourceDualMisses;
      }
      prior.valid = true;
      prior.scale = mAuditScale;
      prior.criterion = crit;
      prior.owner = E;
      prior.sourceY = sourceY;
      prior.time = tNow;
      prior.degree = degree;
      prior.digest = digest;
    }
#endif
    for (Size k = 0; k < degree; ++k) {
      ++mCounters.edgeScans;
#ifdef SCPHASE
      if ((mCounters.edgeScans & (((unsigned long long)1 << 20) - 1)) == 0) {
        std::fprintf(stderr,
          "[SCPHASE]      seed edge-scans=%llu node=%lu leaf=%lu/%lu clock=%.3fs\n",
          mCounters.edgeScans, (unsigned long)E, (unsigned long)ui,
          (unsigned long)
#ifdef SCALING_LEAF_NOVEC
          scalingLfCount(E),
#else
          Everts.size(),
#endif
          (double)clock()/CLOCKS_PER_SEC);
        std::fflush(stderr);
      }
#endif
      Size w = targets[k];
      LL edgeW2 = weightMultiplier * weights[k];
      scheduleEdgeFromEven(u, w, E, crit, tNow, edgeW2,
                           dedupFreshEven, k, sourceY);
      if (mLastStatus != STATUS_SUCCESS) return;
    }
  }
}

// Section 5.2 schedule(u), matched-edge half. W is already a live INNER root
// node with a valid parent/root path. Its base's mate can become OUTER now, be
// scheduled for a later GROW, or already be INNER (cancel GROW and schedule the
// matched blossom event). This same operation is required after ordinary GROW
// and whenever dissolution exposes an INNER child.
inline bool ScalingMatching::scheduleMatchedEdgeFromInner(
    Size W, int crit, LL tNow) {
  if (W >= mNumNodes || mNLabel[W] != LBL_ODD || !treeLive(W)) return false;
  Size bw = mBBase[W];
  if (bw >= mN) return false;
  Size m = mMate[bw];
  if (m == scNull() || m >= mN) return false;
  Size M = outerOfUF(m);
  if (M == W) return true;
  LL matchedSlack = yOfOuter(bw, W) + yOfOuter(m, M) +
                    oldCommonZ(bw, m) - matchedEdgeW2(bw, m);
  bool matchedEligible = crit == CRIT1 ? matchedSlack == 0 :
                         crit == CRIT2 ? matchedSlack == 0 :
                         (matchedSlack == 0 || matchedSlack == -2);
#ifdef MATCHED_GROW_AUDIT
  if (crit == CRIT2 || crit == CRIT3) {
    ++mCounters.matchedGrowTransitions;
    if (!matchedEligible) {
      ++mCounters.matchedGrowIneligible;
      LL auditTarget = (crit == CRIT3 && matchedSlack <= -2) ? -2 : 0;
      LL auditDelay = auditTarget - matchedSlack;
      if (auditDelay == 1) ++mCounters.matchedGrowDelayOne;
      else if (auditDelay == 2) ++mCounters.matchedGrowDelayTwo;
      else ++mCounters.matchedGrowDelayOther;
    }
  }
#endif
  if (matchedEligible && mNLabel[M] == LBL_FREE_OF_TREE) {
    if (M < mN) lazySetVertexRole(M, LBL_EVEN);
    else lazySetNodeRole(M, LBL_EVEN);
    mNLabel[M] = LBL_EVEN; mNParent[M] = W;
    rememberCompactEventLabel(M);
    mNParEdgeU[M] = m; mNParEdgeV[M] = bw; mNRoot[M] = mNRoot[W];
#ifdef EVENT_S5_TARGET_KEY_SHADOW
    s5TargetKeyTrackTreeNode(M);
#endif
    scheduleNodeEdges(M, crit, tNow);
    return mLastStatus == STATUS_SUCCESS;
  }
  if (matchedEligible && mNLabel[M] == LBL_ODD) {
    schedulePracticalEvent(tNow, (long long)bw, (long long)m,
                           8 /*matched inner-inner blossom*/, 0);
    return mLastStatus == STATUS_SUCCESS;
  }
  if (matchedEligible)
    return mNLabel[M] == LBL_EVEN;

  // The mate may already have become OUTER through another alternating-tree
  // path before this matched edge reaches its criterion target. In that case
  // there is no future GROW transition to schedule: the requested destination
  // role already holds. Treat the delayed matched edge as cancelled, just as
  // Section 5.2 cancels a stale grow whose structural effect was superseded.
  if (mNLabel[M] == LBL_EVEN) return true;

  LL target = (crit == CRIT3 && matchedSlack <= -2) ? -2 : 0;
  LL gap = target - matchedSlack, fire = 0;
  if (gap < 0) return false;
  if (mNLabel[M] == LBL_FREE_OF_TREE) {
    if (!checkedAddLL(tNow, gap, &fire)) {
      mLastStatus = STATUS_UNSUPPORTED_OR_OVERFLOWING_INPUT;
      return false;
    }
    schedulePracticalEvent(fire, (long long)bw, (long long)m,
                           7 /*matched grow*/, (long long)W);
    return mLastStatus == STATUS_SUCCESS;
  }
  if (mNLabel[M] == LBL_ODD) {
    if (!checkedAddLL(tNow, (gap + 1) / 2, &fire)) {
      mLastStatus = STATUS_UNSUPPORTED_OR_OVERFLOWING_INPUT;
      return false;
    }
    schedulePracticalEvent(fire, (long long)bw, (long long)m,
                           8 /*matched inner-inner blossom*/,
                           gap & 1 /*half-integral correction*/);
    return mLastStatus == STATUS_SUCCESS;
  }
  return false;
}

inline void ScalingMatching::scheduleEdgeFromEven(
    Size u, Size w, Size E, int crit, LL tNow, LL edgeW2,
    bool dedupFreshEven, Size edgeIndex, LL sourceY) {
#ifdef EVENT_OFFER_AUDIT
  ++mCounters.eventEdgeOffers;
  if (dedupFreshEven) ++mCounters.eventEdgeFreshOffers;
  else ++mCounters.eventEdgeTransitionOffers;
#endif
  if (mShellActive && (w >= mN || !mInShell[w])) {
#ifdef EVENT_OFFER_AUDIT
    ++mCounters.eventEdgeShellRejects;
#endif
    return;
  }
  // Event initialization can visit every edge while the current forest is
  // deeply nested.  Resolve through the maintained outer-blossom union-find;
  // outerOfUF retains the exact laminar walk as a representation slow path.
  Size W = outerOfUF(w);
  if (W == E) {
#ifdef EVENT_OFFER_AUDIT
    ++mCounters.eventEdgeInternalRejects;
#endif
    return;
  }
  int lblW = mNLabel[W];
  if (lblW != LBL_FREE_OF_TREE && !treeLive(W)) lblW = LBL_FREE_OF_TREE;
  // An edge into a live ODD node has no grow/formation event in this forest.
  // Reject it before the substantially more expensive blossom-aware slack.
#if !defined(EVENT_S5_TARGET_KEY_SHADOW) && \
    (!defined(EVENT_SFM_SHADOW) || \
    (defined(EVENT_SFM_BOUNDARY_FRONTIER) && defined(EVENT_SFM_REUSE) && \
     defined(EVENT_SFM_LAZY_ACTIVE)))
  if (lblW == LBL_ODD) {
#if defined(EVENT_SFM_SHADOW) && defined(EVENT_SFM_BOUNDARY_FRONTIER) && \
    defined(EVENT_SFM_REUSE) && defined(EVENT_SFM_LAZY_ACTIVE)
    if (mEventStructuresActive && W >= mN) {
      // Only a real INNER blossom can later produce a split-findmin
      // dissolution query. An ODD singleton has no represented boundary to
      // split, so feeding its incident edges is pure search overhead.
    } else
#endif
    {
#ifdef EVENT_OFFER_AUDIT
      ++mCounters.eventEdgeOddRejects;
#endif
      return;
    }
  }
#endif
  // During a fresh-root initialization both ends of every EVEN--EVEN edge are
  // scanned in the same complete pass. Keep one canonical orientation. This
  // is deliberately not used for later grow/shrink offers: an endpoint that
  // has just become EVEN may need an event that did not exist in the old role.
  if (dedupFreshEven && lblW == LBL_EVEN && u > w) {
#ifdef EVENT_OFFER_AUDIT
    ++mCounters.eventEdgeFreshDedupRejects;
#endif
    return;
  }
  // The caller is already walking this exact adjacency slot. Graph::Set keeps
  // one edge per endpoint pair, so using its weight is identical to w2max but
  // avoids rescanning u's complete row for every scheduled edge.
  // W and E are distinct outer nodes (the same-node case returned above), so
  // no current blossom can enclose both endpoints.  Computing slack through
  // slackZW would repeat the two union-find lookups merely to prove that fact
  // again for every scheduled edge.  Inherited blossoms are a separate,
  // immutable forest and must still contribute through oldCommonZ.
  const LL yu = sourceY == std::numeric_limits<LL>::max()
                  ? yOfOuter(u, E) : sourceY;
  LL s = yu + yOfOuter(w, W) + oldCommonZ(u, w) - edgeW2;
#ifdef EVENT_SFM_SHADOW
  // Section 5.1 shadow feed: while w belongs to a live INNER root, retain the
  // minimum slack of an incident edge from the OUTER forest at w's leaf. An
  // OUTER--INNER edge's slack is invariant under the current dual clock, so
  // this key remains meaningful until a structural role transition. The
  // production queue remains authoritative until dissolution queries match
  // the exact boundary scan on every gate.
  if (lblW == LBL_ODD
#ifdef EVENT_SFM_PRESEED_FREE
      || lblW == LBL_FREE_OF_TREE
#endif
#ifdef EVENT_SFM_BOUNDARY_FRONTIER
      || (lblW == LBL_FREE_OF_TREE && mEventSfmFeedFreeBoundary)
#endif
      ) {
    if (mEventStructuresActive && w < mN && (Val)s < mInnerSfm.key(w)) {
      mInnerSfm.decreasekey(w, (Val)s);
      mSfmKeySource[w] = u;
      mSfmKeyOwner[w] = E;
    }
    if (lblW == LBL_ODD) {
#ifdef EVENT_OFFER_AUDIT
      ++mCounters.eventEdgeOddRejects;
#endif
      return;
    }
  }
#endif
#ifdef EVENT_S5_TARGET_KEY_SHADOW
  bool s5TargetKeyDecreased = false;
  Size s5SourceToken = scNull();
  if (mS5TargetKeyActive && u < mN && E < mNRoot.size() &&
      mNRoot[E] != scNull() && mNLabel[E] == LBL_EVEN && treeLive(E)) {
    Size ro = outerOf(mNRoot[E]);
    if (ro < mBBase.size() && mBBase[ro] < mN) {
      s5SourceToken = mBBase[ro];
      if (mS5TargetKeyTrackedVertexStamp[u] != mS5TargetKeyTouchClock ||
          mS5TargetKeyTrackedVertexToken[u] != s5SourceToken) {
        mS5TargetKeyTrackedVertexToken[u] = s5SourceToken;
        mS5TargetKeyTrackedVertexStamp[u] = mS5TargetKeyTouchClock;
        s5TargetKeyRememberToken(s5SourceToken);
        mS5TargetKeyTreeVertices[s5SourceToken].push_back(u);
      }
    }
  }
  if (mS5TargetKeyActive && mNLabel[E] == LBL_EVEN && treeLive(E) &&
      lblW == LBL_FREE_OF_TREE && w < mN) {
    // An OUTER--FREE edge's slack decreases once per dual unit. Store its
    // criterion-specific absolute eligibility time, not its raw offer slack.
    // This is the clock-normalized Section 5 slack* key that survives forest
    // transitions and SearchOne boundaries.
    LL gap = 0;
    if (s5TargetKeyGap(s, mMate[u] == w, mS5TargetKeyCriterion, &gap)) {
      LL absoluteNow = 0, fire = 0;
      if (!checkedAddLL(mS5TargetKeyClockBase, tNow, &absoluteNow) ||
          !checkedAddLL(absoluteNow, gap, &fire)) {
        mLastStatus = STATUS_UNSUPPORTED_OR_OVERFLOWING_INPUT;
        return;
      }
      ++mCounters.s5TargetKeyOffers;
      if (mS5TargetKeySeeding) ++mCounters.s5TargetKeySeedOffers;
      const Val oldKey = s5TargetSfmKey(w);
      if ((Val)fire < oldKey) {
        if (oldKey == SplitFindmin::INF) s5TargetKeyRememberTarget(w);
        s5TargetSfmDecreaseKey(w, (Val)fire);
        mS5TargetKeySource[w] = u;
        s5TargetKeyReindexPayload(w, E, W, s5SourceToken,
                                  oldKey != SplitFindmin::INF, true);
        mS5TargetKeyEdgeIndex[w] = edgeIndex;
        mS5TargetKeyMatched[w] = mMate[u] == w;
        ++mCounters.s5TargetKeyDecreases;
        if (mS5TargetKeySeeding) ++mCounters.s5TargetKeySeedDecreases;
        s5TargetKeyDecreased = true;
      }
    }
  }
#endif
  // Criteria 1--3 use different slack* targets.  In particular Criterion 3
  // accepts BOTH slack 0 and -2, so an edge already at either value fires at
  // tNow; scheduling every Criterion-3 edge toward -2 incorrectly delayed a
  // slack-0 augmentation and changed the scaling dual trajectory.  Criterion 2
  // similarly targets 0 on matched edges and -2 on unmatched edges.
  LL target = 0;
  bool eligibleNow = (crit == CRIT1) ? (s == 0) :
                     (crit == CRIT3) ? (s == 0 || s == -2) :
                     ((mMate[u] == w) ? (s == 0) : (s == -2));
  if (!eligibleNow) {
    if (crit == CRIT2 && mMate[u] != w) target = -2;
    else if (crit == CRIT3 && s < 0) target = -2;
  }
  LL gap = s - target;
  LL d = -1; int kind = -1;
  // For internal shell edges the C*/D* translation has net zero effect on yz:
  // its endpoint y increases are cancelled by the enclosing inherited-z
  // decreases. Their event rates are therefore the ordinary Edmonds rates.
  if (eligibleNow) {
    d = 0;
    kind = (lblW == LBL_EVEN) ? 1 :
           (lblW == LBL_FREE_OF_TREE ? 0 : -1);
  } else if (gap >= 0 && lblW == LBL_EVEN) {
    // EVEN--EVEN slack changes by two per unit.  An odd gap never reaches the
    // requested exact slack value, so it has no edge event in this forest.
    if ((gap & 1) == 0) d = gap / 2;
    kind = 1;
  } else if (gap >= 0 && lblW == LBL_FREE_OF_TREE) {
    d = gap;
    kind = 0;
  }
  if (d < 0 || kind < 0) {
#ifdef EVENT_S5_TARGET_KEY_DRIVER_AUDIT
    if (s5TargetKeyDecreased) scheduleS5TargetMinEvent(w);
#endif
#ifdef EVENT_OFFER_AUDIT
    ++mCounters.eventEdgeNoEventRejects;
#endif
    return;
  }
#ifdef EVENT_OFFER_AUDIT
  ++mCounters.eventEdgeCandidates;
  if (d == 0) ++mCounters.eventEdgeImmediateCandidates;
#endif
  LL fire = 0;
  if (!checkedAddLL(tNow, d, &fire)) {
    mLastStatus = STATUS_UNSUPPORTED_OR_OVERFLOWING_INPUT; return;
  }
  // Indexed edge events retain the source adjacency slot in the same compact
  // payload. Pop validation then recovers both endpoint and exact weight in
  // O(1), and the slot is a stable identity for later event ownership/minimum
  // maintenance. Kinds 5/6 are indexed forms of edge kinds 0/1.
  // An indexed caller is iterating this exact adjacency row and passes its
  // current slot. Re-reading the vector bounds and target for every candidate
  // repeated work tens of millions of times. Unindexed structural-repair
  // callers retain the endpoint payload form.
  bool indexed = edgeIndex != scNull();
  int queueKind = indexed ? kind + 5 : kind;
  long long queueB = indexed ? (long long)edgeIndex : (long long)w;
#if defined(EVENT_S5_TARGET_KEY_DRIVER_NO_ORACLE)
  // OUTER--FREE scheduling is supplied by the criterion-aware retained list
  // minimum below. EVEN--EVEN formation, dissolution, and shell-control events
  // remain in their ordinary Section-5 scheduler roles.
#ifdef EVENT_S5_FULL_KEY_DRIVER_NO_ORACLE
  // Fixed-dual Figure-7 shell MCM deliberately bypasses the global target-key
  // layer, so its OUTER--FREE event stays in this local exact queue.
  if (kind != 0 || mShellActive)
#else
  if (!(crit == CRIT2 && kind == 0))
#endif
#endif
  schedulePracticalEvent(fire, (long long)u, queueB, queueKind,
                         (long long)E);
#ifdef EVENT_S5_TARGET_KEY_DRIVER_AUDIT
  // Queue the retained minimum after its ordinary oracle payload.  Both queues
  // are FIFO for equal timestamps, so an unchanged ordinary event remains the
  // authoritative structural transition; reveal-only minima can still prove
  // that the target representation is sufficient to drive a missing event.
  if (s5TargetKeyDecreased && mLastStatus == STATUS_SUCCESS)
    scheduleS5TargetMinEvent(w);
#endif
}

inline void ScalingMatching::scheduleShellBoundaryEvent(LL tNow) {
  if (!mShellActive || mShellCstar == scNull()) return;
  if (mShellCstar >= mZ.size() || mZ[mShellCstar] < 0 ||
      (mZ[mShellCstar] & 1) != 0 ||
      (mShellDstar != scNull() &&
       (mShellDstar >= mZ.size() || mZ[mShellDstar] < 0 ||
        (mZ[mShellDstar] & 1) != 0))) {
    mLastStatus = STATUS_INVARIANT_VIOLATION;
    return;
  }
  LL units = mZ[mShellCstar] / 2;
  if (mShellDstar != scNull() && mZ[mShellDstar] / 2 < units)
    units = mZ[mShellDstar] / 2;
  if (units < 0) {
    mLastStatus = STATUS_INVARIANT_VIOLATION;
    return;
  }
  LL fire = 0;
  if (!checkedAddLL(tNow, units, &fire)) {
    mLastStatus = STATUS_UNSUPPORTED_OR_OVERFLOWING_INPUT;
    return;
  }
  schedulePracticalEvent(fire, (long long)mShellCstar,
                         (long long)mShellDstar, 3 /*old boundary*/);
}

// Advance materialized duals from the current time to absolute time newTime,
// applying delta = newTime - *tNow over all labeled outer nodes' vertices/z.
inline void ScalingMatching::applyDualTo(LL newTime, LL* tNow) {
  if (newTime < *tNow) {
    mLastStatus = STATUS_INVARIANT_VIOLATION;
    return;
  }
  LL delta = newTime - *tNow;
  if (delta <= 0) { *tNow = newTime; return; }
  if (mShellActive) {
    ScalingNanosecondGuard shellDualTimer(
      &mCounters.shellDualAdjustmentNanoseconds);
    // Materialized Figure 7 reference: ordinary shell-confined Edmonds dual
    // adjustment followed by one C*/D* translation per unit. Lazy timestamps
    // are intentionally disabled for shell searches until this path is exact.
    if (mShellCstar == scNull() || mShellCstar >= mZ.size() ||
        (mShellDstar != scNull() && mShellDstar >= mZ.size())) {
      mLastStatus = STATUS_INVARIANT_VIOLATION;
      return;
    }
    LL twiceDelta = 0;
    if (!checkedMulLL(2, delta, &twiceDelta)) {
      mLastStatus = STATUS_UNSUPPORTED_OR_OVERFLOWING_INPUT;
      return;
    }
    // Preflight every update before mutating canonical state. This makes an
    // overflowing translation fail atomically instead of leaving half-applied
    // y/z values behind.
    mShellNextY.resize(mShellMarked.size());
    mCounters.shellDualVertexVisits += mShellMarked.size();
    for (Size i = 0; i < mShellMarked.size(); ++i) {
      Size v = mShellMarked[i];
      if (v >= mN) { mShellNextY[i] = 0; continue; }
      Size V = outerOf(v); LL next = mY[v];
      if ((mNLabel[V] == LBL_EVEN && !checkedAddLL(next, -delta, &next)) ||
          (mNLabel[V] == LBL_ODD && !checkedAddLL(next, delta, &next)) ||
          (mInDstar[v] && !checkedAddLL(next, twiceDelta, &next)) ||
          (!mInDstar[v] && mInShell[v] &&
           !checkedAddLL(next, delta, &next))) {
        mLastStatus = STATUS_UNSUPPORTED_OR_OVERFLOWING_INPUT;
        return;
      }
      mShellNextY[i] = next;
    }
    // beginCompactEventLabels/rememberCompactEventLabel maintain the exact
    // set of nodes that have a role in this search.  Historical blossom ids
    // are monotone, so sweeping [mN,mNumNodes) twice here made a small shell
    // pay for every blossom ever allocated.  Build the checked z transaction
    // from the compact labeled-root frontier instead.
    mShellNextZNodes.clear();
    mShellNextZ.clear();
    mShellNextZNodes.reserve(mEventTouchedNodes.size());
    mShellNextZ.reserve(mEventTouchedNodes.size());
#ifdef SCALING_INVARIANTS
    std::vector<char> compactShellDualRoot(mNumNodes, 0);
#endif
    for (Size i = 0; i < mEventTouchedNodes.size(); ++i) {
      Size b = mEventTouchedNodes[i];
      if (b < mN || b >= mNumNodes || mBParent[b] != scNull() ||
          b == mShellCstar || b == mShellDstar) continue;
      int role = mNLabel[b];
      if (role != LBL_EVEN && role != LBL_ODD) continue;
      LL next = mZ[b];
      if (!checkedAddLL(next, role == LBL_EVEN ? twiceDelta : -twiceDelta,
                        &next)) {
        mLastStatus = STATUS_UNSUPPORTED_OR_OVERFLOWING_INPUT;
        return;
      }
      mShellNextZNodes.push_back(b);
      mShellNextZ.push_back(next);
#ifdef SCALING_INVARIANTS
      compactShellDualRoot[b] = 1;
#endif
    }
#ifdef SCALING_INVARIANTS
    // Reviewer-facing semantic guard: prove the compact frontier is exactly
    // the exhaustive root definition before committing any dual change.
    for (Size b = mN; b < mNumNodes; ++b) {
      bool exhaustive = mBParent[b] == scNull() &&
        b != mShellCstar && b != mShellDstar &&
        (mNLabel[b] == LBL_EVEN || mNLabel[b] == LBL_ODD);
      if (exhaustive != (compactShellDualRoot[b] != 0)) {
        mLastStatus = STATUS_INVARIANT_VIOLATION;
        return;
      }
    }
#endif
    LL nextC = 0, nextD = 0;
    if (!checkedAddLL(mZ[mShellCstar], -twiceDelta, &nextC) || nextC < 0 ||
        (mShellDstar != scNull() &&
         (!checkedAddLL(mZ[mShellDstar], -twiceDelta, &nextD) || nextD < 0))) {
      mLastStatus = STATUS_INVARIANT_VIOLATION;
      return;
    }
    for (Size i = 0; i < mShellNextZNodes.size(); ++i) {
      Size b = mShellNextZNodes[i];
      mZ[b] = mShellNextZ[i];
      if (mNLabel[b] == LBL_EVEN && b < mBEverPos.size()) mBEverPos[b] = 1;
    }
    mZ[mShellCstar] = nextC;
    if (mShellDstar != scNull()) mZ[mShellDstar] = nextD;
    mCounters.shellDualVertexVisits += mShellMarked.size();
    for (Size i = 0; i < mShellMarked.size(); ++i) {
      Size v = mShellMarked[i];
      if (v >= mN) continue;
      mY[v] = mShellNextY[i];
    }
    mTNow = newTime;
    *tNow = newTime;
    return;
  }
  if (mLazyDuals) {
    mTNow = newTime;
    *tNow = newTime;
    return;
  }
  if (mEventCompactInitialization) {
    // Compact initialization already maintains the exact set of nodes that
    // received a search role. Current root blossoms partition their leaves,
    // so applying the adjustment through those roots avoids an O(n) vertex
    // sweep and its outerOf() lookup at every SearchOne boundary.
    for (Size i = 0; i < mEventTouchedNodes.size(); ++i) {
      Size node = mEventTouchedNodes[i];
      if (node >= mNumNodes || mBParent[node] != scNull()) continue;
      int role = mNLabel[node];
      if (role != LBL_EVEN && role != LBL_ODD) continue;
      LL vertexDelta = role == LBL_EVEN ? -delta : delta;
#ifdef SCALING_LEAF_NOVEC
      Size leafEnd = node < mLfTailX.size() ? mLfTailX[node] : scNull();
      for (Size v = node < mLfHeadX.size() ? mLfHeadX[node] : scNull(),
                next = scNull();
           v != scNull(); v = next) {
        next = (v == leafEnd) ? scNull() : mLfNextX[v];
        mY[v] += vertexDelta;
      }
#else
      const std::vector<Size>& leaves = mNodeLeaves[node];
      for (Size li = 0; li < leaves.size(); ++li)
        mY[leaves[li]] += vertexDelta;
#endif
      if (node >= mN)
        mZ[node] += role == LBL_EVEN ? 2 * delta : -2 * delta;
    }
    mTNow = newTime;
    *tNow = newTime;
    return;
  }
  for (Size v = 0; v < mN; ++v) {
    Size V = outerOf(v);
    if (mNLabel[V] == LBL_EVEN) mY[v] -= delta;
    else if (mNLabel[V] == LBL_ODD) mY[v] += delta;
  }
  for (Size b = mN; b < mNumNodes; ++b)
    if (mBParent[b] == scNull()) {
      if (mNLabel[b] == LBL_EVEN) mZ[b] += 2 * delta;
      else if (mNLabel[b] == LBL_ODD) mZ[b] -= 2 * delta;
    }
  mTNow = newTime;
  *tNow = newTime;
}

// After a structural change that frees sub-nodes (dissolve), relabel the forest
// from the current matching and re-seed the event PQ from scratch. tNow is reset
// to 0 (a fresh time origin; materialized duals already reflect all past adjusts).
inline void ScalingMatching::rememberCompactEventLabel(Size node) {
  if (!mEventCompactInitialization || node >= mNumNodes) return;
  if (mEventTouchedStamp.size() <= node)
    mEventTouchedStamp.resize(node + 1, 0);
  if (mEventDormantFlag.size() <= node)
    mEventDormantFlag.resize(node + 1, 0);
  // An exposed/relinked node is live again. Its stale entry, if any, remains
  // harmless in mEventDormantNodes until the next scale-boundary compaction.
  mEventDormantFlag[node] = 0;
  if (mEventTouchedStamp[node] == mEventTouchedClock) return;
  mEventTouchedStamp[node] = mEventTouchedClock;
  mEventTouchedNodes.push_back(node);
}

inline void ScalingMatching::beginCompactEventLabels() {
  mCounters.eventInitializationNodeVisits += mEventTouchedNodes.size();
  for (Size j = 0; j < mEventTouchedNodes.size(); ++j) {
    Size node = mEventTouchedNodes[j];
    if (node >= mNLabel.size()) continue;
    // Match the reference initializer exactly: labels on contracted children
    // are dormant state owned by their live parent and are rewritten if a
    // later dissolution exposes them.  Clearing them here changes tied event
    // trajectories even though final optimum weights usually remain equal.
    if (node < mBParent.size() && mBParent[node] != scNull()) {
      ++mCounters.eventInitializationDormantDeferrals;
      if (mEventDormantFlag.size() <= node)
        mEventDormantFlag.resize(node + 1, 0);
      if (!mEventDormantFlag[node]) {
        mEventDormantFlag[node] = 1;
        mEventDormantNodes.push_back(node);
      }
      continue;
    }
    mNLabel[node] = LBL_FREE_OF_TREE;
    mNParent[node] = scNull();
    mNParEdgeU[node] = scNull(); mNParEdgeV[node] = scNull();
    mNRoot[node] = scNull();
    ++mCounters.eventInitializationNodesReset;
  }
  mEventTouchedNodes.clear();
  mEventInitialRoots.clear();
  if (++mEventTouchedClock == 0) {
    std::fill(mEventTouchedStamp.begin(), mEventTouchedStamp.end(), 0);
    ++mEventTouchedClock;
  }
}

inline bool ScalingMatching::edmondsSearchEvent(
    int crit, const std::vector<Size>* explicitRoots,
    LL maxDualUnits, bool exhaustDualBudget, bool batchSameTimeAugment,
    LL maxAdjustments) {
  // Figure 3 requires the paths selected by one SearchOne augmentation phase
  // to be vertex-disjoint.  Mark exactly the contracted path nodes traversed
  // by augmentFrom/flipUpToRoot; unlike the old whole-tree retirement rule,
  // unused branches of a consumed alternating tree remain available to other
  // live roots in the same phase.
  const bool disjointSearchOnePhase = batchSameTimeAugment && crit == CRIT2 &&
                                      !mShellActive;
  if (disjointSearchOnePhase) {
    if (mUsedV.size() < mN) mUsedV.assign(mN, 0);
    else std::fill(mUsedV.begin(), mUsedV.begin() + mN, (char)0);
  }
  bool reuseBoundary = mEventBoundaryReuse && mEventBoundaryReuseReady &&
    crit == CRIT2 && explicitRoots == NULL && !mShellActive &&
    maxAdjustments == 1 && batchSameTimeAugment &&
    !mEventBoundaryReuseEvents.empty();
#ifdef EVENT_S5_TARGET_KEY_SHADOW
  if (!reuseBoundary) mS5TargetKeyActive = false;
#endif
  bool previousCriterion2AdjustmentOnly =
    crit == CRIT2 && mPreviousCriterion2AdjustmentOnly;
  mPreviousCriterion2AdjustmentOnly = false;
  if (!reuseBoundary) ++mCounters.searchInitializations;
  if (crit == CRIT1) ++mCounters.eventCriterion1Searches;
  else if (crit == CRIT2) ++mCounters.eventCriterion2Searches;
  else if (crit == CRIT3) ++mCounters.eventCriterion3Searches;
#ifdef SCPHASE
  std::fprintf(stderr,
    "[SCPHASE]    search=%llu enter crit=%d nodes=%lu clock=%.3fs\n",
    mCounters.searchInitializations, crit, (unsigned long)mNumNodes,
    (double)clock()/CLOCKS_PER_SEC);
  std::fflush(stderr);
#endif
  ScalingNanosecondGuard searchTimer(&mCounters.eventSearchNanoseconds);
  unsigned long long forestInitStart = scalingMonotonicNanoseconds();
  mLastSearchDualUnits = 0;
  mEventOneAdjustment = (maxAdjustments == 1);
#ifdef EVENT_OWNER_MIN_AUDIT
  if (mEventOwnerMinEpoch.size() < mNumNodes) {
    mEventOwnerMinEpoch.resize(mNumNodes, 0);
    mEventOwnerMinTime.resize(mNumNodes, 0);
    mEventOwnerMinTies.resize(mNumNodes, 0);
    mEventOwnerMinVersion.resize(mNumNodes, 0);
    mEventOwnerMinPrimary.resize(mNumNodes);
  }
  if (++mEventOwnerMinClock == 0) {
    std::fill(mEventOwnerMinEpoch.begin(), mEventOwnerMinEpoch.end(), 0);
    ++mEventOwnerMinClock;
  }
  mEventOwnerMinCurrentPayloads = 0;
  mEventOwnerMinTouched.clear();
  mEventOwnerMinTieLog.clear();
  mEventOwnerMatchedFuture.clear();
  mEventOwnerControlFuture.clear();
#endif
  mEventDualBudget = maxDualUnits;
  mEventEpochDualBase = 0;
  Size rootEpoch = 0;
  if (explicitRoots && !mShellActive) {
    if (mSearchRootEpoch.size() < mN) mSearchRootEpoch.resize(mN, 0);
    if (++mSearchRootClock == 0) {
      std::fill(mSearchRootEpoch.begin(), mSearchRootEpoch.end(), 0);
      ++mSearchRootClock;
    }
    rootEpoch = mSearchRootClock;
    for (Size i = 0; i < explicitRoots->size(); ++i)
      if ((*explicitRoots)[i] < mN)
        mSearchRootEpoch[(*explicitRoots)[i]] = rootEpoch;
  }
  std::vector<BucketPQ::Event> reuseEvents;
  if (reuseBoundary) reuseEvents.swap(mEventBoundaryReuseEvents);
  else mEventBoundaryReuseEvents.clear();
  mEventBoundaryReuseReady = false;
  LL tNow = 0;
  if (!reuseBoundary) {
  unsigned long long initPhaseStart = scalingMonotonicNanoseconds();
  // Relabel once (fresh roots), seed events for all EVEN nodes + ODD dissolutions,
  // then process events in time order. Persistent across augments (treeLive).
  mEventFreeCandidates.clear();
  // A shell search is confined to mShellMarked. Scanning every graph vertex
  // here made each exact-MCM certification pass O(n), even for a three-vertex
  // atomic shell. The shell mask is already the authoritative region and is
  // deduplicated while the boundary is built, so seed roots directly from it.
  const bool useCriterion2FreeFrontier =
#ifdef EVENT_CRIT2_FREE_FRONTIER_DISABLE
    false;
#else
    crit == CRIT2 && !mShellActive && explicitRoots == NULL && mAuditScale != 0;
#endif
  if (mShellActive && explicitRoots != NULL) {
    // Exact shell-MCM phases only remove exposed roots: an augmentation
    // matches two roots and never creates an exposed vertex. Consume the
    // driver's monotonically shrinking frontier directly instead of scanning
    // every vertex in the shell and filtering it through an epoch stamp.
    mCounters.eventFreeCandidateVerticesExamined += explicitRoots->size();
    for (Size i = 0; i < explicitRoots->size(); ++i) {
      Size v = (*explicitRoots)[i];
      if (v >= mN || !mInShell[v] || mMate[v] != scNull()) continue;
      Size N = outerOf(v);
      if (N >= mN && mBBase[N] != v) {
        ++mCounters.algorithmFallbackAttempts;
        ++mCounters.malformedFreeRootAttempts;
        mLastStatus = STATUS_INVARIANT_VIOLATION;
        return false;
      }
      mEventFreeCandidates.push_back(v);
    }
  } else if (explicitRoots != NULL) {
    // The request already names the complete non-shell root set. The former
    // path stamped these vertices, scanned all mN vertices, and selected the
    // same set through the stamp—catastrophic for singleton stage-2 PQSearch.
    // Sort/deduplicate to retain the exact ascending order of that old scan.
    mCounters.eventFreeCandidateVerticesExamined += explicitRoots->size();
    for (Size i = 0; i < explicitRoots->size(); ++i) {
      Size v = (*explicitRoots)[i];
      if (v >= mN || mMate[v] != scNull()) continue;
      Size N = outerOf(v);
      if (N >= mN && mBBase[N] != v) {
        ++mCounters.algorithmFallbackAttempts;
        ++mCounters.malformedFreeRootAttempts;
        mLastStatus = STATUS_INVARIANT_VIOLATION;
        return false;
      }
      mEventFreeCandidates.push_back(v);
    }
    std::sort(mEventFreeCandidates.begin(), mEventFreeCandidates.end());
    mEventFreeCandidates.erase(
      std::unique(mEventFreeCandidates.begin(), mEventFreeCandidates.end()),
      mEventFreeCandidates.end());
  } else if (useCriterion2FreeFrontier) {
    if (mEventCriterion2FreeScale != mAuditScale ||
        mEventCriterion2FreeUniverse != mN) {
      mEventCriterion2Free.clear();
      for (Size v = 0; v < mN; ++v)
        if (mMate[v] == scNull()) mEventCriterion2Free.push_back(v);
      mEventCriterion2FreeScale = mAuditScale;
      mEventCriterion2FreeUniverse = mN;
      ++mCounters.eventFreeCandidateFullScans;
      mCounters.eventFreeCandidateVerticesExamined += mN;
    } else {
      ++mCounters.eventFreeCandidateCacheHits;
    }
    Size retained = 0;
    for (Size i = 0; i < mEventCriterion2Free.size(); ++i) {
      Size v = mEventCriterion2Free[i];
      ++mCounters.eventFreeCandidateVerticesExamined;
      if (v >= mN || mMate[v] != scNull()) continue;
      mEventCriterion2Free[retained++] = v;
      Size N = outerOf(v);
      if (N >= mN && mBBase[N] != v) {
        // A free vertex enclosed by a current blossom must be its base. Do not
        // liquidate/flatten a malformed state to make the search continue.
        ++mCounters.algorithmFallbackAttempts;
        ++mCounters.malformedFreeRootAttempts;
        mLastStatus = STATUS_INVARIANT_VIOLATION;
        return false;
      }
      mEventFreeCandidates.push_back(v);
    }
    mEventCriterion2Free.resize(retained);
  } else {
    Size freeScanCount = mShellActive ? mShellMarked.size() : mN;
    mCounters.eventFreeCandidateVerticesExamined += freeScanCount;
    for (Size freeScanIndex = 0; freeScanIndex < freeScanCount; ++freeScanIndex) {
      Size v = mShellActive ? mShellMarked[freeScanIndex] : freeScanIndex;
      if (v >= mN || (mShellActive && !mInShell[v])) continue;
      if (mMate[v] == scNull()) {
        if (rootEpoch != 0 && mSearchRootEpoch[v] != rootEpoch) continue;
        Size N = outerOf(v);
        if (N >= mN && mBBase[N] != v) {
          ++mCounters.algorithmFallbackAttempts;
          ++mCounters.malformedFreeRootAttempts;
          mLastStatus = STATUS_INVARIANT_VIOLATION;
          return false;
        }
        mEventFreeCandidates.push_back(v);
      }
    }
  }
  mCounters.eventInitializationFreeCandidates += mEventFreeCandidates.size();
  mCounters.eventFreeDiscoveryNanoseconds +=
    scalingMonotonicNanoseconds() - initPhaseStart;
#ifdef SCPHASE
  std::fprintf(stderr, "[SCPHASE]    search=%llu flatten-done clock=%.3fs\n",
               mCounters.searchInitializations, (double)clock()/CLOCKS_PER_SEC);
  std::fflush(stderr);
#endif
  initPhaseStart = scalingMonotonicNanoseconds();
  if (mEventCompactInitialization) {
    beginCompactEventLabels();
  } else {
    mEventInitialRoots.clear();
    for (Size i = 0; i < mNumNodes; ++i) {
      ++mCounters.eventInitializationNodeVisits;
      if (mBParent[i] != scNull()) continue;
      mNLabel[i] = LBL_FREE_OF_TREE;
      mNParent[i] = scNull(); mNParEdgeU[i] = scNull(); mNParEdgeV[i] = scNull();
      mNRoot[i] = scNull();
      ++mCounters.eventInitializationNodesReset;
    }
  }
  mCounters.eventLabelResetNanoseconds +=
    scalingMonotonicNanoseconds() - initPhaseStart;
  initPhaseStart = scalingMonotonicNanoseconds();
  bool anyRoot = false;
  for (Size freeIndex = 0; freeIndex < mEventFreeCandidates.size(); ++freeIndex) {
      Size v = mEventFreeCandidates[freeIndex];
      Size N = outerOf(v);
      if (mNLabel[N] == LBL_FREE_OF_TREE) {
        mNLabel[N]=LBL_EVEN; mNRoot[N]=N; anyRoot=true;
        rememberCompactEventLabel(N);
        mEventInitialRoots.push_back(N);
        ++mCounters.eventInitializationRoots;
      }
  }
  if (!anyRoot) return false;
  mCounters.eventRootInitializationNanoseconds +=
    scalingMonotonicNanoseconds() - initPhaseStart;
#ifdef SCPHASE
  std::fprintf(stderr, "[SCPHASE]    search=%llu roots-done clock=%.3fs\n",
               mCounters.searchInitializations, (double)clock()/CLOCKS_PER_SEC);
  std::fflush(stderr);
#endif

  // The primary split-findmin channel owns incremental INNER-dissolution
  // boundaries. Canonical Hybrid normally requests its scan-free frontier.
  // Independently, the shared secondary channel below always requests the same
  // current-Omega list topology and uses key2/findmin2 to drive ordinary
  // OUTER-to-non-outer grow events under Criteria 1--3. Keep the two requests
  // separate: disabling the primary boundary channel must not silently remove
  // the canonical grow representation.
  const bool integratedHybridEvent =
    mHybridCanonical && !mHybridReference && mEventSearch;
  bool eventSfmRequested = mIncrDissolve && !integratedHybridEvent;
#ifdef EVENT_SFM_SHADOW
  if (mIncrDissolve && integratedHybridEvent)
    eventSfmRequested = true;
#endif
  // The shared target-key channel cannot borrow a stale topology merely
  // because this particular search has no incremental-dissolution request.
  // It still needs one authoritative current-Omega list partition for its
  // key2/findmin2 operations.
  bool eventSfmTopologyRequested = eventSfmRequested;
#ifdef EVENT_S5_SHARED_SFM_TOPOLOGY
  // Figure-7 MCM runs at fixed duals. Its zero-time Criterion-1 queue is an
  // ordinary exact cardinality search and does not need the dual-moving
  // Section-5 target-key topology. Rebuilding/resetting that global layer for
  // thousands of tiny dynamic shells dominated real-graph execution while
  // contributing no event that the direct eligible-edge queue cannot express.
  if (!mShellActive) eventSfmTopologyRequested = true;
#endif
#if defined(EVENT_SFM_BOUNDARY_FRONTIER) && defined(EVENT_SFM_REUSE) && \
    defined(EVENT_SFM_LAZY_ACTIVE)
  mEventStructuresActive = false;
#else
  mEventStructuresActive = eventSfmRequested;
#endif
  // Figure 6 stage-2 FVR is bounded BucketSearch even under Criterion 3.
  // Unbounded Criterion-3 Finalization is PQSearch and therefore uses the
  // binary heap. DismantlePath's singleton Criterion-1 call is explicitly a
  // PQSearch in the paper, even though its execution is capped at T units.
  mEventUseBinary = !mShellActive && (crit == CRIT1 || maxDualUnits < 0);
  initPhaseStart = scalingMonotonicNanoseconds();
  if (eventSfmTopologyRequested
#ifdef EVENT_SFM_REUSE
      && (!mEventSfmInitialized || mEventSfmScale != mAuditScale)
#endif
      ) initializeEventStructures();
  mCounters.eventSfmInitializationNanoseconds +=
    scalingMonotonicNanoseconds() - initPhaseStart;
#ifdef EVENT_S5_TARGET_KEY_SHADOW
  // Start the criterion-aware secondary-key epoch. Persistent topology and
  // touched payloads are reused across eligible SearchOne boundaries.
  initPhaseStart = scalingMonotonicNanoseconds();
#ifdef EVENT_S5_FULL_KEY_DRIVER_AUDIT
  if (!mShellActive) initializeS5TargetKeyShadow(crit);
  else mS5TargetKeyActive = false;
#else
  if (crit == CRIT2) initializeS5TargetKeyShadow(crit);
#endif
  mCounters.eventTargetKeyInitializationNanoseconds +=
    scalingMonotonicNanoseconds() - initPhaseStart;
#endif

  mTNow = 0;
  // SearchOne stops after its first positive adjustment, so initializing and
  // later materializing a full timestamp epoch cannot amortize there. Keep the
  // Section 5.2 clock for multi-adjustment BucketSearch/PQSearch, where event
  // processing spans more than one dual boundary.
  initPhaseStart = scalingMonotonicNanoseconds();
  if (mEventLazyClock && !mShellActive && !mEventOneAdjustment) tsInit();
  mCounters.eventLazyInitializationNanoseconds +=
    scalingMonotonicNanoseconds() - initPhaseStart;
  mCounters.eventForestInitializationNanoseconds +=
    scalingMonotonicNanoseconds() - forestInitStart;
  resetPracticalEventQueue();
  if (mEventUseBinary) mUnboundedPQ.beginBatch();
  // Seed exactly the roots discovered above.  The former full node sweep could
  // only find these nodes after the fresh relabel and therefore did redundant
  // O(number-of-allocated-blossoms) work on every SearchOne.
  mEventDeferringSeedBoundary = mEventOneAdjustment &&
                                mEventBoundaryReuse && !mEventUseBinary;
#ifdef EVENT_S5_TARGET_KEY_SHADOW
  mS5TargetKeySeeding = true;
#endif
  for (Size i = 0; i < mEventInitialRoots.size(); ++i)
    scheduleNodeEdges(mEventInitialRoots[i], crit, tNow, true);
#ifdef EVENT_S5_TARGET_KEY_SHADOW
  mS5TargetKeySeeding = false;
#endif
  mEventDeferringSeedBoundary = false;
  if (mEventOneAdjustment && mEventBoundaryReuse && !mEventUseBinary &&
      mEventEarliestDeferredEdge != std::numeric_limits<LL>::max())
    schedulePracticalEvent(mEventEarliestDeferredEdge, -1, -1,
                           4 /*deferred edge boundary*/);
#ifdef SCPHASE
  std::fprintf(stderr,
    "[SCPHASE]    search=%llu seed-done queue=%lu scans=%llu clock=%.3fs\n",
    mCounters.searchInitializations,
    (unsigned long)(mEventUseBinary ? mUnboundedPQ.size() : mPQ.size()),
    mCounters.edgeScans, (double)clock()/CLOCKS_PER_SEC);
  std::fflush(stderr);
#endif
  // Initialize dissolution events only for INNER roots in the requested
  // region.  The former shell path swept [mN,mNumNodes) on every exact-MCM
  // pass.  Since blossom ids are historical and monotone across scales, that
  // made a tiny fixed shell repeatedly inspect tens of thousands of dead or
  // unrelated records (de2010: ~19k passes over ~67k historical blossoms in
  // one inherited phase).  Every live shell root contains a marked shell
  // vertex, so outerOf(marked vertex) is a complete region-local enumeration.
  const auto scheduleInitialDissolution = [&](Size b) -> bool {
    if (b < mN || b >= mNumNodes || mBParent[b] != scNull() ||
        mNLabel[b] != LBL_ODD || b == mShellCstar || b == mShellDstar)
      return true;
    // A blossom can already be INNER when a fresh search begins.  Previously
    // only blossoms that entered the INNER role later in this search called
    // sfmShadowBeginInner(), leaving these initial roots at INF and allowing a
    // dissolution to expose a FREE child without its required grow event.
#ifdef EVENT_SFM_SHADOW
    if (!(mOldForestSeparated && b >= mAuditOldBegin && b < mAuditOldEnd))
      sfmShadowBeginInner(b);
#endif
    LL fire = 0;
    if (!checkedAddLL(tNow, zOf(b) / 2, &fire)) {
      mLastStatus = STATUS_UNSUPPORTED_OR_OVERFLOWING_INPUT;
      return false;
    }
    schedulePracticalEvent(fire, (long long)b, -1, 2 /*dissolve*/);
    return mLastStatus == STATUS_SUCCESS;
  };
  if (mShellActive) {
    if (mShellDissolveStamp.size() < mNumNodes)
      mShellDissolveStamp.resize(mNumNodes, 0);
    if (++mShellDissolveClock == 0) {
      std::fill(mShellDissolveStamp.begin(), mShellDissolveStamp.end(), 0);
      ++mShellDissolveClock;
    }
    for (Size i = 0; i < mShellMarked.size(); ++i) {
      Size v = mShellMarked[i];
      if (v >= mN || !mInShell[v]) continue;
      Size b = outerOf(v);
      if (b < mN || b >= mNumNodes ||
          mShellDissolveStamp[b] == mShellDissolveClock) continue;
      mShellDissolveStamp[b] = mShellDissolveClock;
      if (!scheduleInitialDissolution(b)) {
        materializeLazyDuals();
        mEventStructuresActive = false;
        return false;
      }
    }
  } else {
    for (Size b = mN; b < mNumNodes; ++b)
      if (!scheduleInitialDissolution(b)) {
        materializeLazyDuals();
        mEventStructuresActive = false;
        return false;
      }
  }
  scheduleShellBoundaryEvent(tNow);
  if (mEventUseBinary) mUnboundedPQ.endBatch();
  if (mLastStatus != STATUS_SUCCESS) {
    materializeLazyDuals();
    mEventStructuresActive = false;
    return false;
  }
  } else {
    // The preceding Criterion-2 SearchOne changed no matching edge and stopped
    // exactly after a dual adjustment.  Its surviving forest is therefore
    // still valid at this boundary.  Replay only the tied edge events that the
    // one-adjustment queue deliberately deferred; if they do not augment at
    // time zero, select the request's ordinary fresh initialization below.
    mEventStructuresActive = false;
    mEventUseBinary = false;
    mTNow = 0;
    resetPracticalEventQueue();
    for (Size i = 0; i < reuseEvents.size(); ++i) {
      const BucketPQ::Event& saved = reuseEvents[i];
      schedulePracticalEvent(0, saved.a, saved.b, saved.kind, saved.owner);
      if (mLastStatus != STATUS_SUCCESS) break;
    }
    ++mCounters.eventBoundaryReuseHits;
    mCounters.eventForestInitializationNanoseconds +=
      scalingMonotonicNanoseconds() - forestInitStart;
    if (mLastStatus != STATUS_SUCCESS) return false;
  }

  bool augmentedAny = false;
  LL adjustmentEvents = 0;
  bool adjustmentBoundary = false;
  bool boundaryReuseRequiresFreshInit = false;
  LL guard = 0;
  LL guardCap = 8 * (LL)mN * (LL)mN + 8 * (LL)mN * (mMaxW + 1) + 1000;
  BucketPQ::Event ev;
  while (popPracticalEvent(&ev)) {
#ifdef SCALING_EVENT_INVARIANTS_EVERY_EVENT
    // The previous transition is now complete and no next event has changed
    // state yet. Validate the paper's dual/blossom contract at this exact
    // between-events boundary.
    if (!validateEventInvariantState(crit)) break;
#endif
    if (++guard > guardCap) {
      mLastStatus = STATUS_UNEXPECTED_NO_PROGRESS;
      break;
    }
#ifdef SCPHASE
    if ((guard & (((LL)1 << 20) - 1)) == 0) {
      std::fprintf(stderr,
        "[SCPHASE]    event-search pops=%lld time=%lld kind=%d queue=%lu clock=%.3fs\n",
        (long long)guard, (long long)ev.time, ev.kind,
        (unsigned long)(mEventUseBinary ? mUnboundedPQ.size() : mPQ.size()),
        (double)clock()/CLOCKS_PER_SEC);
      std::fflush(stderr);
    }
#endif
    LL T = ev.time;
    if (reuseBoundary && T > 0) {
      boundaryReuseRequiresFreshInit = true;
      break;
    }
    // A bounded SearchOne returns immediately after its adjustment and the
    // resulting zero-z inner dissolutions. Edge events that became eligible
    // exactly at that boundary belong to the next SearchOne call.
    if (adjustmentBoundary && T > tNow) break;
    // Experimental phase batching stops before the first later dual level.
    // Thus every additional augmentation uses the unchanged duals and only
    // surviving alternating trees from the phase initialized above.
    if (batchSameTimeAugment && augmentedAny && T > tNow) break;
    // Advance only to the caller's exact budget boundary. Event timestamps are
    // absolute within the current reseed epoch; mLastSearchDualUnits persists
    // across event-clock epochs.
    LL delta = T - tNow;
    if (delta < 0) {
#ifdef EVENTTRACE
      std::fprintf(stderr, "[EVENT-INVARIANT] past event T=%lld now=%lld kind=%d\n",
                   (long long)T, (long long)tNow, ev.kind);
#endif
      mLastStatus = STATUS_INVARIANT_VIOLATION; break;
    }
    if (delta > 0 && maxAdjustments >= 0 &&
        adjustmentEvents >= maxAdjustments) break;
    if (maxDualUnits >= 0) {
      LL remaining = maxDualUnits - mLastSearchDualUnits;
      // A closed budget still admits events already eligible at the current
      // time. This is required both by Figure 6's boundary-closing SearchOne
      // and by Figure 7's zero-adjustment eligible-shell MCM phase.
      if (remaining < 0 || (remaining == 0 && delta > 0)) break;
      if (delta > remaining) {
        LL boundary = 0;
        if (!checkedAddLL(tNow, remaining, &boundary)) {
          mLastStatus = STATUS_UNSUPPORTED_OR_OVERFLOWING_INPUT;
          break;
        }
        applyDualTo(boundary, &tNow);
        if (mLastStatus != STATUS_SUCCESS) break;
        mLastSearchDualUnits += remaining;
        mCounters.dualUnits += (unsigned long long)remaining;
        if (remaining > 0) ++adjustmentEvents;
        break;
      }
    }
    // Most queue events are transitions already eligible at the current
    // clock. Avoid routing those millions of zero-time pops through the full
    // dual-adjustment entry point; applyDualTo has no side effect for delta 0.
    if (delta > 0) {
      applyDualTo(T, &tNow);
      if (mLastStatus != STATUS_SUCCESS) break;
    }
    mLastSearchDualUnits += delta;
    mCounters.dualUnits += (unsigned long long)delta;
    if (delta > 0) {
      ++adjustmentEvents;
      if ((maxAdjustments >= 0 && adjustmentEvents >= maxAdjustments) ||
          (maxDualUnits >= 0 && mLastSearchDualUnits >= maxDualUnits))
        adjustmentBoundary = true;
    }
    if (mShellActive && mShellCstar != scNull()) {
      if (mZ[mShellCstar] < 0 ||
          (mShellDstar != scNull() && mZ[mShellDstar] < 0)) {
        mLastStatus = STATUS_INVARIANT_VIOLATION;
        break;
      }
      // Figure 7 dissolves zero-z current roots before removing/updating an
      // inherited boundary. The wrapper performs the old-boundary update.
      if (mZ[mShellCstar] == 0 ||
          (mShellDstar != scNull() && mZ[mShellDstar] == 0)) {
        dissolveZeroCurrentRoots();
        break;
      }
    }
#ifdef SCDBG
    if (guard % 50 == 0 || guard > guardCap - 5) {
      Size cc=0; for(Size v=0;v<mN;++v) if(mMate[v]!=scNull())++cc;
      size_t queueSize = mEventUseBinary ? mUnboundedPQ.size() : mPQ.size();
      std::fprintf(stderr,"  [EVENT g=%ld] T=%ld tNow=%ld kind=%d a=%lld b=%lld matched=%lu nodes=%lu pq=%lu\n",
        (long)guard,(long)T,(long)tNow,ev.kind,ev.a,ev.b,(unsigned long)cc,(unsigned long)mNumNodes,(unsigned long)queueSize);
      std::fflush(stderr);
    }
#endif

#ifdef EVENT_S5_TARGET_KEY_DRIVER_AUDIT
    if (ev.kind == 10) {
      ++mCounters.s5TargetMinEventsPopped;
      Size source = ev.a >= 0 ? (Size)ev.a : scNull();
      Size target = ev.b >= 0 ? (Size)ev.b : scNull();
      unsigned int frontierVersion = ev.owner >= 0
        ? (unsigned int)ev.owner : 0;
      LL absoluteEvent = 0;
      bool eventClockOk = checkedAddLL(mS5TargetKeyClockBase, ev.time,
                                       &absoluteEvent);
      Val currentListKey = SplitFindmin::INF;
      Size currentListAt = scNull();
      Size currentListId = scNull();
      if (target < mN)
        s5TargetSfmFindMinWithList(target, &currentListKey, &currentListAt,
                                   &currentListId);
      bool currentFrontier = frontierVersion == 0;
      if (!currentFrontier && currentListId < mS5TargetFrontier.size()) {
        const S5TargetFrontierEntry& frontier =
          mS5TargetFrontier[currentListId];
        currentFrontier = frontier.version == frontierVersion &&
          frontier.winner == target && frontier.source == source;
      }
      // The logical frontier is now consumed. Any stale-source recovery or
      // role transition below may install its successor generation.
      if (currentFrontier && frontierVersion != 0)
        mS5TargetFrontier[currentListId].version = 0;
      bool currentPayload = target < mN && source < mN && eventClockOk &&
        mS5TargetKeyActive && currentFrontier &&
        currentListAt == target && currentListKey == (Val)absoluteEvent &&
        s5TargetSfmKey(target) == (Val)absoluteEvent &&
        mS5TargetKeySource[target] == source;
      Size currentTargetOwner = target < mN ? outerOfUF(target) : scNull();
      int currentTargetLabel = currentTargetOwner < mNLabel.size()
        ? mNLabel[currentTargetOwner] : LBL_ODD;
      bool consumedBatchTarget = currentTargetOwner < mNumNodes &&
        batchSameTimeAugment && augmentedAny &&
        currentTargetLabel != LBL_FREE_OF_TREE &&
        !treeLive(currentTargetOwner) &&
        (!disjointSearchOnePhase ||
         (target < mUsedV.size() && mUsedV[target]));
      if (currentTargetOwner < mNumNodes &&
          currentTargetLabel != LBL_FREE_OF_TREE &&
          !treeLive(currentTargetOwner))
        currentTargetLabel = LBL_FREE_OF_TREE;
      bool exact = crit == mS5TargetKeyCriterion &&
        target < mN && source < mN &&
        currentPayload && !consumedBatchTarget &&
        mS5TargetKeyTargetOwner[target] == currentTargetOwner &&
        currentTargetLabel == LBL_FREE_OF_TREE &&
        mS5TargetKeyOwner[target] == outerOfUF(source) &&
        mNLabel[outerOfUF(source)] == LBL_EVEN &&
        treeLive(outerOfUF(source)) && outerOfUF(source) != outerOfUF(target) &&
        (bool)mS5TargetKeyMatched[target] == (mMate[source] == target);
      if (exact && mS5TargetKeySourceToken[target] < mN) {
        Size sourceOwner = outerOfUF(source);
        Size root = sourceOwner < mNRoot.size() ? mNRoot[sourceOwner] : scNull();
        Size rootOuter = root < mNumNodes ? outerOf(root) : scNull();
        exact = rootOuter < mBBase.size() &&
                mBBase[rootOuter] == mS5TargetKeySourceToken[target];
      }
      if (!exact) {
        ++mCounters.s5TargetMinEventsStale;
#ifdef EVENT_S5_TARGET_KEY_DRIVER_NO_ORACLE
        // This was the currently retained minimum, not an older decrease-key
        // payload, but its source ceased to be a valid live OUTER owner. The
        // ordinary FREE-edge queue is absent in this mode, so reveal the next
        // live incident source at this exact structural transition.
        if (currentPayload && !consumedBatchTarget) {
          if (currentTargetLabel == LBL_FREE_OF_TREE)
            s5TargetKeyRevealFree(target, tNow);
        }
#endif
        continue;
      }
      ++mCounters.s5TargetMinEventsExact;
      ++mCounters.s5TargetMinEventsDriven;
      // Reuse the ordinary raw FREE-edge validation and transition below.  The
      // retained target payload has just supplied its own cancellable owner.
      ev.kind = 0;
      ev.owner = (long long)outerOfUF(source);
    }
#endif

    if (adjustmentBoundary && ev.kind != 2 &&
        !(ev.kind == 8 && ev.owner == 1)) continue;
    if (ev.kind == 3) {
      // A non-stale boundary marker must have been handled immediately after
      // applyDualTo above. If not, its boundary changed and the marker is stale.
      ++mCounters.eventControlPops;
      continue;
    }
    if (ev.kind == 2) {
      ++mCounters.eventControlPops;
      // dissolution: ODD root blossom b should have z==0 now.
      Size b = (Size)ev.a;
      if (b < mN || mBParent[b] != scNull() || mNLabel[b] != LBL_ODD) continue; // stale
      if (batchSameTimeAugment && augmentedAny && !treeLive(b)) continue;
      if (zOf(b) != 0) { // not yet zero (rescheduled by a later event) -> reschedule
        LL fire = 0; if (!checkedAddLL(tNow, zOf(b) / 2, &fire)) {
          mLastStatus = STATUS_UNSUPPORTED_OR_OVERFLOWING_INPUT; break; }
        schedulePracticalEvent(fire, (long long)b, -1, 2); continue;
      }
      bool dynamicDissolutionReconstruction = false;
#ifdef EVENT_SFM_SHADOW
      const bool sfmMaintainedAtDissolve = mEventStructuresActive;
#ifdef EVENT_SFM_CONSUME
      // Split-findmin has no merge operation.  A current blossom formed after
      // initialization is therefore topology-neutral in eventSfmOnShrink(),
      // and eventSfmOnDissolve() can represent its exposed children only when
      // the stable split pivots for this exact cycle exist.  Otherwise the
      // dissolution boundary must be reconstructed from adjacency; consuming
      // a singleton-list INF here silently loses required grow events.
      const bool sfmHasDissolutionTopology =
        b < mSfmSplitPivots.size() &&
        mSfmSplitPivots[b].size() + 1 == mBCycle[b].size();
      const bool sfmConsumesDissolution = sfmMaintainedAtDissolve &&
        sfmHasDissolutionTopology &&
        !(mOldForestSeparated && b >= mAuditOldBegin && b < mAuditOldEnd);
      dynamicDissolutionReconstruction = sfmMaintainedAtDissolve &&
        !sfmHasDissolutionTopology &&
        !(mOldForestSeparated && b >= mAuditOldBegin && b < mAuditOldEnd);
      if (dynamicDissolutionReconstruction)
        ++mCounters.dynamicDissolutionReconstructions;
#endif
#endif
      // Expand one laminar level in place, including nested members, and
      // schedule only newly exposed nodes. Failure means the maintained tree
      // is not a valid Section-5 dissolution state; never switch to relabel +
      // complete queue reseeding.
      if (mIncrDissolve && relinkOddCycle(b)) {
        std::vector<Size> children = mBCycle[b];
        // Close B's INNER interval and start each exposed child's new role at
        // the current clock. This is the Section 5.2 dissolution transition;
        // no O(n)-wide y/z materialization or timestamp-epoch restart.
        tsOnDissolveInner(b);
        if (mLastStatus != STATUS_SUCCESS) break;
        dissolve(b);
        // B is no longer a live outer node. Leaving its old ODD label behind
        // makes later materialized adjustments decrement the dissolved z below
        // zero (the n=6 Hybrid incremental-dissolution z=-4 failure).
        mNLabel[b] = LBL_FREE_OF_TREE;
        mNParent[b] = scNull(); mNParEdgeU[b] = scNull();
        mNParEdgeV[b] = scNull(); mNRoot[b] = scNull();
#ifdef EVENT_SFM_BOUNDARY_FRONTIER
        // Splitting an INNER blossom exposes every child simultaneously. Feed
        // all newly EVEN sibling rows before asking split-findmin for any FREE
        // child's boundary minimum; otherwise cycle order can make a query
        // miss a smaller candidate owned by a later sibling.
        mEventSfmFeedFreeBoundary = true;
        for (Size i = 0; i < children.size(); ++i) {
          Size child = children[i];
          if (mBParent[child] == scNull() && mNLabel[child] == LBL_EVEN) {
            scheduleNodeEdges(child, crit, tNow);
            if (mLastStatus != STATUS_SUCCESS) break;
          }
        }
        if (mLastStatus != STATUS_SUCCESS) {
          mEventSfmFeedFreeBoundary = false;
          break;
        }
#endif
        for (Size i = 0; i < children.size(); ++i) {
          Size child = children[i];
          if (mBParent[child] != scNull()) continue;
          if (mNLabel[child] == LBL_EVEN) {
#ifndef EVENT_SFM_BOUNDARY_FRONTIER
            scheduleNodeEdges(child, crit, tNow);
#endif
          }
          else if (mNLabel[child] == LBL_FREE_OF_TREE) {
            // Before expansion these edges ended in ODD container B, so no
            // event existed. Offer the newly exposed FREE boundary from the
            // reverse orientation of every unaffected live EVEN source.
#ifdef SCALING_LEAF_NOVEC
            std::vector<Size> leafScratch; scalingLfCollect(child, &leafScratch);
            const std::vector<Size>& leaves = leafScratch;
#else
            const std::vector<Size>& leaves = mNodeLeaves[child];
#endif
            unsigned long long dynamicReconstructionStart = 0;
            if (dynamicDissolutionReconstruction) {
              ++mCounters.dynamicDissolutionReconstructionChildren;
              mCounters.dynamicDissolutionReconstructionVertices += leaves.size();
              dynamicReconstructionStart = scalingMonotonicNanoseconds();
            }
#ifdef EVENT_SFM_SHADOW
            LL shadowScanMin = std::numeric_limits<LL>::max();
#endif
#if defined(EVENT_SFM_CONSUME) && defined(EVENT_SFM_NO_ORACLE)
            if (!sfmConsumesDissolution) {
#endif
            for (Size li = 0; li < leaves.size(); ++li) {
              Size v = leaves[li];
              for (Size ai = 0; ai < mAdj[v].size(); ++ai) {
                ++mCounters.edgeScans;
                ++mCounters.dissolutionBoundaryScans;
                if (dynamicDissolutionReconstruction)
                  ++mCounters.dynamicDissolutionReconstructionScans;
                Size u = mAdj[v][ai];
                if (u >= mN || outerOf(u) == child) continue;
                Size U = outerOf(u);
                if (mNLabel[U] != LBL_EVEN || !treeLive(U)) continue;
                LL edgeW2 = mUseScaleW ? mW[v][ai] : 2 * mWhat[v][ai];
#ifdef EVENT_SFM_SHADOW
                LL boundarySlack = yOfOuter(u, U) + yOfOuter(v, child) +
                                   oldCommonZ(u, v) - edgeW2;
                if (boundarySlack < shadowScanMin) {
                  shadowScanMin = boundarySlack;
                }
#endif
#ifdef EVENT_SFM_CONSUME
                if (!sfmConsumesDissolution)
#endif
                  scheduleEdgeFromEven(u, v, U, crit, tNow, edgeW2);
                if (mLastStatus != STATUS_SUCCESS) break;
              }
              if (mLastStatus != STATUS_SUCCESS) break;
            }
#if defined(EVENT_SFM_CONSUME) && defined(EVENT_SFM_NO_ORACLE)
            }
#endif
            if (dynamicDissolutionReconstruction)
              mCounters.dynamicDissolutionReconstructionNanoseconds +=
                scalingMonotonicNanoseconds() - dynamicReconstructionStart;
#ifdef EVENT_SFM_SHADOW
            if (mLastStatus == STATUS_SUCCESS && !leaves.empty() &&
                sfmMaintainedAtDissolve &&
                !(mOldForestSeparated && b >= mAuditOldBegin &&
                  b < mAuditOldEnd)) {
#ifdef EVENT_SFM_TOPOLOGY_AUDIT
              if (!sfmShadowListMatches(child)) {
                ++mCounters.sfmSplitTopologyFailures;
              }
#endif
              Val shadowKey = SplitFindmin::INF;
              Size shadowVertex = scNull();
              mInnerSfm.findmin(leaves[0], &shadowKey, &shadowVertex);
              ++mCounters.sfmBoundaryQueries;
#ifdef EVENT_SFM_LAZY_REPAIR
              // A retained vertex minimum can become stale when its source is
              // absorbed, relabelled, or moved to the same outer component.
              // Canonical Hybrid lazily deletes an invalid retained minimum,
              // just as it discards a stale event-queue entry. It never scans
              // the target row to reconstruct a successor: live successors
              // are maintained by the full Section-5 key channel.
              Size repairGuard = 0;
              bool repairValidated = shadowKey == SplitFindmin::INF;
              while (shadowKey != SplitFindmin::INF &&
                     repairGuard++ <= leaves.size()) {
                bool shadowInChild = shadowVertex < mN &&
                  outerOf(shadowVertex) == child;
                bool keyValid = false;
                if (shadowInChild && shadowVertex < mSfmKeySource.size()) {
                  Size source = mSfmKeySource[shadowVertex];
                  Size owner = mSfmKeyOwner[shadowVertex];
                  if (source < mN && owner < mNumNodes) {
                    Size currentSource = outerOf(source);
                    Size currentTarget = outerOf(shadowVertex);
                    if (currentSource != currentTarget &&
                        currentSource == outerOf(owner) &&
                        mNLabel[currentSource] == LBL_EVEN &&
                        treeLive(currentSource)) {
                      LL currentSlack = yOfOuter(source, currentSource) +
                        yOfOuter(shadowVertex, currentTarget) +
                        oldCommonZ(source, shadowVertex) -
                        w2max(source, shadowVertex);
                      keyValid = (Val)currentSlack == shadowKey;
                    }
                  }
                }
                if (keyValid) {
                  repairValidated = true;
                  break;
                }

                if (mHybridCanonical) {
#ifdef EVENTTRACE
                  Size traceSource = shadowVertex < mSfmKeySource.size()
                    ? mSfmKeySource[shadowVertex] : scNull();
                  Size traceOwner = shadowVertex < mSfmKeyOwner.size()
                    ? mSfmKeyOwner[shadowVertex] : scNull();
                  Size traceSourceOuter = traceSource < mN
                    ? outerOf(traceSource) : scNull();
                  Size traceTargetOuter = shadowVertex < mN
                    ? outerOf(shadowVertex) : scNull();
                  std::fprintf(stderr,
                    "[EVENT-INVARIANT] stale split-findmin key child=%lu "
                    "vertex=%ld key=%lld source=%ld/source-outer=%ld "
                    "owner=%ld target-outer=%ld labels=%d/%d live=%d t=%lld\n",
                    (unsigned long)child,
                    (long)(shadowVertex == scNull() ? -1 : shadowVertex),
                    (long long)shadowKey,
                    (long)(traceSource == scNull() ? -1 : traceSource),
                    (long)(traceSourceOuter == scNull() ? -1 : traceSourceOuter),
                    (long)(traceOwner == scNull() ? -1 : traceOwner),
                    (long)(traceTargetOuter == scNull() ? -1 : traceTargetOuter),
                    traceSourceOuter < mNLabel.size()
                      ? mNLabel[traceSourceOuter] : -9,
                    traceTargetOuter < mNLabel.size()
                      ? mNLabel[traceTargetOuter] : -9,
                    traceSourceOuter < mNumNodes && treeLive(traceSourceOuter)
                      ? 1 : 0,
                    (long long)tNow);
#endif
                  if (shadowVertex >= mN) {
                    ++mCounters.algorithmFallbackAttempts;
                    ++mCounters.sfmBoundaryRepairFailures;
                    mLastStatus = STATUS_INVARIANT_VIOLATION;
                    break;
                  }
                  mInnerSfm.setkey(shadowVertex, SplitFindmin::INF);
                  mSfmKeySource[shadowVertex] = scNull();
                  mSfmKeyOwner[shadowVertex] = scNull();
                  mInnerSfm.findmin(leaves[0], &shadowKey, &shadowVertex);
                  if (shadowKey == SplitFindmin::INF)
                    repairValidated = true;
                  continue;
                }

                ++mCounters.sfmBoundaryRepairs;
                if (shadowVertex >= mN) {
                  ++mCounters.sfmBoundaryRepairFailures;
                  break;
                }
                mInnerSfm.setkey(shadowVertex, SplitFindmin::INF);
                mSfmKeySource[shadowVertex] = scNull();
                mSfmKeyOwner[shadowVertex] = scNull();
                if (shadowInChild) {
                  for (Size repairAi = 0;
                       repairAi < mAdj[shadowVertex].size(); ++repairAi) {
                    ++mCounters.sfmBoundaryRepairScans;
                    Size source = mAdj[shadowVertex][repairAi];
                    if (source >= mN) continue;
                    Size currentSource = outerOf(source);
                    if (currentSource == child ||
                        mNLabel[currentSource] != LBL_EVEN ||
                        !treeLive(currentSource)) continue;
                    LL edgeW2 = mUseScaleW ?
                      mW[shadowVertex][repairAi] :
                      2 * mWhat[shadowVertex][repairAi];
                    LL currentSlack = yOfOuter(source, currentSource) +
                      yOfOuter(shadowVertex, child) +
                      oldCommonZ(source, shadowVertex) - edgeW2;
                    if ((Val)currentSlack < mInnerSfm.key(shadowVertex)) {
                      mInnerSfm.decreasekey(shadowVertex, (Val)currentSlack);
                      mSfmKeySource[shadowVertex] = source;
                      mSfmKeyOwner[shadowVertex] = currentSource;
                    }
                  }
                }
                mInnerSfm.findmin(leaves[0], &shadowKey, &shadowVertex);
                if (shadowKey == SplitFindmin::INF)
                  repairValidated = true;
              }
              if (!repairValidated)
                ++mCounters.sfmBoundaryRepairFailures;
#endif
              bool shadowInfinite = shadowKey == SplitFindmin::INF;
#ifndef EVENT_SFM_NO_ORACLE
              bool scanInfinite =
                shadowScanMin == std::numeric_limits<LL>::max();
              bool shadowOwnerValid = false;
              if (!shadowInfinite && shadowVertex < mSfmKeySource.size()) {
                Size shadowSource = mSfmKeySource[shadowVertex];
                Size shadowOwner = mSfmKeyOwner[shadowVertex];
                if (shadowSource < mN && shadowOwner < mNumNodes) {
                  Size currentSource = outerOf(shadowSource);
                  shadowOwnerValid = currentSource == outerOf(shadowOwner) &&
                    mNLabel[currentSource] == LBL_EVEN &&
                    treeLive(currentSource);
                }
              }
              if ((scanInfinite && shadowInfinite) ||
                  (!scanInfinite && !shadowInfinite &&
                   shadowKey == (Val)shadowScanMin))
                ++mCounters.sfmBoundaryExact;
              else {
                ++mCounters.sfmBoundaryStale;
                if (shadowInfinite ||
                    (!scanInfinite && shadowKey > (Val)shadowScanMin))
                  ++mCounters.sfmBoundaryTooHigh;
                else
                  ++mCounters.sfmBoundaryTooLow;
                if (!shadowInfinite && !shadowOwnerValid)
                  ++mCounters.sfmBoundaryInvalidOwner;
                if (!shadowInfinite && !scanInfinite &&
                    shadowKey < (Val)shadowScanMin &&
                    shadowVertex < mSfmKeySource.size()) {
                  Size shadowSource = mSfmKeySource[shadowVertex];
                  bool shadowInChild = false;
                  for (Size shadowLi = 0; shadowLi < leaves.size(); ++shadowLi)
                    if (leaves[shadowLi] == shadowVertex) {
                      shadowInChild = true; break;
                    }
                  if (!shadowInChild) {
                    ++mCounters.sfmBoundaryTopologyLeak;
                  } else if (shadowSource >= mN) {
                    ++mCounters.sfmBoundaryMissingSource;
                  } else {
                    Size currentSource = outerOf(shadowSource);
                    Size currentTarget = outerOf(shadowVertex);
                    if (currentSource == currentTarget) {
                      ++mCounters.sfmBoundaryInternalKey;
                    } else {
                      LL currentSlack = yOfOuter(shadowSource, currentSource) +
                        yOfOuter(shadowVertex, currentTarget) +
                        oldCommonZ(shadowSource, shadowVertex) -
                        w2max(shadowSource, shadowVertex);
                      if ((Val)currentSlack != shadowKey)
                        ++mCounters.sfmBoundaryKeyDrift;
                      else
                        ++mCounters.sfmBoundaryUnexplainedLow;
                    }
                  }
                }
              }
#endif
#ifdef EVENT_SFM_CONSUME
              if (!shadowInfinite && shadowVertex < mN) {
                Size shadowSource = mSfmKeySource[shadowVertex];
                Size shadowOwner = mSfmKeyOwner[shadowVertex];
                if (shadowSource >= mN || shadowOwner >= mNumNodes ||
                    outerOf(shadowSource) != outerOf(shadowOwner)) {
#ifdef EVENTTRACE
                  std::fprintf(stderr, "[EVENT-INVARIANT] stale SFM source/owner\n");
#endif
                  mLastStatus = STATUS_INVARIANT_VIOLATION;
                } else {
                  scheduleEdgeFromEven(shadowSource, shadowVertex,
                    outerOf(shadowSource), crit, tNow,
                    w2max(shadowSource, shadowVertex));
                }
              }
#endif
            }
#endif
          }
          else if (child >= mN && mNLabel[child] == LBL_ODD) {
#ifdef EVENT_SFM_SHADOW
            sfmShadowBeginInner(child);
#endif
            LL fire = 0; if (!checkedAddLL(tNow, zOf(child) / 2, &fire)) {
              mLastStatus = STATUS_UNSUPPORTED_OR_OVERFLOWING_INPUT; break; }
            schedulePracticalEvent(fire, (long long)child, -1, 2,
                                   (long long)child);
          }
          if (mLastStatus == STATUS_SUCCESS &&
              mBParent[child] == scNull() && mNLabel[child] == LBL_ODD &&
              !scheduleMatchedEdgeFromInner(child, crit, tNow)) {
            if (mLastStatus == STATUS_SUCCESS) {
#ifdef EVENTTRACE
              std::fprintf(stderr,
                "[EVENT-INVARIANT] exposed INNER matched scheduling child=%lu\n",
                (unsigned long)child);
#endif
              mLastStatus = STATUS_INVARIANT_VIOLATION;
            }
            break;
          }
        }
#ifdef EVENT_SFM_BOUNDARY_FRONTIER
        mEventSfmFeedFreeBoundary = false;
#endif
      } else {
#ifdef EVENTTRACE
        std::fprintf(stderr, "[EVENT-INVARIANT] dissolution relink failed B=%lu k=%lu par=%ld children:",
          (unsigned long)b, (unsigned long)mBCycle[b].size(),
          (long)(mNParent[b] == scNull() ? -1 : mNParent[b]));
        for (Size ti = 0; ti < mBCycle[b].size(); ++ti)
          std::fprintf(stderr, " %lu(z=%lld,l=%d)",
            (unsigned long)mBCycle[b][ti],
            (long long)(mBCycle[b][ti] < mZ.size() ? mZ[mBCycle[b][ti]] : 0),
            mBCycle[b][ti] < mNLabel.size() ? mNLabel[mBCycle[b][ti]] : -9);
        std::fprintf(stderr, "\n");
#endif
        ++mCounters.algorithmFallbackAttempts;
        ++mCounters.dissolutionRelinkFailures;
        mLastStatus = STATUS_INVARIANT_VIOLATION;
        break;
      }
      continue;
    }

    // Section 5.2 deferred matched-edge GROW.  Payload endpoints are the
    // matched base vertices, so blossom contractions can be resolved lazily
    // at pop time.  If the mate was reached independently, this event is the
    // canceled grow described by schedule(u) and is simply stale.
    if (ev.kind == 7) {
      ++mCounters.eventControlPops;
      Size innerBase = (Size)ev.a, mateBase = (Size)ev.b;
      if (innerBase >= mN || mateBase >= mN ||
          mMate[innerBase] != mateBase || mMate[mateBase] != innerBase)
        continue;
      Size W = outerOfUF(innerBase), M = outerOfUF(mateBase);
      if (W == M || mNLabel[W] != LBL_ODD ||
          mNLabel[M] != LBL_FREE_OF_TREE || !treeLive(W))
        continue;
      LL matchedSlack = yOfOuter(innerBase, W) + yOfOuter(mateBase, M) +
                        oldCommonZ(innerBase, mateBase) -
                        matchedEdgeW2(innerBase, mateBase);
      bool matchedEligible = crit == CRIT2 ? matchedSlack == 0 :
                             crit == CRIT3 ?
                               (matchedSlack == 0 || matchedSlack == -2) :
                               matchedSlack == 0;
      if (!matchedEligible) {
        ++mCounters.eventIneligiblePops;
        continue;
      }
      if (M < mN) lazySetVertexRole(M, LBL_EVEN);
      else lazySetNodeRole(M, LBL_EVEN);
      mNLabel[M] = LBL_EVEN; mNParent[M] = W;
      rememberCompactEventLabel(M);
      mNParEdgeU[M] = mateBase; mNParEdgeV[M] = innerBase;
      mNRoot[M] = mNRoot[W];
#ifdef EVENT_S5_TARGET_KEY_SHADOW
      s5TargetKeyTrackTreeNode(M);
#endif
      scheduleNodeEdges(M, crit, tNow);
      continue;
    }

    // If both ends of a matched edge became INNER before it was eligible,
    // Section 5.2 cancels the grow and schedules blossom(e) at half slack*.
    // Orient endpoint A as the OUTER mate of B.  A's former incoming unmatched
    // edge is then an ordinary OUTER--OUTER edge; re-offer it at the same dual
    // time and let the normal formation/augmentation handler decide whether
    // the two trees coincide.  This is the local alternating-tree transition,
    // and avoids rebuilding (and deterministically repeating) the tied forest.
    if (ev.kind == 8) {
      ++mCounters.eventControlPops;
      Size a = (Size)ev.a, b = (Size)ev.b;
      if (a >= mN || b >= mN || mMate[a] != b || mMate[b] != a)
        continue;
      Size A = outerOfUF(a), B = outerOfUF(b);
      if (A == B || mNLabel[A] != LBL_ODD || mNLabel[B] != LBL_ODD ||
          !treeLive(A) || !treeLive(B))
        continue;
      LL matchedSlack = yOfOuter(a, A) + yOfOuter(b, B) +
                        oldCommonZ(a, b) - matchedEdgeW2(a, b);
      LL eventSlack = ev.owner == 1 ? matchedSlack - 1 : matchedSlack;
#ifdef EVENTTRACE
      std::fprintf(stderr,
        "[EVENT-MATCHED-BLOSSOM] crit=%d time=%lld edge=%lu-%lu slack=%lld event-slack=%lld half=%lld\n",
        crit, (long long)tNow, (unsigned long)a, (unsigned long)b,
        (long long)matchedSlack, (long long)eventSlack,
        (long long)ev.owner);
#endif
      bool matchedEligible = crit == CRIT2 ? eventSlack == 0 :
                             crit == CRIT3 ?
                               (eventSlack == 0 || eventSlack == -2) :
                               eventSlack == 0;
      if (!matchedEligible) {
        ++mCounters.eventIneligiblePops;
        continue;
      }
      // blossom(e) can occur at a half-integral time when slack*(e) is odd.
      // The queue clock stores whole dual units.  Such an event is represented
      // at ceil(fire): A was INNER for the first half and OUTER for the second,
      // whereas applyDualTo just charged it as INNER for the whole unit.  Undo
      // exactly that one-unit discrepancy before changing A's role.  All
      // materialized duals remain integral and equal the paper state at the
      // next observable SearchOne boundary.
      if (ev.owner == 1) {
#ifdef SCALING_LEAF_NOVEC
        std::vector<Size> halfLeaves;
        scalingLfCollect(A, &halfLeaves);
        const std::vector<Size>& leaves = halfLeaves;
#else
        const std::vector<Size>& leaves = mNodeLeaves[A];
#endif
        for (Size li = 0; li < leaves.size(); ++li) {
          Size v = leaves[li];
          if (mLazyDuals) { tsEnsureVertex(v); --mYAcc[v]; }
          else --mY[v];
        }
        if (A >= mN) {
          if (mLazyDuals) { tsEnsureNode(A); mZAcc[A] += 2; }
          else mZ[A] += 2;
        }
      }
      Size incomingTarget = mNParEdgeU[A];
      Size incomingSource = mNParEdgeV[A];
      if (mNParent[A] == scNull() || incomingTarget >= mN ||
          incomingSource >= mN) {
#ifdef EVENTTRACE
        std::fprintf(stderr,
          "[EVENT-INVARIANT] matched blossom missing incoming edge A=%lu\n",
          (unsigned long)A);
#endif
        mLastStatus = STATUS_INVARIANT_VIOLATION;
        break;
      }
      if (A < mN) lazySetVertexRole(A, LBL_EVEN);
      else lazySetNodeRole(A, LBL_EVEN);
      mNLabel[A] = LBL_EVEN; mNParent[A] = B;
      rememberCompactEventLabel(A);
      mNParEdgeU[A] = a; mNParEdgeV[A] = b; mNRoot[A] = mNRoot[B];
#ifdef EVENT_S5_TARGET_KEY_SHADOW
      s5TargetKeyTrackTreeNode(A);
#endif
      scheduleNodeEdges(A, crit, tNow);
      Size incomingOwner = outerOfUF(incomingSource);
      if (mNLabel[incomingOwner] != LBL_EVEN || !treeLive(incomingOwner)) {
#ifdef EVENTTRACE
        std::fprintf(stderr,
          "[EVENT-INVARIANT] matched blossom incoming owner stale owner=%lu label=%d\n",
          (unsigned long)incomingOwner, mNLabel[incomingOwner]);
#endif
        mLastStatus = STATUS_INVARIANT_VIOLATION;
        break;
      }
      scheduleEdgeFromEven(incomingSource, incomingTarget, incomingOwner,
                           crit, tNow, w2max(incomingSource, incomingTarget));
      if (mLastStatus != STATUS_SUCCESS) break;
      continue;
    }

    // Edge events use either a raw target (kinds 0/1) or the source adjacency
    // slot (indexed kinds 5/6). The indexed form avoids rescanning u's row to
    // rediscover the weight and gives future ownership a stable edge identity.
    Size u = (Size)ev.a;
    bool indexedEdge = ev.kind == 5 || ev.kind == 6;
    Size edgeIndex = indexedEdge ? (Size)ev.b : scNull();
    if (u >= mAdj.size() ||
        (indexedEdge && edgeIndex >= mAdj[u].size())) {
#ifdef EVENTTRACE
      std::fprintf(stderr,
        "[EVENT-INVARIANT] corrupt indexed payload u=%lu k=%lu\n",
        (unsigned long)u, (unsigned long)edgeIndex);
#endif
      mLastStatus = STATUS_INVARIANT_VIOLATION;
      break;
    }
    Size w = indexedEdge ? mAdj[u][edgeIndex] : (Size)ev.b;
    if (disjointSearchOnePhase && (mUsedV[u] || mUsedV[w])) {
      ++mCounters.eventStaleRolePops;
      continue;
    }
    // The dissolution split restores UF ownership from authoritative laminar
    // ancestry before clearing parent links, so popped raw endpoints resolve
    // through the maintained Section 5 outer-blossom partition in O(alpha).
    Size E = outerOfUF(u), W = outerOfUF(w);
#ifdef UFCHECK_EVENT_POP
    {
      Size exactE = outerOf(u), exactW = outerOf(w);
      if (exactE != E || exactW != W) {
        Size repU = ufFind(u), repW = ufFind(w);
        Size baseUfW = W < mBBase.size() ? mBBase[W] : scNull();
        Size baseOuter = baseUfW < mN ? outerOf(baseUfW) : scNull();
        std::fprintf(stderr,
          "[UFCHECK-EVENT-POP] t=%lld kind=%d u=%lu outer=%lu uf=%lu "
          "rep=%lu w=%lu outer=%lu uf=%lu rep=%lu uf-active=%d "
          "uf-base=%lu base-outer=%lu owner=%lld nodes=%lu\n",
          (long long)tNow, ev.kind, (unsigned long)u, (unsigned long)exactE,
          (unsigned long)E, (unsigned long)repU, (unsigned long)w,
          (unsigned long)exactW, (unsigned long)W, (unsigned long)repW,
          (W < mBActive.size() && mBActive[W]) ? 1 : 0,
          (unsigned long)baseUfW, (unsigned long)baseOuter,
          ev.owner, (unsigned long)mNumNodes);
        std::abort();
      }
    }
#endif
    if (ev.owner >= 0) {
      Size owner = (Size)ev.owner;
      if (owner >= mNumNodes || E != outerOf(owner)) {
        ++mCounters.eventStaleOwnerPops;
        continue;
      }
    }
    if (mNLabel[E] != LBL_EVEN || !treeLive(E)) {
      ++mCounters.eventStaleRolePops;
      continue; // E no longer a live EVEN node
    }
    if (E == W) { ++mCounters.eventInternalPops; continue; }
    // Revalidate through the criterion itself: Criterion 3 admits both slack
    // 0 and -2, whereas Criterion 2 distinguishes matched from unmatched.
    // A separate single-target slack check here used to discard valid
    // Criterion-3 slack-0 events after the scheduler had emitted them.
    LL eventW2 = indexedEdge
      ? (mUseScaleW ? mW[u][edgeIndex] : 2 * mWhat[u][edgeIndex])
      : w2max(u, w);
    LL eventSlack = yOfOuter(u, E) + yOfOuter(w, W) +
                    oldCommonZ(u, w) - eventW2;
    bool eventEligible = (crit == CRIT1) ? (eventSlack == 0) :
                         (crit == CRIT3) ?
                           (eventSlack == 0 || eventSlack == -2) :
                         ((mMate[u] == w) ? (eventSlack == 0) :
                                                (eventSlack == -2));
    if (!eventEligible) {
      ++mCounters.eventIneligiblePops;
      continue;
    }
    ++mCounters.eventUsefulPops;

    int lblW = mNLabel[W];
    if (lblW != LBL_FREE_OF_TREE && !treeLive(W)) {
      // Do not locally reinterpret a stale consumed-tree label: its parent
      // chain belongs to the old matching. The maximal-set wrapper performs a
      // fresh zero-dual event initialization with exact path vertices masked,
      // which safely exposes every unused branch.
      if (batchSameTimeAugment && augmentedAny) continue;
      lblW = LBL_FREE_OF_TREE;
    }

    // SHRINK transfers ownership through the outer UF/split-findmin discipline.
    // AUGMENTATION terminates this request immediately unless same-time shell
    // saturation is requested. Pure GROW changes only the maintained forest.
    if (lblW == LBL_EVEN) {
      if (lcaOuter(E, W) != scNull()) {
        Size formed = mNumNodes;
        shrink(E, W, u, w);
        tsOnFormOuter(formed);
        // Events owned by the cycle's formerly EVEN nodes remain valid and
        // now resolve through outerOf(owner)==formed. Only formerly ODD nodes
        // have newly become outer and need their incident rows initialized.
        const std::vector<Size>& cycle = mBCycle[formed];
        for (Size ci = 0; ci < cycle.size(); ++ci) {
          Size child = cycle[ci];
          if (mNLabel[child] != LBL_ODD) continue;
#ifdef SCALING_LEAF_NOVEC
          std::vector<Size> leafScratch; scalingLfCollect(child, &leafScratch);
          const std::vector<Size>& leaves = leafScratch;
#else
          const std::vector<Size>& leaves = mNodeLeaves[child];
#endif
          for (Size li = 0; li < leaves.size(); ++li) {
            Size source = leaves[li];
            const LL sourceY = yOfOuter(source, formed);
            for (Size ai = 0; ai < mAdj[source].size(); ++ai) {
              ++mCounters.edgeScans;
              LL edgeW2 = mUseScaleW ? mW[source][ai] : 2 * mWhat[source][ai];
              scheduleEdgeFromEven(source, mAdj[source][ai], formed, crit, tNow,
                                   edgeW2, false, ai, sourceY);
              if (mLastStatus != STATUS_SUCCESS) break;
            }
            if (mLastStatus != STATUS_SUCCESS) break;
          }
          if (mLastStatus != STATUS_SUCCESS) break;
        }
      } else {
#ifdef EVENT_S5_TARGET_KEY_SHADOW
        Size retireTokenE = scNull(), retireTokenW = scNull();
        if (mNRoot[E] != scNull()) {
          Size ro = outerOf(mNRoot[E]);
          if (ro < mBBase.size()) retireTokenE = mBBase[ro];
        }
        if (mNRoot[W] != scNull()) {
          Size ro = outerOf(mNRoot[W]);
          if (ro < mBBase.size()) retireTokenW = mBBase[ro];
        }
#endif
        materializeLazyDuals();
        augmentFrom(E, W, u, w);
#ifdef EVENT_S5_TARGET_KEY_SHADOW
        if (retireTokenE < mN) s5TargetKeyRetireTree(retireTokenE, tNow);
        if (retireTokenW < mN && retireTokenW != retireTokenE)
          s5TargetKeyRetireTree(retireTokenW, tNow);
#endif
        if (!batchSameTimeAugment) {
          mEventStructuresActive = false;
#ifdef SCALING_EVENT_INVARIANTS_EVERY_EVENT
          if (!validateEventInvariantState(crit)) return false;
#endif
          return true;
        }
        if (augmentedAny) ++mCounters.eventBatchedAugmentations;
        augmentedAny = true;
      }
    } else if (lblW == LBL_FREE_OF_TREE) {
      Size bw = mBBase[W];
      Size m = mMate[bw];
      if (m == scNull()) {
#ifdef EVENT_S5_TARGET_KEY_SHADOW
        Size retireTokenE = scNull();
        if (mNRoot[E] != scNull()) {
          Size ro = outerOf(mNRoot[E]);
          if (ro < mBBase.size()) retireTokenE = mBBase[ro];
        }
#endif
        materializeLazyDuals();
        augmentFrom(E, W, u, w);
#ifdef EVENT_S5_TARGET_KEY_SHADOW
        if (retireTokenE < mN) s5TargetKeyRetireTree(retireTokenE, tNow);
#endif
        if (!batchSameTimeAugment) {
          mEventStructuresActive = false;
#ifdef SCALING_EVENT_INVARIANTS_EVERY_EVENT
          if (!validateEventInvariantState(crit)) return false;
#endif
          return true;
        }
        if (augmentedAny) ++mCounters.eventBatchedAugmentations;
        augmentedAny = true;
      } else {
        // Section 5.2 GROW first makes W INNER.  Its matched continuation may
        // become eligible only after 0, 1, or 2 adjustments under Criteria 2/3;
        // do not apply Criterion 1's immediate-mate shortcut in those cases.
        if (W < mN) lazySetVertexRole(W, LBL_ODD);
        else lazySetNodeRole(W, LBL_ODD);
        mNLabel[W] = LBL_ODD; mNParent[W] = E;
        rememberCompactEventLabel(W);
        mNParEdgeU[W] = w; mNParEdgeV[W] = u; mNRoot[W] = mNRoot[E];
#ifdef EVENT_S5_TARGET_KEY_SHADOW
        s5TargetKeyTrackTreeNode(W);
#endif
#ifdef EVENT_SFM_SHADOW
        if (W >= mN) sfmShadowBeginInner(W);
#endif
        if (!scheduleMatchedEdgeFromInner(W, crit, tNow)) {
          if (mLastStatus == STATUS_SUCCESS) {
#ifdef EVENTTRACE
            std::fprintf(stderr,
              "[EVENT-INVARIANT] grow matched scheduling W=%lu\n",
              (unsigned long)W);
#endif
            mLastStatus = STATUS_INVARIANT_VIOLATION;
          }
          break;
        }
        if (W >= mN) { LL fire = 0;
          if (!checkedAddLL(tNow, zOf(W) / 2, &fire)) {
            mLastStatus = STATUS_UNSUPPORTED_OR_OVERFLOWING_INPUT; break; }
          schedulePracticalEvent(fire, (long long)W, -1, 2); }
      }
    }
  }
  if (boundaryReuseRequiresFreshInit ||
      (reuseBoundary && !augmentedAny && mLastSearchDualUnits == 0)) {
    materializeLazyDuals();
    mEventStructuresActive = false;
    ++mCounters.eventBoundaryReuseMisses;
    // The recursive full initialization owns its own timer. Close this reuse
    // reuse attempt explicitly so the fresh initialization is not counted twice.
    if (searchTimer.mCounter) {
      *searchTimer.mCounter += scalingMonotonicNanoseconds() - searchTimer.mStart;
      searchTimer.mCounter = NULL;
    }
    std::vector<BucketPQ::Event>().swap(reuseEvents);
    return edmondsSearchEvent(crit, explicitRoots, maxDualUnits,
                              exhaustDualBudget, batchSameTimeAugment,
                              maxAdjustments);
  }
  if (mLastStatus == STATUS_SUCCESS && exhaustDualBudget &&
      maxDualUnits >= 0 && mLastSearchDualUnits < maxDualUnits &&
      !augmentedAny &&
      (maxAdjustments < 0 || adjustmentEvents < maxAdjustments)) {
    LL remaining = maxDualUnits - mLastSearchDualUnits, target = 0;
    if (!checkedAddLL(tNow, remaining, &target))
      mLastStatus = STATUS_UNSUPPORTED_OR_OVERFLOWING_INPUT;
    else {
      applyDualTo(target, &tNow);
      mLastSearchDualUnits += remaining;
      mCounters.dualUnits += (unsigned long long)remaining;
      if (remaining > 0) ++adjustmentEvents;
    }
  }
#ifdef EVENT_OWNER_MIN_AUDIT
  // Validate the exact retained payload set at the SearchOne boundary and
  // classify which records can carry without a reveal-next repair.
  unsigned long long retainedPayloads = 0;
  for (Size oi = 0; oi < mEventOwnerMinTouched.size(); ++oi) {
    Size owner = mEventOwnerMinTouched[oi];
    if (owner >= mEventOwnerMinEpoch.size() ||
        mEventOwnerMinEpoch[owner] != mEventOwnerMinClock) {
      ++mCounters.eventOwnerPayloadFailures;
      continue;
    }
    const BucketPQ::Event& payload = mEventOwnerMinPrimary[owner];
    ++mCounters.eventOwnerPayloadChecks;
    ++retainedPayloads;
    if (payload.time != mEventOwnerMinTime[owner] ||
        payload.owner != (long long)owner)
      ++mCounters.eventOwnerPayloadFailures;
    Size currentOwner = outerOfUF(owner);
    if (payload.time < tNow) ++mCounters.eventOwnerBoundaryExpired;
    else if (outerOfUF((Size)payload.a) != currentOwner)
      ++mCounters.eventOwnerBoundaryStaleOwner;
    else if (mNLabel[currentOwner] != LBL_EVEN || !treeLive(currentOwner))
      ++mCounters.eventOwnerBoundaryStaleRole;
    else
      ++mCounters.eventOwnerBoundaryPayloads;
  }
  for (Size ti = 0; ti < mEventOwnerMinTieLog.size(); ++ti) {
    const EventOwnerMinTie& tie = mEventOwnerMinTieLog[ti];
    if (tie.owner >= mEventOwnerMinEpoch.size() ||
        mEventOwnerMinEpoch[tie.owner] != mEventOwnerMinClock ||
        tie.version != mEventOwnerMinVersion[tie.owner] ||
        tie.event.time != mEventOwnerMinTime[tie.owner])
      continue;
    ++mCounters.eventOwnerPayloadChecks;
    ++retainedPayloads;
    if (tie.event.owner != (long long)tie.owner)
      ++mCounters.eventOwnerPayloadFailures;
    Size currentOwner = outerOfUF(tie.owner);
    if (tie.event.time < tNow) ++mCounters.eventOwnerBoundaryExpired;
    else if (outerOfUF((Size)tie.event.a) != currentOwner)
      ++mCounters.eventOwnerBoundaryStaleOwner;
    else if (mNLabel[currentOwner] != LBL_EVEN || !treeLive(currentOwner))
      ++mCounters.eventOwnerBoundaryStaleRole;
    else
      ++mCounters.eventOwnerBoundaryPayloads;
  }
  if (retainedPayloads != mEventOwnerMinCurrentPayloads)
    ++mCounters.eventOwnerPayloadFailures;
  for (Size mi = 0; mi < mEventOwnerMatchedFuture.size(); ++mi) {
    if (mEventOwnerMatchedFuture[mi].time < tNow)
      ++mCounters.eventOwnerBoundaryMatchedExpired;
    else
      ++mCounters.eventOwnerBoundaryMatched;
  }
  for (Size ci = 0; ci < mEventOwnerControlFuture.size(); ++ci) {
    if (mEventOwnerControlFuture[ci].time < tNow)
      ++mCounters.eventOwnerBoundaryControlExpired;
    else
      ++mCounters.eventOwnerBoundaryControl;
  }
#endif
#ifdef EVENT_S5_TARGET_KEY_SHADOW
#ifndef EVENT_S5_TARGET_KEY_SKIP_EXACT_AUDIT
  auditS5TargetKeys(crit, tNow);
#endif
#endif
  materializeLazyDuals();
#ifdef SCALING_EVENT_INVARIANTS_EVERY_EVENT
  if (mLastStatus == STATUS_SUCCESS)
    validateEventInvariantState(crit);
#endif
  mEventStructuresActive = false;
  mEventBoundaryReuseReady = false;
  mEventBoundaryReuseEvents.clear();
  if (mEventBoundaryReuse && mLastStatus == STATUS_SUCCESS && crit == CRIT2 &&
      !augmentedAny && mLastSearchDualUnits > 0 && maxAdjustments == 1) {
    for (Size i = 0; i < mEventDeferredEdgeEvents.size(); ++i) {
      BucketPQ::Event saved = mEventDeferredEdgeEvents[i];
      if (saved.time != tNow) continue;
      saved.time = 0;
      mEventBoundaryReuseEvents.push_back(saved);
    }
    mEventBoundaryReuseReady = !mEventBoundaryReuseEvents.empty();
    mCounters.eventBoundaryEventsSaved += mEventBoundaryReuseEvents.size();
    if (mEventBoundaryReuseEvents.size() > mCounters.eventBoundaryMaxEventsSaved)
      mCounters.eventBoundaryMaxEventsSaved = mEventBoundaryReuseEvents.size();
  }
#ifdef EVENT_S5_TARGET_KEY_SHADOW
#ifdef EVENT_S5_TARGET_KEY_DRIVER_AUDIT
  if (mS5TargetKeyActive && mEventBoundaryReuseReady) {
    LL nextBase = 0;
    if (!checkedAddLL(mS5TargetKeyClockBase, tNow, &nextBase)) {
      mLastStatus = STATUS_UNSUPPORTED_OR_OVERFLOWING_INPUT;
      mS5TargetKeyActive = false;
      mEventBoundaryReuseReady = false;
      mEventBoundaryReuseEvents.clear();
    } else {
      mS5TargetKeyClockBase = nextBase;
    }
  } else {
    mS5TargetKeyActive = false;
  }
#else
  mS5TargetKeyActive = false;
#endif
#endif
  if (mLastStatus == STATUS_SUCCESS && crit == CRIT2) {
    if (!augmentedAny && mLastSearchDualUnits > 0) {
      ++mCounters.eventCriterion2AdjustmentOnlyReturns;
      mPreviousCriterion2AdjustmentOnly = true;
    } else if (augmentedAny && mLastSearchDualUnits == 0) {
      ++mCounters.eventCriterion2ZeroDualAugmentReturns;
      if (previousCriterion2AdjustmentOnly)
        ++mCounters.eventCriterion2AdjustmentToZeroAugmentPairs;
    }
  }
  if (disjointSearchOnePhase)
    std::fill(mUsedV.begin(), mUsedV.begin() + mN, (char)0);
  return augmentedAny;
}

inline bool ScalingMatching::eventShellMcmAtCurrentDual() {
  ScalingNanosecondGuard shellMcmTimer(&mCounters.shellMcmPhaseNanoseconds);
  ++mCounters.shellMcmCalls;
#ifdef SCPHASE
  const unsigned long long shellPassesAtEntry = mCounters.shellMcmPasses;
  const unsigned long long shellAugPhasesAtEntry =
    mCounters.shellMcmAugmentingPhases;
#endif
  if (!mShellActive ||
      (mShellCstar == scNull() && mShellMarked.empty())) {
    mLastStatus = STATUS_INVARIANT_VIOLATION;
    return false;
  }
  Size initialCard = 0, freeVertices = 0;
  std::vector<Size> phaseRoots;
  phaseRoots.reserve(mShellMarked.size());
  bool initialSymmetric = true;
  const unsigned long long initialScanStart = scalingMonotonicNanoseconds();
  for (Size i = 0; i < mShellMarked.size(); ++i) {
    Size v = mShellMarked[i];
    if (v >= mN || !mInShell[v]) continue;
    Size w = mMate[v];
    if (w == scNull()) {
      ++freeVertices;
      phaseRoots.push_back(v);
      continue;
    }
    if (w >= mN || !mInShell[w] || mMate[w] != v) initialSymmetric = false;
    if (v < w) ++initialCard;
  }
  mCounters.shellMcmInitialScanNanoseconds +=
    scalingMonotonicNanoseconds() - initialScanStart;
  if (!initialSymmetric) {
    mLastStatus = STATUS_INVARIANT_VIOLATION;
    return false;
  }
  const unsigned long long maxAdditional = freeVertices / 2;
  unsigned long long totalAdded = 0;
  bool augmentedAny = false;

  // Exact phase 0: take a maximal set of length-one augmenting paths in the
  // fixed eligible shell.  A free--free eligible edge is a shortest possible
  // augmenting path, so all such vertex-disjoint edges belong to one valid
  // cardinality phase.  The old driver made the general event forest discover
  // these paths amid thousands of later phases even when M entered the shell
  // empty.  Apply them directly to canonical M; longer paths and the terminal
  // no-path certificate remain the responsibility of the unchanged event
  // phases below.
  const unsigned long long greedyStart = scalingMonotonicNanoseconds();
  const unsigned long long beforeGreedy = mCounters.augmentations;
  for (Size i = 0; i < mShellMarked.size(); ++i) {
    Size u = mShellMarked[i];
    if (u >= mN || !mInShell[u] || mMate[u] != scNull()) continue;
    Size E = outerOf(u);
    if (E >= mNumNodes || mBBase[E] != u) continue;
    for (Size k = 0; k < mAdj[u].size(); ++k) {
      Size v = mAdj[u][k];
      if (v >= mN || v == u || !mInShell[v] || mMate[v] != scNull())
        continue;
      Size W = outerOf(v);
      if (W == E || W >= mNumNodes || mBBase[W] != v) continue;
      if (slackZ(u, v) != 0) continue;
      // Present the two free outer nodes as singleton alternating trees so
      // augmentFrom also performs the required blossom-base rotations.
      mNLabel[E] = LBL_EVEN; mNRoot[E] = E; mNParent[E] = scNull();
      mNParEdgeU[E] = scNull(); mNParEdgeV[E] = scNull();
      mNLabel[W] = LBL_EVEN; mNRoot[W] = W; mNParent[W] = scNull();
      mNParEdgeU[W] = scNull(); mNParEdgeV[W] = scNull();
      rememberCompactEventLabel(E); rememberCompactEventLabel(W);
      augmentFrom(E, W, u, v);
      break;
    }
  }
  const unsigned long long greedyAdded =
    mCounters.augmentations - beforeGreedy;
  mCounters.shellMcmGreedyNanoseconds +=
    scalingMonotonicNanoseconds() - greedyStart;
  if (greedyAdded > maxAdditional) {
    mLastStatus = STATUS_INVARIANT_VIOLATION;
    return false;
  }
  if (greedyAdded != 0) {
    ++mCounters.shellMcmPasses;
    ++mCounters.shellMcmAugmentingPhases;
    mCounters.shellMcmPhaseAugmentations += greedyAdded;
    if (greedyAdded > mCounters.shellMcmMaxPhaseAugmentations)
      mCounters.shellMcmMaxPhaseAugmentations = greedyAdded;
    totalAdded = greedyAdded;
    augmentedAny = true;
  }
  if (greedyAdded != 0) {
    const unsigned long long rootFilterStart = scalingMonotonicNanoseconds();
    Size retained = 0;
    for (Size i = 0; i < phaseRoots.size(); ++i) {
      Size v = phaseRoots[i];
      if (v < mN && mMate[v] == scNull()) phaseRoots[retained++] = v;
    }
    phaseRoots.resize(retained);
    mCounters.shellMcmRootFilterNanoseconds +=
      scalingMonotonicNanoseconds() - rootFilterStart;
  }
  while (totalAdded < maxAdditional) {
    ++mCounters.shellMcmPasses;
    const unsigned long long before = mCounters.augmentations;
    const unsigned long long searchStart = scalingMonotonicNanoseconds();
    bool augmented = edmondsSearchEvent(
      CRIT1, &phaseRoots, 0, false, true, 0);
    const unsigned long long searchElapsed =
      scalingMonotonicNanoseconds() - searchStart;
    const unsigned long long added = mCounters.augmentations - before;
    if (mLastStatus != STATUS_SUCCESS) return augmentedAny;
    if (!augmented) {
      mCounters.shellMcmCertificationNanoseconds += searchElapsed;
      ++mCounters.shellMcmCertificationPhases;
      if (added != 0) mLastStatus = STATUS_INVARIANT_VIOLATION;
      break;
    }
    mCounters.shellMcmAugmentSearchNanoseconds += searchElapsed;
    if (added == 0 || added > maxAdditional - totalAdded) {
      mLastStatus = STATUS_INVARIANT_VIOLATION;
      return augmentedAny;
    }
    ++mCounters.shellMcmAugmentingPhases;
    mCounters.shellMcmPhaseAugmentations += added;
    if (added > mCounters.shellMcmMaxPhaseAugmentations)
      mCounters.shellMcmMaxPhaseAugmentations = added;
    totalAdded += added;
    augmentedAny = true;
    const unsigned long long rootFilterStart = scalingMonotonicNanoseconds();
    Size retained = 0;
    for (Size i = 0; i < phaseRoots.size(); ++i) {
      Size v = phaseRoots[i];
      if (v < mN && mMate[v] == scNull()) phaseRoots[retained++] = v;
    }
    phaseRoots.resize(retained);
    mCounters.shellMcmRootFilterNanoseconds +=
      scalingMonotonicNanoseconds() - rootFilterStart;
  }

  // One final O(|shell|) audit replaces the former full-shell scan after every
  // augmenting phase. Cardinality must increase by exactly the number of
  // augmentFrom transitions counted inside those phases.
  Size finalCard = 0;
  bool finalSymmetric = true;
  const unsigned long long auditStart = scalingMonotonicNanoseconds();
  for (Size i = 0; i < mShellMarked.size(); ++i) {
    Size v = mShellMarked[i];
    if (v >= mN || !mInShell[v]) continue;
    Size w = mMate[v];
    if (w != scNull() &&
        (w >= mN || !mInShell[w] || mMate[w] != v)) finalSymmetric = false;
    if (w != scNull() && v < w) ++finalCard;
  }
  mCounters.shellMcmAuditNanoseconds +=
    scalingMonotonicNanoseconds() - auditStart;
  if (!finalSymmetric || finalCard < initialCard ||
      (unsigned long long)(finalCard - initialCard) != totalAdded)
    mLastStatus = STATUS_INVARIANT_VIOLATION;
#ifdef SCPHASE
  if ((mCounters.shellMcmCalls & 2047ULL) == 0) {
    std::fprintf(stderr,
      "[SCPHASE] shell-mcm calls=%llu shell=%lu free=%lu added=%llu "
      "call-passes=%llu call-aug-phases=%llu total-passes=%llu "
      "total-aug=%llu max-phase=%llu mask-s=%.3f shell-search-s=%.3f "
      "phase-mcm-s=%.3f forest-s=%.3f cleanup-s=%.3f "
      "dism-s=%.3f chain-s=%.3f free-s=%.3f shell-s=%.3f progress-s=%.3f "
      "dism-calls=%llu dism-iters=%llu cleanup-calls=%llu "
      "cleanup-visits=%llu init-visits=%llu clock=%.3fs\n",
      mCounters.shellMcmCalls, (unsigned long)mShellMarked.size(),
      (unsigned long)freeVertices, totalAdded,
      mCounters.shellMcmPasses - shellPassesAtEntry,
      mCounters.shellMcmAugmentingPhases - shellAugPhasesAtEntry,
      mCounters.shellMcmPasses, mCounters.shellMcmPhaseAugmentations,
      mCounters.shellMcmMaxPhaseAugmentations,
      mCounters.shellMaskNanoseconds / 1e9,
      mCounters.shellSearchNanoseconds / 1e9,
      mCounters.shellMcmPhaseNanoseconds / 1e9,
      mCounters.eventForestInitializationNanoseconds / 1e9,
      mCounters.zeroRootCleanupNanoseconds / 1e9,
      mCounters.dismantleNanoseconds / 1e9,
      mCounters.dismantleChainNanoseconds / 1e9,
      mCounters.dismantleFreeCountNanoseconds / 1e9,
      mCounters.dismantleShellNanoseconds / 1e9,
      mCounters.dismantleProgressNanoseconds / 1e9,
      mCounters.dismantleCalls, mCounters.dismantleIterations,
      mCounters.zeroRootCleanupCalls, mCounters.zeroRootCleanupNodeVisits,
      mCounters.eventInitializationNodeVisits,
      (double)clock()/CLOCKS_PER_SEC);
    std::fflush(stderr);
  }
#endif
  return augmentedAny;
}

// ===========================================================================
// PHASE 2/3 — Liquidationist scaling loop (paper Fig 5). Contract:
// docs/Scaling_Paper_Exact_Spec.md. Works on doubled-even scale weights mW.
// ===========================================================================

// leaves (original vertices) under node B, in no particular order.
inline void ScalingMatching::collectBlossomVerts(Size B, std::vector<Size>* out) const {
  if (B < mN) { out->push_back(B); return; }
#ifdef CYCDBG
  if (out->size() > mNumNodes + 8) {
    std::fprintf(stderr, "[CYCDBG] collectBlossomVerts RUNAWAY at B=%lu (out=%lu > nodes=%lu) -- cyclic mBCycle!\n",
                 (unsigned long)B, (unsigned long)out->size(), (unsigned long)mNumNodes); std::fflush(stderr);
    return;
  }
#endif
  for (Size i = 0; i < mBCycle[B].size(); ++i)
    collectBlossomVerts(mBCycle[B][i], out);
}

inline Size ScalingMatching::blossomSize(Size B) const {
  if (B < mN) return 1;
  Size s = 0;
#ifdef CYCDBG
  for (Size i = 0; i < mBCycle[B].size(); ++i) {
    s += blossomSize(mBCycle[B][i]);
    if (s > mNumNodes + 8) { std::fprintf(stderr, "[CYCDBG] blossomSize RUNAWAY at B=%lu (s=%lu > nodes=%lu) -- cyclic mBCycle!\n", (unsigned long)B, (unsigned long)s, (unsigned long)mNumNodes); std::fflush(stderr); return s; }
  }
#else
  for (Size i = 0; i < mBCycle[B].size(); ++i) s += blossomSize(mBCycle[B][i]);
#endif
  return s;
}

// ---- HYBRID §4.2 major-path decomposition (pure structure; reuses mNodeLeaves for O(1)
//      sizes instead of recursing blossomSize). A child B1 of B is a MAJOR CHILD iff
//      |B1| > |B|/2; since two children can't both exceed half, there is at most one. ----
inline Size ScalingMatching::majorChild(Size B) const {
  if (B < mN) return scNull();                       // a vertex (trivial) has no children
  Size half = nodeLeafCount(B) / 2;                  // |B|/2  (integer; |B| is odd)
  const std::vector<Size>& cyc = mBCycle[B];
  for (Size i = 0; i < cyc.size(); ++i) {
    Size c = cyc[i];
    Size sz = (c < mN) ? 1 : nodeLeafCount(c);
    if (sz > half) return c;                          // strictly more than half => major
  }
  return scNull();
}

// R is a MAJOR PATH ROOT iff it is NOT the major child of its parent (so it heads a path).
// The dissolution subtree's own root is always a major-path root (no parent within it).
inline bool ScalingMatching::isMajorPathRoot(Size R) const {
  Size p = (mOldForestSeparated && R < mOldParent.size())
             ? mOldParent[R] : mBParent[R];
  // Zero-z records are retained only as immutable membership containers in the
  // separated representation; they are dissolved and therefore absent from
  // Omega_0.  Find R's parent in the LIVE old-blossom forest.
  if (mOldForestSeparated) {
    Size guard = 0;
    while (p != scNull() && p >= mN && p < mZ.size() && mZ[p] == 0 &&
           ++guard <= mNumNodes + 2)
      p = p < mOldParent.size() ? mOldParent[p] : scNull();
  }
  if (p == scNull() || p >= mNumNodes) return true;   // root of the (sub)forest
  if (mOldForestSeparated) {
    Size rsz = nodeLeafCount(R);
    Size psz = nodeLeafCount(p);
    return rsz <= psz / 2;
  }
  return majorChild(p) != R;
}

// Postorder list of major-path roots in B's subtree. Postorder so DismantlePath(R) is
// called on descendants before ancestors (the §4.2 precondition: when DismantlePath(R)
// runs, the only undissolved blossoms left in R are those on P(R) itself).
inline void ScalingMatching::collectMajorPathRoots(Size B, std::vector<Size>* out) const {
  if (B < mN) return;                                 // vertices have no sub-blossoms
  // An inherited blossom chain may be linear in n. Use an explicit DFS stack
  // so the paper's postorder DismantlePath schedule does not depend on the
  // process call-stack limit at million-vertex scale.
  std::vector<Size> nodes;
  std::vector<Size> nextChild;
  nodes.push_back(B);
  nextChild.push_back(0);
  while (!nodes.empty()) {
    Size node = nodes.back();
    Size& next = nextChild.back();
    const std::vector<Size>& cyc = mBCycle[node];
    if (next < cyc.size()) {
      Size child = cyc[next++];
      if (child >= mN) {
        nodes.push_back(child);
        nextChild.push_back(0);
      }
      continue;
    }
    if ((!mOldForestSeparated || (node < mZ.size() && mZ[node] != 0)) &&
        isMajorPathRoot(node)) out->push_back(node);
    nodes.pop_back();
    nextChild.pop_back();
  }
}

// The undissolved blossoms on major path P(R): walk R -> majorChild -> ... collecting
// each still-undissolved blossom node. Order: R first (outermost), down to the deepest.
// Caller (stage 2) wants B1 subset ... subset Bl; we push outermost-first then the caller
// uses the set (sum of z is order-independent). A node is undissolved iff id>=mN and z!=0.
inline void ScalingMatching::undissolvedPathChain(Size R, std::vector<Size>* out) const {
  Size cur = R, guard = 0;
  while (cur != scNull() && cur >= mN) {
    if (++guard > mNumNodes + 8) break;             // cycle guard (matches the slackZ guards)
    if (mZ[cur] != 0) out->push_back(cur);          // still an undissolved (z>0) blossom
    cur = majorChild(cur);                          // descend the major path
  }
}

// §4.2 DismantlePath STAGE 2 (the |F|=1 case): liquidate the remaining blossoms on P(R),
// then run a bounded Criterion-1 search from omega for exactly T dual adjustments. Per
// Lemma 4.4/4.6 this restores y(omega)/yz(R) and finds no augmenting path. We reuse the
// existing liquidate() (folds z into member y) and edmondsSearch(CRIT1, maxAdjust=T).
inline bool ScalingMatching::dismantleStage2(Size R) {
  std::vector<Size> chain;
  undissolvedPathChain(R, &chain);
  if (chain.empty()) {
#ifdef HYBPROF
    ++gHybStage2Empty;
#endif
    return true;                                    // nothing undissolved on P(R)
  }
  // Stage 2 is legal only when the LIVE undissolved major-path region contains
  // one unmatched vertex.  R itself may already have reached z=0 in stage 1;
  // in that case chain.front(), not the retained historical container R, is
  // the outermost live blossom on P(R).
  Size omega = scNull(), freeCount = 0;
  Size liveRegion = chain.front();
  if (liveRegion >= mNodeLeaves.size()) {
    if (mHybridCanonical) mLastStatus = STATUS_INVARIANT_VIOLATION;
    return false;
  }
  std::vector<Size> liveLeaves;
  collectNodeLeaves(liveRegion, &liveLeaves);
  for (Size i = 0; i < liveLeaves.size(); ++i) {
    Size v = liveLeaves[i];
    if (v < mN && mMate[v] == scNull()) { omega = v; ++freeCount; }
  }
#ifdef HYBPROF
  if (freeCount == 0) ++gHybStage2Free0;
  else if (freeCount == 1) ++gHybStage2Free1;
  else ++gHybStage2FreeMany;
  if (freeCount != 1 && std::getenv("MBX_TRACE_DPATH")) {
    std::fprintf(stderr, "[DPATH-STAGE2] scale=%lu R=%lu free=%lu chain=",
      (unsigned long)mAuditScale, (unsigned long)R, (unsigned long)freeCount);
    for (Size i = 0; i < chain.size(); ++i)
      std::fprintf(stderr, " %lu(z=%lld)", (unsigned long)chain[i],
                   (long long)mZ[chain[i]]);
    std::fprintf(stderr, " leaves=");
    for (Size i = 0; i < liveLeaves.size(); ++i) {
      Size q = liveLeaves[i];
      std::fprintf(stderr, " %lu(m=%ld)", (unsigned long)q,
        (long)(q < mN && mMate[q] != scNull() ? mMate[q] : -1));
    }
    std::fprintf(stderr, "\n");
  }
#endif
  if (freeCount != 1 || omega == scNull()) {
    if (mHybridCanonical) mLastStatus = STATUS_INVARIANT_VIOLATION;
    return false;
  }
#ifdef HYBPROF
  ++gHybStage2Calls;
#endif
  ++mCounters.dismantleStage2Calls;
  LL T = 0;
  for (Size i = 0; i < chain.size(); ++i) {
    LL next = 0;
    if (!checkedAddLL(T, mZ[chain[i]] / 2, &next)) {
      mLastStatus = STATUS_UNSUPPORTED_OR_OVERFLOWING_INPUT;
      return false;
    }
    T = next;
  }
  // Lemma 4.6: liquidating the live inherited chain raises y(omega) by T,
  // and the exact-T singleton PQSearch must lower it by exactly T again.
  // Snapshot the pre-liquidation value so the implementation checks this
  // restoration directly instead of inferring it from the search budget.
  LL omegaYBeforeStage2 = mY[omega];
  // Liquidate innermost-first is not required (liquidate folds z into member y and
  // dissolves; doing all of them makes the path vertices bare). Process in collected order.
  for (Size i = 0; i < chain.size(); ++i)
    if (chain[i] >= mN && mZ[chain[i]] != 0) {
      if (mOldForestSeparated) liquidateInherited(chain[i]);
      else liquidate(chain[i]);
    }
#ifdef HYBPROF
  if (std::getenv("MBX_TRACE_DPATH")) {
    std::fprintf(stderr, "[STAGE2-BEFORE-PQ] R=%lu omega=%lu T=%lld",
      (unsigned long)R, (unsigned long)omega, (long long)T);
    for (Size u = 0; u < mN; ++u) for (Size k = 0; k < mAdj[u].size(); ++k) {
      Size v = mAdj[u][k];
      if (u < v && v < mN && slackZ(u, v) <= 0)
        std::fprintf(stderr, " e%lu-%lu:s%lld:y%lld/%lld",
          (unsigned long)u, (unsigned long)v, (long long)slackZ(u, v),
          (long long)mY[u], (long long)mY[v]);
    }
    std::fprintf(stderr, "\n");
  }
#endif
  // PQSearch({omega}) Criterion 1, exactly T units of dual adjustment. Other unmatched
  // vertices may exist elsewhere in the graph, so pass the singleton root set
  // explicitly rather than relying on global-free inference. Event-free intervals
  // are delta-batched but charged in full, including the final interval with no event.
  std::vector<Size> roots(1, omega);
  bool augmented = false;
  mLastSearchDualUnits = 0;
  if (T > 0) {
    SearchRequest request;
    request.roots = &roots; request.criterion = CRIT1;
    request.dualBudget = T; request.maxAdjustments = T;
    request.batchEventFreeTime = true; request.exhaustBudget = true;
    request.stoppingRule = STOP_AT_DUAL_BUDGET;
    SearchResult result = executeSearch(
      request, mEventSearch ? BACKEND_EVENT : BACKEND_SCAN);
    augmented = result.augmented;
  }
  if (mOldForestSeparated) dissolveZeroCurrentRoots();
#ifdef HYBPROF
  if (std::getenv("MBX_TRACE_DPATH")) {
    std::fprintf(stderr, "[STAGE2-AFTER-PQ] R=%lu spent=%lld",
      (unsigned long)R, (long long)mLastSearchDualUnits);
    for (Size u = 0; u < mN; ++u) for (Size k = 0; k < mAdj[u].size(); ++k) {
      Size v = mAdj[u][k];
      if (u < v && v < mN && slackZ(u, v) < 0)
        std::fprintf(stderr, " e%lu-%lu:s%lld:y%lld/%lld",
          (unsigned long)u, (unsigned long)v, (long long)slackZ(u, v),
          (long long)mY[u], (long long)mY[v]);
    }
    std::fprintf(stderr, "\n");
  }
#endif
#ifdef HYBPROF
  if (augmented && std::getenv("MBX_TRACE_STAGE2_AUG")) {
    std::fprintf(stderr,
      "[STAGE2-AUG] scale=%lu R=%lu omega=%lu T=%lld spent=%lld old=[%lu,%lu)\n",
      (unsigned long)mAuditScale, (unsigned long)R, (unsigned long)omega,
      (long long)T, (long long)mLastSearchDualUnits,
      (unsigned long)mAuditOldBegin, (unsigned long)mAuditOldEnd);
    for (Size v = 0; v < mN; ++v)
      if (mMate[v] != scNull() || v == omega)
        std::fprintf(stderr, "  v=%lu mate=%ld y=%lld outer=%lu\n",
          (unsigned long)v, (long)(mMate[v] == scNull() ? -1 : mMate[v]),
          (long long)mY[v], (unsigned long)outerOf(v));
  }
#endif
  // Lemma 4.6 says this search does not find an augmenting path.
  if (augmented || mLastSearchDualUnits != T ||
      mY[omega] != omegaYBeforeStage2) {
    if (mHybridCanonical) mLastStatus = STATUS_INVARIANT_VIOLATION;
    return false;
  }
  return true;
}

// §4.2 ShellSearch(C*,D*): run ONE Criterion-1 Edmonds search CONFINED to the atomic shell
// G(C*,D*) = leaves(C*)\leaves(D*), with the per-dual-adjust C*/D* unit translation. Sets the
// shell masks (mInShell = search region C*\D*; mInDstar = D*), runs a single bounded confined
// edmondsSearch, then clears the masks. NON-persistent (mShellActive needs the rebuild path's
// root-filter; persist path roots differently) -- so force the from-scratch search here.
inline Size ScalingMatching::nextLiveOldParent(Size B, Size pathRoot) const {
  if (B == scNull() || B >= mNumNodes) return scNull();
  Size p = scNull();
  if (mOldForestSeparated && B < mOldParent.size()) p = mOldParent[B];
  else if (B < mBParent.size()) p = mBParent[B];
  Size guard = 0;
  while (p != scNull() && p >= mN && ++guard <= mNumNodes + 2) {
    if (p < mZ.size() && mZ[p] != 0) return p;
    if (p == pathRoot) return scNull();
    if (mOldForestSeparated && p < mOldParent.size()) p = mOldParent[p];
    else if (p < mBParent.size()) p = mBParent[p];
    else p = scNull();
  }
  return scNull();
}

inline Size ScalingMatching::nextLiveOldMajorChild(Size B) const {
  if (B == scNull() || B < mN || B >= mNumNodes) return scNull();
  Size c = majorChild(B), guard = 0;
  while (c != scNull() && c >= mN && ++guard <= mNumNodes + 2) {
    if (c < mZ.size() && mZ[c] != 0) return c;
    c = majorChild(c);
  }
  return scNull();
}

inline bool ScalingMatching::resolveOldShellBounds(Size C, Size D, Size pathRoot,
                                                   Size* Cstar, Size* Dstar) const {
  if (!Cstar || !Dstar || C == scNull() || C < mN) return false;
  Size co = (C < mZ.size() && mZ[C] != 0)
              ? C : ((C == pathRoot) ? scNull()
                                      : nextLiveOldParent(C, pathRoot));
  Size di = scNull();
  if (D != scNull() && D >= mN)
    di = (D < mZ.size() && mZ[D] != 0) ? D : nextLiveOldMajorChild(D);
  if (co == scNull() || co < mN || co == di) return false;
  *Cstar = co; *Dstar = di;
  return true;
}

inline void ScalingMatching::dissolveZeroCurrentRoots() {
  ScalingNanosecondGuard cleanupTimer(&mCounters.zeroRootCleanupNanoseconds);
  ++mCounters.zeroRootCleanupCalls;
  bool changed = true; Size guard = 0;
  while (changed && ++guard <= mNumNodes + 2) {
    changed = false;
    for (Size b = mAuditOldEnd; b < mNumNodes; ++b) {
      ++mCounters.zeroRootCleanupNodeVisits;
      if (b >= mBParent.size() || b >= mZ.size() || b >= mBCycle.size()) continue;
      if (b < mBActive.size() && mBActive[b] &&
          mBParent[b] == scNull() && mZ[b] == 0 && !mBCycle[b].empty()) {
#ifdef HYBPROF
        if (std::getenv("MBX_TRACE_DPATH"))
          std::fprintf(stderr, "[ZERO-CURRENT-DISSOLVE] B=%lu children=",
            (unsigned long)b);
        if (std::getenv("MBX_TRACE_DPATH")) {
          for (Size ci = 0; ci < mBCycle[b].size(); ++ci)
            std::fprintf(stderr, " %lu(z=%lld)",
              (unsigned long)mBCycle[b][ci],
              (long long)(mBCycle[b][ci] < mZ.size() ? mZ[mBCycle[b][ci]] : 0));
          std::fprintf(stderr, "\n");
        }
#endif
        dissolve(b);
        // Retire the historical record.  dissolve() intentionally preserves
        // cycle data for generic search bookkeeping, but this root-cleanup
        // scans by structure.  Leaving the cycle populated lets the same stale
        // record dissolve again after its former leaves have joined a newer
        // blossom, tearing those leaves out of the new positive-z root.
        mBCycle[b].clear();
        if (b < mBCycEdgeU.size()) mBCycEdgeU[b].clear();
        if (b < mBCycEdgeV.size()) mBCycEdgeV[b].clear();
        if (b < mNLabel.size()) mNLabel[b] = LBL_FREE_OF_TREE;
        changed = true;
      }
    }
  }
  if (changed && mHybridCanonical)
    mLastStatus = STATUS_UNEXPECTED_NO_PROGRESS;
}

inline bool ScalingMatching::shellSearch(Size Cstar, Size Dstar, Size pathRoot,
                                         Size iterEpoch, Size callOwner) {
  if (Cstar == scNull() || Cstar < mN) return false;     // C* must be an undissolved blossom
  if ((Size)mInShell.size() < mN) mInShell.assign(mN, 0);
  if ((Size)mInDstar.size() < mN) mInDstar.assign(mN, 0);
  if (!mDynamicShell) {
  unsigned long long shellMaskStart = scalingMonotonicNanoseconds();
  mShellMarked.clear();
  std::vector<Size> cLeaves, dLeaves;
  collectNodeLeaves(Cstar, &cLeaves);
  if (Dstar != scNull() && Dstar >= mN) collectNodeLeaves(Dstar, &dLeaves);
  // mark D* leaves
  if (Dstar != scNull() && Dstar >= mN)
    for (Size i = 0; i < dLeaves.size(); ++i) {
      Size v = dLeaves[i];
      if (v < mN) { mInDstar[v] = 1; mShellMarked.push_back(v); }
    }
  // mark C*\D* leaves (the searchable shell)
  for (Size i = 0; i < cLeaves.size(); ++i) {
    Size v = cLeaves[i];
    if (v < mN && !mInDstar[v]) {
      mInShell[v] = 1; mShellMarked.push_back(v);
    }
  }
  ++mCounters.shellBoundaryConfigurations;
  mCounters.shellMaskVertices += mShellMarked.size();
  mCounters.shellMaskNanoseconds += scalingMonotonicNanoseconds() - shellMaskStart;
#ifdef HYBPROF
  ++gHybShellSearches;
#endif
  mShellActive = true; mShellCstar = Cstar; mShellDstar = Dstar;
  // Confined Criterion-1 search. NON-persistent rebuild path (persistAcrossAdjust=false) so the
  // shell root-filter applies each iteration; maximalAugment=true (MCM in the shell). ADJUSTMENT
  // CAP: per Lemma 4.7/Cor 4.8 #dual adjustments in ShellSearch = O(n(C*\D*)); z(C*) starts at
  // <= the +6 slack bound so the unit translation (z-=2/adjust) dissolves C* within z(C*)/2
  // adjustments. Cap = z(C*)/2 + n(shell) + 2 (paper-justified upper bound, NOT a shortcut) so a
  // logic bug in the translation/dissolution cannot spin forever. NOTE: if this cap is ever hit
  // in practice it indicates a real bug (the search should self-halt first) -- gate must verify.
  LL cap = mZ[Cstar] / 2 + (LL)nodeLeafCount(Cstar) + 2;
  SearchRequest request;
  request.criterion = CRIT1; request.region = REGION_SHELL;
  request.maximalAugment = true; request.maxAdjustments = cap;
  request.stoppingRule = STOP_AT_SHELL_EVENT; request.shellTranslation = true;
  unsigned long long shellStart = scalingMonotonicNanoseconds();
  ++mCounters.shellSearchCalls;
  SearchResult shellResult = executeSearch(
    request, mEventSearch ? BACKEND_EVENT : BACKEND_SCAN);
  const unsigned long long shellElapsed =
    scalingMonotonicNanoseconds() - shellStart;
  mCounters.shellSearchNanoseconds += shellElapsed;
  mCounters.shellMcmNanoseconds += shellElapsed;
  bool aug = shellResult.augmented;
#ifdef HYBPROF
  if (aug) ++gHybShellAugCalls; else ++gHybShellNoAugCalls;
#endif
  // clear shell state
  for (Size i = 0; i < cLeaves.size(); ++i) { Size v = cLeaves[i]; if (v < mN) { mInShell[v]=0; mInDstar[v]=0; } }
  if (Dstar != scNull() && Dstar >= mN)
    for (Size i = 0; i < dLeaves.size(); ++i) { Size v = dLeaves[i]; if (v < mN) mInDstar[v]=0; }
  mShellMarked.clear();
  mShellActive = false; mShellCstar = scNull(); mShellDstar = scNull();
  if (shellResult.status != STATUS_SUCCESS) return false;
  return aug;
  }

  // Paper-faithful dynamic shell. A search that exhausts C* or D* returns to
  // this wrapper; we physically dissolve that old boundary, update C*/D*,
  // rebuild G(C*,D*), check same-iteration absorption, and continue.
  if (pathRoot == scNull()) pathRoot = Cstar;
  if (iterEpoch == 0) {
    if (mDismSeenEpoch.size() < mN) {
      mDismSeenEpoch.resize(mN, 0); mDismSeenOwner.resize(mN, 0);
    }
    if (++mDismEpochClock == 0) {
      std::fill(mDismSeenEpoch.begin(), mDismSeenEpoch.end(), 0);
      ++mDismEpochClock;
    }
    iterEpoch = mDismEpochClock; callOwner = 1;
  }
#ifdef HYBPROF
  ++gHybShellSearches;
#endif
  Size Ccur = scNull(), Dcur = scNull();
  if (!resolveOldShellBounds(Cstar, Dstar, pathRoot, &Ccur, &Dcur))
    return false;
  bool augmented = false;
  bool terminated = false;
  for (Size transition = 0; transition <= mNumNodes + 2; ++transition) {
    unsigned long long shellMaskStart = scalingMonotonicNanoseconds();
    // Clear only vertices touched by the previous boundary configuration.
    for (Size i = 0; i < mShellMarked.size(); ++i) {
      Size v = mShellMarked[i];
      if (v < mN) { mInShell[v] = 0; mInDstar[v] = 0; }
    }
    mShellMarked.clear();
    if (Ccur == scNull() || Ccur < mN || mZ[Ccur] == 0) {
      terminated = true; break;
    }

    std::vector<Size> cLeaves, dLeaves;
    collectNodeLeaves(Ccur, &cLeaves);
    if (Dcur != scNull() && Dcur >= mN && mZ[Dcur] != 0)
      collectNodeLeaves(Dcur, &dLeaves);
    if (Dcur != scNull() && Dcur >= mN && mZ[Dcur] != 0)
      for (Size i = 0; i < dLeaves.size(); ++i) {
        Size v = dLeaves[i];
        if (v < mN) { mInDstar[v] = 1; mShellMarked.push_back(v); }
      }
    bool absorbed = false;
    for (Size i = 0; i < cLeaves.size(); ++i) {
      Size v = cLeaves[i];
      if (v >= mN || mInDstar[v]) continue;
      mInShell[v] = 1; mShellMarked.push_back(v);
      if (v < mDismSeenEpoch.size() && mDismSeenEpoch[v] == iterEpoch &&
          mDismSeenOwner[v] != callOwner) absorbed = true;
    }
    // Figure 7 halting condition 2: the expanded shell absorbed vertices
    // searched by an earlier ShellSearch in this DismantlePath iteration.
    if (absorbed) {
      ++mCounters.absorptionHalts;
#ifdef HYBPROF
      ++gHybShellAbsorbs;
#endif
      terminated = true; break;
    }
    for (Size i = 0; i < mShellMarked.size(); ++i) {
      Size v = mShellMarked[i];
      if (v < mN && mInShell[v]) {
        mDismSeenEpoch[v] = iterEpoch; mDismSeenOwner[v] = callOwner;
      }
    }
    ++mCounters.shellBoundaryConfigurations;
    mCounters.shellMaskVertices += mShellMarked.size();
    mCounters.shellMaskNanoseconds += scalingMonotonicNanoseconds() - shellMaskStart;

    mShellActive = true; mShellCstar = Ccur; mShellDstar = Dcur;
#ifdef HYBPROF
    if (std::getenv("MBX_TRACE_DPATH")) {
      std::fprintf(stderr, "[SHELL-BEFORE] C=%lu D=%ld zC=%lld",
        (unsigned long)Ccur, (long)(Dcur == scNull() ? -1 : Dcur),
        (long long)mZ[Ccur]);
      for (Size uu = 0; uu < mN; ++uu) if (mInShell[uu])
        for (Size kk = 0; kk < mAdj[uu].size(); ++kk) {
          Size vv = mAdj[uu][kk];
          if (uu < vv && vv < mN && mInShell[vv] && slackZ(uu, vv) <= 2)
            std::fprintf(stderr, " e%lu-%lu:s%lld",
              (unsigned long)uu, (unsigned long)vv, (long long)slackZ(uu, vv));
        }
      std::fprintf(stderr, "\n");
    }
#endif
    LL cap = mZ[Ccur] / 2 + (LL)nodeLeafCount(Ccur) + 2;
    SearchRequest request;
    request.criterion = CRIT1; request.region = REGION_SHELL;
    request.maximalAugment = true; request.maxAdjustments = cap;
    request.stoppingRule = STOP_AT_SHELL_EVENT;
    request.shellTranslation = true;
    unsigned long long shellStart = scalingMonotonicNanoseconds();
    ++mCounters.shellSearchCalls;
    SearchResult shellResult = executeSearch(
      request, mEventSearch ? BACKEND_EVENT : BACKEND_SCAN);
    const unsigned long long shellElapsed =
      scalingMonotonicNanoseconds() - shellStart;
    mCounters.shellSearchNanoseconds += shellElapsed;
    mCounters.shellMcmNanoseconds += shellElapsed;
    bool aug = shellResult.augmented;
    mShellActive = false; mShellCstar = scNull(); mShellDstar = scNull();
    if (shellResult.status != STATUS_SUCCESS) { terminated = true; break; }
#ifdef HYBPROF
    if (std::getenv("MBX_TRACE_DPATH")) {
      std::fprintf(stderr, "[SHELL-AFTER] C=%lu zC=%lld",
        (unsigned long)Ccur, (long long)mZ[Ccur]);
      for (Size uu = 0; uu < mN; ++uu) if (mInShell[uu])
        for (Size kk = 0; kk < mAdj[uu].size(); ++kk) {
          Size vv = mAdj[uu][kk];
          if (uu < vv && vv < mN && mInShell[vv] && slackZ(uu, vv) <= 2)
            std::fprintf(stderr, " e%lu-%lu:s%lld",
              (unsigned long)uu, (unsigned long)vv, (long long)slackZ(uu, vv));
        }
      std::fprintf(stderr, "\n");
      if (mN > 6)
        std::fprintf(stderr,
          "[SHELL-AFTER-26] slack=%lld outer=%lu/%lu uf=%lu/%lu parents=%ld/%ld\n",
          (long long)slackZ(2, 6), (unsigned long)outerOf(2),
          (unsigned long)outerOf(6), (unsigned long)ufFind(2),
          (unsigned long)ufFind(6),
          (long)(mBParent[2] == scNull() ? -1 : mBParent[2]),
          (long)(mBParent[6] == scNull() ? -1 : mBParent[6]));
    }
#endif
    if (aug) {
      ++mCounters.shellAugmentationHalts;
      augmented = true; terminated = true; break;
    } // Figure 7 condition 1; M is an MCM

    // Figure 7 dissolves zero-z current roots before changing the inherited
    // shell boundaries or testing an absorption halt.
    dissolveZeroCurrentRoots();
    if (mLastStatus != STATUS_SUCCESS) { terminated = true; break; }

    bool cZero = mZ[Ccur] == 0;
    bool dZero = Dcur != scNull() && Dcur >= mN && mZ[Dcur] == 0;
    if (!cZero && !dZero) {
      if (mHybridCanonical) mLastStatus = STATUS_UNEXPECTED_NO_PROGRESS;
      break;
    }

    // Compute successor boundaries before dissolve() changes child-parent links.
    Size nextC = Ccur;
    if (cZero) {
      ++mCounters.boundaryChanges;
#ifdef HYBPROF
      ++gHybShellCUpdates;
#endif
      nextC = (Ccur == pathRoot) ? scNull()
                                 : nextLiveOldParent(Ccur, pathRoot);
    }
    Size nextD = Dcur;
    if (dZero) {
      ++mCounters.boundaryChanges;
#ifdef HYBPROF
      ++gHybShellDUpdates;
#endif
      nextD = nextLiveOldMajorChild(Dcur);
    }
    // Keep the zero-z old container structurally lazy here. The established engine
    // represents removal from Omega0 by z=0 while retaining its membership links;
    // eagerly calling dissolve() at this point invalidates current-scale blossom/tree
    // relationships (the dynamic A/B gate exposed seven aborts). Boundary selection
    // below skips zero-z old containers, so C*/D* still update semantically.

    // Figure 7 condition 3: the outermost old blossom on P(R) dissolved.
    if (cZero && nextC == scNull()) {
      ++mCounters.shellOutermostDissolutionHalts;
      terminated = true; break;
    }
    if (cZero) Ccur = nextC;
    if (dZero) Dcur = nextD;
    if (Ccur == Dcur) { terminated = true; break; }
  }
  if (!terminated && mLastStatus == STATUS_SUCCESS)
    mLastStatus = STATUS_UNEXPECTED_NO_PROGRESS;
#ifdef HYBPROF
  if (augmented) ++gHybShellAugCalls; else ++gHybShellNoAugCalls;
#endif
  for (Size i = 0; i < mShellMarked.size(); ++i) {
    Size v = mShellMarked[i];
    if (v < mN) { mInShell[v] = 0; mInDstar[v] = 0; }
  }
  mShellMarked.clear();
  mShellActive = false; mShellCstar = scNull(); mShellDstar = scNull();
  return augmented;
}

// #free vertices in shell G(Cstar,Dstar) = unmatched leaves of C* not in D*.
inline Size ScalingMatching::shellFreeCount(Size Cstar, Size Dstar) const {
  if (Cstar == scNull() || Cstar < mN) return 0;
  std::vector<Size> cLeaves, dLeaves;
  collectNodeLeaves(Cstar, &cLeaves);
  if (Dstar != scNull() && Dstar >= mN) collectNodeLeaves(Dstar, &dLeaves);
  // Dstar is the immediately inner inherited boundary on one laminar major
  // path, hence leaves(Dstar) is a subset of leaves(Cstar). The free vertices
  // of the atomic shell Cstar\Dstar are therefore counted exactly by
  // free(Cstar)-free(Dstar). The former implementation allocated and cleared
  // an n-byte membership bitmap for every count, even for a three-vertex
  // shell, which dominated real DismantlePath workloads with thousands of
  // distinct boundary configurations.
  Size freeC = 0, freeD = 0;
  for (Size i = 0; i < cLeaves.size(); ++i) {
    Size v = cLeaves[i];
    if (v < mN && mMate[v] == scNull()) ++freeC;
  }
  for (Size i = 0; i < dLeaves.size(); ++i) {
    Size v = dLeaves[i];
    if (v < mN && mMate[v] == scNull()) ++freeD;
  }
  return freeC >= freeD ? freeC - freeD : 0;
}

// §4.2 DismantlePath(R): STAGE 1 (atomic-shell ShellSearch loop while |F|>=2), then STAGE 2.
inline bool ScalingMatching::dismantlePath(Size R) {
  ScalingNanosecondGuard dismantleTimer(&mCounters.dismantleNanoseconds);
  ++mCounters.dismantleCalls;
#ifdef HYBPROF
  ++gHybDPathEntries;
#endif
  // STAGE 1: while P(R) has undissolved blossoms and total free in P(R) >= 2.
  Size guard = 0;
  for (;;) {
    ++mCounters.dismantleIterations;
    if (++guard > mNumNodes + 8) {
#ifdef HYBPROF
      ++gHybStage1Stalls;
#endif
      if (mHybridCanonical) mLastStatus = STATUS_UNEXPECTED_NO_PROGRESS;
      return false;                                 // unexplained nontermination is not stage 2
    }
    // current undissolved chain on P(R): innermost..outermost (B1 subset ... subset Bk).
    std::vector<Size> chain;
    unsigned long long chainStart = scalingMonotonicNanoseconds();
    undissolvedPathChain(R, &chain);               // outermost-first (R down)
    mCounters.dismantleChainNanoseconds +=
      scalingMonotonicNanoseconds() - chainStart;
    if (chain.empty()) {
#ifdef HYBPROF
      ++gHybDPathTopEmpty;
#endif
      ++mCounters.dismantleStage1AllDissolvedExits;
      return true;                                  // all dissolved
    }
    Size freeTot = 0;
    // Build atomic shells from the chain. undissolvedPathChain gives OUTERMOST-first; reverse
    // to innermost-first B1..Bk so shells are G(B_{i}, B_{i-1}).
    std::vector<Size> inner(chain.rbegin(), chain.rend());  // B1 (innermost) .. Bk (outermost)
    // shells: G(B1, emptyset) [odd], G(B2,B1), G(B3,B2), ... G(Bk, B_{k-1}).
    struct Shell { Size C, D, freeCnt; };
    std::vector<Shell> shells;
    unsigned long long freeCountStart = scalingMonotonicNanoseconds();
    {
      Size D = scNull();
      for (Size i = 0; i < inner.size(); ++i) {
        Size C = inner[i];
        Size fc = shellFreeCount(C, D);
        Shell s; s.C = C; s.D = D; s.freeCnt = fc; shells.push_back(s);
        D = C;
      }
    }
    // total free across shells = the |F| of P(R).
    freeTot = 0; for (Size i = 0; i < shells.size(); ++i) freeTot += shells[i].freeCnt;
    mCounters.dismantleFreeCountNanoseconds +=
      scalingMonotonicNanoseconds() - freeCountStart;
    if (freeTot < 2) {
      ++mCounters.dismantleStage1SingletonExits;
      break;                                        // -> STAGE 2
    }
    // sort shells NON-INCREASING by freeCnt; drop those with <2 free.
    for (Size a = 0; a < shells.size(); ++a)
      for (Size b = a + 1; b < shells.size(); ++b)
        if (shells[b].freeCnt > shells[a].freeCnt) { Shell t = shells[a]; shells[a] = shells[b]; shells[b] = t; }
    Size iterEpoch = 0;
    if (mDynamicShell) {
      if (mDismSeenEpoch.size() < mN) {
        mDismSeenEpoch.resize(mN, 0); mDismSeenOwner.resize(mN, 0);
      }
      if (++mDismEpochClock == 0) {
        std::fill(mDismSeenEpoch.begin(), mDismSeenEpoch.end(), 0);
        ++mDismEpochClock;
      }
      iterEpoch = mDismEpochClock;
    }
    Size chainBefore = chain.size(), freeBefore = freeTot;
    bool anyRun = false;
    unsigned long long shellCallsStart = scalingMonotonicNanoseconds();
    for (Size i = 0; i < shells.size(); ++i) {
      if (shells[i].freeCnt < 2) continue;
      shellSearch(shells[i].C, shells[i].D, R, iterEpoch, i + 1);
      if (mLastStatus != STATUS_SUCCESS) return false;
      anyRun = true;
    }
    mCounters.dismantleShellNanoseconds +=
      scalingMonotonicNanoseconds() - shellCallsStart;
    if (!anyRun) {
#ifdef HYBPROF
      ++gHybStage1Stalls;
#endif
      if (mHybridCanonical) mLastStatus = STATUS_INVARIANT_VIOLATION;
      return false;                                 // |F|>=2 but no legal shell: precondition failure
    }
    // PROGRESS GUARD: a full shell pass must dissolve >=1 blossom OR reduce free count;
    // else we'd loop forever on an unchanged state. (The shellSearch cap guarantees each
    // call returns; this guarantees the OUTER loop terminates.) Recompute post-pass.
    unsigned long long progressStart = scalingMonotonicNanoseconds();
    std::vector<Size> chainAfter; undissolvedPathChain(R, &chainAfter);
    Size freeAfter = 0;
    { Size D = scNull(); std::vector<Size> ia(chainAfter.rbegin(), chainAfter.rend());
      for (Size i = 0; i < ia.size(); ++i) { freeAfter += shellFreeCount(ia[i], D); D = ia[i]; } }
    mCounters.dismantleProgressNanoseconds +=
      scalingMonotonicNanoseconds() - progressStart;
#ifdef HYBPROF
    if (chainAfter.empty()) ++gHybPostShellEmpty;
    else if (freeAfter == 1) ++gHybPostShellFree1;
    else ++gHybPostShellOther;
#endif
    if (chainAfter.size() >= chainBefore && freeAfter >= freeBefore) {
#ifdef HYBPROF
      ++gHybStage1Stalls;
#endif
      if (mHybridCanonical) mLastStatus = STATUS_UNEXPECTED_NO_PROGRESS;
      return false;                                  // do not misroute a stalled stage 1 into stage 2
    }
  }
  // STAGE 2: |F|=1 cleanup (liquidate remaining + bounded PQSearch from omega).
  return dismantleStage2(R);
}

// gabowDissolveSmallBlossom: §4.2 entry -- dissolve all old sub-blossoms in B's subtree by
// calling DismantlePath on each major-path root in POSTORDER (descendants first, so when
// DismantlePath(R) runs the only undissolved blossoms left in R are those on P(R)).
inline bool ScalingMatching::gabowDissolveSmallBlossom(Size B) {
#ifdef HYBPROF
  ++gHybDissolveCalls;
  // budget for Lemma 4.10: m(B) * n(B)^{3/4}. m(B) = #edges with >=1 endpoint in B (approx via
  // sum of degrees of B's leaves); n(B) = #leaves.
  { double nB = (double)nodeLeafCount(B); double mB = 0;
    std::vector<Size> budgetLeaves; collectNodeLeaves(B, &budgetLeaves);
    for (Size i = 0; i < budgetLeaves.size(); ++i) { Size v = budgetLeaves[i]; if (v < mN) mB += (double)mAdj[v].size(); }
    gHybBudget += mB * pow(nB, 0.75);
    gHybInDissolve = true; }
#endif
  std::vector<Size> roots;
  collectMajorPathRoots(B, &roots);     // postorder
  // B is the root of THIS dissolution subtree even when its immutable old
  // parent is a liquidated large blossom.  isMajorPathRoot() sees the full old
  // forest and may therefore omit B when it was that parent's major child.
  if (B < mZ.size() && mZ[B] != 0 &&
      (roots.empty() || roots.back() != B)) roots.push_back(B);
#ifdef HYBPROF
  if (std::getenv("MBX_TRACE_DPATH")) {
    std::fprintf(stderr, "[GABOW] scale=%lu B=%lu z=%lld leaves=",
      (unsigned long)mAuditScale, (unsigned long)B, (long long)mZ[B]);
    std::vector<Size> traceLeaves; collectNodeLeaves(B, &traceLeaves);
    for (Size i = 0; i < traceLeaves.size(); ++i)
      std::fprintf(stderr, " %lu", (unsigned long)traceLeaves[i]);
    std::fprintf(stderr, " roots=");
    for (Size i = 0; i < roots.size(); ++i)
      std::fprintf(stderr, " %lu(z=%lld)", (unsigned long)roots[i],
                   (long long)mZ[roots[i]]);
    std::fprintf(stderr, "\n");
  }
#endif
  for (Size i = 0; i < roots.size(); ++i)
    if (!dismantlePath(roots[i])) {
#ifdef HYBPROF
      gHybInDissolve = false;
#endif
      return false;
    }
#ifdef HYBPROF
  gHybInDissolve = false;
#endif
  return true;
}

// Strip dummy edges (endpoint >= mN0) from each real vertex IN PLACE; clear dummy adjacency.
inline void ScalingMatching::stripDummyAdj(Size mN0) {
  for (Size u = 0; u < mN0; ++u) {
    std::vector<Size>& a = mAdj[u]; std::vector<LL>& wv = mW[u]; std::vector<LL>& hv = mWhat[u]; std::vector<LL>& bv = mWbar[u];
    Size j = 0;
    for (Size i = 0; i < a.size(); ++i)
      if (a[i] < mN0) {                       // keep real-real edge; compact toward front
        if (j != i) { a[j] = a[i]; wv[j] = wv[i]; hv[j] = hv[i]; bv[j] = bv[i]; }
        ++j;
      }
    a.resize(j); wv.resize(j); hv.resize(j); bv.resize(j);   // drop the stripped tail (no realloc)
  }
  for (Size d = mN0; d < 2 * mN0; ++d) { mAdj[d].clear(); mW[d].clear(); mWhat[d].clear(); }
}

// Liquidate(B): y(u) += z(B)/2 for all u in B; z(B)=0; dissolve B (free its subs).
inline void ScalingMatching::liquidate(Size B) {
  CRUMB("liquidate", B, 0);
  if (B < mN) return;
  LL half = mZ[B] / 2;
  std::vector<Size> verts;
  collectBlossomVerts(B, &verts);
  for (Size i = 0; i < verts.size(); ++i) mY[verts[i]] += half;
  mZ[B] = 0;
  dissolve(B);   // sub-nodes become outer (mBParent -> scNull)
}

// Split the two blossom families required by §4.2. The inherited forest Omega_0
// remains available through mOldParent/mOldContaining for major paths and yz,
// while the ordinary mBParent/UF forest is reset to bare vertices so Edmonds
// search can build an independent current-scale Omega.
inline void ScalingMatching::separateInheritedForest(Size oldBegin, Size oldEnd) {
#ifdef EVENT_SFM_REUSE
  // Omega changes identity at the scale boundary: the preceding current
  // forest becomes immutable Omega', so its lists cannot be reused.
  mEventSfmInitialized = false;
#endif
#ifdef PS5_LEAF_SPLICE
  scalingLfSnapshotOld(oldBegin, oldEnd);
#endif
  mOldForestSeparated = true;
  mShellLocalResetReady = false;
  mOldParent = mBParent;
  mOldContaining.assign(mN, std::vector<Size>());
  for (Size v = 0; v < mN; ++v) {
    Size p = v < mOldParent.size() ? mOldParent[v] : scNull();
    Size guard = 0;
    while (p != scNull() && p < oldEnd && ++guard <= oldEnd + 2) {
      if (p >= oldBegin) mOldContaining[v].push_back(p);
      p = p < mOldParent.size() ? mOldParent[p] : scNull();
    }
  }
  // Detach only inherited parent links. Old nodes/cycles remain immutable
  // records; new shrink operations allocate ids >= oldEnd in mBParent.
  for (Size i = 0; i < oldEnd && i < mBParent.size(); ++i)
    mBParent[i] = scNull();
  // Scale separation detaches inherited children without calling dissolve().
  // Wake every still-dormant label once so the next fresh search resets the
  // now-outer records. This is O(number of dormant records) once per scale,
  // rather than once per search.
  for (Size i = 0; i < mEventDormantNodes.size(); ++i) {
    Size node = mEventDormantNodes[i];
    if (node >= oldEnd || node >= mEventDormantFlag.size() ||
        !mEventDormantFlag[node]) continue;
    rememberCompactEventLabel(node);
    ++mCounters.eventScaleResetCandidates;
  }
  mEventDormantNodes.clear();
  if (!mUf.empty())
    for (Size v = 0; v < mN && v < mUf.size(); ++v) mUf[v] = v;
  if (!mBaseToNode.empty())
    for (Size v = 0; v < mN && v < mBaseToNode.size(); ++v) mBaseToNode[v] = v;

}

inline void ScalingMatching::finishInheritedForest() {
#ifdef EVENT_SFM_REUSE
#ifndef EVENT_SFM_BOUNDARY_FRONTIER_FAST
  // The within-scale transition from separated Omega'|Omega to Omega alone
  // is a second list-ownership epoch. Rebuild once at the next search.
  mEventSfmInitialized = false;
#else
  // The fast frontier never inserted inherited Omega' into split-findmin.
  // Current Omega merges/splits were maintained throughout dissolution, so
  // removing the immutable inherited forest does not change any represented
  // list. Preserve that topology until the next scale boundary.
#endif
#endif
#ifdef HYBPROF
  if (std::getenv("MBX_TRACE_STAGE2_AUG")) {
    Size live = 0;
    for (Size b = mAuditOldBegin; b < mAuditOldEnd && b < mZ.size(); ++b)
      if (mZ[b] != 0) ++live;
    if (live)
      std::fprintf(stderr, "[OLD-FINISH] scale=%lu inherited-z-positive=%lu\n",
        (unsigned long)mAuditScale, (unsigned long)live);
  }
#endif
  // Omega' records are immutable history after dissolution.  Their cycles are
  // deliberately retained for diagnostics, but they must not participate in
  // any live-current-blossom structure (notably the Section-5 split-findmin
  // initializer) once the inherited phase ends.  Current-scale blossoms were
  // allocated at ids >= mAuditOldEnd and are unaffected.
  for (Size b = mAuditOldBegin; b < mAuditOldEnd && b < mBActive.size(); ++b)
    mBActive[b] = 0;
  mOldForestSeparated = false;
  mShellLocalResetReady = false;
  mOldParent.clear();
  mOldContaining.clear();
}

inline bool ScalingMatching::checkCanonicalPostDissolution(Size realN) const {
  for (Size u = 0; u < realN; ++u)
    for (Size k = 0; k < mAdj[u].size(); ++k) {
      Size v = mAdj[u][k];
      if (v < realN && u < v && slackZ(u, v) < 0) {
#ifdef HYBPROF
        std::fprintf(stderr,
          "[HYBRID-INVARIANT] scale=%lu post-dissolve domination edge=(%lu,%lu) slack=%lld w=%lld y=%lld/%lld oldSep=%d oldCommon=%lld uf=%lu/%lu outer=%lu/%lu\n",
          (unsigned long)mAuditScale, (unsigned long)u, (unsigned long)v,
          (long long)slackZ(u, v), (long long)mW[u][k],
          (long long)mY[u], (long long)mY[v], mOldForestSeparated ? 1 : 0,
          (long long)oldCommonZ(u, v), (unsigned long)ufFind(u),
          (unsigned long)ufFind(v), (unsigned long)outerOf(u),
          (unsigned long)outerOf(v));
#endif
        return false;
      }
    }
  int parity = -1;
  for (Size v = 0; v < realN; ++v) if (mMate[v] == scNull()) {
    int p = (int)(mY[v] & 1);
    if (parity < 0) parity = p;
    else if (p != parity) return false;
  }
  return true;
}

inline bool ScalingMatching::canonicalDualObjective(Size realN, LL* objective) const {
  if (!objective) return false;
  LL total = 0;
  for (Size v = 0; v < realN; ++v) {
    LL next = 0;
    if (!checkedAddLL(total, mY[v], &next)) return false;
    total = next;
  }
  for (Size b = mAuditOldBegin; b < mNumNodes && b < mZ.size(); ++b) {
    if (mZ[b] == 0 || b >= mNodeLeaves.size()) continue;
    LL term = 0, next = 0;
    if (!checkedMulLL(mZ[b], (LL)(nodeLeafCount(b) / 2), &term) ||
        !checkedAddLL(total, term, &next)) return false;
    total = next;
  }
  *objective = total;
  return true;
}

inline bool ScalingMatching::checkCanonicalSeparatedForests(Size realN) const {
  if (!mOldForestSeparated) return false;
  // Omega' is immutable and detached from the mutable current parent forest.
  for (Size b = mAuditOldBegin; b < mAuditOldEnd; ++b) {
    if (b >= mOldParent.size() || b >= mBParent.size() || mBParent[b] != scNull())
      return false;
    Size p = mOldParent[b], guard = 0;
    while (p != scNull()) {
      if (p < mAuditOldBegin || p >= mAuditOldEnd || ++guard > mAuditOldEnd + 2)
        return false;
      p = p < mOldParent.size() ? mOldParent[p] : scNull();
    }
  }

  // No current blossom may contain an inherited blossom node.
  for (Size b = mAuditOldEnd; b < mNumNodes; ++b)
    if (b < mBCycle.size())
      for (Size i = 0; i < mBCycle[b].size(); ++i) {
        Size child = mBCycle[b][i];
        if (child >= mAuditOldBegin && child < mAuditOldEnd) return false;
      }

  std::vector<char> in(realN, 0);
  // No matched edge crosses a live inherited blossom.
  for (Size b = mAuditOldBegin; b < mAuditOldEnd; ++b) {
    if (b >= mZ.size() || mZ[b] == 0 || b >= mNodeLeaves.size()) continue;
    std::fill(in.begin(), in.end(), (char)0);
    std::vector<Size> leaves; collectNodeLeaves(b, &leaves);
    for (Size i = 0; i < leaves.size(); ++i)
      if (leaves[i] < realN) in[leaves[i]] = 1;
    for (Size u = 0; u < realN; ++u) if (mMate[u] != scNull() && mMate[u] < realN)
      if (in[u] != in[mMate[u]]) return false;
  }

  // Every active positive-z current blossom is full. Every active current root
  // surviving zero-root cleanup has positive z.
  for (Size b = mAuditOldEnd; b < mNumNodes; ++b) {
    if (b >= mBActive.size() || !mBActive[b]) continue;
    if (mBParent[b] == scNull() && (b >= mZ.size() || mZ[b] <= 0)) return false;
    if (b >= mZ.size() || mZ[b] <= 0 || b >= mNodeLeaves.size()) continue;
    std::fill(in.begin(), in.end(), (char)0);
    std::vector<Size> leaves; collectNodeLeaves(b, &leaves);
    for (Size i = 0; i < leaves.size(); ++i)
      if (leaves[i] < realN) in[leaves[i]] = 1;
    Size internalMatchedTwice = 0;
    for (Size u = 0; u < realN; ++u)
      if (in[u] && mMate[u] != scNull() && mMate[u] < realN && in[mMate[u]])
        ++internalMatchedTwice;
    if (internalMatchedTwice / 2 != nodeLeafCount(b) / 2) return false;
  }
  return true;
}

// Liquidation of an Omega_0 record must not call dissolve(): its children may
// now belong to current-scale Omega blossoms. Only the old dual is transferred.
inline void ScalingMatching::liquidateInherited(Size B) {
  if (B < mN || B >= mZ.size()) return;
  LL half = mZ[B] / 2;
  if (B < mNodeLeaves.size()) {
    std::vector<Size> leaves; collectNodeLeaves(B, &leaves);
    for (Size i = 0; i < leaves.size(); ++i) {
      Size v = leaves[i];
      if (v < mY.size()) mY[v] += half;
    }
  }
  mZ[B] = 0;
}

// Figure 5, step 6.  The member sets were captured before Omega' was
// liquidated. Each PQSearch declares one maximal former small blossom as an
// induced vertex-set region. executeSearch uses the same proven mask machinery
// internally, without Figure 7 boundary translation.
inline bool ScalingMatching::runLiquidationistLeveling(
    const std::vector<std::vector<Size> >& maximalOldSmall, Size realN) {
  for (Size bi = 0; bi < maximalOldSmall.size(); ++bi) {
    const std::vector<Size>& members = maximalOldSmall[bi];
    if (members.empty()) {
      mLastStatus = STATUS_INVARIANT_VIOLATION;
      return false;
    }
    ++mCounters.liquidationistStep6Regions;

    while (true) {
      LL top = 0, next = 0;
      Size freeBefore = 0;
      for (Size i = 0; i < members.size(); ++i) {
        Size v = members[i];
        if (v >= realN || mMate[v] != scNull()) continue;
        ++freeBefore;
        if (mY[v] > top) top = mY[v];
      }
      if (top <= 0) break;

      std::vector<Size> roots;
      for (Size i = 0; i < members.size(); ++i) {
        Size v = members[i];
        if (v >= realN || mMate[v] != scNull()) continue;
        if (mY[v] == top) roots.push_back(v);
        else if (mY[v] > next) next = mY[v];
      }
      if (roots.empty() || next < 0 || next >= top) {
        mLastStatus = STATUS_INVARIANT_VIOLATION;
        return false;
      }

      SearchRequest request;
      request.roots = &roots;
      request.regionVertices = &members;
      request.criterion = CRIT1;
      request.region = REGION_VERTEX_SET;
      request.dualBudget = top - next;
      request.maxAdjustments = -1;
      request.batchEventFreeTime = true;       // PQSearch may jump to the next event
      request.stoppingRule = STOP_AT_DUAL_BUDGET;
      SearchResult result = executeSearch(request, BACKEND_SCAN);
      ++mCounters.liquidationistStep6Searches;
      if (result.augmented) ++mCounters.liquidationistStep6Augmentations;
      mCounters.liquidationistStep6DualUnits +=
          (unsigned long long)(result.dualUnits >= 0 ? result.dualUnits : 0);
      if (result.status != STATUS_SUCCESS) return false;

      LL newTop = 0;
      Size freeAfter = 0;
      for (Size i = 0; i < members.size(); ++i) {
        Size v = members[i];
        if (v >= realN || mMate[v] != scNull()) continue;
        ++freeAfter;
        if (mY[v] > newTop) newTop = mY[v];
      }
      // The paper gives the progress proof: a call either augments or spends
      // exactly top-next dual units.  Treat any violation as a failed invariant,
      // never as a behavior-changing progress cap.
      if ((!result.augmented && result.dualUnits != top - next) ||
          (freeAfter >= freeBefore && newTop >= top)) {
        mLastStatus = STATUS_UNEXPECTED_NO_PROGRESS;
        return false;
      }
    }
  }
  return true;
}

// Carried scaling loop shared by Hybrid and Liquidationist
// branch. Fills mMate over the REAL vertices 0..mN0-1. Liquidationist
// materializes step-4 reduced costs; Hybrid omits that step as Figure 6 requires.
inline bool ScalingMatching::runScaling() {
  Size mN0 = mN;                       // real vertices
  if (mN0 == 0) { mNoPerfect = false; mLastStatus = STATUS_SUCCESS; return true; }
  if ((mN0 & 1) != 0 && !mDropPerfectExperimental) {
    mNoPerfect = true; mLastStatus = STATUS_NO_PERFECT_MATCHING; return false;
  }

  // Working node set can grow: blossoms (>=mN0 indexing into mB* arrays) and
  // DUMMY vertices appended at the END of the vertex range. To keep vertices
  // and blossoms cleanly separated we reserve dummy vertex ids in [mN0, 2*mN0).
  // So blossom ids start at 2*mN0. We therefore treat "mN" (the vertex/blossom
  // boundary used everywhere as `>= mN => blossom`) as 2*mN0 during scaling.
  Size dummyBase = mN0;                // dummy d(u) = dummyBase + u
  Size blossomBase = 2 * mN0;          // blossom ids start here
  mN = blossomBase;                    // vertices are [0, 2*mN0): real + dummy slots

  // adjacency for dummies: each real u has a 0-weight edge to its dummy d(u).
  // We APPEND dummy adjacency lazily in Perfection; pre-size arrays now.
  Size maxNodes = 4 * mN0 + 4;
  mAdj.resize(2 * mN0);
  mW.resize(2 * mN0);
  mWhat.resize(2 * mN0);
  for (Size d = mN0; d < 2 * mN0; ++d) { mAdj[d].clear(); mW[d].clear(); mWhat[d].clear(); }

  mY.assign(2 * mN0, 0);
  mZ.assign(maxNodes, 0);
  mMate.assign(2 * mN0, scNull());
  mBParent.assign(maxNodes, scNull());
  mBBase.assign(maxNodes, scNull());
  mBCycle.assign(maxNodes, std::vector<Size>());
  mBActive.assign(maxNodes, 0);
  mBCycEdgeU.assign(maxNodes, std::vector<Size>());
  mBCycEdgeV.assign(maxNodes, std::vector<Size>());
  mNodeLeaves.assign(maxNodes, std::vector<Size>());
  mNLabel.assign(maxNodes, LBL_FREE_OF_TREE);
  mNParent.assign(maxNodes, scNull());
  mNParEdgeU.assign(maxNodes, scNull());
  mNParEdgeV.assign(maxNodes, scNull());
  mNRoot.assign(maxNodes, scNull());
  mUf.assign(maxNodes, 0);             // union-find: each vertex its own rep initially
  mBEverPos.assign(maxNodes, 0);       // no blossom has been z-positive yet
  mBaseToNode.assign(maxNodes, scNull());
  for (Size v = 0; v < 2 * mN0; ++v) { mBBase[v] = v; mNodeLeaves[v].assign(1, v); mUf[v] = v; mBaseToNode[v] = v; }
#ifdef PS5_LEAF_SPLICE
  scalingLfInit(2 * mN0, maxNodes);
#endif
  mNumNodes = blossomBase;             // active node ids; blossoms appended >= here
  mUseScaleW = true;

  // current scale even weight per directed edge slot, starts at 0.
  for (Size v = 0; v < mN0; ++v) { mW[v].assign(mAdj[v].size(), 0); }
  // faithful-reweight: accumulated per-edge offset to recover the true weight at the end.
  if (mFaithfulReweight) { mWoff.assign(2 * mN0, std::vector<LL>()); for (Size v = 0; v < mN0; ++v) mWoff[v].assign(mAdj[v].size(), 0); }

  // w'(prev scale) per slot; y',z' carried via mY/mZ between scales (in place).
  std::vector<std::vector<LL> > wPrev(mN0);
  for (Size v = 0; v < mN0; ++v) wPrev[v].assign(mAdj[v].size(), 0);

  Size nbits = mNumScales;             // inclusive extended-weight bit length
#ifdef SCPHASE
  std::fprintf(stderr, "[SCPHASE] scale-loop START clock=%.3fs nbits=%lu\n", (double)clock()/CLOCKS_PER_SEC, (unsigned long)nbits);
  double phLast = (double)clock()/CLOCKS_PER_SEC;
#endif
  for (Size scale = 1; scale <= nbits; ++scale) {
    // Blossom ids are monotone. Before this scale performs any search, every
    // allocated nontrivial node is inherited from the preceding scale.
    mAuditScale = scale;
    mAuditOldBegin = blossomBase;
    mAuditOldEnd = mNumNodes;
#ifdef SCMAG
    { LL ys=0, zs=0; for (Size v=0; v<mN0; ++v) ys+=mY[v]; for (Size b=2*mN0; b<mNumNodes; ++b) zs+=mZ[b];
      Size matched=0; for (Size v=0; v<mN0; ++v) if (mMate[v]!=scNull()) ++matched;
      std::fprintf(stderr,"[CARRY-IN scale %lu] sumY=%lld sumZ=%lld matchedReals=%lu\n",
        (unsigned long)scale,(long long)ys,(long long)zs,(unsigned long)matched);
      if (std::getenv("MBX_DUMPSCALE") && scale == (Size)atoi(std::getenv("MBX_DUMPSCALE"))) {
        std::fprintf(stderr,"  [DUMP scale %lu] y:", (unsigned long)scale);
        for (Size v=0; v<mN0; ++v) std::fprintf(stderr," %lld", (long long)mY[v]);
        std::fprintf(stderr,"\n  [DUMP] mate:");
        for (Size v=0; v<mN0; ++v) std::fprintf(stderr," %ld", (long)(mMate[v]==scNull()?-1:(long)mMate[v]));
        std::fprintf(stderr,"\n  [DUMP] blossom z (b>=2mN0, z!=0):");
        for (Size b=2*mN0; b<mNumNodes; ++b) if (mBParent[b]==scNull() && mZ[b]!=0) std::fprintf(stderr," b%lu:z%lld", (unsigned long)b,(long long)mZ[b]);
        std::fprintf(stderr,"\n");
      } }
#endif
#ifdef SCPHASE
    { double t=(double)clock()/CLOCKS_PER_SEC;
#ifdef SCPROF
      static long long lastAdj=0,lastIt=0,lastRel=0,lastSh=0,lastDis=0,lastCalls=0; long long rl=gProfRelabels-lastRel, sh=gRelShrink-lastSh, dis=gRelDissolve-lastDis, calls=gProfCalls-lastCalls; std::fprintf(stderr,"[SCPHASE]  scale %lu start +%.3fs  (relabels=%lld [shrink=%lld dissolve=%lld first/call=%lld] adj=%lld iters=%lld calls=%lld)\n",(unsigned long)scale,t-phLast,rl,sh,dis,rl-sh-dis,gProfScaleAdj-lastAdj,gProfIters-lastIt,calls); lastAdj=gProfScaleAdj; lastIt=gProfIters; lastRel=gProfRelabels; lastSh=gRelShrink; lastDis=gRelDissolve; lastCalls=gProfCalls;
#else
      std::fprintf(stderr,"[SCPHASE]  scale %lu start +%.3fs\n",(unsigned long)scale,t-phLast);
#endif
      phLast=t; }
#ifdef OPCNT
    { static long long lE=0,lO=0,lS=0,lU=0,lM=0,lD=0,lH=0,l42=0,lFv=0,lEl=0; extern long long gDisIn42,gDisFvr; long long sc=gOpEdgeScans-lE, el=gOpEligible-lEl; std::fprintf(stderr,"        [OPCNT scale %lu] edgeScans=%lld ELIGIBLE=%lld (scan/elig=%.0f) shrinks=%lld dissolves=%lld outerSteps=%lld ufSlowPath=%lld\n",(unsigned long)scale,sc,el,el?(double)sc/el:0.0,gOpShrinks-lH,gOpDissolves-lD,gOpOuterSteps-lO,gOpUfSlowPath-lU); lE=gOpEdgeScans; lO=gOpOuterSteps; lS=gOpSlackSteps; lU=gOpUfSlowPath; lM=gOpUfMiss; lD=gOpDissolves; lH=gOpShrinks; l42=gDisIn42; lFv=gDisFvr; lEl=gOpEligible; }
#endif
#endif
    // --- Scaling step: w(e) <- 2(w'(e) + i-th bit of wbar(e)); y<-2y'+3; z<-2z' ---
    // i-th bit counting MSB first: bit index = nbits - scale.
    Size bitpos = nbits - scale;
    for (Size v = 0; v < mN0; ++v) {
      // mAdj[v] may have GROWN since wPrev was sized (dummy parks from a prior scale's
      // Perfection append edges). Keep wPrev[v] in sync (new slots = 0; a dummy edge's
      // w' is 0). Without this, wPrev[v][i] read below overflows (was a heap bug at :2175).
      if (wPrev[v].size() < mAdj[v].size()) wPrev[v].resize(mAdj[v].size(), 0);
      for (Size i = 0; i < mAdj[v].size(); ++i) {
        LL bit = (mWbar[v][i] >> bitpos) & 1;
        LL sum = 0;
        if (!checkedAddLL(wPrev[v][i], bit, &sum) ||
            !checkedMulLL(2, sum, &mW[v][i])) {
          mLastStatus = STATUS_UNSUPPORTED_OR_OVERFLOWING_INPUT;
          return false;
        }
      }
    }
    // Rescale duals of ALL vertices (real AND dummy -- a dummy parked in an
    // earlier scale must rescale too, else its 0-weight edge slack relative to
    // the inflated real y's blows up and the parked real can never be un-parked
    // == the multi-scale carry bug). z of all live blossoms doubles.
    for (Size v = 0; v < 2 * mN0; ++v) {
      LL twice = 0;
      if (!checkedMulLL(2, mY[v], &twice) ||
          !checkedAddLL(twice, 3, &mY[v])) {
        mLastStatus = STATUS_UNSUPPORTED_OR_OVERFLOWING_INPUT;
        return false;
      }
    }
    for (Size b = blossomBase; b < mNumNodes; ++b)
      if (!checkedMulLL(2, mZ[b], &mZ[b])) {
        mLastStatus = STATUS_UNSUPPORTED_OR_OVERFLOWING_INPUT;
        return false;
      }
#ifdef SCPHASE
    std::fprintf(stderr, "[SCPHASE]  scale %lu rescale-done clock=%.3fs\n",
                 (unsigned long)scale, (double)clock()/CLOCKS_PER_SEC);
    std::fflush(stderr);
#endif
#ifdef HYBPROF
    if (mHybridCanonical && std::getenv("MBX_TRACE_DPATH")) {
      for (Size u = 0; u < mN0; ++u) for (Size k = 0; k < mAdj[u].size(); ++k) {
        Size v = mAdj[u][k];
        if (u < v && v < mN0 && slackZ(u, v) < 0)
          std::fprintf(stderr,
            "[HYBRID-INVARIANT] scale=%lu pre-dissolve edge=(%lu,%lu) slack=%lld w=%lld y=%lld/%lld oldCommon=%lld\n",
            (unsigned long)scale, (unsigned long)u, (unsigned long)v,
            (long long)slackZ(u, v), (long long)mW[u][k],
            (long long)mY[u], (long long)mY[v], (long long)oldCommonZ(u, v));
      }
    }
#endif
#ifdef SCMAG
    if (scale==2){ std::fprintf(stderr,"[LOOP scale 2 POST-RESCALE pre-search] y:"); for(Size v=0;v<mN0;++v) std::fprintf(stderr," %lld",(long long)mY[v]); std::fprintf(stderr,"\n"); }
#endif

    // --- Inherited-blossom step. LIQUIDATIONIST (default): liquidate ALL inherited root
    // blossoms (z folds into member y's; domination preserved). HYBRID (mHybrid): liquidate
    // LARGE blossoms (|B|>=tau) but DISSOLVE small ones via Gabow's §4.2 O(m n^{3/4}) routine
    // (gabowDissolveSmallBlossom) -- the paper's Thm-4.3 path. Snapshot the root list first
    // (dissolution mutates mNumNodes/mBParent). ---
    if (mHybrid) {
      // PAPER-FAITHFUL ORDER (Fig 6 / Fig 5 step 1): the matching is EMPTIED (M<-emptyset)
      // at scale start, BEFORE small-blossom dissolution -- so §4.2's precondition "B
      // contains NO matched edges" holds and stage-1 ShellSearch (|F|>=2) is reachable.
      // (Our default loop reorders M<-empty to after dissolution; for Hybrid we must empty
      // FIRST.) Un-park carried reals + strip dummy adjacency (in place), then clear M.
      stripDummyAdj(mN0);
      mMate.assign(2 * mN0, scNull());      // M <- empty (spec step 1) BEFORE dissolution
      const bool splitOldForest = mHybridCanonical ||
                                  std::getenv("MBX_SEPARATE_OLD_FOREST") != 0;
      if (splitOldForest) {
        separateInheritedForest(blossomBase, mAuditOldEnd);
        // Figure 6 treats every inherited large blossom by liquidation, not
        // merely the outer roots.  The old forest is immutable in this mode,
        // so transfer all large z values first.
        for (Size b = blossomBase; b < mAuditOldEnd; ++b)
          if (b < mZ.size() && mZ[b] != 0 &&
              b < mNodeLeaves.size() && nodeLeafCount(b) >= mTau)
            liquidateInherited(b);
#ifdef SCALING_INVARIANTS
        if (!checkCanonicalSeparatedForests(mN0)) {
          finishInheritedForest();
          mLastStatus = STATUS_INVARIANT_VIOLATION;
          return false;
        }
#endif

        // What remains is a forest of small positive-z blossoms.  Dismantle
        // each maximal remaining subtree; zero-z large ancestors do not own
        // that subtree anymore.
        for (Size b = blossomBase; b < mAuditOldEnd; ++b) {
          if (b >= mZ.size() || mZ[b] == 0 || b >= mNodeLeaves.size() ||
              nodeLeafCount(b) >= mTau) continue;
          bool positiveAncestor = false;
          Size p = b < mOldParent.size() ? mOldParent[b] : scNull();
          Size guard = 0;
          while (p != scNull() && p < mAuditOldEnd && ++guard <= mAuditOldEnd + 2) {
            if (p >= blossomBase && p < mZ.size() && mZ[p] != 0) {
              positiveAncestor = true; break;
            }
            p = p < mOldParent.size() ? mOldParent[p] : scNull();
          }
          if (!positiveAncestor) {
#ifdef SCALING_INVARIANTS
            LL objectiveBefore = 0, objectiveAfter = 0;
            if (!canonicalDualObjective(mN0, &objectiveBefore)) {
              finishInheritedForest();
              mLastStatus = STATUS_UNSUPPORTED_OR_OVERFLOWING_INPUT;
              return false;
            }
#endif
            if (!gabowDissolveSmallBlossom(b)) {
              finishInheritedForest();
              if (mLastStatus == STATUS_SUCCESS)
                mLastStatus = STATUS_INVARIANT_VIOLATION;
              return false;
            }
#ifdef SCALING_INVARIANTS
            if (!canonicalDualObjective(mN0, &objectiveAfter)) {
              finishInheritedForest();
              mLastStatus = STATUS_UNSUPPORTED_OR_OVERFLOWING_INPUT;
              return false;
            }
            bool forestOk = checkCanonicalSeparatedForests(mN0);
            if (objectiveAfter > objectiveBefore || !forestOk) {
#ifdef HYBPROF
              std::fprintf(stderr,
                "[HYBRID-INVARIANT] post-Gabow B=%lu objective=%lld->%lld forestOk=%d\n",
                (unsigned long)b, (long long)objectiveBefore,
                (long long)objectiveAfter, forestOk ? 1 : 0);
#endif
              finishInheritedForest();
              mLastStatus = STATUS_INVARIANT_VIOLATION;
              return false;
            }
#endif
#ifdef HYBPROF
            if (!checkCanonicalPostDissolution(mN0))
              std::fprintf(stderr,
                "[HYBRID-INVARIANT] first visible after Gabow subtree B=%lu scale=%lu\n",
                (unsigned long)b, (unsigned long)scale);
#endif
          }
        }
        // Omega' must be empty before FVR begins.  Retained immutable records
        // may still carry their cycles/leaves, but no inherited dual may remain
        // live after all maximal small subtrees have been dismantled.
        for (Size b = blossomBase; b < mAuditOldEnd && b < mZ.size(); ++b)
          if (mZ[b] != 0) {
            finishInheritedForest();
            mLastStatus = STATUS_INVARIANT_VIOLATION;
            return false;
        }
        finishInheritedForest();
#ifdef SCALING_INVARIANTS
        if (!checkCanonicalPostDissolution(mN0)) {
          mLastStatus = STATUS_INVARIANT_VIOLATION;
          return false;
        }
#endif
      } else {
      std::vector<Size> roots;
      for (Size b = blossomBase; b < mNumNodes; ++b)
        if (mBParent[b] == scNull() &&
            mZ[b] != 0) roots.push_back(b);
      for (Size i = 0; i < roots.size(); ++i) {
        Size b = roots[i];
        if (mBParent[b] != scNull() || mZ[b] == 0) continue;   // already gone
        if (nodeLeafCount(b) >= mTau) {
          liquidate(b);                                       // LARGE -> liquidate
        }
        else {
#ifdef HYBPROF
          { Size fr = 0; std::vector<Size> leaves; collectNodeLeaves(b, &leaves); for (Size j = 0; j < leaves.size(); ++j) { Size v = leaves[j]; if (v < mN0 && mMate[v] == scNull()) ++fr; }
            if (fr > gHybMaxFreeInSmall) gHybMaxFreeInSmall = fr; }
#endif
          if (!gabowDissolveSmallBlossom(b)) {
            return false;
          }                                                    // SMALL -> Gabow §4.2
        }
      }
      }
#ifdef SCMAG
      if (std::getenv("MBX_TRACE_FVR_PARITY")) {
        std::fprintf(stderr, "[SCMAG]   scale%lu-postDissolve free-y:", (unsigned long)scale);
        for (Size v=0; v<mN0; ++v) if (mMate[v]==scNull())
          std::fprintf(stderr, " v%lu=%lld", (unsigned long)v, (long long)mY[v]);
        std::fprintf(stderr, "\n");
      }
#endif
#ifdef HYBPROF
      // POST-DISSOLUTION Property-1 check: domination yz(e)>=w(e) on real edges + worst slack +
      // whether all small old blossoms actually dissolved. domViol>0 => §4.2 left duals INVALID.
      { Size dom=0, undissolved=0; LL worst=0;
        for (Size u = 0; u < mN0; ++u)
          for (Size i = 0; i < mAdj[u].size(); ++i) {
            Size v = mAdj[u][i]; if (v >= mN0 || v <= u) continue;
            LL sk = slackZ(u, v);
            if (sk < 0) { ++dom; if (sk < worst) worst = sk; }
          }
        // Current-scale Omega blossoms are allowed (and expected) here.  Only
        // inherited Omega' ids belong to the dissolution postcondition.
        for (Size b = blossomBase; b < mAuditOldEnd && b < mZ.size(); ++b)
          if (mZ[b] != 0 && b < mNodeLeaves.size() &&
              nodeLeafCount(b) < mTau) ++undissolved;
        if (dom || undissolved) {
          std::fprintf(stderr, "[HYBPROF] post-dissolve: domViol=%lu WORST=%lld undissolvedSmall=%lu\n",(unsigned long)dom,(long long)worst,(unsigned long)undissolved);
          // dump up to 3 violating edges: endpoints' y, the edge w, the slack
          Size shown=0;
          for (Size u = 0; u < mN0 && shown < 3; ++u)
            for (Size i = 0; i < mAdj[u].size() && shown < 3; ++i) {
              Size v = mAdj[u][i]; if (v >= mN0 || v <= u) continue;
              LL sk = slackZ(u, v);
              if (sk < 0) { std::fprintf(stderr, "    viol e(%lu,%lu) w=%lld y[u]=%lld y[v]=%lld slack=%lld\n",(unsigned long)u,(unsigned long)v,(long long)mW[u][i],(long long)mY[u],(long long)mY[v],(long long)sk); ++shown; }
            }
        }
      }
#endif
    } else {
      // Figure 5 steps 1--6.  Preserve Omega' long enough to classify maximal
      // old small blossoms, then transfer every inherited dual (including
      // nested descendants).  A single root-only dissolve pass is insufficient:
      // children exposed after their parent's dissolution would be skipped.
      stripDummyAdj(mN0);
      mMate.assign(2 * mN0, scNull());                 // step 1: M <- empty

      std::vector<std::vector<Size> > maximalOldSmall;
      for (Size b = blossomBase; b < mAuditOldEnd; ++b) {
        if (b >= mZ.size() || mZ[b] <= 0 ||
            b >= mNodeLeaves.size() || nodeLeafCount(b) >= mTau)
          continue;
        bool hasSmallAncestor = false;
        Size p = b < mBParent.size() ? mBParent[b] : scNull();
        Size guard = 0;
        while (p != scNull() && p < mAuditOldEnd && ++guard <= mAuditOldEnd + 2) {
          if (p >= blossomBase && p < mZ.size() && mZ[p] > 0 &&
              p < mNodeLeaves.size() && nodeLeafCount(p) < mTau) {
            hasSmallAncestor = true;
            break;
          }
          p = p < mBParent.size() ? mBParent[p] : scNull();
        }
        if (!hasSmallAncestor) {
          std::vector<Size> leaves;
          collectNodeLeaves(b, &leaves);
          std::vector<Size> realLeaves;
          for (Size i = 0; i < leaves.size(); ++i)
            if (leaves[i] < mN0) realLeaves.push_back(leaves[i]);
          if (!realLeaves.empty()) maximalOldSmall.push_back(realLeaves);
        }
      }

      separateInheritedForest(blossomBase, mAuditOldEnd);
      // Steps 3 and 5, with the materialized reduced-cost step 4 between them.
      for (Size b = blossomBase; b < mAuditOldEnd; ++b)
        if (b < mZ.size() && mZ[b] != 0 && b < mNodeLeaves.size() &&
            nodeLeafCount(b) >= mTau)
          liquidateInherited(b);

      // Step 4, materialized in reduced-cost coordinates.  This cannot be
      // postponed until scale end: step 6 compares vertex duals with the
      // absolute level zero, and step 8 relies on every residual free real
      // vertex reaching y=-tau.  mWoff records the transformation so final
      // reported weights remain in the caller's original coordinates.
      if (mFaithfulReweight) {
        for (Size v = 0; v < mN0; ++v)
          for (Size i = 0; i < mAdj[v].size(); ++i) {
            Size w = mAdj[v][i];
            if (w >= mN0) continue;
            LL shift = 0;
            if (!checkedAddLL(mY[v], mY[w], &shift) ||
                !checkedAddLL(mW[v][i], -shift, &mW[v][i]) ||
                !checkedAddLL(mWoff[v][i], shift, &mWoff[v][i])) {
              finishInheritedForest();
              mLastStatus = STATUS_UNSUPPORTED_OR_OVERFLOWING_INPUT;
              return false;
            }
          }
        for (Size v = 0; v < 2 * mN0; ++v) mY[v] = 0;
      }

      for (Size b = blossomBase; b < mAuditOldEnd; ++b)
        if (b < mZ.size() && mZ[b] != 0 && b < mNodeLeaves.size() &&
            nodeLeafCount(b) < mTau)
          liquidateInherited(b);
      for (Size b = blossomBase; b < mAuditOldEnd && b < mZ.size(); ++b)
        if (mZ[b] != 0) {
          finishInheritedForest();
          mLastStatus = STATUS_INVARIANT_VIOLATION;
          return false;
        }
      finishInheritedForest();
      if (!runLiquidationistLeveling(maximalOldSmall, mN0)) return false;
#ifdef HYBPROF
      // SAME check on the NON-Hybrid path: do transient post-LIQUIDATION domination violations
      // also occur here (and get fixed by the later search)? If YES, the violations are normal
      // and the real Hybrid bug is the search FAILING to fix them; if NO, §4.2 is the culprit.
      { Size dom=0; LL worst=0;
        for (Size u = 0; u < mN0; ++u)
          for (Size i = 0; i < mAdj[u].size(); ++i) {
            Size v = mAdj[u][i]; if (v >= mN0 || v <= u) continue;
            LL sk = slackZ(u, v); if (sk < 0) { ++dom; if (sk<worst) worst=sk; }
          }
        if (dom) std::fprintf(stderr, "[HYBPROF] (NON-HYBRID) post-liquidate: domViol=%lu WORST=%lld\n",(unsigned long)dom,(long long)worst);
      }
#endif
    }

#ifdef SCPHASE
    std::fprintf(stderr, "[SCPHASE]  scale %lu inherited-done clock=%.3fs\n",
                 (unsigned long)scale, (double)clock()/CLOCKS_PER_SEC);
    std::fflush(stderr);
#endif

    // --- Perfection cleanup: UN-PARK reals from their dummies so this scale's
    // search can re-route them under the refined weights. ---
    for (Size u = 0; u < mN0; ++u) {
      Size m = mMate[u];
      if (m != scNull() && m >= mN0) { mMate[u] = scNull(); mMate[m] = scNull(); }
    }

    // Figure 5 step 4 was applied above, between the large and small passes.
    // Figure 6 deliberately omits it.  Nothing remains to do at this point.

    // --- RE-OPTIMIZATION under the refined scale weights (steps 5-8 collapsed to
    // a CORRECTNESS-FIRST realization). The paper restores near-tightness via the
    // bounded step-6 small-blossom searches + tau-bounded SearchOne + dummies, all
    // tuned so each scale costs O(m sqrt n). DIAGNOSTIC PROOF (UNMATCH_ALL, 50/50):
    // if instead we UN-MATCH all real-real edges and run Criterion-1 to a perfect
    // matching under the (warm-started, rescaled) duals, the result is EXACTLY
    // optimal every scale -- the rescaled duals + the (n/2+1) extended weight make
    // the last scale's exact Crit-1 solve THE MWPM (Lemma 2.1/2.2). This is the
    // paper's algorithm with the bounded searches replaced by run-to-completion --
    // correct, just not yet the paper's per-scale speed (that's Layer-1's job: the
    // fast event engine + the warm start make each such solve cheap). Keeping the
    // stale matching + only tau gentle adjustments was the PLATEAU bug (21/50).
    if (mFaithfulLoop && mBoundedSearch) {
      // --- WEIGHT-INDEPENDENT per-scale step (spec Fig-5 steps 5-8, faithful): the
      // re-optimization is the TAU-BOUNDED SearchOne(Crit2) on the REAL graph with FREE reals
      // as roots (NOT a perfect-via-dummies matching -- with everything matched there are no
      // free roots, so the search would be a no-op: the bug found this session). SearchOne =
      // maximal augment + exactly ONE dual adjustment; we call it tau=O(sqrt n) times (spec
      // step 7), so the duals grow by only O(tau) and the carried reduced cost w' stays
      // O(sqrt n)-bounded => weight-INDEPENDENT (the whole point; see
      // scaling-weight-blowup-diagnosed). We do NOT grow to a real perfect matching (the
      // spec MISTAKE-TO-AVOID, and what over-grew y -> blew up w'); free reals are PARKED on
      // dummies in step 8 AFTER the bounded searches; Finalization deletes the dummies and
      // runs Crit3 to the true optimum.
      // Step 1/5: Liquidationist enters FVR with M empty.  Hybrid already set M<-empty
      // BEFORE §4.2 dissolution, whose ShellSearch calls then construct the matching
      // that Figure 6 passes into FVR.  Do not erase that matching a second time: doing
      // so leaves current blossoms contracted while all of their leaves are exposed.
      // In both cases strip dummy adjacency so FVR runs on real edges only.
      stripDummyAdj(mN0);
      // EXPERIMENT (mCarryM via MBX_CARRYM env): do NOT empty M -- carry the real-real matches across
      // scales so the bounded search only REPAIRS instead of rebuilding from empty every scale (the
      // measured cliff = ~n/2 augments x log(nN) scales, the decorative full-rebuild). Free only
      // dummy-parked reals. If this holds 840 + kills the cliff, the full-rebuild WAS the cliff.
      if (mHybrid || mAlgorithmIdentity == ALGORITHM_LIQUIDATIONIST_FIGURE5) {
        // Old dummy slots are outside the real FVR universe.  Their incident edges
        // were stripped above; clear only those slots and retain the real matching
        // produced by DismantlePath/ShellSearch (Hybrid) or Figure-5 step 6
        // (Liquidationist).
        for (Size d = mN0; d < 2 * mN0; ++d) mMate[d] = scNull();
      } else if (mCarryM) {
        for (Size u = 0; u < mN0; ++u) { Size mm = mMate[u]; if (mm != scNull() && mm >= mN0) { mMate[u] = scNull(); mMate[mm] = scNull(); } }
        for (Size d = mN0; d < 2 * mN0; ++d) mMate[d] = scNull();   // clear dummy slots only
      } else {
        mMate.assign(2 * mN0, scNull());      // M <- empty (spec step 1)
      }
      // Step 7/Figure 6 FVR on the real graph (free reals root the search).  The
      // paper budget is measured in unit dual adjustments. Event-free intervals may
      // be batched, but their full delta is charged to the same unit budget.
      // Restrict the universe to real vertices for rooting (mN = mN0).
      { Size sv = mN; mN = mN0;
#ifdef SCMAG
        Size freeStart = 0; for (Size v = 0; v < mN0; ++v) if (mMate[v] == scNull()) ++freeStart;
        Size searchAug = 0;
#endif
        // HYBRID two-stage FVR (Fig 6 steps 2-3): stage 1 = sqrt(n) units of
        // SearchOne(Crit2); the batched implementation is checked against literal calls.
        // stage 2 = BucketSearch(Crit3) until delta=tau or M perfect. Because Hybrid's tau is
        // LARGER (sqrt n log n), literal unit calls would be too slow -- the second
        // stage uses delta-batched Crit3 searches for the remaining tau-sqrt(n) units.
        // (Liquidationist: single stage, tau=sqrt(n) Criterion-2 units.) stage1Cnt is
        // the unit boundary where Hybrid changes from Criterion 2 to Criterion 3.
        Size stage1Cnt = mHybrid ? mSqrtN : mTau;
#ifdef SCMAG
        Size tRan = 0;
#endif
        // D PATH (a): run the whole FVR via the S5 engine bridge (sqrt(n) per-level batching) instead
        // of the tau-loop of edmondsSearchEx calls. One state-exchange per scale. Crit3 (covers
        // Crit2's eligibility as a superset for the grow). Gated OFF by default.
        // DEBUG bisect: MBX_BRIDGE_S1 => bridge ONLY scale 1 (empty forest), tau-loop for scale 2+.
        // If correctness restored under S1-only, the bug is the CARRIED-FOREST load (scale 2+).
        bool useBridge = mBridgeFvr && (!std::getenv("MBX_BRIDGE_S1") || scale == 1);
        LL fvrUnits = 0;
        Size fvrCalls = 0;
        const bool fvrBatch = !mFvrUnitReference &&
                              (mHybridCanonical || !std::getenv("MBX_FVR_UNIT"));
        if (useBridge) runFvrBridged(CRIT3);
        else {
        while (fvrUnits < (LL)mTau) {
#ifdef SCPHASE
          if (fvrCalls < 64 || (fvrCalls & 31) == 0) {
            std::fprintf(stderr,
              "[SCPHASE]  scale %lu FVR call=%lu begin units=%lld queue=%llu/%llu scans=%llu clock=%.3fs\n",
              (unsigned long)scale, (unsigned long)fvrCalls,
              (long long)fvrUnits, mCounters.queueInsertions,
              mCounters.queuePops, mCounters.edgeScans,
              (double)clock()/CLOCKS_PER_SEC);
            std::fflush(stderr);
          }
#endif
#ifdef SCMAG
          tRan = fvrCalls + 1;
#endif
          // EARLY EXIT: FVR's purpose is to REDUCE free vertices; once none are free, further
          // SearchOne/BucketSearch calls only churn duals on a fully-matched graph (pure waste --
          // profiling showed scale 3 doing ALL 684 searches with freeReals already 0, = 3.45s of
          // the 3.8s at n=3200). The paper runs "tau times" as an UPPER bound on the work needed
          // to drive free below O(n/tau); finishing early is correct (the remaining adjustments
          // would change nothing structural). Re-checked each iteration: O(n) scan, dwarfed by a
          // search. (Park-on-dummies in step 8 handles any residual free; Finalization optimizes.)
          { bool anyFree = false; for (Size v = 0; v < mN0; ++v) if (mMate[v] == scNull()) { anyFree = true; break; }
            if (!anyFree) break; }
          bool useCrit3 = mHybrid && (fvrUnits >= (LL)stage1Cnt); // switch after sqrt(n) DUAL UNITS
          LL phaseEnd = useCrit3 ? (LL)mTau : (LL)stage1Cnt;
          LL unitBudget = phaseEnd - fvrUnits;
          bool a = false;
          SearchRequest searchRequest;
          searchRequest.criterion = useCrit3 ? CRIT3 : CRIT2;
          searchRequest.dualBudget = unitBudget;
          searchRequest.maximalAugment = true;
          searchRequest.batchEventFreeTime = fvrBatch;
          searchRequest.persistAcrossAdjust = mPersistSearch;
          searchRequest.persistAcrossAugment = mPersistSearch;
          searchRequest.exhaustBudget = true;
          searchRequest.stoppingRule = STOP_AT_DUAL_BUDGET;
          if (useCrit3) {
            // Figure 6 stage 2 is one BucketSearch that keeps its alternating
            // forest until it augments or consumes the remaining tau budget;
            // elapsed delta is charged exactly. The event backend expresses
            // that directly. Canonical's conservative scan oracle may rebuild
            // and re-enter after each adjustment without changing M/y/z.
            searchRequest.maxAdjustments = mEventSearch ? -1 : 1;
            a = executeSearch(searchRequest,
                              mEventSearch ? BACKEND_EVENT : BACKEND_SCAN).augmented;
#ifdef SCMAG
            if (a) ++searchAug;
#endif
          } else {
          // SearchOne = maximal augment + ONE dual adjustment. The "one adjustment" is one
          // UNIT of paper time.  We may batch an event-free interval, but the jump distance
          // is charged against the sqrt(n)/tau budget by mLastSearchDualUnits; it is never
          // counted as one unit merely because it was one implementation call.
          // EDM WIN (opt-in via setPersistSearch): the PERSISTENT-forest variant
          // (persistAcrossAdjust+Augment=true) skips the from-scratch relabel+regrow
          // (O(n^2+nm)/search). Profiling: iters grew O(n) + time super-quadratic while the
          // dual-adjust count stayed flat -- the cost was the per-search REBUILD. Persistent
          // gives 3.7x-22x (n=800: 0.17->0.05s; n=3200: 6.5->0.3s), correctness 180/180. BUT it
          // INFINITE-LOOPS on some instances (a persistent-path no-progress spin, e.g. one n=1600
          // seed, while n=3200 finishes in 0.3s) -- so it is GATED OFF by default. The safe
          // default uses the non-persistent search (always terminates, just slower at scale).
          // EVENT-SEARCH (opt-in, mEventSearch): use the Stage-3 event-PQ search (paper §5
          // BucketSearch, validated 150/150) -- it pops the next tight-edge EVENT instead of
          // rescanning all m edges per round, which is the O(m)-per-iteration cost that made
          // the bounded loop super-linear in n (n=6400: 35s -> event 3.3s). It advances duals
          // by delta-event (not unit step) so it stays weight-robust; inside the scaling loop
          // the per-scale weights are small (post-reweight), so its weight-dependence does not
          // bite. NOTE: event search relabels fresh per call (no cross-call persist), so it does
          // its own full search to an augment each call -- still O(events) per call, not O(m*iters).
          if (mPersistSearch) {
            // maxAdjust: how many delta-batched events ONE rooted search applies before returning.
            // =1 means the tau-loop re-roots+rescans (O(m)) per dual jump -> ~900 wasted non-augmenting
            // searches/scale in pathological seeds (MEASURED: late-scale edgeScans ~40M = 500x m). A
            // LARGER maxAdjust amortizes the O(m) root over many dual jumps within one persistent
            // search. MBX_MAXADJ overrides (default 1 = old behavior; experiment with larger).
            LL maxAdj = 1;
            if (!mHybridCanonical) {
              const char* e = std::getenv("MBX_MAXADJ");
              if (e) maxAdj = atoll(e);
            }
            searchRequest.maxAdjustments = maxAdj;
            a = executeSearch(searchRequest,
                              mEventSearch ? BACKEND_EVENT : BACKEND_SCAN).augmented;
          } else {
            searchRequest.maxAdjustments = 1;
            a = executeSearch(searchRequest,
                              mEventSearch ? BACKEND_EVENT : BACKEND_SCAN).augmented;
          }
#ifdef SCMAG
          if (a) ++searchAug;
#endif
          }
          ++fvrCalls;
          LL spent = mLastSearchDualUnits;
          if (spent > unitBudget) spent = unitBudget; // defensive; the search clamps internally
          fvrUnits += spent;
#ifdef SCPHASE
          if (fvrCalls <= 64 || (fvrCalls & 31) == 0) {
            std::fprintf(stderr,
              "[SCPHASE]  scale %lu FVR call=%lu end aug=%d spent=%lld units=%lld queue=%llu/%llu scans=%llu clock=%.3fs\n",
              (unsigned long)scale, (unsigned long)fvrCalls, a ? 1 : 0,
              (long long)spent, (long long)fvrUnits,
              mCounters.queueInsertions, mCounters.queuePops,
              mCounters.edgeScans, (double)clock()/CLOCKS_PER_SEC);
            std::fflush(stderr);
          }
#endif
#ifdef SCMAG
          if (std::getenv("MBX_TRACE_FVR_CALLS")) {
            std::fprintf(stderr,
              "[SCMAG]   scale%lu FVR-call%lu crit=%d aug=%d spent=%lld total=%lld free-y:",
              (unsigned long)scale, (unsigned long)fvrCalls,
              useCrit3 ? 3 : 2, a ? 1 : 0, (long long)spent, (long long)fvrUnits);
            for (Size v=0; v<mN0; ++v) if (mMate[v]==scNull())
              std::fprintf(stderr," v%lu=%lld(o%lu,l%d)",(unsigned long)v,
                (long long)mY[v],(unsigned long)outerOf(v),mNLabel[outerOf(v)]);
            std::fprintf(stderr,"\n");
          }
#endif

          if (mLastStatus != STATUS_SUCCESS) break;
          // A maximal same-time search may complete the matching while its
          // boolean reports no additional path from the terminal forest.  The
          // Figure-6 stop condition is the canonical matching state itself.
          bool fvrPerfect = true;
          for (Size v = 0; v < mN0; ++v)
            if (mMate[v] == scNull()) { fvrPerfect = false; break; }
          if (fvrPerfect) break;
          // With free roots and remaining FVR budget, Figure 6 permits neither a
          // silent halt nor an implementation progress cap.  The search must
          // augment or advance the dual clock; otherwise its state is invalid.
          if (!a && spent == 0) {
            mLastStatus = STATUS_UNEXPECTED_NO_PROGRESS;
            break;
          }
        }
        // A delta jump may make Criterion-3 events eligible exactly at the
        // delta=tau boundary.  BucketSearch processes events at the current
        // time before testing whether another dual adjustment is permitted.
        // Close that boundary without charging an additional dual unit.
        if (mHybrid && fvrUnits == (LL)mTau) {
          bool anyFree = false;
          for (Size v = 0; v < mN0; ++v)
            if (mMate[v] == scNull()) { anyFree = true; break; }
          if (anyFree) {
            SearchRequest boundary;
            boundary.criterion = CRIT3; boundary.maximalAugment = true;
            boundary.maxAdjustments = 0; boundary.dualBudget = 0;
            boundary.batchEventFreeTime = true;
            boundary.stoppingRule = STOP_AT_DUAL_BUDGET;
            executeSearch(boundary,
                          mEventSearch ? BACKEND_EVENT : BACKEND_SCAN);
          }
        }
        }
#ifdef SCMAG
        Size freeEnd = 0; for (Size v = 0; v < mN0; ++v) if (mMate[v] == scNull()) ++freeEnd;
        std::fprintf(stderr, "[SCMAG]   step7: tauCap=%lu calls=%lu dualUnits=%lld, freeReals %lu->%lu, augmenting-searches=%lu\n",
                     (unsigned long)mTau, (unsigned long)tRan, (long long)fvrUnits,
                     (unsigned long)freeStart, (unsigned long)freeEnd, (unsigned long)searchAug);
        if (scale == 1 || std::getenv("MBX_TRACE_FVR_PARITY")) {
          std::fprintf(stderr, "[SCMAG]   scale%lu-FVR y:", (unsigned long)scale);
          for (Size v = 0; v < mN0; ++v) std::fprintf(stderr, " %lld", (long long)mY[v]);
          std::fprintf(stderr, "\n[SCMAG]   scale%lu-FVR mate:", (unsigned long)scale);
          for (Size v = 0; v < mN0; ++v) std::fprintf(stderr, " %ld", mMate[v] == scNull() ? -1L : (long)mMate[v]);
          std::fprintf(stderr, "\n");
        }
#endif
        mN = sv;
        if (mLastStatus != STATUS_SUCCESS) return false;
      }
      // Step 8 (Perfection): park each still-free real on its dummy (0-weight matched edge),
      // y(dummy)=tau. These parked slots are re-optimizable in Finalization.
      for (Size u = 0; u < mN0; ++u) {
        if (mMate[u] != scNull()) continue;
        if (mAlgorithmIdentity == ALGORITHM_LIQUIDATIONIST_FIGURE5 &&
            mY[u] != -(LL)mTau) {
          mLastStatus = STATUS_INVARIANT_VIOLATION;
          return false;
        }
        Size d = mN0 + u;
        ++mCounters.perfectionParkedVertices;
        // mWbar[u] MUST track mAdj[u] (the next scale's Scaling step reads mWbar[u][i]
        // for i<mAdj[u].size()); dummy edges have wbar=0. Failing to grow mWbar here was
        // a heap-buffer-overflow at the scaling step (ScalingMatching.h:2174).
        mAdj[u].push_back(d); mW[u].push_back(0); mWhat[u].push_back(0); mWbar[u].push_back(0);
        mAdj[d].push_back(u); mW[d].push_back(0); mWhat[d].push_back(0);
        // Property 2 requires the newly matched zero-weight pendant edge to be
        // tight: y(u)+y(d)=0 (it shares no blossom). Liquidationist has
        // y(u)=-tau here, yielding the paper's y(d)=tau. Hybrid can retain a
        // nonzero pre-FVR free dual, so use the equivalent invariant form while
        // the faithful Omega_0/Omega prototype is enabled.
        // Liquidationist has y(u)=-tau here and therefore writes y(d)=tau.
        // Hybrid's no-reweighting path guarantees common free-dual parity, not
        // y(u)=-tau (Lemma 4.1 and the y6 term in Lemma 3.7).  "Perfection as
        // in Liquidationist" must preserve the operative condition: the new
        // zero-weight matched edge is tight.
        LL dummyY = (LL)mTau;
        if (mHybridCanonical ||
            (mHybrid && std::getenv("MBX_SEPARATE_OLD_FOREST"))) {
          // Figure 5's literal y(dummy)=tau is inconsistent with reachable
          // no-reweighting Hybrid states (and with its own proof that the new
          // zero edge is tight). Preserve that proof invariant explicitly.
          if (mY[u] == std::numeric_limits<LL>::min()) {
            mLastStatus = STATUS_UNSUPPORTED_OR_OVERFLOWING_INPUT;
            return false;
          }
          dummyY = -mY[u];
          if (mY[u] != -(LL)mTau)
            ++mCounters.perfectionNonTauCorrections;
        }
        mY[d] = dummyY;
        LL pendantDual = 0;
        if (!checkedAddLL(mY[u], mY[d], &pendantDual) || pendantDual != 0) {
          ++mCounters.perfectionTightnessFailures;
          mLastStatus = STATUS_INVARIANT_VIOLATION;
          return false;
        }
        mMate[u] = d; mMate[d] = u;
      }
    } else if (mFaithfulLoop) {
      // --- FAITHFUL per-scale step (steps 5-8, using the PROVEN single-scale mechanism:
      // real-edge search + park-on-dummies, NO eager dummies). Spec Fig-5: M is rebuilt
      // each scale via the search under the CARRIED+rescaled duals; reals that can't match
      // on real tight edges get PARKED on dummies (step 8). Finalization (after the LAST
      // scale) un-parks + re-optimizes to the true optimum. We do NOT drive each scale to
      // a real perfect matching (the spec's MISTAKE-TO-AVOID) and do NOT use dummies in the
      // search (the n=8-rep17 premature-park bug). Matching is reset M<-empty (spec line 69)
      // then rebuilt from the carried duals. ---
      // strip any dummy edges/parks from the previous scale's adjacency (Perfection added
      // them lazily); the new scale searches on real edges only.
      stripDummyAdj(mN0);
      mMate.assign(2 * mN0, scNull());      // M <- empty (spec step 1)
      // Phase A: real-edge search under carried duals, to a real PM if possible.
      // NOTE: this uses unbounded Criterion 1 (correct, 180/180 at field-norm weights).
      // The spec's step-7 BOUNDED SearchOne(Crit2) is what would give weight-magnitude
      // INDEPENDENCE (Crit2's Property-2 bounds all slacks -> reweight bounds magnitude),
      // but swapping it in (edmondsSearch CRIT2 maxAdjust=1, tau times) did NOT bound the
      // magnitude in practice here (mW still doubled / hung at large N) -- the gap between
      // our edmondsSearch and the spec's exact Property-2 maintenance is unclosed. So we
      // keep the CORRECT Crit1 loop (exact MWPM for weights up to the field norm ~1e6);
      // full weight-independence to 1e9+ remains the open piece.
#ifdef SCMAG
      if (scale==2){ std::fprintf(stderr,"  [LOOP scale2 edge w2 (=mW)]:"); for(Size u=0;u<mN0;++u) for(Size k=0;k<mAdj[u].size();++k) if(mAdj[u][k]>u && mAdj[u][k]<mN0) std::fprintf(stderr," (%lu-%lu:%lld)",(unsigned long)u,(unsigned long)mAdj[u][k],(long long)mW[u][k]); std::fprintf(stderr,"\n"); }
#endif
      { Size sv = mN; mN = mN0;
        for (Size gg = 0; gg <= 4 * mN0 + 8; ++gg) {
          bool any = false;
          for (Size v = 0; v < mN0; ++v) if (mMate[v] == scNull()) { any = true; break; }
          if (!any) break;
          if (!edmondsSearch(CRIT1, /*maximalAugment=*/true, /*maxAdjust=*/-1)) break;
        }
        mN = sv;
      }
#ifdef SCMAG
      if (scale <= 3) { LL ys=0; for (Size v=0; v<mN0; ++v) ys+=mY[v];
        std::fprintf(stderr,"[LOOP scale %lu post-search] sumY=%lld y:", (unsigned long)scale,(long long)ys);
        for (Size v=0; v<mN0; ++v) std::fprintf(stderr," %lld",(long long)mY[v]);
        std::fprintf(stderr,"\n  mate:"); for (Size v=0;v<mN0;++v) std::fprintf(stderr," %ld",(long)(mMate[v]==scNull()?-1:(long)mMate[v]));
        std::fprintf(stderr,"\n"); }
#endif
      // Phase B (Perfection, step 8): park each still-free real on its dummy (lazily create
      // the 0-weight matched edge), y(dummy)=tau. These parked slots are re-optimizable.
      for (Size u = 0; u < mN0; ++u) {
        if (mMate[u] != scNull()) continue;
        Size d = mN0 + u;
        mAdj[u].push_back(d); mW[u].push_back(0); mWhat[u].push_back(0); mWbar[u].push_back(0);  // mWbar tracks mAdj
        mAdj[d].push_back(u); mW[d].push_back(0); mWhat[d].push_back(0);
        mY[d] = (LL)mTau;
        mMate[u] = d; mMate[d] = u;
      }
    } else {
    for (Size u = 0; u < mN0; ++u)
      if (mMate[u] != scNull() && mMate[u] < mN0) { mMate[mMate[u]] = scNull(); mMate[u] = scNull(); }
    if (mUseFastEngine) {
      // LAYER 2: re-solve this scale on the FAST event engine (UF + lazy duals +
      // event PQ + sqrt(n) multi-augment) instead of the old rebuild edmondsSearch.
      // Each scale is a full MWPM solve under the scale weights mW (the correctness-
      // first realization re-solves per scale); the fast engine makes that cheap.
      ScalingFastMWPM fe; fe.setUseEvent(true); fe.init(mN0);
      for (Size v = 0; v < mN0; ++v)
        for (Size i = 0; i < mAdj[v].size(); ++i) {
          Size w = mAdj[v][i];
          if (w < mN0 && w > v) fe.addEdge(v, w, mW[v][i] / 2);   // mW is doubled
        }
      std::vector<Size> fmate;
      LL fg = fe.solve(&fmate);
      if (fg >= 0) for (Size v = 0; v < mN0; ++v) mMate[v] = fmate[v];
    } else {
      Size sv = mN; mN = mN0;            // search over real vertices only
      for (Size gg = 0; gg <= 2 * mN0 + 4; ++gg) {
        bool perfect = true;
        for (Size v = 0; v < mN0; ++v) if (mMate[v] == scNull()) { perfect = false; break; }
        if (perfect) break;
        if (!edmondsSearch(CRIT1, /*maximalAugment=*/false, /*maxAdjust=*/-1)) break;
      }
      mN = sv;
    }
    }

    // Figure 5's materialized step-4 reweighting already occurred between the
    // large- and small-blossom passes. Hybrid intentionally omits reweighting.
#ifdef SCMAG
    {
      // Property-2 check: slack = yz - w should be in [-2, 0] on ALL edges if the
      // search maintained near-domination (>= -2) and near-tightness (<= 0). If
      // slacks are large-positive, the reduced cost w-yz is large-negative -> w
      // can't be bounded by reweight. This pinpoints whether Property 2 holds.
      LL maxS = 0, minS = 0; bool first = true;
      LL maxMatchedS = 0; LL nBig = 0; LL nEdges = 0;   // matched-edge slack + count of |slack|>4
      for (Size v = 0; v < mN0; ++v)
        for (Size i = 0; i < mAdj[v].size(); ++i) {
          Size w = mAdj[v][i]; if (w >= mN0) continue;
          LL sk = slackZ(v, w);
          ++nEdges; if (sk > 4 || sk < -4) ++nBig;
          if (mMate[v] == w) { LL a = sk < 0 ? -sk : sk; if (a > maxMatchedS) maxMatchedS = a; }
          if (first) { maxS = minS = sk; first = false; }
          else { if (sk > maxS) maxS = sk; if (sk < minS) minS = sk; }
        }
      std::fprintf(stderr, "[SCMAG]   pre-reweight slack range [%lld, %lld]  matched|slack|max=%lld  bigEdges(|s|>4)=%lld/%lld\n",
                   minS, maxS, maxMatchedS, nBig, nEdges);
    }
#endif
    // Carry w' <- w (the REWEIGHTED current w) and y' <- y (post-search y, carried
    // in place) to the next scale. The paper's recurrence works in yz-space: the
    // next scale does w<-2(w'+bit), y<-2y'+3, so the reweight offset rides along in
    // y and stays consistent. (Carrying the un-reweighted w instead breaks
    // domination because y reflects the reweighted world.)
    for (Size v = 0; v < mN0; ++v)
      for (Size i = 0; i < mAdj[v].size() && i < wPrev[v].size(); ++i)
        wPrev[v][i] = mW[v][i];
#ifdef SCMAG
    {
      LL maxW = 0, maxY = 0;
      for (Size v = 0; v < mN0; ++v) {
        if (mY[v] > maxY) maxY = mY[v]; if (-mY[v] > maxY) maxY = -mY[v];
        for (Size i = 0; i < mAdj[v].size(); ++i) {
          LL ww = mW[v][i]; if (ww > maxW) maxW = ww; if (-ww > maxW) maxW = -ww;
        }
      }
      std::fprintf(stderr, "[SCMAG] scale %lu/%lu  max|w|=%lld  max|y|=%lld\n",
                   (unsigned long)scale, (unsigned long)nbits, maxW, maxY);
    }
#endif
  }

#ifdef SCPHASE
  std::fprintf(stderr, "[SCPHASE] scale-loop done at clock=%.3fs\n", (double)clock()/CLOCKS_PER_SEC);
#endif
#ifdef SCPROF
  { extern long long gRelinkOK, gRelinkBail, gRelinkB[8], gRelabelShell, gRelabelMain;
    std::fprintf(stderr, "[RELINK] ok=%lld bail=%lld  reasons[oddlen=%lld noEntry=%lld noAttach=%lld multiChild=%lld]\n",
                 gRelinkOK, gRelinkBail, gRelinkB[0], gRelinkB[1], gRelinkB[2], gRelinkB[3]);
    std::fprintf(stderr, "[RELABEL-SPLIT] inShellSearch=%lld inMainFVR=%lld\n", gRelabelShell, gRelabelMain); }
#endif
  // --- Finalization: delete all dummy vertices; PQSearch(F) Criterion 3 until
  // every real vertex is matched on REAL edges. ---
  for (Size u = 0; u < mN0; ++u) {
    Size m = mMate[u];
    if (m != scNull() && m >= mN0) { mMate[u] = scNull(); mMate[m] = scNull(); }
  }
  // Drop dummy adjacency so the search can't use 0-weight dummy edges (in place).
  stripDummyAdj(mN0);
  for (Size d = mN0; d < 2 * mN0; ++d) mMate[d] = scNull();   // (stripDummyAdj cleared their adj)
  if (mSkipFinalizationExperimental) {
    bool partial = false;
    for (Size v = 0; v < mN0; ++v)
      if (mMate[v] == scNull()) { partial = true; break; }
    mN = mN0;
    mNoPerfect = partial;
    mLastStatus = partial ? STATUS_SUCCESS_EXPERIMENTAL_PARTIAL : STATUS_SUCCESS;
    return true;
  }
  // Now run Criterion-3 searches to perfection on real vertices only.
  // (Restrict the vertex universe back to real vertices for rooting.)
  // pqDelta: the bounded (weight-independent) path uses delta-BATCHED dual adjustment so
  // Finalization stays O(events)=O(m) per search instead of O(maxW) unit steps -- without
  // this, Finalization is the pseudo-polynomial spike (some seeds at W=1e4 needed ~87k unit
  // adjustments vs ~60 with batching). The non-bounded faithful path keeps the proven
  // unit-step behavior (pqDelta=false) so its 180/180 gate is untouched.
  bool finalPQ = mBoundedSearch;
  // Finalization now KEEPS incr-shrink ON (it is EXACT -- 840/840 -- so no need to force the slow
  // relabel path here). Disabling it forced the O(depth)+full-relabel path during Finalization,
  // which at n=25600 was ~7.7s of 9.8s (the dominant cost, missed because per-scale profiling
  // stops at scale-loop end). Keeping incr-shrink = the no-relabel + O(alpha) fast path here too.
  bool savedIncr = mIncrShrink;   /* keep mIncrShrink as-is (ON) for Finalization */
#ifdef SCPROF
  gProfInFinal = true;
#endif
  Size savedN = mN; mN = mN0;          // search roots only real free vertices now
  for (Size guard = 0; guard <= 2 * mN0 + 4; ++guard) {
    bool perfect = true;
    for (Size v = 0; v < mN0; ++v) if (mMate[v] == scNull()) { perfect = false; break; }
    if (perfect) { mN = savedN; mNoPerfect = false;
#ifdef SCPROF
      gProfInFinal = false;
#endif
      return mLastStatus == STATUS_SUCCESS; }
#ifdef SCPROF
    ++gProfFinalSrch;
#endif
#ifdef SCMAG
    { Size fr = 0; for (Size v = 0; v < mN0; ++v) if (mMate[v] == scNull()) ++fr;
      std::fprintf(stderr, "[SCMAG]   FINAL search %lu: freeReals=%lu\n",
                   (unsigned long)guard, (unsigned long)fr);
      for (Size v=0; v<mN0; ++v) if (mMate[v]==scNull())
        { std::fprintf(stderr,"    free v%lu outer=%lu uf=%lu y=%lld\n",(unsigned long)v,
          (unsigned long)outerOf(v),(unsigned long)ufFind(v),(long long)mY[v]);
          for (Size k=0;k<mAdj[v].size();++k) { Size w=mAdj[v][k]; if(w<mN0)
            std::fprintf(stderr,"      e(%lu,%lu) slack=%lld outerW=%lu ufW=%lu\n",
              (unsigned long)v,(unsigned long)w,(long long)slackZ(v,w),
              (unsigned long)outerOf(w),(unsigned long)ufFind(w)); } }
      for (Size b=2*mN0; b<mNumNodes; ++b)
        if (mBParent[b]==scNull() && mZ[b]!=0)
          std::fprintf(stderr,"    root B%lu base=%lu z=%lld leaves=%lu\n",
            (unsigned long)b,(unsigned long)mBBase[b],(long long)mZ[b],
            (unsigned long)nodeLeafCount(b)); }
#endif
    SearchRequest finalRequest;
    finalRequest.criterion = CRIT3;
    finalRequest.maxAdjustments = -1;
    finalRequest.batchEventFreeTime = finalPQ;
    finalRequest.stoppingRule = STOP_AFTER_AUGMENTATION;
#ifdef PS5_FIN_MAXAUG
    finalRequest.maximalAugment = true;
#endif
    if (!executeSearch(finalRequest,
                       mEventSearch ? BACKEND_EVENT : BACKEND_SCAN).augmented) break;
  }
#ifdef SCPROF
  gProfInFinal = false;
#endif
  mIncrShrink = savedIncr;   // restore (experiment: incr-shrink off only during Finalization)
  mN = savedN;
  // Check perfection on real vertices.  The named drop-perfect experiment is
  // the sole exception: when the ordinary Finalization search has exhausted
  // augmenting paths, publish that partial matching for quality/runtime
  // measurement.  Do not silently grant this behavior to canonical Hybrid.
  for (Size v = 0; v < mN0; ++v) if (mMate[v] == scNull()) {
    mNoPerfect = true;
    if (mDropPerfectExperimental && mLastStatus == STATUS_SUCCESS) {
      mLastStatus = STATUS_SUCCESS_EXPERIMENTAL_PARTIAL;
      return true;
    }
    mLastStatus = STATUS_NO_PERFECT_MATCHING;
    return false;
  }
  mNoPerfect = false;
  return mLastStatus == STATUS_SUCCESS;
}

// ISOLATION TEST for the re-optimization mechanism. ONE scale at full weight
// w=2*what (skips the bit-by-bit scaling; w-bar's role is only the optimality
// margin, unneeded for a single full-weight solve). Builds M leaving free reals
// PARKED on 0-weight dummies, then DELETES dummies and runs Criterion-3 search
// to the optimum. If this matches LEMON, the dummy+dual-adjust+Finalization
// re-optimization works and the multi-scale wrapper is the remaining glue.
inline bool ScalingMatching::runFaithfulSingleScale() {
  Size mN0 = mN;
  mNoPerfect = false;
  if (mN0 == 0) return true;
  if (mN0 % 2 != 0) { mNoPerfect = true; return false; }

  Size blossomBase = 2 * mN0;
  Size maxNodes = 4 * mN0 + 4;
  mN = blossomBase;                          // vertices [0,2mN0): real + dummy slots
  mAdj.resize(2 * mN0); mW.resize(2 * mN0); mWhat.resize(2 * mN0);
  for (Size d = mN0; d < 2 * mN0; ++d) { mAdj[d].clear(); mW[d].clear(); mWhat[d].clear(); }
  // full weight w = 2*what on real edges. DUMMY edges are NOT added here: per the spec
  // (Fig-5 steps 7/8 + the MISTAKE-TO-AVOID note), dummies are used ONLY in Perfection
  // (step 8), AFTER the main search -- NOT eagerly during it. Adding them before Phase A
  // let a real vertex get PARKED on its dummy mid-search (when its y hit 0 the 0-weight
  // dummy edge went tight and Phase A augmented onto it), freezing it out of the search
  // and blocking re-optimization (the n=8 rep17 -302-domination-violation bug). So we add
  // each dummy edge lazily, only when Phase B actually parks that vertex.
  for (Size v = 0; v < mN0; ++v) {
    mW[v].assign(mAdj[v].size(), 0);
    for (Size i = 0; i < mAdj[v].size(); ++i) mW[v][i] = 2 * mWhat[v][i];
  }
  mUseScaleW = true;

  mY.assign(2 * mN0, 0);
  for (Size v = 0; v < mN0; ++v) mY[v] = mMaxW;   // domination for w=2*what
  // dummy y: 0 so dummy edge slack = mY[real]+0-0 = mY[real] >= 0 (dominated).
  mZ.assign(maxNodes, 0);
  mMate.assign(2 * mN0, scNull());
  mBParent.assign(maxNodes, scNull());
  mBBase.assign(maxNodes, scNull());
  mBCycle.assign(maxNodes, std::vector<Size>());
  mBCycEdgeU.assign(maxNodes, std::vector<Size>());
  mBCycEdgeV.assign(maxNodes, std::vector<Size>());
  mNodeLeaves.assign(maxNodes, std::vector<Size>());
  mNLabel.assign(maxNodes, LBL_FREE_OF_TREE);
  mNParent.assign(maxNodes, scNull());
  mNParEdgeU.assign(maxNodes, scNull());
  mNParEdgeV.assign(maxNodes, scNull());
  mNRoot.assign(maxNodes, scNull());
  for (Size v = 0; v < 2 * mN0; ++v) { mBBase[v] = v; mNodeLeaves[v].assign(1, v); }
#ifdef PS5_LEAF_SPLICE
  scalingLfInit(2 * mN0, maxNodes);
#endif
  mNumNodes = blossomBase;

  // Phase A: search to a perfect matching using CRITERION 1 (exact tightness).
  // Criterion 1 maintains PROPERTY 1, so a perfect matching it produces is the
  // exact MWPM (Lemma 2.1). [Criterion 2 only maintains Property 2 = within-n of
  // optimal -- that relaxation is for the SCALING LOOP's speed, not for a single
  // full-weight solve. For the standalone exact reference we use Criterion 1.]
  // (The +2 matched-edge drift bug from the z-clamp in dual adjust is now fixed,
  // so Criterion 1 reaches a true Property-1 perfect matching.)
  for (Size t = 0; t < 4 * mN0 + 8; ++t) {
    bool any = false;
    for (Size v = 0; v < mN0; ++v) if (mMate[v] == scNull()) { any = true; break; }
    if (!any) break;
    if (!edmondsSearch(CRIT1, /*maximalAugment=*/true, /*maxAdjust=*/-1, /*pqDelta=*/false)) break;
  }

  // Phase B: Perfection -- park any still-free real on its dummy (0-weight match).
  // Now that dummies were NOT in the search, create the dummy edge HERE, only for the
  // reals that remained free after the real-edge search. y(dummy)=tau (spec step 8).
  for (Size u = 0; u < mN0; ++u) {
    if (mMate[u] != scNull()) continue;
    Size d = mN0 + u;
    mAdj[u].push_back(d); mW[u].push_back(0); mWhat[u].push_back(0); mWbar[u].push_back(0);  // mWbar tracks mAdj
    mAdj[d].push_back(u); mW[d].push_back(0); mWhat[d].push_back(0);
    mY[d] = (LL)mTau;
    mMate[u] = d; mMate[d] = u;
  }

  // Phase C: Finalization -- DELETE dummies (free their parked reals), strip dummy
  // edges, then Criterion-3 search to perfection on real vertices -> optimum.
  for (Size u = 0; u < mN0; ++u) {
    Size m = mMate[u];
    if (m != scNull() && m >= mN0) { mMate[u] = scNull(); mMate[m] = scNull(); }
  }
  stripDummyAdj(mN0);                        // strip dummy edges in place
  for (Size d = mN0; d < 2 * mN0; ++d) mMate[d] = scNull();
  mN = mN0;                                  // root only real free vertices
  for (Size guard = 0; guard <= 2 * mN0 + 4; ++guard) {
    bool perfect = true;
    for (Size v = 0; v < mN0; ++v) if (mMate[v] == scNull()) { perfect = false; break; }
    if (perfect) break;
    if (!edmondsSearch(CRIT3, /*maximalAugment=*/false, /*maxAdjust=*/-1, /*pqDelta=*/false)) break;
  }
  for (Size v = 0; v < mN0; ++v) if (mMate[v] == scNull()) { mNoPerfect = true; return false; }
  mNoPerfect = false;
  return true;
}

// ===========================================================================
// S5 PERSISTENT-FOREST engine (paper S5; docs/Scaling_S5_DataStructures.md).
// Same MWPM result as runPerfectBlossom but NO per-iteration relabel. The forest
// PERSISTS across dual adjustments; a dual adjustment is a DELTA-JUMP applied
// incrementally to labeled nodes only (O(#labeled), not O(n)). Relabel happens
// ONLY after an augment (tree consumed) or a dissolve (sub-nodes freed) -- both
// O(#augments+#dissolves) total, not per dual adjustment. This removes the
// O(mNumNodes)-per-dual-adjustment rebuild that made scaling blow up (n=12800:
// 20s). Reuses shrink/augmentFrom/rotateToBase/lcaOuter/eligible/slackZ unchanged
// (they read mY/mZ, which we maintain incrementally here).
//
// STEP 2 of the S5 plan: persistent forest + delta-jump incremental duals (still
// scans EVEN edges to find the next delta + newly-tight edges -- O(m) per round,
// O(#events)=O(n) rounds => O(mn), already a big win over the rebuild). STEP 3
// swaps the O(m) delta-scan for the BucketPQ; STEP 4 adds union-find/split-findmin.
// ===========================================================================

inline bool ScalingMatching::edmondsSearchS5(int crit) {
  bool augmentedAny = false;
  bool needRelabel = true;
  LL guardAdj = 0;
  LL guardCap = 16 * (LL)mN * (LL)mN + 8 * (LL)mN * (mMaxW + 1) + 1000;
#ifdef SCDBG
  LL dbgLoop = 0;
#endif
  for (;;) {
#ifdef SCDBG
    if (++dbgLoop % 200 == 0) { Size cc=0; for(Size v=0;v<mN;++v) if(mMate[v]!=scNull())++cc; std::fprintf(stderr,"  [S5 loop %ld] matched=%lu nodes=%lu guardAdj=%ld relabel=%d\n",(long)dbgLoop,(unsigned long)cc,(unsigned long)mNumNodes,(long)guardAdj,(int)needRelabel); std::fflush(stderr); if(dbgLoop>3000) return augmentedAny; }
#endif
    if (needRelabel) {
      // repair malformed blossoms (free non-base member) then relabel forest.
      for (Size v = 0; v < mN; ++v)
        if (mMate[v] == scNull()) {
          Size N = outerOf(v);
          LL rg = 0;
          while (N >= mN && mBBase[N] != v) {
#ifdef SCDBG
            if (++rg > 100) { std::fprintf(stderr,"  [S5 REPAIR STUCK] v=%lu outer=%lu base=%lu mBParent[N]=%ld cyc=%lu\n",(unsigned long)v,(unsigned long)N,(unsigned long)mBBase[N],(long)mBParent[N],(unsigned long)mBCycle[N].size()); std::fflush(stderr); break; }
#endif
            flattenNode(N); N = outerOf(v);
          }
        }
      for (Size i = 0; i < mNumNodes; ++i) {
        if (mBParent[i] != scNull()) continue;
        mNLabel[i] = LBL_FREE_OF_TREE;
        mNParent[i] = scNull(); mNParEdgeU[i] = scNull(); mNParEdgeV[i] = scNull();
        mNRoot[i] = scNull();
      }
      bool anyRoot = false;
      for (Size v = 0; v < mN; ++v)
        if (mMate[v] == scNull()) {
          Size N = outerOf(v);
          if (mNLabel[N] == LBL_FREE_OF_TREE) { mNLabel[N]=LBL_EVEN; mNRoot[N]=N; anyRoot=true; }
        }
      if (!anyRoot) return augmentedAny;
      needRelabel = false;
    }

    // ---- GROW to fixpoint (worklist of EVEN outer nodes). Labels PERSIST. ----
    bool didAugment = false, didShrink = false;
    std::vector<Size> wq;
    for (Size i = 0; i < mNumNodes; ++i)
      if (mBParent[i] == scNull() && mNLabel[i] == LBL_EVEN) wq.push_back(i);
    Size wh = 0;
    while (wh < wq.size() && !didShrink && !didAugment) {
      Size E = wq[wh++];
      if (mBParent[E] != scNull() || mNLabel[E] != LBL_EVEN) continue;
      std::vector<Size> Everts; collectNodeLeaves(E, &Everts);
      for (Size ui = 0; ui < Everts.size() && !didShrink && !didAugment; ++ui) {
        Size u = Everts[ui];
        for (Size k = 0; k < mAdj[u].size() && !didShrink && !didAugment; ++k) {
          Size w = mAdj[u][k];
          Size W = outerOf(w);
          if (W == E) continue;
          if (!eligible(u, w, crit)) continue;
          if (mNLabel[W] == LBL_EVEN) {
            if (lcaOuter(E, W) != scNull()) { shrink(E, W, u, w); didShrink = true; break; }
            else { augmentFrom(E, W, u, w); didAugment = true; augmentedAny = true; break; }
          } else if (mNLabel[W] == LBL_ODD) {
            continue;
          } else {
            Size bw = mBBase[W];
            Size m = mMate[bw];
            if (m == scNull()) { augmentFrom(E, W, u, w); didAugment = true; augmentedAny = true; break; }
            Size M = outerOf(m);
            // GROW only if BOTH W and its matched partner M are currently OUTSIDE
            // the tree. In the PERSISTENT forest, M may already be labeled (e.g.
            // M==E, or M is W's ancestor) -- blindly re-labeling M then created a
            // PARENT CYCLE (0<->blossom20) and lcaOuter spun forever. If M is
            // already in the structure, this edge yields no clean grow; skip it.
            if (mNLabel[M] != LBL_FREE_OF_TREE) continue;
            if (W == E || M == W || M == E) continue;  // never self/2-cycle the forest
            mNLabel[W] = LBL_ODD; mNParent[W] = E;
            mNParEdgeU[W] = w; mNParEdgeV[W] = u;
            mNRoot[W] = mNRoot[E];
            mNLabel[M] = LBL_EVEN; mNParent[M] = W;
            mNParEdgeU[M] = m; mNParEdgeV[M] = bw;
            mNRoot[M] = mNRoot[E];
            wq.push_back(M);
          }
        }
      }
    }
    if (didShrink) continue;                 // labels preserved -> re-grow, no relabel
    if (didAugment) {                        // tree consumed -> relabel for more paths
      if (crit != CRIT2 && crit != CRIT3) return true;  // CRIT1: one augment per call
      // CRIT2/3 (finalization) keep going to perfection within this call.
      needRelabel = true; continue;
    }

    // ---- DUAL ADJUSTMENT: delta-jump (incremental, labels persist) ----
    LL unmTarget = (crit == CRIT1) ? 0 : -2;
    LL best = -1;
    for (Size u = 0; u < mN; ++u) {
      Size U = outerOf(u);
      if (mNLabel[U] != LBL_EVEN) continue;
      for (Size k = 0; k < mAdj[u].size(); ++k) {
        Size w = mAdj[u][k]; Size W = outerOf(w);
        if (W == U) continue;
        LL s = slackZ(u, w);
        int lw = mNLabel[W];
        if (lw == LBL_EVEN) { LL d=(s-unmTarget)/2; if (d>=0 && (best<0||d<best)) best=d; }
        else if (lw == LBL_FREE_OF_TREE) { LL d=s-unmTarget; if (d>=0 && (best<0||d<best)) best=d; }
      }
    }
    for (Size b = mN; b < mNumNodes; ++b)
      if (mBParent[b] == scNull() && mNLabel[b] == LBL_ODD) {
        LL d = mZ[b] / 2; if (d>=0 && (best<0||d<best)) best=d;
      }
    if (best < 0) return augmentedAny;       // no event -> can't progress

    // pre-dissolve zero-z ODD blossoms (z-clamp fix: don't push z below 0).
    bool dissolvedAny = false;
    for (Size b = mN; b < mNumNodes; ++b)
      if (mBParent[b] == scNull() && mNLabel[b] == LBL_ODD && mZ[b] == 0)
        { dissolve(b); dissolvedAny = true; }
    if (dissolvedAny) { needRelabel = true; continue; }

    LL delta = (best == 0 ? 1 : best);       // DELTA-JUMP (>=1); best>=1 unless a slack-0
                                             // edge slipped the grow -> unit step is safe.
    // apply incrementally to LABELED nodes only.
    for (Size v = 0; v < mN; ++v) {
      Size V = outerOf(v);
      if (mNLabel[V] == LBL_EVEN) mY[v] -= delta;
      else if (mNLabel[V] == LBL_ODD) mY[v] += delta;
    }
    for (Size b = mN; b < mNumNodes; ++b)
      if (mBParent[b] == scNull()) {
        if (mNLabel[b] == LBL_EVEN) mZ[b] += 2 * delta;
        else if (mNLabel[b] == LBL_ODD) mZ[b] -= 2 * delta;
      }
    mTNow += delta;
    // dissolve any ODD blossom whose z hit 0.
    bool dis2 = false;
    for (Size b = mN; b < mNumNodes; ++b)
      if (mBParent[b] == scNull() && mNLabel[b] == LBL_ODD && mZ[b] == 0)
        { dissolve(b); dis2 = true; }
    if (dis2) needRelabel = true;
    guardAdj += delta;
    if (guardAdj > guardCap) return augmentedAny;  // safety
    // loop: re-grow with new duals (labels persist unless needRelabel set above).
  }
}

inline bool ScalingMatching::runPerfectS5() {
  mNoPerfect = false;
  if (mN == 0) return true;
  if (mN % 2 != 0) { mNoPerfect = true; return false; }

  Size cap = 2 * mN + 1;
  mBParent.assign(cap, scNull());
  mBBase.assign(cap, scNull());
  mBCycle.assign(cap, std::vector<Size>());
  mBCycEdgeU.assign(cap, std::vector<Size>());
  mBCycEdgeV.assign(cap, std::vector<Size>());
  mNodeLeaves.assign(cap, std::vector<Size>());
  mZ.assign(cap, 0);
  mY.assign(mN, mMaxW);
  for (Size v = 0; v < mN; ++v) { mBBase[v] = v; mNodeLeaves[v].assign(1, v); }
#ifdef PS5_LEAF_SPLICE
  scalingLfInit(mN, cap);
#endif
  mMate.assign(mN, scNull());
  mNumNodes = mN;
  mNLabel.assign(cap, LBL_FREE_OF_TREE);
  mNParent.assign(cap, scNull());
  mNParEdgeU.assign(cap, scNull());
  mNParEdgeV.assign(cap, scNull());
  mNRoot.assign(cap, scNull());
  mTNow = 0;

  for (Size guard = 0; guard <= 2 * mN + 2; ++guard) {
    bool perfect = true;
    for (Size v = 0; v < mN; ++v) if (mMate[v] == scNull()) { perfect = false; break; }
    if (perfect) return true;
    if (!edmondsSearchS5(CRIT1)) { mNoPerfect = true; return false; }
  }
  mNoPerfect = true;
  return false;
}

} // namespace Matchbox

#endif // SCALING_MATCHING_H
