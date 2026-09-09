// Smoke tests for FastMatchingEngineV2 in kGabow2018 backend mode.
// Mirrors tests/testMatchingEngineFast.cc but exercises the Gabow path.
// Compile from project root:
//   g++ -O2 -Wall -ansi -Iinclude tests/testGabowBackend.cc \
//       lib/optimized/libmatchbox.a -o tests/testGabow

#include <iostream>
#include <vector>
#include "MatchingEngine.h"

using namespace Matchbox;

#define CHECK(cond) do { if (!(cond)) { \
  std::cerr << "FAIL: " << #cond << " at line " << __LINE__ << std::endl; \
  return 1; } } while(0)

static int testEmpty() {
  Graph g;
  FastMatchingEngineV2 eng(FastMatchingEngineV2::kGabow2018);
  std::vector<Size> mate;
  Size card = 99;
  Val w = 99;
  eng.compute(g, &mate, &card, &w);
  CHECK(card == 0);
  CHECK(w == 0);
  std::cout << "testEmpty OK" << std::endl;
  return 0;
}

static int testSingleEdge() {
  // 2 vertices, 1 edge, weight 5.
  std::vector<Size> sv;  sv.push_back(0);
  std::vector<Size> tv;  tv.push_back(1);
  std::vector<Val>  wv;  wv.push_back(5.0);
  Graph g(2, 1, sv, tv, wv);
  Err initErr;
  g.GetInitErr(&initErr);
  CHECK(initErr == eErrNone);

  FastMatchingEngineV2 eng(FastMatchingEngineV2::kGabow2018);
  std::vector<Size> mate;
  Size card = 0;
  Val w = 0;
  eng.compute(g, &mate, &card, &w);

  std::cout << "  card=" << card << " weight=" << w
            << " mate[0]=" << (mate[0] == cNullItm ? -1 : (int)mate[0])
            << " mate[1]=" << (mate[1] == cNullItm ? -1 : (int)mate[1]) << std::endl;
  CHECK(card == 1);
  CHECK(w == 5.0);
  CHECK(mate[0] == 1);
  CHECK(mate[1] == 0);
  std::cout << "testSingleEdge OK" << std::endl;
  return 0;
}

static int testTwoDisjointEdges() {
  // 4 vertices, 2 disjoint edges (0-1) w=5, (2-3) w=3.
  std::vector<Size> sv;  sv.push_back(0); sv.push_back(2);
  std::vector<Size> tv;  tv.push_back(1); tv.push_back(3);
  std::vector<Val>  wv;  wv.push_back(5.0); wv.push_back(3.0);
  Graph g(4, 2, sv, tv, wv);

  FastMatchingEngineV2 eng(FastMatchingEngineV2::kGabow2018);
  std::vector<Size> mate;
  Size card = 0;
  Val w = 0;
  eng.compute(g, &mate, &card, &w);

  std::cout << "  card=" << card << " weight=" << w << std::endl;
  CHECK(card == 2);
  CHECK(w == 8.0);
  std::cout << "testTwoDisjointEdges OK" << std::endl;
  return 0;
}

static int testDisconnectedWithIsolates() {
  // Two components plus two isolated vertices.  The unmatched sentinels must
  // remain valid and must not affect the exact objective.
  std::vector<Size> sv; sv.push_back(0); sv.push_back(2); sv.push_back(3);
  std::vector<Size> tv; tv.push_back(1); tv.push_back(3); tv.push_back(4);
  std::vector<Val> wv; wv.push_back(9.0); wv.push_back(4.0); wv.push_back(7.0);
  Graph g(7, 3, sv, tv, wv);

  FastMatchingEngineV2 eng(FastMatchingEngineV2::kGabow2018);
  std::vector<Size> mate; Size card = 0; Val w = 0;
  eng.compute(g, &mate, &card, &w);
  CHECK(card == 2);
  CHECK(w == 16.0);
  CHECK(mate[0] == 1 && mate[1] == 0);
  CHECK(mate[3] == 4 && mate[4] == 3);
  CHECK(mate[5] == cNullItm && mate[6] == cNullItm);
  std::cout << "testDisconnectedWithIsolates OK" << std::endl;
  return 0;
}

static int testParallelEdges() {
  // Graph normalizes parallel input records to one simple edge carrying the
  // maximum weight.  Exercise both endpoint orientations and put the maximum
  // first, middle, and last across distinct pairs.
  std::vector<Size> sv, tv;
  std::vector<Val> wv;
  sv.push_back(0); tv.push_back(1); wv.push_back(11.0);
  sv.push_back(1); tv.push_back(0); wv.push_back(2.0);
  sv.push_back(0); tv.push_back(1); wv.push_back(7.0);
  sv.push_back(2); tv.push_back(3); wv.push_back(5.0);
  sv.push_back(3); tv.push_back(2); wv.push_back(13.0);
  sv.push_back(2); tv.push_back(3); wv.push_back(4.0);
  sv.push_back(4); tv.push_back(5); wv.push_back(1.0);
  sv.push_back(5); tv.push_back(4); wv.push_back(6.0);
  sv.push_back(4); tv.push_back(5); wv.push_back(9.0);
  Graph g(6, 9, sv, tv, wv);

  Size graphEdges = 0;
  g.GetNumEdgs(&graphEdges);
  CHECK(graphEdges == 3);
  Size degree = 0;
  std::vector<Size> adj;
  std::vector<Val> weights;
  CHECK(g.GetAdj(0, &degree, &adj, &weights) == eErrNone);
  CHECK(degree == 1 && adj[0] == 1 && weights[0] == 11.0);
  CHECK(g.GetAdj(3, &degree, &adj, &weights) == eErrNone);
  CHECK(degree == 1 && adj[0] == 2 && weights[0] == 13.0);
  CHECK(g.GetAdj(5, &degree, &adj, &weights) == eErrNone);
  CHECK(degree == 1 && adj[0] == 4 && weights[0] == 9.0);

  FastMatchingEngineV2 eng(FastMatchingEngineV2::kGabow2018);
  std::vector<Size> mate; Size card = 0; Val w = 0;
  eng.compute(g, &mate, &card, &w);
  CHECK(card == 3);
  CHECK(w == 33.0);
  CHECK(mate[0] == 1 && mate[1] == 0);
  CHECK(mate[2] == 3 && mate[3] == 2);
  CHECK(mate[4] == 5 && mate[5] == 4);
  std::cout << "testParallelEdges OK" << std::endl;
  return 0;
}

// 5-cycle with a chord forcing a blossom shrink-then-augment path.
// Vertices 0-1-2-3-4-0 form a 5-cycle, plus a chord (0,2).
// All edges weight 1. Max matching of a 5-cycle has cardinality 2 weight 2;
// the chord lets one extra edge live without changing card here.
static int testFiveCycleWithChord() {
  std::vector<Size> sv;
  std::vector<Size> tv;
  std::vector<Val>  wv;
  sv.push_back(0); tv.push_back(1); wv.push_back(1.0);
  sv.push_back(1); tv.push_back(2); wv.push_back(1.0);
  sv.push_back(2); tv.push_back(3); wv.push_back(1.0);
  sv.push_back(3); tv.push_back(4); wv.push_back(1.0);
  sv.push_back(4); tv.push_back(0); wv.push_back(1.0);
  sv.push_back(0); tv.push_back(2); wv.push_back(1.0);
  Graph g(5, 6, sv, tv, wv);

  FastMatchingEngineV2 eng(FastMatchingEngineV2::kGabow2018);
  std::vector<Size> mate;
  Size card = 0;
  Val w = 0;
  eng.compute(g, &mate, &card, &w);

  std::cout << "  card=" << card << " weight=" << w << std::endl;
  CHECK(card == 2);
  CHECK(w == 2.0);
  std::cout << "testFiveCycleWithChord OK" << std::endl;
  return 0;
}

// Same problem on two engines, default (GMG'86) vs Gabow: cards and weights
// must match.
static int testTriangleParity() {
  std::vector<Size> sv;
  std::vector<Size> tv;
  std::vector<Val>  wv;
  sv.push_back(0); tv.push_back(1); wv.push_back(3.0);
  sv.push_back(1); tv.push_back(2); wv.push_back(2.0);
  sv.push_back(2); tv.push_back(0); wv.push_back(4.0);
  Graph g(3, 3, sv, tv, wv);

  FastMatchingEngineV2 engGmg(FastMatchingEngineV2::kGMG1986);
  std::vector<Size> mGmg; Size cGmg = 0; Val wGmg = 0;
  engGmg.compute(g, &mGmg, &cGmg, &wGmg);

  FastMatchingEngineV2 engGab(FastMatchingEngineV2::kGabow2018);
  std::vector<Size> mGab; Size cGab = 0; Val wGab = 0;
  engGab.compute(g, &mGab, &cGab, &wGab);

  std::cout << "  gmg: card=" << cGmg << " weight=" << wGmg
            << "  gabow: card=" << cGab << " weight=" << wGab << std::endl;
  CHECK(cGmg == cGab);
  CHECK(wGmg == wGab);
  std::cout << "testTriangleParity OK" << std::endl;
  return 0;
}

int main() {
  if (testEmpty()) return 1;
  if (testSingleEdge()) return 1;
  if (testTwoDisjointEdges()) return 1;
  if (testDisconnectedWithIsolates()) return 1;
  if (testParallelEdges()) return 1;
  if (testFiveCycleWithChord()) return 1;
  if (testTriangleParity()) return 1;
  std::cout << "ALL TESTS PASSED" << std::endl;
  return 0;
}
