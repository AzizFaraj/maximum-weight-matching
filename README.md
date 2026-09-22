# Experimental Evaluation: Maximum Weighted Matching Algorithms

This repository contains the three exact implementations evaluated in the
paper:

- **Edmonds**: unrestricted maximum-weight matching;
- **Gabow**: unrestricted maximum-weight matching; and
- **Scaling**: maximum-weight perfect matching using the practical Figure-6
  configuration used for the paper's headline experiments.

The public API contains only these three entry points. Supporting data
structures remain private so that users do not have to navigate development
variants or experiment-campaign scripts.

## Where the algorithms are

| Algorithm | Public entry point | Implementation |
|---|---|---|
| Edmonds | `Matchbox::Edmonds` | `src/classic/MatchingEngine.h` (`ComputeMaxEdgWghtMatchingEdmondsFastV2`) |
| Gabow | `Matchbox::Gabow` | `src/classic/MatchingEngine.h` (`ComputeMaxEdgWghtMatchingGabow`) |
| Scaling | `Matchbox::Scaling` | `src/scaling/ScalingMatching.h` (`computeHybrid`) |

Edmonds and Gabow intentionally share one source file because their optimized
implementations reuse the same primal-dual blossom machinery. Scaling has its
own private implementation directory. `src/MatchingAlgorithms.cc` is the thin
public adapter.

Scaling is compiled with `SCALING_S5_FULL_KEY_DRIVER_DISABLE=1`. This freezes
the practical configuration used in the experiments; the slower full-key
research path is not exposed as a selectable public algorithm.

## Build and test

Requirements: CMake 3.20+, a C++17 compiler, and OpenMP.

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

For the native CPU optimization used by the HPC timing builds:

```bash
cmake -S . -B build-paper \
  -DCMAKE_BUILD_TYPE=Release \
  -DMWM_NATIVE_OPTIMIZATION=ON
cmake --build build-paper --parallel
```

## Run

```bash
build/bin/mwm_run edmonds graph.mtx
build/bin/mwm_run gabow graph.mtx
build/bin/mwm_run scaling graph.mtx
```

Inputs are weighted, square Matrix Market graphs. Scaling returns
`no-perfect-matching` when its MWPM objective is not defined for the input.
Edmonds and Gabow do not require a perfect matching.

The paper therefore uses two experimental suites: Suite A compares Edmonds and
Gabow on unrestricted MWM instances; Suite B compares all three on identical
stored-weight graphs independently screened so that the MWM and MWPM optimum
weights coincide. No algorithm-specific weight shift or graph reduction is
used in a three-way timing run.

## Scope

This is the algorithm release, not the experiment-data archive. It deliberately
omits generated graphs, raw HPC results, paper sources, and campaign-specific
submission scripts. The code derives from the GPL-3.0 Matchbox project and is
distributed under the same license; see `LICENSE`, `NOTICE.md`, and
`PROVENANCE.md`.
