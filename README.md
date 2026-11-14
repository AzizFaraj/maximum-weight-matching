# Maximum-weight matching: Edmonds, Gabow and Scaling

Three exact research implementations, with one small public API:

| Program | Public function | Implementation |
|---|---|---|
| Edmonds | `Matchbox::Edmonds` | [`src/classic/MatchingEngine.h`](src/classic/MatchingEngine.h), `ComputeMaxEdgWghtMatchingEdmondsFastV2` |
| Gabow | `Matchbox::Gabow` | [`src/classic/MatchingEngine.h`](src/classic/MatchingEngine.h), `ComputeMaxEdgWghtMatchingGabow` |
| Scaling | `Matchbox::Scaling` | [`src/scaling/ScalingMatching.h`](src/scaling/ScalingMatching.h), `computeHybrid` |

Edmonds and Gabow solve unrestricted maximum-weight matching (MWM). Scaling
solves maximum-weight **perfect** matching (MWPM) with practical Hybrid:
full-key driver disabled, with both dissolution frontiers enabled. No
alternative Scaling configuration is exposed by the public API.

The classical programs share blossom/search machinery and a tie-breaking
rule; their candidate-management backends differ. Supporting headers are
private dependencies, not additional public algorithms.
[`src/MatchingAlgorithms.cc`](src/MatchingAlgorithms.cc) is the adapter.

## Validated baseline and results

The experimental kernels, including the classical membership repairs, are
preserved at tag `hpc-experiments-20261104-v3` (commit `72fd9a7`). This branch
adds one post-campaign correctness repair: valid pointer construction for
empty classical CSR rows. The reported timings were **not** produced by this
repaired revision. See [PROVENANCE.md](PROVENANCE.md), the unchanged experimental
[source-identity manifest](tests/kernel_source_manifest.json), and the narrow
[repair record](tests/post_campaign_repairs.json).

[Results and scope](results/README.md) include all **807 timing runs** on
**102 instances**, with [readable tables](results/TABLES.md), per-trial data
and a lightweight Python check. Original and deliberately targeted structural
panels are kept separate. These are not universal speedup claims or an
external-solver ranking.

## Build and test

Requirements: CMake 3.20+, a C++17 compiler and OpenMP. GNU/Linux is the CI
target; Cygwin is also used for local verification. Python 3 runs the
source-identity/results checks, not the matching algorithms.

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel 2
ctest --test-dir build --output-on-failure
python3 results/check_results.py
```

Portable compilation is the default. To select the fixed campaign CPU target
(requires an x86-64-v3-capable CPU):

```bash
cmake -S . -B build-paper -DCMAKE_BUILD_TYPE=Release \
  -DMWM_PAPER_CPU_TARGET=ON
cmake --build build-paper --parallel 2
ctest --test-dir build-paper --output-on-failure
```

`MWM_NATIVE_OPTIMIZATION=ON` instead uses `-march=native`, which is
machine-specific and is **not** the campaign's fixed target. Do not combine
these options. This adapter/CLI is not the exact HPC benchmark executable;
frontend, compiler and validation differences still matter.

On Cygwin, if another Windows Make is on PATH, configure with
`-G "Unix Makefiles" -DCMAKE_MAKE_PROGRAM=/usr/bin/make`.

Local verification passes with the portable Cygwin build. The fixed-CPU-target
Cygwin build currently fails two tiny tests; see PROVENANCE.md. The Linux CI
fixed-target gate must pass before that distribution build is promoted.

## Run

```bash
build/bin/mwm_run edmonds graph.mtx
build/bin/mwm_run gabow graph.mtx
build/bin/mwm_run scaling graph.mtx
```

The CLI accepts weighted square Matrix Market graphs. Scaling reports
`no-perfect-matching` when its perfect-matching objective is infeasible;
Edmonds and Gabow do not require perfection. Exact integer-weight inputs
within the checked arithmetic range are the scope of the experiments.

For a fair three-way comparison, MWM and MWPM optimum weights must coincide.
Eligible Suite B inputs were screened independently with LEMON. The new
structural additions compare only Edmonds and Gabow; no MWPM eligibility is
assumed. No algorithm-specific weight shift, supplied initial matching or
supplied duals are used.

## Repository scope

This tree contains the three implementations, required data structures, a
CLI, small tests and compact results. It omits large datasets, raw matching
exports, development binaries, manuscript drafts and historical campaign
directories. Exact cluster scripts and generators remain in the separate
frozen snapshot identified in PROVENANCE.md.

These are engineered research implementations, not claims to realize every
original data-structure or asymptotic guarantee. See PROVENANCE.md for source
identity and evidence limits.

Derived from GPL-3.0 Matchbox; see [LICENSE](LICENSE) and [NOTICE.md](NOTICE.md).
LEMON is an external experimental oracle, not a library dependency.
