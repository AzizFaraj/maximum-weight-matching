# Experimental Evaluation: Maximum Weighted Matching Algorithms

This repository contains the research software used to compare three exact
weighted-matching implementations on general undirected graphs:

- **Edmonds** — the repository's engineered GMG/Edmonds implementation;
- **Gabow** — the repository's Gabow 2018 backend; and
- **Scaling** — the Figure 6 hybrid scaling implementation.

The names **Edmonds**, **Gabow**, and **Scaling** are used consistently in the
software, experiment manifests, and accompanying paper. See
[`docs/ALGORITHM_IDENTITIES.md`](docs/ALGORITHM_IDENTITIES.md) for the exact
entry points and the important distinction between the practical and full-key
Scaling configurations.

> **Pre-release status:** this branch is being prepared while the final HPC
> campaigns are still running. It is suitable for local build and correctness
> work, but it must not be advertised as the final paper artifact until every
> item in [`RELEASE_CHECKLIST.md`](RELEASE_CHECKLIST.md) is closed.

## Supported environment

The currently tested platform is GNU/Linux or Cygwin with:

- a C++17 compiler;
- CMake 3.20 or newer; and
- OpenMP.

Native Windows/MinGW is not currently a supported target because the inherited
timing layer uses POSIX interfaces.

## Build

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

On a Cygwin installation whose `PATH` also exposes a Windows build of GNU
Make, force the native Cygwin tool explicitly:

```bash
cmake -S . -B build -G "Unix Makefiles" \
  -DCMAKE_MAKE_PROGRAM=/usr/bin/make \
  -DCMAKE_BUILD_TYPE=Release
```

Portable builds are the default. To reproduce the native CPU optimization used
for paper timing binaries:

```bash
cmake -S . -B build-paper \
  -DCMAKE_BUILD_TYPE=Release \
  -DMWM_NATIVE_OPTIMIZATION=ON
cmake --build build-paper --parallel
```

The benchmark executable is `build/bin/mwm_benchmark`.

## Quick check

```bash
build/bin/mwm_benchmark edmonds tests/fixtures/perfect_k4.mtx general file 1000 1 1 1
build/bin/mwm_benchmark gabow   tests/fixtures/perfect_k4.mtx general file 1000 1 1 1
build/bin/mwm_benchmark scaling tests/fixtures/perfect_k4.mtx general file 1000 1 1 1
```

Each command must report `status=success` and `exact-weight=19`.

## Experimental scope

Edmonds and Gabow solve unrestricted maximum-weight matching. Scaling solves
maximum-weight **perfect** matching and therefore has a narrower input domain.
The paper uses two suites rather than silently changing either objective:

- **Suite A** compares Edmonds and Gabow on unrestricted MWM instances.
- **Suite B** compares all three on the same stored-weight instances after an
  independent eligibility check establishes that the MWM and MWPM optimum
  objective values coincide.

No algorithm-specific weight offset, vertex deletion, or graph reduction is
used in the three-way runs. Input normalization is documented in
[`docs/INPUT_FORMAT.md`](docs/INPUT_FORMAT.md), and the evolving reproduction
record is in [`repro/README.md`](repro/README.md).

## Repository contents

- `include/`, `src/`: the minimal implementation closure used by the paper;
- `benchmarks/drivers/`: the instrumented three-way driver and optional LEMON
  oracle;
- `benchmarks/generators/`: deterministic Suite B graph generation;
- `tests/`: focused correctness, identity, and input-contract gates;
- `repro/`: experiment protocol and HPC release material; and
- `docs/`: precise algorithm and input identities.

Large graph datasets, raw multi-gigabyte result archives, binaries, profiling
outputs, and internal development backups are deliberately excluded.

## License and provenance

The software is distributed under the GNU General Public License version 3.0
only; see [`LICENSE`](LICENSE). This is a modified research fork of Matchbox.
Upstream revision and human-author metadata must be verified before public
release; see [`NOTICE.md`](NOTICE.md) and the release checklist. LEMON is an
optional external validation dependency and is not included in this repository.
