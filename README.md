# Maximum-Weight Matching

C++17 research implementations of [Edmonds][edmonds-paper], [Gabow][gabow-paper]
and [Scaling][scaling-paper] for weighted undirected graphs.

Edmonds and Gabow solve maximum-weight matching. Scaling solves maximum-weight
perfect matching: it implements the Hybrid algorithm of Ran Duan, Seth Pettie
and Hsin-Hao Su, using a practical search backend with the full-key driver disabled.

## Code

| Algorithm | Public function | Implementation |
|---|---|---|
| [Edmonds][edmonds-paper] | `Matchbox::Edmonds` | [MatchingEngine.h](src/classic/MatchingEngine.h) |
| [Gabow][gabow-paper] | `Matchbox::Gabow` | [MatchingEngine.h](src/classic/MatchingEngine.h) |
| [Scaling (Hybrid)][scaling-paper] | `Matchbox::Scaling` | [ScalingMatching.h](src/scaling/ScalingMatching.h) |

The [public API](include/MatchingAlgorithms.h) is implemented by a small
[adapter](src/MatchingAlgorithms.cc). Other headers provide required data structures.

## Build and test

Requires CMake 3.20+, a C++17 compiler and OpenMP. Python 3.9+ runs the
source and results checks. GNU/Linux is the primary supported platform.

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel 2
ctest --test-dir build --output-on-failure
```

## Run

Input: a weighted square Matrix Market graph.

```bash
build/bin/mwm_run edmonds graph.mtx
build/bin/mwm_run gabow graph.mtx
build/bin/mwm_run scaling graph.mtx
```

Scaling reports `no-perfect-matching` when no perfect matching exists.

## Results and reproducibility

[Results](results/README.md) contain tables and per-trial data for 807 timing
runs on 102 instances. Large graph files and raw matching exports are not included.

Reported timings use the experimental code preserved at tag
`hpc-experiments-20261104-v3`. The current distribution includes a documented
empty-row pointer repair made after those experiments. See
[PROVENANCE.md](PROVENANCE.md) for source identity, validation, build settings
and limitations. LEMON was an independent oracle, not a library dependency.

## License

Derived from [Matchbox](https://github.com/CSCsw/matchbox) under GPL-3.0.
See [LICENSE](LICENSE) and [NOTICE.md](NOTICE.md).

[edmonds-paper]: https://nvlpubs.nist.gov/nistpubs/jres/69B/jresv69Bn1-2p125_A1b.pdf
[gabow-paper]: https://arxiv.org/abs/1611.07541
[scaling-paper]: https://arxiv.org/abs/1411.1919
