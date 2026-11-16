# Source and experiment identity

## Verified HPC baseline

The implementation kernels come from the submitted package
`hpc_mwm_full_peredge_20261104_v3.tar.gz`, not the evolving research workspace:

```text
base commit: e3f3b6710393c25e2512812beec6aa1d024f17e5
classical source state:
e3f3b6710393c25e2512812beec6aa1d024f17e5+shared-membership-reuse-e33a6f4c4886
package SHA-256:
7ccb4a3546b01b54c328875524cae96ff09f7240261b129708311ff1e86d2864
```

The exact submission and preservation index are retained separately in the
local archival branch `codex/hpc-20261104-v3`, commit `2eeb7a4`, unchanged.
That archival branch is not yet published by this clean-distribution update;
its package includes original cluster identifiers and needs a separate
publication review. Its pending/running language describes when it was
created. Completed results are documented in this branch's
[results](results/README.md).

## Source mapping

The program named `Edmonds` follows the weighted blossom-search organization
of [Galil, Micali and Gabow (1986)](https://doi.org/10.1137/0215009), with the
engineering adaptations documented here and in the study.
[Edmonds (1965)](https://nvlpubs.nist.gov/nistpubs/jres/69B/jresv69Bn1-2p125_A1b.pdf)
is the foundational weighted-matching reference, not the implementation-specific link.

`src/classic/MatchingEngine.h`, `ConcatenableQueue.h`, `BlossomMerging.h` and
required candidate helpers come from the submitted `candidate/` overrides.
Other private classic/scaling supports come from `frozen/include/`.
`src/scaling/ScalingMatching.h` is the submitted practical Scaling kernel;
no new Scaling optimization was made for this distribution. The classical
header has the explicitly documented post-campaign repair described below.

Public graph/type/utility headers and `src/Graph.cc` retain submitted bytes.
The unchanged [manifest](tests/kernel_source_manifest.json) records the 35
experimental source paths and hashes. Of these, 34 still match byte-for-byte;
one has the allowlisted two-expression repair in
[post_campaign_repairs.json](tests/post_campaign_repairs.json).
`python3 tests/test_source_identity.py` requires its repaired hash and reverses
exactly those two replacements to recover the original byte hash. Other
changes are rejected. Targeted `.gitattributes` rules preserve historical
kernel line endings. The adapter, CLI and CMake layout are distribution code,
not the experiment frontend.

## Post-campaign empty-row repair (2026-11-14)

The exact pre-repair experimental distribution is preserved at annotated tag
`hpc-experiments-20261104-v3`, commit
`5bd1fabc0ceb9965f841596d99727da35fb63f5e`. Original HPC packages, campaign
results, stored weights, and timing tables were not changed by this repair.

The mutable and const classical `FlatJagged` accessors formed a row pointer
using `&mData[mOff[u]]`. For an empty trailing row backed by a nonempty vector,
the offset equals `mData.size()`, violating `vector::operator[]`'s bounds
contract. The repaired expressions use `mData.data() + mOff[u]`, retaining the
existing empty-storage guard. This permits the intended one-past pointer
without changing CSR layout, matching policy, weights, or search algorithms.
Bounds checks are not suppressed. A regression covers isolated vertices and
completely empty adjacency storage in both classical backends.

The reported campaign timings belong to the preserved experimental code,
not this repaired revision. Tiny-case checks are correctness/regression
evidence, not a new performance campaign or a claim that all possible inputs
are unaffected. The repaired code passed the Linux checks recorded below;
successful finite tests do not prove safety for every possible input.

Local repair validation passed all six CTest checks in both Debug (including
`_GLIBCXX_ASSERTIONS`) and Release with native Windows GCC 13.2. The run used
the existing empty MinGW `sys/times.h` compatibility include outside this
repository. This is local regression evidence, not Linux sanitizer evidence.
The original exact, randomized-brute-force, and 676-vertex tests were retained;
the new public regression adds ten empty-row fixtures for both classical
backends. Negative provenance checks reject undeclared or additional edits.

## Frozen behavior and build

The public build fixes shared classical D3 tie ordering, Edmonds publication
ownership, practical Scaling with the full-key driver disabled, and both
initial-dissolution and zero-root-cleanup frontiers. Hot profiling and
diagnostic output are disabled. No initial matching, learned duals or
algorithm-specific weight offset is supplied.

The measured executable used Linux, GCC 11.4.1, C++17/OpenMP,
`-O3 -DNDEBUG -march=x86-64-v3 -mtune=generic`, one thread and one inherited
CPU on exclusive EPYC 9334 nodes. Portable compilation is the default here.
`MWM_PAPER_CPU_TARGET=ON` selects the fixed instruction target;
`MWM_NATIVE_OPTIMIZATION=ON` is machine-specific and not equivalent.

Identical kernels do not make this library's CLI identical to the HPC
benchmark: input compatibility readers, validation, timing scopes and output
exports differ. The archived cluster package is the reproduction authority
for reported timings.

## Completed evidence

Jobs 2370437-2370441 cover setup, generation, timing, diagnostics and new
LEMON references. The result archive was received and reviewed on 2026-11-13.
Gzip integrity, 807 raw clean records, 1,011 matching-export hashes, 204
diagnostic checkpoint files, source manifests and 24 new oracle outputs
passed their recorded consistency checks. Raw matching files are not in Git.
Original large graphs and HPC executable bytes were absent from the archive;
the local review did not repeat every input-edge weight scan or independently
hash the actual executable bytes.

Derived results omit absolute workstation/HPC paths, process IDs and raw
hostnames. Stable node-group pseudonyms retain observed associations. This
does not imply that all hardware effects are controlled. Checkpoint agreement
is probabilistic, not a directly compared complete event stream.

## Distribution checks (2026-11-13)

Portable Release with Cygwin GCC 13.4 passes all three C++ tests, including
randomized brute-force checks and the 676-vertex membership regression.
Standalone Python verifies all 35 frozen source files and the exported
807-row dataset. The fixed x86-64-v3 Cygwin build compiles, but its randomized
and membership tests exit with signal 11. This is an unresolved distribution
check, not a finding that the separately validated Linux/GCC 11.4.1 HPC runs
failed. No kernel was changed to suppress it. The cause of the Cygwin failure
was not established and that platform-specific issue has not been reassessed
by the later Linux checks.

## Linux validation (2026-11-14)

GitHub Actions CI run #6 for repair commit `867b63f` (`28d7fa5` before the e-mail change below) completed successfully.
All six jobs passed: GCC Debug/Release, Clang Debug/Release, the fixed
x86-64-v3 paper CPU target, and AddressSanitizer/UndefinedBehaviorSanitizer.
The run's commit and status were confirmed from the supplied GitHub summary.
The original randomized test and the added empty-row regressions remain
enabled; no bounds or sanitizer checks were suppressed.

This validates that code revision's regression suite. Later documentation
or merge commits have their own CI runs and are not represented as already
tested merely because their kernel bytes match `867b63f`. None of these
checks reruns or replaces the paper's HPC experiments.

## Commit IDs (2026-11-16)

On 2026-11-16 the author, committer and tagger e-mail of every commit and of the annotated tag
was changed. No file changed, but every commit ID did:

| Before | After | Commit |
|---|---|---|
| `e21a48c` | `a81dfbb` | Initial research artifact |
| `7707c57` | `13a847a` | Create production-only solver layout |
| `369a296` | `ef148e4` | Publish clean Edmonds Gabow and Scaling implementation |
| `72fd9a7` | `5bd1fab` | Prepare verified HPC kernels and compact November results |
| `28d7fa5` | `867b63f` | Fix empty CSR row pointer construction without changing experimental snapshot |
| `4cee257` | `1547e44` | Promote verified artifact while preserving original main history |
| `eea846f` | `e2e604a` | Simplify README and record successful Linux validation |
| `6e695ee` | `353d3a4` | Link algorithm papers and identify Hybrid scaling |
| `e4c1d87` | `ebea56e` | Reference GMG for the engineered Edmonds implementation |
| `60bdb2d` | `0bfbfba` | Simplify the Scaling introduction in README |

## Attribution

Derived from the GPL-3.0 [CSCsw Matchbox project](https://github.com/CSCsw/matchbox).
The local upstream directory was named `matchbox-10.12`; its exact upstream
commit was not recorded and none is invented here. See LICENSE and NOTICE.md.
LEMON 1.3.1 was an independent oracle, not linked or redistributed by this
algorithm library.
