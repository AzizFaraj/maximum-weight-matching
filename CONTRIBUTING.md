# Contributing

This repository is a research artifact whose algorithm identities and timing
configurations must remain reproducible.

Before proposing a change:

1. build from a clean directory with CMake;
2. run `ctest --test-dir build --output-on-failure`;
3. state whether the change affects Edmonds, Gabow, Scaling, input
   normalization, result verification, or timing boundaries;
4. add a deterministic regression test for correctness changes; and
5. never relabel an experimental or sensitivity configuration as canonical.

Do not commit generated graphs, raw benchmark archives, binaries, profiler
outputs, or downloaded papers. Performance claims require the exact compiler,
flags, source revision, graph digest, seed, and trial protocol.

Large refactors of the matching engines should be separated from experimental
paper freezes because even cosmetic structural changes can affect correctness
and timing.
