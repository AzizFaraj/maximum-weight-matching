# Source snapshot

The initial public-release preparation tree was assembled on 2026-09-09 from
the research branch based at:

```text
0513299ba14b6d7ada7f2330aab82ac90f3aa88f
```

The selected implementation files match the Suite B HPC package identified as:

```text
0513299ba14b6d7ada7f2330aab82ac90f3aa88f+tree-3fe93aad8224
```

That compound identity matters because `ScalingMatching.h` contained validated
but not-yet-committed changes relative to the base commit. Once the release
branch is committed and verified, the release commit itself becomes the
authoritative source identity. Future algorithm updates must be imported as
reviewable commits and rerun through the release tests.
