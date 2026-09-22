# Source provenance

The exact Edmonds, Gabow, and Scaling implementation sources in this tree were
copied from the research/HPC source state based on commit:

```text
e3f3b6710393c25e2512812beec6aa1d024f17e5
```

The Scaling adapter permanently selects the practical paper configuration by
defining `SCALING_S5_FULL_KEY_DRIVER_DISABLE=1` before including
`ScalingMatching.h`.

This repository is derived from the GPL-3.0 Matchbox codebase at
<https://github.com/CSCsw/matchbox>. The local upstream source directory was
named `matchbox-10.12`. The exact upstream commit represented by that local
archive was not recorded, so this file does not claim one.

After this clean tree receives its first verified commit, that commit becomes
the authoritative identity for public builds.
