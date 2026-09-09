# Algorithm identities

The public names in this project are **Edmonds**, **Gabow**, and **Scaling**.
They have equal status in the experimental study, but they do not all solve the
same optimization problem on arbitrary inputs.

| Public name | Exact entry point | Backend/configuration | Objective |
|---|---|---|---|
| Edmonds | `MatchingEngine::ComputeMaxEdgWghtMatchingEdmondsFastV2()` | `gmg1986` | unrestricted MWM |
| Gabow | `MatchingEngine::ComputeMaxEdgWghtMatchingGabow()` | `gabow2018` | unrestricted MWM |
| Scaling | `ScalingMatching::computeHybrid()` | Figure 6 hybrid | MWPM |

## Practical Scaling used in headline timings

The benchmark target defines `SCALING_S5_FULL_KEY_DRIVER_DISABLE`. It retains
the Figure 6 hybrid entry point and integrated event search, while selecting
the faster practical frontier instead of the substantially slower full-key
Section 5 driver. Output records therefore identify it as
`configuration=practical` and `full-key-driver=disabled`.

This configuration must be called **practical Scaling** whenever the backend
distinction matters. Its timings must not be presented as establishing the
exact implementation bound of the strongest theoretical Section 5
realization.

The full-key implementation remains available for conformance and sensitivity
work by building a target without `SCALING_S5_FULL_KEY_DRIVER_DISABLE`; it is
not the headline timing configuration.

## Objective-safe comparison

On an arbitrary graph, a maximum-weight matching need not be perfect. Merely
removing a final perfectness check from Scaling does not convert its MWPM state
machine into an unrestricted-MWM algorithm.

Accordingly:

- Suite A reports Edmonds and Gabow on unrestricted MWM inputs.
- Suite B admits a graph to the three-way comparison only after an independent
  oracle shows that its MWM and MWPM optimum objective values are equal.
- All three Suite B runs consume the same graph bytes and stored weights.
