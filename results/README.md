# Verified November 2026 results

Compact derived records from timing job **2370439**, cross-checked against the
complete raw result archive on 2026-11-13. The 10.04 GB archive, full matching
exports and large graph files are deliberately not committed here.

- [All graph tables](TABLES.md): medians, paired ratios and peak memory.
- [`campaign-2370439.json`](campaign-2370439.json): all 807 per-trial rows,
  instance/configuration metadata, reference objectives and diagnostic counts.
- [`check_results.py`](check_results.py): coverage/consistency checks and
  deterministic table generation, using Python's standard library.
- [`export_from_summary.py`](export_from_summary.py): export from the original
  collector summary and passed audit receipts, which are preserved separately.

```bash
python3 results/check_results.py
python3 results/check_results.py --write-tables
```

These verify the small distributed dataset, not a new solver run or a new
check of all graph edges. Exact cluster material remains in the frozen
snapshot identified in [PROVENANCE.md](../PROVENANCE.md).

## Inventory

| Panel | Instances | Clean solves | Scope |
|---|---:|---:|---|
| Original Suite A | 25 | 186 | MWM; six designated instances retain six trials |
| Eligible Suite B | 50 | 450 | Identical inputs with independently equal MWM/MWPM optima |
| BA different-objective controls | 3 | 27 | Excluded from same-objective Scaling speedup rankings |
| Targeted structural additions | 24 | 144 | 12 construction variants, two topology seeds, E/G only |
| Total | 102 | 807 | 324 pair/triple timing tasks |

Every other graph/program pair has three trials. Separate instrumented runs
provide 102 classical diagnostic pairs; 24 new LEMON references and 12
topology-seed nonisomorphism checks passed. The additions use roughly one
million vertices, topology seeds 41001/41002, fixed weight seed 51001, and the
submitted per-edge hash law with weights in 1..1000. Topology seeds change
topology, not merely labels or weights.

## Findings and limits

- Edmonds wins 24/25 original Suite A instances. Gabow's narrow advantage on
  `spal_004` repeats in every paired trial.
- Scaling wins 36/50 eligible Suite B instances; Edmonds wins 14. Scaling was
  not run on the targeted MWM-only additions.
- Gabow wins all 24 additions in every paired trial, with median time
  reductions of about 4.1-53.1% and peak RSS about 1.49-1.90 times Edmonds.

The added structures were selected after exploratory work to probe when
candidate maintenance pays off. They are not an unbiased population sample,
12 unrelated invariant graph families, or a generic Gabow win rate. Keep their
aggregates separate from the original panel; no original losses are omitted.

`solver_s` excludes loading/validation; `process_s` does not. Memory is
whole-process peak RSS in KiB, including loading. Tables use the largest peak
across trials divided by 2^20 for GiB. G/E is a ratio of medians; paired ranges
describe these observations, not population confidence intervals.

Node-associated variation remains on some large Suite A inputs and affects
the programs differently. All observations are retained. Node groups are
stable pseudonyms, not causal diagnoses or different processor classes.

Classical mate hashes agree across engines and trials. Logical-event
checkpoint fingerprints agree too, but are probabilistic rather than a direct
complete sequence comparison. Old-build phase profiles must not be used as a
direct decomposition of new-build runtimes. These diagnostics record counts
and traces, not a new phase-timing campaign.

The comparison studies these research implementations, not a fastest-solver
leaderboard. LEMON supplies optimum references, not a timed competitor.
