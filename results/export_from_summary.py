"""Export compact, path-free campaign evidence; no solvers are executed."""
import argparse
import hashlib
import json
from pathlib import Path


def read(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--summary", required=True)
    ap.add_argument("--raw-audit", required=True)
    ap.add_argument("--metadata-audit", required=True)
    ap.add_argument("--output", default=str(Path(__file__).with_name("campaign-2370439.json")))
    args = ap.parse_args()
    d, audit, meta = map(read, (args.summary, args.raw_audit, args.metadata_audit))
    if not (audit["passed"] and meta["passed"] and audit["gzip_complete_crc_valid"]):
        raise ValueError("A complete passed raw archive audit is required")
    tasks = {r["task"]: r for r in audit["hardware_rows"]}
    nodes = {node: f"host-{i:02d}" for i, node in enumerate(sorted(meta["hostnames"]), 1)}
    oracles = {r["graph_id"]: r["record"] for r in d["oracles"]}
    cases = {}
    diagnostics = []
    for item in d["diagnostics"]:
        r, c = item["record"], item["record"]["case"]
        gid = item["graph_id"]
        fields = ("graph_id", "suite", "graph_kind", "n", "m", "digest", "seed", "W",
                  "objective_class", "trials", "expected_mwm_weight", "expected_mwpm_weight")
        case = {k: c[k] for k in fields}
        if c["suite"] == "N":
            case.update({k: c[k] for k in ("family", "topology_seed", "weight_seed", "target_n", "weight_model")})
            case["planned_specification"] = c["planned_specification"]
            case["expected_mwm_weight"] = int(oracles[gid]["oracle_result"]["exact-weight"])
            case["independent_reference"] = "LEMON-1.3.1-new-MWM"
        else:
            case["independent_reference"] = "archived-LEMON-MWM-and-where-applicable-MWPM"
        case["engines"] = c["engines"]
        case["common_objective_comparison"] = c["suite"] == "B" and c["objective_class"] == "eligible"
        cases[gid] = case
        diagnostics.append(dict(graph_id=gid, reference_verified=item["reference_verified"],
            full_classical_mates_equal=r["full_mates_equal"],
            logical_comparison=r["logical_comparison"],
            counters={run["engine"]: run["counters"] for run in r["runs"]}))
    rows = []
    for r in d["timing_rows"]:
        keys = ("task", "graph_id", "suite", "trial", "engine", "status", "validated",
                "solver_s", "process_s", "peak_rss_KiB", "exact_weight", "cardinality", "mate_sha256")
        row = {k: r[k] for k in keys}
        row["node_group"] = nodes[tasks[r["task"]]["hostname"]]
        rows.append(row)
    result = dict(schema="mwm-paper-results-v1", campaign_date="2026-11-04",
        completed_results_review_date="2026-11-13", timing_job=2370439,
        source_state=meta["expected_source_state"],
        snapshot_branch="codex/hpc-20261104-v3", snapshot_commit="2eeb7a4",
        source_summary_sha256=hashlib.sha256(Path(args.summary).read_bytes()).hexdigest(),
        protocol=dict(compiler="GCC 11.4.1", language="C++17", optimization="-O3 -march=x86-64-v3 -mtune=generic",
            cpu="AMD EPYC 9334", omp_threads=1, same_task_sequential_engines=True,
            process_cap_seconds=54000, allocation_GiB=240,
            solver_time_excludes_loading_and_validation=True,
            memory="Per-process whole-process getrusage peak RSS in KiB; divide by 2^20 for GiB",
            shared_classical_tie_order=True, edmonds_publication_ownership=True,
            scaling="practical Hybrid; full-key disabled; both dissolution frontiers enabled",
            node_group_scope="Stable pseudonyms within this dataset only; not processor classes or causal labels"),
        audit=dict(clean_runs=807, timing_tasks=324, diagnostic_pairs=102,
            new_independent_oracles=24, topology_seed_pairs=12,
            complete_archive_gzip_crc_valid=True, matching_exports_checked=1011,
            all_clean_runs_validated=d["all_clean_runs_validated"],
            all_topology_pairs_verified=d["all_topology_pairs_verified"],
            scope="Raw receipts, source identities and mate-file hashes checked. Original large input files and executable bytes are not redistributed or independently reread by this export."),
        qualifications=[
            "The 24 additions are a targeted structural panel selected after exploratory work, not an unbiased population sample.",
            "Three BA graphs are different-objective controls and must not enter same-objective Scaling rankings.",
            "The new structures compare only Edmonds and Gabow; no MWPM eligibility is claimed.",
            "Checkpoint fingerprints are probabilistic evidence, not a direct complete logical-event comparison.",
            "Node-associated Suite A variation is retained; no node cause or causal old/new code speedup is established.",
            "No old-build phase timing is a direct decomposition of these new-build solver times."],
        instances=list(cases.values()), timing_rows=rows, diagnostics=diagnostics)
    Path(args.output).write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print("EXPORTED", len(cases), "instances", len(rows), "runs", args.output)


if __name__ == "__main__":
    main()
