"""Validate the small published export and regenerate its readable tables."""
import argparse
import collections
import json
import math
from pathlib import Path
import statistics


def require(value, message):
    if not value:
        raise ValueError(message)


def analyze(path):
    data = json.loads(path.read_text(encoding="utf-8"))
    cases = {c["graph_id"]: c for c in data["instances"]}
    rows = data["timing_rows"]
    require(len(cases) == 102 and len(rows) == 807, "Wrong campaign counts")
    keys = {(r["graph_id"], r["engine"], r["trial"]) for r in rows}
    expected = {(c["graph_id"], engine, trial) for c in cases.values()
                for engine in c["engines"] for trial in range(1, c["trials"] + 1)}
    require(len(keys) == len(rows) and keys == expected, "Duplicate/missing run cells")
    groups = collections.defaultdict(lambda: collections.defaultdict(list))
    for r in rows:
        case = cases[r["graph_id"]]
        require(r["validated"] and r["status"] == "success", "Failed run")
        require(math.isfinite(r["solver_s"]) and r["solver_s"] > 0 and r["peak_rss_KiB"] > 0, "Bad metric")
        opt = case["expected_mwpm_weight"] if r["engine"] == "scaling" else case["expected_mwm_weight"]
        require(r["exact_weight"] == opt, "Independent optimum mismatch")
        groups[r["graph_id"]][r["engine"]].append(r)
    graphs = []
    for gid, engines in groups.items():
        medians = {e: statistics.median(r["solver_s"] for r in runs) for e, runs in engines.items()}
        memory = {e: max(r["peak_rss_KiB"] for r in runs) / 2**20 for e, runs in engines.items()}
        paired = {e: {r["trial"]: r for r in runs} for e, runs in engines.items()}
        ratios = [paired["gabow"][t]["solver_s"] / paired["edmonds"][t]["solver_s"] for t in paired["edmonds"]]
        require(len({r["mate_sha256"] for e in ("edmonds", "gabow") for r in engines[e]}) == 1,
                "Classical mate-vector hashes differ")
        graphs.append(dict(case=cases[gid], times=medians, memory=memory,
            g_over_e=medians["gabow"] / medians["edmonds"], paired_min=min(ratios), paired_max=max(ratios)))
    eligible = [g for g in graphs if g["case"]["common_objective_comparison"]]
    winners = collections.Counter(min(g["times"], key=g["times"].get) for g in eligible)
    require(len(eligible) == 50 and winners == {"scaling": 36, "edmonds": 14}, "Ranking counts differ")
    additions = [g for g in graphs if g["case"]["suite"] == "N"]
    require(len(additions) == 24 and all(g["paired_max"] < 1 for g in additions), "Added-panel result differs")
    require(len(data["diagnostics"]) == 102 and all(d["reference_verified"] and d["full_classical_mates_equal"]
            and d["logical_comparison"]["status"] == "checkpoint-equal" for d in data["diagnostics"]), "Diagnostic summary differs")
    require(not any(s in path.read_text() for s in ("/home/", "C:\\Users", '"pid"', '"hostname"')),
            "Raw workstation/host identifiers leaked into public data")
    return data, graphs


def table(graphs):
    lines = ["# Verified campaign result tables", "",
        "Generated from `campaign-2370439.json` using `python3 results/check_results.py --write-tables`.", "",
        "Times are medians of solver-only seconds; RSS is the largest per-process peak across trials, in GiB. G/E is the ratio of medians, not a population confidence statement. Original and targeted panels are kept separate. BA different-objective controls are explicitly labelled.", ""]
    for suite, label in (("A", "Original Suite A"), ("B", "Original Suite B (including labelled controls)"), ("N", "Targeted per-edge structural additions")):
        lines += ["## " + label, "", "| Graph | E s | G s | G/E | Paired G/E range | S s | E GiB | G GiB | S GiB | Objective |",
            "|---|---:|---:|---:|---:|---:|---:|---:|---:|---|"]
        for g in graphs:
            c, t, m = g["case"], g["times"], g["memory"]
            if c["suite"] != suite:
                continue
            st = f"{t['scaling']:.3f}" if "scaling" in t else "-"
            sm = f"{m['scaling']:.3f}" if "scaling" in m else "-"
            lines.append(f"| {c['graph_id']} | {t['edmonds']:.3f} | {t['gabow']:.3f} | {g['g_over_e']:.4f} | {g['paired_min']:.4f}-{g['paired_max']:.4f} | {st} | {m['edmonds']:.3f} | {m['gabow']:.3f} | {sm} | {c['objective_class']} |")
        lines.append("")
    return "\n".join(lines).rstrip() + "\n"


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data", type=Path, default=Path(__file__).with_name("campaign-2370439.json"))
    ap.add_argument("--write-tables", action="store_true")
    args = ap.parse_args()
    data, graphs = analyze(args.data)
    if args.write_tables:
        args.data.with_name("TABLES.md").write_text(table(graphs), encoding="utf-8")
    print("RESULTS_CHECK_PASS instances=102 tasks=324 clean_runs=807 eligible=50 additions=24")


if __name__ == "__main__":
    main()
