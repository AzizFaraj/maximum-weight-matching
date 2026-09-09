#!/usr/bin/env python3
"""Collect per-graph JSON sidecars into a stable TSV ledger."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path


COLUMNS = [
    "schema", "graph_id", "family", "n", "m", "parameter",
    "topology_seed", "actual_topology_seed", "weight_seed", "weight_max",
    "min_weight", "max_weight", "average_degree", "components", "connected",
    "sha256", "bytes", "generation_seconds", "generated_at_utc", "notes",
]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--graph-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--verify-hash", action="store_true")
    args = parser.parse_args()

    rows = []
    for metadata_path in sorted(args.graph_dir.glob("*.json")):
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
        graph_path = args.graph_dir / f"{metadata['graph_id']}.mtx"
        if not graph_path.is_file():
            raise FileNotFoundError(graph_path)
        if graph_path.stat().st_size != metadata["bytes"]:
            raise ValueError(f"size mismatch: {graph_path}")
        if args.verify_hash and sha256_file(graph_path) != metadata["sha256"]:
            raise ValueError(f"SHA-256 mismatch: {graph_path}")
        rows.append(metadata)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=COLUMNS, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    print(f"LEDGER graphs={len(rows)} output={args.output}")


if __name__ == "__main__":
    main()
