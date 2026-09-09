#!/usr/bin/env python3
"""Deterministic Matrix Market generators for the three-way HPC campaign.

Topology and edge weights deliberately use independent pseudo-random streams.
Changing weight_seed or weight_max therefore cannot change the edge sequence.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import random
import re
import sys
import time
from array import array
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable


SCHEMA = "three-way-graph-v1"
ID_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")
PROGRESS_EDGES = 5_000_000


@dataclass(frozen=True)
class Spec:
    physical_row: int
    graph_id: str
    family: str
    n: int
    parameter: int
    topology_seed: int
    weight_seed: int
    weight_max: int
    require_connected: bool
    enabled: bool
    notes: str


class DSU:
    def __init__(self, n: int) -> None:
        self.parent = array("i", range(n))
        self.rank = bytearray(n)
        self.components = n

    def find(self, x: int) -> int:
        parent = self.parent
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(self, a: int, b: int) -> None:
        a = self.find(a)
        b = self.find(b)
        if a == b:
            return
        rank = self.rank
        if rank[a] < rank[b]:
            a, b = b, a
        self.parent[b] = a
        if rank[a] == rank[b]:
            rank[a] += 1
        self.components -= 1


def parse_bool(value: str) -> bool:
    normalized = value.strip().lower()
    if normalized in {"1", "true", "yes", "y"}:
        return True
    if normalized in {"0", "false", "no", "n"}:
        return False
    raise ValueError(f"invalid boolean: {value!r}")


def load_manifest(path: Path) -> list[Spec]:
    specs: list[Spec] = []
    seen: set[str] = set()
    with path.open("r", encoding="utf-8", newline="") as stream:
        for row_number, raw in enumerate(stream, start=1):
            if not raw.strip() or raw.lstrip().startswith("#"):
                continue
            fields = next(csv.reader([raw], delimiter="\t"))
            if len(fields) != 10:
                raise ValueError(
                    f"{path}:{row_number}: expected 10 tab-separated fields, got {len(fields)}"
                )
            spec = Spec(
                physical_row=row_number,
                graph_id=fields[0],
                family=fields[1],
                n=int(fields[2]),
                parameter=int(fields[3]),
                topology_seed=int(fields[4]),
                weight_seed=int(fields[5]),
                weight_max=int(fields[6]),
                require_connected=parse_bool(fields[7]),
                enabled=parse_bool(fields[8]),
                notes=fields[9],
            )
            validate_spec(spec)
            if spec.graph_id in seen:
                raise ValueError(f"duplicate graph_id: {spec.graph_id}")
            seen.add(spec.graph_id)
            specs.append(spec)
    return specs


def validate_spec(spec: Spec) -> None:
    if not ID_RE.fullmatch(spec.graph_id):
        raise ValueError(f"unsafe graph_id: {spec.graph_id!r}")
    if spec.family not in {"er", "ba", "clique_blocks", "biclique_blocks"}:
        raise ValueError(f"{spec.graph_id}: unsupported family {spec.family!r}")
    if spec.n <= 1 or spec.n % 2:
        raise ValueError(f"{spec.graph_id}: n must be positive and even")
    if spec.parameter <= 0:
        raise ValueError(f"{spec.graph_id}: parameter must be positive")
    if spec.weight_max <= 0:
        raise ValueError(f"{spec.graph_id}: weight_max must be positive")
    if spec.family == "er" and spec.parameter >= spec.n:
        raise ValueError(f"{spec.graph_id}: ER expected degree must be below n")
    if spec.family == "ba" and spec.parameter >= spec.n:
        raise ValueError(f"{spec.graph_id}: BA m must be below n")
    if spec.family == "clique_blocks":
        if spec.parameter < 2 or spec.parameter % 2 or spec.n % spec.parameter:
            raise ValueError(
                f"{spec.graph_id}: clique block size must be even and divide n"
            )
    if spec.family == "biclique_blocks":
        if spec.n % (2 * spec.parameter):
            raise ValueError(
                f"{spec.graph_id}: 2 * biclique side size must divide n"
            )


class Writer:
    def __init__(self, path: Path, n: int, weight_rng: random.Random, weight_max: int,
                 dsu: DSU | None) -> None:
        self.path = path
        self.n = n
        self.weight_rng = weight_rng
        self.weight_max = weight_max
        self.dsu = dsu
        self.edges = 0
        self.min_weight = weight_max
        self.max_weight = 1
        self.stream = path.open("w", encoding="ascii", newline="\n", buffering=1024 * 1024)
        self.stream.write("%%MatrixMarket matrix coordinate integer symmetric\n")
        self.header_pos = self.stream.tell()
        self.stream.write(f"{n} {n} {0:20d}\n")

    def edge(self, u: int, v: int) -> None:
        if u == v:
            raise AssertionError("self-loop generated")
        if u > v:
            u, v = v, u
        weight = self.weight_rng.randint(1, self.weight_max)
        self.stream.write(f"{u + 1} {v + 1} {weight}\n")
        self.edges += 1
        self.min_weight = min(self.min_weight, weight)
        self.max_weight = max(self.max_weight, weight)
        if self.dsu is not None:
            self.dsu.union(u, v)
        if self.edges % PROGRESS_EDGES == 0:
            print(f"PROGRESS edges={self.edges}", file=sys.stderr, flush=True)

    def finish(self) -> None:
        self.stream.seek(self.header_pos)
        self.stream.write(f"{self.n} {self.n} {self.edges:20d}\n")
        self.stream.close()

    def abort(self) -> None:
        self.stream.close()
        self.path.unlink(missing_ok=True)


def generate_er(spec: Spec, path: Path, actual_seed: int) -> tuple[int, int, int, int]:
    topology_rng = random.Random(actual_seed)
    weight_rng = random.Random(spec.weight_seed)
    p = spec.parameter / (spec.n - 1)
    log_failure = math.log1p(-p)
    dsu = DSU(spec.n) if spec.require_connected else None
    writer = Writer(path, spec.n, weight_rng, spec.weight_max, dsu)
    try:
        for u in range(spec.n - 1):
            v = u + 1
            while v < spec.n:
                failures = int(math.log1p(-topology_rng.random()) / log_failure)
                v += failures
                if v >= spec.n:
                    break
                writer.edge(u, v)
                v += 1
        writer.finish()
    except BaseException:
        writer.abort()
        raise
    components = dsu.components if dsu is not None else -1
    return writer.edges, writer.min_weight, writer.max_weight, components


def generate_ba(spec: Spec, path: Path) -> tuple[int, int, int, int]:
    topology_rng = random.Random(spec.topology_seed)
    weight_rng = random.Random(spec.weight_seed)
    repeated = array("I")
    writer = Writer(path, spec.n, weight_rng, spec.weight_max, None)
    m = spec.parameter
    try:
        for u in range(m + 1):
            for v in range(u + 1, m + 1):
                writer.edge(u, v)
                repeated.append(u)
                repeated.append(v)
        for u in range(m + 1, spec.n):
            chosen: set[int] = set()
            while len(chosen) < m:
                chosen.add(repeated[topology_rng.randrange(len(repeated))])
            for v in sorted(chosen):
                writer.edge(v, u)
                repeated.append(v)
                repeated.append(u)
        writer.finish()
    except BaseException:
        writer.abort()
        raise
    return writer.edges, writer.min_weight, writer.max_weight, 1


def generate_clique_blocks(spec: Spec, path: Path) -> tuple[int, int, int, int]:
    writer = Writer(path, spec.n, random.Random(spec.weight_seed), spec.weight_max, None)
    k = spec.parameter
    try:
        for first in range(0, spec.n, k):
            for u in range(first, first + k):
                for v in range(u + 1, first + k):
                    writer.edge(u, v)
        writer.finish()
    except BaseException:
        writer.abort()
        raise
    return writer.edges, writer.min_weight, writer.max_weight, spec.n // k


def generate_biclique_blocks(spec: Spec, path: Path) -> tuple[int, int, int, int]:
    writer = Writer(path, spec.n, random.Random(spec.weight_seed), spec.weight_max, None)
    side = spec.parameter
    block = 2 * side
    try:
        for first in range(0, spec.n, block):
            for u in range(first, first + side):
                for v in range(first + side, first + block):
                    writer.edge(u, v)
        writer.finish()
    except BaseException:
        writer.abort()
        raise
    return writer.edges, writer.min_weight, writer.max_weight, spec.n // block


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def validate_matrix_market(path: Path, expected: Spec | None = None) -> dict[str, object]:
    digest = hashlib.sha256()
    count = 0
    minimum: int | None = None
    maximum: int | None = None
    with path.open("rb") as raw:
        first = raw.readline()
        digest.update(first)
        if first.rstrip() != b"%%MatrixMarket matrix coordinate integer symmetric":
            raise ValueError(f"{path}: unexpected Matrix Market banner")
        line = raw.readline()
        digest.update(line)
        while line.startswith(b"%"):
            line = raw.readline()
            digest.update(line)
        dimensions = line.split()
        if len(dimensions) != 3:
            raise ValueError(f"{path}: invalid dimension line")
        rows, cols, declared = map(int, dimensions)
        if rows != cols:
            raise ValueError(f"{path}: graph matrix is not square")
        for line_number, line in enumerate(raw, start=3):
            digest.update(line)
            fields = line.split()
            if len(fields) != 3:
                raise ValueError(f"{path}:{line_number}: expected u v weight")
            u, v, weight = map(int, fields)
            if not (1 <= u < v <= rows):
                raise ValueError(f"{path}:{line_number}: endpoints are not canonical")
            if weight <= 0:
                raise ValueError(f"{path}:{line_number}: nonpositive weight")
            if expected is not None and weight > expected.weight_max:
                raise ValueError(f"{path}:{line_number}: weight exceeds manifest maximum")
            count += 1
            minimum = weight if minimum is None else min(minimum, weight)
            maximum = weight if maximum is None else max(maximum, weight)
    if count != declared:
        raise ValueError(f"{path}: declared {declared} edges but read {count}")
    if expected is not None and rows != expected.n:
        raise ValueError(f"{path}: expected n={expected.n}, read n={rows}")
    return {
        "n": rows,
        "m": count,
        "min_weight": minimum,
        "max_weight": maximum,
        "sha256": digest.hexdigest(),
        "bytes": path.stat().st_size,
    }


def generate_one(spec: Spec, output_dir: Path, overwrite: bool, resume: bool) -> dict[str, object]:
    graph_path = output_dir / f"{spec.graph_id}.mtx"
    metadata_path = output_dir / f"{spec.graph_id}.json"
    if graph_path.exists() and not overwrite:
        if not resume:
            raise FileExistsError(f"{graph_path} exists; use --resume or --overwrite")
        validation = validate_matrix_market(graph_path, spec)
        if not metadata_path.is_file():
            raise FileNotFoundError(
                f"{metadata_path} is missing; inspect the existing graph and use --overwrite explicitly"
            )
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
        if metadata.get("sha256") != validation["sha256"]:
            raise ValueError(f"metadata SHA-256 mismatch: {graph_path}")
        print(f"VALID_EXISTING id={spec.graph_id} sha256={validation['sha256']}")
        return metadata

    partial = output_dir / f".{spec.graph_id}.partial.{os.getpid()}"
    partial.unlink(missing_ok=True)
    started = time.monotonic()
    actual_topology_seed = spec.topology_seed
    if spec.family == "er":
        for attempt in range(8):
            actual_topology_seed = spec.topology_seed + attempt
            result = generate_er(spec, partial, actual_topology_seed)
            if not spec.require_connected or result[3] == 1:
                break
            print(
                f"REJECT_DISCONNECTED id={spec.graph_id} attempt={attempt + 1} "
                f"seed={actual_topology_seed} components={result[3]}",
                file=sys.stderr,
            )
            partial.unlink(missing_ok=True)
        else:
            raise RuntimeError(f"{spec.graph_id}: no connected ER realization in 8 attempts")
    elif spec.family == "ba":
        result = generate_ba(spec, partial)
    elif spec.family == "clique_blocks":
        result = generate_clique_blocks(spec, partial)
    else:
        result = generate_biclique_blocks(spec, partial)

    validation = validate_matrix_market(partial, spec)
    if validation["m"] != result[0]:
        raise AssertionError("post-generation edge count changed")
    os.replace(partial, graph_path)
    elapsed = time.monotonic() - started
    metadata: dict[str, object] = {
        "schema": SCHEMA,
        **asdict(spec),
        "actual_topology_seed": actual_topology_seed,
        "components": result[3],
        "connected": result[3] == 1,
        "average_degree": 2.0 * int(validation["m"]) / spec.n,
        "duplicate_policy": "unique_by-construction",
        "self_loop_policy": "none-generated",
        "edge_order": "generator-canonical",
        "generation_seconds": elapsed,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "generator_sha256": sha256_file(Path(__file__)),
        **validation,
    }
    temporary_metadata = metadata_path.with_suffix(".json.partial")
    temporary_metadata.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(temporary_metadata, metadata_path)
    print(
        f"GENERATED id={spec.graph_id} family={spec.family} n={spec.n} "
        f"m={validation['m']} sha256={validation['sha256']} seconds={elapsed:.3f}"
    )
    return metadata


def selected_specs(specs: Iterable[Spec], row: int | None, graph_id: str | None) -> list[Spec]:
    if row is not None:
        result = [spec for spec in specs if spec.physical_row == row]
    elif graph_id is not None:
        result = [spec for spec in specs if spec.graph_id == graph_id]
    else:
        result = [spec for spec in specs if spec.enabled]
    if not result:
        raise ValueError("selection matched no manifest row")
    return result


def expected_edges(spec: Spec) -> float | int:
    if spec.family == "er":
        return spec.n * spec.parameter / 2.0
    if spec.family == "ba":
        m = spec.parameter
        return m * (m + 1) // 2 + (spec.n - m - 1) * m
    if spec.family == "clique_blocks":
        return spec.n * (spec.parameter - 1) // 2
    return spec.n * spec.parameter // 2


def self_test() -> None:
    import tempfile

    with tempfile.TemporaryDirectory(prefix="three_way_graph_test_") as temporary:
        root = Path(temporary)
        specs = [
            Spec(1, "er", "er", 100, 12, 11, 21, 1000, True, True, ""),
            Spec(2, "ba", "ba", 100, 4, 12, 22, 1000, True, True, ""),
            Spec(3, "clique", "clique_blocks", 96, 16, 0, 23, 1000, False, True, ""),
            Spec(4, "biclique", "biclique_blocks", 96, 8, 0, 24, 1000, False, True, ""),
        ]
        for spec in specs:
            validate_spec(spec)
            generate_one(spec, root, False, False)
        # Weight changes must not alter ER topology.
        a = Spec(5, "er_w10", "er", 100, 12, 91, 31, 10, True, True, "")
        b = Spec(6, "er_w1m", "er", 100, 12, 91, 32, 1_000_000, True, True, "")
        generate_one(a, root, False, False)
        generate_one(b, root, False, False)
        def endpoints(path: Path) -> list[tuple[str, str]]:
            lines = path.read_text(encoding="ascii").splitlines()[2:]
            return [(line.split()[0], line.split()[1]) for line in lines]
        if endpoints(root / "er_w10.mtx") != endpoints(root / "er_w1m.mtx"):
            raise AssertionError("weight configuration changed ER topology")
        clone = Spec(7, "er_clone", "er", 100, 12, 91, 31, 10, True, True, "")
        generate_one(clone, root, False, False)
        if (root / "er_w10.mtx").read_bytes() != (root / "er_clone.mtx").read_bytes():
            raise AssertionError("identical seeds did not reproduce identical bytes")
    print("SELF_TEST_OK")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--row", type=int, help="physical TSV line number (for SLURM arrays)")
    parser.add_argument("--graph-id")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--resume", action="store_true", help="validate and retain an existing graph")
    parser.add_argument("--validate-only", action="store_true")
    parser.add_argument("--list", action="store_true", help="list manifest rows without generating")
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test()
        return
    if args.manifest is None or args.output_dir is None:
        parser.error("--manifest and --output-dir are required")
    if args.row is not None and args.graph_id is not None:
        parser.error("--row and --graph-id are mutually exclusive")
    manifest_specs = load_manifest(args.manifest)
    if args.list and args.row is None and args.graph_id is None:
        specs = manifest_specs
    else:
        specs = selected_specs(manifest_specs, args.row, args.graph_id)
    if args.list:
        for spec in specs:
            print(
                f"row={spec.physical_row} enabled={int(spec.enabled)} id={spec.graph_id} "
                f"family={spec.family} n={spec.n} expected_m={expected_edges(spec):.0f}"
            )
        return
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for spec in specs:
        path = args.output_dir / f"{spec.graph_id}.mtx"
        if args.validate_only:
            result = validate_matrix_market(path, spec)
            print(f"VALID id={spec.graph_id} n={result['n']} m={result['m']} sha256={result['sha256']}")
        else:
            generate_one(spec, args.output_dir, args.overwrite, args.resume)


if __name__ == "__main__":
    main()
