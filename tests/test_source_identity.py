"""Verify the byte-exact private source closure of the validated HPC kernel."""
import hashlib
import json
from pathlib import Path


def main():
    root = Path(__file__).resolve().parents[1]
    manifest = json.loads((root / "tests/kernel_source_manifest.json").read_text())
    failures = []
    for row in manifest["files"]:
        source = root / row["path"]
        if not source.is_file() or hashlib.sha256(source.read_bytes()).hexdigest() != row["sha256"]:
            failures.append(row["path"])
    if failures:
        raise SystemExit("Frozen source identity mismatch: " + ", ".join(failures))
    print(str(len(manifest["files"])) + " frozen HPC source files match byte-for-byte")


if __name__ == "__main__":
    main()
