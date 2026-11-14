"""Verify the frozen HPC baseline plus the one declared distribution repair."""
import hashlib
import json
from pathlib import Path


REPAIR_PATH = "src/classic/MatchingEngine.h"
REPAIR_ID = "empty-csr-row-pointer-v1"
BEFORE = b"&mData[mOff[u]]"
AFTER = b"mData.data() + mOff[u]"


def digest(raw):
    return hashlib.sha256(raw).hexdigest()


def repair_map(manifest, document):
    if document.get("schema") != "mwm-post-campaign-repairs-v1":
        raise ValueError("Unknown repair schema")
    repairs = document.get("repairs")
    if not isinstance(repairs, list) or len(repairs) != 1:
        raise ValueError("Exactly one declared repair is supported")
    repair = repairs[0]
    if repair.get("id") != REPAIR_ID or repair.get("path") != REPAIR_PATH:
        raise ValueError("Undeclared repair")
    expected = {row["path"]: row["sha256"] for row in manifest["files"]}
    if len(expected) != len(manifest["files"]):
        raise ValueError("Duplicate baseline path")
    if repair.get("baseline_sha256") != expected.get(REPAIR_PATH):
        raise ValueError("Wrong repair baseline hash")
    if repair.get("before") != BEFORE.decode() or repair.get("after") != AFTER.decode() or repair.get("count") != 2:
        raise ValueError("Unsupported repair recipe")
    return {REPAIR_PATH: repair}


def verify_bytes(raw, baseline_sha256, repair=None):
    if repair is None:
        if digest(raw) != baseline_sha256:
            raise ValueError("Undeclared source change")
        return
    if digest(raw) != repair["repaired_sha256"]:
        raise ValueError("Repaired source hash mismatch")
    if raw.count(AFTER) != 2 or BEFORE in raw:
        raise ValueError("Repair must contain exactly the two declared replacements")
    if digest(raw.replace(AFTER, BEFORE)) != baseline_sha256:
        raise ValueError("Reversing the repair does not restore the frozen baseline")


def verify(root):
    root = Path(root)
    manifest = json.loads((root / "tests/kernel_source_manifest.json").read_text())
    repairs = repair_map(manifest, json.loads((root / "tests/post_campaign_repairs.json").read_text()))
    for row in manifest["files"]:
        path = (root / row["path"]).resolve()
        if not path.is_relative_to(root.resolve()) or not path.is_file():
            raise ValueError("Missing/outside source: " + row["path"])
        try:
            verify_bytes(path.read_bytes(), row["sha256"], repairs.get(row["path"]))
        except ValueError as error:
            raise ValueError(row["path"] + ": " + str(error)) from error
    return len(manifest["files"]) - len(repairs), len(repairs)


def main():
    exact, repaired = verify(Path(__file__).resolve().parents[1])
    print(str(exact) + " frozen HPC source files match byte-for-byte; " + str(repaired) +
          " declared post-campaign repair reverses to the frozen baseline")


if __name__ == "__main__":
    main()
