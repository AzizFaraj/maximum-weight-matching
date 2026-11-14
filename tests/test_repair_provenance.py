"""Negative checks for the narrowly allowlisted post-campaign repair."""
import copy
import json
from pathlib import Path
import tempfile
import unittest

import test_source_identity as identity

ROOT = Path(__file__).resolve().parents[1]


class RepairProvenanceTests(unittest.TestCase):
    def setUp(self):
        self.manifest = json.loads((ROOT / "tests/kernel_source_manifest.json").read_text())
        self.document = json.loads((ROOT / "tests/post_campaign_repairs.json").read_text())
        self.repair = self.document["repairs"][0]
        self.raw = (ROOT / self.repair["path"]).read_bytes()

    def test_complete_current_source_closure(self):
        self.assertEqual(identity.verify(ROOT), (34, 1))

    def test_extra_repaired_byte_rejected(self):
        with self.assertRaises(ValueError):
            identity.verify_bytes(self.raw + b"\n", self.repair["baseline_sha256"], self.repair)

    def test_extra_edit_rejected_even_if_repaired_hash_is_updated(self):
        raw = self.raw + b"\n"
        repair = dict(self.repair, repaired_sha256=identity.digest(raw))
        with self.assertRaisesRegex(ValueError, "restore the frozen"):
            identity.verify_bytes(raw, repair["baseline_sha256"], repair)

    def test_unrepaired_source_not_accepted_as_the_repaired_revision(self):
        raw = self.raw.replace(identity.AFTER, identity.BEFORE)
        with self.assertRaises(ValueError):
            identity.verify_bytes(raw, self.repair["baseline_sha256"], self.repair)

    def test_one_missing_replacement_rejected(self):
        raw = self.raw.replace(identity.AFTER, identity.BEFORE, 1)
        repair = dict(self.repair, repaired_sha256=identity.digest(raw))
        with self.assertRaisesRegex(ValueError, "exactly the two"):
            identity.verify_bytes(raw, repair["baseline_sha256"], repair)

    def test_undeclared_other_file_change_rejected(self):
        with self.assertRaises(ValueError):
            identity.verify_bytes(b"changed", identity.digest(b"original"))

    def test_unknown_duplicate_missing_and_invalid_repairs_rejected(self):
        for change in ("unknown", "duplicate", "missing", "baseline", "recipe", "schema"):
            with self.subTest(change=change):
                document = copy.deepcopy(self.document)
                if change == "unknown":
                    document["repairs"][0]["path"] = "src/scaling/ScalingMatching.h"
                elif change == "duplicate":
                    document["repairs"] *= 2
                elif change == "missing":
                    document["repairs"] = []
                elif change == "baseline":
                    document["repairs"][0]["baseline_sha256"] = "0" * 64
                elif change == "recipe":
                    document["repairs"][0]["count"] = 1
                else:
                    document["schema"] = "unknown"
                with self.assertRaises(ValueError):
                    identity.repair_map(self.manifest, document)

    def test_missing_source_rejected(self):
        with tempfile.TemporaryDirectory(prefix="repair-provenance-", dir=ROOT / "tests") as directory:
            root = Path(directory)
            (root / "tests").mkdir()
            (root / "tests/kernel_source_manifest.json").write_text(json.dumps(self.manifest))
            (root / "tests/post_campaign_repairs.json").write_text(json.dumps(self.document))
            with self.assertRaisesRegex(ValueError, "Missing/outside"):
                identity.verify(root)


if __name__ == "__main__":
    unittest.main(verbosity=2)
