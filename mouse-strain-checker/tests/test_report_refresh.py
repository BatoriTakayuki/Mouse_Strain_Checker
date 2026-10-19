"""Interpretation refresh must preserve saved calls and all expensive inputs."""
from dataclasses import replace
import gzip
import hashlib
import json
from pathlib import Path
import tempfile
import unittest

from mouse_strain.core import Rules
from mouse_strain.panel import build_panel, sha256
from mouse_strain.pipeline import infer, write_json
from mouse_strain.reports import refresh_reports
from test_checker import fasta_fixture, vcf_fixture, pileup


class ReportRefreshTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.fasta = fasta_fixture(self.root)
        self.panel = self.root / "panel"
        build_panel(vcf_fixture(self.root), self.fasta, self.panel, "GRCm39", spacing=1)
        self.output = self.root / "result"
        self.rules = replace(Rules(), min_sites=2, min_window_sites=2, min_diagnostic=2, window_size=200)
        self.rows = infer(None, self.fasta, self.panel, self.output, self.rules,
                          pileup_lines=[pileup(10, "A" * 10), pileup(30, "A" * 10),
                                        pileup(210, "G" * 10), pileup(230, "G" * 10)])
        write_json(self.output / "complete.json", dict(accession="SRX1", bam={"path": "original.bam"}))

    def test_refresh_preserves_exact_calls_and_counts_and_keeps_backup(self):
        counts = sha256(self.output / "sites.tsv.gz")
        original_csv = (self.output / "chromosomes.csv").read_bytes()
        result = refresh_reports(self.panel, self.output)
        self.assertEqual(result["used_sites"], 4)
        self.assertEqual(sha256(self.output / "sites.tsv.gz"), counts)
        self.assertEqual((Path(result["backup"]) / "chromosomes.csv").read_bytes(), original_csv)
        diagnosis = json.loads((self.output / "diagnosis.json").read_text(encoding="utf-8"))
        self.assertEqual(diagnosis["summary"]["used_sites"], 4)
        completion = json.loads((self.output / "complete.json").read_text(encoding="utf-8"))
        self.assertEqual(completion["accession"], "SRX1")
        self.assertEqual(completion["bam"], {"path": "original.bam"})
        self.assertEqual(completion["report_schema"], 3)
        self.assertIn("background_reason", (self.output / "chromosomes.csv").read_text(encoding="utf-8"))
        self.assertFalse(list(self.output.glob(".report-refresh-*")))

    def test_panel_mismatch_leaves_report_and_marker_untouched(self):
        before = (self.output / "complete.json").read_bytes()
        evidence = json.loads((self.output / "evidence.json").read_text(encoding="utf-8"))
        evidence["panel"]["strains"] = ["wrong", "names"]
        write_json(self.output / "evidence.json", evidence)
        with self.assertRaisesRegex(ValueError, "Panel metadata differs"):
            refresh_reports(self.panel, self.output)
        self.assertEqual((self.output / "complete.json").read_bytes(), before)
        self.assertFalse((self.output / "reports_backup").exists())

    def test_intentionally_cleaned_bam_keeps_saved_sex_through_repeated_refresh(self):
        sex = dict(schema=1, status="male_likely", label="オス候補（RNA推定）",
                   reason="削除前のBAMでYマーカーを観測しました。", evidence=dict(y_fragments=123))
        write_json(self.output / "sex.json", sex)
        digest = hashlib.sha256((self.output / "sex.json").read_bytes()).hexdigest()
        write_json(self.output / "cleanup.json", dict(schema=1, scope="cache", sex_sha256=digest))
        for _ in range(2):
            result = refresh_reports(self.panel, self.output)
            self.assertEqual(result["sex"], sex)
            self.assertEqual(json.loads((self.output / "sex.json").read_text(encoding="utf-8")), sex)
    def test_changed_counts_cannot_silently_change_exact_calls(self):
        before = (self.output / "chromosomes.csv").read_bytes()
        with gzip.open(self.output / "sites.tsv.gz", "rt", encoding="utf-8") as handle:
            rows = handle.read()
        rows = rows.replace("10\t0\t0\t1", "0\t10\t0\t1")
        with gzip.open(self.output / "sites.tsv.gz", "wt", encoding="utf-8") as handle:
            handle.write(rows)
        with self.assertRaisesRegex(ValueError, "do not reproduce"):
            refresh_reports(self.panel, self.output)
        self.assertEqual((self.output / "chromosomes.csv").read_bytes(), before)
        self.assertFalse((self.output / "reports_backup").exists())
        self.assertFalse(list(self.output.glob(".report-refresh-*")))

    def test_original_exact_reason_and_regional_background_are_explained(self):
        self.assertEqual(self.rows[0]["status"], "regional_discordance")
        self.assertIn("一致しない", self.rows[0]["reason"])
        self.assertEqual(self.rows[0]["background_status"], "regional_discordance")
        self.assertEqual(self.rows[0]["background_call"], "")


if __name__ == "__main__":
    unittest.main()
