"""Sex evidence must use the original BAM and remain optional for strain reports."""
import json
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

from mouse_strain.core import Rules
from mouse_strain.pipeline import identity, render_html, sex_report, write_json
from mouse_strain.reports import refresh_reports
import test_report_refresh


class SexReportTests(unittest.TestCase):
    def test_collection_failure_preserves_an_explicit_unavailable_result(self):
        meta = {"assembly": "GRCm39", "contigs": {"X": 100}}
        with patch("mouse_strain.pipeline.collect_sex_evidence", side_effect=subprocess.TimeoutExpired("samtools", 60)):
            result = sex_report("sample.bam", meta, Rules())
        self.assertEqual(result["status"], "unavailable")
        self.assertIn("失敗", result["reason"])

    def test_html_contains_alternatives_and_inferred_sex_with_escaping(self):
        with tempfile.TemporaryDirectory() as directory:
            render_html(directory, [{"chromosome": "1", "call": "", "candidate_options": "A + B or A + C"}], [],
                        ["A", "B", "C"], {"label": "オス候補（推定）", "reason": "<script>",
                        "markers": [{"gene": "Uty", "fragments": 20}]},
                        {"stages": ["prefetch", "fastp"], "quality_notes": ["SRA Lite <品質>"]})
            report = (Path(directory) / "report.html").read_text(encoding="utf-8")
            self.assertIn("A + B or A + C", report)
            self.assertIn("オス候補（推定）", report)
            self.assertIn("Uty=20", report)
            self.assertNotIn("<script>", report)
            self.assertIn("&lt;script&gt;", report)
            self.assertIn("SRA Lite &lt;品質&gt;", report)

    def test_refresh_retains_lite_quality_note_and_processing_evidence(self):
        fixture = self._fixture()
        workflow = {"stages": ["prefetch", "fasterq-dump", "fastp"], "quality_notes": ["SRA Liteの品質値は簡略化されています。"]}
        completion = json.loads((fixture.output / "complete.json").read_text(encoding="utf-8"))
        completion["workflow"] = workflow
        write_json(fixture.output / "complete.json", completion)
        refresh_reports(fixture.panel, fixture.output)
        self.assertIn(workflow["quality_notes"][0], (fixture.output / "report.html").read_text(encoding="utf-8"))
        self.assertEqual(json.loads((fixture.output / "evidence.json").read_text(encoding="utf-8"))["workflow"], workflow)

    def _fixture(self):
        fixture = test_report_refresh.ReportRefreshTests(methodName="test_original_exact_reason_and_regional_background_are_explained")
        fixture.setUp()
        self.addCleanup(fixture.doCleanups)
        return fixture

    def test_missing_saved_bam_is_unavailable_and_old_sex_is_backed_up(self):
        fixture = self._fixture()
        write_json(fixture.output / "sex.json", {"label": "旧結果"})
        with patch("mouse_strain.reports.sex_report") as collector:
            result = refresh_reports(fixture.panel, fixture.output)
        collector.assert_not_called()
        self.assertEqual(result["sex"]["status"], "unavailable")
        self.assertEqual(json.loads((Path(result["backup"]) / "sex.json").read_text(encoding="utf-8"))["label"], "旧結果")

    def test_regional_hold_explains_that_options_pool_observed_sites(self):
        fixture = self._fixture()
        self.assertEqual(fixture.rows[0]["status"], "regional_discordance")
        self.assertIn("合算", fixture.rows[0]["candidate_options_reason"])
        self.assertEqual(fixture.rows[0]["call"], "")

    def test_modified_saved_bam_is_not_used_for_sex(self):
        fixture = self._fixture()
        bam = fixture.root / "original.bam"
        bam.write_bytes(b"original")
        write_json(fixture.output / "complete.json", {"accession": "SRX1", "bam": identity(bam)})
        bam.write_bytes(b"replaced-with-different-content")
        with patch("mouse_strain.reports.sex_report") as collector:
            result = refresh_reports(fixture.panel, fixture.output)
        collector.assert_not_called()
        self.assertIn("異なる", result["sex"]["reason"])

    def test_original_saved_bam_supplies_sex_without_changing_strain_calls(self):
        fixture = self._fixture()
        bam = fixture.root / "original.bam"
        bam.write_bytes(b"original")
        write_json(fixture.output / "complete.json", {"accession": "SRX1", "bam": identity(bam)})
        expected = {"schema": 1, "status": "inconclusive", "label": "判定保留", "reason": "不足", "markers": []}
        with patch("mouse_strain.reports.sex_report", return_value=expected) as collector:
            result = refresh_reports(fixture.panel, fixture.output)
        collector.assert_called_once()
        self.assertEqual(result["sex"], expected)
        self.assertEqual(json.loads((fixture.output / "sex.json").read_text(encoding="utf-8")), expected)


if __name__ == "__main__":
    unittest.main()
