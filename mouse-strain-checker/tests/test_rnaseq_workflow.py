"""Public RNA-seq acquisition/preprocessing must feed the mapper and retain provenance."""
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch

from mouse_strain.__main__ import parser
from mouse_strain.panel import sha256
from mouse_strain.pipeline import analyze, tool_version, write_json


class RnaSeqWorkflowTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.fasta = self.root / "reference.fa"
        self.fasta.write_text(">1\nAAAA\n")
        self.meta = dict(assembly="GRCm39", fasta_sha256=sha256(self.fasta), contigs={"1": 4}, strains=["A", "B"])
        self.row = dict(run_accession="SRR1", sample_accession="SRS1", library_layout="SINGLE", fastq_ftp="", fastq_md5="")
        self.output = self.root / "result"
        self.events = []
        self.input_quality = None

    def acquire(self, row, directory, threads, method):
        self.events.append(("get_reads", method))
        directory.mkdir(parents=True)
        read = directory / "SRR1.fastq"
        read.write_text("@read\nACGT\n+\nIIII\n")
        acquisition = {"method": "sra", "schema": 3}
        if self.input_quality is not None:
            acquisition["quality"] = self.input_quality
        write_json(directory / "sra_reads.json", acquisition)
        return [read]

    def trim(self, paths, row, directory, threads):
        self.events.append(("fastp", paths[0].name))
        directory.mkdir(parents=True)
        read = directory / "SRR1.trim.fastq.gz"
        read.write_bytes(b"trimmed fixture")
        write_json(directory / "complete.json", {"provenance": {"method": "fastp"}, "summary": {"reads_after": 1}})
        return [read]

    def map(self, paths, row, prefix, index, bam, threads):
        self.events.append(("mapping", paths[0].name))
        bam.write_bytes(b"bam fixture")
        return bam

    def execute(self, **kwargs):
        with patch("mouse_strain.pipeline.load_panel", return_value=(self.meta, Mock())), \
             patch("mouse_strain.pipeline.index_identity", return_value={}), \
             patch("mouse_strain.pipeline.require_tools"), \
             patch("mouse_strain.pipeline.resolve_srx", return_value=[self.row]), \
             patch("mouse_strain.pipeline.get_reads", side_effect=self.acquire), \
             patch("mouse_strain.pipeline.trim_reads", side_effect=self.trim) as trimming, \
             patch("mouse_strain.pipeline.map_run", side_effect=self.map), \
             patch("mouse_strain.pipeline.infer", return_value=[]) as inference, \
             patch("mouse_strain.pipeline.shutil.which", return_value="/bin/tool"), \
             patch("mouse_strain.pipeline.tool_version", return_value="fixture version"):
            analyze("SRX1", self.fasta, "panel", "index", self.output, **kwargs)
        return trimming, inference, json.loads((self.output / "complete.json").read_text(encoding="utf-8"))

    def test_default_is_sra_then_fastp_then_trimmed_mapping(self):
        _, _, completion = self.execute()
        self.assertEqual(self.events, [("get_reads", "sra"), ("fastp", "SRR1.fastq"), ("mapping", "SRR1.trim.fastq.gz")])
        flow = json.loads((self.output / "workflow.json").read_text(encoding="utf-8"))
        self.assertEqual(completion["workflow"], flow)
        self.assertEqual(flow["name"], "sra_fastp_v1")
        self.assertEqual(flow["stages"][:4], ["prefetch", "vdb-validate", "fasterq-dump", "fastp"])
        self.assertEqual(flow["runs"][0]["preprocessing"]["summary"]["reads_after"], 1)
        self.assertEqual(flow["mapping"]["flags"], ["--dta", "-k", "5"])

    def test_explicit_ena_and_skip_fastp_are_recorded(self):
        trimming, _, completion = self.execute(download_method="ena", preprocess=False)
        trimming.assert_not_called()
        self.assertEqual(self.events, [("get_reads", "ena"), ("mapping", "SRR1.fastq")])
        flow = completion["workflow"]
        self.assertEqual(flow["name"], "ena_raw_v1")
        self.assertEqual(flow["stages"][0], "ENA HTTPS")
        self.assertFalse(flow["fastp_enabled"])

    def test_sra_lite_quality_remains_visible_in_completed_workflow_and_report_input(self):
        self.input_quality = {"format": "sra_lite", "base_quality": "simplified", "note": "SRA Liteの品質値は簡略化されています。"}
        _, inference, completion = self.execute()
        workflow = completion["workflow"]
        self.assertEqual(workflow["runs"][0]["acquisition"]["quality"], self.input_quality)
        self.assertEqual(workflow["quality_notes"], ["SRR1: " + self.input_quality["note"]])
        self.assertEqual(inference.call_args.kwargs["workflow"], workflow)

    def test_external_bam_does_not_claim_to_have_downloaded_or_trimmed(self):
        bam = self.root / "input.bam"
        bam.write_bytes(b"input fixture")
        trimming, _, completion = self.execute(bam=bam, bam_assembly="GRCm39")
        trimming.assert_not_called()
        self.assertEqual(self.events, [])
        self.assertEqual(completion["workflow"]["source"], "external_bam")
        self.assertIsNone(completion["workflow"]["download_method"])
        self.assertFalse(completion["workflow"]["fastp_enabled"])

    def test_cli_default_and_opt_out_are_explicit(self):
        args = parser().parse_args(["run", "--srx", "SRX1", "--out", "result"])
        self.assertEqual(args.download_method, "sra")
        self.assertFalse(args.skip_fastp)
        args = parser().parse_args(["run", "--srx", "SRX1", "--out", "result", "--download-method", "ena", "--skip-fastp"])
        self.assertEqual(args.download_method, "ena")
        self.assertTrue(args.skip_fastp)

    def test_version_ignores_toolkit_blank_lines_and_accepts_stderr(self):
        for stdout, stderr in [("\nprefetch : 3.4.1\n", ""), ("", "\nfastp 0.23.4\n")]:
            with self.subTest(stdout=stdout), patch("mouse_strain.pipeline.subprocess.run", return_value=SimpleNamespace(stdout=stdout, stderr=stderr)):
                self.assertEqual(tool_version("tool"), (stdout + stderr).strip())


if __name__ == "__main__":
    unittest.main()
