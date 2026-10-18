"""Regression coverage for metadata validation and interrupted pipeline stages."""
import http.client
import io
import json
from pathlib import Path
import subprocess
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch

from mouse_strain.pipeline import (
    analyze, build_index, download, get_reads, index_identity, infer, map_run,
    render_html, resolve_srx, write_csv, write_json,
)


class Response(io.BytesIO):
    def __init__(self, data, length=None):
        super().__init__(data)
        self.headers = {} if length is None else {"Content-Length": str(length)}


class PipelineReviewTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)

    def tearDown(self):
        self.temporary.cleanup()

    def ena_response(self, **changes):
        row = dict(run_accession="SRR1", experiment_accession="SRX1", sample_accession="SRS1",
                   tax_id="10090", scientific_name="Mus musculus", library_strategy="RNA-Seq",
                   library_source="TRANSCRIPTOMIC", library_layout="PAIRED",
                   instrument_platform="ILLUMINA", fastq_ftp="", fastq_md5="")
        row.update(changes)
        return Response(("\t".join(row) + "\n" + "\t".join(row.values()) + "\n").encode())

    def reads_row(self, names):
        return dict(run_accession="SRR1", library_layout="PAIRED",
                    fastq_ftp=";".join("ftp.sra.ebi.ac.uk/vol1/" + name for name in names),
                    fastq_md5=";".join(["a" * 32] * len(names)))

    def test_download_replaces_corrupt_cache_record(self):
        destination = self.root / "reference.gz"
        destination.write_bytes(b"old")
        Path(str(destination) + ".download.json").write_text("{", encoding="utf-8")
        with patch("urllib.request.urlopen", return_value=Response(b"new", 3)) as remote:
            download("https://example.org/reference.gz", destination)
        self.assertEqual(destination.read_bytes(), b"new")
        self.assertEqual(remote.call_count, 1)
        with patch("urllib.request.urlopen", side_effect=AssertionError("cache should be reused")):
            download("https://example.org/reference.gz", destination)

    def test_truncated_download_is_never_marked_complete(self):
        destination = self.root / "reference.gz"
        with patch("urllib.request.urlopen", side_effect=lambda *a, **k: Response(b"short", 100)), \
             patch("mouse_strain.pipeline.time.sleep"):
            with self.assertRaisesRegex(ValueError, "Content-Length"):
                download("https://example.org/reference.gz", destination)
        self.assertFalse(destination.exists())
        self.assertFalse(Path(str(destination) + ".partial").exists())
        self.assertFalse(Path(str(destination) + ".download.json").exists())

    def test_incomplete_http_read_is_retried(self):
        response = Mock()
        response.__enter__ = Mock(return_value=response)
        response.__exit__ = Mock(return_value=None)
        response.read.side_effect = http.client.IncompleteRead(b"partial", 100)
        destination = self.root / "reference.gz"
        with patch("urllib.request.urlopen", side_effect=[response, Response(b"complete", 8)]) as remote, \
             patch("mouse_strain.pipeline.time.sleep"):
            download("https://example.org/reference.gz", destination)
        self.assertEqual(remote.call_count, 2)
        self.assertEqual(destination.read_bytes(), b"complete")

    def test_exhausted_http_retry_reports_download_context(self):
        destination = self.root / "reference.gz"
        with patch("urllib.request.urlopen", side_effect=http.client.IncompleteRead(b"partial", 100)) as remote, \
             patch("mouse_strain.pipeline.time.sleep"):
            with self.assertRaisesRegex(RuntimeError, r"reference\.gz.*https://example\.org/reference\.gz") as caught:
                download("https://example.org/reference.gz", destination)
        self.assertIsInstance(caught.exception.__cause__, http.client.IncompleteRead)
        self.assertEqual(remote.call_count, 3)
        self.assertFalse(destination.exists())
        self.assertFalse(Path(str(destination) + ".download.json").exists())

    def test_malformed_ena_columns_fail_with_clear_error(self):
        for body in (b"run_accession\nSRR1\n", self.ena_response().getvalue().rsplit(b"\t", 1)[0]):
            with self.subTest(body=body), patch("urllib.request.urlopen", return_value=Response(body)):
                with self.assertRaisesRegex(ValueError, "Malformed ENA metadata"):
                    resolve_srx("SRX1")

    def test_missing_sample_cannot_be_pooled(self):
        with patch("urllib.request.urlopen", return_value=self.ena_response(sample_accession="")):
            with self.assertRaisesRegex(ValueError, "sample accession is missing"):
                resolve_srx("SRX1")

    def test_inconsistent_or_prefix_only_mouse_metadata_is_rejected(self):
        for changes in (dict(tax_id="10090", scientific_name="Homo sapiens"),
                        dict(tax_id="9606", scientific_name="Mus musculus"),
                        dict(tax_id="10090", scientific_name="Mus musculus_other")):
            with self.subTest(changes=changes), \
                 patch("urllib.request.urlopen", return_value=self.ena_response(**changes)):
                with self.assertRaisesRegex(ValueError, "Non-mouse"):
                    resolve_srx("SRX1")

    def test_mouse_subspecies_are_still_supported(self):
        with patch("urllib.request.urlopen", return_value=self.ena_response(
                tax_id="10092", scientific_name="Mus musculus domesticus")):
            self.assertEqual(resolve_srx("SRX1")[0]["run_accession"], "SRR1")

    def test_duplicate_or_unrelated_ena_mates_are_rejected(self):
        for names in (["SRR1_1.fastq.gz", "SRR1_1.fastq.gz"],
                      ["SRR1_1.fastq.gz", "SRR2_2.fastq.gz"],
                      ["SRR1.fastq.gz", "SRR1_1.fastq.gz"]):
            with self.subTest(names=names), patch("mouse_strain.pipeline.download") as remote:
                with self.assertRaisesRegex(ValueError, "complete FASTQ/checksum"):
                    get_reads(self.reads_row(names), self.root / "reads", 1, "ena")
                remote.assert_not_called()

    def test_complete_ena_pair_is_ordered_and_orphans_excluded(self):
        names = ["SRR1_2.fastq.gz", "SRR1.fastq.gz", "SRR1_1.fastq.gz"]
        with patch("mouse_strain.pipeline.download") as remote:
            paths = get_reads(self.reads_row(names), self.root / "reads", 1, "ena")
        self.assertEqual([p.name for p in paths], ["SRR1_1.fastq.gz", "SRR1_2.fastq.gz"])
        self.assertEqual(remote.call_count, 2)

    def test_corrupt_sra_cache_record_triggers_extraction(self):
        directory = self.root / "reads"
        directory.mkdir()
        for mate in (1, 2):
            (directory / f"SRR1_{mate}.fastq").write_text("old", encoding="utf-8")
        (directory / "sra_reads.json").write_text("{", encoding="utf-8")
        def fake_run(command):
            if command[0] == "fasterq-dump":
                for mate in (1, 2):
                    (directory / f"SRR1_{mate}.fastq").write_text("new", encoding="utf-8")
        with patch("mouse_strain.pipeline.require_tools"), \
             patch("mouse_strain.pipeline.run", side_effect=fake_run) as external:
            paths = get_reads(self.reads_row([]), directory, 1, "sra")
        self.assertEqual(external.call_count, 2)
        self.assertTrue(all(path.read_text(encoding="utf-8") == "new" for path in paths))

    def test_invalid_pipeline_inputs_do_not_create_outputs(self):
        for kwargs in (dict(threads=0), dict(threads=1.5), dict(download_method="invalid")):
            with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
                analyze("SRX1", "missing", "missing", "missing", self.root / "out", **kwargs)
        with self.assertRaises(ValueError):
            analyze("../SRX1", "missing", "missing", "missing", self.root / "out")
        self.assertFalse((self.root / "out").exists())
        with self.assertRaises(ValueError):
            get_reads(self.reads_row([]), self.root / "reads", 1, "invalid")
        self.assertFalse((self.root / "reads").exists())

    def test_failed_index_build_can_be_retried(self):
        fasta = self.root / "ref.fa"
        fasta.write_text(">1\nACGT\n", encoding="utf-8")
        prefix = self.root / "index" / "genome"
        def fake_failure(command):
            Path(str(prefix) + ".1.ht2").write_bytes(b"incomplete")
            raise subprocess.CalledProcessError(1, command)
        with patch("mouse_strain.pipeline.run", side_effect=fake_failure):
            with self.assertRaises(subprocess.CalledProcessError):
                build_index(fasta, prefix, "GRCm39", 1)
        self.assertFalse(Path(str(prefix) + ".1.ht2").exists())
        self.assertFalse(Path(str(prefix) + ".reference.json").exists())
        def fake_success(command):
            for i in range(1, 9):
                Path(f"{prefix}.{i}.ht2").write_bytes(b"index")
        with patch("mouse_strain.pipeline.run", side_effect=fake_success):
            build_index(fasta, prefix, "GRCm39", 1)
        self.assertTrue(Path(str(prefix) + ".reference.json").exists())

    def test_successful_exit_without_complete_index_is_rejected(self):
        fasta = self.root / "ref.fa"
        fasta.write_text(">1\nACGT\n", encoding="utf-8")
        prefix = self.root / "index"
        with patch("mouse_strain.pipeline.run"):
            with self.assertRaisesRegex(ValueError, "8 nonempty"):
                build_index(fasta, prefix, "GRCm39", 1)
        self.assertFalse(Path(str(prefix) + ".reference.json").exists())

    def test_empty_index_parts_and_malformed_provenance_are_rejected(self):
        prefix = self.root / "index"
        meta = dict(assembly="GRCm39", fasta_sha256="fixture")
        write_json(str(prefix) + ".reference.json", meta)
        for i in range(1, 9):
            Path(f"{prefix}.{i}.ht2").touch()
        with self.assertRaisesRegex(ValueError, "Incomplete HISAT2"):
            index_identity(prefix, meta)
        write_json(str(prefix) + ".reference.json", {})
        with self.assertRaisesRegex(ValueError, "provenance"):
            index_identity(prefix, meta)
        Path(str(prefix) + ".reference.json").write_text("{", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "provenance sidecar"):
            index_identity(prefix, meta)

    def test_corrupt_bam_marker_triggers_mapping(self):
        reads = self.root / "reads.fastq"
        reads.write_bytes(b"reads")
        bam = self.root / "mapped.bam"
        bam.write_bytes(b"old")
        Path(str(bam) + ".done.json").write_text("{", encoding="utf-8")
        mapper = SimpleNamespace(stdout=io.BytesIO(), poll=lambda: 0, wait=lambda: 0)
        def fake_sort(command, **kwargs):
            Path(command[command.index("-o") + 1]).write_bytes(b"new")
            return SimpleNamespace(returncode=0)
        with patch("mouse_strain.pipeline.tool_version", return_value="fixture"), \
             patch("mouse_strain.pipeline.subprocess.Popen", return_value=mapper), \
             patch("mouse_strain.pipeline.subprocess.run", side_effect=fake_sort), \
             patch("mouse_strain.pipeline.run"):
            map_run([reads], dict(run_accession="SRR1", sample_accession="SRS1"),
                    "index", {}, bam, 1)
        self.assertEqual(bam.read_bytes(), b"new")
        self.assertTrue(json.loads(Path(str(bam) + ".done.json").read_text(encoding="utf-8")))

    def test_duplicate_mapping_inputs_are_rejected(self):
        with self.assertRaisesRegex(ValueError, "distinct mate"):
            map_run(["same.fastq", "same.fastq"], {}, "index", {}, self.root / "mapped.bam", 1)

    def test_bam_check_failure_closes_panel_database(self):
        database = Mock()
        meta = dict(strains=["A", "B"], contigs={"1": 100})
        with patch("mouse_strain.pipeline.load_panel", return_value=(meta, database)), \
             patch("mouse_strain.pipeline.check_bam", side_effect=ValueError("invalid BAM")):
            with self.assertRaisesRegex(ValueError, "invalid BAM"):
                infer("input.bam", "ref.fa", "panel", self.root / "result")
        database.close.assert_called_once()

    def test_mpileup_start_failure_closes_database_and_log(self):
        database = Mock()
        meta = dict(strains=["A", "B"], contigs={"1": 100})
        logs = []
        def fail_to_start(*args, **kwargs):
            logs.append(kwargs["stderr"])
            raise FileNotFoundError("samtools")
        with patch("mouse_strain.pipeline.load_panel", return_value=(meta, database)), \
             patch("mouse_strain.pipeline.check_bam"), \
             patch("mouse_strain.pipeline.subprocess.Popen", side_effect=fail_to_start):
            with self.assertRaises(FileNotFoundError):
                infer("input.bam", "ref.fa", "panel", self.root / "result")
        database.close.assert_called_once()
        self.assertTrue(logs[0].closed)

    def test_reports_use_utf8_for_japanese_content(self):
        write_json(self.root / "data.json", {"strain": "系統A"})
        self.assertEqual(json.loads((self.root / "data.json").read_text(encoding="utf-8")),
                         {"strain": "系統A"})
        write_csv(self.root / "data.csv", [dict(strain="系統A")])
        self.assertIn("系統A", (self.root / "data.csv").read_text(encoding="utf-8"))
        render_html(self.root, [dict(chromosome="1", call="系統A")],
                    [dict(chromosome="1", start=1, end=10, status="candidate_supported",
                          call="系統A", n_sites=10)], ["系統A"])
        self.assertIn("系統A", (self.root / "report.html").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
