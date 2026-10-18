"""Regression coverage for panel quality checks and incomplete panel handling."""
import json
from pathlib import Path
import sqlite3
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch

from mouse_strain.panel import Fasta, build_panel, load_panel, vcf_samples


class PanelReviewTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.fasta = self.root / "ref.fa"
        # The FAI indexes bytes, so fixture offsets must not depend on OS newlines.
        self.fasta.write_bytes(b">1\n" + b"A" * 100 + b"\n")
        Path(str(self.fasta) + ".fai").write_bytes(b"1\t100\t3\t100\t101\n")
        self.vcf = self.root / "strains.vcf"
        self.output = self.root / "panel"

    def write_vcf(self, records=None, samples=("A", "B")):
        header = ("##fileformat=VCFv4.2\n##reference=GRCm39\n"
                  "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
                  + "\t".join(samples) + "\n")
        if records is None:
            records = ["1\t10\t.\tA\tG\t60\tPASS\t.\tGT:FI\t0/0:1\t1/1:1"]
        self.vcf.write_text(header + "\n".join(records) + "\n")

    def build(self, **kwargs):
        return build_panel(self.vcf, self.fasta, self.output, "GRCm39", spacing=1, **kwargs)

    def metadata(self, **changes):
        path = self.output / "metadata.json"
        meta = json.loads(path.read_text())
        meta.update(changes)
        path.write_text(json.dumps(meta))

    def test_nonfinite_gq_cannot_pass_quality_filter(self):
        records = [f"1\t{pos}\t.\tA\tG\t60\tPASS\t.\tGT:FI:GQ\t0/0:1:{gq}\t1/1:1:60"
                   for pos, gq in enumerate(("NaN", "Inf", "-Inf", ".", "19", "20"), 1)]
        self.write_vcf(records)
        meta = self.build()
        self.assertEqual(meta["stats"]["missing_het_or_low_quality"], 5)
        self.assertEqual((self.output / "sites.bed").read_text(), "1\t5\t6\n")

    def test_duplicate_sample_header_is_rejected_even_with_explicit_subset(self):
        self.write_vcf(samples=("A", "A", "B"))
        with self.assertRaisesRegex(ValueError, "sample names"):
            self.build(strains=["A", "B"])
        self.assertFalse(self.output.exists())

    def test_sites_only_vcf_has_no_samples(self):
        self.vcf.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        self.assertEqual(vcf_samples(self.vcf), [])

    def test_duplicate_format_keys_are_rejected(self):
        self.write_vcf(["1\t10\t.\tA\tG\t60\tPASS\t.\tGT:GT:FI\t1/1:0/0:1\t0/0:1/1:1"])
        with self.assertRaisesRegex(ValueError, "Duplicate VCF FORMAT"):
            self.build()
        self.assertFalse((self.output / "metadata.json").exists())

    def test_explicit_empty_strain_selection_does_not_select_every_sample(self):
        self.write_vcf()
        with self.assertRaisesRegex(ValueError, "At least two"):
            self.build(strains=[])
        self.assertFalse(self.output.exists())

    def test_nonfinite_threshold_is_rejected_before_output(self):
        self.write_vcf()
        with self.assertRaises(ValueError):
            self.build(min_gq=float("nan"))
        self.assertFalse(self.output.exists())

    def test_database_initialization_failure_closes_fasta(self):
        self.write_vcf()
        reference = SimpleNamespace(index={"1": (100, 3, 100, 101)}, close=Mock())
        with patch("mouse_strain.panel.Fasta", return_value=reference), \
             patch("mouse_strain.panel.sqlite3.connect", side_effect=OSError("disk full")):
            with self.assertRaisesRegex(OSError, "disk full"):
                self.build()
        reference.close.assert_called_once()
        self.assertFalse((self.output / "metadata.json").exists())

    def test_duplicate_fasta_contigs_are_rejected(self):
        Path(str(self.fasta) + ".fai").write_bytes(b"1\t100\t3\t100\t101\n" * 2)
        with self.assertRaisesRegex(ValueError, "Duplicate FASTA contig"):
            Fasta(self.fasta)

    def test_invalid_fasta_index_is_rejected(self):
        Path(str(self.fasta) + ".fai").write_bytes(b"1\t100\t3\t0\t1\n")
        with self.assertRaisesRegex(ValueError, "Invalid FASTA index"):
            Fasta(self.fasta)

    def test_provenance_failure_does_not_publish_panel(self):
        self.write_vcf()
        with patch("mouse_strain.panel.sha256", side_effect=OSError("read failed")):
            with self.assertRaisesRegex(OSError, "read failed"):
                self.build()
        self.assertFalse((self.output / "metadata.json").exists())
        self.assertFalse((self.output / "panel.sqlite").exists())

    def test_invalid_strain_metadata_is_rejected(self):
        self.write_vcf()
        self.build()
        for strains in (["A"], ["A", "A"], ["A", ""], ["A", []], "AB"):
            with self.subTest(strains=strains):
                self.metadata(strains=strains)
                with self.assertRaisesRegex(ValueError, "Panel strains"):
                    load_panel(self.output)

    def test_invalid_contig_metadata_is_rejected(self):
        self.write_vcf()
        self.build()
        for contigs in ({}, {"1": 0}, {"1": True}, {"1": 100, "chr1": 100}):
            with self.subTest(contigs=contigs):
                self.metadata(contigs=contigs)
                with self.assertRaisesRegex(ValueError, "Panel contigs"):
                    load_panel(self.output)

    def test_panel_missing_bed_is_incomplete(self):
        self.write_vcf()
        self.build()
        (self.output / "sites.bed").unlink()
        with self.assertRaisesRegex(ValueError, "Incomplete"):
            load_panel(self.output)

    def test_missing_or_invalid_fasta_checksum_is_rejected(self):
        self.write_vcf()
        self.build()
        for digest in (None, "", "incorrect", 123):
            with self.subTest(digest=digest):
                self.metadata(fasta_sha256=digest)
                with self.assertRaisesRegex(ValueError, "FASTA SHA-256"):
                    load_panel(self.output)

    def test_invalid_database_schema_is_rejected(self):
        self.write_vcf()
        self.build()
        with sqlite3.connect(self.output / "panel.sqlite") as db:
            db.execute("DROP TABLE snps")
        with self.assertRaisesRegex(ValueError, "Invalid panel database"):
            load_panel(self.output)


if __name__ == "__main__":
    unittest.main()
