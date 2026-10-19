"""Abstention, fragment accounting and reference safeguards for RNA sex markers."""
from pathlib import Path
from types import SimpleNamespace
import subprocess
import tempfile
import unittest
from unittest.mock import patch

from mouse_strain.sex import (ASSEMBLY_LENGTHS, MARKERS, classify_sex,
                              collect_sex_evidence, count_marker_alignments,
                              unavailable_sex)


def marker(gene):
    return next(m for m in MARKERS if m["gene"] == gene)


def sam(name, gene="Ddx3y", position=None, cigar="25M", flag=0, mapq=60,
        nh="NH:i:1", rg="run1", chrom=None):
    m = marker(gene)
    fields = [name, str(flag), chrom or m["chromosome"],
              str(position if position is not None else m["exons"][0][0]),
              str(mapq), cigar, "*", "0", "0", "A" * 25, "I" * 25]
    if nh:
        fields.append(nh)
    if rg:
        fields.append("RG:Z:" + rg)
    return "\t".join(fields) + "\n"


def records(gene, count, prefix=None):
    start = marker(gene)["exons"][0][0]
    return [sam(f"{prefix or gene}-{i}", gene, position=start + i) for i in range(count)]


def xist_records():
    lines = records("Xist", 20)
    exons = marker("Xist")["exons"]
    for i, offset in enumerate((0, 1, 2, 0)):
        left, right = (exons[0], exons[1]) if i < 3 else (exons[1], exons[2])
        anchor = 10 - offset
        lines.append(sam(f"xist-junction-{i}", "Xist", position=left[1] - anchor + 1,
                         cigar=f"{anchor}M{right[0] - left[1] - 1}N10M"))
    return lines


class SexMarkerTests(unittest.TestCase):
    def test_no_y_never_implies_female(self):
        result = classify_sex(count_marker_alignments(records("Actb", 60) + records("Gapdh", 60)))
        self.assertEqual(result["status"], "inconclusive")
        self.assertIn("Yの未検出", result["reason"])

    def test_no_expression_is_inconclusive(self):
        result = classify_sex(count_marker_alignments([]))
        self.assertEqual(result["status"], "inconclusive")
        self.assertEqual(result["evidence"]["y_fragments"], 0)

    def test_two_y_genes_independent_fragments_support_male(self):
        result = classify_sex(count_marker_alignments(records("Ddx3y", 10) + records("Uty", 10)))
        self.assertEqual(result["status"], "male_likely")
        self.assertEqual(result["evidence"]["y_genes_supported"], 2)
        self.assertEqual(result["evidence"]["y_fragments"], 20)
        self.assertIn("RNA", result["label"])

    def test_one_y_gene_is_not_enough(self):
        result = classify_sex(count_marker_alignments(records("Ddx3y", 50)))
        self.assertEqual(result["status"], "inconclusive")

    def test_repeated_start_positions_are_not_independent_evidence(self):
        lines = [sam(f"dup-{gene}-{i}", gene) for gene in ("Ddx3y", "Uty") for i in range(30)]
        result = classify_sex(count_marker_alignments(lines))
        self.assertEqual(result["status"], "inconclusive")
        self.assertEqual(result["evidence"]["y_genes_supported"], 0)

    def test_female_requires_controls_and_xist_specific_splicing(self):
        controls = records("Actb", 60) + records("Gapdh", 60)
        result = classify_sex(count_marker_alignments(controls + xist_records()))
        self.assertEqual(result["status"], "female_likely")
        self.assertEqual(result["evidence"]["xist_junctions"], 2)
        self.assertEqual(result["evidence"]["xist_spliced_fragments"], 4)
        # The same locus-only evidence cannot distinguish unstranded Tsix.
        locus_only = classify_sex(count_marker_alignments(controls + records("Xist", 30)))
        self.assertEqual(locus_only["status"], "inconclusive")
        self.assertIn("Tsix", locus_only["reason"])
        no_controls = classify_sex(count_marker_alignments(xist_records()))
        self.assertEqual(no_controls["status"], "inconclusive")

    def test_any_y_fragment_prevents_female_call(self):
        lines = records("Actb", 60) + records("Gapdh", 60) + xist_records() + records("Ddx3y", 1)
        self.assertEqual(classify_sex(count_marker_alignments(lines))["status"], "inconclusive")

    def test_strong_xist_and_y_are_conflicting_evidence(self):
        lines = records("Ddx3y", 10) + records("Uty", 10) + xist_records()
        result = classify_sex(count_marker_alignments(lines))
        self.assertEqual(result["status"], "inconclusive")
        self.assertIn("併存", result["reason"])

    def test_locus_only_signal_is_not_confirmed_xist_conflicting_with_y(self):
        lines = records("Ddx3y", 10) + records("Uty", 10) + records("Xist", 30)
        result = classify_sex(count_marker_alignments(lines))
        self.assertEqual(result["status"], "male_likely")
        self.assertIn("Tsixを区別できません", result["reason"])

    def test_read_pairs_deduplicate_but_read_groups_remain_distinct(self):
        lines = [sam("same", flag=65), sam("same", flag=129, position=1260800),
                 sam("same", rg="run2")]
        evidence = next(m for m in count_marker_alignments(lines) if m["gene"] == "Ddx3y")
        self.assertEqual(evidence["fragments"], 2)

    def test_low_quality_nonunique_and_excluded_flags_are_ignored(self):
        lines = [sam("good"), sam("mapq", mapq=29), sam("no-nh", nh=""),
                 sam("multi", nh="NH:i:2"), sam("nh-wrong", nh="NH:Z:1")]
        lines += [sam(f"flag-{flag}", flag=flag) for flag in (4, 256, 512, 1024, 2048)]
        evidence = next(m for m in count_marker_alignments(lines) if m["gene"] == "Ddx3y")
        self.assertEqual(evidence["fragments"], 1)

    def test_skipped_or_deleted_exonic_region_does_not_count(self):
        lines = [sam("skipped", position=1260000, cigar="10M30000N10M"),
                 sam("deleted", position=1260000, cigar="10M30000D10M"),
                 sam("intronic", position=1263400)]
        evidence = next(m for m in count_marker_alignments(lines) if m["gene"] == "Ddx3y")
        self.assertEqual(evidence["fragments"], 0)

    def test_short_anchors_do_not_count_as_specific_xist_splicing(self):
        end, begin = marker("Xist")["exons"][0][1], marker("Xist")["exons"][1][0]
        line = sam("weak-splice", "Xist", position=end - 3, cigar=f"4M{begin-end-1}N25M")
        evidence = next(m for m in count_marker_alignments([line]) if m["gene"] == "Xist")
        self.assertEqual(evidence["fragments"], 1)
        self.assertEqual(evidence["spliced_fragments"], 0)

    def test_expression_cap_does_not_hide_later_splice_evidence(self):
        # Coordinate order encounters first-exon expression before the splice
        # donors. The independent splice counter and junction set must keep
        # processing even after their preceding counters have reached a cap.
        lines = sorted(xist_records(), key=lambda line: int(line.split("\t")[3]))
        with patch("mouse_strain.sex.MAX_FRAGMENTS", 3):
            evidence = next(m for m in count_marker_alignments(lines) if m["gene"] == "Xist")
        self.assertEqual(evidence["fragments"], 3)
        self.assertEqual(evidence["distinct_starts"], 3)
        self.assertTrue(evidence["count_capped"])
        self.assertEqual(evidence["spliced_fragments"], 3)
        self.assertEqual(evidence["splice_junctions"], 2)

    def test_chr_prefix_and_sources_are_supported(self):
        evidence = count_marker_alignments([sam("prefix", chrom="chrY")])
        result = classify_sex(evidence)
        self.assertEqual(result["evidence"]["y_fragments"], 1)
        self.assertTrue(all(m["source"].startswith("https://rest.ensembl.org/") for m in result["markers"]))
        self.assertTrue(result["notes"])

    def test_bad_reference_is_rejected_before_any_subprocess(self):
        with patch("mouse_strain.sex.subprocess.run") as run:
            with self.assertRaisesRegex(ValueError, "GRCm39 only"):
                collect_sex_evidence("none.bam", "GRCm38", ASSEMBLY_LENGTHS)
            with self.assertRaisesRegex(ValueError, "complete GRCm39"):
                collect_sex_evidence("none.bam", "GRCm39", {"X": 100, "Y": 100})
            run.assert_not_called()

    def test_collector_checks_bam_header_then_uses_indexed_bounded_command(self):
        with tempfile.TemporaryDirectory() as temporary:
            bam = Path(temporary) / "sample.bam"
            bam.touch()
            header = "\n".join(f"@SQ\tSN:chr{k}\tLN:{v}" for k, v in ASSEMBLY_LENGTHS.items())
            def fake_run(command, **kwargs):
                if "-H" in command:
                    return SimpleNamespace(stdout=header)
                self.assertIn("-M", command)
                self.assertEqual(kwargs["timeout"], 180)
                self.assertIn("chrY:1260715-1286663", command)
                kwargs["stdout"].write("".join(records("Ddx3y", 10) + records("Uty", 10)))
                return SimpleNamespace(returncode=0)
            with patch("mouse_strain.sex.subprocess.run", side_effect=fake_run) as run:
                self.assertEqual(collect_sex_evidence(bam, "GRCm39", ASSEMBLY_LENGTHS)["status"], "male_likely")
                self.assertEqual(run.call_count, 2)

    def test_bad_bam_header_and_subprocess_failure_never_become_female(self):
        with tempfile.TemporaryDirectory() as temporary:
            bam = Path(temporary) / "sample.bam"
            bam.touch()
            with patch("mouse_strain.sex.subprocess.run", return_value=SimpleNamespace(stdout="@SQ\tSN:X\tLN:100")):
                with self.assertRaisesRegex(ValueError, "complete GRCm39"):
                    collect_sex_evidence(bam, "GRCm39", ASSEMBLY_LENGTHS)
            with patch("mouse_strain.sex.subprocess.run", side_effect=subprocess.TimeoutExpired("samtools", 30)):
                with self.assertRaises(subprocess.TimeoutExpired):
                    collect_sex_evidence(bam, "GRCm39", ASSEMBLY_LENGTHS)
            result = unavailable_sex("BAM unavailable", "GRCm39")
            self.assertEqual(result["status"], "unavailable")
            self.assertEqual(result["markers"], [])


if __name__ == "__main__":
    unittest.main()
