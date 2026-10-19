import csv
from dataclasses import replace
import gzip
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from types import SimpleNamespace

from mouse_strain.core import Rules, Summary, count_fragments, parse_bases, site_scores
from mouse_strain.panel import build_panel, load_panel, sha256
from mouse_strain.pipeline import infer, resolve_srx, index_identity, analyze, map_run, get_reads
from mouse_strain.__main__ import main


def fasta_fixture(root):
    fasta = root / "ref.fa"
    fasta.write_text(">1\n" + "A" * 1000 + "\n>2\n" + "A" * 300 + "\n>X\n" + "A" * 100 + "\n")
    Path(str(fasta) + ".fai").write_text("1\t1000\t3\t1000\t1001\n2\t300\t1007\t300\t301\nX\t100\t1311\t100\t101\n")
    return fasta


def vcf_fixture(root, records=None, assembly="GRCm39"):
    path = root / "panel.vcf"
    header = ("##fileformat=VCFv4.2\n##reference=" + assembly + "\n"
              "##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>\n"
              "##FORMAT=<ID=FI,Number=1,Type=Integer,Description=Confidence>\n"
              "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\n")
    if records is None:
        records = [f"1\t{p}\t.\tA\tG\t60\tPASS\t.\tGT:FI\t0/0:1\t1/1:1" for p in (10, 30, 210, 230)]
    path.write_text(header + "\n".join(records) + "\n")
    return path


def pileup(pos, bases, names=None, nh=None, groups=None, chrom="1"):
    n = len(bases)
    return "\t".join([chrom, str(pos), "A", str(n), bases, "I" * n,
                     ",".join(names or [f"r{i}" for i in range(n)]),
                     ",".join(groups or ["run"] * n), ",".join(nh or ["1"] * n)]) + "\n"


class ScoringTests(unittest.TestCase):
    def setUp(self):
        self.rules = replace(Rules(), min_sites=2, min_window_sites=2, min_diagnostic=2)

    def test_pure_and_mixed_candidates(self):
        for counts, expected in [((10, 0), "A"), ((0, 10), "B"), ((5, 5), "A + B")]:
            summary = Summary(["A", "B"])
            for i in range(30):
                summary.add(i + 1, "01", *counts, self.rules)
            result = summary.result(self.rules)
            self.assertEqual(result["call"], expected)

    def test_deep_single_site_cannot_pass_site_requirement(self):
        summary = Summary(["A", "B"])
        summary.add(1, "01", 0, 100000, self.rules)
        self.assertEqual(summary.result(self.rules)["status"], "insufficient_sites")

    def test_observed_indistinguishable_strains_no_call(self):
        summary = Summary(["A", "B", "C"], pairs=False)
        for pos in range(20):
            summary.add(pos + 1, "001", 10, 0, self.rules)
        result = summary.result(self.rules)
        self.assertEqual(result["status"], "ambiguous")
        self.assertEqual(result["call"], "")
        self.assertEqual(result["diagnostic_sites"], 0)

    def test_poor_fit_is_unknown(self):
        summary = Summary(["A", "B"], pairs=False)
        for pos in range(20):
            summary.add(pos + 1, "01", 5, 5, self.rules)
        self.assertEqual(summary.result(self.rules)["status"], "poor_fit")

    def test_depth_cap(self):
        self.assertEqual(site_scores(10, 10, self.rules), site_scores(10000, 10000, self.rules))

    def test_no_reads(self):
        result = Summary(["A", "B"]).result(self.rules)
        self.assertEqual(result["status"], "no_observed_snps")
        self.assertEqual(result["candidate"], "")

    def test_invalid_rules(self):
        for kwargs in [dict(min_depth=0), dict(error_rate=0), dict(window_size=0), dict(min_baseq=-1)]:
            with self.assertRaises(ValueError):
                Rules(**kwargs)


class PileupTests(unittest.TestCase):
    def test_splice_indels_and_start_end_tokens(self):
        self.assertEqual(parse_bases("^].,+2agG-1t$*<>", "A"), ["A", "A", "G", "*", "<", ">"])

    def test_matching_mates_once_conflicting_mates_discarded(self):
        fields = pileup(10, "AAGGAG", names=["a", "a", "b", "b", "c", "c"]).strip().split("\t")
        self.assertEqual(count_fragments(fields, Rules()), {"A": 1, "G": 1})

    def test_read_groups_separate_runs(self):
        fields = pileup(10, "AG", names=["r", "r"], groups=["run1", "run2"]).strip().split("\t")
        self.assertEqual(count_fragments(fields, Rules()), {"A": 1, "G": 1})

    def test_nh_multimappers_and_missing_tags_are_excluded(self):
        fields = pileup(10, "AAG", nh=["1", "2", "*"]).strip().split("\t")
        self.assertEqual(count_fragments(fields, Rules()), {"A": 1})

    def test_low_base_quality(self):
        fields = pileup(10, "AG").strip().split("\t")
        fields[5] = "!I"
        self.assertEqual(count_fragments(fields, Rules()), {"G": 1})

    def test_malformed_output_fails(self):
        fields = pileup(10, "AG").strip().split("\t")
        fields[-1] = "1"
        with self.assertRaises(ValueError):
            count_fragments(fields, Rules())
        with self.assertRaises(ValueError):
            parse_bases("A+3ac", "A")


class PanelAndIntegrationTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.fasta = fasta_fixture(self.root)

    def tearDown(self):
        self.temp.cleanup()

    def prepare(self, records=None, **kwargs):
        vcf = vcf_fixture(self.root, records)
        output = self.root / "panel"
        meta = build_panel(vcf, self.fasta, output, "GRCm39", spacing=1, **kwargs)
        return output, meta

    def test_panel_provenance_and_bed_coordinates(self):
        output, meta = self.prepare()
        self.assertEqual(meta["sites_per_chromosome"], {"1": 4})
        self.assertEqual(meta["fasta_sha256"], sha256(self.fasta))
        self.assertEqual((output / "sites.bed").read_text().splitlines()[0], "1\t9\t10")
        loaded, db = load_panel(output)
        self.assertEqual(db.execute("SELECT gt FROM snps LIMIT 1").fetchone()[0], "01")
        db.close()

    def test_missing_low_fi_and_heterozygous_calls_not_reference(self):
        records = [
            "1\t10\t.\tA\tG\t60\tPASS\t.\tGT:FI\t./.:0\t1/1:1",
            "1\t20\t.\tA\tG\t60\tPASS\t.\tGT:FI\t0/0:1\t1/1:0",
            "1\t30\t.\tA\tG\t60\tPASS\t.\tGT:FI\t0/1:1\t1/1:1",
            "1\t40\t.\tA\tG\t60\tPASS\t.\tGT:FI\t0/0:1\t1/1:1"]
        output, meta = self.prepare(records)
        self.assertEqual(meta["sites_per_chromosome"], {"1": 1})
        self.assertEqual(meta["stats"]["missing_het_or_low_quality"], 3)

    def test_ref_mismatch_rejected(self):
        with self.assertRaisesRegex(ValueError, "REF mismatch"):
            self.prepare(["1\t10\t.\tC\tG\t60\tPASS\t.\tGT:FI\t0/0:1\t1/1:1"])
        self.assertFalse((self.root / "panel" / "metadata.json").exists())

    def test_assembly_mismatch_rejected(self):
        vcf = vcf_fixture(self.root, assembly="GRCm38")
        with self.assertRaisesRegex(ValueError, "assembly"):
            build_panel(vcf, self.fasta, self.root / "panel", "GRCm39")

    def test_unknown_strain_name_rejected(self):
        with self.assertRaisesRegex(ValueError, "absent strains"):
            self.prepare(strains=["not_in_vcf", "B"])

    def test_chr_prefix_aliases(self):
        output, _ = self.prepare(["chr1\t10\t.\tA\tG\t60\tPASS\t.\tGT:FI\t0/0:1\t1/1:1"])
        self.assertIn("1\t9\t10", (output / "sites.bed").read_text())

    def test_inference_regional_discordance_and_no_coverage(self):
        panel, _ = self.prepare()
        lines = [pileup(10, "AAAAAAAAAA"), pileup(30, "AAAAAAAAAA"),
                 pileup(210, "GGGGGGGGGG"), pileup(230, "GGGGGGGGGG")]
        rules = replace(Rules(), window_size=200, min_sites=2, min_window_sites=2, min_diagnostic=2)
        out = self.root / "result"
        rows = infer(None, self.fasta, panel, out, rules, pileup_lines=lines)
        self.assertEqual(rows[0]["status"], "regional_discordance")
        self.assertEqual(rows[0]["call"], "")
        self.assertEqual(rows[1]["status"], "no_observed_snps")
        self.assertEqual(rows[2]["status"], "no_observed_snps")
        with open(out / "windows.csv") as f:
            windows = list(csv.DictReader(f))
        self.assertEqual([r["call"] for r in windows[:2]], ["A", "B"])
        self.assertEqual(windows[2]["status"], "no_observed_snps")
        self.assertTrue((out / "report.html").exists())
        with gzip.open(out / "sites.tsv.gz", "rt") as f:
            self.assertEqual(len(f.readlines()), 5)

    def test_other_bases_and_insufficient_depth(self):
        panel, _ = self.prepare()
        lines = [pileup(10, "AAACCC"), pileup(30, "GG")]
        out = self.root / "result"
        rows = infer(None, self.fasta, panel, out, pileup_lines=lines)
        self.assertEqual(rows[0]["n_sites"], 0)
        self.assertEqual(json.loads((out / "evidence.json").read_text())["qc"]["filtered_panel_positions"], 2)

    def test_pileup_depth_limit_excluded(self):
        panel, _ = self.prepare()
        rows = infer(None, self.fasta, panel, self.root / "result", replace(Rules(), max_pileup_depth=10),
                     pileup_lines=[pileup(10, "A" * 10)])
        self.assertEqual(rows[0]["n_sites"], 0)

    def test_localized_reads_cannot_label_whole_chromosome(self):
        panel, _ = self.prepare()
        rules = replace(Rules(), window_size=200, min_sites=2, min_window_sites=2, min_diagnostic=2)
        rows = infer(None, self.fasta, panel, self.root / "result", rules,
                     pileup_lines=[pileup(10, "A" * 10), pileup(30, "A" * 10)])
        self.assertEqual(rows[0]["status"], "localized_evidence")
        self.assertEqual(rows[0]["candidate"], "A")
        self.assertEqual(rows[0]["call"], "")

    def test_sra_accession_directory_and_fastq_paths(self):
        directory = self.root / "reads"
        commands = []
        def fake_run(command, **kwargs):
            commands.append(list(map(str, command)))
            if command[0] == "prefetch":
                accession_directory = Path(command[-1])
                accession_directory.mkdir(parents=True)
                (accession_directory / "SRR1.sra").write_text("sra")
            if command[0] == "fasterq-dump":
                output = Path(command[command.index("-O") + 1])
                (output / "SRR1_1.fastq").write_text("read1")
                (output / "SRR1_2.fastq").write_text("read2")
        row = dict(run_accession="SRR1", fastq_ftp="", fastq_md5="", library_layout="PAIRED")
        with patch("mouse_strain.pipeline.run", side_effect=fake_run), \
             patch("mouse_strain.pipeline.require_tools"), \
             patch("mouse_strain.pipeline.tool_version", return_value="fixture"):
            paths = get_reads(row, directory, 4)
        self.assertEqual(Path(commands[0][-1]).name, "SRR1")
        self.assertEqual(commands[1], ["vdb-validate", str(directory / "sra/SRR1/SRR1.sra")])
        self.assertEqual(commands[2][1], str(directory / "sra/SRR1"))
        self.assertIn("--split-3", commands[2])
        self.assertEqual([p.name for p in paths], ["SRR1_1.fastq", "SRR1_2.fastq"])

    def test_mapper_failure_not_hidden_by_successful_sort(self):
        import io
        reads = self.root / "read.fastq"
        reads.write_text("fixture")
        mapper = SimpleNamespace(stdout=io.BytesIO(), wait=lambda: 1, poll=lambda: 1)
        output = self.root / "SRR1.bam"
        with patch("mouse_strain.pipeline.tool_version", return_value="fixture"), \
             patch("subprocess.Popen", return_value=mapper), \
             patch("subprocess.run", return_value=SimpleNamespace(returncode=0)):
            with self.assertRaisesRegex(RuntimeError, "HISAT2=1"):
                map_run([reads], dict(run_accession="SRR1", sample_accession="SRS1"),
                        "index", {}, output, 2)
        self.assertFalse(output.exists())
        self.assertFalse(Path(str(output) + ".done.json").exists())

    def test_external_bam_assembly_must_be_explicit(self):
        panel, _ = self.prepare()
        with self.assertRaisesRegex(ValueError, "bam-assembly"):
            analyze("local", self.fasta, panel, None, self.root / "result", bam="x.bam")

    def test_index_provenance_guard(self):
        panel, meta = self.prepare()
        prefix = self.root / "index"
        with self.assertRaisesRegex(ValueError, "reference.json"):
            index_identity(prefix, meta)
        Path(str(prefix) + ".reference.json").write_text(json.dumps(dict(assembly="GRCm38", fasta_sha256="wrong")))
        with self.assertRaisesRegex(ValueError, "provenance"):
            index_identity(prefix, meta)

    def test_invalid_cli_accession_fails_before_any_output(self):
        self.assertEqual(main(["run", "--srx", "../bad", "--out", str(self.root / "out"),
                               "--fasta", "ref", "--panel", "panel", "--index", "index"]), 1)
        self.assertFalse((self.root / "out").exists())

    def test_failed_rerun_invalidates_old_completion(self):
        panel, _ = self.prepare()
        out = self.root / "out"
        out.mkdir()
        (out / "complete.json").write_text("{}")
        (out / "report.html").write_text("old")
        (out / "sex.json").write_text('{"label":"old"}')
        with self.assertRaises(ValueError):
            analyze("SRX1", self.fasta, panel, "missing", out)
        self.assertFalse((out / "complete.json").exists())
        self.assertFalse((out / "report.html").exists())
        self.assertFalse((out / "sex.json").exists())


class EnaTests(unittest.TestCase):
    def test_input_path_injection(self):
        for value in ("../SRX1", "SRR1", "SRX1;echo", "SRX1/x"):
            with self.assertRaises(ValueError):
                resolve_srx(value)

    def test_species_guard(self):
        class Response:
            def __enter__(self): return self
            def __exit__(self, *args): return None
            def read(self):
                return ("run_accession\texperiment_accession\tsample_accession\ttax_id\tscientific_name\tlibrary_strategy\tlibrary_source\tlibrary_layout\tinstrument_platform\tfastq_ftp\tfastq_md5\n"
                        "SRR1\tSRX1\tSRS1\t9606\tHomo sapiens\tRNA-Seq\tTRANSCRIPTOMIC\tPAIRED\tILLUMINA\t\t\n").encode()
        with patch("urllib.request.urlopen", return_value=Response()):
            with self.assertRaisesRegex(ValueError, "Non-mouse"):
                resolve_srx("SRX1")


if __name__ == "__main__":
    unittest.main()
