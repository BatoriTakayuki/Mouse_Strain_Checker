"""fastp provenance/restart coverage, plus offline tests with the real executable."""
import gzip
import json
from pathlib import Path
import random
import shutil
import subprocess
import tempfile
import unittest
from unittest.mock import patch

from mouse_strain.preprocess import trim_reads


def fastq(records):
    return "".join(f"@{name}\n{sequence}\n+\n{quality}\n" for name, sequence, quality in records)


class PreprocessTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.paths = [self.root / "SRR1_1.fastq", self.root / "SRR1_2.fastq"]
        for path in self.paths:
            path.write_text(fastq([("read", "ACGT" * 20, "I" * 80)]), encoding="ascii")
        self.row = dict(run_accession="SRR1", library_layout="PAIRED")
        self.output = self.root / "preprocessed" / "SRR1"
        self.executions = []
        self.version = "fastp 0.23.4"
        self.after = 2
        self.corrupt = False

    def tearDown(self):
        self.temporary.cleanup()

    def fake_run(self, command, **kwargs):
        if command == ["fastp", "--version"]:
            return subprocess.CompletedProcess(command, 0, stdout="", stderr=self.version + "\n")
        self.executions.append(command)
        for flag in ("-o", "-O"):
            if flag in command:
                destination = Path(command[command.index(flag) + 1])
                if self.corrupt:
                    destination.write_bytes(b"not gzip")
                else:
                    with gzip.open(destination, "wt", encoding="ascii") as handle:
                        handle.write(fastq([("trimmed", "ACGT" * 10, "I" * 40)]))
        Path(command[command.index("--html") + 1]).write_text("<html>QC report</html>")
        report = {"summary": {"before_filtering": {"total_reads": 4},
                              "after_filtering": {"total_reads": self.after,
                                                  "total_bases": self.after * 40}}}
        Path(command[command.index("--json") + 1]).write_text(json.dumps(report))
        return subprocess.CompletedProcess(command, 0)

    def run_trim(self, paths=None, row=None, threads=30):
        with patch("mouse_strain.preprocess.subprocess.run", side_effect=self.fake_run):
            return trim_reads(self.paths if paths is None else paths,
                              self.row if row is None else row, self.output, threads)

    def test_paired_command_reports_provenance_and_reuses_cache(self):
        originals = [path.read_bytes() for path in self.paths]
        reads = self.run_trim()
        self.assertEqual([path.name for path in reads], ["SRR1_1.trim.fastq.gz", "SRR1_2.trim.fastq.gz"])
        command = self.executions[0]
        self.assertIn("--detect_adapter_for_pe", command)
        self.assertEqual(command[command.index("--thread") + 1], "16")
        marker = json.loads((self.output / "complete.json").read_text())
        self.assertEqual(marker["provenance"]["tool"]["version"], self.version)
        self.assertEqual(marker["provenance"]["options"]["threads"], 16)
        self.assertEqual(marker["summary"], dict(reads_before=4, reads_after=2, bases_after=80,
                                                  fragments_after=1))
        self.assertEqual(len(marker["outputs"]), 4)
        self.assertEqual(self.run_trim(), reads)
        self.assertEqual(len(self.executions), 1)
        self.assertEqual([path.read_bytes() for path in self.paths], originals)
        self.assertFalse(list(self.output.glob(".fastp-*")))

    def test_single_end_uses_one_gzip_output(self):
        reads = self.run_trim(self.paths[:1], dict(run_accession="SRR1", library_layout="SINGLE"), 3)
        command = self.executions[0]
        self.assertEqual([path.name for path in reads], ["SRR1.trim.fastq.gz"])
        self.assertNotIn("--detect_adapter_for_pe", command)
        self.assertNotIn("-I", command)
        self.assertEqual(command[command.index("--thread") + 1], "3")

    def test_changed_raw_input_invalidates_cache(self):
        self.run_trim()
        self.paths[0].write_text(self.paths[0].read_text() + "\n")
        self.run_trim()
        self.assertEqual(len(self.executions), 2)

    def test_changed_tool_or_effective_thread_options_invalidates_cache(self):
        self.run_trim()
        self.version = "fastp 0.23.5"
        self.run_trim()
        self.run_trim(threads=4)
        self.assertEqual(len(self.executions), 3)

    def test_missing_or_changed_outputs_and_reports_invalidate_cache(self):
        reads = self.run_trim()
        for file in [reads[1], self.output / "fastp.json", self.output / "fastp.html"]:
            with self.subTest(file=file.name):
                file.unlink()
                self.run_trim()
        reads[0].write_bytes(b"truncated gzip")
        self.run_trim()
        self.assertEqual(len(self.executions), 5)

    def test_malformed_marker_is_not_a_completion_record(self):
        self.run_trim()
        (self.output / "complete.json").write_text("{")
        self.run_trim()
        self.assertEqual(len(self.executions), 2)

    def test_process_failure_invalidates_old_marker_and_cleans_only_own_staging(self):
        self.run_trim()
        sentinel = self.output / "keep.txt"
        sentinel.write_text("user data")
        self.paths[0].write_text(self.paths[0].read_text() + "\n")
        def fail(command, **kwargs):
            if "--version" in command:
                return self.fake_run(command, **kwargs)
            raise subprocess.CalledProcessError(1, command)
        with patch("mouse_strain.preprocess.subprocess.run", side_effect=fail):
            with self.assertRaises(subprocess.CalledProcessError):
                trim_reads(self.paths, self.row, self.output, 1)
        self.assertFalse((self.output / "complete.json").exists())
        self.assertEqual(sentinel.read_text(), "user data")
        self.assertFalse(list(self.output.glob(".fastp-*")))
        self.assertTrue(all(path.is_file() for path in self.paths))

    def test_input_mutation_during_processing_is_not_published(self):
        def mutate(command, **kwargs):
            result = self.fake_run(command, **kwargs)
            if "--version" not in command:
                self.paths[0].write_text("changed")
            return result
        with patch("mouse_strain.preprocess.subprocess.run", side_effect=mutate):
            with self.assertRaisesRegex(ValueError, "changed while"):
                trim_reads(self.paths, self.row, self.output, 1)
        self.assertFalse((self.output / "complete.json").exists())
        self.assertFalse(list(self.output.glob("*.fastq.gz")))

    def test_empty_or_inconsistent_filtering_summary_is_not_published(self):
        for value in (0, 1, 5, True):
            with self.subTest(reads_after=value):
                self.after = value
                with self.assertRaisesRegex(ValueError, "no usable reads|inconsistent"):
                    self.run_trim()
                self.assertFalse((self.output / "complete.json").exists())

    def test_invalid_gzip_is_not_published(self):
        self.corrupt = True
        with self.assertRaisesRegex(ValueError, "Invalid trimmed FASTQ"):
            self.run_trim()
        self.assertFalse((self.output / "complete.json").exists())

    def test_bad_inputs_rejected_before_running_tool(self):
        cases = [dict(threads=0), dict(threads=True), dict(threads=1.5),
                 dict(paths=self.paths[:1]), dict(paths=self.paths * 2),
                 dict(paths=[self.paths[0], self.paths[0]]),
                 dict(row=dict(run_accession="../SRR1", library_layout="PAIRED")),
                 dict(row=dict(run_accession="SRR1", library_layout="unknown"))]
        for arguments in cases:
            with self.subTest(arguments=arguments):
                with patch("mouse_strain.preprocess.subprocess.run") as run:
                    with self.assertRaises(ValueError):
                        trim_reads(arguments.get("paths", self.paths), arguments.get("row", self.row),
                                   self.output, arguments.get("threads", 1))
                    run.assert_not_called()

    def test_preprocessing_cannot_overwrite_raw_fastq(self):
        self.output.mkdir(parents=True)
        raw = self.output / "SRR1.trim.fastq.gz"
        raw.write_bytes(b"raw data")
        with self.assertRaisesRegex(ValueError, "overwrite input"):
            self.run_trim([raw], dict(run_accession="SRR1", library_layout="SINGLE"))
        self.assertEqual(raw.read_bytes(), b"raw data")


@unittest.skipUnless(shutil.which("fastp"), "fastp is required; run inside the supplied Docker image")
class RealFastpTests(unittest.TestCase):
    def test_real_fastp_filters_bad_reads_and_retains_pairing_for_se_and_pe(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            randomizer = random.Random(24269)
            good = ["".join(randomizer.choices("ACGT", k=100)) for _ in range(2)]
            for paired in (False, True):
                with self.subTest(paired=paired):
                    paths = []
                    for mate in range(2 if paired else 1):
                        path = root / f"input_{paired}_{mate}.fastq"
                        records = [(f"good{i}/{mate+1}", sequence, "I" * 100)
                                   for i, sequence in enumerate(good)]
                        records += [(f"bad/{mate+1}", "ACGT" * 25, "!" * 100),
                                    (f"short/{mate+1}", "ACGTACGT", "I" * 8)]
                        path.write_text(fastq(records), encoding="ascii")
                        paths.append(path)
                    output = root / ("paired" if paired else "single")
                    row = dict(run_accession="SRR1", library_layout="PAIRED" if paired else "SINGLE")
                    reads = trim_reads(paths, row, output, 1)
                    marker = json.loads((output / "complete.json").read_text())
                    self.assertEqual(marker["summary"]["reads_before"], 8 if paired else 4)
                    self.assertEqual(marker["summary"]["fragments_after"], 2)
                    for path in reads:
                        with gzip.open(path, "rt") as handle:
                            content = handle.read()
                        self.assertIn("@good0/", content)
                        self.assertIn("@good1/", content)
                        self.assertNotIn("@bad/", content)
                        self.assertNotIn("@short/", content)
                    self.assertEqual(trim_reads(paths, row, output, 1), reads)


if __name__ == "__main__":
    unittest.main()
