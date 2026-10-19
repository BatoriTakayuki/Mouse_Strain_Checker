"""Offline checks for prefetch, validated conversion, and interrupted SRA caches."""
import json
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

from mouse_strain.pipeline import get_reads, identity


class SRAWorkflowTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.directory = self.root / "reads"
        self.row = dict(run_accession="SRR1", library_layout="PAIRED", fastq_ftp="", fastq_md5="")
        self.commands = []
        self.payload_name = "SRR1.sra"
        self.payload_subdirectory = False
        self.output_names = ["SRR1_1.fastq", "SRR1_2.fastq"]
        self.fail_conversion = False
        self.fail_validation = False
        self.addCleanup(patch.stopall)
        patch("mouse_strain.pipeline.require_tools").start()
        patch("mouse_strain.pipeline.tool_version", side_effect=lambda name: name + " fixture").start()
        patch("mouse_strain.pipeline.run", side_effect=self.fake_run).start()

    def fake_run(self, command, **kwargs):
        self.commands.append(list(command))
        if command[0] == "prefetch":
            target = Path(command[-1])
            if self.payload_subdirectory:
                target = target / "SRR1"
            target.mkdir(parents=True, exist_ok=True)
            if self.payload_name:
                payload = target / self.payload_name
                if not payload.exists():
                    payload.write_bytes(b"complete-sra")
        elif command[0] == "vdb-validate" and self.fail_validation:
            raise subprocess.CalledProcessError(3, command)
        elif command[0] == "fasterq-dump":
            output = Path(command[command.index("-O") + 1])
            for name in self.output_names:
                (output / name).write_bytes(b"@read\nACGT\n+\nIIII\n")
            if self.fail_conversion:
                raise subprocess.CalledProcessError(1, command)

    def reads(self, method="sra"):
        return get_reads(self.row, self.directory, 4, method)

    def marker(self):
        return json.loads((self.directory / "sra_reads.json").read_text(encoding="utf-8"))

    def write_marker(self, marker):
        (self.directory / "sra_reads.json").write_text(json.dumps(marker), encoding="utf-8")

    def test_default_and_auto_use_prefetch_even_with_ena_urls(self):
        self.row.update(fastq_ftp="ftp.sra.ebi.ac.uk/SRR1_1.fastq.gz;ftp.sra.ebi.ac.uk/SRR1_2.fastq.gz",
                        fastq_md5=";".join(["a" * 32] * 2))
        with patch("mouse_strain.pipeline.download") as download:
            get_reads(self.row, self.directory, 4)
            self.assertEqual([c[0] for c in self.commands], ["prefetch", "vdb-validate", "fasterq-dump"])
            self.commands.clear()
            get_reads(self.row, self.root / "auto", 4, "auto")
            self.assertEqual(self.commands[0][0], "prefetch")
            download.assert_not_called()

    def test_commands_and_provenance_allow_verified_cache_reuse(self):
        paths = self.reads()
        prefetch, validate, convert = self.commands
        accession_directory = self.directory / "sra" / "SRR1"
        self.assertEqual(prefetch, ["prefetch", "SRR1", "--max-size", "u",
                                   "--output-directory", accession_directory])
        self.assertEqual(validate, ["vdb-validate", accession_directory / "SRR1.sra"])
        self.assertEqual(convert[1], accession_directory)
        self.assertIn("--split-3", convert)
        self.assertIn("--skip-technical", convert)
        self.assertEqual(convert[convert.index("-e") + 1], 4)
        self.assertNotEqual(convert[convert.index("-O") + 1], self.directory)
        self.assertNotEqual(convert[convert.index("-O") + 1], convert[convert.index("-t") + 1])
        marker = self.marker()
        self.assertEqual(marker["schema"], 3)
        self.assertEqual(marker["options"], dict(split_mode="split-3", skip_technical=True))
        self.assertEqual(marker["method"], "sra")
        self.assertEqual(marker["accession"], "SRR1")
        self.assertEqual(marker["layout"], "PAIRED")
        self.assertEqual(marker["reads"], [identity(p) for p in paths])
        self.assertEqual(set(marker["tools"]), {"prefetch", "fasterq-dump", "vdb-validate"})
        self.commands.clear()
        self.assertEqual(self.reads("auto"), paths)
        self.assertEqual(self.commands, [])

    def test_exact_no_suffix_payload_is_supported(self):
        self.payload_name = "SRR1"
        self.reads()
        self.assertEqual(Path(self.commands[1][1]).name, "SRR1")
        self.assertEqual(Path(self.marker()["sra"]["path"]).name, "SRR1")

    def test_prefetch_nested_directory_is_used_for_validation_conversion_and_cache(self):
        # Actual prefetch 3.4.1 appends SRR1 to its requested output directory.
        self.payload_subdirectory = True
        for payload_name in ("SRR1.sra", "SRR1.sralite", "SRR1"):
            with self.subTest(payload_name=payload_name):
                self.directory = self.root / payload_name
                self.payload_name = payload_name
                self.commands.clear()
                paths = self.reads()
                actual = self.directory / "sra" / "SRR1" / "SRR1"
                self.assertEqual(self.commands[1], ["vdb-validate", actual / payload_name])
                self.assertEqual(self.commands[2][1], actual)
                self.assertEqual(self.marker()["sra"], identity(actual / payload_name))
                self.commands.clear()
                self.assertEqual(self.reads(), paths)
                self.assertEqual(self.commands, [])

    def test_sra_lite_is_validated_converted_and_quality_is_explicit(self):
        self.payload_name = "SRR1.sralite"
        self.reads()
        source = self.directory / "sra" / "SRR1" / "SRR1.sralite"
        self.assertEqual(self.commands[1], ["vdb-validate", source])
        self.assertEqual(self.commands[2][1], source.parent)
        quality = self.marker()["quality"]
        self.assertEqual(quality["format"], "sra_lite")
        self.assertEqual(quality["base_quality"], "simplified")
        self.assertIn("簡略化", quality["note"])
        self.commands.clear()
        self.reads()
        self.assertEqual(self.commands, [])

    def test_normalized_sra_is_preferred_when_lite_is_also_cached(self):
        folder = self.directory / "sra" / "SRR1"
        folder.mkdir(parents=True)
        (folder / "SRR1.sralite").write_bytes(b"lite")
        self.reads()
        self.assertEqual(self.commands[1], ["vdb-validate", folder / "SRR1.sra"])
        self.assertEqual(self.marker()["quality"]["base_quality"], "full")

    def test_flat_and_nested_completed_payloads_are_ambiguous(self):
        folder = self.directory / "sra" / "SRR1" / "SRR1"
        folder.mkdir(parents=True)
        (folder / "SRR1.sra").write_bytes(b"nested")
        with self.assertRaisesRegex(ValueError, "found 2"):
            self.reads()
        self.assertEqual([c[0] for c in self.commands], ["prefetch"])

    def test_nested_partial_or_unrelated_files_are_not_payloads(self):
        self.payload_name = None
        self.payload_subdirectory = True
        folder = self.directory / "sra" / "SRR1" / "SRR1"
        folder.mkdir(parents=True)
        for name in ("SRR1.sra.partial", "SRR1.sra.lock", "SRR999.sra"):
            (folder / name).write_bytes(b"unfinished")
        with self.assertRaisesRegex(ValueError, "found 0"):
            self.reads()
        self.assertEqual([c[0] for c in self.commands], ["prefetch"])

    def test_partial_and_lock_files_are_never_chosen_as_source(self):
        self.payload_name = None
        folder = self.directory / "sra" / "SRR1"
        folder.mkdir(parents=True)
        for name in ("SRR1.sra.partial", "SRR1.sra.lock", "SRR1.sra.tmp", "SRR999.sra",
                     "SRR1.sralite.partial", "SRR1.sralite.lock", "SRR999.sralite"):
            (folder / name).write_bytes(b"unfinished")
        with self.assertRaisesRegex(ValueError, "exactly one nonempty SRA"):
            self.reads()
        self.assertEqual([c[0] for c in self.commands], ["prefetch"])

    def test_two_exact_payloads_are_rejected(self):
        folder = self.directory / "sra" / "SRR1"
        folder.mkdir(parents=True)
        (folder / "SRR1").write_bytes(b"ambiguous")
        with self.assertRaisesRegex(ValueError, "exactly one"):
            self.reads()
        self.assertEqual([c[0] for c in self.commands], ["prefetch"])

    def test_validation_failure_preserves_existing_fastq_files(self):
        self.directory.mkdir()
        old = self.directory / "SRR1_1.fastq"
        old.write_bytes(b"prior-read")
        self.fail_validation = True
        with self.assertRaises(subprocess.CalledProcessError):
            self.reads()
        self.assertEqual(old.read_bytes(), b"prior-read")
        self.assertFalse((self.directory / "sra_reads.json").exists())
        self.assertEqual([c[0] for c in self.commands], ["prefetch", "vdb-validate"])

    def test_failed_conversion_preserves_old_cache_and_retry_uses_new_stage(self):
        paths = self.reads()
        old_reads = [p.read_bytes() for p in paths]
        marker = self.marker()
        marker["schema"] = 1
        self.write_marker(marker)
        old_marker = (self.directory / "sra_reads.json").read_bytes()
        self.fail_conversion = True
        self.output_names = ["SRR1_1.fastq"]
        with self.assertRaises(subprocess.CalledProcessError):
            self.reads()
        failed = self.commands[-1]
        self.assertEqual([p.read_bytes() for p in paths], old_reads)
        self.assertEqual((self.directory / "sra_reads.json").read_bytes(), old_marker)
        self.assertFalse(Path(failed[failed.index("-O") + 1]).exists())
        self.fail_conversion = False
        self.output_names = ["SRR1_1.fastq", "SRR1_2.fastq"]
        self.reads()
        retry = self.commands[-1]
        self.assertNotEqual(failed[failed.index("-O") + 1], retry[retry.index("-O") + 1])
        self.assertEqual(self.marker()["schema"], 3)

    def test_incomplete_pair_cannot_borrow_stale_final_or_interrupted_stage_mate(self):
        self.directory.mkdir()
        (self.directory / "SRR1_2.fastq").write_bytes(b"old-final")
        stale = self.directory / ".SRR1.fasterq-interrupted" / "reads"
        stale.mkdir(parents=True)
        (stale / "SRR1_2.fastq").write_bytes(b"old-staged")
        self.output_names = ["SRR1_1.fastq"]
        with self.assertRaisesRegex(ValueError, "declared layout"):
            self.reads()
        self.assertEqual((self.directory / "SRR1_2.fastq").read_bytes(), b"old-final")
        self.assertEqual((stale / "SRR1_2.fastq").read_bytes(), b"old-staged")
        self.assertFalse((self.directory / "sra_reads.json").exists())

    def test_empty_fastq_is_rejected(self):
        original_run = self.fake_run
        def empty_output(command, **kwargs):
            original_run(command, **kwargs)
            if command[0] == "fasterq-dump":
                output = Path(command[command.index("-O") + 1])
                (output / "SRR1_2.fastq").write_bytes(b"")
        with patch("mouse_strain.pipeline.run", side_effect=empty_output):
            with self.assertRaisesRegex(ValueError, "declared layout"):
                self.reads()

    def test_paired_orphans_are_excluded_from_published_fastq_and_cache(self):
        self.output_names = ["SRR1_1.fastq", "SRR1_2.fastq", "SRR1.fastq"]
        with patch("mouse_strain.pipeline.log") as log:
            paths = self.reads()
        self.assertEqual([p.name for p in paths], ["SRR1_1.fastq", "SRR1_2.fastq"])
        self.assertFalse((self.directory / "SRR1.fastq").exists())
        self.assertEqual([Path(p["path"]).name for p in self.marker()["reads"]],
                         ["SRR1_1.fastq", "SRR1_2.fastq"])
        log.assert_any_call("SRA orphan reads excluded; paired reads only")

    def test_paired_unrelated_output_is_rejected_even_when_both_mates_exist(self):
        self.output_names = ["SRR1_1.fastq", "SRR1_2.fastq", "SRR2.fastq"]
        with self.assertRaisesRegex(ValueError, "declared layout"):
            self.reads()
        self.assertFalse((self.directory / "sra_reads.json").exists())

    def test_schema_two_split_files_cache_is_regenerated_with_split_three(self):
        self.reads()
        marker = self.marker()
        marker.update(schema=2, options=dict(split_files=True, skip_technical=True))
        self.write_marker(marker)
        self.commands.clear()
        self.reads()
        self.assertEqual([c[0] for c in self.commands], ["prefetch", "vdb-validate", "fasterq-dump"])
        self.assertIn("--split-3", self.commands[-1])
        self.assertEqual(self.marker()["schema"], 3)

    def test_single_read_suffix_is_normalized_and_unexpected_second_mate_rejected(self):
        self.row["library_layout"] = "SINGLE"
        self.output_names = ["SRR1_1.fastq"]
        self.assertEqual([p.name for p in self.reads()], ["SRR1.fastq"])
        marker = self.marker()
        marker["schema"] = 1
        self.write_marker(marker)
        self.output_names = ["SRR1_1.fastq", "SRR1_2.fastq"]
        with self.assertRaisesRegex(ValueError, "declared layout"):
            self.reads()

    def test_old_marker_and_mismatched_provenance_are_cache_misses(self):
        self.reads()
        valid = self.marker()
        variants = [[identity(self.directory / f"SRR1_{mate}.fastq") for mate in (1, 2)]]
        for key, value in (("schema", 1), ("method", "ena"), ("accession", "SRR2"), ("layout", "SINGLE")):
            variants.append(dict(valid, **{key: value}))
        variants.append(dict(valid, tools={"prefetch": "old"}))
        for marker in variants:
            with self.subTest(marker=marker):
                self.write_marker(marker)
                self.commands.clear()
                self.reads()
                self.assertEqual([c[0] for c in self.commands], ["prefetch", "vdb-validate", "fasterq-dump"])


if __name__ == "__main__":
    unittest.main()
