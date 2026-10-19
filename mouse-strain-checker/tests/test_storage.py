"""Destructive operations use synthetic jobs and never touch real analysis data."""
import json
from pathlib import Path
import tempfile
import threading
import unittest
from unittest.mock import patch

from mouse_strain.storage import retained_sex
from mouse_strain.web_jobs import Jobs, atomic_json
import test_web_server as fixtures


class StorageTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.base = Path(temporary.name)
        self.jobs = Jobs(self.base / "jobs", self.base / "resources")
        self.addCleanup(self.jobs.close)
        self.storage = self.jobs.storage
        self.jid = "a" * 32
        self.sample = self.make_job(self.jid)

    def put(self, path, data=b"synthetic data"):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(data)
        return path

    def make_job(self, jid, status="completed"):
        sample = self.jobs.root / jid / "output" / "SRX1"
        fixtures.sample(sample, "SRX1")
        job = dict(id=jid, status=status, created_at="2026-10-19", request=dict(srx=["SRX1"]))
        self.jobs.items[jid] = job
        self.jobs._events[jid] = threading.Event()
        atomic_json(self.jobs.root / jid / "job.json", job)
        self.put(sample / "cache" / "SRR1.bam")
        return sample

    def ready_references(self):
        root = self.jobs.resources
        fixtures.ready(root)
        required = [root / "reference.fa", root / "reference.fa.fai", root / "panel" / "metadata.json",
                    root / "panel" / "panel.sqlite", root / "panel" / "sites.bed"]
        required += [root / "hisat2" / f"genome.{i}.ht2" for i in range(1, 9)]
        for path in required:
            self.put(path)
        return required

    def cleanup(self, scope="cache", jid=None):
        preview = self.storage.preview(scope, self.jid if jid is None and scope != "reference_downloads" else jid)
        return self.storage.cleanup(preview["token"], True)

    def test_preview_is_read_only_and_capacity_matches_selected_files(self):
        raw = self.put(self.sample / "cache" / "SRR1" / "raw.fastq", b"A" * 200)
        before = {p: p.read_bytes() for p in self.jobs.root.rglob("*") if p.is_file()}
        overview = self.storage.overview()["jobs"][0]
        preview = self.storage.preview("cache", self.jid)
        self.assertEqual(preview["bytes"], 200 + len(b"synthetic data"))
        self.assertEqual(overview["cache_bytes"], preview["bytes"])
        self.assertEqual(preview["files"], 2)
        self.assertEqual(len(preview["token"]), 64)
        self.assertIn(str(raw.relative_to(self.jobs.root)), preview["paths"])
        self.assertEqual(before, {p: p.read_bytes() for p in self.jobs.root.rglob("*") if p.is_file()})

    def test_cache_cleanup_keeps_results_qc_and_sex_but_removes_all_intermediates(self):
        cache = self.sample / "cache"
        for name in ("raw/SRR1_1.fastq", "sra/SRR1/SRR1.sralite", ".orphan/temp.bin", "SRR1.bam.bai",
                     "preprocess/SRR1/paired_1.fastq.gz", "preprocess/.fastp-staging/fastp.json"):
            self.put(cache / name)
        retained = [self.put(cache / "preprocess/SRR1/fastp.html", b"QC report"),
                    self.put(cache / "preprocess/SRR1/fastp.json", b"{\"valid\":true}")]
        sex = dict(schema=1, status="male_likely", label="オス候補（RNA推定）", reason="Yマーカーを観測")
        atomic_json(self.sample / "sex.json", sex)
        saved = {p: p.read_bytes() for p in self.sample.iterdir() if p.is_file()}
        saved.update({p: p.read_bytes() for p in retained})
        result = self.cleanup()
        self.assertEqual(result["files"], 7)
        self.assertEqual(set(p for p in cache.rglob("*") if p.is_file()), set(retained))
        self.assertEqual(saved, {p: p.read_bytes() for p in saved})
        self.assertTrue(self.jobs.sample_complete(self.jid, "SRX1"))
        self.assertEqual(retained_sex(self.sample), sex)
        self.assertEqual(self.storage.overview()["jobs"][0]["cache_bytes"], 0)
        self.assertIn("last_cleanup", self.jobs.get(self.jid))

    def test_saved_sex_requires_matching_cleanup_hash(self):
        atomic_json(self.sample / "sex.json", dict(schema=1, status="male_likely", label="オス", reason="Y"))
        self.assertIsNone(retained_sex(self.sample))
        self.cleanup()
        self.assertIsNotNone(retained_sex(self.sample))
        atomic_json(self.sample / "sex.json", dict(schema=1, status="female_likely", label="メス", reason="Xist"))
        self.assertIsNone(retained_sex(self.sample))

    def test_full_job_cleanup_removes_history_and_only_selected_job(self):
        other = self.make_job("b" * 32)
        source = self.put(self.jobs.resources / "reference.fa")
        self.jobs._generation[self.jid] = 3
        self.cleanup("job")
        self.assertFalse((self.jobs.root / self.jid).exists())
        self.assertIsNone(self.jobs.get(self.jid))
        self.assertNotIn(self.jid, self.jobs._generation)
        self.assertTrue((other / "cache/SRR1.bam").is_file())
        self.assertEqual(source.read_bytes(), b"synthetic data")

    def test_every_active_status_is_rejected_for_both_scopes(self):
        for status in ("queued", "running", "pausing"):
            for scope in ("cache", "job"):
                with self.subTest(status=status, scope=scope):
                    self.jobs.items[self.jid]["status"] = status
                    with self.assertRaises(ValueError):
                        self.storage.preview(scope, self.jid)
        self.assertTrue((self.sample / "cache/SRR1.bam").is_file())

    def test_status_rechecked_after_preview_and_child_process_is_also_protected(self):
        for condition in ("queued", "running", "pausing", "process"):
            with self.subTest(condition=condition):
                self.jobs.items[self.jid]["status"] = "paused"
                preview = self.storage.preview("cache", self.jid)
                if condition == "process":
                    self.jobs._processes[self.jid] = object()
                else:
                    self.jobs.items[self.jid]["status"] = condition
                with self.assertRaises(ValueError):
                    self.storage.cleanup(preview["token"], True)
                self.jobs._processes.clear()
        self.jobs.items[self.jid]["status"] = "completed"
        self.assertTrue((self.sample / "cache/SRR1.bam").is_file())

    def test_inactive_job_can_be_cleaned_while_another_job_runs(self):
        self.make_job("b" * 32, "running")
        self.cleanup()
        self.assertTrue((self.jobs.root / ("b" * 32) / "output/SRX1/cache/SRR1.bam").is_file())

    def test_reference_download_cleanup_preserves_all_prepared_resources(self):
        required = self.ready_references()
        sources = [self.put(self.jobs.resources / name) for name in (
            "reference.fa.gz", "strains.vcf.gz", "strains.vcf.gz.download.json", "strains.vcf.gz.partial",
            "strains.vcf.gz.partial.json")]
        unrelated = self.put(self.jobs.resources / "user.fasta")
        result = self.cleanup("reference_downloads")
        self.assertEqual(result["files"], len(sources))
        self.assertTrue(all(not p.exists() for p in sources))
        self.assertTrue(all(p.read_bytes() == b"synthetic data" for p in required + [unrelated]))
        self.assertTrue((self.jobs.resources / "config.json").is_file())

    def test_reference_cleanup_requires_prepared_resources_and_idle_queue(self):
        source = self.put(self.jobs.resources / "strains.vcf.gz")
        with self.assertRaises((OSError, ValueError)):
            self.storage.preview("reference_downloads")
        self.ready_references()
        for status in ("queued", "running", "pausing"):
            self.jobs.items[self.jid]["status"] = "completed"
            preview = self.storage.preview("reference_downloads")
            self.jobs.items[self.jid]["status"] = status
            with self.assertRaises(ValueError):
                self.storage.cleanup(preview["token"], True)
        self.jobs.items[self.jid]["status"] = "completed"
        (self.jobs.resources / "hisat2/genome.8.ht2").unlink()
        with self.assertRaises(ValueError):
            self.storage.preview("reference_downloads")
        self.assertTrue(source.is_file())

    def test_confirmation_and_token_are_required_and_token_is_single_use(self):
        preview = self.storage.preview("cache", self.jid)
        for confirmed in (False, "true", 1, None):
            with self.assertRaises(ValueError):
                self.storage.cleanup(preview["token"], confirmed)
        with self.assertRaises(ValueError):
            self.storage.cleanup("guessed-token", True)
        self.storage.cleanup(preview["token"], True)
        with self.assertRaises(ValueError):
            self.storage.cleanup(preview["token"], True)

    def test_expired_or_changed_preview_never_deletes_files(self):
        preview = self.storage.preview("cache", self.jid)
        self.storage.previews[preview["token"]]["expires"] = 0
        with self.assertRaises(ValueError):
            self.storage.cleanup(preview["token"], True)
        preview = self.storage.preview("cache", self.jid)
        raw = self.put(self.sample / "cache/new.fastq")
        with self.assertRaises(ValueError):
            self.storage.cleanup(preview["token"], True)
        self.assertTrue(raw.is_file())
        preview = self.storage.preview("cache", self.jid)
        raw.write_bytes(b"changed size")
        with self.assertRaises(ValueError):
            self.storage.cleanup(preview["token"], True)
        self.assertTrue((self.sample / "cache/SRR1.bam").is_file())

    def test_invalid_scope_and_paths_never_choose_arbitrary_files(self):
        for scope, jid in (("cache", "../a"), ("job", str(self.base)), ("all", self.jid),
                           ("cache", "a" * 31), ("reference_downloads", self.jid), ([], self.jid)):
            with self.subTest(scope=scope, jid=jid), self.assertRaises(ValueError):
                self.storage.preview(scope, jid)
        self.assertTrue((self.sample / "cache/SRR1.bam").is_file())

    def test_links_inside_cache_or_replacing_job_are_rejected(self):
        external = self.put(self.base / "outside/keep.bam")
        link = self.sample / "cache/link"
        try:
            link.symlink_to(external.parent, target_is_directory=True)
        except OSError:
            self.skipTest("Symlinks require platform support")
        for scope in ("cache", "job"):
            with self.assertRaises(ValueError):
                self.storage.preview(scope, self.jid)
        self.assertEqual(external.read_bytes(), b"synthetic data")
        link.unlink()
        directory = self.jobs.root / self.jid
        moved = self.jobs.root / "saved-job"
        directory.rename(moved)
        directory.symlink_to(moved, target_is_directory=True)
        with self.assertRaises(ValueError):
            self.storage.preview("job", self.jid)

    def test_link_inserted_after_preview_cannot_escape_scope(self):
        preview = self.storage.preview("cache", self.jid)
        external = self.put(self.base / "outside/keep.bam")
        target = self.sample / "cache/SRR1.bam"
        target.unlink()
        try:
            target.symlink_to(external)
        except OSError:
            self.skipTest("Symlinks require platform support")
        with self.assertRaises(ValueError):
            self.storage.cleanup(preview["token"], True)
        self.assertEqual(external.read_bytes(), b"synthetic data")

    def test_io_failure_keeps_persistent_history_available_for_retry(self):
        preview = self.storage.preview("job", self.jid)
        original = Path.rmdir
        def fail_cache(path):
            if path == self.sample / "cache":
                raise PermissionError("injected fixture directory failure")
            return original(path)
        with patch.object(Path, "rmdir", fail_cache), self.assertRaises(OSError):
            self.storage.cleanup(preview["token"], True)
        self.assertTrue((self.jobs.root / self.jid / "job.json").is_file())
        self.assertIsNotNone(self.jobs.get(self.jid))
        self.cleanup("job")
        self.assertIsNone(self.jobs.get(self.jid))

    def test_failed_root_removal_restores_history(self):
        preview = self.storage.preview("job", self.jid)
        original = Path.rmdir
        def fail_root(path):
            if path == self.jobs.root / self.jid:
                raise PermissionError("injected fixture root failure")
            return original(path)
        with patch.object(Path, "rmdir", fail_root), self.assertRaises(OSError):
            self.storage.cleanup(preview["token"], True)
        persisted = json.loads((self.jobs.root / self.jid / "job.json").read_text(encoding="utf-8"))
        self.assertEqual(persisted["id"], self.jid)
        self.cleanup("job")


class StorageHTTPTests(unittest.TestCase):
    setUp = fixtures.HTTPTests.setUp
    tearDown = fixtures.HTTPTests.tearDown
    request = fixtures.HTTPTests.request

    def test_preview_confirm_and_deleted_history_roundtrip(self):
        overview = json.loads(self.request("GET", "/api/storage")[1])
        self.assertEqual(overview["jobs"][0]["id"], self.jid)
        status, data, _ = self.request("POST", "/api/storage/preview", dict(scope="job", job_id=self.jid))
        self.assertEqual(status, 200)
        preview = json.loads(data)
        self.assertTrue(self.directory.is_dir())
        self.assertEqual(self.request("POST", "/api/storage/cleanup", dict(token=preview["token"], confirm=False))[0], 400)
        self.assertTrue(self.directory.is_dir())
        self.assertEqual(self.request("POST", "/api/storage/cleanup", dict(token=preview["token"], confirm=True))[0], 200)
        self.assertFalse(self.directory.exists())
        self.assertEqual(json.loads(self.request("GET", "/api/jobs")[1]), [])
        self.assertEqual(self.request("GET", f"/api/jobs/{self.jid}")[0], 404)

    def test_cross_origin_schema_and_arbitrary_paths_are_rejected(self):
        for endpoint in ("preview", "cleanup"):
            self.assertEqual(self.request("POST", "/api/storage/" + endpoint, {},
                                          {"Origin": "http://attacker.example"})[0], 403)
            self.assertEqual(self.request("POST", "/api/storage/" + endpoint, {},
                                          {"Content-Type": "text/plain"})[0], 415)
        for payload in ([], dict(scope="job", job_id="../"), dict(scope="job", job_id=self.jid, path="/etc"),
                        dict(scope=["job"], job_id=self.jid)):
            self.assertEqual(self.request("POST", "/api/storage/preview", payload)[0], 400)
        for payload in ({}, dict(confirm=True, path=str(self.directory)), dict(token="x", confirm=True)):
            self.assertEqual(self.request("POST", "/api/storage/cleanup", payload)[0], 400)
        self.assertTrue(self.directory.is_dir())


if __name__ == "__main__":
    unittest.main()
