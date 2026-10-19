"""Persistent queue, pause/resume and local HTTP API regression tests."""
import csv
from http.client import HTTPConnection
from http.server import ThreadingHTTPServer
import json
from pathlib import Path
import tempfile
import threading
import time
import unittest
from unittest.mock import patch

from mouse_strain.web_jobs import Jobs, REQUIRED_RESULTS, atomic_json, validate_request
from mouse_strain.webapp import handler


def ready(directory):
    atomic_json(directory / "config.json", dict(fasta="reference.fa", panel="panel", index="hisat2/genome"))


def sample(directory, accession):
    directory.mkdir(parents=True, exist_ok=True)
    for name in REQUIRED_RESULTS:
        if name.endswith(".csv"):
            (directory / name).write_text("chromosome,model\n1,C57BL_6J\n", encoding="utf-8")
        else:
            (directory / name).write_text("{}" if name.endswith(".json") else "<html>結果</html>", encoding="utf-8")
    atomic_json(directory / "complete.json", dict(accession=accession))


def wait_status(jobs, jid, status):
    deadline = time.monotonic() + 5
    while time.monotonic() < deadline:
        job = jobs.get(jid)
        if job["status"] == status:
            return job
        time.sleep(0.01)
    raise AssertionError(f"Expected {status}, got {jobs.get(jid)}")


class FakeProcess:
    def __init__(self, callback, gate=None):
        self.returncode = None
        self.callback, self.gate, self.stopped = callback, gate, False

    def poll(self):
        if self.returncode is None and (self.gate is None or self.gate.is_set()):
            self.callback()
            self.returncode = 0
        return self.returncode


class QueueTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.resources = self.root / "resources"
        self.calls, self.processes = [], []
        self.gates = {}
        self.launcher = patch("mouse_strain.web_jobs.launch", side_effect=self.launch).start()
        self.stopper = patch("mouse_strain.web_jobs.stop", side_effect=self.stop).start()
        self.jobs = Jobs(self.root / "jobs", self.resources, 3)

    def tearDown(self):
        self.jobs.close()
        patch.stopall()
        self.temporary.cleanup()

    def launch(self, command, log_path, **options):
        self.calls.append((command, log_path, options))
        if command[3] == "setup":
            callback, gate = lambda: ready(self.resources), self.gates.get("setup")
        else:
            accession = command[command.index("--srx") + 1]
            output = Path(command[command.index("--out") + 1]) / accession
            callback, gate = lambda: sample(output, accession), self.gates.get(accession)
        process = FakeProcess(callback, gate)
        self.processes.append(process)
        return process

    def stop(self, process):
        process.stopped, process.returncode = True, -15

    def test_validates_and_deduplicates_accessions_without_paths(self):
        request = validate_request(dict(srx="SRX1, ERX2\nSRX1"), 3)
        self.assertEqual(request["srx"], ["SRX1", "ERX2"])
        self.assertEqual(request["threads"], 3)
        for payload in (dict(srx="../SRX1"), dict(srx="SRX１２"), dict(srx="SRX1", threads=True),
                        dict(srx="SRX1", resources="/etc"), dict(srx="SRX1", single_only="false"),
                        dict(srx=[]), dict(srx=[f"SRX{i}" for i in range(31)]), []):
            with self.subTest(payload=payload), self.assertRaises(ValueError):
                validate_request(payload)

    def test_automatic_setup_pause_and_resume_once(self):
        job = self.jobs.submit(dict(srx="SRX1", pause_after_reference=True, single_only=True))
        paused = wait_status(self.jobs, job["id"], "paused")
        self.assertTrue(paused["reference_ready"])
        self.assertEqual(len(self.calls), 1)
        self.assertIn("--resume", self.calls[0][0])
        self.jobs.resume(job["id"])
        final = wait_status(self.jobs, job["id"], "completed")
        self.assertEqual(final["completed_accessions"], ["SRX1"])
        self.assertEqual(len(self.calls), 2)
        self.assertIn("--single-only", self.calls[-1][0])
        command = self.calls[-1][0]
        self.assertEqual(command[command.index("--download-method") + 1], "sra")
        self.assertTrue(self.calls[-1][2]["append"])
        self.assertEqual(self.calls[-1][2]["env"]["MOUSE_STRAIN_PROGRESS_FILE"],
                         str(self.jobs.root / job["id"] / "progress.json"))
        self.assertEqual(self.jobs.get(job["id"])["id"], job["id"])

    def test_setup_only_completes_without_analysis(self):
        job = self.jobs.submit(dict(mode="setup"))
        wait_status(self.jobs, job["id"], "completed")
        self.assertEqual([call[0][3] for call in self.calls], ["setup"])

    def test_running_pause_stops_child_and_resume_skips_finished_sample(self):
        ready(self.resources)
        self.gates["SRX2"] = threading.Event()
        job = self.jobs.submit(dict(srx=["SRX1", "SRX2"]))
        deadline = time.monotonic() + 5
        while len(self.calls) < 2 and time.monotonic() < deadline:
            time.sleep(0.01)
        self.assertEqual(len(self.calls), 2)
        self.assertEqual(self.jobs.pause(job["id"])["status"], "pausing")
        paused = wait_status(self.jobs, job["id"], "paused")
        self.assertTrue(self.processes[-1].stopped)
        self.assertFalse(paused["pause_requested"])
        self.gates["SRX2"].set()
        self.jobs.resume(job["id"])
        wait_status(self.jobs, job["id"], "completed")
        self.assertEqual([call[0][-1] for call in self.calls], ["SRX1", "SRX2", "SRX2"])

    def test_old_queue_tokens_cannot_execute_a_resumed_job_twice(self):
        ready(self.resources)
        self.gates["SRX1"] = threading.Event()
        first = self.jobs.submit(dict(srx="SRX1"))
        wait_status(self.jobs, first["id"], "running")
        second = self.jobs.submit(dict(srx="SRX2"))
        for _ in range(3):
            self.jobs.pause(second["id"])
            self.jobs.resume(second["id"])
        self.gates["SRX1"].set()
        wait_status(self.jobs, first["id"], "completed")
        wait_status(self.jobs, second["id"], "completed")
        self.jobs.queue.join()
        self.assertEqual([call[0][-1] for call in self.calls], ["SRX1", "SRX2"])

    def test_progress_is_structured_and_cannot_overwrite_status(self):
        ready(self.resources)
        self.gates["SRX1"] = threading.Event()
        job = self.jobs.submit(dict(srx="SRX1"))
        wait_status(self.jobs, job["id"], "running")
        deadline = time.monotonic() + 5
        while not self.calls and time.monotonic() < deadline:
            time.sleep(0.01)
        atomic_json(self.jobs.root / job["id"] / "progress.json",
                    dict(stage="mapping", message="処理中", status="completed", id="bad",
                         transfer={"done_bytes": 10, "total_bytes": 100}, runs_done=1))
        self.jobs._read_progress(job["id"])
        result = self.jobs.get(job["id"])
        self.assertEqual(result["status"], "running")
        self.assertEqual(result["stage"], "mapping")
        self.assertEqual(result["id"], job["id"])
        self.assertEqual(result["transfer"]["done_bytes"], 10)
        self.assertEqual(result["runs_done"], 1)

    def test_close_stops_active_children_and_persists_paused(self):
        ready(self.resources)
        self.gates["SRX1"] = threading.Event()
        job = self.jobs.submit(dict(srx="SRX1"))
        deadline = time.monotonic() + 5
        while not self.calls and time.monotonic() < deadline:
            time.sleep(0.01)
        self.jobs.close()
        self.assertFalse(self.jobs._worker.is_alive())
        self.assertTrue(self.processes[-1].stopped)
        self.assertEqual(self.jobs.get(job["id"])["status"], "paused")
        with self.assertRaises(ValueError):
            self.jobs.submit(dict(srx="SRX2"))

    def test_restart_recovers_running_and_pausing_jobs_as_paused(self):
        self.jobs.close()
        for index, status in enumerate(("queued", "running", "pausing")):
            jid = f"{index:032x}"
            atomic_json(self.jobs.root / jid / "job.json",
                        dict(id=jid, status=status, request=dict(srx="SRX1"), created_at="2026-01-01"))
        restored = Jobs(self.jobs.root, self.resources)
        try:
            self.assertEqual(len(restored.list()), 3)
            self.assertTrue(all(job["status"] == "paused" for job in restored.list()))
        finally:
            restored.close()

    def test_failed_job_can_resume_with_same_id(self):
        ready(self.resources)
        with patch.object(self.jobs, "_run", side_effect=RuntimeError("fixture failure")):
            job = self.jobs.submit(dict(srx="SRX1"))
            wait_status(self.jobs, job["id"], "failed")
        self.jobs.resume(job["id"])
        result = wait_status(self.jobs, job["id"], "completed")
        self.assertEqual(result["id"], job["id"])


class HTTPTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.resources = self.root / "resources"
        ready(self.resources)
        self.jobs = Jobs(self.root / "jobs", self.resources)
        # A completed persisted job avoids external bioinformatics tools in API tests.
        self.jid = "a" * 32
        self.directory = self.jobs.root / self.jid
        sample(self.directory / "output" / "SRX1", "SRX1")
        (self.directory / "pipeline.log").write_text("解析ログ\n", encoding="utf-8")
        self.jobs.items[self.jid] = dict(id=self.jid, status="completed", created_at="2026-01-01",
                                       request=dict(srx=["SRX1"]))
        self.server = ThreadingHTTPServer(("127.0.0.1", 0), handler(self.jobs))
        self.thread = threading.Thread(target=self.server.serve_forever, kwargs={"poll_interval": 0.01})
        self.thread.start()

    def tearDown(self):
        self.server.shutdown()
        self.server.server_close()
        self.thread.join()
        self.jobs.close()
        self.temporary.cleanup()

    def request(self, method, path, payload=None, headers=None):
        connection = HTTPConnection("127.0.0.1", self.server.server_port, timeout=5)
        body = None if payload is None else json.dumps(payload).encode("utf-8")
        combined = {"Content-Type": "application/json"}
        combined.update(headers or {})
        connection.request(method, path, body, combined)
        response = connection.getresponse()
        status, data, response_headers = response.status, response.read(), dict(response.getheaders())
        connection.close()
        return status, data, response_headers

    def test_health_and_persisted_jobs(self):
        status, data, _ = self.request("GET", "/api/health")
        self.assertEqual(status, 200)
        health = json.loads(data)
        self.assertEqual(health["app"], "mouse-strain-checker")
        self.assertTrue(health["reference_ready"])
        self.assertIn("hisat2", health["tools"])
        self.assertEqual(json.loads(self.request("GET", "/api/jobs")[1])[0]["id"], self.jid)

    def test_cross_origin_and_non_json_requests_are_rejected(self):
        self.assertEqual(self.request("POST", "/api/jobs", dict(srx="SRX1"),
                                      {"Origin": "http://attacker.example"})[0], 403)
        self.assertEqual(self.request("POST", "/api/jobs", dict(srx="SRX1"),
                                      {"Sec-Fetch-Site": "cross-site"})[0], 403)
        self.assertEqual(self.request("POST", "/api/jobs", dict(srx="SRX1"),
                                      {"Content-Type": "text/plain"})[0], 415)
        self.assertEqual(self.request("POST", "/api/jobs", dict(srx="x" * 17000))[0], 413)

    def test_same_origin_json_submission_and_actions(self):
        origin = f"http://127.0.0.1:{self.server.server_port}"
        with patch.object(self.jobs, "submit", return_value={"id": self.jid}) as submit:
            status, data, _ = self.request("POST", "/api/jobs", dict(srx="SRX1"), {"Origin": origin})
            self.assertEqual(status, 202)
            self.assertEqual(json.loads(data)["id"], self.jid)
            submit.assert_called_once_with(dict(srx="SRX1"))
        with patch.object(self.jobs, "pause", return_value={"status": "paused"}) as pause:
            status, data, _ = self.request("POST", f"/api/jobs/{self.jid}/pause", {})
            self.assertEqual(status, 202)
            pause.assert_called_once_with(self.jid)
        self.assertEqual(self.request("POST", "/api/jobs/" + "b" * 32 + "/resume", {})[0], 404)

    def test_results_and_downloads_are_limited_to_completed_allowlisted_files(self):
        path = f"/api/jobs/{self.jid}"
        status, data, _ = self.request("GET", path + "/results")
        self.assertEqual(status, 200)
        self.assertEqual(json.loads(data)["samples"][0]["accession"], "SRX1")
        status, data, headers = self.request("GET", path + "/files/SRX1/chromosomes.csv")
        self.assertEqual(status, 200)
        self.assertIn(b"C57BL_6J", data)
        self.assertIn("attachment", headers["Content-Disposition"])
        for suffix in ("/files/SRX1/cache/bam", "/files/SRX1/../job.json", "/files/SRX2/complete.json",
                       "/files/SRX1/%2e%2e%2fjob.json", "/files/job.json"):
            self.assertEqual(self.request("GET", path + suffix)[0], 404)
        self.jobs.items[self.jid]["status"] = "paused"
        self.assertEqual(self.request("GET", path + "/results")[0], 409)
        self.assertEqual(self.request("GET", path + "/files/SRX1/chromosomes.csv")[0], 404)
        self.assertEqual(self.request("GET", path + "/files/pipeline.log")[0], 200)

    def test_table_response_is_bounded_but_full_csv_can_be_downloaded(self):
        file = self.directory / "output" / "SRX1" / "windows.csv"
        with file.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(["chromosome", "start"])
            writer.writerows(("1", index) for index in range(1005))
        path = f"/api/jobs/{self.jid}"
        result = json.loads(self.request("GET", path + "/results")[1])
        self.assertTrue(result["truncated"])
        self.assertEqual(result["samples"][0]["counts"]["windows"], 1005)
        self.assertEqual(len(result["samples"][0]["windows"]), 1000)
        data = self.request("GET", path + "/files/SRX1/windows.csv")[1]
        self.assertEqual(len(data.splitlines()), 1006)

    def test_log_tail_and_html_are_served_as_safe_downloads(self):
        path = f"/api/jobs/{self.jid}"
        self.assertIn("解析ログ", json.loads(self.request("GET", path + "/log")[1])["text"])
        status, data, headers = self.request("GET", path + "/files/SRX1/report.html")
        self.assertEqual(status, 200)
        self.assertEqual(headers["Content-Type"], "application/octet-stream")
        self.assertIn("style-src-attr 'unsafe-inline'", headers["Content-Security-Policy"])


if __name__ == "__main__":
    unittest.main()
