"""Verify cancellation against real CLI child process groups on Linux."""
import os
from pathlib import Path
import signal
import sys
import tempfile
import time
import unittest
from unittest.mock import patch

from mouse_strain.web_jobs import Jobs


@unittest.skipUnless(sys.platform == "linux", "Linux process-group integration test")
class RealProcessTests(unittest.TestCase):
    def test_pause_kills_even_a_descendant_ignoring_sigterm(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            marker = root / "child.pid"
            child = ("import os,signal,time; from pathlib import Path; "
                     "signal.signal(signal.SIGTERM,signal.SIG_IGN); "
                     f"Path({str(marker)!r}).write_text(str(os.getpid())); time.sleep(60)")
            parent = f"import subprocess,sys,time; subprocess.Popen([sys.executable,'-c',{child!r}]); time.sleep(60)"

            def run_fixture(manager, jid):
                manager._command(jid, [sys.executable, "-c", parent], "mapping", "fixture")

            pid = None
            with patch.object(Jobs, "_run", run_fixture):
                jobs = Jobs(root / "jobs", root / "resources")
                try:
                    job = jobs.submit(dict(srx="SRX1"))
                    deadline = time.monotonic() + 5
                    while not marker.exists() and time.monotonic() < deadline:
                        time.sleep(0.02)
                    self.assertTrue(marker.exists(), "child did not start")
                    pid = int(marker.read_text())
                    jobs.pause(job["id"])
                    while jobs.get(job["id"])["status"] != "paused" and time.monotonic() < deadline:
                        time.sleep(0.02)
                    self.assertEqual(jobs.get(job["id"])["status"], "paused")
                    process_state = Path(f"/proc/{pid}/stat")
                    exit_deadline = time.monotonic() + 2
                    while True:
                        try:
                            state = process_state.read_text().split()[2]
                        except FileNotFoundError:
                            break
                        if state == "Z":
                            break
                        self.assertLess(time.monotonic(), exit_deadline, "descendant is still running")
                        time.sleep(0.01)
                finally:
                    jobs.close()
                    if pid:
                        try:
                            os.kill(pid, signal.SIGKILL)
                        except ProcessLookupError:
                            pass


if __name__ == "__main__":
    unittest.main()
