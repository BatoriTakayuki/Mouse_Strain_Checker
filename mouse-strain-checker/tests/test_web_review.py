"""UI process management and result visibility regression tests."""
import os
from pathlib import Path
import signal
import subprocess
import tempfile
import unittest
from unittest.mock import Mock, patch

from mouse_strain.jobs import completion_signature, launch, stop

try:
    from streamlit.testing.v1 import AppTest
except ImportError:
    AppTest = None


class JobTests(unittest.TestCase):
    def test_launch_uses_utf8_logs_and_separate_process_group(self):
        with tempfile.TemporaryDirectory() as temporary:
            with patch("mouse_strain.jobs.subprocess.Popen") as popen:
                launch(["python", "-m", "mouse_strain"], Path(temporary) / "job.log")
            options = popen.call_args.kwargs
            self.assertEqual(options["env"]["PYTHONIOENCODING"], "utf-8")
            self.assertEqual(options["env"]["PYTHONUNBUFFERED"], "1")
            self.assertTrue(options.get("start_new_session") or options.get("creationflags"))

    @unittest.skipUnless(os.name == "posix", "POSIX process groups required")
    def test_stop_escalates_after_timeout(self):
        process = Mock(pid=9876)
        process.poll.return_value = None
        process.wait.side_effect = [subprocess.TimeoutExpired("fixture", 10), 0]
        with patch("mouse_strain.jobs.os.killpg") as kill:
            stop(process)
        self.assertEqual(kill.call_args_list[0].args, (9876, signal.SIGTERM))
        self.assertEqual(kill.call_args_list[1].args, (9876, signal.SIGKILL))

    def test_windows_stop_targets_child_processes(self):
        process = Mock(pid=9876)
        process.poll.return_value = None
        with patch("mouse_strain.jobs.os.name", "nt"), \
             patch("mouse_strain.jobs.subprocess.run", return_value=Mock(returncode=0)) as run:
            stop(process)
        self.assertEqual(run.call_args.args[0], ["taskkill", "/PID", "9876", "/T", "/F"])
        process.wait.assert_called_once_with(timeout=10)

    @unittest.skipUnless(os.name == "posix", "POSIX process groups required")
    def test_stop_finishes_children_after_parent_exits(self):
        process = Mock(pid=9876)
        process.poll.return_value = None
        process.wait.return_value = 0
        with patch("mouse_strain.jobs.os.killpg") as kill:
            stop(process)
        self.assertEqual(kill.call_args_list[-1].args, (9876, signal.SIGKILL))


@unittest.skipIf(AppTest is None, "Streamlit is optional; install requirements-web.txt for UI tests")
class WebReviewTests(unittest.TestCase):
    def app(self):
        return AppTest.from_file(str(Path(__file__).resolve().parents[1] / "mouse_strain/app.py"))

    def test_launch_deduplicates_and_hides_previous_completion(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            resources = root / "resources"
            resources.mkdir()
            (resources / "config.json").write_text("{}", encoding="utf-8")
            output = root / "results"
            existing = output / "SRX1"
            existing.mkdir(parents=True)
            (existing / "complete.json").write_text("{}", encoding="utf-8")
            process = Mock()
            process.poll.return_value = None
            with patch("mouse_strain.jobs.launch", return_value=process) as start:
                app = self.app().run()
                app.text_input[0].set_value(str(resources))
                app.text_input[1].set_value(str(output))
                app.text_area[0].set_value("SRX1, SRX1")
                app.button[0].click().run()
            self.assertFalse(app.exception)
            command = start.call_args.args[0]
            self.assertEqual(command.count("SRX1"), 1)
            self.assertEqual(app.session_state["job"]["previous"]["SRX1"], completion_signature(existing))
            self.assertEqual(len(app.subheader), 0)
            self.assertEqual(len(app.dataframe), 0)

    def test_launch_error_is_shown_without_crashing(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "config.json").write_text("{}", encoding="utf-8")
            with patch("mouse_strain.jobs.launch", side_effect=OSError("fixture failure")):
                app = self.app().run()
                app.text_input[0].set_value(str(root))
                app.text_input[1].set_value(str(root / "out"))
                app.text_area[0].set_value("SRX1")
                app.button[0].click().run()
            self.assertFalse(app.exception)
            self.assertIn("fixture failure", app.error[0].value)


if __name__ == "__main__":
    unittest.main()
