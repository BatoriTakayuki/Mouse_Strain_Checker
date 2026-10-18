"""Local HTTP transfers and bounded recovery of interrupted setup artifacts."""
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
import json
from pathlib import Path
import tempfile
import threading
import unittest
from unittest.mock import patch

from mouse_strain.__main__ import _clear_unfinished_index, _clear_unfinished_panel, main
from mouse_strain.pipeline import download


class Paused(Exception):
    pass


class ResumeTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.target = Path(self.temp.name) / "reference.gz"
        self.data = b"ACGT" * 200000
        self.etag = '"v1"'
        self.ignore_range = False
        self.bad_range = False
        self.truncate_once = False
        self.requests = []
        owner = self

        class Handler(BaseHTTPRequestHandler):
            def log_message(self, *args):
                pass

            def do_GET(self):
                owner.requests.append(dict(self.headers))
                requested = self.headers.get("Range")
                offset = int(requested.split("=")[1].split("-")[0]) if requested else 0
                ranged = requested and not owner.ignore_range and self.headers.get("If-Range") == owner.etag
                if ranged and offset >= len(owner.data):
                    self.send_response(416)
                    self.end_headers()
                    return
                start = offset if ranged else 0
                body = owner.data[start:]
                self.send_response(206 if ranged else 200)
                self.send_header("Content-Length", str(len(body)))
                self.send_header("ETag", owner.etag)
                if ranged:
                    actual = 0 if owner.bad_range else start
                    self.send_header("Content-Range", f"bytes {actual}-{len(owner.data) - 1}/{len(owner.data)}")
                self.end_headers()
                if owner.truncate_once:
                    owner.truncate_once = False
                    body = body[:len(body) // 2]
                try:
                    self.wfile.write(body)
                except (BrokenPipeError, ConnectionResetError):
                    pass

        self.server = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
        self.thread = threading.Thread(target=self.server.serve_forever, daemon=True)
        self.thread.start()
        self.addCleanup(self.close_server)
        self.url = f"http://127.0.0.1:{self.server.server_port}/reference.gz"

    def close_server(self):
        self.server.shutdown()
        self.server.server_close()
        self.thread.join()

    def interrupt(self):
        def progress(**details):
            if details["done_bytes"] >= 128000:
                raise Paused()
        with patch("mouse_strain.pipeline.transfer", side_effect=progress):
            with self.assertRaises(Paused):
                download(self.url, self.target)
        self.assertFalse(self.target.exists())
        return Path(str(self.target) + ".partial").stat().st_size

    def test_interrupted_download_resumes_with_validator_and_range(self):
        offset = self.interrupt()
        download(self.url, self.target)
        self.assertEqual(self.target.read_bytes(), self.data)
        self.assertEqual(self.requests[-1]["Range"], f"bytes={offset}-")
        self.assertEqual(self.requests[-1]["If-Range"], self.etag)
        self.assertFalse(Path(str(self.target) + ".partial.json").exists())

    def test_range_ignored_replaces_old_prefix(self):
        self.interrupt()
        self.ignore_range = True
        download(self.url, self.target)
        self.assertEqual(self.target.read_bytes(), self.data)

    def test_changed_representation_is_retrieved_fresh(self):
        self.interrupt()
        self.data = b"TGCA" * 150000
        self.etag = '"v2"'
        download(self.url, self.target)
        self.assertEqual(self.target.read_bytes(), self.data)

    def test_bad_content_range_retries_without_appending(self):
        self.interrupt()
        self.bad_range = True
        with patch("mouse_strain.pipeline.time.sleep"):
            download(self.url, self.target)
        self.assertEqual(self.target.read_bytes(), self.data)
        self.assertNotIn("Range", self.requests[-1])

    def test_short_http_transfer_retries_from_validated_prefix(self):
        self.truncate_once = True
        with patch("mouse_strain.pipeline.time.sleep"):
            download(self.url, self.target)
        self.assertEqual(self.target.read_bytes(), self.data)
        self.assertIn("Range", self.requests[-1])

    def test_416_refetches_without_trusting_partial_completion(self):
        self.interrupt()
        partial = Path(str(self.target) + ".partial")
        partial.write_bytes(self.data)
        with patch("mouse_strain.pipeline.time.sleep"):
            download(self.url, self.target)
        self.assertEqual(self.target.read_bytes(), self.data)
        self.assertNotIn("Range", self.requests[-1])


class SetupRecoveryTests(unittest.TestCase):
    def test_unfinished_panel_cleanup_preserves_unknown_files(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary).resolve()
            panel = root / "panel"
            panel.mkdir()
            (panel / "panel.sqlite.partial").write_bytes(b"partial")
            (panel / "notes.txt").write_text("keep", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "unrecognized"):
                _clear_unfinished_panel(root)
            self.assertTrue((panel / "panel.sqlite.partial").exists())
            self.assertEqual((panel / "notes.txt").read_text(), "keep")
            (panel / "notes.txt").unlink()
            _clear_unfinished_panel(root)
            self.assertFalse(panel.exists())

    def test_unfinished_index_cleanup_is_limited_to_generated_parts(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary).resolve()
            index = root / "hisat2"
            index.mkdir()
            (index / "genome.1.ht2").write_bytes(b"partial")
            (index / "custom.1.ht2").write_bytes(b"keep")
            _clear_unfinished_index(root)
            self.assertFalse((index / "genome.1.ht2").exists())
            self.assertEqual((index / "custom.1.ht2").read_bytes(), b"keep")

    def test_setup_resume_rejects_changed_recorded_inputs(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "setup_inputs.json").write_text(json.dumps({"changed": True}), encoding="utf-8")
            with patch("mouse_strain.__main__.require_tools"), patch("mouse_strain.__main__.download") as fetch:
                self.assertEqual(main(["setup", "--resume", "--out", str(root)]), 1)
            fetch.assert_not_called()


if __name__ == "__main__":
    unittest.main()
