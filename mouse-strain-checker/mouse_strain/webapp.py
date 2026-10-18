"""Local browser application, implemented with the Python standard library."""
from __future__ import annotations

import csv
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
import json
import mimetypes
import os
from pathlib import Path
import re
import shutil
import signal
import threading
from urllib.parse import urlsplit

from .web_jobs import RESULT_FILES, Jobs, reference_ready

STATIC = Path(__file__).with_name("static")
MAX_BODY = 16384
TABLE_LIMIT = 1000


def result_data(jobs, job):
    samples, truncated = [], False
    for accession in job["request"]["srx"]:
        if not jobs.sample_complete(job["id"], accession):
            raise ValueError(f"{accession}の結果が揃っていません。")
        directory = jobs.root / job["id"] / "output" / accession
        sample = dict(accession=accession, counts={})
        for key, filename in (("chromosomes", "chromosomes.csv"), ("windows", "windows.csv"),
                              ("scores", "candidate_scores.csv")):
            rows, count = [], 0
            file = directory / filename
            if not file.resolve().is_relative_to(directory.resolve()):
                raise ValueError("結果ファイルのパスが不正です。")
            with file.open(encoding="utf-8-sig", newline="") as handle:
                for row in csv.DictReader(handle):
                    count += 1
                    if len(rows) < TABLE_LIMIT:
                        rows.append(row)
            sample[key], sample["counts"][key] = rows, count
            truncated = truncated or count > TABLE_LIMIT
        sample["files"] = sorted(name for name in RESULT_FILES if (directory / name).is_file()
                                 and (directory / name).resolve().is_relative_to(directory.resolve()))
        samples.append(sample)
    return dict(samples=samples, truncated=truncated)


def handler(jobs):
    class Handler(BaseHTTPRequestHandler):
        server_version = "MouseStrainChecker"

        def setup(self):
            super().setup()
            self.connection.settimeout(15)

        def log_message(self, format, *args):
            # Avoid logging accessions and request parameters to the server console.
            pass

        def _headers(self, status, content_type, length, filename=None):
            self.send_response(status)
            self.send_header("Content-Type", content_type)
            self.send_header("Content-Length", str(length))
            self.send_header("Cache-Control", "no-store")
            self.send_header("X-Content-Type-Options", "nosniff")
            self.send_header("Referrer-Policy", "same-origin")
            self.send_header("Content-Security-Policy", "default-src 'self'; script-src 'self'; "
                             "style-src 'self'; style-src-attr 'unsafe-inline'; img-src 'self' data:; "
                             "object-src 'none'; base-uri 'none'; frame-ancestors 'none'")
            if filename:
                self.send_header("Content-Disposition", f'attachment; filename="{filename}"')
            self.end_headers()

        def send(self, status, data, content_type="application/json; charset=utf-8"):
            if not isinstance(data, bytes):
                data = json.dumps(data, ensure_ascii=False, allow_nan=False).encode("utf-8")
            self._headers(status, content_type, len(data))
            try:
                self.wfile.write(data)
            except (BrokenPipeError, ConnectionResetError):
                pass

        def _file(self, file, base, filename):
            if not file.resolve().is_relative_to(base.resolve()) or not file.is_file():
                return self.send(404, {"error": "ファイルが見つかりません。"})
            try:
                with file.open("rb") as handle:
                    self._headers(200, "application/octet-stream", os.fstat(handle.fileno()).st_size, filename)
                    shutil.copyfileobj(handle, self.wfile, length=1024 * 1024)
            except (BrokenPipeError, ConnectionResetError):
                pass

        def do_POST(self):
            origin = self.headers.get("Origin")
            if origin:
                try:
                    parsed_origin = urlsplit(origin)
                except ValueError:
                    return self.send(403, {"error": "同じ画面から操作してください。"})
                if parsed_origin.scheme != "http" or parsed_origin.netloc != self.headers.get("Host"):
                    return self.send(403, {"error": "同じ画面から操作してください。"})
            if self.headers.get("Sec-Fetch-Site") == "cross-site":
                return self.send(403, {"error": "同じ画面から操作してください。"})
            path = urlsplit(self.path).path
            action = re.fullmatch(r"/api/jobs/([a-f0-9]{32})/(pause|resume)", path)
            if path != "/api/jobs" and not action:
                return self.send(404, {"error": "Not found"})
            if self.headers.get("Content-Type", "").split(";")[0].strip().lower() != "application/json":
                return self.send(415, {"error": "application/jsonを指定してください。"})
            try:
                length = int(self.headers.get("Content-Length", "0"))
                if not 0 < length <= MAX_BODY or self.headers.get("Transfer-Encoding"):
                    return self.send(413, {"error": "入力サイズが不正です。"})
                payload = json.loads(self.rfile.read(length))
                if action:
                    jid, operation = action.groups()
                    if not jobs.get(jid):
                        return self.send(404, {"error": "ジョブが見つかりません。"})
                    if payload != {}:
                        raise ValueError("操作の入力は空のJSONオブジェクトを指定してください。")
                    job = jobs.pause(jid) if operation == "pause" else jobs.resume(jid)
                else:
                    job = jobs.submit(payload)
            except (ValueError, UnicodeError, TypeError) as error:
                return self.send(400, {"error": str(error)})
            except OSError:
                return self.send(500, {"error": "ジョブ情報を保存できませんでした。空き容量を確認してください。"})
            return self.send(202, job)

        def do_GET(self):
            try:
                self._get()
            except (OSError, UnicodeError, csv.Error, ValueError):
                self.send(500, {"error": "保存ファイルを読み取れませんでした。ログを確認してください。"})

        def _get(self):
            path = urlsplit(self.path).path
            if path in {"/", "/app.js", "/style.css"}:
                file = STATIC / ("index.html" if path == "/" else path[1:])
                if not file.is_file():
                    return self.send(404, {"error": "画面ファイルが見つかりません。"})
                mime = mimetypes.guess_type(file.name)[0] or "application/octet-stream"
                return self.send(200, file.read_bytes(), mime + "; charset=utf-8")
            if path in {"/api/health", "/_stcore/health"}:
                tools = ("prefetch", "vdb-validate", "fasterq-dump", "hisat2", "hisat2-build", "samtools")
                return self.send(200, dict(app="mouse-strain-checker", tools={name: bool(shutil.which(name)) for name in tools},
                                          default_threads=jobs.threads, logical_cpus=os.cpu_count(),
                                          disk_free_gb=round(shutil.disk_usage(jobs.root).free / 1024**3, 1),
                                          reference_ready=bool(reference_ready(jobs.resources))))
            if path == "/api/jobs":
                return self.send(200, jobs.list())
            match = re.fullmatch(r"/api/jobs/([a-f0-9]{32})(?:/(results|log|files/.+))?", path)
            if not match:
                return self.send(404, {"error": "Not found"})
            jid, action = match.groups()
            job = jobs.get(jid)
            if not job:
                return self.send(404, {"error": "ジョブが見つかりません。"})
            directory = jobs.root / jid
            if not action:
                return self.send(200, job)
            if action == "log":
                file = directory / "pipeline.log"
                if not file.is_file() or not file.resolve().is_relative_to(directory):
                    return self.send(200, {"text": ""})
                with file.open("rb") as handle:
                    handle.seek(max(0, file.stat().st_size - 48000))
                    return self.send(200, {"text": handle.read().decode("utf-8", errors="replace")})
            if action == "results":
                if job["status"] != "completed":
                    return self.send(409, {"error": "解析はまだ完了していません。"})
                try:
                    return self.send(200, result_data(jobs, job))
                except ValueError as error:
                    return self.send(409, {"error": str(error)})
            if action == "files/pipeline.log":
                return self._file(directory / "pipeline.log", directory, "pipeline.log")
            file_match = re.fullmatch(r"files/([SED]RX[0-9]+)/([^/]+)", action)
            if not file_match:
                return self.send(404, {"error": "ファイルが見つかりません。"})
            accession, name = file_match.groups()
            if (job["status"] != "completed" or accession not in job["request"]["srx"]
                    or name not in RESULT_FILES or not jobs.sample_complete(jid, accession)):
                return self.send(404, {"error": "ファイルが見つかりません。"})
            return self._file(directory / "output" / accession / name, directory, accession + "_" + name)

    return Handler


def serve(host="127.0.0.1", port=8501, root="data/jobs", resources="data/resources", threads=4):
    jobs = Jobs(root, resources, threads)
    server = None
    previous_term = None
    try:
        server = ThreadingHTTPServer((host, port), handler(jobs))
        server.daemon_threads = True
        if threading.current_thread() is threading.main_thread():
            def terminate(signum, frame):
                raise KeyboardInterrupt
            previous_term = signal.signal(signal.SIGTERM, terminate)
        print(f"Mouse Strain Checker: http://{host}:{server.server_port}", flush=True)
        server.serve_forever(poll_interval=0.2)
    except KeyboardInterrupt:
        pass
    finally:
        if previous_term is not None:
            signal.signal(signal.SIGTERM, previous_term)
        if server:
            server.server_close()
        jobs.close()
