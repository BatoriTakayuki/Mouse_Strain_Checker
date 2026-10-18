"""Persistent, serial CLI jobs used by the local browser application."""
from __future__ import annotations

from copy import deepcopy
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import queue
import re
import sys
import threading
import traceback
import uuid

from .jobs import launch, stop

JOB_ID = re.compile(r"[a-f0-9]{32}")
ACCESSION = re.compile(r"[SED]RX[0-9]+")
RESULT_FILES = frozenset({"report.html", "chromosomes.csv", "windows.csv",
                          "candidate_scores.csv", "evidence.json", "sites.tsv.gz",
                          "complete.json", "ena_runs.json"})
REQUIRED_RESULTS = ("report.html", "chromosomes.csv", "windows.csv",
                    "candidate_scores.csv", "evidence.json")


def now():
    return datetime.now(timezone.utc).isoformat()


def atomic_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    temporary.replace(path)


def read_json(path):
    try:
        return json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, ValueError, UnicodeError):
        return None


def reference_ready(resources):
    """Setup writes this configuration only after validating all reference assets."""
    config = read_json(Path(resources) / "config.json")
    return isinstance(config, dict) and all(
        isinstance(config.get(key), str) and config[key].strip()
        for key in ("fasta", "panel", "index"))


def validate_request(payload, default_threads=4):
    if not isinstance(payload, dict):
        raise ValueError("入力はJSONオブジェクトで指定してください。")
    allowed = {"mode", "srx", "threads", "single_only", "pause_after_reference"}
    if set(payload) - allowed:
        raise ValueError("使用できない設定が含まれています。")
    mode = payload.get("mode", "analyze")
    if mode not in ("analyze", "setup"):
        raise ValueError("modeはanalyzeまたはsetupを指定してください。")
    threads = payload.get("threads", default_threads)
    if type(threads) is not int or not 1 <= threads <= 128:
        raise ValueError("CPUスレッド数は1～128の整数で指定してください。")
    for name in ("single_only", "pause_after_reference"):
        if name in payload and type(payload[name]) is not bool:
            raise ValueError(f"{name}はtrueまたはfalseで指定してください。")
    raw = payload.get("srx", [])
    if isinstance(raw, str):
        raw = [value for value in re.split(r"[\s,]+", raw.strip()) if value]
    if not isinstance(raw, list) or any(not isinstance(value, str) for value in raw):
        raise ValueError("SRX番号は文字列または文字列の配列で指定してください。")
    accessions = list(dict.fromkeys(raw))
    if (mode == "analyze" and not accessions) or len(accessions) > 30:
        raise ValueError("解析するSRX/ERX/DRX番号を1～30件指定してください。")
    if any(not ACCESSION.fullmatch(value) or len(value) > 32 for value in accessions):
        raise ValueError("SRX/ERX/DRXに半角数字が続く番号を指定してください。")
    return dict(mode=mode, srx=accessions, threads=threads,
                single_only=payload.get("single_only", False),
                pause_after_reference=payload.get("pause_after_reference", False))


class _Paused(Exception):
    pass


class Jobs:
    """Run one child at a time; all queue and lifecycle changes hold one lock."""

    def __init__(self, root, resources, threads=4):
        if type(threads) is not int or not 1 <= threads <= 128:
            raise ValueError("threads must be an integer between 1 and 128")
        self.root, self.resources = Path(root).resolve(), Path(resources).resolve()
        self.root.mkdir(parents=True, exist_ok=True)
        self.threads = threads
        self.lock = threading.RLock()
        self.queue = queue.Queue()
        self.items, self._events, self._generation, self._processes = {}, {}, {}, {}
        self._closing = False
        for path in self.root.glob("*/job.json"):
            if not JOB_ID.fullmatch(path.parent.name) or not path.resolve().is_relative_to(self.root):
                continue
            job = read_json(path)
            if not isinstance(job, dict) or job.get("id") != path.parent.name:
                continue
            try:
                job["request"] = validate_request(job["request"], threads)
            except (KeyError, ValueError, TypeError):
                continue
            if job.get("status") not in {"queued", "running", "pausing", "paused", "failed", "completed"}:
                continue
            if job["status"] in {"queued", "running", "pausing"}:
                job.update(status="paused", pause_requested=False, paused_at=now(), updated_at=now(),
                           message="サーバーの再起動により中断しました。「再開」で保存済みデータから続けられます。")
                atomic_json(path, job)
            self.items[job["id"]] = job
            self._events[job["id"]] = threading.Event()
            self._generation[job["id"]] = 0
        self._worker = threading.Thread(target=self._work, name="strain-jobs", daemon=True)
        self._worker.start()

    def get(self, jid):
        with self.lock:
            return deepcopy(self.items.get(jid))

    def list(self):
        with self.lock:
            return sorted((deepcopy(job) for job in self.items.values()),
                          key=lambda job: job.get("created_at", ""), reverse=True)

    def _update(self, jid, **changes):
        with self.lock:
            job = self.items[jid]
            if "stage" in changes and changes["stage"] != job.get("stage"):
                job.update(stage_started_at=now(), transfer=None, activity="")
            job.update(changes, updated_at=now())
            atomic_json(self.root / jid / "job.json", job)

    def _enqueue(self, jid):
        self._generation[jid] = self._generation.get(jid, 0) + 1
        self.queue.put((jid, self._generation[jid]))

    def _capacity(self):
        if self._closing:
            raise ValueError("サーバーを終了しています。")
        if sum(job["status"] in {"queued", "running", "pausing"} for job in self.items.values()) >= 10:
            raise ValueError("待機中のジョブが多いため、完了後にもう一度実行してください。")

    def submit(self, payload):
        request = validate_request(payload, self.threads)
        with self.lock:
            self._capacity()
            jid = uuid.uuid4().hex
            job = dict(id=jid, request=request, accessions=request["srx"], threads=request["threads"],
                       status="queued", stage="queued", message="実行待ち", progress=0,
                       pause_requested=False, reference_ready=False, reference_pause_done=False,
                       completed_accessions=[], created_at=now(), updated_at=now())
            atomic_json(self.root / jid / "job.json", job)
            self.items[jid] = job
            self._events[jid] = threading.Event()
            self._enqueue(jid)
            return deepcopy(job)

    def pause(self, jid):
        with self.lock:
            job = self.items[jid]
            if job["status"] in {"paused", "pausing"}:
                return deepcopy(job)
            if job["status"] not in {"queued", "running"}:
                raise ValueError("実行中・待機中のジョブだけを中断できます。")
            self._events[jid].set()
            if job["status"] == "queued":
                # Invalidate the old queue token before allowing another resume.
                self._generation[jid] += 1
                self._update(jid, status="paused", pause_requested=False, paused_at=now(),
                             message="実行前に中断しました。")
            else:
                self._update(jid, status="pausing", pause_requested=True,
                             message="処理を停止しています。中断完了までお待ちください。")
            return self.get(jid)

    def resume(self, jid):
        with self.lock:
            if self.items[jid]["status"] not in {"paused", "failed"}:
                raise ValueError("中断中・失敗したジョブだけを再開できます。")
            self._capacity()
            self._events[jid] = threading.Event()
            self._update(jid, status="queued", pause_requested=False, error=None,
                         message="再開待ち。完了済みの解析と検証済みキャッシュを再利用します。")
            self._enqueue(jid)
            return self.get(jid)

    def _check_pause(self, jid):
        if self._closing or self._events[jid].is_set():
            raise _Paused("中断しました。「再開」で保存済みデータから続けられます。")

    def _read_progress(self, jid):
        value = read_json(self.root / jid / "progress.json")
        if not isinstance(value, dict):
            return
        allowed = {"stage", "message", "progress", "activity", "active_command", "transfer",
                   "storage_bytes", "storage_label", "accession", "current", "total",
                   "current_accession", "runs_done", "runs_total"}
        changes = {key: value[key] for key in allowed if key in value}
        # A malformed/stale file must never change the persisted job's identity or status.
        if changes and changes != {key: self.get(jid).get(key) for key in changes}:
            self._update(jid, **changes)

    def _command(self, jid, command, stage, message):
        directory = self.root / jid
        progress = directory / "progress.json"
        with self.lock:
            self._check_pause(jid)
            progress.unlink(missing_ok=True)
            self._update(jid, stage=stage, message=message, progress=None)
            with (directory / "pipeline.log").open("a", encoding="utf-8") as handle:
                handle.write(f"\n[{now()}] {' '.join(command)}\n")
            process = launch(command, directory / "pipeline.log", append=True,
                             env={"MOUSE_STRAIN_PROGRESS_FILE": str(progress)})
            self._processes[jid] = process
        try:
            while process.poll() is None:
                self._read_progress(jid)
                if self._closing or self._events[jid].wait(0.2):
                    stop(process)
                    self._check_pause(jid)
            self._read_progress(jid)
            self._check_pause(jid)
            if process.returncode:
                raise RuntimeError(f"処理が終了コード{process.returncode}で失敗しました。ログを確認してください。")
        finally:
            try:
                if process.poll() is None:
                    stop(process)
            finally:
                with self.lock:
                    self._processes.pop(jid, None)

    def sample_complete(self, jid, accession):
        if not ACCESSION.fullmatch(accession):
            return False
        directory = self.root / jid / "output" / accession
        if not directory.resolve().is_relative_to(self.root / jid):
            return False
        complete = read_json(directory / "complete.json")
        return (isinstance(complete, dict) and complete.get("accession") == accession
                and all((directory / name).is_file() for name in REQUIRED_RESULTS))

    def _run(self, jid):
        request = self.get(jid)["request"]
        threads = str(request["threads"])
        if not reference_ready(self.resources):
            self._command(jid, [sys.executable, "-m", "mouse_strain", "setup", "--out",
                               str(self.resources), "--threads", threads, "--resume"],
                          "setup", "参照データを準備しています。初回はダウンロードと計算に時間がかかります。")
            if not reference_ready(self.resources):
                raise RuntimeError("参照データの準備完了を確認できませんでした。")
        self._check_pause(jid)
        self._update(jid, reference_ready=True)
        if request["mode"] == "setup":
            self._update(jid, status="completed", stage="completed", progress=100,
                         message="参照データの準備が完了しました。", completed_at=now())
            return
        if request["pause_after_reference"] and not self.get(jid).get("reference_pause_done"):
            self._update(jid, reference_pause_done=True)
            raise _Paused("参照データの準備が完了しました。「再開」でSRXの解析を開始します。")
        completed = []
        for accession in request["srx"]:
            self._check_pause(jid)
            if not self.sample_complete(jid, accession):
                self._update(jid, accession=accession, completed_accessions=completed)
                command = [sys.executable, "-m", "mouse_strain", "run", "--resources",
                           str(self.resources), "--out", str(self.root / jid / "output"),
                           "--threads", threads, "--srx", accession]
                if request["single_only"]:
                    command.append("--single-only")
                self._command(jid, command, "analyze", f"{accession}を解析しています。")
                if not self.sample_complete(jid, accession):
                    raise RuntimeError(f"{accession}の解析完了ファイルを確認できませんでした。")
            completed.append(accession)
            self._update(jid, completed_accessions=list(completed))
        self._check_pause(jid)
        self._update(jid, status="completed", stage="completed", message="解析が完了しました。",
                     progress=100, pause_requested=False, completed_at=now())

    def _work(self):
        while True:
            token = self.queue.get()
            try:
                if token is None:
                    return
                jid, generation = token
                with self.lock:
                    if (self._closing or self._generation.get(jid) != generation
                            or self.items[jid]["status"] != "queued"):
                        continue
                    self._update(jid, status="running", attempt_started_at=now(), paused_at=None)
                try:
                    self._run(jid)
                except _Paused as error:
                    self._update(jid, status="paused", pause_requested=False, paused_at=now(),
                                 message=str(error))
                except Exception as error:
                    with (self.root / jid / "pipeline.log").open("a", encoding="utf-8") as handle:
                        traceback.print_exc(file=handle)
                    self._update(jid, status="failed", stage="failed", pause_requested=False,
                                 message=str(error), error=str(error))
            finally:
                self.queue.task_done()

    def close(self):
        with self.lock:
            if not self._closing:
                self._closing = True
                for jid, job in self.items.items():
                    if job["status"] in {"queued", "running", "pausing"}:
                        self.pause(jid)
                self.queue.put(None)
        if threading.current_thread() is not self._worker:
            self._worker.join(timeout=30)
            if self._worker.is_alive():
                raise RuntimeError("解析プロセスの停止を確認できませんでした。")
