"""Scoped storage cleanup with fresh previews and no traversal through links."""
from datetime import datetime, timezone
import errno
import hashlib
import json
import os
from pathlib import Path
import re
import secrets
import shutil
import stat
import time


ACTIVE = {"queued", "running", "pausing"}
TERMINAL = {"completed", "paused", "failed"}
SOURCE_NAMES = tuple(name + suffix for name in ("reference.fa.gz", "strains.vcf.gz")
                     for suffix in ("", ".download.json", ".partial", ".partial.json"))
MAX_ENTRIES = 100_000


def _guard(base, path):
    base, path = Path(base).resolve(), Path(path).absolute()
    if path == base or not path.is_relative_to(base) or not path.resolve().is_relative_to(base):
        raise ValueError("削除対象が許可された保存先の外にあります。")
    current = path
    while current != base:
        try:
            info = current.lstat()
        except FileNotFoundError:
            pass
        else:
            if (stat.S_ISLNK(info.st_mode) or getattr(info, "st_file_attributes", 0) & 0x400
                    or current.is_mount()):
                raise ValueError("リンク・ジャンクション・別のマウント先を含むため削除できません。")
        current = current.parent
    return path


def _key(info):
    return [stat.S_IFMT(info.st_mode), info.st_size, info.st_mtime_ns, info.st_ino, info.st_dev]


def _scan(base, target):
    target = _guard(base, target)
    if not target.exists():
        return []
    entries, pending = [], [target]
    while pending:
        # Descendants are reached only through checked directories. Rechecking
        # every ancestor for every file is costly on Windows bind mounts;
        # cleanup still checks the complete path immediately before unlinking.
        path = pending.pop()
        info = path.lstat()
        directory = stat.S_ISDIR(info.st_mode)
        if (stat.S_ISLNK(info.st_mode) or getattr(info, "st_file_attributes", 0) & 0x400
                or directory and path.is_mount()):
            raise ValueError("リンク・ジャンクション・別のマウント先を含むため削除できません。")
        if not directory and not stat.S_ISREG(info.st_mode):
            raise ValueError("通常のファイル以外を含むため削除できません。")
        entries.append(dict(path=path, directory=directory, key=_key(info)))
        if len(entries) > MAX_ENTRIES:
            raise ValueError("ファイルが多すぎます。保存先を手動で確認してください。")
        if directory:
            with os.scandir(path) as children:
                pending.extend(Path(child.path) for child in children)
    return entries


def _bytes(entries):
    return sum(entry["key"][1] for entry in entries if not entry["directory"])


def _atomic_json(base, path, data):
    path = _guard(base, path)
    temporary = _guard(base, path.with_name(path.name + ".partial"))
    temporary.write_text(json.dumps(data, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    temporary.replace(path)


def retained_sex(output):
    """Only reuse saved sex evidence tied to an intentional cache cleanup."""
    output = Path(output)
    try:
        cleanup, sex = output / "cleanup.json", output / "sex.json"
        _guard(output, cleanup)
        _guard(output, sex)
        if cleanup.stat().st_size > 65536 or sex.stat().st_size > 1024 * 1024:
            return None
        record = json.loads(cleanup.read_text(encoding="utf-8"))
        blob = sex.read_bytes()
        result = json.loads(blob)
        if (record.get("scope") != "cache" or record.get("sex_sha256") != hashlib.sha256(blob).hexdigest()
                or not isinstance(result, dict) or result.get("schema") != 1
                or result.get("status") not in {"male_likely", "female_likely", "inconclusive", "unavailable"}
                or not isinstance(result.get("label"), str) or not isinstance(result.get("reason"), str)):
            return None
        json.dumps(result, allow_nan=False)
        return result
    except (OSError, ValueError, TypeError, AttributeError):
        return None


class Storage:
    def __init__(self, jobs):
        self.jobs = jobs
        self.previews = {}

    def _job(self, jid):
        if not isinstance(jid, str) or not re.fullmatch(r"[a-f0-9]{32}", jid) or jid not in self.jobs.items:
            raise ValueError("解析が見つかりません。容量表示を更新してください。")
        job = self.jobs.items[jid]
        if job["status"] not in TERMINAL or jid in self.jobs._processes:
            raise ValueError("実行中・待機中・中断処理中の解析は削除できません。")
        return job

    def _build(self, scope, jid=None, scanned_job=None):
        if self.jobs._closing:
            raise ValueError("サーバーを終了しています。")
        if not isinstance(scope, str) or scope not in {"cache", "job", "reference_downloads"}:
            raise ValueError("削除方法が不正です。")
        roots, entries, selected, samples = [], [], [], []
        if scope == "reference_downloads":
            if jid is not None or any(j["status"] in ACTIVE for j in self.jobs.items.values()):
                raise ValueError("参照の元ファイルは、実行中・待機中の解析がないときに削除できます。")
            base = self.jobs.resources
            config = json.loads(_guard(base, base / "config.json").read_text(encoding="utf-8"))
            if not isinstance(config, dict) or any(not isinstance(config.get(k), str) or not config[k] for k in ("fasta", "panel", "index")):
                raise ValueError("参照データの設定が不正です。")
            fasta = _guard(base, base / config["fasta"])
            panel = _guard(base, base / config["panel"])
            index = _guard(base, base / config["index"])
            required = [fasta, Path(str(fasta) + ".fai"), panel / "metadata.json", panel / "panel.sqlite", panel / "sites.bed"]
            required = [_guard(base, p) for p in required]
            complete_index = any(all(_guard(base, Path(f"{index}.{i}.{ext}")).is_file()
                                     and Path(f"{index}.{i}.{ext}").stat().st_size
                                     for i in range(1, 9)) for ext in ("ht2", "ht2l"))
            if not complete_index or not all(p.is_file() and p.stat().st_size for p in required):
                raise ValueError("参照の準備が完了していないため、元ファイルを削除できません。")
            roots = [base / name for name in SOURCE_NAMES]
            if any(path.resolve() in {p.resolve() for p in required} for path in roots):
                raise ValueError("解析に必要な参照ファイルは削除できません。")
            for root in roots:
                part = _scan(base, root)
                if any(entry["directory"] for entry in part):
                    raise ValueError("参照の元ファイルと同名のフォルダがあるため削除できません。")
                entries.extend(part)
            selected = [entry for entry in entries if not entry["directory"]]
            title = "参照作成に使ったダウンロード元"
            effect = "VCF.gzと参照FASTA.gz、その転送記録を削除します。解析に必要な展開済みFASTA・SNPパネル・索引は保持します。参照を作り直す場合は再ダウンロードが必要です。"
        else:
            job = self._job(jid)
            base = self.jobs.root
            directory = _guard(base, base / jid)
            title = ", ".join(job["request"]["srx"]) or "参照準備の履歴"
            if scope == "job":
                roots = [directory]
                entries = _scan(base, directory)
                selected = [entry for entry in entries if not entry["directory"]]
                effect = "この解析の結果・ログ・中間ファイル・履歴をすべて削除します。元に戻せません。共有の参照データは保持します。"
            else:
                for accession in job["request"]["srx"]:
                    if not re.fullmatch(r"[SED]RX[0-9]+", accession):
                        raise ValueError("解析番号が不正です。")
                    sample = _guard(base, directory / "output" / accession)
                    cache = sample / "cache"
                    roots.append(cache)
                    part = (_scan(base, cache) if scanned_job is None else
                            [entry for entry in scanned_job if entry["path"].is_relative_to(cache)])
                    entries.extend(part)
                    if len(entries) > MAX_ENTRIES:
                        raise ValueError("ファイルが多すぎます。保存先を手動で確認してください。")
                    removable = []
                    for entry in part:
                        relative = entry["path"].relative_to(cache).parts
                        qc = (len(relative) == 3 and relative[0] == "preprocess"
                              and re.fullmatch(r"[SED]RR[0-9]+", relative[1])
                              and relative[2] in {"fastp.html", "fastp.json"})
                        if not entry["directory"] and not qc:
                            removable.append(entry)
                    selected.extend(removable)
                    if removable:
                        samples.append((sample, _bytes(removable)))
                effect = "FASTQ・SRA・BAM・一時ファイルを削除します。結果CSV/HTML・性別結果・観測SNP・fastpレポートは保持します。再開・再解析に必要なデータは再取得・再作成されます。"
        signature = [[str(entry["path"].relative_to(base)), entry["key"]] for entry in entries]
        digest = hashlib.sha256(json.dumps([scope, jid, sorted(map(str, roots)), sorted(signature)], ensure_ascii=False).encode()).hexdigest()
        return dict(scope=scope, job_id=jid, title=title, effect=effect, base=base, entries=entries,
                    selected=selected, samples=samples, fingerprint=digest, bytes=_bytes(selected), files=len(selected))

    def overview(self):
        with self.jobs.lock:
            rows = []
            for job in self.jobs.list():
                row = dict(id=job["id"], accessions=job["request"]["srx"], status=job["status"],
                           total_bytes=None, cache_bytes=None, deletable=False, reason="")
                try:
                    scanned = _scan(self.jobs.root, self.jobs.root / job["id"])
                    row["total_bytes"] = _bytes(scanned)
                    plan = self._build("cache", job["id"], scanned_job=scanned)
                    row.update(cache_bytes=plan["bytes"], deletable=True)
                except (OSError, ValueError, KeyError, TypeError) as error:
                    row["reason"] = str(error)
                rows.append(row)
            reference = dict(bytes=0, deletable=False, reason="")
            try:
                plan = self._build("reference_downloads")
                reference.update(bytes=plan["bytes"], deletable=bool(plan["files"]))
            except (OSError, ValueError, KeyError, TypeError) as error:
                reference["reason"] = str(error)
            usage = shutil.disk_usage(self.jobs.root)
            return dict(jobs=rows, reference=reference, free_bytes=usage.free, total_bytes=usage.total)

    def preview(self, scope, jid=None):
        with self.jobs.lock:
            self.previews = {k: v for k, v in self.previews.items() if v["expires"] > time.monotonic()}
            if len(self.previews) >= 128:
                self.previews.pop(next(iter(self.previews)))
            plan = self._build(scope, jid)
            if not plan["files"]:
                raise ValueError("削除対象のファイルはありません。")
            token = secrets.token_hex(32)
            self.previews[token] = dict(scope=scope, jid=jid, fingerprint=plan["fingerprint"], expires=time.monotonic() + 600)
            return dict(token=token, scope=scope, job_id=jid, title=plan["title"], effect=plan["effect"],
                        bytes=plan["bytes"], files=plan["files"],
                        paths=[str(entry["path"].relative_to(plan["base"])) for entry in plan["selected"][:15]])

    def cleanup(self, token, confirmed=False):
        with self.jobs.lock:
            if confirmed is not True or not isinstance(token, str):
                raise ValueError("削除対象を確認してから実行してください。")
            saved = self.previews.pop(token, None)
            if not saved or saved["expires"] <= time.monotonic():
                raise ValueError("削除確認の有効期限が切れました。容量を再確認してください。")
            plan = self._build(saved["scope"], saved["jid"])
            if plan["fingerprint"] != saved["fingerprint"]:
                raise ValueError("保存内容が変わりました。容量と削除対象を再確認してください。")
            stamp = datetime.now(timezone.utc).isoformat()
            for sample, size in plan["samples"]:
                sex = _guard(plan["base"], sample / "sex.json")
                digest = hashlib.sha256(sex.read_bytes()).hexdigest() if sex.is_file() and sex.stat().st_size <= 1024 * 1024 else None
                _atomic_json(plan["base"], sample / "cleanup.json", dict(schema=1, scope="cache", at=stamp,
                    bytes_selected=size, sex_sha256=digest, note="中間ファイルを削除しました。結果と性別評価は保存済みの解析結果です。"))
            # Keep the history record until both files and child directories
            # have been removed, so an I/O failure remains visible for retrying.
            jid = plan["job_id"]
            history = plan["base"] / jid / "job.json" if plan["scope"] == "job" else None
            files = sorted(plan["selected"], key=lambda e: str(e["path"]))
            for entry in files:
                if entry["path"] == history:
                    continue
                path = _guard(plan["base"], entry["path"])
                if _key(path.lstat()) != entry["key"]:
                    raise ValueError("削除中にファイルが変わりました。残りの対象を再確認してください。")
                path.unlink()
            directories = sorted([entry["path"] for entry in plan["entries"] if entry["directory"]], key=lambda p: len(p.parts), reverse=True)
            for path in directories:
                if history is not None and path == history.parent:
                    continue
                try:
                    _guard(plan["base"], path).rmdir()
                except OSError as error:
                    if plan["scope"] != "cache" or error.errno not in (errno.ENOTEMPTY, errno.EEXIST):
                        raise
            if plan["scope"] == "job":
                if history.exists():
                    entry = next(e for e in files if e["path"] == history)
                    if _key(_guard(plan["base"], history).lstat()) != entry["key"]:
                        raise ValueError("削除中に履歴が変わりました。残りの対象を再確認してください。")
                    history.unlink()
                try:
                    _guard(plan["base"], history.parent).rmdir()
                except OSError:
                    _atomic_json(plan["base"], history, self.jobs.items[jid])
                    raise
                self.jobs.items.pop(jid)
                self.jobs._events.pop(jid, None)
                self.jobs._generation.pop(jid, None)
            elif jid:
                self.jobs._update(jid, last_cleanup=dict(scope="cache", at=stamp, bytes=plan["bytes"]))
            return dict(scope=plan["scope"], job_id=jid, bytes=plan["bytes"], files=plan["files"],
                        free_bytes=shutil.disk_usage(self.jobs.root).free,
                        message="削除が完了しました。" if plan["scope"] == "job" else "対象ファイルを削除しました。解析結果は保持しています。")
