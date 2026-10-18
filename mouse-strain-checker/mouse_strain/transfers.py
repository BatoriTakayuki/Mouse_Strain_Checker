"""One validated HTTP transfer attempt, retaining resumable partial files."""
import json
from pathlib import Path
import re
import time
import urllib.error
import urllib.request


def _header(response, name):
    value = response.headers.get(name)
    return value if isinstance(value, str) else None


def fetch(url, partial, notify):
    partial = Path(partial)
    metadata = Path(str(partial) + ".json")
    try:
        saved = json.loads(metadata.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        saved = {}
    if not isinstance(saved, dict):
        saved = {}
    validator = saved.get("validator")
    offset = partial.stat().st_size if partial.is_file() else 0
    if (saved.get("url") != url or not isinstance(validator, str) or not validator
            or saved.get("validator_header") not in ("ETag", "Last-Modified")):
        offset = 0
    headers = {"User-Agent": "Mouse-Strain-Checker/0.1", "Accept-Encoding": "identity"}
    if offset:
        headers.update(Range=f"bytes={offset}-", **{"If-Range": validator})
    try:
        response = urllib.request.urlopen(urllib.request.Request(url, headers=headers), timeout=30)
    except urllib.error.HTTPError as error:
        if error.code == 416:
            metadata.unlink(missing_ok=True)
        raise
    with response as source:
        status = getattr(source, "status", 200)
        status = status if isinstance(status, int) else 200
        total = None
        length = _header(source, "Content-Length")
        if length is not None:
            total = int(length)
        if status == 206:
            match = re.fullmatch(r"bytes ([0-9]+)-([0-9]+)/([0-9]+)", _header(source, "Content-Range") or "")
            if (not offset or not match or int(match[1]) != offset
                    or int(match[2]) != int(match[3]) - 1
                    or (total is not None and total != int(match[3]) - offset)
                    or (saved.get("total") is not None and int(match[3]) != saved["total"])):
                metadata.unlink(missing_ok=True)
                raise ValueError("Invalid Content-Range; refusing to append mismatched data")
            total = int(match[3])
            returned = _header(source, saved["validator_header"])
            if returned and returned != validator:
                metadata.unlink(missing_ok=True)
                raise ValueError("Download changed; refusing to append mismatched data")
        elif status == 200:
            offset = 0
        else:
            raise ValueError(f"Unexpected download HTTP status: {status}")
        etag = _header(source, "ETag")
        modified = _header(source, "Last-Modified")
        if etag and not etag.startswith("W/"):
            validator, kind = etag, "ETag"
        elif modified:
            validator, kind = modified, "Last-Modified"
        elif offset:
            kind = saved["validator_header"]
        else:
            validator, kind = None, None
        # Truncate first, then publish the validator: an old prefix must never
        # acquire metadata for a new remote representation after a crash.
        with partial.open("ab" if offset else "wb") as target:
            state = dict(url=url, validator=validator, validator_header=kind, total=total)
            temporary = Path(str(metadata) + ".partial")
            temporary.write_text(json.dumps(state), encoding="utf-8")
            temporary.replace(metadata)
            done, started, notified = offset, time.monotonic(), 0.0
            name = partial.name.removesuffix(".partial")
            notify(name=name, done_bytes=done, total_bytes=total, rate_bps=0, resumed_from=offset)
            while True:
                block = source.read(256 * 1024)
                if not block:
                    break
                target.write(block)
                target.flush()
                done += len(block)
                current = time.monotonic()
                if current - notified >= 0.5:
                    notify(name=name, done_bytes=done, total_bytes=total,
                           rate_bps=(done - offset) / max(current - started, 0.001), resumed_from=offset)
                    notified = current
            if not done or (total is not None and done != total):
                raise ValueError(f"Incomplete download: Content-Length mismatch ({done}/{total})")
            notify(name=name, done_bytes=done, total_bytes=total,
                   rate_bps=(done - offset) / max(time.monotonic() - started, 0.001), resumed_from=offset)


def discard_unvalidated_partial(partial):
    """Retain bytes only when their representation has an HTTP validator."""
    metadata = Path(str(partial) + ".json")
    try:
        saved = json.loads(metadata.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        saved = {}
    if not isinstance(saved, dict) or not saved.get("validator"):
        Path(partial).unlink(missing_ok=True)
        metadata.unlink(missing_ok=True)
