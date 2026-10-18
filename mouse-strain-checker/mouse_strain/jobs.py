"""Process lifecycle helpers for the local web UI (no Streamlit dependency)."""
import os
from pathlib import Path
import signal
import subprocess


def launch(command, log_path, *, env=None, append=False):
    options = {"creationflags": subprocess.CREATE_NEW_PROCESS_GROUP} if os.name == "nt" else {"start_new_session": True}
    environment = dict(os.environ, PYTHONUNBUFFERED="1", PYTHONIOENCODING="utf-8")
    if env:
        environment.update(env)
    with open(log_path, "ab" if append else "wb") as handle:
        return subprocess.Popen(command, stdout=handle, stderr=subprocess.STDOUT,
                                env=environment, **options)


def stop(process, timeout=10):
    """Stop the analysis and its mapping/download subprocesses within a bounded time."""
    if process.poll() is not None:
        return
    if os.name == "nt":
        result = subprocess.run(["taskkill", "/PID", str(process.pid), "/T", "/F"],
                                capture_output=True, timeout=timeout)
        if result.returncode and process.poll() is None:
            raise RuntimeError("Could not stop the analysis process tree")
    else:
        try:
            os.killpg(process.pid, signal.SIGTERM)
        except ProcessLookupError:
            pass
        try:
            process.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            pass
        finally:
            # The CLI can exit before a child that ignores SIGTERM. Finish the
            # entire group even when waiting for the group leader succeeded.
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
    process.wait(timeout=timeout)


def completion_signature(directory):
    """Distinguish a new completion marker from results left by an earlier job."""
    try:
        status = (Path(directory) / "complete.json").stat()
    except FileNotFoundError:
        return None
    return status.st_mtime_ns, status.st_ctime_ns, status.st_size, status.st_ino
