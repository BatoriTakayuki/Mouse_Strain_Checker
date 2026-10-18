"""Optional structured progress for CLI workers managed by the web server."""
from datetime import datetime, timezone
import json
import os
from pathlib import Path

_state = {}


def emit(stage, message, **details):
    destination = os.environ.get("MOUSE_STRAIN_PROGRESS_FILE")
    if not destination:
        return
    if stage != _state.get("stage"):
        _state.update(transfer=None, activity="")
    _state.update(stage=stage, message=message, **details,
                  updated_at=datetime.now(timezone.utc).isoformat())
    path = Path(destination)
    temporary = path.with_name(path.name + ".partial")
    temporary.write_text(json.dumps(_state, ensure_ascii=False) + "\n", encoding="utf-8")
    temporary.replace(path)


def transfer(**details):
    emit(_state.get("stage", "download"), _state.get("message", "データを取得中"), transfer=details)
