"""Per-run fastp preprocessing with validated, restartable output publication."""
import gzip
import json
from pathlib import Path
import re
import subprocess
import tempfile

from .progress import emit


def _identity(path):
    path = Path(path).resolve()
    stat = path.stat()
    return dict(path=str(path), size=stat.st_size, mtime_ns=stat.st_mtime_ns)


def _version():
    result = subprocess.run(["fastp", "--version"], check=True, capture_output=True,
                            text=True, encoding="utf-8", errors="replace")
    lines = (result.stdout + "\n" + result.stderr).strip().splitlines()
    if not lines:
        raise ValueError("fastp did not report its version")
    return lines[0].strip()


def _first_read(path):
    """Check a gzip stream and its first FASTQ record without scanning all reads."""
    if not path.is_file() or not path.stat().st_size:
        raise ValueError(f"Missing or empty trimmed FASTQ: {path.name}")
    try:
        with gzip.open(path, "rt", encoding="ascii") as handle:
            header, sequence, plus, quality = [handle.readline().rstrip("\r\n") for _ in range(4)]
    except (OSError, EOFError, UnicodeError) as error:
        raise ValueError(f"Invalid trimmed FASTQ: {path.name}") from error
    if (not header.startswith("@") or not sequence or not plus.startswith("+")
            or len(sequence) != len(quality)):
        raise ValueError(f"Invalid trimmed FASTQ record: {path.name}")


def _report_summary(path, paired):
    try:
        report = json.loads(path.read_text(encoding="utf-8"))
        before = report["summary"]["before_filtering"]["total_reads"]
        after = report["summary"]["after_filtering"]["total_reads"]
        bases = report["summary"]["after_filtering"]["total_bases"]
    except (OSError, ValueError, KeyError, TypeError) as error:
        raise ValueError("Missing or invalid fastp JSON summary") from error
    if (any(not isinstance(value, int) or isinstance(value, bool) or value < 0
            for value in (before, after, bases)) or before < after
            or after == 0 or bases == 0 or (paired and after % 2)):
        raise ValueError("fastp produced no usable reads or an inconsistent read summary")
    return dict(reads_before=before, reads_after=after, bases_after=bases,
                fragments_after=after // (2 if paired else 1))


def _validate_outputs(reads, html, report, paired):
    for path in reads:
        _first_read(path)
    if not html.is_file() or not html.stat().st_size:
        raise ValueError("Missing or empty fastp HTML report")
    return _report_summary(report, paired)


def trim_reads(paths, row, output, threads):
    """Return SE/PE gzip FASTQs; ``output`` is this run's preprocessing folder.

    ``complete.json`` records the exact inputs, tool, options and published files.
    Output is reusable only when all identities and the lightweight read/report
    checks still match. Raw input FASTQs are always retained.
    """
    if not isinstance(threads, int) or isinstance(threads, bool) or threads < 1:
        raise ValueError("threads must be a positive integer")
    accession = row.get("run_accession", "")
    layout = row.get("library_layout")
    if not re.fullmatch(r"[SED]RR\d+", accession):
        raise ValueError("Invalid run accession for preprocessing")
    if layout not in ("SINGLE", "PAIRED"):
        raise ValueError("Unknown library layout for preprocessing")
    paired = layout == "PAIRED"
    paths = [Path(path).resolve() for path in paths]
    if len(paths) != (2 if paired else 1) or len(set(paths)) != len(paths):
        raise ValueError("FASTQ inputs do not match the library layout")
    for path in paths:
        if not path.is_file() or not path.stat().st_size:
            raise ValueError(f"Missing or empty input FASTQ: {path}")
    output = Path(output).resolve()
    output.mkdir(parents=True, exist_ok=True)
    names = ([f"{accession}_1.trim.fastq.gz", f"{accession}_2.trim.fastq.gz"] if paired
             else [f"{accession}.trim.fastq.gz"])
    reads = [output / name for name in names]
    html, report, marker = [output / name for name in ("fastp.html", "fastp.json", "complete.json")]
    published = reads + [html, report, marker]
    if any(path.is_symlink() or path.resolve() in paths for path in published):
        raise ValueError("Preprocessing outputs must not overwrite input FASTQs or symbolic links")
    # The supplied Ubuntu fastp 0.23.x explicitly caps worker threads at 16.
    effective_threads = min(threads, 16)
    expected = dict(schema=1, method="fastp", run_accession=accession,
                    library_layout=layout, inputs=[_identity(path) for path in paths],
                    tool=dict(name="fastp", version=_version()),
                    options=dict(threads=effective_threads, detect_adapter_for_pe=paired,
                                 defaults=True, gzip_output=True))
    emit("read_preprocess", f"{accession}: fastpで品質管理・アダプター除去を実行しています。",
         run=accession, activity=f"fastp: {effective_threads} threads",
         requested_threads=threads, effective_threads=effective_threads, transfer=None)
    try:
        record = json.loads(marker.read_text(encoding="utf-8"))
        if (isinstance(record, dict) and record.get("provenance") == expected
                and record.get("outputs") == [_identity(path) for path in reads + [html, report]]
                and record.get("summary") == _validate_outputs(reads, html, report, paired)):
            print(f"Reuse verified fastp reads: {accession}", flush=True)
            return reads
    except (OSError, ValueError, EOFError, TypeError):
        pass
    # Invalidate first: a stopped fastp process can never leave a valid marker.
    marker.unlink(missing_ok=True)
    with tempfile.TemporaryDirectory(prefix=".fastp-", dir=output) as staging:
        staged = Path(staging)
        staged_reads = [staged / name for name in names]
        staged_html, staged_report = staged / html.name, staged / report.name
        command = ["fastp", "-i", str(paths[0]), "-o", str(staged_reads[0])]
        if paired:
            command += ["-I", str(paths[1]), "-O", str(staged_reads[1]), "--detect_adapter_for_pe"]
        command += ["--html", str(staged_html), "--json", str(staged_report),
                    "--thread", str(effective_threads)]
        print("$ " + " ".join(command), flush=True)
        subprocess.run(command, check=True)
        summary = _validate_outputs(staged_reads, staged_html, staged_report, paired)
        if [_identity(path) for path in paths] != expected["inputs"]:
            raise ValueError("Input FASTQs changed while fastp was running")
        for source, destination in zip(staged_reads + [staged_html, staged_report], reads + [html, report]):
            source.replace(destination)
        record = dict(provenance=expected, summary=summary,
                      outputs=[_identity(path) for path in reads + [html, report]])
        staged_marker = staged / marker.name
        staged_marker.write_text(json.dumps(record, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        staged_marker.replace(marker)
    return reads
