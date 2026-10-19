"""Refresh interpretation from saved allele counts, without mapping or downloads."""
from collections import Counter
from datetime import datetime, timezone
import csv
import gzip
import json
import math
from pathlib import Path
import shutil
import uuid

from .core import Rules, Summary
from .panel import canonical, load_panel
from .pipeline import identity, sex_report, write_json, write_reports
from .sex import unavailable_sex
from .storage import retained_sex


REPORT_FILES = ("chromosomes.csv", "windows.csv", "candidate_scores.csv",
                "report.html", "evidence.json", "diagnosis.json", "sex.json")


def _compare_exact_calls(old, new):
    """A display refresh must reproduce every previously saved exact call."""
    with old.open(encoding="utf-8-sig", newline="") as handle:
        before = list(csv.DictReader(handle))
    with new.open(encoding="utf-8-sig", newline="") as handle:
        after = list(csv.DictReader(handle))
    if len(before) != len(after):
        raise ValueError("Saved counts do not reproduce the original report row count")
    exact = ("chromosome", "start", "end", "status", "call", "candidate", "runner_up", "n_sites",
             "diagnostic_sites", "first_snp", "last_snp", "called_windows", "observed_windows")
    for previous, refreshed in zip(before, after):
        if any(previous.get(key) != refreshed.get(key) for key in exact if key in previous):
            raise ValueError("Saved counts do not reproduce the original exact calls; no report replaced")
        for key in ("margin_per_site", "fit_loss_per_site"):
            if previous.get(key):
                if not math.isclose(float(previous[key]), float(refreshed[key]), rel_tol=1e-10, abs_tol=1e-10):
                    raise ValueError("Saved counts do not reproduce original model scores; no report replaced")


def refresh_reports(panel, output):
    output = Path(output).resolve()
    evidence = json.loads((output / "evidence.json").read_text(encoding="utf-8"))
    completion = json.loads((output / "complete.json").read_text(encoding="utf-8"))
    if not isinstance(evidence, dict) or not isinstance(completion, dict):
        raise ValueError("Completed evidence and completion records are required")
    if not isinstance(evidence.get("rules"), dict) or type(evidence.get("pairs")) is not bool:
        raise ValueError("Saved rules and candidate-model settings are required")
    try:
        rules = Rules(**evidence["rules"])
    except TypeError as error:
        raise ValueError("Unsupported saved rule settings") from error
    pairs = evidence["pairs"]
    meta, database = load_panel(panel)
    stats = Counter(evidence.get("qc", {}))
    chromosome = {}
    windows = {}
    try:
        if evidence.get("panel") != meta:
            raise ValueError("Panel metadata differs from the original analysis")
        chromosome = {name: Summary(meta["strains"], pairs and canonical(name) not in ("X", "Y", "MT"))
                      for name in meta["contigs"]}
        positions, covered, used = {}, 0, 0
        with gzip.open(output / "sites.tsv.gz", "rt", encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                covered += 1
                if row["used"] == "0":
                    continue
                if row["used"] != "1":
                    raise ValueError("Invalid saved site-use flag")
                chrom, pos = row["chromosome"], int(row["position"])
                if chrom not in chromosome or pos <= positions.get(chrom, 0):
                    raise ValueError("Saved SNP coordinates are invalid or duplicated")
                positions[chrom] = pos
                site = database.execute("SELECT ref,alt,gt FROM snps WHERE chrom=? AND pos=?", (chrom, pos)).fetchone()
                if not site or site[:2] != (row["ref"], row["alt"]):
                    raise ValueError("Saved SNP does not match the original panel")
                nr, na = int(row["ref_fragments"]), int(row["alt_fragments"])
                chromosome[chrom].add(pos, site[2], nr, na, rules)
                key = chrom, (pos - 1) // rules.window_size
                if key not in windows:
                    windows[key] = Summary(meta["strains"], pairs and canonical(chrom) not in ("X", "Y", "MT"))
                windows[key].add(pos, site[2], nr, na, rules)
                used += 1
        if used != stats.get("used_panel_positions", 0) or covered != stats.get("covered_panel_positions", 0):
            raise ValueError("Saved site counts do not match the original QC totals")
    finally:
        database.close()

    staging = output / (".report-refresh-" + uuid.uuid4().hex)
    staging.mkdir()
    try:
        # Use only the exact original BAM. Replaced/missing files cannot supply sex evidence.
        sex = unavailable_sex("元のBAMが保存されていないため、性別は未評価です。", meta.get("assembly"))
        preserved = retained_sex(output)
        if preserved is not None:
            sex = preserved
        recorded_bam = completion.get("bam")
        if isinstance(recorded_bam, dict) and isinstance(recorded_bam.get("path"), str):
            bam = Path(recorded_bam["path"])
            if bam.is_file():
                if recorded_bam == identity(bam):
                    sex = sex_report(bam, meta, rules)
                else:
                    sex = unavailable_sex("保存BAMが元の解析時と異なるため、性別は未評価です。", meta.get("assembly"))
        workflow = completion.get("workflow")
        rows = write_reports(staging, meta, chromosome, windows, rules, pairs, stats, sex,
                             workflow if isinstance(workflow, dict) else None)
        for name in ("chromosomes.csv", "windows.csv"):
            _compare_exact_calls(output / name, staging / name)
        stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
        backup = output / "reports_backup" / stamp
        backup.mkdir(parents=True)
        for name in (*REPORT_FILES, "complete.json"):
            if (output / name).is_file():
                shutil.copy2(output / name, backup / name)
        for name in REPORT_FILES:
            (staging / name).replace(output / name)
        completion.update(report_schema=3, reports_updated_at=datetime.now(timezone.utc).isoformat(),
                          previous_reports=str(backup.relative_to(output)))
        write_json(output / "complete.json", completion)
        return dict(backup=str(backup), used_sites=used, sex=sex,
                    background_status_counts=dict(Counter(row["background_status"] for row in rows)),
                    background_calls={row["chromosome"]: row["background_call"] for row in rows if row["background_call"]})
    finally:
        for name in REPORT_FILES:
            (staging / name).unlink(missing_ok=True)
            (staging / (name + ".partial")).unlink(missing_ok=True)
        staging.rmdir()
