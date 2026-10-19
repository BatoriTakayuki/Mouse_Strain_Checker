"""SRX -> FASTQ -> HISAT2 -> expressed-SNP evidence and reports."""
from dataclasses import asdict
from collections import Counter
import csv
import hashlib
import html
import http.client
import io
import json
from pathlib import Path
import re
import shutil
import subprocess
import time
import urllib.parse
import urllib.request

from . import __version__
from .core import Rules, Summary, count_fragments
from .background import background_result
from .panel import load_panel, sha256, canonical
from .progress import emit, transfer
from .transfers import discard_unvalidated_partial, fetch
from .sex import collect_sex_evidence, unavailable_sex
from .preprocess import trim_reads


HISAT2_FLAGS = ("--dta", "-k", "5")


def log(message):
    print(message, flush=True)


def run(command, **kwargs):
    log("$ " + " ".join(map(str, command)))
    return subprocess.run(list(map(str, command)), check=True, **kwargs)


def write_json(path, value):
    path = Path(path)
    tmp = path.with_name(path.name + ".partial")
    tmp.write_text(json.dumps(value, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    tmp.replace(path)


def identity(path):
    p = Path(path).resolve()
    s = p.stat()
    return dict(path=str(p), size=s.st_size, mtime_ns=s.st_mtime_ns)


def tool_version(name):
    result = subprocess.run([name, "--version"], capture_output=True, text=True,
                            encoding="utf-8", errors="replace", check=True)
    lines = [line.strip() for line in (result.stdout + "\n" + result.stderr).splitlines() if line.strip()]
    if not lines:
        raise ValueError(f"{name} did not report its version")
    return lines[0]


def require_tools(names):
    missing = [name for name in names if not shutil.which(name)]
    if missing:
        raise ValueError("Missing executables: " + ", ".join(missing) + "; use Docker or install them first")


def _cache_matches(path, expected):
    """Interrupted or unreadable cache records are misses, never completion markers."""
    try:
        return json.loads(Path(path).read_text(encoding="utf-8")) == expected
    except (OSError, ValueError):
        return False


def _validate_threads(threads):
    if not isinstance(threads, int) or isinstance(threads, bool) or threads < 1:
        raise ValueError("threads must be a positive integer")


def download(url, dest, expected_md5=None):
    if expected_md5 is not None and not re.fullmatch(r"[0-9a-fA-F]{32}", expected_md5):
        raise ValueError("Expected MD5 must contain 32 hexadecimal digits")
    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and expected_md5:
        if file_md5(dest) == expected_md5.lower():
            log(f"Reuse verified download: {dest.name}")
            return
    elif dest.is_file() and dest.stat().st_size:
        if _cache_matches(str(dest) + ".download.json", dict(url=url, file=identity(dest))):
            return
    temp = dest.with_name(dest.name + ".partial")
    for attempt in range(3):
        try:
            log(f"Download: {url}")
            fetch(url, temp, transfer)
            if not temp.stat().st_size:
                raise ValueError("Empty download")
            if expected_md5 and file_md5(temp) != expected_md5.lower():
                temp.unlink(missing_ok=True)
                Path(str(temp) + ".json").unlink(missing_ok=True)
                raise ValueError(f"MD5 mismatch: {dest.name}")
            temp.replace(dest)
            Path(str(temp) + ".json").unlink(missing_ok=True)
            write_json(str(dest) + ".download.json", dict(url=url, file=identity(dest)))
            return
        except (OSError, ValueError, http.client.HTTPException) as error:
            discard_unvalidated_partial(temp)
            if attempt == 2:
                if isinstance(error, http.client.HTTPException):
                    raise RuntimeError(f"Download failed for {dest.name} ({url}): {error}") from error
                raise
            time.sleep(attempt + 1)


def file_md5(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def resolve_srx(accession):
    if not re.fullmatch(r"[SED]RX\d+", accession):
        raise ValueError("Experiment accession must be SRX/ERX/DRX followed by digits")
    fields = ("run_accession,experiment_accession,sample_accession,tax_id,scientific_name,"
              "library_strategy,library_source,library_layout,instrument_platform,fastq_ftp,fastq_md5")
    url = "https://www.ebi.ac.uk/ena/portal/api/filereport?" + urllib.parse.urlencode(
        dict(accession=accession, result="read_run", fields=fields, format="tsv"))
    last = None
    for attempt in range(3):
        try:
            with urllib.request.urlopen(url, timeout=90) as response:
                reader = csv.DictReader(io.StringIO(response.read().decode("utf-8")), delimiter="\t")
                if not reader.fieldnames or not set(fields.split(",")).issubset(reader.fieldnames):
                    raise ValueError("Malformed ENA metadata: required columns are missing")
                rows = list(reader)
            break
        except (OSError, http.client.HTTPException) as error:
            last = error
            if attempt == 2:
                raise RuntimeError(f"ENA lookup failed for {accession}") from last
            time.sleep(attempt + 1)
    if not rows:
        raise ValueError("No public ENA runs found for experiment")
    for row in rows:
        if None in row or any(row.get(field) is None for field in fields.split(",")):
            raise ValueError("Malformed ENA metadata: inconsistent column count")
        if row["experiment_accession"] != accession:
            raise ValueError("ENA returned a different experiment")
        if not re.fullmatch(r"[SED]RR\d+", row["run_accession"]):
            raise ValueError("Invalid ENA run accession")
        if not row["sample_accession"].strip():
            raise ValueError("ENA sample accession is missing; refusing to pool unresolved samples")
        name = row["scientific_name"]
        if not (name == "Mus musculus" or name.startswith("Mus musculus ")):
            raise ValueError("Non-mouse or unresolved species metadata")
        if not row["tax_id"].isdigit() or (name == "Mus musculus" and row["tax_id"] != "10090"):
            raise ValueError("Non-mouse or unresolved species metadata")
        if row["library_strategy"].upper() != "RNA-SEQ":
            raise ValueError("Experiment is not annotated RNA-Seq")
        if not row["library_source"].upper().startswith("TRANSCRIPTOMIC"):
            raise ValueError("Experiment is not a transcriptomic library")
        if row["instrument_platform"].upper() != "ILLUMINA":
            raise ValueError("v0.1 supports Illumina short-read RNA-seq only")
        if row["library_layout"] not in ("SINGLE", "PAIRED"):
            raise ValueError("Unknown library layout")
    if len({r["sample_accession"] for r in rows}) != 1:
        raise ValueError("Experiment maps to multiple biological samples; refusing to pool")
    if len({r["run_accession"] for r in rows}) != len(rows):
        raise ValueError("ENA returned duplicate runs")
    return sorted(rows, key=lambda x: x["run_accession"])


def _sra_payload(accession_directory, accession):
    """Select the completed run itself, never prefetch locks or partial files."""
    # prefetch can append its own accession folder to --output-directory.
    # Accept only exact run names in the requested directory or that one child.
    # Never recursively pick a partial file or an unrelated run/dependency.
    roots = [accession_directory]
    child = accession_directory / accession
    if child.is_dir() and not child.is_symlink():
        roots.append(child)
    candidates = [root / name for root in roots for name in (f"{accession}.sra", f"{accession}.sralite", accession)]
    completed = [p for p in candidates if not p.is_symlink() and p.is_file() and p.stat().st_size
                 and p.resolve().is_relative_to(accession_directory.resolve())]
    # A run can have both distributions cached. Prefer the normalized source;
    # retaining Lite alongside it must not make a valid normalized run ambiguous.
    normalized = [p for p in completed if p.suffix != ".sralite"]
    if normalized:
        completed = normalized
    if len(completed) != 1:
        raise ValueError("Prefetch must produce exactly one nonempty SRA payload for " + accession
                         + f" in {accession_directory} or its accession subdirectory; found {len(completed)}")
    return completed[0]


def _sra_quality(payload):
    """Record the quality representation instead of treating Lite as full BQ."""
    if payload.suffix == ".sralite":
        return dict(format="sra_lite", base_quality="simplified", note=
            "SRA Lite由来です。塩基品質値はリードごとに簡略化され、元の塩基別品質値ではありません。fastpとBQフィルターは簡略化された品質値に適用しています。",
            source="https://www.ncbi.nlm.nih.gov/sra/docs/sra-data-formats/")
    return dict(format="sra_normalized" if payload.suffix == ".sra" else "unspecified",
                base_quality="full" if payload.suffix == ".sra" else "unspecified", note="")


def _sra_fastqs(staging, accession, layout):
    """Check this conversion's isolated output before publishing any cached reads."""
    names = ([f"{accession}_1.fastq", f"{accession}_2.fastq"]
             if layout == "PAIRED" else [f"{accession}.fastq"])
    if layout == "SINGLE" and not (staging / names[0]).exists():
        names = [f"{accession}_1.fastq"]
    paths = [staging / name for name in names]
    actual = {p.name for p in staging.glob("*.fastq")}
    allowed = set(names)
    if layout == "PAIRED":
        # split-3 keeps reads lacking their mate in a separate orphan file.
        allowed.add(f"{accession}.fastq")
    if (not set(names).issubset(actual) or not actual.issubset(allowed)
            or not all(not p.is_symlink() and p.is_file() and p.stat().st_size for p in paths)):
        raise ValueError("SRA FASTQ output does not match declared layout")
    if layout == "PAIRED" and f"{accession}.fastq" in actual:
        log("SRA orphan reads excluded; paired reads only")
    return paths


def get_reads(row, directory, threads, method="sra"):
    _validate_threads(threads)
    if method not in ("auto", "ena", "sra"):
        raise ValueError("download method must be auto, ena, or sra")
    if not re.fullmatch(r"[SED]RR\d+", row["run_accession"]):
        raise ValueError("Invalid ENA run accession")
    if row["library_layout"] not in ("SINGLE", "PAIRED"):
        raise ValueError("Unknown library layout")
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    urls = row["fastq_ftp"].split(";") if row["fastq_ftp"] else []
    checksums = row["fastq_md5"].split(";") if row["fastq_md5"] else []
    valid_ena = len(urls) == len(checksums) and all(re.fullmatch(r"[0-9a-fA-F]{32}", x) for x in checksums)
    selected = list(zip(urls, checksums))
    if row["library_layout"] == "PAIRED":
        mates = []
        for item in selected:
            match = re.fullmatch(r"(.+)_([12])(\.f(?:ast)?q\.gz)", item[0])
            if match:
                mates.append((match.group(1), match.group(3), match.group(2), item))
        selected = [mate[3] for mate in sorted(mates)]
        valid_pair = (len(mates) == 2 and {mate[2] for mate in mates} == {"1", "2"}
                      and len({mate[:2] for mate in mates}) == 1)
        valid_ena = valid_ena and valid_pair
        if valid_pair and len(urls) != 2:
            log("ENA orphan FASTQ excluded; only paired files are used")
    expected = 2 if row["library_layout"] == "PAIRED" else 1
    valid_ena = valid_ena and len(selected) == expected
    if method == "ena" and not valid_ena:
        raise ValueError("ENA lacks a complete FASTQ/checksum set; use --download-method sra")
    if method == "ena" and valid_ena:
        emit("download", f"{row['run_accession']}: FASTQを取得・検証中")
        paths = []
        for remote, checksum in sorted(selected):
            url = remote if "://" in remote else "https://" + remote
            parsed = urllib.parse.urlparse(url)
            if parsed.scheme != "https" or parsed.hostname not in ("ftp.sra.ebi.ac.uk", "ftp.ebi.ac.uk"):
                raise ValueError("Unexpected ENA FASTQ host")
            dest = directory / Path(parsed.path).name
            download(url, dest, checksum)
            paths.append(dest)
        return paths
    accession = row["run_accession"]
    tools = ["prefetch", "vdb-validate", "fasterq-dump"]
    require_tools(tools)
    versions = {name: tool_version(name) for name in tools}
    marker = directory / "sra_reads.json"
    paths = ([directory / f"{accession}_1.fastq", directory / f"{accession}_2.fastq"]
             if expected == 2 else [directory / f"{accession}.fastq"])
    accession_directory = directory / "sra" / accession
    provenance = dict(schema=3, method="sra", accession=accession, layout=row["library_layout"],
                      options=dict(split_mode="split-3", skip_technical=True), tools=versions)
    if all(not p.is_symlink() and p.is_file() and p.stat().st_size for p in paths):
        try:
            cached_payload = _sra_payload(accession_directory, accession)
            cached = dict(provenance, sra=identity(cached_payload), quality=_sra_quality(cached_payload),
                          reads=[identity(p) for p in paths])
        except (OSError, ValueError):
            cached = None
        if cached is not None and _cache_matches(marker, cached):
            log(f"Reuse verified SRA FASTQ: {accession}")
            return paths
    emit("sra_download", f"{accession}: SRAデータを取得中")
    # Repeating prefetch resumes its own incomplete download and fetches dependencies.
    # Keep this location stable so existing prefetched data is reused.
    run(["prefetch", accession, "--max-size", "u", "--output-directory", accession_directory])
    payload = _sra_payload(accession_directory, accession)
    quality = _sra_quality(payload)
    if quality["note"]:
        log(f"{accession}: {quality['note']}")
    emit("sra_validate", f"{accession}: SRAデータの整合性を確認中")
    run(["vdb-validate", payload])
    emit("fastq_convert", f"{accession}: SRAをFASTQへ変換中")
    # A killed or failed conversion cannot turn partial FASTQ files into a cache hit.
    # Keep the prior cache until all replacement outputs pass validation.
    import tempfile
    with tempfile.TemporaryDirectory(prefix=f".{accession}.fasterq-", dir=directory) as temporary:
        staging = Path(temporary)
        output = staging / "reads"
        scratch = staging / "tmp"
        output.mkdir()
        scratch.mkdir()
        # Use the accession directory so reference-compressed runs find dependencies.
        run(["fasterq-dump", payload.parent, "--split-3", "--skip-technical",
             "-e", threads, "-O", output, "-t", scratch])
        extracted = _sra_fastqs(output, accession, row["library_layout"])
        for source, destination in zip(extracted, paths):
            source.replace(destination)
    write_json(marker, dict(provenance, sra=identity(payload), quality=quality, reads=[identity(p) for p in paths]))
    return paths


def build_index(fasta, prefix, assembly, threads=4, snp=None):
    _validate_threads(threads)
    emit("hisat2_index", "HISAT2のゲノム索引を作成中（初回は時間がかかります）")
    prefix = Path(prefix).resolve()
    prefix.parent.mkdir(parents=True, exist_ok=True)
    if any(prefix.parent.glob(prefix.name + ".*.ht2*")):
        raise ValueError("Index already exists; use another prefix")
    provenance = dict(assembly=assembly, fasta_sha256=sha256(fasta),
                      alignment_reference="snp_graph" if snp else "linear",
                      snp_sha256=sha256(snp) if snp else None)
    sidecar = Path(str(prefix) + ".reference.json")
    sidecar.unlink(missing_ok=True)
    command = ["hisat2-build", "-p", threads]
    if snp:
        command += ["--snp", str(Path(snp).resolve())]
    command += [str(Path(fasta).resolve()), str(prefix)]
    try:
        run(command)
        if not any(all(Path(f"{prefix}.{i}.{ext}").is_file()
                       and Path(f"{prefix}.{i}.{ext}").stat().st_size for i in range(1, 9))
                   for ext in ("ht2", "ht2l")):
            raise ValueError("HISAT2 build did not produce 8 nonempty index parts")
        write_json(sidecar, provenance)
    except BaseException:
        # No index parts existed before this invocation, so these belong to the failed build.
        for ext in ("ht2", "ht2l"):
            for i in range(1, 9):
                Path(f"{prefix}.{i}.{ext}").unlink(missing_ok=True)
        sidecar.unlink(missing_ok=True)
        raise


def index_identity(prefix, meta):
    prefix = Path(prefix).resolve()
    sidecar = Path(str(prefix) + ".reference.json")
    if not sidecar.exists():
        raise ValueError("Index needs .reference.json; build it with the build-index command")
    try:
        provenance = json.loads(sidecar.read_text(encoding="utf-8"))
    except (OSError, ValueError) as error:
        raise ValueError("Invalid HISAT2 index provenance sidecar; rebuild the index") from error
    if (not isinstance(provenance, dict) or provenance.get("assembly") != meta["assembly"]
            or provenance.get("fasta_sha256") != meta["fasta_sha256"]):
        raise ValueError("HISAT2 index provenance does not match panel/FASTA")
    for ext in ("ht2", "ht2l"):
        parts = [Path(f"{prefix}.{i}.{ext}") for i in range(1, 9)]
        if all(p.is_file() and p.stat().st_size for p in parts):
            return dict(provenance=provenance, parts=[identity(p) for p in parts])
    raise ValueError("Incomplete HISAT2 index (8 parts required)")


def map_run(paths, row, prefix, index_meta, output, threads):
    _validate_threads(threads)
    if len(paths) not in (1, 2) or len({Path(p).resolve() for p in paths}) != len(paths):
        raise ValueError("Mapping requires one FASTQ or two distinct mate FASTQs")
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    marker = Path(str(output) + ".done.json")
    key = dict(reads=[identity(p) for p in paths], run=row["run_accession"],
               index=index_meta, flags=list(HISAT2_FLAGS), version=__version__,
               hisat2=tool_version("hisat2"), samtools=tool_version("samtools"))
    if output.is_file():
        if _cache_matches(marker, dict(inputs=key, bam=identity(output))):
            run(["samtools", "quickcheck", output])
            log(f"Reuse HISAT2 BAM: {output.name}")
            return output
    marker.unlink(missing_ok=True)
    emit("mapping", f"{row['run_accession']}: HISAT2マッピング・BAMを作成中")
    temporary = output.with_name(output.stem + ".partial.bam")
    command = ["hisat2", "-p", str(threads), *HISAT2_FLAGS, "-x", str(prefix),
               "--rg-id", row["run_accession"], "--rg", "SM:" + row["sample_accession"]]
    command += ["-1", str(paths[0]), "-2", str(paths[1])] if len(paths) == 2 else ["-U", str(paths[0])]
    log("$ " + " ".join(command) + " | samtools sort")
    with open(str(output) + ".hisat2.log", "w", encoding="utf-8") as err:
        mapper = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=err)
        try:
            sorter = subprocess.run(["samtools", "sort", "-@", str(max(0, threads - 1)),
                                      "-m", "512M", "-o", str(temporary), "-"], stdin=mapper.stdout)
            mapper.stdout.close()
            mapper_rc = mapper.wait()
            if mapper_rc or sorter.returncode:
                raise RuntimeError(f"Mapping failed: HISAT2={mapper_rc}, sort={sorter.returncode}; see {output}.hisat2.log")
        except BaseException:
            mapper.stdout.close()
            if mapper.poll() is None:
                mapper.terminate()
                mapper.wait()
            temporary.unlink(missing_ok=True)
            raise
    run(["samtools", "quickcheck", temporary])
    temporary.replace(output)
    run(["samtools", "index", "-@", threads, output])
    write_json(marker, dict(inputs=key, bam=identity(output)))
    return output


def check_bam(bam, contigs):
    run(["samtools", "quickcheck", bam])
    header = run(["samtools", "view", "-H", bam], capture_output=True, text=True).stdout
    found, sort_order = {}, None
    for line in header.splitlines():
        values = dict(x.split(":", 1) for x in line.split("\t")[1:] if ":" in x)
        if line.startswith("@SQ\t"):
            found[values["SN"]] = int(values["LN"])
        elif line.startswith("@HD\t"):
            sort_order = values.get("SO")
    if sort_order != "coordinate":
        raise ValueError("BAM must be coordinate sorted")
    if any(found.get(chrom) != length for chrom, length in contigs.items()):
        raise ValueError("BAM chromosome names/lengths do not match panel reference; check assembly")


def write_csv(path, rows):
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def sex_report(bam, meta, rules):
    """Sex inference is optional evidence; failure must not erase strain results."""
    if bam is None:
        return unavailable_sex("性別推定に必要なBAMがないため、未評価です。", meta.get("assembly"))
    if meta.get("assembly") != "GRCm39":
        return unavailable_sex("性別マーカーはGRCm39の座標にのみ対応しているため、この参照では未評価です。", meta.get("assembly"))
    try:
        return collect_sex_evidence(bam, meta["assembly"], meta["contigs"], rules.min_mapq)
    except (OSError, ValueError, RuntimeError, subprocess.SubprocessError) as error:
        log(f"Sex marker evaluation unavailable: {error}")
        return unavailable_sex("性別マーカーの取得に失敗したため判定できません。解析ログを確認してください。", meta.get("assembly"))


def render_html(directory, chromosomes, windows, strains, sex=None, workflow=None):
    colors = ["#2563eb", "#059669", "#c026d3", "#d97706", "#dc2626", "#0891b2", "#7c3aed", "#65a30d"]
    palette = {s: colors[i % len(colors)] for i, s in enumerate(strains)}
    bars = []
    for row in chromosomes:
        cells = []
        for w in windows:
            if w["chromosome"] != row["chromosome"]:
                continue
            call = w["call"]
            background = w.get("background_call", "")
            color = (palette.get(call, "#737373") if call else
                     "repeating-linear-gradient(135deg,#c7dacc 0,#c7dacc 5px,#eef4ef 5px,#eef4ef 10px)" if background else "#d1d5db")
            label = f'{w["start"]}-{w["end"]}: {w["status"]}; {call}; 候補={w.get("candidate_options", "")}; 参考背景={background}; SNPs={w["n_sites"]}; {w.get("reason", "")}'
            cells.append(f'<span style="flex:{w["end"]-w["start"]+1};background:{color}" title="{html.escape(label)}"></span>')
        bars.append(f'<div class="chr"><b>{html.escape(row["chromosome"])}</b><div class="bar">{"".join(cells)}</div></div>')
    headers = list(chromosomes[0])
    table = "<tr><th>結果（候補は未確定）</th>" + "".join(f"<th>{html.escape(h)}</th>" for h in headers) + "</tr>"
    for row in chromosomes:
        finding = row.get("call") or ("候補：" + row["candidate_options"] + "（未確定）" if row.get("candidate_options") else "候補を絞れません")
        table += "<tr><td>" + html.escape(finding) + "</td>" + "".join(f"<td>{html.escape(str(row[h]) if row[h] is not None else '')}</td>" for h in headers) + "</tr>"
    legend = " ".join(f'<span style="color:{c}">● {html.escape(s)}</span>' for s, c in palette.items())
    document = '''<!doctype html><html lang="ja"><meta charset="utf-8"><meta name="viewport" content="width=device-width">
    <title>Mouse Strain Checker</title><style>body{font:15px system-ui;margin:36px;color:#182238;max-width:1400px}h1{font-size:28px}p{line-height:1.8}.chr{display:flex;gap:15px;margin:7px 0}.chr b{width:65px}.bar{display:flex;flex:1;height:20px;gap:1px}table{border-collapse:collapse;font-size:12px}td,th{padding:9px;border-bottom:1px solid #ddd;text-align:left}.scroll{overflow:auto}.notice{background:#fff4d7;padding:18px;border-radius:10px}</style>
    <h1>Mouse Strain Checker</h1><p class="notice">探索用の暫定結果です。RNAで読めた既知SNPにおける候補系統との適合を示します。染色体全体の由来、DNA上の系統、交雑率、F1を確定するものではありません。灰色の区間は判定不能、濃い灰色は2系統モデルです。</p>'''
    sex = sex or dict(label="未評価", reason="性別マーカーを評価していません。")
    document += '<h2>性別（RNAからの推定）</h2><p><b>' + html.escape(str(sex.get("label", "判定保留"))) + "</b><br>" + html.escape(str(sex.get("reason", ""))) + "</p>"
    markers = sex.get("markers", [])
    if markers:
        document += "<p>マーカーごとのユニーク断片数：" + " / ".join(html.escape(str(m.get("gene", ""))) + "=" + ("≥" if m.get("count_capped") else "") + html.escape(str(m.get("fragments", 0))) for m in markers) + "</p>"
    document += "<p>候補の or は、今回のデータで区別しきれない比較対象モデルを表します。+ は2系統モデルで、F1や混在の確定ではありません。情報不足や適合不良の場合は候補を絞りません。</p>"
    document += "<p>" + legend + "</p><p>淡い斜線は亜系統をまとめた参考背景です。厳密な判定やF1の確定ではありません。</p>" + "".join(bars)
    if isinstance(workflow, dict):
        stages = [str(stage) for stage in workflow.get("stages", [])]
        document += "<h2>RNA-seqの処理</h2><p>" + html.escape(" → ".join(stages)) + "</p>"
        for note in workflow.get("quality_notes", []):
            document += '<p class="notice">' + html.escape(str(note)) + "</p>"
    document += '<p>各バーは同じ幅に正規化しています。区間にカーソルを合わせると詳細を表示します。全区間の数値は windows.csv に保存されています。</p><div class="scroll"><table>' + table + "</table></div></html>"
    (Path(directory) / "report.html").write_text(document, encoding="utf-8")


def infer(bam, fasta, panel, output, rules=None, pairs=True, pileup_lines=None, workflow=None):
    """Optional pileup_lines supports reproducible fixture-based validation."""
    rules = rules or Rules()
    has_bam = pileup_lines is None
    emit("inference", "既知SNPと照合し、染色体・区間別の候補を集計中")
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    meta, db = load_panel(panel)
    chromosome = {c: Summary(meta["strains"], pairs and canonical(c) not in ("X", "Y", "MT")) for c in meta["contigs"]}
    windows = {}
    stats = Counter()
    process = None
    error_handle = None
    try:
        if pileup_lines is None:
            check_bam(bam, meta["contigs"])
            command = ["samtools", "mpileup", "-B", "-q", str(rules.min_mapq),
                       "-Q", "0", "-d", str(rules.max_pileup_depth),
                       "--ff", "UNMAP,SECONDARY,QCFAIL,DUP,SUPPLEMENTARY",
                       "-x", "--output-extra", "QNAME,RG,NH",
                       "-f", str(fasta), "-l", str(Path(panel) / "sites.bed"), str(bam)]
            # -Q 0 preserves entry/quality alignment; BQ filtering and mate collapse are Python-side.
            log("$ " + " ".join(command))
            error_handle = open(output / "mpileup.log", "w", encoding="utf-8")
            process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=error_handle,
                                       text=True, encoding="utf-8")
            pileup_lines = process.stdout
        with gzip_counts(output / "sites.tsv.gz") as counts_file:
            counts_file.write("chromosome\tposition\tref\talt\tref_fragments\talt_fragments\tother_fragments\tused\n")
            for line in pileup_lines:
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 4:
                    raise ValueError("Malformed mpileup line")
                chrom, pos = fields[0], int(fields[1])
                snp = db.execute("SELECT ref,alt,gt FROM snps WHERE chrom=? AND pos=?", (chrom, pos)).fetchone()
                if not snp:
                    raise ValueError("Pileup emitted a position outside panel")
                r, alt, gt = snp
                if fields[2].upper() != r:
                    raise ValueError("Pileup reference base disagrees with panel")
                counts = count_fragments(fields, rules, stats)
                nr, na = counts.get(r, 0), counts.get(alt, 0)
                other = sum(counts.values()) - nr - na
                total = nr + na + other
                stats["covered_panel_positions"] += 1
                saturated = int(fields[3]) >= rules.max_pileup_depth
                if saturated:
                    stats["pileup_depth_limit_positions"] += 1
                use = not saturated and nr + na >= rules.min_depth and other / max(1, total) <= rules.max_other_fraction
                counts_file.write(f"{chrom}\t{pos}\t{r}\t{alt}\t{nr}\t{na}\t{other}\t{int(use)}\n")
                if not use:
                    stats["filtered_panel_positions"] += 1
                    continue
                stats["used_panel_positions"] += 1
                chromosome[chrom].add(pos, gt, nr, na, rules)
                key = (chrom, (pos - 1) // rules.window_size)
                if key not in windows:
                    windows[key] = Summary(meta["strains"], pairs and canonical(chrom) not in ("X", "Y", "MT"))
                windows[key].add(pos, gt, nr, na, rules)
        if process and process.wait():
            raise RuntimeError("samtools mpileup failed; see mpileup.log")
    finally:
        db.close()
        if process:
            if process.poll() is None:
                process.terminate()
                process.wait()
            process.stdout.close()
        if error_handle:
            error_handle.close()
    emit("sex_inference", "XistとY染色体マーカーから性別を推定中")
    sex = sex_report(bam if has_bam else None, meta, rules)
    return write_reports(output, meta, chromosome, windows, rules, pairs, stats, sex, workflow)


def write_reports(output, meta, chromosome, windows, rules, pairs, stats, sex=None, workflow=None):
    """Write exact calls and a separate, explicitly coarser background summary."""
    output = Path(output)
    if sex is None:
        sex = unavailable_sex("性別マーカーを評価していないため、未評価です。", meta.get("assembly"))
    chrom_rows, window_rows, score_rows = [], [], []
    for chrom, length in meta["contigs"].items():
        local = []
        for index in range((length - 1) // rules.window_size + 1):
            summary = windows.get((chrom, index)) or Summary(meta["strains"], pairs and canonical(chrom) not in ("X", "Y", "MT"))
            row = dict(chromosome=chrom, start=index * rules.window_size + 1,
                       end=min(length, (index + 1) * rules.window_size), **summary.result(rules, window=True),
                       **background_result(summary, rules, window=True))
            window_rows.append(row)
            local.append(row)
        result = chromosome[chrom].result(rules)
        result.update(background_result(chromosome[chrom], rules))
        called = [r for r in local if r["call"]]
        unique_calls = sorted({r["call"] for r in called})
        # A single chromosome label must never override conflicting local ancestry evidence.
        if len(unique_calls) > 1:
            result.update(status="regional_discordance", call="", reason="区間ごとの判定が一致しないため、染色体を単一の判定で代表できません。")
        elif result["call"] and canonical(chrom) not in ("X", "Y", "MT") and len(called) < rules.min_chrom_windows:
            result.update(status="localized_evidence", call="", reason=f"判定できた区間が{len(called)}個で、染色体判定に必要な{rules.min_chrom_windows}個未満です。")
        if result["status"] in ("regional_discordance", "localized_evidence"):
            result["candidate_options_reason"] = result["reason"] + "候補欄は染色体内の観測SNPを合算した比較です。区間別結果も確認してください。"
        background_called = [r for r in local if r["background_call"]]
        background_calls = sorted({r["background_call"] for r in background_called})
        if len(background_calls) > 1:
            result.update(background_status="regional_discordance", background_call="",
                          background_reason="区間ごとの参考背景が一致しないため、染色体全体の参考背景を保留します。")
        elif result["background_call"] and canonical(chrom) not in ("X", "Y", "MT") and len(background_called) < rules.min_chrom_windows:
            result.update(background_status="localized_evidence", background_call="",
                          background_reason=f"参考背景を示せた区間が{len(background_called)}個で、必要な{rules.min_chrom_windows}個未満です。")
        result.update(background_called_windows=len(background_called),
                      background_window_calls="; ".join(background_calls))
        result.update(observed_windows=sum(r["n_sites"] > 0 for r in local),
                      called_windows=len(called), total_windows=len(local),
                      window_calls="; ".join(unique_calls),
                      note="sex_chromosome_or_mitochondrial_single_candidate_only" if canonical(chrom) in ("X", "Y", "MT") else "")
        chrom_rows.append(dict(chromosome=chrom, **result))
        score_rows.extend(dict(chromosome=chrom, **s) for s in chromosome[chrom].rankings())
    write_csv(output / "chromosomes.csv", chrom_rows)
    write_csv(output / "windows.csv", window_rows)
    if score_rows:
        write_csv(output / "candidate_scores.csv", score_rows)
    else:
        (output / "candidate_scores.csv").write_text("chromosome,rank,model,n_sites,mean_log_score,fit_loss_per_site\n", encoding="utf-8")
    render_html(output, chrom_rows, window_rows, meta["strains"], sex, workflow)
    write_json(output / "sex.json", sex)
    write_json(output / "evidence.json", dict(rules=asdict(rules), pairs=pairs,
        qc=dict(stats), panel=meta, exploratory=True, version=__version__, workflow=workflow))
    pair_diagnostics = []
    names = meta["strains"]
    if "C57BL_6J" in names and "C57BL_6NJ" in names:
        a, b = names.index("C57BL_6J"), names.index("C57BL_6NJ")
        count = sum(n for summary in chromosome.values() for gt, n in summary.patterns.items() if gt[a] != gt[b])
        if count < rules.min_diagnostic:
            pair_diagnostics.append(dict(strains=["C57BL_6J", "C57BL_6NJ"], diagnostic_sites=count,
                                         minimum_per_region=rules.min_diagnostic))
    write_json(output / "diagnosis.json", dict(schema=1, summary=dict(
        used_sites=sum(summary.n for summary in chromosome.values()),
        chromosome_status_counts=dict(Counter(row["status"] for row in chrom_rows)),
        background_status_counts=dict(Counter(row["background_status"] for row in chrom_rows)),
        background_called_chromosomes=sum(bool(row["background_call"]) for row in chrom_rows),
        background_called_windows=sum(bool(row["background_call"]) for row in window_rows),
        candidate_options_chromosomes=sum(bool(row.get("candidate_options")) and not row["call"] for row in chrom_rows),
        unresolved_substrain_pairs=pair_diagnostics),
        reference_assumption=meta.get("reference_assumption"), added_reference=meta.get("added_reference"), notes=[
            "厳密判定と既定の閾値は維持しています。参考背景はC57BL_6JとC57BL_6NJだけを同じ表示群として比較します。",
            "参考背景は群内で最上位のモデルを、群外の全モデルに対して既定基準で比較した結果です。群内の亜系統は確定しません。",
            "2系統の参考モデルはF1や混在、汚染を区別するものではありません。",
            "染色体内の区間で参考背景が異なる場合は、染色体全体の参考背景も保留します。",
            "候補欄のorは、比較したモデルのうち既定の識別SNP数またはスコア差で最上位から区別しきれない候補です。情報不足・適合不良では絞り込みません。",
        ]))
    if stats.get("missing_nh_entries", 0):
        log(f"Warning: {stats['missing_nh_entries']} SNP-overlapping entries lacked NH tags and were excluded")
    return chrom_rows


def gzip_counts(path):
    import gzip
    return gzip.open(path, "wt", encoding="utf-8")


def analyze(accession, fasta, panel, index, output, threads=4, rules=None,
            pairs=True, bam=None, bam_assembly=None, download_method="sra", preprocess=True):
    _validate_threads(threads)
    if download_method not in ("auto", "ena", "sra"):
        raise ValueError("download method must be auto, ena, or sra")
    if type(preprocess) is not bool:
        raise ValueError("preprocess must be a boolean")
    if not bam and not re.fullmatch(r"[SED]RX\d+", accession):
        raise ValueError("Experiment accession must be SRX/ERX/DRX followed by digits")
    emit("metadata", f"{accession}: 入力と参照データを確認中", current_accession=accession,
         runs_done=0, runs_total=0)
    output = Path(output).resolve()
    output.mkdir(parents=True, exist_ok=True)
    # Invalidate prior reports immediately so a failed rerun cannot present old results as new.
    for name in ("complete.json", "report.html", "chromosomes.csv", "windows.csv",
                 "candidate_scores.csv", "evidence.json", "diagnosis.json", "sex.json", "workflow.json", "cleanup.json", "sites.tsv.gz", "ena_runs.json", "mpileup.log"):
        (output / name).unlink(missing_ok=True)
    meta, db = load_panel(panel)
    db.close()
    if sha256(fasta) != meta["fasta_sha256"]:
        raise ValueError("FASTA checksum does not match panel")
    rows, idx = [], None
    method = "ena" if download_method == "ena" else "sra"
    workflow = dict(schema=1, name="external_bam_v1" if bam else f"{method}_{'fastp' if preprocess else 'raw'}_v1",
                    source="external_bam" if bam else "public_srx", assembly=meta["assembly"],
                    download_method=None if bam else method,
                    requested_download_method=None if bam else download_method,
                    fastp_enabled=bool(preprocess and not bam),
                    mapping=dict(tool="external" if bam else "hisat2", flags=[] if bam else list(HISAT2_FLAGS)),
                    stages=["既存BAM", "SNP・性別推定"] if bam else
                        (["ENA HTTPS"] if method == "ena" else ["prefetch", "vdb-validate", "fasterq-dump"])
                        + (["fastp"] if preprocess else []) + ["HISAT2", "samtools sort/index", "SNP・性別推定"],
                    runs=[], quality_notes=[])
    if bam:
        if bam_assembly != meta["assembly"]:
            raise ValueError("Existing BAM requires --bam-assembly matching the panel")
        bam = Path(bam).resolve()
    else:
        idx = index_identity(index, meta)
        require_tools(["hisat2", "samtools"] + (["prefetch", "vdb-validate", "fasterq-dump"] if method == "sra" else [])
                      + (["fastp"] if preprocess else []))
        rows = resolve_srx(accession)
        emit("metadata", f"{accession}: {len(rows)} runを確認しました", runs_total=len(rows))
        write_json(output / "ena_runs.json", rows)
        cache = output / "cache"
        cache.mkdir(exist_ok=True)
        bams = []
        for row in rows:
            paths = get_reads(row, cache / row["run_accession"], threads, download_method)
            read_inputs = [identity(path) for path in paths]
            acquisition = (json.loads((cache / row["run_accession"] / "sra_reads.json").read_text(encoding="utf-8"))
                           if method == "sra" else dict(method="ena", fastq_ftp=row["fastq_ftp"], fastq_md5=row["fastq_md5"]))
            note = acquisition.get("quality", {}).get("note")
            if note:
                workflow["quality_notes"].append(row["run_accession"] + ": " + note)
            preprocessing = None
            if preprocess:
                processed = cache / "preprocess" / row["run_accession"]
                paths = trim_reads(paths, row, processed, threads)
                preprocessing = json.loads((processed / "complete.json").read_text(encoding="utf-8"))
            workflow["runs"].append(dict(run_accession=row["run_accession"], library_layout=row["library_layout"],
                                         raw_reads=read_inputs, mapping_reads=[identity(path) for path in paths],
                                         acquisition=acquisition, preprocessing=preprocessing))
            bams.append(map_run(paths, row, index, idx, cache / (row["run_accession"] + ".bam"), threads))
            emit("mapping", f"{accession}: {len(bams)}/{len(rows)} runのBAMが完成しました", runs_done=len(bams))
        if len(bams) == 1:
            bam = bams[0]
        else:
            bam = cache / "merged.bam"
            emit("merge", f"{accession}: 複数runのBAMを統合中")
            run(["samtools", "merge", "-f", "-@", threads, bam, *bams])
            run(["samtools", "index", "-@", threads, bam])
    versions = {}
    for name in (("samtools",) if workflow["source"] == "external_bam" else
                 ("hisat2", "samtools", "prefetch", "vdb-validate", "fasterq-dump", "fastp")):
        if shutil.which(name):
            versions[name] = tool_version(name)
    write_json(output / "workflow.json", workflow)
    result = infer(bam, fasta, panel, output, rules, pairs, workflow=workflow)
    write_json(output / "complete.json", dict(accession=accession, bam=identity(bam),
        assembly=meta["assembly"], index=idx, external_bam_assembly_attestation=bam_assembly,
        source_runs=rows, tools=versions, version=__version__, report_schema=3, workflow=workflow))
    log(f"Complete: {output / 'report.html'}")
    emit("sample_completed", f"{accession}: 系統候補の出力が完了しました")
    return result
