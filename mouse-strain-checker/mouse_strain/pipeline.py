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
from .panel import load_panel, sha256, canonical
from .progress import emit, transfer
from .transfers import discard_unvalidated_partial, fetch


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
    return result.stdout.splitlines()[0] if result.stdout else result.stderr[:200]


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


def get_reads(row, directory, threads, method="auto"):
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
    if method != "sra" and valid_ena:
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
    require_tools(["prefetch", "fasterq-dump"])
    marker = directory / "sra_reads.json"
    paths = ([directory / f"{accession}_1.fastq", directory / f"{accession}_2.fastq"]
             if expected == 2 else [directory / f"{accession}.fastq"])
    if all(p.is_file() and p.stat().st_size for p in paths):
        if _cache_matches(marker, [identity(p) for p in paths]):
            return paths
    for p in directory.glob(f"{accession}*.fastq"):
        p.unlink()
    sra = directory / "sra"
    emit("sra_download", f"{accession}: SRAデータを取得中")
    run(["prefetch", accession, "--max-size", "u", "-O", sra / accession])
    emit("fastq_convert", f"{accession}: SRAをFASTQへ変換中")
    run(["fasterq-dump", sra / accession / f"{accession}.sra", "--split-3",
         "--skip-technical", "-e", threads, "-O", directory, "-t", directory / "tmp"])
    if expected == 1 and not paths[0].exists() and (directory / f"{accession}_1.fastq").exists():
        (directory / f"{accession}_1.fastq").rename(paths[0])
    if not all(p.exists() and p.stat().st_size for p in paths):
        raise ValueError("SRA FASTQ output does not match declared layout")
    if expected == 2 and (directory / f"{accession}.fastq").exists():
        log("SRA orphan reads excluded; paired reads only")
    write_json(marker, [identity(p) for p in paths])
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
               index=index_meta, flags=["--dta", "-k", "5"], version=__version__,
               hisat2=tool_version("hisat2"), samtools=tool_version("samtools"))
    if output.is_file():
        if _cache_matches(marker, dict(inputs=key, bam=identity(output))):
            run(["samtools", "quickcheck", output])
            log(f"Reuse HISAT2 BAM: {output.name}")
            return output
    marker.unlink(missing_ok=True)
    emit("mapping", f"{row['run_accession']}: HISAT2マッピング・BAMを作成中")
    temporary = output.with_name(output.stem + ".partial.bam")
    command = ["hisat2", "-p", str(threads), "--dta", "-k", "5", "-x", str(prefix),
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


def render_html(directory, chromosomes, windows, strains):
    colors = ["#2563eb", "#059669", "#c026d3", "#d97706", "#dc2626", "#0891b2", "#7c3aed", "#65a30d"]
    palette = {s: colors[i % len(colors)] for i, s in enumerate(strains)}
    bars = []
    for row in chromosomes:
        cells = []
        for w in windows:
            if w["chromosome"] != row["chromosome"]:
                continue
            call = w["call"]
            color = palette.get(call, "#737373" if call else "#d1d5db")
            label = f'{w["start"]}-{w["end"]}: {w["status"]}; {call}; SNPs={w["n_sites"]}'
            cells.append(f'<span style="flex:{w["end"]-w["start"]+1};background:{color}" title="{html.escape(label)}"></span>')
        bars.append(f'<div class="chr"><b>{html.escape(row["chromosome"])}</b><div class="bar">{"".join(cells)}</div></div>')
    headers = list(chromosomes[0])
    table = "<tr>" + "".join(f"<th>{html.escape(h)}</th>" for h in headers) + "</tr>"
    for row in chromosomes:
        table += "<tr>" + "".join(f"<td>{html.escape(str(row[h]) if row[h] is not None else '')}</td>" for h in headers) + "</tr>"
    legend = " ".join(f'<span style="color:{c}">● {html.escape(s)}</span>' for s, c in palette.items())
    document = '''<!doctype html><html lang="ja"><meta charset="utf-8"><meta name="viewport" content="width=device-width">
    <title>Mouse Strain Checker</title><style>body{font:15px system-ui;margin:36px;color:#182238;max-width:1400px}h1{font-size:28px}p{line-height:1.8}.chr{display:flex;gap:15px;margin:7px 0}.chr b{width:65px}.bar{display:flex;flex:1;height:20px;gap:1px}table{border-collapse:collapse;font-size:12px}td,th{padding:9px;border-bottom:1px solid #ddd;text-align:left}.scroll{overflow:auto}.notice{background:#fff4d7;padding:18px;border-radius:10px}</style>
    <h1>Mouse Strain Checker</h1><p class="notice">探索用の暫定結果です。RNAで読めた既知SNPにおける候補系統との適合を示します。染色体全体の由来、DNA上の系統、交雑率、F1を確定するものではありません。灰色の区間は判定不能、濃い灰色は2系統モデルです。</p>'''
    document += "<p>" + legend + "</p>" + "".join(bars)
    document += '<p>各バーは同じ幅に正規化しています。区間にカーソルを合わせると詳細を表示します。全区間の数値は windows.csv に保存されています。</p><div class="scroll"><table>' + table + "</table></div></html>"
    (Path(directory) / "report.html").write_text(document, encoding="utf-8")


def infer(bam, fasta, panel, output, rules=None, pairs=True, pileup_lines=None):
    """Optional pileup_lines supports reproducible fixture-based validation."""
    rules = rules or Rules()
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
    chrom_rows, window_rows, score_rows = [], [], []
    for chrom, length in meta["contigs"].items():
        local = []
        for index in range((length - 1) // rules.window_size + 1):
            summary = windows.get((chrom, index)) or Summary(meta["strains"], pairs and canonical(chrom) not in ("X", "Y", "MT"))
            row = dict(chromosome=chrom, start=index * rules.window_size + 1,
                       end=min(length, (index + 1) * rules.window_size), **summary.result(rules, window=True))
            window_rows.append(row)
            local.append(row)
        result = chromosome[chrom].result(rules)
        called = [r for r in local if r["call"]]
        unique_calls = sorted({r["call"] for r in called})
        # A single chromosome label must never override conflicting local ancestry evidence.
        if len(unique_calls) > 1:
            result.update(status="regional_discordance", call="")
        elif result["call"] and canonical(chrom) not in ("X", "Y", "MT") and len(called) < rules.min_chrom_windows:
            result.update(status="localized_evidence", call="")
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
    render_html(output, chrom_rows, window_rows, meta["strains"])
    write_json(output / "evidence.json", dict(rules=asdict(rules), pairs=pairs,
        qc=dict(stats), panel=meta, exploratory=True, version=__version__))
    if stats["missing_nh_entries"]:
        log(f"Warning: {stats['missing_nh_entries']} SNP-overlapping entries lacked NH tags and were excluded")
    return chrom_rows


def gzip_counts(path):
    import gzip
    return gzip.open(path, "wt", encoding="utf-8")


def analyze(accession, fasta, panel, index, output, threads=4, rules=None,
            pairs=True, bam=None, bam_assembly=None, download_method="auto"):
    _validate_threads(threads)
    if download_method not in ("auto", "ena", "sra"):
        raise ValueError("download method must be auto, ena, or sra")
    if not bam and not re.fullmatch(r"[SED]RX\d+", accession):
        raise ValueError("Experiment accession must be SRX/ERX/DRX followed by digits")
    emit("metadata", f"{accession}: 入力と参照データを確認中", current_accession=accession,
         runs_done=0, runs_total=0)
    output = Path(output).resolve()
    output.mkdir(parents=True, exist_ok=True)
    # Invalidate prior reports immediately so a failed rerun cannot present old results as new.
    for name in ("complete.json", "report.html", "chromosomes.csv", "windows.csv",
                 "candidate_scores.csv", "evidence.json", "sites.tsv.gz", "ena_runs.json", "mpileup.log"):
        (output / name).unlink(missing_ok=True)
    meta, db = load_panel(panel)
    db.close()
    if sha256(fasta) != meta["fasta_sha256"]:
        raise ValueError("FASTA checksum does not match panel")
    rows, idx = [], None
    if bam:
        if bam_assembly != meta["assembly"]:
            raise ValueError("Existing BAM requires --bam-assembly matching the panel")
        bam = Path(bam).resolve()
    else:
        idx = index_identity(index, meta)
        rows = resolve_srx(accession)
        emit("metadata", f"{accession}: {len(rows)} runを確認しました", runs_total=len(rows))
        write_json(output / "ena_runs.json", rows)
        cache = output / "cache"
        cache.mkdir(exist_ok=True)
        bams = []
        for row in rows:
            paths = get_reads(row, cache / row["run_accession"], threads, download_method)
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
    for name in ("hisat2", "samtools"):
        if shutil.which(name):
            result = subprocess.run([name, "--version"], capture_output=True, text=True,
                                    encoding="utf-8", errors="replace")
            versions[name] = result.stdout.splitlines()[0] if result.stdout else result.stderr[:200]
    result = infer(bam, fasta, panel, output, rules, pairs)
    write_json(output / "complete.json", dict(accession=accession, bam=identity(bam),
        assembly=meta["assembly"], index=idx, external_bam_assembly_attestation=bam_assembly,
        source_runs=rows, tools=versions, version=__version__))
    log(f"Complete: {output / 'report.html'}")
    emit("sample_completed", f"{accession}: 系統候補の出力が完了しました")
    return result
