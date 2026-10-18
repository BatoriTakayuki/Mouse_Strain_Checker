import argparse
from dataclasses import fields
import gzip
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

from .core import Rules
from .panel import build_panel, load_panel, sha256, vcf_samples
from .pipeline import analyze, build_index, download, index_identity, write_json, require_tools
from .progress import emit


FASTA_URL = "https://ftp.ensembl.org/pub/release-113/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz"
VCF_URL = "https://ftp.ebi.ac.uk/pub/databases/mousegenomes/REL-2112-v8-SNPs_Indels/mgp_REL2021_snps.vcf.gz"
DEFAULT_STRAINS = "129S1_SvImJ,C57BL_6NJ,BALB_cJ,DBA_2J,FVB_NJ,CAST_EiJ,PWK_PhJ,WSB_EiJ"


def resources_paths(directory):
    root = Path(directory).resolve()
    config = json.loads((root / "config.json").read_text(encoding="utf-8"))
    if not isinstance(config, dict) or any(
        not isinstance(config.get(key), str) or not config[key].strip()
        for key in ("fasta", "panel", "index")
    ):
        raise ValueError("Invalid resources config: fasta, panel, and index paths are required")
    return {k: str(root / config[k]) for k in ("fasta", "panel", "index")}


def _clear_unfinished_panel(root):
    directory = root / "panel"
    if directory.resolve().parent != root or directory.is_symlink():
        raise ValueError("Refusing to replace a panel outside the setup directory")
    allowed = {"panel.sqlite.partial", "panel.sqlite.partial-journal", "panel.sqlite.partial-wal",
               "panel.sqlite.partial-shm", "panel.sqlite", "sites.bed", "metadata.json.partial"}
    children = list(directory.iterdir())
    if any(child.name not in allowed or not child.is_file() or child.is_symlink() for child in children):
        raise ValueError("Incomplete panel contains unrecognized files; move it aside before retrying setup")
    for child in children:
        child.unlink()
    directory.rmdir()


def _clear_unfinished_index(root):
    directory = root / "hisat2"
    if directory.resolve().parent != root or directory.is_symlink():
        raise ValueError("Refusing to replace an index outside the setup directory")
    names = [f"genome.{i}.{ext}" for ext in ("ht2", "ht2l") for i in range(1, 9)]
    names += ["genome.reference.json", "genome.reference.json.partial"]
    for name in names:
        path = directory / name
        if path.is_symlink() or (path.exists() and not path.is_file()):
            raise ValueError("Unexpected file type in unfinished HISAT2 index")
    for name in names:
        (directory / name).unlink(missing_ok=True)


def parser():
    p = argparse.ArgumentParser(description="Exploratory mouse strain inference from RNA-seq using HISAT2")
    sub = p.add_subparsers(dest="command", required=True)
    setup = sub.add_parser("setup", help="Download GRCm39/MGP v8, prepare panel, and build HISAT2 index once")
    setup.add_argument("--out", required=True)
    setup.add_argument("--threads", type=int, default=4)
    setup.add_argument("--strains", default=DEFAULT_STRAINS)
    setup.add_argument("--fasta-url", default=FASTA_URL)
    setup.add_argument("--vcf-url", default=VCF_URL)
    setup.add_argument("--resume", action="store_true", help="Rebuild unfinished setup artifacts from the same recorded inputs")
    sample = sub.add_parser("samples", help="List exact VCF sample names")
    sample.add_argument("--vcf", required=True)
    panel = sub.add_parser("panel", help="Create a panel from your own assembly-matched FASTA and VCF")
    panel.add_argument("--vcf", required=True)
    panel.add_argument("--fasta", required=True)
    panel.add_argument("--out", required=True)
    panel.add_argument("--assembly", choices=("GRCm38", "GRCm39"), required=True)
    panel.add_argument("--strains", help="Comma-separated exact VCF sample names; default all")
    panel.add_argument("--add-reference", help="Explicitly assume this added strain is homozygous REF")
    panel.add_argument("--spacing", type=int, default=1000)
    panel.add_argument("--min-gq", type=int, default=20)
    panel.add_argument("--allow-no-fi", action="store_true")
    idx = sub.add_parser("build-index", help="Build HISAT2 index with a reference provenance sidecar")
    idx.add_argument("--fasta", required=True)
    idx.add_argument("--out", required=True)
    idx.add_argument("--assembly", choices=("GRCm38", "GRCm39"), required=True)
    idx.add_argument("--threads", type=int, default=4)
    idx.add_argument("--snp", help="Optional HISAT2 SNP-format file in matching coordinates (graph index)")
    r = sub.add_parser("run", help="Analyze SRX experiments independently, or reuse a sorted HISAT2 BAM")
    r.add_argument("--srx", nargs="+", help="One or more SRX/ERX/DRX experiment accessions")
    r.add_argument("--resources", help="Directory created by setup")
    r.add_argument("--fasta")
    r.add_argument("--panel")
    r.add_argument("--index")
    r.add_argument("--out", required=True)
    r.add_argument("--threads", type=int, default=4)
    r.add_argument("--bam", help="Existing sorted BAM (HISAT2 NH tags required)")
    r.add_argument("--bam-assembly", choices=("GRCm38", "GRCm39"))
    r.add_argument("--download-method", choices=("auto", "ena", "sra"), default="auto")
    r.add_argument("--single-only", action="store_true", help="Disable two-strain models on autosomes")
    for field in fields(Rules):
        default = getattr(Rules(), field.name)
        r.add_argument("--" + field.name.replace("_", "-"), default=default, type=type(default))
    web = sub.add_parser("web", help="Launch the persistent local browser UI")
    web.add_argument("--host", default="127.0.0.1")
    web.add_argument("--port", type=int, default=8501)
    web.add_argument("--data", default=os.environ.get("STRAIN_DATA", "./data"))
    web.add_argument("--resources", help="Prepared reference directory; default DATA/resources")
    web.add_argument("--jobs", help="Persistent analysis history; default DATA/results/web")
    web.add_argument("--threads", type=int, default=max(1, min(4, os.cpu_count() or 1)))
    web.add_argument("--streamlit", action="store_true", help="Use the earlier optional Streamlit UI")
    return p


def main(argv=None):
    p = parser()
    args = p.parse_args(argv)
    try:
        if args.command == "samples":
            print("\n".join(vcf_samples(args.vcf)))
        elif args.command == "panel":
            meta = build_panel(args.vcf, args.fasta, args.out, args.assembly,
                               args.strains.split(",") if args.strains else None,
                               args.add_reference, args.spacing, args.min_gq, args.allow_no_fi)
            print(json.dumps(meta["sites_per_chromosome"], indent=2))
        elif args.command == "build-index":
            if args.threads < 1:
                raise ValueError("threads must be positive")
            require_tools(["hisat2-build"])
            build_index(args.fasta, args.out, args.assembly, args.threads, args.snp)
        elif args.command == "setup":
            if args.threads < 1:
                raise ValueError("threads must be positive")
            require_tools(["samtools", "hisat2-build"])
            root = Path(args.out).resolve()
            root.mkdir(parents=True, exist_ok=True)
            if (root / "config.json").exists() and not args.resume:
                raise ValueError("Setup already complete; use run --resources or choose a new directory")
            setup_key = dict(fasta_url=args.fasta_url, vcf_url=args.vcf_url, strains=args.strains,
                             assembly="GRCm39", reference="C57BL_6J")
            setup_record = root / "setup_inputs.json"
            recorded = setup_record.exists()
            if recorded and json.loads(setup_record.read_text(encoding="utf-8")) != setup_key:
                raise ValueError("Setup inputs changed; use a new resources directory")
            if (root / "config.json").exists():
                paths = resources_paths(root)
                meta, db = load_panel(paths["panel"])
                db.close()
                if meta["fasta_sha256"] != sha256(paths["fasta"]):
                    raise ValueError("Prepared FASTA does not match the panel")
                index_identity(paths["index"], meta)
                emit("reference_ready", "参照データは準備済みです")
                return 0
            write_json(setup_record, setup_key)
            emit("reference_download", "参照FASTAを取得・検証中")
            download(args.fasta_url, root / "reference.fa.gz")
            emit("reference_download", "マウス系統のSNP VCFを取得・検証中")
            download(args.vcf_url, root / "strains.vcf.gz")
            if not (root / "reference.fa").exists():
                emit("reference_extract", "参照FASTAを展開中")
                temp = root / "reference.fa.partial"
                with gzip.open(root / "reference.fa.gz", "rb") as src, open(temp, "wb") as dst:
                    shutil.copyfileobj(src, dst)
                temp.replace(root / "reference.fa")
            if not (root / "panel" / "metadata.json").exists():
                if (root / "panel").exists():
                    if not args.resume or not recorded:
                        raise ValueError("Incomplete panel found; use setup --resume for the same recorded inputs")
                    _clear_unfinished_panel(root)
                if args.resume and recorded:
                    (root / "reference.fa.fai").unlink(missing_ok=True)
                emit("panel", "FASTA索引と系統識別用SNPパネルを作成中")
                build_panel(root / "strains.vcf.gz", root / "reference.fa", root / "panel", "GRCm39",
                            args.strains.split(","), add_reference="C57BL_6J")
            meta, db = load_panel(root / "panel")
            db.close()
            if meta["assembly"] != "GRCm39" or meta["fasta_sha256"] != sha256(root / "reference.fa"):
                raise ValueError("Existing panel does not match setup FASTA/assembly; use a new directory")
            prefix = root / "hisat2" / "genome"
            if args.resume and recorded:
                try:
                    index_identity(prefix, meta)
                except (OSError, ValueError):
                    _clear_unfinished_index(root)
            if not Path(str(prefix) + ".reference.json").exists():
                build_index(root / "reference.fa", root / "hisat2" / "genome", "GRCm39", args.threads)
            index_identity(root / "hisat2" / "genome", meta)
            write_json(root / "config.json", dict(fasta="reference.fa", panel="panel", index="hisat2/genome"))
            emit("reference_ready", "参照FASTA・SNPパネル・HISAT2索引の準備が完了しました")
            print(f"Ready: run --resources {root} --srx SRX... --out results")
        elif args.command == "run":
            if args.threads < 1:
                raise ValueError("threads must be positive")
            paths = resources_paths(args.resources) if args.resources else {}
            paths.update({key: getattr(args, key) for key in ("fasta", "panel", "index") if getattr(args, key)})
            if "fasta" not in paths or "panel" not in paths or (not args.bam and "index" not in paths):
                raise ValueError("Use --resources or provide --fasta --panel --index")
            if args.bam and args.srx:
                raise ValueError("Use either --bam or --srx")
            if not args.bam and not args.srx:
                raise ValueError("--srx is required when --bam is absent")
            rules = Rules(**{f.name: getattr(args, f.name) for f in fields(Rules)})
            accessions = args.srx or ["local_bam"]
            # Validate before constructing output paths.
            import re
            if args.srx and any(not re.fullmatch(r"[SED]RX\d+", x) for x in accessions):
                raise ValueError("Invalid experiment accession")
            require_tools(["samtools"] if args.bam else ["hisat2", "samtools"])
            for acc in dict.fromkeys(accessions):
                analyze(acc, paths["fasta"], paths["panel"], paths.get("index"),
                        Path(args.out) / acc, args.threads, rules, not args.single_only,
                        args.bam, args.bam_assembly, args.download_method)
        else:
            if not 1 <= args.threads <= 128:
                raise ValueError("web threads must be between 1 and 128")
            if args.streamlit:
                script = str(Path(__file__).with_name("app.py"))
                return subprocess.call([sys.executable, "-m", "streamlit", "run", script,
                                        "--server.address", args.host, "--server.port", str(args.port)])
            from .webapp import serve
            root = Path(args.data).resolve()
            serve(args.host, args.port, args.jobs or root / "results" / "web",
                  args.resources or root / "resources", args.threads)
    except (OSError, ValueError, RuntimeError, subprocess.CalledProcessError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
