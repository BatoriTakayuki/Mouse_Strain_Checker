import argparse
from dataclasses import fields
import gzip
import json
from pathlib import Path
import shutil
import subprocess
import sys

from .core import Rules
from .panel import build_panel, vcf_samples
from .pipeline import analyze, build_index, download, write_json, require_tools


FASTA_URL = "https://ftp.ensembl.org/pub/release-113/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz"
VCF_URL = "https://ftp.ebi.ac.uk/pub/databases/mousegenomes/REL-2112-v8-SNPs_Indels/mgp_REL2021_snps.vcf.gz"
DEFAULT_STRAINS = "129S1_SvImJ,C57BL_6NJ,BALB_cJ,DBA_2J,FVB_NJ,CAST_EiJ,PWK_PhJ,WSB_EiJ"


def resources_paths(directory):
    root = Path(directory).resolve()
    config = json.loads((root / "config.json").read_text())
    return {k: str(root / config[k]) for k in ("fasta", "panel", "index")}


def parser():
    p = argparse.ArgumentParser(description="Exploratory mouse strain inference from RNA-seq using HISAT2")
    sub = p.add_subparsers(dest="command", required=True)
    setup = sub.add_parser("setup", help="Download GRCm39/MGP v8, prepare panel, and build HISAT2 index once")
    setup.add_argument("--out", required=True)
    setup.add_argument("--threads", type=int, default=4)
    setup.add_argument("--strains", default=DEFAULT_STRAINS)
    setup.add_argument("--fasta-url", default=FASTA_URL)
    setup.add_argument("--vcf-url", default=VCF_URL)
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
    web = sub.add_parser("web", help="Launch the browser UI (requires streamlit)")
    web.add_argument("--host", default="127.0.0.1")
    web.add_argument("--port", type=int, default=8501)
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
            if (root / "config.json").exists():
                raise ValueError("Setup already complete; use run --resources or choose a new directory")
            setup_key = dict(fasta_url=args.fasta_url, vcf_url=args.vcf_url, strains=args.strains,
                             assembly="GRCm39", reference="C57BL_6J")
            setup_record = root / "setup_inputs.json"
            if setup_record.exists() and json.loads(setup_record.read_text()) != setup_key:
                raise ValueError("Setup inputs changed; use a new resources directory")
            write_json(setup_record, setup_key)
            download(args.fasta_url, root / "reference.fa.gz")
            download(args.vcf_url, root / "strains.vcf.gz")
            if not (root / "reference.fa").exists():
                temp = root / "reference.fa.partial"
                with gzip.open(root / "reference.fa.gz", "rb") as src, open(temp, "wb") as dst:
                    shutil.copyfileobj(src, dst)
                temp.replace(root / "reference.fa")
            if not (root / "panel" / "metadata.json").exists():
                if (root / "panel").exists():
                    raise ValueError("Incomplete panel found; move panel directory aside before retrying setup")
                build_panel(root / "strains.vcf.gz", root / "reference.fa", root / "panel", "GRCm39",
                            args.strains.split(","), add_reference="C57BL_6J")
            if not Path(str(root / "hisat2" / "genome") + ".reference.json").exists():
                build_index(root / "reference.fa", root / "hisat2" / "genome", "GRCm39", args.threads)
            write_json(root / "config.json", dict(fasta="reference.fa", panel="panel", index="hisat2/genome"))
            print(f"Ready: run --resources {root} --srx SRX... --out results")
        elif args.command == "run":
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
            script = str(Path(__file__).with_name("app.py"))
            return subprocess.call([sys.executable, "-m", "streamlit", "run", script,
                                    "--server.address", args.host, "--server.port", str(args.port)])
    except (OSError, ValueError, RuntimeError, subprocess.CalledProcessError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
