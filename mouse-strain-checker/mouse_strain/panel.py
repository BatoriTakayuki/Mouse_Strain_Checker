"""Build a compact panel from a multi-strain, reference-matched VCF."""
from collections import Counter
import gzip
import hashlib
import json
import math
from pathlib import Path
import re
import sqlite3
import subprocess


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def open_text(path):
    return gzip.open(path, "rt", encoding="utf-8") if str(path).endswith(".gz") else open(path, encoding="utf-8")


def canonical(name):
    name = name.removeprefix("chr")
    return "MT" if name in ("M", "MT") else name


def primary(name):
    return canonical(name) in [str(i) for i in range(1, 20)] + ["X", "Y", "MT"]


class Fasta:
    def __init__(self, path):
        self.path = Path(path).resolve()
        if str(path).endswith(".gz"):
            raise ValueError("FASTA must be uncompressed (gunzip it first)")
        if not Path(str(self.path) + ".fai").exists():
            subprocess.run(["samtools", "faidx", str(self.path)], check=True)
        self.index = {}
        with open(str(self.path) + ".fai", encoding="utf-8") as f:
            for line in f:
                name, length, offset, bases, width, *_ = line.rstrip().split("\t")
                entry = tuple(map(int, (length, offset, bases, width)))
                if name in self.index:
                    raise ValueError(f"Duplicate FASTA contig name: {name}")
                if entry[0] < 1 or entry[1] < 0 or entry[2] < 1 or entry[3] < entry[2]:
                    raise ValueError(f"Invalid FASTA index entry: {name}")
                self.index[name] = entry
        self.handle = open(self.path, "rb")

    def base(self, chrom, pos):
        length, offset, bases, width = self.index[chrom]
        if not 1 <= pos <= length:
            raise ValueError(f"VCF coordinate outside FASTA: {chrom}:{pos}")
        zero = pos - 1
        self.handle.seek(offset + zero // bases * width + zero % bases)
        return self.handle.read(1).decode().upper()

    def close(self):
        self.handle.close()


def vcf_samples(path):
    with open_text(path) as f:
        for line in f:
            if line.startswith("#CHROM\t"):
                fields = line.rstrip("\r\n").split("\t")
                fixed = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
                if fields[:8] != fixed or (len(fields) > 8 and (fields[8] != "FORMAT" or len(fields) < 10)):
                    raise ValueError("VCF sample header must include the fixed columns and samples")
                samples = fields[9:]
                if any(not s or any(c.isspace() for c in s) for s in samples) or len(samples) != len(set(samples)):
                    raise ValueError("VCF sample names must be nonempty and unique")
                return samples
    raise ValueError("VCF #CHROM header not found")


def build_panel(vcf, fasta, output, assembly, strains=None, add_reference=None,
                spacing=1000, min_gq=20, allow_no_fi=False):
    if not math.isfinite(spacing) or not math.isfinite(min_gq) or spacing < 1 or min_gq < 0:
        raise ValueError("spacing must be positive and min_gq nonnegative")
    if assembly not in ("GRCm38", "GRCm39"):
        raise ValueError("assembly must be GRCm38 or GRCm39")
    if add_reference is not None and (not isinstance(add_reference, str) or not add_reference.strip()):
        raise ValueError("Added reference name must be nonempty")
    samples = vcf_samples(vcf)
    strains = list(samples if strains is None else strains)
    if len(strains) != len(set(strains)) or set(strains) - set(samples):
        raise ValueError(f"Duplicate or absent strains. VCF samples: {samples}")
    if add_reference in strains:
        raise ValueError("Added reference already exists in selected VCF samples")
    selected = [samples.index(s) for s in strains]
    if add_reference:
        strains.append(add_reference)
    if len(strains) < 2:
        raise ValueError("At least two candidate strains required")
    out = Path(output)
    if out.exists():
        raise ValueError("Panel output already exists; choose a new directory")
    ref = Fasta(fasta)
    db = None
    created_output = False
    try:
        contigs = {n: x[0] for n, x in ref.index.items() if primary(n)}
        aliases = {}
        for n in contigs:
            if canonical(n) in aliases:
                raise ValueError("Ambiguous primary FASTA contig aliases")
            aliases[canonical(n)] = n
        if not contigs:
            raise ValueError("FASTA has no primary mouse chromosomes")
        out.mkdir(parents=True)
        created_output = True
        db = sqlite3.connect(out / "panel.sqlite.partial")
        db.execute("CREATE TABLE snps(chrom TEXT,pos INTEGER,ref TEXT,alt TEXT,gt TEXT,PRIMARY KEY(chrom,pos)) WITHOUT ROWID")
        stats, per_chr, last, seen = Counter(), Counter(), {}, {}
        headers = []
        with open_text(vcf) as source:
            for line in source:
                if line.startswith("##"):
                    if line.startswith(("##reference=", "##contig=", "##VEP=")):
                        headers.append(line.rstrip())
                    # Reject explicit known build conflicts, even if some REF bases happen to match.
                    if line.startswith(("##reference=", "##VEP=")):
                        builds = set(re.findall(r"GRCm(?:38|39)|mm(?:10|39)", line))
                        normalized = {"mm10": "GRCm38", "mm39": "GRCm39"}
                        if any(normalized.get(b, b) != assembly for b in builds):
                            raise ValueError("VCF header assembly disagrees with --assembly")
                    continue
                if line.startswith("#"):
                    continue
                fields = line.rstrip().split("\t")
                stats["records"] += 1
                if len(fields) != 9 + len(samples):
                    raise ValueError("VCF record column count disagrees with sample header")
                chrom = aliases.get(canonical(fields[0]))
                if not chrom:
                    stats["nonprimary"] += 1
                    continue
                pos = int(fields[1])
                if pos <= seen.get(chrom, 0):
                    raise ValueError("VCF must contain unique, position-sorted records per chromosome")
                seen[chrom] = pos
                r, alt = fields[3:5]
                if fields[6] != "PASS" or len(r) != 1 or len(alt) != 1 or r not in "ACGT" or alt not in "ACGT" or r == alt:
                    stats["not_pass_biallelic_snp"] += 1
                    continue
                keys = fields[8].split(":")
                if len(keys) != len(set(keys)):
                    raise ValueError(f"Duplicate VCF FORMAT keys at {chrom}:{pos}")
                if "GT" not in keys:
                    stats["no_gt"] += 1
                    continue
                if "FI" not in keys and not allow_no_fi:
                    raise ValueError("VCF lacks FORMAT/FI. Use --allow-no-fi only for an independently quality-filtered VCF")
                gt = []
                for idx in selected:
                    values = dict(zip(keys, fields[9 + idx].split(":")))
                    alleles = values.get("GT", ".").replace("|", "/").split("/")
                    if not alleles or alleles[0] not in ("0", "1") or any(a != alleles[0] for a in alleles):
                        break
                    if "FI" in keys and values.get("FI") != "1":
                        break
                    if "GQ" in keys:
                        quality = values.get("GQ", ".")
                        if quality == ".":
                            break
                        quality = float(quality)
                        if not math.isfinite(quality) or quality < min_gq:
                            break
                    gt.append(alleles[0])
                if len(gt) != len(selected):
                    stats["missing_het_or_low_quality"] += 1
                    continue
                if add_reference:
                    gt.append("0")
                if len(set(gt)) < 2:
                    stats["uninformative"] += 1
                    continue
                # Check REF even for subsequently thinned sites.
                if ref.base(chrom, pos) != r:
                    raise ValueError(f"REF mismatch at {chrom}:{pos}; check FASTA/VCF assembly")
                if pos - last.get(chrom, -spacing) < spacing:
                    stats["thinned"] += 1
                    continue
                db.execute("INSERT INTO snps VALUES(?,?,?,?,?)", (chrom, pos, r, alt, "".join(gt)))
                last[chrom] = pos
                per_chr[chrom] += 1
        if not sum(per_chr.values()):
            raise ValueError("No informative high-confidence SNPs survived panel filtering")
        db.commit()
        with open(out / "sites.bed", "w", encoding="utf-8") as bed:
            for chrom in contigs:
                for (pos,) in db.execute("SELECT pos FROM snps WHERE chrom=? ORDER BY pos", (chrom,)):
                    bed.write(f"{chrom}\t{pos-1}\t{pos}\n")
        meta = dict(schema=1, assembly=assembly, strains=strains, contigs=contigs,
                    fasta_sha256=sha256(fasta), vcf_sha256=sha256(vcf),
                    vcf_headers=headers, spacing=spacing, min_gq=min_gq,
                    allow_no_fi=allow_no_fi, added_reference=add_reference,
                    reference_assumption="Added strain is assumed REF at selected sites; this is not a measured genotype" if add_reference else None,
                    stats=dict(stats), sites_per_chromosome=dict(per_chr))
        db.close()
        (out / "panel.sqlite.partial").rename(out / "panel.sqlite")
        (out / "metadata.json.partial").write_text(json.dumps(meta, indent=2) + "\n", encoding="utf-8")
        (out / "metadata.json.partial").rename(out / "metadata.json")
        return meta
    except BaseException:
        if db is not None:
            db.close()
        # A failed panel is never accepted by the run command.
        if created_output:
            (out / "metadata.json").unlink(missing_ok=True)
        raise
    finally:
        ref.close()


def load_panel(path):
    p = Path(path).resolve()
    meta = json.loads((p / "metadata.json").read_text(encoding="utf-8"))
    if (not isinstance(meta, dict) or type(meta.get("schema")) is not int or meta["schema"] != 1
            or not (p / "panel.sqlite").is_file() or not (p / "sites.bed").is_file()):
        raise ValueError("Incomplete or unsupported panel")
    strains = meta.get("strains")
    if (not isinstance(strains, list) or len(strains) < 2
            or any(not isinstance(s, str) or not s.strip() for s in strains)
            or len(strains) != len(set(strains))):
        raise ValueError("Panel strains must contain at least two nonempty, unique names")
    contigs = meta.get("contigs")
    if (not isinstance(contigs, dict) or not contigs
            or any(not primary(c) or type(n) is not int or n < 1 for c, n in contigs.items())
            or len({canonical(c) for c in contigs}) != len(contigs)):
        raise ValueError("Panel contigs must have unique primary chromosomes and positive lengths")
    if meta.get("assembly") not in ("GRCm38", "GRCm39"):
        raise ValueError("Panel assembly must be GRCm38 or GRCm39")
    digest = meta.get("fasta_sha256")
    if not isinstance(digest, str) or not re.fullmatch(r"[0-9a-f]{64}", digest):
        raise ValueError("Panel requires a valid FASTA SHA-256 checksum")
    db = sqlite3.connect((p / "panel.sqlite").as_uri() + "?mode=ro", uri=True)
    try:
        db.execute("SELECT chrom,pos,ref,alt,gt FROM snps LIMIT 0")
        if db.execute("SELECT 1 FROM snps LIMIT 1").fetchone() is None:
            raise ValueError("Panel contains no SNPs")
    except (sqlite3.DatabaseError, ValueError) as error:
        db.close()
        raise ValueError(f"Invalid panel database: {error}") from error
    return meta, db
