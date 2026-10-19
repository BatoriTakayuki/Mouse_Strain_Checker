"""Conservative, explicitly uncalibrated sex-associated RNA marker evidence.

Coordinates are 1-based closed GRCm39 canonical exons, obtained from the
official Ensembl REST lookup/symbol endpoint on 2026-10-19.  No annotation
download or network request is needed during analysis.  These are expression
signals, not a karyotype, gonadal sex, or a contamination estimate.
"""
from pathlib import Path
import re
import subprocess
import tempfile


ASSEMBLY_LENGTHS = {"5": 151758149, "6": 149588044,
                    "X": 169476592, "Y": 91455967}
MAX_FRAGMENTS = 10000
THRESHOLDS = dict(y_gene_fragments=5, y_gene_starts=3, y_genes=2,
                  y_total_fragments=20, xist_fragments=20, xist_starts=10,
                  xist_spliced_fragments=3, xist_junctions=2,
                  control_gene_fragments=20, control_genes=2,
                  control_total_fragments=100, splice_anchor_bases=8)
SOURCES = [
    "https://rest.ensembl.org/documentation/info/lookup",
    "https://rest.ensembl.org/info/assembly/mus_musculus?content-type=application%2Fjson",
    "https://pmc.ncbi.nlm.nih.gov/articles/PMC4058436/",
    "https://pubmed.ncbi.nlm.nih.gov/10192391/",
    "https://pubmed.ncbi.nlm.nih.gov/26739568/",
]
NOTES = [
    "RNA発現からの参考推定です。性染色体構成や個体の性別を確定する検査ではありません。",
    "Yマーカーが見えないだけではメスと判定しません。組織・培養細胞・遺伝子操作・混合試料では判定できないことがあります。",
    "Xist候補は対照遺伝子の発現とXist固有のスプライス接合を併用します。Xist領域のリードだけでは反対鎖のTsixと区別できません。",
    "閾値は保守的な実装規則で、精度や確率を校正したものではありません。",
]


def _marker(gene, chrom, start, end, gene_id, transcript, exons, role):
    return dict(gene=gene, chromosome=chrom, start=start, end=end,
                gene_id=gene_id, transcript=transcript, exons=exons, role=role,
                source="https://rest.ensembl.org/lookup/symbol/mus_musculus/"
                       + gene + "?content-type=application%2Fjson&expand=1")


MARKERS = [
    _marker("Ddx3y", "Y", 1260715, 1286663, "ENSMUSG00000069045", "ENSMUST00000091190.13",
            [(1260771,1263367),(1263570,1263700),(1263796,1263949),(1264841,1264958),
             (1265028,1265209),(1265575,1265719),(1265829,1265973),(1266359,1266519),
             (1266595,1266693),(1267103,1267188),(1269705,1269840),(1278897,1278996),
             (1279423,1279578),(1280310,1280439),(1282901,1282948),(1284322,1284382),
             (1286479,1286663)], "Y"),
    _marker("Uty", "Y", 1096855, 1245759, "ENSMUSG00000068457", "ENSMUST00000069309.15",
            [(1096855,1098573),(1099582,1099752),(1102765,1102891),(1129870,1130011),
             (1131036,1131223),(1134872,1134986),(1136848,1136996),(1137658,1137732),
             (1137807,1137871),(1151948,1152153),(1152769,1152874),(1154155,1154284),
             (1157907,1158634),(1168090,1168185),(1169123,1169257),(1169965,1170184),
             (1170904,1171002),(1174728,1174854),(1176491,1176584),(1186774,1186808),
             (1189336,1189390),(1197206,1197326),(1199455,1199513),(1224497,1224546),
             (1239711,1239819),(1245156,1245219),(1245407,1245718)], "Y"),
    _marker("Kdm5d", "Y", 897569, 956786, "ENSMUSG00000056673", "ENSMUST00000055032.15",
            [(897788,898115),(898550,898627),(899428,899550),(899787,899957),
             (900479,900613),(910448,910586),(910686,910849),(914025,914183),
             (916641,916760),(916846,917004),(921402,921583),(921674,921836),
             (927007,927126),(927245,927439),(927710,927891),(927970,928094),
             (928237,928384),(936876,936972),(937813,938171),(939819,939957),
             (940599,940778),(940905,941042),(941227,941823),(942231,942312),
             (942442,942641),(942824,946316)], "Y"),
    _marker("Eif2s3y", "Y", 1008409, 1047462, "ENSMUSG00000069049", "ENSMUST00000091197.5",
            [(1010543,1010685),(1011260,1011323),(1011840,1011967),(1012048,1012169),
             (1013379,1013473),(1014634,1014792),(1016057,1016191),(1016590,1016684),
             (1017364,1017508),(1019986,1020155),(1023302,1023474),(1028188,1031318)], "Y"),
    _marker("Xist", "X", 102503972, 102537573, "ENSMUSG00000086503", "ENSMUST00000127786.4",
            [(102503972,102511642),(102512424,102512578),(102512906,102513052),
             (102513197,102513407),(102514156,102514287),(102514437,102514527),
             (102517322,102526860)], "Xist"),
    _marker("Actb", "5", 142888870, 142892509, "ENSMUSG00000029580", "ENSMUST00000100497.11",
            [(142888870,142889696),(142889822,142890003),(142890099,142890537),
             (142890992,142891231),(142891319,142891447),(142892407,142892509)], "control"),
    _marker("Gapdh", "6", 125135289, 125143518, "ENSMUSG00000057666", "ENSMUST00000118875.8",
            [(125138817,125139064),(125139175,125139356),(125139441,125139671),
             (125139806,125140003),(125140102,125140399),(125142234,125142274),
             (125142515,125142736)], "control"),
    _marker("Rplp0", "5", 115694515, 115701793, "ENSMUSG00000067274", "ENSMUST00000086519.13",
            [(115697526,115697598),(115697846,115697947),(115698800,115699063),
             (115699158,115699304),(115699392,115699577),(115700485,115700625),
             (115701342,115701793)], "control"),
]


def unavailable_sex(reason, assembly=None):
    return dict(schema=1, status="unavailable", label="判定できません（性別推定未実施）",
                reason=str(reason), method="GRCm39 RNA marker evidence v1", assembly=assembly,
                markers=[], evidence={}, thresholds=dict(THRESHOLDS),
                sources=list(SOURCES), notes=list(NOTES))


def _canonical(name):
    return name[3:] if name.startswith("chr") else name


def _valid_contigs(contigs):
    if not isinstance(contigs, dict):
        raise ValueError("Sex inference requires reference chromosome lengths")
    normalized = {}
    for name, length in contigs.items():
        key = _canonical(name)
        if key in normalized or type(length) is not int or length < 1:
            raise ValueError("Invalid or duplicate reference chromosomes")
        normalized[key] = (name, length)
    for name, length in ASSEMBLY_LENGTHS.items():
        if name not in normalized or normalized[name][1] != length:
            raise ValueError("Sex markers require complete GRCm39 chromosomes 5, 6, X and Y")
    return normalized


def _alignment_parts(position, cigar):
    """Return covered reference blocks and well-anchored skipped junctions."""
    tokens = re.findall(r"([1-9][0-9]*)([MIDNSHP=X])", cigar)
    if not tokens or "".join(n + op for n, op in tokens) != cigar:
        raise ValueError("Invalid SAM CIGAR in sex marker evidence")
    blocks, junctions = [], []
    for i, (size, op) in enumerate(tokens):
        size = int(size)
        if op in "M=X":
            blocks.append((position, position + size - 1))
            position += size
        elif op in "DN":
            if (op == "N" and 0 < i < len(tokens) - 1
                    and tokens[i - 1][1] in "M=X" and tokens[i + 1][1] in "M=X"
                    and min(int(tokens[i - 1][0]), int(tokens[i + 1][0]))
                    >= THRESHOLDS["splice_anchor_bases"]):
                junctions.append((position - 1, position + size))
            position += size
    return blocks, junctions


def count_marker_alignments(lines, min_mapq=30):
    """Count unique RG+QNAME fragments overlapping annotated exonic bases."""
    state = {m["gene"]: dict(fragments=set(), starts=set(), spliced=set(), junctions=set())
             for m in MARKERS}
    for line in lines:
        if not line.strip() or line.startswith("@"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 11:
            raise ValueError("Malformed SAM in sex marker evidence")
        flag, position, mapq = int(fields[1]), int(fields[3]), int(fields[4])
        if flag & 0xF04 or mapq < min_mapq or position < 1:
            continue
        tags = fields[11:]
        # HISAT2 uses NH; MAPQ alone is not evidence of a unique alignment.
        nh = [t for t in tags if t.startswith("NH:")]
        if nh != ["NH:i:1"]:
            continue
        readgroups = [t[5:] for t in tags if t.startswith("RG:Z:")]
        if len(readgroups) > 1:
            raise ValueError("Duplicate SAM read-group tags")
        fragment = (readgroups[0] if readgroups else "", fields[0])
        blocks, junctions = _alignment_parts(position, fields[5])
        chrom = _canonical(fields[2])
        for marker in MARKERS:
            if chrom != marker["chromosome"]:
                continue
            exons = marker["exons"]
            if not any(a <= d and c <= b for a, b in blocks for c, d in exons):
                continue
            evidence = state[marker["gene"]]
            # Bound each counter independently. Coordinate-sorted output can
            # reach the expression cap in an early exon before any splice
            # junction appears; later alignments must still supply junctions.
            if len(evidence["fragments"]) < MAX_FRAGMENTS:
                evidence["fragments"].add(fragment)
            if len(evidence["starts"]) < MAX_FRAGMENTS:
                evidence["starts"].add((position, bool(flag & 16)))
            known = {(exons[i][1], exons[i + 1][0]) for i in range(len(exons) - 1)}
            matched = known.intersection(junctions)
            if matched:
                if len(evidence["spliced"]) < MAX_FRAGMENTS:
                    evidence["spliced"].add(fragment)
                evidence["junctions"].update(matched)
    return [dict(marker, fragments=len(state[marker["gene"]]["fragments"]),
                 distinct_starts=len(state[marker["gene"]]["starts"]),
                 spliced_fragments=len(state[marker["gene"]]["spliced"]),
                 splice_junctions=len(state[marker["gene"]]["junctions"]),
                 count_capped=len(state[marker["gene"]]["fragments"]) >= MAX_FRAGMENTS)
            for marker in MARKERS]


def classify_sex(markers, assembly="GRCm39", min_mapq=30):
    """Classify marker patterns; abstain on conflicting or absent evidence."""
    t = THRESHOLDS
    y = [m for m in markers if m["role"] == "Y"]
    controls = [m for m in markers if m["role"] == "control"]
    xist = next(m for m in markers if m["gene"] == "Xist")
    supported_y = sum(m["fragments"] >= t["y_gene_fragments"]
                      and m["distinct_starts"] >= t["y_gene_starts"] for m in y)
    y_total = sum(m["fragments"] for m in y)
    control_genes = sum(m["fragments"] >= t["control_gene_fragments"] for m in controls)
    control_total = sum(m["fragments"] for m in controls)
    adequate = control_genes >= t["control_genes"] and control_total >= t["control_total_fragments"]
    xist_locus = xist["fragments"] >= t["xist_fragments"] and xist["distinct_starts"] >= t["xist_starts"]
    xist_specific = (xist_locus and xist["spliced_fragments"] >= t["xist_spliced_fragments"]
                     and xist["splice_junctions"] >= t["xist_junctions"])
    male = supported_y >= t["y_genes"] and y_total >= t["y_total_fragments"]
    status, label = "inconclusive", "判定保留（性別を推定できません）"
    if male and xist_specific:
        reason = "複数のYマーカーとXist固有のスプライス接合を伴う強い発現が併存しています。混合試料や細胞状態などを区別できません。"
    elif male:
        status, label = "male_likely", "オスの可能性が高い（RNA推定）"
        reason = "非PAR領域の複数のY遺伝子に、独立した開始位置を持つ一意マッピング断片を認めました。"
        if xist_locus:
            reason += "Xist領域にもリードがありますが、固有スプライス接合が不足し、XistとTsixを区別できません。"
    elif adequate and xist_specific and y_total == 0:
        status, label = "female_likely", "メスの可能性が高い（RNA推定）"
        reason = "対照遺伝子の発現を確認し、Xistの発現と複数の固有スプライス接合を検出しました。Yマーカーは検出されませんでした。"
    elif not adequate:
        reason = "対照遺伝子の発現量が不足しています。Yマーカーの未検出から性別を判断できません。"
    elif y_total:
        reason = "Yマーカーの断片はありますが、複数遺伝子・独立開始位置の基準を満たさず、性別を推定できません。"
    elif xist_locus and not xist_specific:
        reason = "Xist領域に発現がありますが、固有スプライス接合が不足し、Tsixなどとの区別が不十分です。"
    else:
        reason = "YマーカーとXistの確かな発現が不足しています。Yの未検出だけではメスと判定できません。"
    return dict(schema=1, status=status, label=label, reason=reason,
                method="GRCm39 RNA marker evidence v1", assembly=assembly,
                markers=markers, min_mapq=min_mapq, min_nh=1, excluded_flags="0xF04",
                fragment_key="RG + QNAME", annotation_checked_at="2026-10-19",
                count_cap_per_marker=MAX_FRAGMENTS,
                evidence=dict(y_genes_supported=supported_y, y_fragments=y_total,
                              xist_fragments=xist["fragments"],
                              xist_spliced_fragments=xist["spliced_fragments"],
                              xist_junctions=xist["splice_junctions"],
                              control_genes_supported=control_genes,
                              control_fragments=control_total),
                thresholds=dict(t), sources=list(SOURCES), notes=list(NOTES))


def collect_sex_evidence(bam, assembly, contigs, min_mapq=30):
    """Read only eight indexed BAM gene regions, with a bounded subprocess.

    Raises on missing tools/index, malformed data, unsupported references or
    failed commands, so callers can record an explicit unavailable result.
    Marker coordinates are never applied to mm38 or an unchecked reference.
    """
    if assembly != "GRCm39":
        raise ValueError("Sex marker inference currently supports GRCm39 only")
    if type(min_mapq) is not int or not 0 <= min_mapq <= 255:
        raise ValueError("min_mapq must be an integer in 0..255")
    _valid_contigs(contigs)
    bam = Path(bam)
    if not bam.is_file():
        raise ValueError("BAM for sex inference is missing")
    header = subprocess.run(["samtools", "view", "-H", str(bam)], check=True,
                            capture_output=True, text=True, encoding="utf-8", timeout=30)
    sequences = {}
    for line in header.stdout.splitlines():
        if line.startswith("@SQ\t"):
            fields = dict(f.split(":", 1) for f in line.split("\t")[1:] if ":" in f)
            name = fields.get("SN", "")
            if name in sequences:
                raise ValueError("Duplicate BAM sequence names")
            sequences[name] = int(fields.get("LN", 0))
    actual = _valid_contigs(sequences)
    regions = [f"{actual[m['chromosome']][0]}:{m['start']}-{m['end']}" for m in MARKERS]
    # Disk streaming avoids retaining unbounded high-expression SAM output.
    with tempfile.TemporaryFile(mode="w+t", encoding="utf-8") as output:
        subprocess.run(["samtools", "view", "-M", "-q", str(min_mapq), "-F", "3844",
                        str(bam), *regions], stdout=output, stderr=subprocess.PIPE,
                       text=True, encoding="utf-8", check=True, timeout=180)
        output.seek(0)
        markers = count_marker_alignments(output, min_mapq)
    return classify_sex(markers, assembly, min_mapq)
