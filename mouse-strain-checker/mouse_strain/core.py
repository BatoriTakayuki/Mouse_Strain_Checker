"""Site-balanced, explicitly uncalibrated candidate scoring (stdlib only)."""
from dataclasses import dataclass
from itertools import combinations
import math
from numbers import Real


@dataclass(frozen=True)
class Rules:
    min_depth: int = 5
    min_mapq: int = 30
    min_baseq: int = 25
    min_sites: int = 30
    min_window_sites: int = 10
    min_chrom_windows: int = 2
    min_diagnostic: int = 10
    min_margin: float = 0.1
    max_fit_loss: float = 0.5
    window_size: int = 10_000_000
    depth_cap: int = 20
    error_rate: float = 0.02
    concentration: float = 10.0
    max_other_fraction: float = 0.1
    max_pileup_depth: int = 10000

    def __post_init__(self):
        for key, value in vars(self).items():
            if isinstance(value, bool) or not isinstance(value, Real) or not math.isfinite(value):
                raise ValueError(f"{key} must be a finite number")
        integer_fields = ("min_depth", "min_mapq", "min_baseq", "min_sites", "min_window_sites",
                          "min_chrom_windows", "min_diagnostic", "window_size", "depth_cap",
                          "max_pileup_depth")
        for key in integer_fields:
            if not isinstance(getattr(self, key), int):
                raise ValueError(f"{key} must be an integer")
        for key in set(integer_fields) - {"min_mapq", "min_baseq"}:
            if getattr(self, key) < 1:
                raise ValueError(f"{key} must be positive")
        if not 0 < self.error_rate < 0.5 or self.concentration <= 0:
            raise ValueError("Invalid likelihood parameters")
        if not 0 <= self.max_other_fraction < 1:
            raise ValueError("Invalid other-base fraction")
        if min(self.min_mapq, self.min_baseq, self.min_margin, self.max_fit_loss) < 0:
            raise ValueError("Quality thresholds must be nonnegative")


def site_scores(ref_count, alt_count, rules):
    """Beta-binomial log scores without the shared binomial coefficient.

    Cap effective depth without rounding; use the same overdispersion for
    homozygous and two-strain models. No independence/posterior claim is made.
    """
    if any(isinstance(count, bool) or not isinstance(count, Real)
           or not math.isfinite(count) or count < 0 for count in (ref_count, alt_count)):
        raise ValueError("Allele counts must be finite and nonnegative")
    depth = ref_count + alt_count
    if depth <= 0 or not math.isfinite(depth):
        raise ValueError("At least one finite allele observation is required")
    scale = min(1.0, rules.depth_cap / depth)
    k, n = alt_count * scale, depth * scale
    result = []
    for p in (rules.error_rate, 0.5, 1 - rules.error_rate):
        a, b = p * rules.concentration, (1 - p) * rules.concentration
        result.append(math.lgamma(k + a) + math.lgamma(n - k + b)
                      - math.lgamma(n + a + b) - math.lgamma(a)
                      - math.lgamma(b) + math.lgamma(a + b))
    return result


class Summary:
    def __init__(self, strains, pairs=True):
        self.strains = tuple(strains)
        if len(self.strains) < 2 or any(not isinstance(s, str) or not s for s in self.strains):
            raise ValueError("At least two named candidate strains are required")
        if len(set(self.strains)) != len(self.strains):
            raise ValueError("Candidate strain names must be unique")
        self.models = [(i, i) for i in range(len(self.strains))]
        if pairs:
            self.models += list(combinations(range(len(self.strains)), 2))
        self.scores = [0.0] * len(self.models)
        self.saturated = 0.0
        self.n = 0
        self.first = self.last = None
        self.patterns = {}

    def add(self, pos, genotypes, ref_count, alt_count, rules):
        if not isinstance(genotypes, str) or len(genotypes) != len(self.strains) or set(genotypes) - {"0", "1"}:
            raise ValueError("Expected one binary genotype per candidate strain")
        likelihoods = site_scores(ref_count, alt_count, rules)
        for i, (a, b) in enumerate(self.models):
            self.scores[i] += likelihoods[int(genotypes[a]) + int(genotypes[b])]
        self.saturated += max(likelihoods)
        self.n += 1
        self.first = pos if self.first is None else min(self.first, pos)
        self.last = pos if self.last is None else max(self.last, pos)
        self.patterns[genotypes] = self.patterns.get(genotypes, 0) + 1

    def name(self, model):
        a, b = self.models[model]
        return self.strains[a] if a == b else self.strains[a] + " + " + self.strains[b]

    def _candidate_options(self, rules, order):
        """Ranked panel models not separated from the best by existing rules.

        A model must fit in its own right. Rejecting it also requires both
        enough observed diagnostic sites and enough score separation from the
        best model. These are reference alternatives, not a confidence set or
        a claim that the true strain is included in the panel.
        """
        top = order[0]
        a, b = self.models[top]
        options = []
        for index in order:
            loss = max(0.0, (self.saturated - self.scores[index]) / self.n)
            if loss > rules.max_fit_loss:
                continue
            c, d = self.models[index]
            diagnostic = sum(n for g, n in self.patterns.items()
                             if int(g[a]) + int(g[b]) != int(g[c]) + int(g[d]))
            margin = (self.scores[top] - self.scores[index]) / self.n
            if index == top or diagnostic < rules.min_diagnostic or margin < rules.min_margin:
                options.append(self.name(index))
        return options

    def result(self, rules, window=False):
        row = dict(status="insufficient_sites", call="", candidate="", runner_up="",
                   n_sites=self.n, diagnostic_sites=0, margin_per_site=None,
                   fit_loss_per_site=None, first_snp=self.first, last_snp=self.last,
                   reason="", candidate_options="", candidate_count=0,
                   candidate_options_reason="")
        if not self.n:
            row["status"] = "no_observed_snps"
            row["reason"] = "RNA上で有効なパネルSNPを観測できませんでした。"
            row["candidate_options_reason"] = "有効なSNPがないため、候補を絞れません。"
            return row
        order = sorted(range(len(self.models)), key=lambda i: -self.scores[i])
        top, second = order[:2]
        a, b = self.models[top]
        c, d = self.models[second]
        diagnostic = sum(n for g, n in self.patterns.items()
                         if int(g[a]) + int(g[b]) != int(g[c]) + int(g[d]))
        margin = (self.scores[top] - self.scores[second]) / self.n
        loss = max(0.0, (self.saturated - self.scores[top]) / self.n)
        row.update(candidate=self.name(top), runner_up=self.name(second),
                   diagnostic_sites=diagnostic, margin_per_site=margin,
                   fit_loss_per_site=loss)
        if self.n < (rules.min_window_sites if window else rules.min_sites):
            row["reason"] = f"有効SNPが{self.n}か所で、基準の{rules.min_window_sites if window else rules.min_sites}か所に届きません。"
            row["candidate_options_reason"] = "有効SNP数が基準未満のため、候補を絞れません。"
            return row
        if loss > rules.max_fit_loss:
            row["status"] = "poor_fit"
            row["reason"] = f"候補モデルの適合損失/SNPが{loss:.4g}で、上限{rules.max_fit_loss:g}を超えています。"
            row["candidate_options_reason"] = "最上位モデルも適合度の基準を満たさないため、候補を絞れません。"
        elif diagnostic < rules.min_diagnostic or margin < rules.min_margin:
            row["status"] = "ambiguous"
            reasons = []
            if diagnostic == 0:
                reasons.append("最上位と次点は観測したSNPでは同型で、区別できません（識別SNP 0か所）")
            elif diagnostic < rules.min_diagnostic:
                reasons.append(f"最上位と次点の識別SNPが{diagnostic}か所で、基準{rules.min_diagnostic}か所未満です")
            if margin < rules.min_margin:
                reasons.append(f"スコア差/SNPが{margin:.4g}で、基準{rules.min_margin:g}未満です")
            row["reason"] = "。".join(reasons) + "。"
        else:
            row["status"] = "candidate_supported" if a == b else "mixed_or_heterozygous"
            row["call"] = self.name(top)
            row["reason"] = "有効SNP数・識別SNP数・スコア差・適合度の既定基準を満たしました。"
        if loss <= rules.max_fit_loss:
            options = self._candidate_options(rules, order)
            row.update(candidate_options=" or ".join(options), candidate_count=len(options),
                       candidate_options_reason="候補パネル内で適合度を満たし、最上位との差または識別SNP数が区別の基準に届かないモデルをスコア順に示します。確率や系統の確定を示すものではありません。")
        return row

    def rankings(self):
        if not self.n:
            return []
        return [dict(rank=rank, model=self.name(i), n_sites=self.n,
                     mean_log_score=self.scores[i] / self.n,
                     fit_loss_per_site=max(0.0, (self.saturated - self.scores[i]) / self.n))
                for rank, i in enumerate(sorted(range(len(self.models)), key=lambda i: -self.scores[i]), 1)]


def parse_bases(bases, ref):
    """One token per pileup entry, including deletions and intron skips."""
    if not isinstance(ref, str) or len(ref) != 1 or ref.upper() not in "ACGTN":
        raise ValueError("Expected one pileup reference base")
    result, i = [], 0
    while i < len(bases):
        if bases[i] == "^":
            if i + 2 >= len(bases):
                raise ValueError("Truncated pileup read-start")
            i += 2
        char = bases[i]
        if char not in ".,ACGTNacgtn*#<>":
            raise ValueError(f"Unsupported pileup token: {char!r}")
        result.append(ref.upper() if char in ".," else char.upper())
        i += 1
        while i < len(bases) and bases[i] in "+-":
            operation = bases[i]
            i += 1
            start = i
            while i < len(bases) and bases[i] in "0123456789":
                i += 1
            if i == start:
                raise ValueError("Invalid pileup indel")
            length = int(bases[start:i])
            if length < 1:
                raise ValueError("Invalid pileup indel length")
            if i + length > len(bases):
                raise ValueError("Truncated pileup indel")
            allowed = "ACGTNacgtn*#" if operation == "+" else "ACGTNacgtn"
            if any(base not in allowed for base in bases[i:i + length]):
                raise ValueError("Invalid pileup indel sequence")
            i += length
        if i < len(bases) and bases[i] == "$":
            i += 1
    return result


def count_fragments(fields, rules, qc=None):
    """Expected fields: CHROM POS REF DP BASES BQ QNAME RG NH.

    Agreeing mates count once. Conflicting mates are both discarded.
    Group by (RG, QNAME) so names repeated across SRRs do not collide.
    """
    if len(fields) != 9:
        raise ValueError("Expected 9 mpileup columns (QNAME,RG,NH); samtools >=1.13 required")
    depth = int(fields[3])
    if depth < 0:
        raise ValueError("Pileup depth must be nonnegative")
    if depth == 0:
        if any(field != "*" for field in fields[4:]):
            raise ValueError("Zero-depth pileup columns must be empty placeholders")
        return {}
    bases = parse_bases(fields[4], fields[2])
    qualities = [ord(c) - 33 for c in fields[5]]
    if any(quality < 0 or quality > 93 for quality in qualities):
        raise ValueError("Invalid pileup base-quality encoding")
    names, groups, nh = [s.split(",") for s in fields[6:9]]
    if not (len(bases) == len(qualities) == len(names) == len(groups) == len(nh) == depth):
        raise ValueError("Pileup column lengths disagree; refusing misassigned base evidence")
    fragments = {}
    for base, quality, name, group, hits in zip(bases, qualities, names, groups, nh):
        if qc is not None:
            qc["pileup_entries"] += 1
            if hits == "*":
                qc["missing_nh_entries"] += 1
            elif hits != "1":
                qc["nonunique_nh_entries"] += 1
        if quality < rules.min_baseq or hits != "1" or base not in "ACGT":
            continue
        key = (group, name)
        if key not in fragments:
            fragments[key] = base
        elif fragments[key] != base:
            fragments[key] = None
    counts = {}
    for base in fragments.values():
        if base:
            counts[base] = counts.get(base, 0) + 1
    return counts
