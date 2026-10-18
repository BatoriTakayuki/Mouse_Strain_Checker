"""Site-balanced, explicitly uncalibrated candidate scoring (stdlib only)."""
from dataclasses import dataclass
from itertools import combinations
import math


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
        if any(not math.isfinite(value) for value in vars(self).values()):
            raise ValueError("Thresholds must be finite")
        for key in ("min_depth", "min_sites", "min_window_sites", "min_chrom_windows", "min_diagnostic",
                    "window_size", "depth_cap", "max_pileup_depth"):
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
    depth = ref_count + alt_count
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
        self.strains = strains
        self.models = [(i, i) for i in range(len(strains))]
        if pairs:
            self.models += list(combinations(range(len(strains)), 2))
        self.scores = [0.0] * len(self.models)
        self.saturated = 0.0
        self.n = 0
        self.first = self.last = None
        self.patterns = {}

    def add(self, pos, genotypes, ref_count, alt_count, rules):
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

    def result(self, rules, window=False):
        row = dict(status="insufficient_sites", call="", candidate="", runner_up="",
                   n_sites=self.n, diagnostic_sites=0, margin_per_site=None,
                   fit_loss_per_site=None, first_snp=self.first, last_snp=self.last)
        if not self.n:
            row["status"] = "no_observed_snps"
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
            return row
        if loss > rules.max_fit_loss:
            row["status"] = "poor_fit"
        elif diagnostic < rules.min_diagnostic or margin < rules.min_margin:
            row["status"] = "ambiguous"
        else:
            row["status"] = "candidate_supported" if a == b else "mixed_or_heterozygous"
            row["call"] = self.name(top)
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
    result, i = [], 0
    while i < len(bases):
        char = bases[i]
        if char == "^":
            if i + 1 >= len(bases):
                raise ValueError("Truncated pileup read-start")
            i += 2
        elif char == "$":
            i += 1
        elif char in "+-":
            i += 1
            start = i
            while i < len(bases) and bases[i].isdigit():
                i += 1
            if i == start:
                raise ValueError("Invalid pileup indel")
            length = int(bases[start:i])
            i += length
            if i > len(bases):
                raise ValueError("Truncated pileup indel")
        elif char in ".,ACGTNacgtn*#<>":
            result.append(ref.upper() if char in ".," else char.upper())
            i += 1
        else:
            raise ValueError(f"Unsupported pileup token: {char!r}")
    return result


def count_fragments(fields, rules, qc=None):
    """Expected fields: CHROM POS REF DP BASES BQ QNAME RG NH.

    Agreeing mates count once. Conflicting mates are both discarded.
    Group by (RG, QNAME) so names repeated across SRRs do not collide.
    """
    if len(fields) != 9:
        raise ValueError("Expected 9 mpileup columns (QNAME,RG,NH); samtools >=1.18 required")
    if int(fields[3]) == 0:
        return {}
    bases = parse_bases(fields[4], fields[2])
    qualities = [ord(c) - 33 for c in fields[5]]
    names, groups, nh = [s.split(",") for s in fields[6:9]]
    if not (len(bases) == len(qualities) == len(names) == len(groups) == len(nh) == int(fields[3])):
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
