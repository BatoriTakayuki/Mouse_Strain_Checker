"""Input-integrity regressions for candidate scoring and pileup evidence."""
from collections import Counter
from dataclasses import replace
import math
import unittest

from mouse_strain.core import Rules, Summary, count_fragments, parse_bases, site_scores


class CoreInputTests(unittest.TestCase):
    def test_rules_reject_fractional_counts_and_nonnumeric_thresholds(self):
        for key in ("min_depth", "min_mapq", "min_baseq", "min_sites", "min_window_sites",
                    "min_chrom_windows", "min_diagnostic", "window_size", "depth_cap", "max_pileup_depth"):
            for value in (1.5, True):
                with self.subTest(key=key, value=value), self.assertRaises(ValueError):
                    replace(Rules(), **{key: value})
        for value in ("0.1", None, float("nan"), float("inf"), True):
            with self.subTest(value=value), self.assertRaises(ValueError):
                replace(Rules(), min_margin=value)

    def test_invalid_counts_cannot_enter_summary(self):
        rules = Rules()
        summary = Summary(["A", "B"])
        for ref, alt in ((0, 0), (-1, 10), (1, float("nan")), (float("inf"), 1)):
            with self.subTest(ref=ref, alt=alt), self.assertRaises(ValueError):
                summary.add(1, "01", ref, alt, rules)
        self.assertEqual(summary.n, 0)
        self.assertEqual(summary.scores, [0.0, 0.0, 0.0])
        self.assertTrue(all(math.isfinite(score) for score in site_scores(3, 7, rules)))

    def test_summary_checks_strain_identity_and_genotype_mapping(self):
        for strains in ([], ["A"], ["A", "A"], ["A", ""], ["A", None]):
            with self.subTest(strains=strains), self.assertRaises(ValueError):
                Summary(strains)
        names = ["A", "B"]
        summary = Summary(names)
        names[0] = "changed"
        self.assertEqual(summary.name(0), "A")
        for genotypes in ("0", "010", "02", (0, 1)):
            with self.subTest(genotypes=genotypes), self.assertRaises(ValueError):
                summary.add(1, genotypes, 10, 0, Rules())
        self.assertEqual(summary.n, 0)


class PileupIntegrityTests(unittest.TestCase):
    def fields(self, bases, depth=1, quality="I"):
        return ["1", "10", "A", str(depth), bases, quality, "read", "run", "1"]

    def test_annotation_must_be_attached_to_a_pileup_entry(self):
        # These formerly produced one apparently valid base despite invalid syntax.
        for bases in ("+1tA", "$A", "A$$", "A^]", "A+0", "A+1$", "A-1*", "A+\u0661t"):
            with self.subTest(bases=bases), self.assertRaises(ValueError):
                count_fragments(self.fields(bases), Rules())

    def test_all_supported_entry_types_keep_quality_and_name_alignment(self):
        bases = "^].+2a*,$G-1t$*#<>n"
        self.assertEqual(parse_bases(bases, "A"), ["A", "A", "G", "*", "#", "<", ">", "N"])
        fields = ["1", "10", "A", "8", bases, "I" * 8,
                  "a,a,b,c,d,e,f,g", ",".join(["run"] * 8), ",".join(["1"] * 8)]
        qc = Counter()
        self.assertEqual(count_fragments(fields, Rules(), qc), {"A": 1, "G": 1})
        self.assertEqual(qc["pileup_entries"], 8)

    def test_zero_depth_requires_documented_empty_columns(self):
        empty = ["1", "10", "A", "0", "*", "*", "*", "*", "*"]
        self.assertEqual(count_fragments(empty, Rules()), {})
        for column in range(4, 9):
            fields = empty.copy()
            fields[column] = "A"
            with self.subTest(column=column), self.assertRaises(ValueError):
                count_fragments(fields, Rules())

    def test_malformed_base_quality_cannot_pass_quality_filter(self):
        for quality in (" ", "\x7f", "\u00e9"):
            with self.subTest(quality=quality), self.assertRaises(ValueError):
                count_fragments(self.fields("A", quality=quality), Rules())


if __name__ == "__main__":
    unittest.main()
