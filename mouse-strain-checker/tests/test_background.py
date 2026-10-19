"""Background groups cannot replace evidence required for exact strain calls."""
from dataclasses import replace
import unittest

from mouse_strain.background import background_result
from mouse_strain.core import Rules, Summary


class BackgroundResultTests(unittest.TestCase):
    def setUp(self):
        self.rules = replace(Rules(), min_sites=20, min_window_sites=10, min_diagnostic=10)

    def summary(self, strains, genotypes, counts=(10, 0), n=30, pairs=True):
        summary = Summary(strains, pairs=pairs)
        for pos in range(n):
            summary.add(pos + 1, genotypes, *counts, self.rules)
        return summary

    def test_exact_b6_tie_can_have_supported_background_group(self):
        summary = self.summary(["C57BL_6J", "C57BL_6NJ", "A"], "001")
        self.assertEqual(summary.result(self.rules)["status"], "ambiguous")
        result = background_result(summary, self.rules)
        self.assertEqual(result["background_status"], "candidate_group_supported")
        self.assertEqual(result["background_call"], "C57BL/6 (6J/6NJ)")
        self.assertEqual(result["background_diagnostic_sites"], 30)
        self.assertEqual(set(result["background_members"].split(" | ")),
                         {"C57BL_6J", "C57BL_6NJ", "C57BL_6J + C57BL_6NJ"})

    def test_mixed_group_is_only_a_group_label(self):
        summary = self.summary(["129", "C57BL_6J", "C57BL_6NJ", "A"], "0110", counts=(5, 5))
        result = background_result(summary, self.rules)
        # 129 + B6 and A + B6 are observationally identical in this fixture.
        self.assertEqual(result["background_status"], "ambiguous")
        self.assertEqual(result["background_call"], "")
        self.assertIn("観測SNP上で", result["background_reason"])

    def test_129_plus_b6_group_keeps_both_exact_pair_members(self):
        summary = self.summary(["129", "C57BL_6J", "C57BL_6NJ"], "011", counts=(5, 5))
        result = background_result(summary, self.rules)
        self.assertEqual(summary.result(self.rules)["call"], "")
        self.assertEqual(result["background_status"], "candidate_group_supported")
        self.assertEqual(result["background_call"], "129 + C57BL/6 (6J/6NJ)")
        self.assertEqual(set(result["background_members"].split(" | ")),
                         {"129 + C57BL_6J", "129 + C57BL_6NJ"})

    def test_background_requires_margin_to_other_groups(self):
        summary = self.summary(["C57BL_6J", "C57BL_6NJ", "A"], "001", counts=(6, 4))
        result = background_result(summary, replace(self.rules, min_margin=100))
        self.assertEqual(result["background_status"], "ambiguous")
        self.assertEqual(result["background_call"], "")
        self.assertIn("スコア差", result["background_reason"])

    def test_unknown_substrains_remain_separate(self):
        summary = self.summary(["C57BL_6J", "C57BL_6N", "A"], "001", pairs=False)
        result = background_result(summary, self.rules)
        self.assertEqual(result["background_status"], "ambiguous")
        self.assertEqual(result["background_runner_up"], "C57BL_6N")
        self.assertEqual(result["background_candidate"], "C57BL_6J")
        self.assertEqual(result["background_members"], "C57BL_6J")

    def test_one_selected_b6_substrain_keeps_exact_label(self):
        for strain in ("C57BL_6J", "C57BL_6NJ"):
            with self.subTest(strain=strain):
                summary = self.summary([strain, "A"], "01", pairs=False)
                result = background_result(summary, self.rules)
                self.assertEqual(result["background_status"], "candidate_group_supported")
                self.assertEqual(result["background_candidate"], strain)
                self.assertEqual(result["background_call"], strain)
                self.assertEqual(result["background_members"], strain)

    def test_pair_with_one_selected_b6_substrain_keeps_exact_label(self):
        summary = self.summary(["129", "C57BL_6NJ"], "01", counts=(5, 5))
        result = background_result(summary, self.rules)
        self.assertEqual(result["background_status"], "candidate_group_supported")
        self.assertEqual(result["background_call"], "129 + C57BL_6NJ")
        self.assertEqual(result["background_members"], "129 + C57BL_6NJ")

    def test_all_outside_models_require_diagnostic_sites(self):
        summary = Summary(["C57BL_6J", "C57BL_6NJ", "A", "B"], pairs=False)
        for pos in range(10):
            summary.add(pos + 1, "0010", 5, 0, self.rules)
        for pos in range(10, 19):
            summary.add(pos + 1, "0001", 20, 0, self.rules)
        for pos in range(19, 30):
            summary.add(pos + 1, "0000", 10, 0, self.rules)
        result = background_result(summary, self.rules)
        self.assertEqual(result["background_runner_up"], "A")
        # The score runner-up has 10 diagnostic SNPs, but B has only nine.
        self.assertEqual(result["background_diagnostic_sites"], 9)
        self.assertEqual(result["background_status"], "ambiguous")
        self.assertIn("SNPが9個", result["background_reason"])

    def test_few_diagnostic_sites_hold_even_with_large_margin(self):
        summary = self.summary(["C57BL_6J", "C57BL_6NJ", "A"], "000", n=25, pairs=False)
        for pos in range(25, 30):
            summary.add(pos + 1, "001", 10, 0, self.rules)
        result = background_result(summary, replace(self.rules, min_margin=0))
        self.assertEqual(result["background_diagnostic_sites"], 5)
        self.assertEqual(result["background_status"], "ambiguous")
        self.assertIn("SNPが5個", result["background_reason"])

    def test_insufficient_sites_and_window_threshold_match_exact_rules(self):
        summary = self.summary(["C57BL_6J", "C57BL_6NJ", "A"], "001", n=15)
        result = background_result(summary, self.rules)
        self.assertEqual(result["background_status"], "insufficient_sites")
        self.assertEqual(result["background_call"], "")
        self.assertEqual(background_result(summary, self.rules, window=True)["background_status"],
                         "candidate_group_supported")

    def test_poor_fit_keeps_background_on_hold(self):
        summary = self.summary(["C57BL_6J", "C57BL_6NJ", "A"], "001", counts=(5, 5), pairs=False)
        result = background_result(summary, self.rules)
        self.assertEqual(result["background_status"], "poor_fit")
        self.assertEqual(result["background_call"], "")

    def test_no_outside_group_cannot_be_supported(self):
        summary = self.summary(["C57BL_6J", "C57BL_6NJ"], "01")
        result = background_result(summary, self.rules)
        self.assertEqual(result["background_status"], "ambiguous")
        self.assertIn("比較できる別の背景候補群がない", result["background_reason"])
        self.assertEqual(result["background_call"], "")
        self.assertIsNone(result["background_margin_per_site"])

    def test_no_observations_has_same_columns_and_no_candidate(self):
        summary = Summary(["C57BL_6J", "C57BL_6NJ", "A"])
        result = background_result(summary, self.rules)
        populated = background_result(self.summary(summary.strains, "001"), self.rules)
        self.assertEqual(result.keys(), populated.keys())
        self.assertEqual(result["background_status"], "no_observed_snps")
        self.assertEqual(result["background_candidate"], "")
        self.assertEqual(result["background_call"], "")


if __name__ == "__main__":
    unittest.main()
