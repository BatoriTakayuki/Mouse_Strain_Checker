"""Explicit alternatives retain all plausible models without forcing a call."""
from dataclasses import replace
import unittest

from mouse_strain.core import Rules, Summary


class CandidateOptionsTests(unittest.TestCase):
    def setUp(self):
        self.rules = Rules()

    def summary(self, strains, genotypes, counts=(10, 0), n=30, pairs=False):
        summary = Summary(strains, pairs=pairs)
        for pos in range(1, n + 1):
            summary.add(pos, genotypes, *counts, self.rules)
        return summary

    def test_three_indistinguishable_models_are_all_reported(self):
        summary = self.summary(["A", "B", "C", "D"], "0001")
        result = summary.result(self.rules)
        self.assertEqual(result["status"], "ambiguous")
        self.assertEqual(result["call"], "")
        self.assertEqual(result["candidate_options"], "A or B or C")
        self.assertEqual(result["candidate_count"], 3)
        self.assertNotIn("D", result["candidate_options"])

    def test_two_strain_model_names_keep_plus_and_options_use_or(self):
        summary = self.summary(["A", "B", "C"], "011", (5, 5), pairs=True)
        result = summary.result(self.rules)
        self.assertEqual(result["candidate_options"], "A + B or A + C")
        self.assertEqual(result["candidate_count"], 2)
        self.assertEqual(result["status"], "ambiguous")

    def test_diagnostic_requirement_can_retain_a_lower_scoring_third_model(self):
        summary = self.summary(["A", "B", "C"], "000", n=29)
        summary.add(30, "011", 10, 0, self.rules)
        result = summary.result(self.rules)
        self.assertGreater(result["margin_per_site"], self.rules.min_margin)
        self.assertEqual(result["diagnostic_sites"], 1)
        self.assertEqual(result["candidate_options"], "A or B or C")
        self.assertEqual(result["candidate_count"], 3)

    def test_small_score_gap_retains_alternative_even_with_enough_diagnostic_sites(self):
        rules = replace(self.rules, min_diagnostic=2)
        summary = self.summary(["A", "B"], "00", n=398)
        for pos in (399, 400):
            summary.add(pos, "01", 10, 0, rules)
        result = summary.result(rules)
        self.assertEqual(result["diagnostic_sites"], rules.min_diagnostic)
        self.assertLess(result["margin_per_site"], rules.min_margin)
        self.assertEqual(result["candidate_options"], "A or B")

    def test_individually_poor_fitting_alternative_is_excluded(self):
        summary = self.summary(["A", "B"], "00", n=29)
        summary.add(30, "01", 20, 0, self.rules)
        result = summary.result(self.rules)
        # One diagnostic SNP does not support the strict call. It also does
        # not justify presenting an individually poor-fitting B as an option.
        self.assertEqual(result["status"], "ambiguous")
        self.assertEqual(result["call"], "")
        self.assertEqual(result["candidate_options"], "A")
        self.assertEqual(result["candidate_count"], 1)

    def test_supported_result_retains_its_call_and_single_option(self):
        summary = self.summary(["A", "B"], "01")
        scores_before = list(summary.scores)
        rankings_before = summary.rankings()
        result = summary.result(self.rules)
        self.assertEqual(result["status"], "candidate_supported")
        self.assertEqual(result["call"], "A")
        self.assertEqual(result["candidate_options"], "A")
        self.assertEqual(result["candidate_count"], 1)
        self.assertEqual(summary.scores, scores_before)
        self.assertEqual(summary.rankings(), rankings_before)

    def test_no_data_insufficient_sites_and_poor_fit_do_not_narrow_candidates(self):
        cases = [
            (Summary(["A", "B"]), "no_observed_snps"),
            (self.summary(["A", "B"], "01", n=29), "insufficient_sites"),
            (self.summary(["A", "B"], "01", (5, 5)), "poor_fit"),
        ]
        for summary, status in cases:
            with self.subTest(status=status):
                result = summary.result(self.rules)
                self.assertEqual(result["status"], status)
                self.assertEqual(result["candidate_options"], "")
                self.assertEqual(result["candidate_count"], 0)
                self.assertTrue(result["candidate_options_reason"])

    def test_window_threshold_controls_whether_options_can_be_shown(self):
        summary = self.summary(["A", "B", "C"], "001", n=15)
        self.assertEqual(summary.result(self.rules)["candidate_count"], 0)
        result = summary.result(self.rules, window=True)
        self.assertEqual(result["candidate_options"], "A or B")
        self.assertEqual(result["candidate_count"], 2)

    def test_empty_and_populated_results_have_the_same_columns(self):
        empty = Summary(["A", "B"]).result(self.rules)
        populated = self.summary(["A", "B"], "01").result(self.rules)
        self.assertEqual(empty.keys(), populated.keys())


if __name__ == "__main__":
    unittest.main()
