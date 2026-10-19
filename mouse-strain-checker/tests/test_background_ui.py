"""Keep exact calls separate from useful, explicitly provisional background evidence."""
import csv
import json
from pathlib import Path
import shutil
import subprocess
import tempfile
from types import SimpleNamespace
import unittest

from mouse_strain.webapp import result_data


STATIC = Path(__file__).resolve().parents[1] / "mouse_strain" / "static"


JS_CHECKS = r"""
const assert = require('node:assert/strict');
const ui = require(process.argv[1]);
const palette = new Map([['129S1_SvImJ', '#112233'], ['C57BL/6 (6J/6NJ)', '#445566']]);
const held = {n_sites:'1247', status:'ambiguous', call:'', candidate:'129S1_SvImJ + C57BL_6NJ',
  background_status:'ambiguous', background_candidate:'129S1_SvImJ + C57BL/6 (6J/6NJ)',
  background_call:'', reason:'6J/6NJを区別するSNPが不足しています。', background_reason:'群間の差が基準未満です。'};
const candidateOnly = ui.windowAppearance(held, palette);
assert.equal(candidateOnly.kind, '観測あり・判定保留');
assert.ok(!candidateOnly.reference);
assert.ok(!candidateOnly.background.includes('#112233'));
const supportedBackground = {...held, background_call:held.background_candidate, background_status:'mixed_or_heterozygous'};
const reference = ui.windowAppearance(supportedBackground, palette);
assert.equal(reference.reference, true);
assert.ok(reference.kind.includes('参考背景'));
assert.ok(reference.background.includes('#ffffff'));
assert.notEqual(reference.background, candidateOnly.background);
const strict = ui.windowAppearance({...supportedBackground, call:'129S1_SvImJ'}, palette);
assert.equal(strict.background, '#112233');
assert.ok(!strict.reference);
assert.equal(ui.windowAppearance({...supportedBackground, n_sites:0}, palette).kind, '未観測');
assert.equal(ui.strictReason(held), held.reason);
assert.equal(ui.strictReason({...held,call:'129S1_SvImJ'}), '—');
const text = ui.sampleExplanation({chromosomes:[supportedBackground,held], diagnosis:{summary:{used_sites:26879,
  unresolved_substrain_pairs:[{strains:['C57BL_6J','C57BL_6NJ'],diagnostic_sites:7}]}}});
assert.ok(text.includes('26,879'));
assert.ok(text.includes('厳密判定 0染色体'));
assert.ok(text.includes('参考背景の基準を満たすもの 1染色体'));
assert.ok(text.includes('6J/6NJは識別不能'));
const detail = ui.windowDescription({...held,chromosome:'10',start:1,end:10000000});
assert.ok(detail.includes('厳密判定：保留'));
assert.ok(detail.includes('参考背景：保留'));
assert.ok(detail.includes(held.reason) && detail.includes(held.background_reason));
assert.doesNotThrow(() => ui.sampleExplanation({chromosomes:[held]}));
console.log('Background UI evidence separation passed.');
"""


class BackgroundResultAPITests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.job = {"id": "a" * 32, "request": {"srx": ["SRX1"]}}
        self.directory = self.root / self.job["id"] / "output" / "SRX1"
        self.directory.mkdir(parents=True)
        self.jobs = SimpleNamespace(root=self.root, sample_complete=lambda *_: True)
        self.row = {"chromosome": "10", "status": "ambiguous", "call": "", "n_sites": "1247",
                    "candidate": "129S1_SvImJ + C57BL_6NJ", "background_status": "ambiguous",
                    "background_call": "", "background_candidate": "129S1_SvImJ + C57BL/6 (6J/6NJ)",
                    "background_margin_per_site": "0.09668", "reason": "亜系統を区別できません。",
                    "background_reason": "群間のスコア差が基準未満です。"}
        for name in ("chromosomes.csv", "windows.csv", "candidate_scores.csv"):
            with (self.directory / name).open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=list(self.row))
                writer.writeheader()
                writer.writerow(self.row)

    def sample(self):
        return result_data(self.jobs, self.job)["samples"][0]

    def test_result_preserves_hold_and_background_evidence(self):
        diagnosis = {"schema": 1, "summary": {"used_sites": 26879},
                     "added_reference": "C57BL_6J", "reference_assumption": "Assumed REF, not measured"}
        (self.directory / "diagnosis.json").write_text(json.dumps(diagnosis), encoding="utf-8")
        sample = self.sample()
        self.assertEqual(sample["chromosomes"][0], self.row)
        self.assertEqual(sample["diagnosis"], diagnosis)
        self.assertIn("diagnosis.json", sample["files"])
        self.assertEqual(sample["chromosomes"][0]["call"], "")
        self.assertEqual(sample["chromosomes"][0]["background_call"], "")

    def test_legacy_results_have_no_optional_diagnosis(self):
        sample = self.sample()
        self.assertNotIn("diagnosis", sample)
        self.assertEqual(sample["chromosomes"][0]["status"], "ambiguous")

    def test_invalid_optional_diagnosis_does_not_hide_results(self):
        for value in ("{broken", "[]", '{"score":NaN}', '{"score":Infinity}', '{"score":1e999}'):
            with self.subTest(value=value):
                (self.directory / "diagnosis.json").write_text(value, encoding="utf-8")
                sample = self.sample()
                self.assertNotIn("diagnosis", sample)
                self.assertEqual(sample["chromosomes"][0], self.row)

    @unittest.skipUnless(shutil.which("node"), "Node is optional; JS_CHECKS can run in an available JS runtime")
    def test_reference_evidence_never_becomes_an_exact_call(self):
        result = subprocess.run([shutil.which("node"), "-e", JS_CHECKS, str(STATIC / "app.js")],
                                text=True, encoding="utf-8", capture_output=True, timeout=15)
        self.assertEqual(result.returncode, 0, result.stderr)


if __name__ == "__main__":
    unittest.main()
