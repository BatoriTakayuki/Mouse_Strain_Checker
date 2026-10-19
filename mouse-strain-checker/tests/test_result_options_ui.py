"""Unresolved strain options and provisional RNA sex stay visible and separate."""
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
const options = '(129S1_SvImJ + C57BL_6J) or (129S1_SvImJ + C57BL_6NJ)';
const held = {chromosome:'2',n_sites:1247,status:'ambiguous',call:'',candidate_options:options,
  candidate_count:2,candidate_options_reason:'6Jと6NJを識別する根拠が不足しています。'};
assert.equal(ui.resultCall(held), '候補：' + options + '（未確定）');
assert.equal(ui.resultCall({...held,call:'129S1_SvImJ'}), '129S1_SvImJ');
assert.equal(ui.resultCall({status:'no_observed_snps',call:'',candidate:'A'}), '保留');
assert.ok(ui.resultCall({...held,candidate_options:'A',candidate_count:1}).includes('未確定'));
assert.ok(ui.sampleExplanation({chromosomes:[held]}).includes('未確定の 1染色体'));
const summary = ui.candidateSummary({chromosomes:[held,{...held,chromosome:'3'},
  {...held,chromosome:'4',call:'A'},{chromosome:'Y',call:'',candidate:'A'}]});
assert.equal(summary, '染色体 2, 3：' + options + '（未確定）');
const detail = ui.windowDescription({...held,start:1,end:100});
assert.ok(detail.includes(options) && detail.includes(held.candidate_options_reason));
assert.equal(ui.sexDescription({}).label, '未評価（旧結果）');
const male = ui.sexDescription({sex:{status:'male_likely',label:'オス候補（RNA推定）',reason:'Yマーカーを観測しました。',
  evidence:{y_fragments:1234,xist_fragments:0},markers:[{gene:'Ddx3y',fragments:1234}]}});
assert.equal(male.label, 'オス候補（RNA推定）');
assert.equal(male.evidence[0][1], '1,234');
assert.equal(male.markers[0].fragments, 1234);
assert.ok(ui.sexDescription({sex:{status:'female_likely'}}).label.includes('推定'));
assert.equal(ui.sexDescription({sex:{status:'inconclusive',reason:'根拠不足'}}).label, '判定保留');
assert.equal(ui.sexDescription({sex:[]}).label, '未評価（旧結果）');
assert.doesNotThrow(() => ui.sexDescription({sex:{evidence:{extra:{a:2},nil:null}}}));
assert.ok(ui.sampleExplanation({workflow:{stages:['prefetch','fasterq-dump','fastp','HISAT2']}}).includes('prefetch → fasterq-dump → fastp → HISAT2'));
assert.ok(ui.sampleExplanation({workflow:{quality_notes:['SRA Liteの品質値は簡略化されています。']}}).includes('SRA Liteの品質値は簡略化されています。'));
console.log('Unresolved options and RNA sex display passed.');
"""


class ResultOptionsAPITests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.job = {"id": "b" * 32, "request": {"srx": ["SRX1"]}}
        self.directory = self.root / self.job["id"] / "output" / "SRX1"
        self.directory.mkdir(parents=True)
        self.jobs = SimpleNamespace(root=self.root, sample_complete=lambda *_: True)
        self.row = {"chromosome": "2", "status": "ambiguous", "call": "", "candidate": "A",
                    "candidate_options": "A or B", "candidate_count": "2",
                    "candidate_options_reason": "候補間の識別根拠が不足しています。"}
        for name in ("chromosomes.csv", "windows.csv", "candidate_scores.csv"):
            with (self.directory / name).open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=list(self.row))
                writer.writeheader()
                writer.writerow(self.row)

    def sample(self):
        return result_data(self.jobs, self.job)["samples"][0]

    def test_options_and_rna_sex_are_preserved_without_turning_hold_into_call(self):
        sex = {"schema": 1, "status": "male_likely", "label": "オス候補（RNA推定）",
               "reason": "Yマーカーを観測しました。", "evidence": {"y_fragments": 60},
               "markers": [{"gene": "Ddx3y", "fragments": 60}]}
        (self.directory / "sex.json").write_text(json.dumps(sex), encoding="utf-8")
        sample = self.sample()
        self.assertEqual(sample["sex"], sex)
        self.assertIn("sex.json", sample["files"])
        self.assertEqual(sample["chromosomes"][0], self.row)
        self.assertEqual(sample["windows"][0]["call"], "")

    def test_legacy_result_does_not_invent_a_sex(self):
        sample = self.sample()
        self.assertNotIn("sex", sample)
        self.assertNotIn("workflow", sample)
        self.assertEqual(sample["chromosomes"][0]["candidate_options"], "A or B")

    def test_workflow_provenance_is_available_without_claiming_reanalysis(self):
        workflow = {"schema": 1, "stages": ["prefetch", "fasterq-dump", "fastp", "HISAT2"], "fastp_enabled": True}
        (self.directory / "workflow.json").write_text(json.dumps(workflow), encoding="utf-8")
        sample = self.sample()
        self.assertEqual(sample["workflow"], workflow)
        self.assertIn("workflow.json", sample["files"])
        self.assertEqual(sample["chromosomes"][0]["call"], "")

    def test_invalid_or_oversized_optional_sex_does_not_hide_strain_result(self):
        for value in ("{broken", "[]", '{"count":NaN}', '{"count":1e999}',
                      json.dumps({"reason": "a" * (1024 * 1024)})):
            with self.subTest(prefix=value[:30]):
                (self.directory / "sex.json").write_text(value, encoding="utf-8")
                sample = self.sample()
                self.assertNotIn("sex", sample)
                self.assertEqual(sample["chromosomes"][0], self.row)

    @unittest.skipUnless(shutil.which("node"), "Node is optional; JS_CHECKS can run in an available JS runtime")
    def test_options_and_sex_display_rules(self):
        result = subprocess.run([shutil.which("node"), "-e", JS_CHECKS, str(STATIC / "app.js")],
                                text=True, encoding="utf-8", capture_output=True, timeout=15)
        self.assertEqual(result.returncode, 0, result.stderr)


if __name__ == "__main__":
    unittest.main()
