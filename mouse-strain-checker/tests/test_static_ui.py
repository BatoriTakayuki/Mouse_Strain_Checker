"""Browser entry-point accessibility and JavaScript input/display regressions."""
from html.parser import HTMLParser
from pathlib import Path
import shutil
import subprocess
import unittest


STATIC = Path(__file__).resolve().parents[1] / "mouse_strain" / "static"


class Document(HTMLParser):
    def __init__(self, text):
        super().__init__()
        self.tags = []
        self.feed(text)

    def handle_starttag(self, tag, attrs):
        self.tags.append((tag, dict(attrs)))


JS_CHECKS = r"""
const assert = require('node:assert/strict');
const ui = require(process.argv[1]);
assert.equal(ui.parseAccessions('srx12, ERX34\nDRX56;SRX12').join(','), 'SRX12,ERX34,DRX56');
assert.equal(ui.parseAccessions('SRX1、SRX2，SRX3').length, 3);
assert.equal(ui.parseAccessions('', false).length, 0);
for (const value of ['', 'SRR123', 'SRX-12', '<script>', 'SRX1/../']) assert.throws(() => ui.parseAccessions(value));
assert.equal(ui.parseAccessions(Array.from({length:30}, (_, i) => 'SRX' + i).join('\n')).length, 30);
assert.throws(() => ui.parseAccessions(Array.from({length:31}, (_, i) => 'SRX' + i).join('\n')));
assert.equal(ui.parseThreads('', 8), 8);
assert.equal(ui.parseThreads(' 128 '), 128);
for (const value of ['0', '129', '1.5', '-1', 'NaN', 'Infinity']) assert.throws(() => ui.parseThreads(value));
const palette = new Map([['A', '#112233'], ['B', '#445566']]);
const unobserved = ui.windowAppearance({n_sites:'0',call:'A'}, palette);
assert.equal(unobserved.kind, '未観測');
const held = ui.windowAppearance({n_sites:'20',candidate:'A',call:''}, palette);
assert.equal(held.kind, '観測あり・判定保留');
assert.notEqual(held.background, '#112233');
assert.notEqual(held.background, unobserved.background);
assert.equal(ui.windowAppearance({n_sites:'20',call:'A'}, palette).background, '#112233');
const mixed = ui.windowAppearance({n_sites:20,call:'A + B'}, palette);
assert.equal(mixed.kind, '2系統モデル');
assert.ok(mixed.background.includes('#112233') && mixed.background.includes('#445566'));
assert.equal(ui.modelParts('129S1/SvImJ').length, 1);
assert.equal(ui.selectedStillCurrent('old', 1, 'new', 2), false);
assert.equal(ui.selectedStillCurrent('same', 1, 'same', 2), false);
assert.equal(ui.selectedStillCurrent('same', 2, 'same', 2), true);
assert.equal(ui.formatBytes(1024), '1 KiB');
assert.equal(ui.formatBytes(null), '—');
console.log('Static UI validation and display rules passed.');
"""


class StaticUITests(unittest.TestCase):
    def test_accessible_self_contained_entry_point(self):
        document = Document((STATIC / "index.html").read_text(encoding="utf-8"))
        ids = [attrs["id"] for _, attrs in document.tags if "id" in attrs]
        self.assertEqual(len(ids), len(set(ids)), "HTML IDs must be unique")
        self.assertIn(("html", {"lang": "ja"}), document.tags)
        for tag, attrs in document.tags:
            if tag == "label" and "for" in attrs:
                self.assertIn(attrs["for"], ids)
            self.assertFalse(any(name.startswith("on") for name in attrs), "Inline event handlers are blocked by CSP")
            if tag == "script":
                self.assertEqual(attrs.get("src"), "/app.js")
            if tag == "link" and attrs.get("rel") == "stylesheet":
                self.assertEqual(attrs.get("href"), "/style.css")

    def test_api_strings_are_rendered_without_html_parsing(self):
        source = (STATIC / "app.js").read_text(encoding="utf-8")
        for unsafe in (".innerHTML", ".outerHTML", "insertAdjacentHTML", "document.write(", "eval("):
            self.assertNotIn(unsafe, source)

    @unittest.skipUnless(shutil.which("node"), "Node is optional; run JS_CHECKS in an available JS runtime")
    def test_javascript_input_and_evidence_display_rules(self):
        result = subprocess.run([shutil.which("node"), "-e", JS_CHECKS, str(STATIC / "app.js")],
                                text=True, encoding="utf-8", capture_output=True, timeout=15)
        self.assertEqual(result.returncode, 0, result.stderr)


if __name__ == "__main__":
    unittest.main()
