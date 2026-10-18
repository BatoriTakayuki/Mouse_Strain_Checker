"""CLI resource validation and setup completion regression tests."""
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import Mock, patch

from mouse_strain.__main__ import main, resources_paths
from mouse_strain.panel import sha256


class CliReviewTests(unittest.TestCase):
    def test_invalid_resource_configs_are_explained(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            for value in (None, [], {}, {"fasta": 3, "panel": "panel", "index": "index"}):
                with self.subTest(value=value):
                    (root / "config.json").write_text(json.dumps(value), encoding="utf-8")
                    with self.assertRaisesRegex(ValueError, "Invalid resources config"):
                        resources_paths(root)
                    self.assertEqual(main(["run", "--resources", str(root), "--srx", "SRX1",
                                           "--out", str(root / "out")]), 1)
                    self.assertFalse((root / "out").exists())

    def test_setup_does_not_publish_config_for_incomplete_index(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fasta = root / "reference.fa"
            fasta.write_bytes(b">1\nAAAA\n")
            panel = root / "panel"
            panel.mkdir()
            (panel / "metadata.json").write_text("{}", encoding="utf-8")
            meta = dict(assembly="GRCm39", fasta_sha256=sha256(fasta))
            (root / "hisat2").mkdir()
            (root / "hisat2/genome.reference.json").write_text(json.dumps(meta), encoding="utf-8")
            database = Mock()
            with patch("mouse_strain.__main__.require_tools"), \
                 patch("mouse_strain.__main__.download"), \
                 patch("mouse_strain.__main__.load_panel", return_value=(meta, database)):
                self.assertEqual(main(["setup", "--out", str(root)]), 1)
            database.close.assert_called_once()
            self.assertFalse((root / "config.json").exists())


if __name__ == "__main__":
    unittest.main()
