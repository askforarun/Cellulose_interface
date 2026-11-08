from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from build_cellulose_slab import NO_RESTRAINT_MARKER, collect_outputs, mark_unrestrained


class BuilderTests(unittest.TestCase):
    def test_collect_outputs_renames_files_and_report_entries(self) -> None:
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            staging = root / "staging"
            output = root / "output"
            staging.mkdir()
            output.mkdir()
            files = {
                "coordinates": "structure.pdb",
                "gromacs_coordinates": "structure.gro",
                "gromacs_topology": "system.top",
            }
            for filename in files.values():
                (staging / filename).write_text(filename + "\n")
            (staging / "build_report.json").write_text(json.dumps({"outputs": files}))
            written = collect_outputs(staging, output)
            self.assertEqual(
                {path.name for path in written},
                {"cellulose.pdb", "cellulose.gro", "cellulose.top", "cellulose_build_report.json"},
            )
            report = json.loads((output / "cellulose_build_report.json").read_text())
            self.assertEqual(report["outputs"]["gromacs_topology"], "cellulose.top")

    def test_no_restraints_marker_is_inserted_before_system(self) -> None:
        with tempfile.TemporaryDirectory() as name:
            topology = Path(name) / "cellulose.top"
            topology.write_text("[ moleculetype ]\nCELL 3\n\n[ system ]\ncellulose\n")
            mark_unrestrained(topology)
            text = topology.read_text()
            self.assertIn(NO_RESTRAINT_MARKER, text)
            self.assertLess(text.index("POSRES_CELLULOSE"), text.index("[ system ]"))


if __name__ == "__main__":
    unittest.main()
