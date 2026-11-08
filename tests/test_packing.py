from __future__ import annotations

import contextlib
import io
import json
import shutil
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pack_liquid_cellulose


SLAB_PDB = """CRYST1   20.000   20.000   10.000  90.00  90.00  90.00 P 1           1
ATOM      1   C1 CEL    1       1.000   1.000   2.000  1.00  0.00           C
ATOM      2   C2 CEL    1       1.000   1.000   8.000  1.00  0.00           C
END
"""


class PackingTests(unittest.TestCase):
    def arguments(self, root: Path, packmol: Path, *extra: str) -> list[str]:
        slab = root / "slab.pdb"
        slab.write_text(SLAB_PDB)
        return [
            "pack-liquid-cellulose",
            "--cellulose-coordinates",
            str(slab),
            "--water-model",
            "spc",
            "--molecules",
            "2",
            "--density",
            "1.0",
            "--output-directory",
            str(root),
            "--packmol",
            str(packmol),
            "--obabel",
            str(packmol),
            *extra,
        ]

    def test_dry_run_reports_liquid_only_growth_and_density_refusal(self) -> None:
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            output = io.StringIO()
            arguments = self.arguments(
                root,
                Path(shutil.which("true") or "/usr/bin/true"),
                "--dry-run",
                "--max-attempts",
                "3",
                "--z-growth-factor",
                "1.1",
                "--max-density-loss",
                "0.05",
            )
            with patch.object(sys, "argv", arguments), contextlib.redirect_stdout(output):
                self.assertEqual(pack_liquid_cellulose.main(), 0)
            text = output.getvalue()
            self.assertIn("Attempt 1:", text)
            self.assertIn("density=1.000000 g/cm^3 (allowed)", text)
            self.assertIn("Attempt 2:", text)
            self.assertIn("density=0.909091 g/cm^3 (refused)", text)

    def test_existing_output_survives_packmol_failure(self) -> None:
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            packmol = root / "packmol-fails"
            packmol.write_text("#!/bin/sh\nexit 1\n")
            packmol.chmod(0o755)
            existing = root / "interface.pdb"
            existing.write_text("previous valid coordinates\n")
            arguments = self.arguments(root, packmol, "--max-attempts", "1")
            with patch.object(sys, "argv", arguments):
                with self.assertRaisesRegex(RuntimeError, "Packmol failed"):
                    pack_liquid_cellulose.main()
            self.assertEqual(existing.read_text(), "previous valid coordinates\n")

    def test_density_bound_refuses_retry_and_preserves_output(self) -> None:
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            packmol = root / "packmol-fails"
            packmol.write_text("#!/bin/sh\nexit 1\n")
            packmol.chmod(0o755)
            existing = root / "interface.pdb"
            existing.write_text("keep me\n")
            arguments = self.arguments(
                root,
                packmol,
                "--max-attempts",
                "3",
                "--z-growth-factor",
                "1.05",
                "--max-density-loss",
                "0",
            )
            with patch.object(sys, "argv", arguments):
                with self.assertRaisesRegex(RuntimeError, "density-loss bound") as caught:
                    pack_liquid_cellulose.main()
            self.assertIn("attempt=2", str(caught.exception))
            self.assertIn("achieved_density=0.952381", str(caught.exception))
            self.assertEqual(existing.read_text(), "keep me\n")

    def test_packmol_version_capture_cannot_block_and_report_records_version(self) -> None:
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            packmol = root / "packmol-fake"
            packmol.write_text(
                "#!" + sys.executable + "\n"
                "import pathlib, re, sys\n"
                "data = sys.stdin.read()\n"
                "print('######## PACKMOL ########')\n"
                "print('Version 21.0.1')\n"
                "if not data: sys.exit(0)\n"
                "output = pathlib.Path(re.search(r'^output (.+)$', data, re.M).group(1))\n"
                "paths = re.findall(r'^structure (.+)$', data, re.M)\n"
                "atoms = []\n"
                "for resid, path in enumerate(paths, 1):\n"
                "    for line in pathlib.Path(path).read_text().splitlines():\n"
                "        if line.startswith(('ATOM  ', 'HETATM')):\n"
                "            atoms.append(line[:6] + f'{len(atoms)+1:5d}' + line[11:22] + f'{resid:4d}' + line[26:])\n"
                "output.write_text('\\n'.join(atoms + ['END', '']))\n"
                "print('Success!')\n"
            )
            packmol.chmod(0o755)
            arguments = self.arguments(root, packmol, "--max-attempts", "1")
            with patch.object(sys, "argv", arguments), contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(pack_liquid_cellulose.main(), 0)
            report = json.loads((root / "interface_packing_report.json").read_text())
            self.assertEqual(report["tool_versions"]["packmol"], "21.0.1")
            self.assertTrue((root / "interface.pdb").is_file())


if __name__ == "__main__":
    unittest.main()
