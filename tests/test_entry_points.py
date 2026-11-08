from __future__ import annotations

import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path


class EntryPointTests(unittest.TestCase):
    def run_bad(self, command: str, arguments: list[str]) -> None:
        executable = shutil.which(command)
        if executable is None:
            self.skipTest(f"installed entry point is unavailable: {command}")
        result = subprocess.run([executable, *arguments], capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertNotIn("Traceback", result.stderr)
        error_lines = [line for line in result.stderr.splitlines() if line.strip()]
        self.assertEqual(len(error_lines), 1, result.stderr)

    def test_installed_entry_points_report_one_line_errors(self) -> None:
        with tempfile.TemporaryDirectory() as name:
            missing = str(Path(name) / "missing")
            false = shutil.which("false") or "/usr/bin/false"
            cases = {
                "build-cellulose-slab": [
                    "--layers", "1", "--width", "10", "--length", "10", "--no-restraints",
                    "--cellulose-forge-executable", false, "--output-directory", name,
                ],
                "parameterize-molecule": [
                    "--smiles", "CCO", "--output-directory", "/dev/null/output",
                    "--antechamber", false, "--parmchk2", false, "--tleap", false,
                ],
                "pack-liquid-cellulose": [
                    "--cellulose-coordinates", missing, "--water-model", "spc", "--molecules", "2",
                    "--density", "1", "--packmol", false, "--obabel", false,
                ],
                "combine-gromacs-topologies": [
                    "--slab", missing, "--liquid", missing, "--combined-coordinates", missing,
                ],
                "prepare-water-model": [
                    "--model", "spc", "--output-directory", "/dev/null/output",
                ],
            }
            for command, arguments in cases.items():
                with self.subTest(command=command):
                    self.run_bad(command, arguments)


if __name__ == "__main__":
    unittest.main()
