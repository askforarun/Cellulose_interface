"""AmberTools integration checks; skipped in the lightweight CI environment."""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from combine_gromacs_topologies import (
    CoordinateAtom,
    molecule_block,
    molecule_coordinate_signature,
    validate_liquid_residues,
)


ROOT = Path(__file__).resolve().parents[1]


def mol2_atom_rows(path: Path) -> list[list[str]]:
    lines = path.read_text().splitlines()
    start = lines.index("@<TRIPOS>ATOM") + 1
    end = lines.index("@<TRIPOS>BOND", start)
    return [line.split() for line in lines[start:end] if line.strip()]


class ParameterizationIntegrationTests(unittest.TestCase):
    @unittest.skipUnless(
        all(shutil.which(command) for command in ("antechamber", "parmchk2", "tleap")),
        "AmberTools executables are unavailable",
    )
    def test_multiple_bonds_keep_gaff2_types_and_template_order(self) -> None:
        try:
            import parmed  # noqa: F401
            import rdkit  # noqa: F401
        except ImportError:
            self.skipTest("ParmEd and RDKit are required")
        cases = (
            ("acetic_acid", "CC(=O)O", {"C2": "c", "O1": "o"}),
            ("acetone", "CC(=O)C", {"C2": "c", "O1": "o"}),
            ("benzene", "c1ccccc1", {f"C{i}": "ca" for i in range(1, 7)}),
            ("ethyl_acetate", "CCOC(=O)C", {"C3": "c", "O2": "o"}),
        )
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for prefix, smiles, expected_types in cases:
                with self.subTest(molecule=prefix):
                    result = subprocess.run(
                        [
                            sys.executable, str(ROOT / "parameterize_molecule.py"),
                            "--smiles", smiles, "--residue-name", "MOL",
                            "--prefix", prefix, "--output-directory", directory,
                        ],
                        cwd=ROOT,
                        env=os.environ.copy(),
                        capture_output=True,
                        text=True,
                    )
                    self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                    atom_types = {row[1]: row[5] for row in mol2_atom_rows(root / f"{prefix}.mol2")}
                    for atom_name, expected_type in expected_types.items():
                        self.assertEqual(atom_types[atom_name], expected_type)
                    topology = (root / f"{prefix}.top").read_text()
                    residue, expected_names = molecule_coordinate_signature(molecule_block(topology))
                    atoms = [
                        CoordinateAtom(int(line[22:26]), line[17:20].strip(), line[12:16].strip(),
                                       (float(line[30:38]), float(line[38:46]), float(line[46:54])))
                        for line in (root / f"{prefix}.pdb").read_text().splitlines()
                        if line.startswith(("ATOM  ", "HETATM"))
                    ]
                    self.assertEqual(validate_liquid_residues(atoms, residue, expected_names), 1)


if __name__ == "__main__":
    unittest.main()
