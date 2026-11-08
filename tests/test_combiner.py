from __future__ import annotations

import contextlib
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import combine_gromacs_topologies as combine


SLAB_TOP = """[ defaults ]
1 2 yes 1.0 1.0
[ atomtypes ]
CS 6 12.0 0.0 A 0.30 0.20
[ moleculetype ]
SLB 3
[ atoms ]
1 CS 1 CEL C1 1 0.0 12.0
[ system ]
slab
[ molecules ]
SLB 1
"""

LIQUID_TOP = """[ defaults ]
1 2 yes 0.5 0.83333333
[ atomtypes ]
A 6 12.0 0.0 A 0.30 0.20
B 1 1.0 0.0 A 0.40 0.80
[ moleculetype ]
MOLECULE_TYPE 3
[ atoms ]
1 A 1 LIQ A1 1 -0.1 12.0
2 A 1 LIQ A2 1  0.0 12.0
3 B 1 LIQ B1 1  0.0 1.0
4 B 1 LIQ B2 1  0.1 1.0
[ pairs ]
1 4 1
[ system ]
liquid
[ molecules ]
MOLECULE_TYPE 1
"""

COORDINATES = """CRYST1   20.000   20.000   30.000  90.00  90.00  90.00 P 1           1
ATOM      1   C1 CEL    1       1.000   1.000  10.000  1.00  0.00           C
HETATM    2   A1 LIQ    2       2.000   2.000  20.000  1.00  0.00           C
HETATM    3   A2 LIQ    2       2.500   2.000  20.000  1.00  0.00           C
HETATM    4   B1 LIQ    2       3.000   2.000  20.000  1.00  0.00           H
HETATM    5   B2 LIQ    2       3.500   2.000  20.000  1.00  0.00           H
END
"""


class CombinerTests(unittest.TestCase):
    def test_type1_to_type2_pair_values(self) -> None:
        block = combine.molecule_block(LIQUID_TOP, "MOLECULE_TYPE")
        converted = combine.explicit_secondary_pairs(
            block,
            combine.defaults(LIQUID_TOP),
            combine.atomtypes(LIQUID_TOP)[1],
        )
        row = next(line.split() for line in converted if line.split()[:3] == ["1", "4", "2"])
        self.assertAlmostEqual(float(row[3]), 0.83333333)
        self.assertAlmostEqual(float(row[4]), -0.1)
        self.assertAlmostEqual(float(row[5]), 0.1)
        self.assertAlmostEqual(float(row[6]), 0.35)
        self.assertAlmostEqual(float(row[7]), 0.2)

    def test_suffix_names_and_include_lines(self) -> None:
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            slab = root / "slab.top"
            liquid = root / "liquid.top"
            coordinates = root / "packed.pdb"
            slab.write_text(SLAB_TOP)
            liquid.write_text(LIQUID_TOP)
            coordinates.write_text(COORDINATES)
            arguments = [
                "combine-gromacs-topologies",
                "--slab", str(slab),
                "--liquid", str(liquid),
                "--combined-coordinates", str(coordinates),
                "--output-directory", str(root),
                "--suffix", "case",
            ]
            with patch.object(sys, "argv", arguments):
                self.assertEqual(combine.main(), 0)
            topology = (root / "system_case.top").read_text()
            self.assertIn('#include "forcefield_case.itp"', topology)
            self.assertIn('#include "slab_case.itp"', topology)
            self.assertIn('#include "liquid_case.itp"', topology)
            for filename in combine.output_names("case").values():
                self.assertTrue((root / filename).is_file())

    def test_residue_name_atom_count_and_order_checks(self) -> None:
        atoms = [
            combine.CoordinateAtom(1, "LIQ", name, (0.0, 0.0, 0.0))
            for name in ("A1", "A2", "B1", "B2")
        ]
        self.assertEqual(combine.validate_liquid_residues(atoms, "LIQ", ["A1", "A2", "B1", "B2"]), 1)
        with self.assertRaisesRegex(ValueError, "atom mismatch"):
            combine.validate_liquid_residues(atoms[:-1], "LIQ", ["A1", "A2", "B1", "B2"])
        with self.assertRaisesRegex(ValueError, "residue name"):
            combine.validate_liquid_residues(atoms, "MOL", ["A1", "A2", "B1", "B2"])

    def test_pairtypes_are_refused(self) -> None:
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            slab = root / "slab.top"
            liquid = root / "liquid.top"
            coordinates = root / "packed.pdb"
            slab.write_text(SLAB_TOP)
            liquid.write_text(LIQUID_TOP.replace("[ moleculetype ]", "[ pairtypes ]\nA B 1 0.3 0.2\n[ moleculetype ]"))
            coordinates.write_text(COORDINATES)
            arguments = [
                "combine-gromacs-topologies", "--slab", str(slab), "--liquid", str(liquid),
                "--combined-coordinates", str(coordinates), "--output-directory", str(root),
            ]
            with patch.object(sys, "argv", arguments):
                with self.assertRaisesRegex(ValueError, "pairtypes"):
                    combine.main()


if __name__ == "__main__":
    unittest.main()
