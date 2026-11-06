from __future__ import annotations

import math
import unittest

from water_models import MODELS, geometry, topology_parts


class WaterModelTests(unittest.TestCase):
    def test_supported_models_have_neutral_charge(self) -> None:
        for name, model in MODELS.items():
            with self.subTest(model=name):
                self.assertAlmostEqual(sum(atom[2] for atom in model.atoms), 0.0, places=7)

    def test_geometry_site_count_matches_topology(self) -> None:
        for name, model in MODELS.items():
            with self.subTest(model=name):
                self.assertEqual(len(geometry(name)), len(model.atoms))

    def test_tip4p_virtual_site_distances(self) -> None:
        expected = {"tip4p": 0.15, "tip4pew": 0.125}
        for name, distance in expected.items():
            with self.subTest(model=name):
                _, x, y, z = geometry(name)[3]
                self.assertAlmostEqual(math.sqrt(x * x + y * y + z * z), distance, places=6)

    def test_tip5p_lone_pair_distances(self) -> None:
        for _, x, y, z in geometry("tip5p")[3:]:
            self.assertAlmostEqual(math.sqrt(x * x + y * y + z * z), 0.7, places=6)

    def test_topology_contains_virtual_sites_only_when_required(self) -> None:
        for name in MODELS:
            with self.subTest(model=name):
                _, molecule = topology_parts(name)
                self.assertEqual("[ virtual_sites3 ]" in molecule, name.startswith("tip4") or name == "tip5p")


if __name__ == "__main__":
    unittest.main()
