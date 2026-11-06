"""Rigid AMBER-family water-model geometries and GROMACS topologies."""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path

WATER_MOLAR_MASS = 18.01528


@dataclass(frozen=True)
class WaterModel:
    name: str
    atoms: tuple[tuple[str, str, float, float], ...]
    oh_nm: float
    hh_nm: float
    virtual_sites: str = ""


def _three_site_geometry(oh_angstrom: float, angle_degrees: float) -> list[tuple[str, float, float, float]]:
    angle = math.radians(angle_degrees)
    return [
        ("OW", 0.0, 0.0, 0.0),
        ("HW1", oh_angstrom, 0.0, 0.0),
        ("HW2", oh_angstrom * math.cos(angle), oh_angstrom * math.sin(angle), 0.0),
    ]


def geometry(model_name: str) -> list[tuple[str, float, float, float]]:
    model = get_water_model(model_name)
    angle = 109.47 if model.name in {"spc", "spce"} else 104.52
    sites = _three_site_geometry(model.oh_nm * 10, angle)
    if model.name in {"tip4p", "tip4pew"}:
        distance = 0.15 if model.name == "tip4p" else 0.125
        _, hx1, hy1, _ = sites[1]
        _, hx2, hy2, _ = sites[2]
        length = math.hypot(hx1 + hx2, hy1 + hy2)
        sites.append(("MW", distance * (hx1 + hx2) / length, distance * (hy1 + hy2) / length, 0.0))
    elif model.name == "tip5p":
        # Apply the GROMACS virtual-site construction in the generated H-O-H
        # frame. Coordinates here are angstrom, so c is converted from nm^-1
        # to angstrom^-1. This gives the required 0.070 nm O-LP distance.
        _, hx1, hy1, hz1 = sites[1]
        _, hx2, hy2, hz2 = sites[2]
        a = -0.344908262
        c = 0.64437903493
        cross = (
            hy1 * hz2 - hz1 * hy2,
            hz1 * hx2 - hx1 * hz2,
            hx1 * hy2 - hy1 * hx2,
        )
        base = (a * (hx1 + hx2), a * (hy1 + hy2), a * (hz1 + hz2))
        sites.extend(
            (
                ("LP1", *(base[index] - c * cross[index] for index in range(3))),
                ("LP2", *(base[index] + c * cross[index] for index in range(3))),
            )
        )
    return sites


MODELS = {
    "spc": WaterModel(
        "spc",
        (("OW", "OW_spc", -0.82, 15.9994), ("HW1", "HW_spc", 0.41, 1.008), ("HW2", "HW_spc", 0.41, 1.008)),
        0.10000,
        0.16330,
    ),
    "spce": WaterModel(
        "spce",
        (("OW", "OW_spc", -0.8476, 15.9994), ("HW1", "HW_spc", 0.4238, 1.008), ("HW2", "HW_spc", 0.4238, 1.008)),
        0.10000,
        0.16330,
    ),
    "tip3p": WaterModel(
        "tip3p",
        (("OW", "OW", -0.834, 16.0), ("HW1", "HW", 0.417, 1.008), ("HW2", "HW", 0.417, 1.008)),
        0.09572,
        0.15139,
    ),
    "tip4p": WaterModel(
        "tip4p",
        (("OW", "OW_tip4p", 0.0, 16.0), ("HW1", "HW_tip4p", 0.52, 1.008), ("HW2", "HW_tip4p", 0.52, 1.008), ("MW", "MW", -1.04, 0.0)),
        0.09572,
        0.15139,
        "[ virtual_sites3 ]\n; site from funct a b\n4 1 2 3 1 0.128012065 0.128012065\n",
    ),
    "tip4pew": WaterModel(
        "tip4pew",
        (("OW", "OW_tip4pew", 0.0, 16.0), ("HW1", "HW_tip4pew", 0.52422, 1.008), ("HW2", "HW_tip4pew", 0.52422, 1.008), ("MW", "MW", -1.04844, 0.0)),
        0.09572,
        0.15139,
        "[ virtual_sites3 ]\n; site from funct a b\n4 1 2 3 1 0.106676721 0.106676721\n",
    ),
    "tip5p": WaterModel(
        "tip5p",
        (("OW", "OW_tip5p", 0.0, 16.0), ("HW1", "HW_tip5p", 0.241, 1.008), ("HW2", "HW_tip5p", 0.241, 1.008), ("LP1", "MW", -0.241, 0.0), ("LP2", "MW", -0.241, 0.0)),
        0.09572,
        0.15139,
        "[ virtual_sites3 ]\n; site from funct a b c\n4 1 2 3 4 -0.344908262 -0.34490826 -6.4437903493\n5 1 2 3 4 -0.344908262 -0.34490826 6.4437903493\n",
    ),
}


ATOMTYPES = {
    "OW": (8, 16.0, "A", 0.315061, 0.636386),
    "HW": (1, 1.008, "A", 0.0, 0.0),
    "OW_spc": (8, 15.9994, "A", 0.316557, 0.650629),
    "HW_spc": (1, 1.008, "A", 0.0, 0.0),
    "OW_tip4p": (8, 16.0, "A", 0.315365, 0.648520),
    "HW_tip4p": (1, 1.008, "A", 0.0, 0.0),
    "OW_tip4pew": (8, 16.0, "A", 0.316435, 0.680946),
    "HW_tip4pew": (1, 1.008, "A", 0.0, 0.0),
    "OW_tip5p": (8, 16.0, "A", 0.312000, 0.669440),
    "HW_tip5p": (1, 1.008, "A", 0.0, 0.0),
    "MW": (0, 0.0, "D", 0.0, 0.0),
}


def get_water_model(name: str) -> WaterModel:
    key = name.lower().replace("-", "")
    aliases = {"spce": "spce", "tip4pew": "tip4pew"}
    key = aliases.get(key, key)
    if key not in MODELS:
        raise ValueError(f"unsupported water model {name!r}; choose: {', '.join(MODELS)}")
    return MODELS[key]


def write_water_pdb(model_name: str, output: Path, residue: str = "SOL") -> None:
    lines = [f"REMARK {model_name.upper()} rigid water template"]
    for index, (name, x, y, z) in enumerate(geometry(model_name), 1):
        element = "O" if name == "OW" else ("H" if name.startswith("H") else "X")
        lines.append(
            f"HETATM{index:5d} {name:>4.4s} {residue:>3.3s} {1:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {element:>2s}"
        )
    lines.extend(("TER", "END"))
    output.write_text("\n".join(lines) + "\n")


def topology_parts(model_name: str) -> tuple[str, str]:
    model = get_water_model(model_name)
    used_types = dict.fromkeys(atom_type for _, atom_type, _, _ in model.atoms)
    atomtype_lines = []
    for name in used_types:
        atomic_number, mass, ptype, sigma, epsilon = ATOMTYPES[name]
        atomtype_lines.append(
            f"{name:<14}{atomic_number:3d}{mass:12.6f}  0.000000  {ptype}  {sigma:.8g}  {epsilon:.8g}"
        )
    atoms = ["[ atoms ]", "; nr type resnr residue atom cgnr charge mass"]
    for index, (atom_name, atom_type, charge, mass) in enumerate(model.atoms, 1):
        atoms.append(
            f"{index:5d} {atom_type:<12} 1 SOL {atom_name:<4} 1 {charge: .8f} {mass:.6f}"
        )
    exclusions = ["[ exclusions ]"]
    indices = list(range(1, len(model.atoms) + 1))
    for index in indices:
        exclusions.append(" ".join(str(value) for value in (index, *[other for other in indices if other != index])))
    molecule = (
        "[ moleculetype ]\nSOL 2\n\n"
        + "\n".join(atoms)
        + f"\n\n[ settles ]\n1 1 {model.oh_nm:.5f} {model.hh_nm:.5f}\n\n"
        + model.virtual_sites
        + "\n"
        + "\n".join(exclusions)
        + "\n"
    )
    return "\n".join(atomtype_lines) + "\n", molecule


def write_water_topology(model_name: str, output: Path) -> None:
    atomtypes, molecule = topology_parts(model_name)
    output.write_text(
        "[ defaults ]\n; nbfunc comb-rule gen-pairs fudgeLJ fudgeQQ\n1 2 yes 1.0 1.0\n\n"
        "[ atomtypes ]\n"
        + atomtypes
        + "\n"
        + molecule
        + f"\n[ system ]\n{model_name.upper()} water\n\n[ molecules ]\nSOL 1\n"
    )
