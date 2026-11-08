#!/usr/bin/env python3
"""Combine GROMACS topologies with different 1-4 scaling conventions."""

from __future__ import annotations

import argparse
import math
import re
import sys
from dataclasses import dataclass
from pathlib import Path

from project_version import project_version

VERSION = project_version()
SECTION = re.compile(r"^\s*\[\s*([^]]+?)\s*]\s*(?:;.*)?$")


@dataclass(frozen=True)
class Defaults:
    nbfunc: int
    combination_rule: int
    generate_pairs: str
    fudge_lj: float
    fudge_qq: float


@dataclass
class CoordinateAtom:
    resid: int
    resname: str
    atomname: str
    xyz: tuple[float, float, float]


def read_coordinates(path: Path) -> tuple[list[CoordinateAtom], list[float]]:
    lines = path.read_text().splitlines()
    if path.suffix.lower() == ".gro":
        count = int(lines[1])
        atoms = [
            CoordinateAtom(
                int(line[:5]), line[5:10].strip(), line[10:15].strip(),
                (float(line[20:28]), float(line[28:36]), float(line[36:44])),
            )
            for line in lines[2 : 2 + count]
        ]
        fields = lines[2 + count].split()
        if len(fields) != 3:
            raise ValueError("only orthorhombic GRO boxes are supported")
        return atoms, [float(value) for value in fields]
    if path.suffix.lower() == ".pdb":
        atoms: list[CoordinateAtom] = []
        box = None
        for line in lines:
            if line.startswith("CRYST1"):
                box = [float(line[6:15]) / 10, float(line[15:24]) / 10, float(line[24:33]) / 10]
            elif line.startswith(("ATOM  ", "HETATM")):
                atoms.append(
                    CoordinateAtom(
                        int(line[22:26]), line[17:20].strip(), line[12:16].strip(),
                        (float(line[30:38]) / 10, float(line[38:46]) / 10, float(line[46:54]) / 10),
                    )
                )
        if box is None:
            raise ValueError(f"PDB coordinate file has no CRYST1 record: {path}")
        return atoms, box
    raise ValueError("coordinate files must use .gro or .pdb")


def residue_count(atoms: list[CoordinateAtom]) -> int:
    return sum(
        index == 0
        or (atom.resid, atom.resname) != (atoms[index - 1].resid, atoms[index - 1].resname)
        for index, atom in enumerate(atoms)
    )


def write_gro(atoms: list[CoordinateAtom], box: list[float], output: Path) -> None:
    lines = ["Combined coordinates written by combine-gromacs-topologies", f"{len(atoms):5d}"]
    for serial, atom in enumerate(atoms, 1):
        lines.append(
            f"{atom.resid % 100000:5d}{atom.resname:<5.5s}{atom.atomname:>5.5s}"
            f"{serial % 100000:5d}{atom.xyz[0]:8.3f}{atom.xyz[1]:8.3f}{atom.xyz[2]:8.3f}"
        )
    lines.append("".join(f"{value:10.5f}" for value in box))
    output.write_text("\n".join(lines) + "\n")


def use_combined_coordinates(
    path: Path,
    output: Path,
    liquid_resname: str,
    expected_atom_names: list[str],
) -> tuple[int, int, list[float]]:
    atoms, box = read_coordinates(path)
    liquid_indices = [
        index for index, atom in enumerate(atoms) if atom.resname == liquid_resname
    ]
    if not liquid_indices:
        raise ValueError(
            f"combined coordinates contain no liquid residue named {liquid_resname}"
        )
    first_liquid = min(liquid_indices)
    liquid_atoms = atoms[first_liquid:]
    if any(atom.resname != liquid_resname for atom in liquid_atoms):
        raise ValueError(
            "combined coordinates must list all slab atoms before liquid molecules"
        )
    liquid_count = validate_liquid_residues(
        liquid_atoms, liquid_resname, expected_atom_names
    )
    write_gro(atoms, box, output)
    return first_liquid, liquid_count, box


def validate_liquid_residues(
    atoms: list[CoordinateAtom], residue_name: str, expected_atom_names: list[str]
) -> int:
    if not expected_atom_names:
        raise ValueError("liquid topology contains no atoms")
    count = 0
    start = 0
    while start < len(atoms):
        atom = atoms[start]
        key = (atom.resid, atom.resname)
        end = start + 1
        while end < len(atoms) and (atoms[end].resid, atoms[end].resname) == key:
            end += 1
        if atom.resname != residue_name:
            raise ValueError(
                f"expected liquid residue name {residue_name}, found {atom.resname}"
            )
        found_names = [entry.atomname for entry in atoms[start:end]]
        if found_names != expected_atom_names:
            raise ValueError(
                f"liquid residue {atom.resid} ({residue_name}) atom mismatch: "
                f"expected {expected_atom_names}, found {found_names}"
            )
        count += 1
        start = end
    return count


def assemble_coordinates(
    primary_path: Path,
    secondary_path: Path,
    output_path: Path,
    axis_name: str,
    gap: float,
    layout: str,
) -> tuple[int, int, list[float]]:
    primary_atoms, primary_box = read_coordinates(primary_path)
    secondary_atoms, secondary_box = read_coordinates(secondary_path)
    axis = "xyz".index(axis_name)
    lateral = [index for index in range(3) if index != axis]
    if any(abs(primary_box[index] - secondary_box[index]) > 1e-4 for index in lateral):
        raise ValueError("slab and liquid boxes must have matching lateral dimensions")
    # Preserve whole molecules. Atom-wise wrapping can split bonded molecules across
    # a periodic boundary and creates very long excluded/bonded distances in grompp.
    primary_xyz = [list(atom.xyz) for atom in primary_atoms]
    secondary_xyz = [list(atom.xyz) for atom in secondary_atoms]
    primary_min = min(xyz[axis] for xyz in primary_xyz)
    primary_max = max(xyz[axis] for xyz in primary_xyz)
    secondary_min = min(xyz[axis] for xyz in secondary_xyz)
    secondary_max = max(xyz[axis] for xyz in secondary_xyz)
    primary_extent = primary_max - primary_min
    secondary_extent = secondary_max - secondary_min
    box = list(primary_box)
    box[axis] = primary_extent + secondary_extent + 2 * gap

    if layout == "stacked":
        for xyz in primary_xyz:
            xyz[axis] -= primary_min
        shift = primary_extent + gap - secondary_min
        for xyz in secondary_xyz:
            xyz[axis] += shift
    elif layout == "centered-slab":
        # Put the slab in the middle of the box. Place the liquid after its upper
        # surface and wrap whole residues through the periodic boundary. This
        # writes an explicit liquid | slab | liquid view without splitting any
        # bonded molecule atom-by-atom.
        primary_start = secondary_extent / 2 + gap
        for xyz in primary_xyz:
            xyz[axis] += primary_start - primary_min

        liquid_shift = primary_start + primary_extent + gap - secondary_min
        start = 0
        while start < len(secondary_atoms):
            key = (secondary_atoms[start].resid, secondary_atoms[start].resname)
            end = start + 1
            while end < len(secondary_atoms) and (
                secondary_atoms[end].resid, secondary_atoms[end].resname
            ) == key:
                end += 1
            molecule = secondary_xyz[start:end]
            center = sum(xyz[axis] + liquid_shift for xyz in molecule) / len(molecule)
            periodic_shift = liquid_shift
            while center >= box[axis]:
                center -= box[axis]
                periodic_shift -= box[axis]
            while center < 0:
                center += box[axis]
                periodic_shift += box[axis]
            for xyz in molecule:
                xyz[axis] += periodic_shift
            start = end
    else:
        raise ValueError(f"unsupported coordinate layout: {layout}")

    lines = ["Interface assembled by combine-gromacs-topologies", f"{len(primary_atoms) + len(secondary_atoms):5d}"]
    maximum_primary_residue = max(atom.resid for atom in primary_atoms)
    all_atoms = list(zip(primary_atoms, primary_xyz)) + [
        (
            CoordinateAtom(
                maximum_primary_residue + atom.resid, atom.resname, atom.atomname, atom.xyz
            ),
            xyz,
        )
        for atom, xyz in zip(secondary_atoms, secondary_xyz)
    ]
    for serial, (atom, xyz) in enumerate(all_atoms, 1):
        lines.append(
            f"{atom.resid % 100000:5d}{atom.resname:<5.5s}{atom.atomname:>5.5s}"
            f"{serial % 100000:5d}{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}"
        )
    lines.append("".join(f"{value:10.5f}" for value in box))
    output_path.write_text("\n".join(lines) + "\n")
    return len(primary_atoms), residue_count(secondary_atoms), box


def sections(text: str) -> list[tuple[str, int, int]]:
    lines = text.splitlines()
    starts: list[tuple[str, int]] = []
    for index, line in enumerate(lines):
        match = SECTION.match(line)
        if match:
            starts.append((match.group(1).strip().lower(), index))
    return [
        (name, start, starts[position + 1][1] if position + 1 < len(starts) else len(lines))
        for position, (name, start) in enumerate(starts)
    ]


def section(text: str, name: str) -> tuple[int, int, list[str]]:
    lines = text.splitlines()
    matches = [(start, end) for found, start, end in sections(text) if found == name]
    if len(matches) != 1:
        raise ValueError(f"expected exactly one [{name}] section, found {len(matches)}")
    start, end = matches[0]
    return start, end, lines[start:end]


def data_rows(lines: list[str]) -> list[list[str]]:
    rows: list[list[str]] = []
    for line in lines[1:]:
        value = line.split(";", 1)[0].strip()
        if value and not value.startswith("#"):
            rows.append(value.split())
    return rows


def defaults(text: str) -> Defaults:
    _, _, lines = section(text, "defaults")
    rows = data_rows(lines)
    if len(rows) != 1 or len(rows[0]) < 5:
        raise ValueError("[defaults] must contain one row with at least five fields")
    row = rows[0]
    return Defaults(int(row[0]), int(row[1]), row[2], float(row[3]), float(row[4]))


def atomtypes(text: str) -> tuple[list[str], dict[str, tuple[float, float]]]:
    _, _, lines = section(text, "atomtypes")
    parsed: dict[str, tuple[float, float]] = {}
    output: list[str] = []
    for line in lines[1:]:
        fields = line.split(";", 1)[0].split()
        if not fields:
            continue
        if fields[0].startswith("#"):
            continue
        if len(fields) < 7:
            raise ValueError(f"unsupported [atomtypes] row: {line}")
        parsed[fields[0]] = (float(fields[-2]), float(fields[-1]))
        output.append(line)
    return output, parsed


def molecule_block(text: str, desired_name: str | None = None) -> list[str]:
    lines = text.splitlines()
    starts = [start for name, start, _ in sections(text) if name == "moleculetype"]
    system_starts = [start for name, start, _ in sections(text) if name == "system"]
    if not starts or len(system_starts) != 1:
        raise ValueError("each input must contain molecule and system definitions")
    blocks = [
        lines[start : starts[index + 1] if index + 1 < len(starts) else system_starts[0]]
        for index, start in enumerate(starts)
    ]
    if desired_name is None:
        if len(blocks) != 1:
            raise ValueError("topology defines multiple molecules but no selection was provided")
        return blocks[0]
    for block in blocks:
        if molecule_name(block) == desired_name:
            return block
    raise ValueError(f"molecule type referenced by [molecules] was not found: {desired_name}")


def molecule_name(block: list[str]) -> str:
    for fields in data_rows(block):
        return fields[0]
    raise ValueError("molecule name is missing")


def atom_data(block: list[str]) -> tuple[dict[int, str], dict[int, float]]:
    text = "\n".join(block)
    _, _, lines = section(text, "atoms")
    types: dict[int, str] = {}
    charges: dict[int, float] = {}
    for row in data_rows(lines):
        if len(row) < 8:
            raise ValueError(f"unsupported [atoms] row: {' '.join(row)}")
        index = int(row[0])
        types[index] = row[1]
        charges[index] = float(row[6])
    return types, charges


def molecule_coordinate_signature(block: list[str]) -> tuple[str, list[str]]:
    """Return the single residue name and ordered atom names for one molecule."""
    _, _, lines = section("\n".join(block), "atoms")
    rows = data_rows(lines)
    if not rows:
        raise ValueError("liquid topology contains no [ atoms ] rows")
    if any(len(row) < 5 for row in rows):
        raise ValueError("unsupported liquid [ atoms ] row")
    residue_names = {row[3] for row in rows}
    if len(residue_names) != 1:
        raise ValueError(
            "liquid topology must use one residue name to infer packed molecules; "
            "pass a one-residue liquid topology"
        )
    return next(iter(residue_names)), [row[4] for row in rows]


def correct_charge_rounding(block: list[str], tolerance: float = 1e-4) -> list[str]:
    """Move a tiny printed-charge residual onto the final atom."""
    output = block.copy()
    text = "\n".join(block)
    start, end, lines = section(text, "atoms")
    atom_rows: list[tuple[int, list[str], str]] = []
    for relative_index, line in enumerate(lines[1:], 1):
        value, separator, comment = line.partition(";")
        fields = value.split()
        if fields:
            atom_rows.append((start + relative_index, fields, separator + comment))
    charge = sum(float(fields[6]) for _, fields, _ in atom_rows)
    target = round(charge)
    residual = target - charge
    if abs(residual) > tolerance:
        raise ValueError(
            f"molecule charge {charge:.8f} is not within {tolerance:g} of an integer"
        )
    if residual and atom_rows:
        line_index, fields, comment = atom_rows[-1]
        fields[6] = f"{float(fields[6]) + residual:.8f}"
        output[line_index] = "    " + " ".join(fields) + (f" {comment}" if comment else "")
    return output


def explicit_secondary_pairs(
    block: list[str],
    secondary_defaults: Defaults,
    type_parameters: dict[str, tuple[float, float]],
) -> list[str]:
    if secondary_defaults.combination_rule != 2:
        raise ValueError("explicit pair conversion currently supports combination rule 2 only")
    atom_types, charges = atom_data(block)
    output = block.copy()
    text = "\n".join(block)
    pair_sections = [(start, end) for name, start, end in sections(text) if name == "pairs"]
    if not pair_sections:
        return block
    if len(pair_sections) != 1:
        raise ValueError(f"expected at most one liquid [pairs] section, found {len(pair_sections)}")
    start, end = pair_sections[0]
    replacement = [block[start], "; ai aj funct fudgeQQ qi qj sigma14 epsilon14"]
    for line in block[start + 1 : end]:
        value = line.split(";", 1)[0].strip()
        if not value:
            continue
        fields = value.split()
        if len(fields) != 3 or fields[2] != "1":
            raise ValueError(f"unsupported liquid pair row: {line}")
        first, second = int(fields[0]), int(fields[1])
        sigma_i, epsilon_i = type_parameters[atom_types[first]]
        sigma_j, epsilon_j = type_parameters[atom_types[second]]
        sigma_14 = (sigma_i + sigma_j) / 2
        epsilon_14 = math.sqrt(epsilon_i * epsilon_j) * secondary_defaults.fudge_lj
        replacement.append(
            f"{first:7d}{second:7d}  2  {secondary_defaults.fudge_qq:.8f}"
            f"  {charges[first]:.8f}  {charges[second]:.8f}"
            f"  {sigma_14:.10g}  {epsilon_14:.10g}"
        )
    return output[:start] + replacement + output[end:]


def selected_molecule(text: str) -> tuple[str, int]:
    rows = molecule_composition(text)
    if len(rows) != 1:
        raise ValueError("topology must contain exactly one [molecules] row here")
    return rows[0]


def molecule_composition(text: str) -> list[tuple[str, int]]:
    _, _, lines = section(text, "molecules")
    rows = data_rows(lines)
    if not rows:
        raise ValueError("topology contains no [molecules] rows")
    return [(row[0], int(row[1])) for row in rows]


def suffix(value: str) -> str:
    if value and not all(character.isalnum() or character in "-_" for character in value):
        raise argparse.ArgumentTypeError("suffix may contain only letters, digits, '-' and '_'")
    return value


def output_names(case_suffix: str) -> dict[str, str]:
    """File names written in --output-directory; a suffix is appended to every stem."""
    tail = f"_{case_suffix}" if case_suffix else ""
    return {
        "coordinates": f"interface{tail}.gro",
        "topology": f"system{tail}.top",
        "forcefield": f"forcefield{tail}.itp",
        "slab": f"slab{tail}.itp",
        "liquid": f"liquid{tail}.itp",
    }


def main() -> int:
    parser = argparse.ArgumentParser(prog="combine-gromacs-topologies")
    parser.add_argument("--version", action="version", version=f"%(prog)s {VERSION}")
    parser.add_argument("--slab", required=True, type=Path)
    parser.add_argument("--liquid", required=True, type=Path)
    parser.add_argument(
        "--liquid-residue-name",
        help="packed-coordinate residue name; defaults to the liquid [ atoms ] residue",
    )
    parser.add_argument("--slab-coordinates", type=Path)
    parser.add_argument("--liquid-coordinates", type=Path)
    parser.add_argument(
        "--combined-coordinates",
        type=Path,
        help=(
            "preassembled coordinates with slab atoms first; copy to GRO and "
            "infer the liquid molecule count"
        ),
    )
    parser.add_argument(
        "--output-directory",
        type=Path,
        default=Path("."),
        help="directory for interface.gro, system.top, forcefield.itp, slab.itp and liquid.itp",
    )
    parser.add_argument(
        "--suffix",
        type=suffix,
        default="",
        help="case suffix: writes interface_SUFFIX.gro, system_SUFFIX.top, forcefield_SUFFIX.itp, ...",
    )
    parser.add_argument("--axis", choices=("x", "y", "z"), default="z")
    parser.add_argument(
        "--interface-gap",
        type=float,
        default=0.0,
        help="optional gap at each periodic interface in nm (default: 0.0)",
    )
    parser.add_argument(
        "--coordinate-layout",
        choices=("centered-slab", "stacked"),
        default="centered-slab",
        help=(
            "coordinate arrangement: center the slab and put liquid on both sides "
            "(default), or place the two boxes sequentially"
        ),
    )
    parser.add_argument("--system-name", default="Combined molecular system")
    args = parser.parse_args()
    if args.interface_gap < 0:
        parser.error("--interface-gap cannot be negative")
    if args.combined_coordinates:
        if args.slab_coordinates or args.liquid_coordinates:
            parser.error(
                "--combined-coordinates cannot be used with --slab-coordinates "
                "or --liquid-coordinates"
            )
    elif not (args.slab_coordinates and args.liquid_coordinates):
        parser.error(
            "provide --combined-coordinates or both --slab-coordinates and "
            "--liquid-coordinates"
        )

    primary_text = args.slab.read_text()
    secondary_text = args.liquid.read_text()
    if any(name == "pairtypes" for name, _, _ in sections(secondary_text)):
        raise ValueError(
            "liquid topology contains [ pairtypes ]; explicit pairtype conversion is unsupported"
        )
    primary_defaults = defaults(primary_text)
    secondary_defaults = defaults(secondary_text)
    if primary_defaults.nbfunc != secondary_defaults.nbfunc:
        raise ValueError("input topologies use different nonbonded functions")
    if primary_defaults.combination_rule != secondary_defaults.combination_rule:
        raise ValueError("input topologies use different combination rules")

    primary_type_lines, primary_types = atomtypes(primary_text)
    secondary_type_lines, secondary_types = atomtypes(secondary_text)
    for name in primary_types.keys() & secondary_types.keys():
        if primary_types[name] != secondary_types[name]:
            raise ValueError(f"conflicting atom type definition: {name}")
    secondary_unique = [
        line for line in secondary_type_lines if line.split(";", 1)[0].split()[0] not in primary_types
    ]

    primary_composition = molecule_composition(primary_text)
    secondary_selected_name, _ = selected_molecule(secondary_text)
    primary_blocks = [
        correct_charge_rounding(molecule_block(primary_text, name))
        for name, _ in primary_composition
    ]
    secondary_source_block = correct_charge_rounding(
        molecule_block(secondary_text, secondary_selected_name)
    )
    topology_residue_name, liquid_atom_names = molecule_coordinate_signature(
        secondary_source_block
    )
    liquid_residue_name = args.liquid_residue_name or topology_residue_name
    secondary_block = explicit_secondary_pairs(
        secondary_source_block,
        secondary_defaults,
        secondary_types,
    )
    secondary_name = molecule_name(secondary_block)
    output = args.output_directory.resolve()
    output.mkdir(parents=True, exist_ok=True)
    names = output_names(args.suffix)
    coordinate_output = output / names["coordinates"]
    if args.combined_coordinates:
        primary_atom_count, secondary_count, combined_box = use_combined_coordinates(
            args.combined_coordinates,
            coordinate_output,
            liquid_residue_name,
            liquid_atom_names,
        )
    else:
        primary_atom_count, secondary_count, combined_box = assemble_coordinates(
            args.slab_coordinates,
            args.liquid_coordinates,
            coordinate_output,
            args.axis,
            args.interface_gap,
            args.coordinate_layout,
        )

    defaults_line = (
        f"{primary_defaults.nbfunc} {primary_defaults.combination_rule} "
        f"{primary_defaults.generate_pairs} {primary_defaults.fudge_lj:g} "
        f"{primary_defaults.fudge_qq:g}"
    )
    (output / names["forcefield"]).write_text(
        "[ defaults ]\n; nbfunc comb-rule gen-pairs fudgeLJ fudgeQQ\n"
        + defaults_line
        + "\n\n[ atomtypes ]\n"
        + "\n".join(primary_type_lines + secondary_unique)
        + "\n"
    )
    (output / names["slab"]).write_text(
        "\n\n".join("\n".join(block).rstrip() for block in primary_blocks) + "\n"
    )
    (output / names["liquid"]).write_text("\n".join(secondary_block).rstrip() + "\n")
    (output / names["topology"]).write_text(
        f'#include "{names["forcefield"]}"\n#include "{names["slab"]}"\n#include "{names["liquid"]}"\n\n'
        f"[ system ]\n{args.system_name}\n\n[ molecules ]\n"
        + "".join(f"{name} {count}\n" for name, count in primary_composition)
        +
        f"{secondary_name} {secondary_count}\n"
    )
    print(f"Wrote combined topology to {output}: " + ", ".join(names.values()))
    print(
        f"Slab scaling: LJ={primary_defaults.fudge_lj:g}, QQ={primary_defaults.fudge_qq:g}"
    )
    print(
        "Liquid explicit-pair scaling: "
        f"LJ={secondary_defaults.fudge_lj:g}, QQ={secondary_defaults.fudge_qq:g}"
    )
    print(
        f"Coordinates: {primary_atom_count} slab atoms, {secondary_count} liquid "
        f"residues; box={' x '.join(f'{value:.5f}' for value in combined_box)} nm"
    )
    if args.combined_coordinates:
        print("Coordinates: using preassembled slab-liquid interface")
    else:
        print(f"Configured separation at each periodic interface: {args.interface_gap:.5f} nm")
        print(f"Coordinate layout: {args.coordinate_layout}")
    return 0


def cli() -> None:
    try:
        sys.exit(main())
    except (OSError, ValueError) as error:
        print(f"Error: {error}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    cli()
