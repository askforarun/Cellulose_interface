#!/usr/bin/env python3
"""Pack a molecular liquid around a fixed, centered cellulose slab."""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

from combine_gromacs_topologies import CoordinateAtom, read_coordinates
from water_models import MODELS, WATER_MOLAR_MASS, write_water_pdb

AVOGADRO_VOLUME_FACTOR = 602.214076  # molecules nm^-3 from mol cm^-3
VERSION = "0.2.0"


def positive_float(value: str) -> float:
    number = float(value)
    if number <= 0:
        raise argparse.ArgumentTypeError("value must be positive")
    return number


def positive_int(value: str) -> int:
    number = int(value)
    if number <= 0:
        raise argparse.ArgumentTypeError("value must be a positive integer")
    return number


def greater_than_one(value: str) -> float:
    number = float(value)
    if number <= 1:
        raise argparse.ArgumentTypeError("value must be greater than 1")
    return number


def executable(value: str) -> str:
    resolved = shutil.which(value)
    if resolved is None:
        raise argparse.ArgumentTypeError(f"executable not found: {value}")
    return resolved


def residue_name(value: str) -> str:
    name = value.upper()
    if not name.isalnum() or not 1 <= len(name) <= 3:
        raise argparse.ArgumentTypeError("residue name must contain 1-3 letters or digits")
    return name


def molar_mass(smiles: str, obabel: str) -> float:
    result = subprocess.run(
        [obabel, f"-:{smiles}", "-osmi", "--append", "MW"],
        check=True,
        capture_output=True,
        text=True,
    )
    fields = result.stdout.split()
    if len(fields) < 2:
        raise RuntimeError("Open Babel did not return a molecular weight")
    return float(fields[-1])


def normalize_template_names(pdb_path: Path, residue: str) -> None:
    """Give each element a unique, topology-friendly name within the molecule."""
    counts: defaultdict[str, int] = defaultdict(int)
    output: list[str] = []
    for line in pdb_path.read_text().splitlines(keepends=True):
        if line.startswith(("ATOM  ", "HETATM")):
            element = line[76:78].strip().upper()
            if not element:
                element = "".join(character for character in line[12:16] if character.isalpha())
                element = element[:2].upper()
            counts[element] += 1
            atom_name = f"{element}{counts[element]}"
            if len(atom_name) > 4:
                raise RuntimeError(f"cannot represent atom name in PDB format: {atom_name}")
            line = line[:12] + f"{atom_name:>4}" + line[16:17] + f"{residue:>3}" + line[20:]
        output.append(line)
    pdb_path.write_text("".join(output))


def element_from_atom_name(name: str) -> str:
    letters = "".join(character for character in name if character.isalpha())
    if not letters:
        return "X"
    return letters[0].upper() + (letters[1].lower() if len(letters) > 1 else "")


def write_centered_cellulose(
    atoms: list[CoordinateAtom], box: list[float], final_z: float, output: Path
) -> tuple[float, float]:
    z_min = min(atom.xyz[2] for atom in atoms)
    z_max = max(atom.xyz[2] for atom in atoms)
    thickness = z_max - z_min
    shift_z = (final_z - thickness) / 2 - z_min
    lines = [
        f"CRYST1{box[0] * 10:9.3f}{box[1] * 10:9.3f}{final_z * 10:9.3f}"
        "  90.00  90.00  90.00 P 1           1"
    ]
    for serial, atom in enumerate(atoms, 1):
        x = (atom.xyz[0] % box[0]) * 10
        y = (atom.xyz[1] % box[1]) * 10
        z = (atom.xyz[2] + shift_z) * 10
        element = element_from_atom_name(atom.atomname)
        lines.append(
            f"ATOM  {serial % 100000:5d} {atom.atomname:>4.4s} "
            f"{atom.resname:>3.3s} {atom.resid % 10000:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {element:>2.2s}"
        )
    lines.extend(("TER", "END"))
    output.write_text("\n".join(lines) + "\n")
    return thickness, shift_z


def add_cryst1(path: Path, box: list[float]) -> None:
    cryst1 = (
        f"CRYST1{box[0] * 10:9.3f}{box[1] * 10:9.3f}{box[2] * 10:9.3f}"
        "  90.00  90.00  90.00 P 1           1\n"
    )
    lines = path.read_text().splitlines(keepends=True)
    path.write_text(cryst1 + "".join(line for line in lines if not line.startswith("CRYST1")))


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(
        prog="pack-liquid-cellulose",
        description=(
            "Center a cellulose slab and use periodic Packmol packing to place a "
            "molecular liquid on both sides. Only combined coordinates are written."
        ),
    )
    result.add_argument("--version", action="version", version=f"%(prog)s {VERSION}")
    result.add_argument("--cellulose-coordinates", required=True, type=Path)
    liquid = result.add_mutually_exclusive_group(required=True)
    liquid.add_argument("--smiles", help='molecular-liquid SMILES, e.g. "CCO"')
    liquid.add_argument(
        "--water-model",
        choices=tuple(MODELS),
        help="rigid water preset (uses model-specific geometry and virtual sites)",
    )
    result.add_argument(
        "--residue-name",
        type=residue_name,
        help="molecular-liquid residue name (default: MOL); water presets use SOL",
    )
    result.add_argument("--molecules", required=True, type=positive_int)
    result.add_argument("--density", required=True, type=positive_float, help="g/cm^3")
    result.add_argument("--output", type=Path, default=Path("interface.pdb"))
    result.add_argument("--tolerance-angstrom", type=positive_float, default=2.0)
    result.add_argument("--maxit", type=positive_int, default=50)
    result.add_argument("--nloop", type=positive_int, default=5)
    result.add_argument("--max-attempts", type=positive_int, default=5)
    result.add_argument("--z-growth-factor", type=greater_than_one, default=1.05)
    result.add_argument("--seed", type=int, default=1234567)
    result.add_argument("--packmol", type=executable, default="packmol")
    result.add_argument("--obabel", type=executable, default="obabel")
    result.add_argument("--dry-run", action="store_true")
    return result


def main() -> int:
    args = parser().parse_args()
    if args.molecules % 2:
        raise ValueError(
            "--molecules must be even so the two cellulose surfaces receive equal counts"
        )
    cellulose_atoms, cellulose_box = read_coordinates(args.cellulose_coordinates)
    if args.water_model:
        mass = WATER_MOLAR_MASS
        liquid_name = args.water_model.upper()
        residue = "SOL"
        if args.residue_name and args.residue_name != residue:
            raise ValueError("water presets require --residue-name SOL")
    else:
        mass = molar_mass(args.smiles, args.obabel)
        liquid_name = args.smiles
        residue = args.residue_name or "MOL"
    slab_z_min = min(atom.xyz[2] for atom in cellulose_atoms)
    slab_z_max = max(atom.xyz[2] for atom in cellulose_atoms)
    slab_thickness = slab_z_max - slab_z_min
    liquid_z = (
        args.molecules * mass
        / (args.density * AVOGADRO_VOLUME_FACTOR * cellulose_box[0] * cellulose_box[1])
    )
    initial_z = slab_thickness + liquid_z
    print(f"Cellulose thickness: {slab_thickness:.6f} nm")
    print(f"Liquid thickness:    {liquid_z:.6f} nm")
    print(f"Initial final box:   {cellulose_box[0]:.6f} x {cellulose_box[1]:.6f} x {initial_z:.6f} nm")
    print(f"Target liquid density: {args.density:.6f} g/cm^3 ({args.molecules} molecules)")
    print(f"Molecules per side:    {args.molecules // 2}")
    print(f"Liquid molar mass:     {mass:.6f} g/mol")
    print(f"Liquid definition:     {liquid_name}")
    if args.dry_run:
        return 0

    output = args.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    final_z = initial_z
    last_result: subprocess.CompletedProcess[str] | None = None
    with tempfile.TemporaryDirectory(prefix="pack_liquid_cellulose_") as temporary_name:
        temporary = Path(temporary_name)
        centered_cellulose = temporary / "cellulose_centered.pdb"
        liquid_template = temporary / "liquid.pdb"
        if args.water_model:
            write_water_pdb(args.water_model, liquid_template)
        else:
            subprocess.run(
                [
                    args.obabel,
                    f"-:{args.smiles}",
                    "--gen3d",
                    "-h",
                    "-opdb",
                    "-O",
                    str(liquid_template),
                ],
                check=True,
            )
            normalize_template_names(liquid_template, residue)
        for attempt in range(1, args.max_attempts + 1):
            final_box = [cellulose_box[0], cellulose_box[1], final_z]
            write_centered_cellulose(cellulose_atoms, cellulose_box, final_z, centered_cellulose)
            liquid_region = final_z - slab_thickness
            lower_surface = liquid_region / 2
            upper_surface = lower_surface + slab_thickness
            molecules_per_side = args.molecules // 2
            packmol_input = temporary / "packmol.inp"
            packmol_input.write_text(
                f"""tolerance {args.tolerance_angstrom:.6f}
filetype pdb
output {output}
seed {args.seed}
nloop {args.nloop}
maxit {args.maxit}
pbc {final_box[0] * 10:.6f} {final_box[1] * 10:.6f} {final_box[2] * 10:.6f}

structure {centered_cellulose}
  number 1
  fixed 0.0 0.0 0.0 0.0 0.0 0.0
  resnumbers 1
end structure

structure {liquid_template}
  number {molecules_per_side}
  inside box 0.0 0.0 0.0 {final_box[0] * 10:.6f} {final_box[1] * 10:.6f} {lower_surface * 10:.6f}
  resnumbers 3
end structure

structure {liquid_template}
  number {molecules_per_side}
  inside box 0.0 0.0 {upper_surface * 10:.6f} {final_box[0] * 10:.6f} {final_box[1] * 10:.6f} {final_box[2] * 10:.6f}
  resnumbers 3
end structure
"""
            )
            print(f"Packmol attempt {attempt}/{args.max_attempts}: z = {final_z:.6f} nm")
            output.unlink(missing_ok=True)
            with packmol_input.open() as input_file:
                last_result = subprocess.run(
                    [args.packmol], stdin=input_file, capture_output=True, text=True
                )
            if last_result.returncode == 0 and "Success!" in last_result.stdout:
                add_cryst1(output, final_box)
                print(f"Packed interface coordinates written to: {output}")
                return 0
            if attempt < args.max_attempts:
                final_z *= args.z_growth_factor
                print("Packing failed; increasing the total z dimension and retrying.")
    assert last_result is not None
    details = (last_result.stderr or last_result.stdout)[-2000:]
    raise RuntimeError(
        f"Packmol failed after {args.max_attempts} attempts. Last output:\n{details}"
    )


def cli() -> None:
    try:
        sys.exit(main())
    except (OSError, RuntimeError, subprocess.CalledProcessError, ValueError) as error:
        print(f"Error: {error}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    cli()
