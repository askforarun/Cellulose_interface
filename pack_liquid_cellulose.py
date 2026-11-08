#!/usr/bin/env python3
"""Pack a molecular liquid around a fixed, centered cellulose slab."""

from __future__ import annotations

import argparse
import json
import re
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

from combine_gromacs_topologies import CoordinateAtom, read_coordinates
from project_version import project_version
from water_models import MODELS, WATER_MOLAR_MASS, write_water_pdb

AVOGADRO_VOLUME_FACTOR = 602.214076  # molecules nm^-3 from mol cm^-3
VERSION = project_version()


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


def fraction(value: str) -> float:
    number = float(value)
    if not 0 <= number < 1:
        raise argparse.ArgumentTypeError("value must be in the interval [0, 1)")
    return number


def liquid_volume_nm3(box_x: float, box_y: float, final_z: float, slab_thickness: float) -> float:
    return box_x * box_y * (final_z - slab_thickness)


def density_from_volume(molecules: int, molar_mass: float, volume_nm3: float) -> float:
    return molecules * molar_mass / (AVOGADRO_VOLUME_FACTOR * volume_nm3)


def command_version(command: list[str], *, packmol: bool = False) -> str:
    try:
        result = subprocess.run(
            command, stdin=subprocess.DEVNULL, capture_output=True, text=True,
            check=False, timeout=10,
        )
    except (OSError, subprocess.TimeoutExpired):
        return "unknown"
    output = result.stdout + "\n" + result.stderr
    if packmol:
        match = re.search(r"Version\s+(\S+)", output, re.IGNORECASE)
        return match.group(1) if match else "unknown"
    return next((line.strip() for line in output.splitlines() if line.strip()), "unknown")


def executable(value: str) -> str:
    resolved = shutil.which(value)
    if resolved is None:
        raise argparse.ArgumentTypeError(f"executable not found: {value}")
    return resolved


def suffix(value: str) -> str:
    if value and not all(character.isalnum() or character in "-_" for character in value):
        raise argparse.ArgumentTypeError("suffix may contain only letters, digits, '-' and '_'")
    return value


def residue_name(value: str) -> str:
    name = value.upper()
    if not name.isalnum() or not 1 <= len(name) <= 3:
        raise argparse.ArgumentTypeError("residue name must contain 1-3 letters or digits")
    return name


def molar_mass(smiles: str) -> float:
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors
    except ImportError as error:
        raise RuntimeError("RDKit is required to calculate molecular mass") from error
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError(f"RDKit could not parse SMILES: {smiles!r}")
    return float(Descriptors.MolWt(molecule))


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
        "--liquid-template",
        type=Path,
        help="parameterized liquid PDB template; with --smiles, preferred over Open Babel",
    )
    result.add_argument(
        "--residue-name",
        type=residue_name,
        help="molecular-liquid residue name (default: MOL); water presets use SOL",
    )
    result.add_argument("--molecules", required=True, type=positive_int)
    result.add_argument("--density", required=True, type=positive_float, help="g/cm^3")
    result.add_argument(
        "--output-directory",
        type=Path,
        default=Path("."),
        help="directory for the packed interface.pdb (default: current directory)",
    )
    result.add_argument(
        "--suffix",
        type=suffix,
        default="",
        help="case suffix: writes interface_SUFFIX.pdb, e.g. --suffix 3000",
    )
    result.add_argument("--tolerance-angstrom", type=positive_float, default=2.0)
    result.add_argument("--maxit", type=positive_int, default=50)
    result.add_argument("--nloop", type=positive_int, default=5)
    result.add_argument("--max-attempts", type=positive_int, default=5)
    result.add_argument("--z-growth-factor", type=greater_than_one, default=1.05)
    result.add_argument(
        "--max-density-loss",
        type=fraction,
        default=0.05,
        help="maximum fractional density loss allowed after retries (default: 0.05)",
    )
    result.add_argument("--seed", type=int, default=1234567)
    result.add_argument("--packmol", type=executable, default="packmol")
    result.add_argument(
        "--obabel",
        default="obabel",
        help="Open Babel executable used only when --liquid-template is omitted",
    )
    result.add_argument("--dry-run", action="store_true")
    return result


def main() -> int:
    args = parser().parse_args()
    if args.liquid_template and args.water_model:
        raise ValueError("--liquid-template cannot be used with --water-model")
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
        mass = molar_mass(args.smiles)
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
        for attempt in range(1, args.max_attempts + 1):
            attempt_liquid_z = liquid_z * args.z_growth_factor ** (attempt - 1)
            attempt_z = slab_thickness + attempt_liquid_z
            volume = liquid_volume_nm3(
                cellulose_box[0], cellulose_box[1], attempt_z, slab_thickness
            )
            achieved = density_from_volume(args.molecules, mass, volume)
            status = "allowed" if achieved >= args.density * (1 - args.max_density_loss) else "refused"
            print(
                f"Attempt {attempt}: liquid thickness={attempt_liquid_z:.6f} nm, "
                f"box Z={attempt_z:.6f} nm, density={achieved:.6f} g/cm^3 ({status})"
            )
            if status == "refused":
                break
        return 0

    output = args.output_directory.resolve() / f"interface{'_' + args.suffix if args.suffix else ''}.pdb"
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        dir=output.parent,
        prefix=f".{output.stem}.",
        suffix=".pdb",
        delete=False,
    ) as candidate_file:
        candidate = Path(candidate_file.name)
    candidate.unlink()
    final_z = initial_z
    attempts_used = 0
    last_numbers: tuple[int, float, float, float] | None = None
    last_result: subprocess.CompletedProcess[str] | None = None
    obabel_used = False
    try:
        with tempfile.TemporaryDirectory(prefix="pack_liquid_cellulose_") as temporary_name:
            temporary = Path(temporary_name)
            centered_cellulose = temporary / "cellulose_centered.pdb"
            liquid_template = temporary / "liquid.pdb"
            if args.water_model:
                write_water_pdb(args.water_model, liquid_template)
            elif args.liquid_template:
                if not args.liquid_template.is_file():
                    raise ValueError(f"liquid template does not exist: {args.liquid_template}")
                shutil.copy2(args.liquid_template, liquid_template)
                template_records = [
                    line
                    for line in liquid_template.read_text().splitlines()
                    if line.startswith(("ATOM  ", "HETATM"))
                ]
                template_names = [line[12:16].strip() for line in template_records]
                if not template_names or any(not name for name in template_names):
                    raise ValueError("liquid template must contain named PDB atoms")
                if len(set(template_names)) != len(template_names):
                    raise ValueError("liquid template atom names must be unique within the molecule")
            else:
                args.obabel = executable(args.obabel)
                obabel_used = True
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
            template_atom_count = sum(
                line.startswith(("ATOM  ", "HETATM"))
                for line in liquid_template.read_text().splitlines()
            )
            if template_atom_count == 0:
                raise RuntimeError(f"liquid template contains no atoms: {liquid_template}")
            expected_atom_count = len(cellulose_atoms) + args.molecules * template_atom_count
            for attempt in range(1, args.max_attempts + 1):
                attempts_used = attempt
                liquid_region = liquid_z * args.z_growth_factor ** (attempt - 1)
                final_z = slab_thickness + liquid_region
                volume = liquid_volume_nm3(
                    cellulose_box[0], cellulose_box[1], final_z, slab_thickness
                )
                achieved_density = density_from_volume(args.molecules, mass, volume)
                last_numbers = (attempt, final_z, volume, achieved_density)
                minimum_density = args.density * (1 - args.max_density_loss)
                if achieved_density < minimum_density:
                    raise RuntimeError(
                        "refusing Packmol attempt because the density-loss bound would be exceeded: "
                        f"attempt={attempt}, final_z={final_z:.6f} nm, "
                        f"liquid_volume={volume:.6f} nm^3, "
                        f"achieved_density={achieved_density:.6f} g/cm^3, "
                        f"minimum_density={minimum_density:.6f} g/cm^3"
                    )
                final_box = [cellulose_box[0], cellulose_box[1], final_z]
                write_centered_cellulose(
                    cellulose_atoms, cellulose_box, final_z, centered_cellulose
                )
                lower_surface = liquid_region / 2
                upper_surface = lower_surface + slab_thickness
                molecules_per_side = args.molecules // 2
                packmol_input = temporary / "packmol.inp"
                packmol_input.write_text(
                    f"""tolerance {args.tolerance_angstrom:.6f}
filetype pdb
output {candidate}
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
                candidate.unlink(missing_ok=True)
                with packmol_input.open() as input_file:
                    last_result = subprocess.run(
                        [args.packmol], stdin=input_file, capture_output=True, text=True
                    )
                if last_result.returncode == 0 and "Success!" in last_result.stdout:
                    if not candidate.is_file():
                        raise RuntimeError("Packmol reported success but wrote no coordinate file")
                    add_cryst1(candidate, final_box)
                    packed_atoms, packed_box = read_coordinates(candidate)
                    if len(packed_atoms) != expected_atom_count:
                        raise RuntimeError(
                            "Packmol atom-count mismatch: "
                            f"expected {expected_atom_count}, found {len(packed_atoms)}"
                        )
                    if any(abs(found - expected) > 1e-4 for found, expected in zip(packed_box, final_box)):
                        raise RuntimeError(
                            f"packed CRYST1 box mismatch: expected {final_box}, found {packed_box}"
                        )
                    report = {
                        "requested_density_g_cm3": args.density,
                        "molecules": args.molecules,
                        "molar_mass_g_mol": mass,
                        "slab_thickness_nm": slab_thickness,
                        "initial_z_nm": initial_z,
                        "final_z_nm": final_z,
                        "attempts_used": attempts_used,
                        "liquid_volume_nm3": volume,
                        "achieved_density_g_cm3": achieved_density,
                        "packmol_seed": args.seed,
                        "packmol_settings": {
                            "tolerance_angstrom": args.tolerance_angstrom,
                            "maxit": args.maxit,
                            "nloop": args.nloop,
                            "max_attempts": args.max_attempts,
                            "liquid_z_growth_factor": args.z_growth_factor,
                            "max_density_loss": args.max_density_loss,
                        },
                        "tool_versions": {
                            "packmol": command_version([args.packmol], packmol=True),
                            "openbabel": (
                                command_version([args.obabel, "-V"])
                                if obabel_used
                                else "not used"
                            ),
                            "pack_liquid_cellulose": VERSION,
                        },
                    }
                    report_path = output.with_name(output.stem + "_packing_report.json")
                    with tempfile.NamedTemporaryFile(
                        mode="w", dir=output.parent, prefix=f".{report_path.stem}.",
                        suffix=".json", delete=False,
                    ) as report_file:
                        report_candidate = Path(report_file.name)
                        try:
                            report_file.write(json.dumps(report, indent=2) + "\n")
                        except BaseException:
                            report_candidate.unlink(missing_ok=True)
                            raise
                    try:
                        report_candidate.replace(report_path)
                    finally:
                        report_candidate.unlink(missing_ok=True)
                    candidate.replace(output)
                    print(f"Requested density:     {args.density:.6f} g/cm^3")
                    print(f"Final box Z:           {final_z:.6f} nm")
                    print(f"Liquid volume:         {volume:.6f} nm^3")
                    print(f"Achieved density:      {achieved_density:.6f} g/cm^3")
                    print(f"Packed interface coordinates written to: {output}")
                    print(f"Packing report written to: {report_path}")
                    return 0
                if attempt < args.max_attempts:
                    print(
                        "Packing failed; increasing only the liquid thickness by "
                        f"a factor of {args.z_growth_factor:g} and retrying."
                    )
    finally:
        candidate.unlink(missing_ok=True)
    assert last_result is not None
    details = (last_result.stderr or last_result.stdout)[-2000:]
    numbers = ""
    if last_numbers is not None:
        attempt, attempted_z, attempted_volume, attempted_density = last_numbers
        numbers = (
            f" Last attempt: attempt={attempt}, final_z={attempted_z:.6f} nm, "
            f"liquid_volume={attempted_volume:.6f} nm^3, "
            f"achieved_density={attempted_density:.6f} g/cm^3."
        )
    raise RuntimeError(
        f"Packmol failed after {attempts_used} attempts.{numbers} Last output:\n{details}"
    )


def cli() -> None:
    try:
        sys.exit(main())
    except (OSError, RuntimeError, subprocess.CalledProcessError, ValueError) as error:
        print(f"Error: {error}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    cli()
