#!/usr/bin/env python3
"""Parameterize one small molecule with GAFF2/AM1-BCC for GROMACS."""

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

from project_version import project_version

from combine_gromacs_topologies import (
    atomtypes,
    defaults,
    explicit_secondary_pairs,
    molecule_block,
    selected_molecule,
)
VERSION = project_version()


def generate_conformer(smiles: str, seed: int, output: Path) -> float:
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except ImportError as error:
        raise RuntimeError("RDKit is required for deterministic conformer generation") from error
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError(f"RDKit could not parse SMILES: {smiles!r}")
    molecule = Chem.AddHs(molecule)
    parameters = AllChem.ETKDGv3()
    parameters.randomSeed = seed
    if AllChem.EmbedMolecule(molecule, parameters) != 0:
        raise RuntimeError("RDKit ETKDG conformer generation failed")
    properties = AllChem.MMFFGetMoleculeProperties(molecule)
    if properties is None:
        raise RuntimeError("RDKit MMFF parameters are unavailable for this molecule")
    status = AllChem.MMFFOptimizeMolecule(molecule, mmffVariant="MMFF94")
    if status != 0:
        raise RuntimeError(f"RDKit MMFF conformer optimization did not converge (status {status})")
    force_field = AllChem.MMFFGetMoleculeForceField(molecule, properties)
    if force_field is None:
        raise RuntimeError("RDKit could not construct the MMFF force field")
    energy = float(force_field.CalcEnergy())
    Chem.MolToPDBFile(molecule, str(output))
    return energy


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


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(
        prog="parameterize-molecule",
        description=(
            "Generate a 3D molecule from SMILES, assign GAFF2 atom types and "
            "AM1-BCC charges, and write Amber and GROMACS topologies."
        ),
    )
    result.add_argument("--version", action="version", version=f"%(prog)s {VERSION}")
    result.add_argument("--smiles", required=True, help='molecular SMILES, e.g. "CCO"')
    result.add_argument("--residue-name", type=residue_name, default="MOL")
    result.add_argument("--net-charge", type=int, default=0)
    result.add_argument("--seed", type=int, default=2026, help="RDKit ETKDG random seed")
    result.add_argument("--prefix", default="molecule", help="output filename prefix")
    result.add_argument("--output-directory", type=Path, default=Path("."))
    result.add_argument("--antechamber", type=executable, default="antechamber")
    result.add_argument("--parmchk2", type=executable, default="parmchk2")
    result.add_argument("--tleap", type=executable, default="tleap")
    return result


def run(args: argparse.Namespace) -> tuple[Path, ...]:
    try:
        import parmed as pmd
    except ImportError as error:
        raise RuntimeError(
            "ParmEd is unavailable; run this command inside the AmberTools environment"
        ) from error

    output_directory = args.output_directory.resolve()
    output_directory.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix="parameterize_molecule_") as temporary_name:
        temporary = Path(temporary_name)
        input_pdb = temporary / "input.pdb"
        mol2 = temporary / "molecule.mol2"
        frcmod = temporary / "molecule.frcmod"
        prmtop = temporary / "molecule.prmtop"
        inpcrd = temporary / "molecule.inpcrd"

        conformer_energy = generate_conformer(args.smiles, args.seed, input_pdb)
        normalize_template_names(input_pdb, args.residue_name)
        subprocess.run(
            [
                args.antechamber,
                "-i",
                str(input_pdb),
                "-fi",
                "pdb",
                "-o",
                str(mol2),
                "-fo",
                "mol2",
                "-at",
                "gaff2",
                "-c",
                "bcc",
                "-nc",
                str(args.net_charge),
                "-rn",
                args.residue_name,
                "-s",
                "2",
            ],
            cwd=temporary,
            check=True,
        )
        subprocess.run(
            [
                args.parmchk2,
                "-i",
                str(mol2),
                "-f",
                "mol2",
                "-o",
                str(frcmod),
                "-s",
                "gaff2",
            ],
            cwd=temporary,
            check=True,
        )
        if "ATTN, need revision" in frcmod.read_text():
            raise RuntimeError("parmchk2 reported 'ATTN, need revision' in the frcmod file")

        leap_input = temporary / "tleap.in"
        leap_input.write_text(
            "\n".join(
                (
                    "source leaprc.gaff2",
                    f"loadamberparams {frcmod}",
                    f"MOL = loadmol2 {mol2}",
                    "check MOL",
                    f"saveamberparm MOL {prmtop} {inpcrd}",
                    "quit",
                    "",
                )
            )
        )
        tleap_result = subprocess.run(
            [args.tleap, "-f", str(leap_input)],
            cwd=temporary,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            check=True,
        )
        if re.search(r"Errors\s*=\s*0\b", tleap_result.stdout) is None:
            raise RuntimeError("tleap output did not report 'Errors = 0'")

        destinations = {
            input_pdb: output_directory / f"{args.prefix}.pdb",
            mol2: output_directory / f"{args.prefix}.mol2",
            frcmod: output_directory / f"{args.prefix}.frcmod",
            prmtop: output_directory / f"{args.prefix}.prmtop",
            inpcrd: output_directory / f"{args.prefix}.inpcrd",
        }
        for source, destination in destinations.items():
            shutil.copy2(source, destination)

        structure = pmd.load_file(str(prmtop), xyz=str(inpcrd))
        gromacs_top = output_directory / f"{args.prefix}.top"
        gromacs_gro = output_directory / f"{args.prefix}.gro"
        structure.save(str(gromacs_top), overwrite=True)
        structure.save(str(gromacs_gro), overwrite=True)

        topology_text = gromacs_top.read_text()
        topology_defaults = defaults(topology_text)
        atomtype_lines, atomtype_parameters = atomtypes(topology_text)
        selected_name, _ = selected_molecule(topology_text)
        explicit_block = explicit_secondary_pairs(
            molecule_block(topology_text, selected_name),
            topology_defaults,
            atomtype_parameters,
        )
        atomtypes_itp = output_directory / f"{args.prefix}_atomtypes.itp"
        molecule_itp = output_directory / f"{args.prefix}.itp"
        atomtypes_itp.write_text(
            "; Generated GAFF2 atom types; include before any [ moleculetype ]\n"
            "[ atomtypes ]\n"
            + "\n".join(atomtype_lines)
            + "\n"
        )
        molecule_itp.write_text(
            "; GAFF2 molecule with explicit type-2 1-4 pairs\n"
            + "\n".join(explicit_block).rstrip()
            + "\n"
        )

        report = output_directory / f"{args.prefix}_parameterization_report.json"
        report.write_text(
            json.dumps(
                {
                    "smiles": args.smiles,
                    "residue_name": args.residue_name,
                    "net_charge": args.net_charge,
                    "rdkit_etkdg_seed": args.seed,
                    "mmff_energy_kcal_mol": conformer_energy,
                    "tleap_errors": 0,
                },
                indent=2,
            )
            + "\n"
        )

    return (
        *destinations.values(),
        gromacs_top,
        gromacs_gro,
        atomtypes_itp,
        molecule_itp,
        report,
    )


def main() -> int:
    args = parser().parse_args()
    outputs = run(args)
    print("Parameterization completed successfully:")
    for output in outputs:
        print(f"  {output}")
    return 0


def cli() -> None:
    try:
        sys.exit(main())
    except (OSError, RuntimeError, subprocess.CalledProcessError, ValueError) as error:
        print(f"Error: {error}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    cli()
