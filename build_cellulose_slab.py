#!/usr/bin/env python3
"""Build a dimension-driven cellulose slab with an installed CelluloseForge 1.1.2."""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path


VERSION = "0.3.0"
DEFAULT_FORGE = "cellulose-forge"
DEFAULT_RESTRAINT_FORCE = (1000.0, 1000.0, 1000.0)
RESTRAINT_DEFINE = "POSRES_CELLULOSE"
NO_RESTRAINT_MARKER = (
    "; No position restraints (built with --no-restraints).\n"
    f"; This empty block keeps MDP files that define {RESTRAINT_DEFINE} valid\n"
    "; without changes; it adds no restraint.\n"
    f"#ifdef {RESTRAINT_DEFINE}\n"
    "#endif\n\n"
)


def mark_unrestrained(topology: Path) -> None:
    """Insert an empty POSRES_CELLULOSE block before [ system ] in the slab topology.

    GROMACS rejects a -D macro from the MDP 'define' field that the topology never
    uses. The supplied MDP files always define POSRES_CELLULOSE, so an unrestrained
    topology carries this empty block to keep them usable unchanged.
    """
    text = topology.read_text()
    if RESTRAINT_DEFINE in text:
        raise RuntimeError(f"{topology} already refers to {RESTRAINT_DEFINE}; expected an unrestrained topology")
    lines = text.splitlines(keepends=True)
    systems = [index for index, line in enumerate(lines) if line.split(";", 1)[0].strip().lower().replace(" ", "") == "[system]"]
    if len(systems) != 1:
        raise RuntimeError(f"expected exactly one [ system ] section in {topology}, found {len(systems)}")
    lines.insert(systems[0], NO_RESTRAINT_MARKER)
    topology.write_text("".join(lines))


def cli() -> int:
    parser = argparse.ArgumentParser(prog="build-cellulose-slab")
    parser.add_argument("--version", action="version", version=f"%(prog)s {VERSION}")
    parser.add_argument("--allomorph", default="I-beta")
    parser.add_argument("--face", nargs=3, type=int, default=(1, 1, 0), metavar=("H", "K", "L"))
    parser.add_argument("--layers", type=int, required=True)
    parser.add_argument("--width", type=float, required=True, help="slab width in angstrom")
    parser.add_argument("--length", type=float, required=True, help="slab length in angstrom")
    parser.add_argument("--chain-boundary", choices=("periodic", "finite"), default="periodic")
    parser.add_argument("--vacuum", type=float, default=0.0, help="vacuum padding in angstrom")
    parser.add_argument(
        "--padding-longitudinal",
        type=float,
        help="total chain-axis padding in angstrom for finite cellulose chains",
    )
    parser.add_argument("--force-field", default="glycam", choices=("none", "glycam", "charmm36"))
    restraints = parser.add_mutually_exclusive_group(required=True)
    restraints.add_argument(
        "--restrain-layers",
        nargs="+",
        type=int,
        help="layer numbers whose heavy atoms are position-restrained",
    )
    restraints.add_argument(
        "--no-restraints",
        action="store_true",
        help="build without position restraints; the supplied MDP files still work unchanged",
    )
    parser.add_argument(
        "--restraint-force",
        nargs=3,
        type=float,
        metavar=("FX", "FY", "FZ"),
        help="restraint force constants in kJ/mol/nm^2 (default: 1000 1000 1000)",
    )
    parser.add_argument("--tleap-executable", type=Path)
    parser.add_argument("--parmed-python", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument(
        "--cellulose-forge-executable",
        default=DEFAULT_FORGE,
        help="cellulose-forge command name or path (default: cellulose-forge on PATH)",
    )
    parser.add_argument("--dry-run", action="store_true", help="print the CelluloseForge command only")
    args = parser.parse_args()

    if args.layers < 1:
        parser.error("--layers must be positive")
    if args.width <= 0 or args.length <= 0:
        parser.error("--width and --length must be positive")
    if args.no_restraints and args.restraint_force is not None:
        parser.error("--restraint-force cannot be used with --no-restraints")
    if args.no_restraints and args.force_field == "none":
        parser.error("--no-restraints needs a force-field topology; use --force-field glycam or charmm36")
    forge = shutil.which(args.cellulose_forge_executable)
    if forge is None:
        parser.error(
            f"CelluloseForge executable {args.cellulose_forge_executable!r} not found; "
            "install CelluloseForge 1.1.2 or pass --cellulose-forge-executable"
        )

    command = [
        forge,
        "build",
        "--allomorph",
        args.allomorph,
        "--shape",
        "slab",
        "--face",
        *(str(value) for value in args.face),
        "--layers",
        str(args.layers),
        "--slab-width",
        str(args.width),
        "--slab-length",
        str(args.length),
        "--chain-boundary",
        args.chain_boundary,
        "--vacuum",
        str(args.vacuum),
        "--force-field",
        args.force_field,
        "--output",
        str(args.output.resolve()),
    ]
    if args.padding_longitudinal is not None:
        command.extend(("--padding-longitudinal", str(args.padding_longitudinal)))
    if args.restrain_layers:
        force = args.restraint_force or DEFAULT_RESTRAINT_FORCE
        command.extend(("--restrain-layers", *(str(value) for value in args.restrain_layers)))
        command.extend(("--restraint-force", *(str(value) for value in force)))
    if args.tleap_executable:
        command.extend(("--tleap-executable", str(args.tleap_executable.resolve())))
    if args.parmed_python:
        command.extend(("--parmed-python", str(args.parmed_python.resolve())))
    print(" ".join(command))
    if args.no_restraints:
        print(f"No position restraints; an empty {RESTRAINT_DEFINE} block will be added to the topology.")
    if args.dry_run:
        return 0
    completed = subprocess.run(command, check=False)
    if completed.returncode != 0 or not args.no_restraints:
        return completed.returncode
    report = json.loads((args.output / "build_report.json").read_text())
    topology = args.output / report["outputs"]["gromacs_topology"]
    mark_unrestrained(topology)
    print(f"Unrestrained topology marked for the supplied MDP files: {topology}")
    return 0


if __name__ == "__main__":
    sys.exit(cli())
