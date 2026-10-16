#!/usr/bin/env python3
"""Build a dimension-driven cellulose slab with a local CelluloseForge checkout."""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path


VERSION = "0.1.0"
DEFAULT_FORGE = Path(__file__).resolve().parents[1] / "cellulose_functionalization"
DEFAULT_UV = shutil.which("uv") or str(Path.home() / ".local" / "bin" / "uv")


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
    parser.add_argument("--force-field", default="glycam", choices=("none", "glycam", "charmm36"))
    parser.add_argument("--restrain-layers", nargs="+", type=int)
    parser.add_argument("--restraint-force", nargs=3, type=float, default=(1000.0, 1000.0, 1000.0))
    parser.add_argument("--tleap-executable", type=Path)
    parser.add_argument("--parmed-python", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--cellulose-forge-directory", type=Path, default=DEFAULT_FORGE)
    parser.add_argument("--uv", type=Path, default=Path(DEFAULT_UV), help="path to the uv executable")
    parser.add_argument("--dry-run", action="store_true", help="print the CelluloseForge command only")
    args = parser.parse_args()

    if args.layers < 1:
        parser.error("--layers must be positive")
    if args.width <= 0 or args.length <= 0:
        parser.error("--width and --length must be positive")
    forge = args.cellulose_forge_directory.resolve()
    if not (forge / "pyproject.toml").is_file():
        parser.error(f"CelluloseForge project not found at {forge}")
    uv = args.uv.resolve()
    if not uv.is_file():
        parser.error(f"uv executable not found at {uv}")

    command = [
        str(uv),
        "run",
        "--project",
        str(forge),
        "cellulose-forge",
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
    if args.restrain_layers:
        command.extend(("--restrain-layers", *(str(value) for value in args.restrain_layers)))
        command.extend(("--restraint-force", *(str(value) for value in args.restraint_force)))
    if args.tleap_executable:
        command.extend(("--tleap-executable", str(args.tleap_executable.resolve())))
    if args.parmed_python:
        command.extend(("--parmed-python", str(args.parmed_python.resolve())))
    print(" ".join(command))
    if args.dry_run:
        return 0
    completed = subprocess.run(command, check=False)
    return completed.returncode


if __name__ == "__main__":
    sys.exit(cli())
