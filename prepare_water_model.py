#!/usr/bin/env python3
"""Write a rigid water-model template and standalone GROMACS topology."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from project_version import project_version
from water_models import MODELS, topology_parts, write_water_pdb, write_water_topology

VERSION = project_version()


def main() -> int:
    parser = argparse.ArgumentParser(prog="prepare-water-model")
    parser.add_argument("--version", action="version", version=f"%(prog)s {VERSION}")
    parser.add_argument("--model", required=True, choices=tuple(MODELS))
    parser.add_argument("--prefix", help="output prefix; default: water_MODEL")
    parser.add_argument("--output-directory", type=Path, default=Path("."))
    args = parser.parse_args()
    output = args.output_directory.resolve()
    output.mkdir(parents=True, exist_ok=True)
    prefix = args.prefix or f"water_{args.model}"
    pdb = output / f"{prefix}.pdb"
    top = output / f"{prefix}.top"
    atomtypes_itp = output / f"{prefix}_atomtypes.itp"
    molecule_itp = output / f"{prefix}.itp"
    write_water_pdb(args.model, pdb)
    write_water_topology(args.model, top)
    atomtypes, molecule = topology_parts(args.model)
    atomtypes_itp.write_text("[ atomtypes ]\n" + atomtypes)
    molecule_itp.write_text(molecule)
    print(f"Prepared {args.model.upper()} water model:")
    for path in (pdb, top, atomtypes_itp, molecule_itp):
        print(f"  {path}")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (OSError, ValueError) as error:
        print(f"Error: {error}", file=sys.stderr)
        sys.exit(1)
