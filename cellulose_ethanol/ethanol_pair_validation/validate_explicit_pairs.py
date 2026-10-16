#!/usr/bin/env python3
"""Validate explicit GROMACS type-2 pairs against ParmEd's native GAFF2 topology."""

from __future__ import annotations

import argparse
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "ethanol_parameters"
PAIR_HEADER = "; ai aj funct fudgeQQ qi qj sigma14 epsilon14"
TERMS = ("Bond", "Angle", "Proper-Dih.", "LJ-14", "Coulomb-14", "Potential")


def run(command: list[str], cwd: Path, stdin: str | None = None) -> str:
    result = subprocess.run(
        command,
        cwd=cwd,
        input=stdin,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    if result.returncode:
        raise RuntimeError(f"Command failed ({result.returncode}): {' '.join(command)}\n{result.stdout}")
    return result.stdout


def topology_parts(text: str) -> tuple[str, str]:
    molecule = text.index("[ moleculetype ]")
    system = text.index("[ system ]")
    return text[:molecule].rstrip(), text[molecule:system].rstrip()


def explicit_forcefield(native_forcefield: str) -> str:
    return re.sub(
        r"(?m)^(1\s+2\s+yes\s+)0\.5\s+0\.83333333\s*$",
        r"\g<1>1.0          1.0",
        native_forcefield,
        count=1,
    )


def make_topologies(work: Path) -> None:
    native_text = (SOURCE / "ethanol.top").read_text()
    native_forcefield, native_molecule = topology_parts(native_text)
    explicit_molecule = (SOURCE / "ethanol.itp").read_text().rstrip()
    if PAIR_HEADER not in explicit_molecule:
        raise ValueError("source ethanol.itp does not contain explicit type-2 pairs")
    tail = "\n\n[ system ]\nEthanol pair-scaling validation\n\n[ molecules ]\nETO 1\n"
    (work / "native.top").write_text(native_forcefield + "\n\n" + native_molecule + tail)
    (work / "explicit.top").write_text(
        explicit_forcefield(native_forcefield) + "\n\n" + explicit_molecule + tail
    )


def read_xvg(path: Path) -> list[float]:
    rows = [line for line in path.read_text().splitlines() if line and line[0] not in "#@"]
    if len(rows) != 1:
        raise ValueError(f"expected one energy frame in {path}, found {len(rows)}")
    return [float(value) for value in rows[0].split()[1:]]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gmx", default=shutil.which("gmx") or "gmx")
    parser.add_argument("--tolerance", type=float, default=1.0e-5)
    args = parser.parse_args()
    work = Path(__file__).resolve().parent
    make_topologies(work)
    shutil.copyfile(SOURCE / "ethanol.gro", work / "source.gro")
    run(
        [args.gmx, "editconf", "-f", "source.gro", "-o", "ethanol_3nm.gro", "-box", "3", "3", "3", "-c"],
        work,
    )
    (work / "single_point.mdp").write_text(
        "integrator = md\n"
        "nsteps = 0\n"
        "dt = 0.001\n"
        "cutoff-scheme = Verlet\n"
        "nstlist = 1\n"
        "rlist = 1.0\n"
        "coulombtype = Cut-off\n"
        "rcoulomb = 1.0\n"
        "vdwtype = Cut-off\n"
        "rvdw = 1.0\n"
        "pbc = xyz\n"
        "nstenergy = 1\n"
        "continuation = no\n"
        "gen-vel = no\n"
    )
    # These are terms 1, 2, 3, 4, 5, and 8 in the GROMACS 2024 energy menu.
    selection = "1 2 3 4 5 8\n0\n"
    logs: dict[str, str] = {}
    energies: dict[str, dict[str, float]] = {}
    for case in ("native", "explicit"):
        logs[f"grompp_{case}"] = run(
            [args.gmx, "grompp", "-f", "single_point.mdp", "-c", "ethanol_3nm.gro", "-p", f"{case}.top", "-o", f"{case}.tpr", "-po", f"{case}.mdout.mdp", "-maxwarn", "0"],
            work,
        )
        logs[f"mdrun_{case}"] = run(
            [args.gmx, "mdrun", "-s", f"{case}.tpr", "-deffnm", case, "-nt", "1"], work
        )
        logs[f"energy_{case}"] = run(
            [args.gmx, "energy", "-f", f"{case}.edr", "-o", f"{case}_energy.xvg", "-xvg", "none"],
            work,
            selection,
        )
        values = read_xvg(work / f"{case}_energy.xvg")
        if len(values) != len(TERMS):
            raise ValueError(f"expected {len(TERMS)} selected terms, found {len(values)}")
        energies[case] = dict(zip(TERMS, values, strict=True))
    differences = {term: energies["explicit"][term] - energies["native"][term] for term in TERMS}
    passed = all(abs(value) <= args.tolerance for value in differences.values())
    report = {
        "test": "ParmEd native GAFF2 type-1 pairs versus explicit GROMACS type-2 pairs",
        "units": "kJ/mol",
        "tolerance": args.tolerance,
        "passed": passed,
        "native": energies["native"],
        "explicit": energies["explicit"],
        "explicit_minus_native": differences,
    }
    (work / "validation_report.json").write_text(json.dumps(report, indent=2) + "\n")
    (work / "commands.log").write_text("\n\n".join(f"## {key}\n{value}" for key, value in logs.items()))
    print(json.dumps(report, indent=2))
    return 0 if passed else 1


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (OSError, RuntimeError, ValueError) as error:
        print(f"Error: {error}", file=sys.stderr)
        sys.exit(2)
