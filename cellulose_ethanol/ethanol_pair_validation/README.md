# Ethanol 1-4 pair validation

This isolated test compares the same GAFF2/AM1-BCC ethanol coordinates and parameters in two GROMACS representations:

1. `native.top`: ParmEd's ordinary type-1 pairs with global `fudgeLJ = 0.5` and `fudgeQQ = 0.83333333`.
2. `explicit.top`: global `fudgeLJ = fudgeQQ = 1.0`, matching the intended GLYCAM host topology, with GAFF2 scaling encoded in every type-2 ethanol pair.

Both systems are evaluated as one ethanol molecule in the same 3 nm periodic box. The script compares bonded energy, proper-dihedral energy, LJ-14, Coulomb-14, and total potential energy at step zero.

Run from this directory or elsewhere:

```console
python validate_explicit_pairs.py
```

The machine-readable result is written to `validation_report.json`; complete GROMACS output is retained in `commands.log`.
