# Cellulose–liquid interface workflow

Build a crystalline cellulose slab with **GLYCAM**, parameterize a small liquid molecule, pack liquid on both surfaces, and prepare a combined GROMACS system. The example below uses cellulose I-beta (110) and ethanol.

The force-field assignment for this project is **GLYCAM06j-1 for cellulose** (`leaprc.GLYCAM_06j-1`) and **GAFF2/AM1-BCC for ethanol**. The cellulose builder defaults to `glycam`, and the build command below explicitly selects it. The saved `cellulose_ethanol/cellulose_build_report.json` records GLYCAM and this release for the provided cellulose structure. The liquid parameterization script generates GAFF2 parameters; it does not parameterize cellulose.

Run commands from the `Cellulose_interface` project root.

## Files and tools

| File or directory | Purpose |
| --- | --- |
| `build_cellulose_slab.py` | Call an external CelluloseForge checkout to generate the slab and its topology |
| `parameterize_molecule.py` | Generate GAFF2/AM1-BCC small-molecule parameters using AmberTools and ParmEd |
| `pack_liquid_cellulose.py` | Use Packmol to place liquid on both sides of the cellulose slab |
| `combine_gromacs_topologies.py` | Combine topologies and convert packed coordinates to GROMACS format |
| `pyproject.toml` | Package metadata and optional installed command-line tools |
| `cellulose_ethanol/` | Existing ethanol/cellulose inputs, MDP templates, and pair-validation artifacts |
| `cellulose_ethanol/run_restrained.sbatch` | Cluster submission script; all ethanol simulation inputs and outputs live directly in `cellulose_ethanol/` |

## Requirements

- Python 3.11 or newer.
- A separate CelluloseForge checkout and `uv` for building cellulose. The default checkout location is `../cellulose_functionalization`; supply `--cellulose-forge-directory` if it is elsewhere. That default checkout is not included in this folder.
- AmberTools: `antechamber`, `parmchk2`, `tleap`, and the GLYCAM parameter files for the cellulose example.
- ParmEd importable by the Python interpreter used for parameterization.
- Open Babel (`obabel`) and Packmol (`packmol`) on PATH. Packmol must support the `pbc` input keyword.
- GROMACS for preprocessing and simulation; Slurm for the supplied batch script.

Activate your AmberTools environment before parameterization. Installing this project's Python package does not install the external simulation tools or force fields.

You can run the Python files directly. Optionally install the command-line entry points:

```bash
python -m pip install -e .
```

The installed names are `build-cellulose-slab`, `parameterize-molecule`, `pack-liquid-cellulose`, and `combine-gromacs-topologies`.

## 1. Build the cellulose slab

Example: four layers, approximately 40 Å wide and 40 Å long, with periodic cellulose chains and restraints on layers 2 and 3.

```bash
mkdir -p cellulose_ethanol
python build_cellulose_slab.py \
  --cellulose-forge-directory /path/to/cellulose_functionalization \
  --allomorph I-beta \
  --face 1 1 0 \
  --layers 4 \
  --width 40 \
  --length 40 \
  --chain-boundary periodic \
  --vacuum 0 \
  --force-field glycam \
  --restrain-layers 2 3 \
  --restraint-force 1000 1000 1000 \
  --output cellulose_ethanol/slab_build
```

Replace the checkout path with your installation. If needed, also pass `--uv /path/to/uv`, `--tleap-executable /path/to/tleap`, or `--parmed-python /path/to/python`.

The wrapper forwards `--output` to CelluloseForge. Confirm the output-path convention for your installed CelluloseForge version with `uv run --project /path/to/cellulose_functionalization cellulose-forge build --help`. The saved build report in this project records outputs named `structure.pdb`, `structure.gro`, and `system.top`; use the actual generated GRO and TOP paths in the next steps.

The commands below assume these are `cellulose_ethanol/slab_build/structure.gro` and `cellulose_ethanol/slab_build/system.top`. This intermediate slab-build directory keeps CelluloseForge outputs separate from the final combined `cellulose_ethanol/system.top`. Cellulose dimensions are rounded to fit the crystal lattice. Review the generated build report for realized dimensions, chain ordering, periodic connectivity, and restrained layers. The existing example has 28 chains, eight glucose units per chain, and 4,704 cellulose atoms.

`--dry-run` prints the CelluloseForge command, but still requires a valid checkout and `uv` path.

## 2. Parameterize the liquid

For neutral ethanol, use SMILES `CCO` and residue name `ETO`:

```bash
python parameterize_molecule.py \
  --smiles 'CCO' \
  --residue-name ETO \
  --net-charge 0 \
  --prefix ethanol \
  --output-directory cellulose_ethanol
```

This generates `ethanol.pdb`, `.mol2`, `.frcmod`, `.prmtop`, `.inpcrd`, `.top`, `.gro`, `.itp`, and `ethanol_atomtypes.itp` in `cellulose_ethanol/`.

Use `ethanol.top` as the secondary topology when combining systems. The `.itp` molecule file separately encodes explicit type-2 1–4 interactions. The script uses GAFF2 atom types and AM1-BCC charges.

## 3. Pack liquid around cellulose

```bash
python pack_liquid_cellulose.py \
  --cellulose-coordinates cellulose_ethanol/slab_build/structure.gro \
  --smiles 'CCO' \
  --residue-name ETO \
  --molecules 1000 \
  --density 0.789 \
  --seed 1234567 \
  --output cellulose_ethanol/interface.pdb
```

Density is in **g/cm³**; `0.789` is an example packing target. Width, length, and Packmol tolerance arguments are in **Å**, while GROMACS box dimensions printed by the packing script are in **nm**.

The script retains the cellulose X/Y box lengths, centers the slab along Z, and places half the liquid molecules on each surface. The molecule count must therefore be even. Use the same SMILES and residue name as in parameterization so atom naming and molecule identification remain consistent.

The initial Z length is calculated from the liquid mass, target density, lateral area, and slab thickness. Failed packing attempts increase the total Z dimension; the resulting initial liquid density can therefore be lower than the requested target. Packing density is an initial construction parameter, not a measurement of equilibrated bulk density.

Useful options include `--tolerance-angstrom 2.0`, `--max-attempts 5`, and `--z-growth-factor 1.05`. `--dry-run` prints the initial box calculation without packing, but still requires Open Babel to determine molecular mass.

## 4. Assemble the combined GROMACS system

```bash
python combine_gromacs_topologies.py \
  --primary cellulose_ethanol/slab_build/system.top \
  --secondary cellulose_ethanol/ethanol.top \
  --combined-coordinates cellulose_ethanol/interface.pdb \
  --output-directory cellulose_ethanol \
  --system-name 'Cellulose I-beta ethanol interface'
```

The resulting directory contains:

```text
cellulose_ethanol/
├── interface.gro
├── system.top
├── forcefield.itp
├── primary.itp
└── secondary.itp
```

Keep these five files together. The packed coordinates must list cellulose atoms first, followed by liquid molecules, with a liquid residue name matching the secondary topology. The packing script produces that ordering. Only orthorhombic boxes are supported by this workflow.

Use the GLYCAM cellulose topology as `--primary`. The combiner retains its global 1–4 scaling (`fudgeLJ = 1`, `fudgeQQ = 1`) and writes the GAFF2 secondary molecule's scaling explicitly into its pairs. This supports the GLYCAM/GAFF2 combination used here. It requires matching nonbonded functions and combination rules, and rejects conflicting atom-type Lennard-Jones parameters. Inspect any warnings from GROMACS preprocessing rather than bypassing them automatically.

The provided example already has these five files directly in `cellulose_ethanol/`. Rebuilding with the commands above replaces the final interface inputs there; run construction only when no simulation is using those files.

## 5. Run restrained minimization, NPT, and NVT

The batch script currently follows this sequence:

| Stage | Settings |
| --- | --- |
| Minimization | Steepest descent, maximum 50,000 steps, force tolerance 500 kJ/mol/nm |
| NPT | 500 ps, 2 fs timestep, 298.15 K, semi-isotropic C-rescale pressure coupling at 1 bar; X/Y fixed and Z free to relax |
| NVT production | 2 ns, 2 fs timestep, 298.15 K; final NPT box and velocities retained |

All stages use `periodic-molecules = yes` and `define = -DPOSRES_CELLULOSE`, which activates the middle-layer heavy-atom restraints in this example. A generic `-DPOSRES` does not activate this topology's restraints.

Fixed lateral area during NPT is implemented with:

```ini
pcoupltype      = semiisotropic
compressibility = 0.0 4.5e-5
refcoord-scaling = all
```

Semi-isotropic coupling alone allows X and Y to change together; the zero first compressibility fixes them. NVT uses the final instantaneous NPT box. Its restraint reference is reset to the final NPT cellulose coordinates.

### Submit the ethanol system

All ethanol inputs, MDP files, submission script, and simulation outputs live directly in `cellulose_ethanol/`. Shared Python tools and package metadata stay in the project root.

```bash
cd cellulose_ethanol
sbatch run_restrained.sbatch
```

The script uses the `interface.gro`, `system.top`, included ITP files, and MDP files in this directory. It runs minimization, 500 ps fixed-X/Y NPT, and 2 ns restrained NVT in sequence. No separate run directory is created. To change the production duration before a fresh run, edit `nsteps` in `nvt_production.mdp`: at 2 fs, 1,000,000 steps is 2 ns.

The construction commands above write the five final interface files directly into `cellulose_ethanol/`, ready for submission. Do not overwrite inputs while a job is running. Existing TPRs and checkpoints are reused on resubmission; changed settings require fresh stage outputs rather than reusing old TPRs.

The cluster script uses `flight`, the `gridware` environment, module `apps/gromacs_cuda/2026.1`, and `gmx_mpi`, with one MPI process and 32 OpenMP threads. Its Slurm `--chdir` points to this checkout's `cellulose_ethanol/` directory; update this and the partition/resources when moving the project to another location or cluster.

Resubmit to resume available checkpoints. Completed stages are skipped when their final `.gro` file exists. Do not submit concurrent jobs into this directory.

## 6. Check progress and equilibration

```bash
squeue -u "$USER"
tail -n 40 cellulose_ethanol/npt.log
tail -n 40 cellulose_ethanol/nvt_production.log
```

Each stage writes its own log and energy file; dynamics also writes compressed trajectories and checkpoints. A final `.gro` indicates that the stage returned its final coordinates. Inspect minimization's actual convergence message, not only file existence.

Use `gmx energy` (or `gmx_mpi energy`) to inspect temperature, potential energy, pressure, and box lengths from `npt.edr` and `nvt_production.edr`. X and Y should remain fixed throughout this NPT protocol. Check Z relaxation, cellulose structure, and the ethanol mass-density profile along Z before judging equilibration. The current batch script advances automatically after 500 ps; it does not decide whether equilibration has been achieved.

## Other liquids and troubleshooting

- Change SMILES, residue name, molecular charge, molecule count, and target packing density consistently for another liquid. This parameterizer specifically uses GAFF2/AM1-BCC; assess whether that model is suitable for the molecule and study.
- For missing executables, activate the relevant environment or supply the scripts' executable-path options. Use `python SCRIPT.py --help` for all supported flags.
- If ParmEd cannot be imported, run parameterization with the Python interpreter in your AmberTools/ParmEd environment.
- If the cellulose builder cannot find CelluloseForge, supply its checkout path explicitly.
- If packing fails, review liquid count, target density, tolerance, and retry settings.
- If restraints are inactive, confirm that the topology contains `#ifdef POSRES_CELLULOSE` and that the processed MDP enables it.
- If GROMACS reports an atom-order/count mismatch, check the coordinate/topology pair and consistent liquid atom/residue naming.
- `cellulose_ethanol/ethanol_pair_validation/` contains a saved isolated ethanol 1–4 energy comparison. Its saved report passed; it is not a full-system equilibration test. The validation script currently expects an `ethanol_parameters/` subdirectory under `cellulose_ethanol/`, whereas the provided ethanol parameter files are directly in `cellulose_ethanol/`; update that source path before rerunning it.
