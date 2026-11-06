# Cellulose–liquid interface workflow

Build a crystalline cellulose slab with **GLYCAM**, parameterize a small liquid molecule, pack liquid on both surfaces, and prepare a combined GROMACS system. The example below uses cellulose I-beta (110) and ethanol.

The force-field assignment for this project is **GLYCAM06j-1 for cellulose** (`leaprc.GLYCAM_06j-1`) and **GAFF2/AM1-BCC for ethanol**. The cellulose builder defaults to `glycam`, and the build command below explicitly selects it. The saved `cellulose_ethanol/cellulose_build_report.json` records GLYCAM for the provided cellulose structure. It was produced by the CelluloseForge development build that became 1.1.0; because it predates the version bump it reports `1.0.0`, but it already uses the 1.1.0 slab frame (surface normal along `-z`) and heavy-atom-only restraint selection. The liquid parameterization script generates GAFF2 parameters; it does not parameterize cellulose.

Run commands from the `Cellulose_interface` project root.

## Files and tools

| File or directory | Purpose |
| --- | --- |
| `build_cellulose_slab.py` | Call the installed `cellulose-forge` command (CelluloseForge 1.1.0) to generate the slab and its topology |
| `parameterize_molecule.py` | Generate GAFF2/AM1-BCC small-molecule parameters using AmberTools and ParmEd |
| `pack_liquid_cellulose.py` | Use Packmol to place liquid on both sides of the cellulose slab |
| `combine_gromacs_topologies.py` | Combine topologies and convert packed coordinates to GROMACS format |
| `pyproject.toml` | Package metadata and optional installed command-line tools |
| `cellulose_ethanol/` | Existing ethanol/cellulose inputs, MDP templates, and pair-validation artifacts |
| `cellulose_ethanol/run_restrained.sbatch` | Cluster submission script; all ethanol simulation inputs and outputs live directly in `cellulose_ethanol/` |

## Requirements

- Python 3.11 or newer.
- [CelluloseForge 1.1.0](https://github.com/askforarun/cellulose-forge/releases/tag/v1.1.0) for building cellulose, providing the `cellulose-forge` command (see below).
- AmberTools: `antechamber`, `parmchk2`, `tleap`, and the GLYCAM parameter files for the cellulose example.
- ParmEd importable by the Python interpreter used for parameterization.
- Open Babel (`obabel`) and Packmol (`packmol`) on PATH. Packmol must support the `pbc` input keyword.
- GROMACS for preprocessing and simulation; Slurm for the supplied batch script.

### Recommended Conda installation

The repository includes a complete `environment.yml` containing AmberTools,
ParmEd, Open Babel, periodic-capable Packmol, GROMACS, analysis libraries,
CelluloseForge 1.0.0, and this project's command-line tools:

```bash
conda env create --file environment.yml
conda activate cellulose-interface
```

To use a different environment name without editing the file:

```bash
conda env create --name MY_ENVIRONMENT --file environment.yml
conda activate MY_ENVIRONMENT
```

Verify the required commands after installation:

```bash
python --version
cellulose-forge --version
antechamber -h >/dev/null
parmchk2 -h >/dev/null
tleap -h >/dev/null
obabel -V
packmol --version
gmx --version
```

The supplied Slurm scripts use a cluster-provided `gmx_mpi` module and Slurm;
those site-specific tools are intentionally not installed by this local Conda
environment.

Packmol is installed from its official PyPI distribution because the
conda-forge AmberTools package intentionally conflicts with conda-forge's
Packmol package. The environment includes a Fortran build toolchain as a
fallback when PyPI does not provide a prebuilt Packmol wheel for the platform.

Activate your AmberTools environment before parameterization. Installing this project's Python package does not install CelluloseForge, the external simulation tools, or force fields.

### Install CelluloseForge 1.1.0

CelluloseForge is not yet on PyPI. Install the tagged 1.1.0 release from GitHub into the environment you use for this workflow:

```bash
python -m pip install "cellulose-forge @ git+https://github.com/askforarun/cellulose-forge.git@v1.1.0"
cellulose-forge --help
```

`build_cellulose_slab.py` runs `cellulose-forge` from your PATH. Pass `--cellulose-forge-executable /path/to/cellulose-forge` if it is installed elsewhere.

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

If needed, also pass `--cellulose-forge-executable /path/to/cellulose-forge`, `--tleap-executable /path/to/tleap`, or `--parmed-python /path/to/python`.

The wrapper forwards `--output` to CelluloseForge. Confirm the output-path convention with `cellulose-forge build --help`. The saved build report in this project records outputs named `structure.pdb`, `structure.gro`, and `system.top`; use the actual generated GRO and TOP paths in the next steps.

The commands below assume these are `cellulose_ethanol/slab_build/structure.gro` and `cellulose_ethanol/slab_build/system.top`. This intermediate slab-build directory keeps CelluloseForge outputs separate from the final combined `cellulose_ethanol/system.top`. Cellulose dimensions are rounded to fit the crystal lattice. Review the generated build report for realized dimensions, chain ordering, periodic connectivity, and restrained layers. The existing example has 28 chains, eight glucose units per chain, and 4,704 cellulose atoms.

`--dry-run` prints the CelluloseForge command, but still requires `cellulose-forge` to be installed.

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

### Encoding each force field’s 1–4 scaling in the combined topology

A 1–4 pair consists of atoms separated by three bonds. Its Lennard-Jones and electrostatic interactions can use force-field-specific scaling. Simply combining two molecule topologies under one set of GROMACS `[ defaults ]` can therefore change the liquid's intramolecular energy.

For this workflow, the scaling conventions are:

| Component | Lennard-Jones 1–4 multiplier | Electrostatic 1–4 multiplier |
| --- | --- | --- |
| GLYCAM cellulose | 1.0 | 1.0 |
| GAFF2 liquid generated by this parameterizer | 0.5 | approximately 0.83333333 |

GROMACS uses one set of global 1–4 scaling factors for the combined system, but GLYCAM and GAFF2 require different values. We keep the global factors at **1.0 for GLYCAM** and specify GAFF2’s scaling directly for each liquid 1–4 pair: **0.5 for Lennard-Jones** and **1/1.2 for electrostatics**. This lets both force fields retain their original 1–4 interactions within the same simulation.

`combine_gromacs_topologies.py` handles this automatically:

1. Retain the primary topology's global defaults, including GLYCAM's unit scaling.
2. Read the secondary topology's native scaling factors, atom charges, and atom-type Lennard-Jones parameters.
3. Convert each secondary type-1 `[ pairs ]` entry into an explicit type-2 entry with fields `ai aj funct fudgeQQ qi qj sigma14 epsilon14`.
4. Encode electrostatic scaling in the pair's `fudgeQQ` field and Lennard-Jones scaling in `epsilon14`.

For the supported combination rule 2, the converter uses:

```text
sigma14   = (sigma_i + sigma_j) / 2
epsilon14 = sqrt(epsilon_i * epsilon_j) * secondary_fudgeLJ
pair_fudgeQQ = secondary_fudgeQQ
```

Although the global scaling factors are set to 1.0 for GLYCAM, the liquid retains its original force-field-specific 1–4 scaling because those factors are specified directly for each liquid 1–4 pair. These per-pair scaling settings apply to intramolecular 1–4 interactions; it does not introduce a separate scaling factor for cellulose–liquid interactions. `parameterize_molecule.py` also writes an explicit-pair molecule ITP, but the combiner should receive the native secondary `.top` so that it can read its original defaults and type-1 pairs. Do not apply the conversion twice.

The current converter supports combination rule 2 and ordinary three-field type-1 secondary pair rows. It rejects other pair formats rather than assuming their meaning. Review the source force-field conventions when introducing another liquid or topology source.

### Checking that the original 1–4 energies are preserved

Validate the conversion with an isolated molecule before relying on a new force-field combination. This method applies to any supported liquid:

1. Prepare a native topology with its original global scaling and type-1 pairs.
2. Prepare an equivalent topology with the host's global defaults and explicit type-2 liquid pairs.
3. Use identical coordinates, periodic box, nonbonded settings, and GROMACS version for both single-point evaluations.
4. Compare bond, angle, proper-dihedral, LJ-14, Coulomb-14, and total potential energies. They should agree within a specified numerical tolerance.
5. Retain the settings and a machine-readable report of the energies, differences, tolerance, and pass/fail result.

The worked example in this repository is a **cellulose–ethanol interface**, using GLYCAM for cellulose and GAFF2/AM1-BCC for ethanol. The accompanying test in [`cellulose_ethanol/ethanol_pair_validation/`](cellulose_ethanol/ethanol_pair_validation/) checks that the ethanol 1–4 interactions retain their original GAFF2 energies when represented under the GLYCAM global defaults used in the combined interface topology. It compares an isolated ethanol molecule using native GAFF2 type-1 pairs with the same molecule using explicit type-2 pairs. The saved `validation_report.json` reports matching energies within **1e-5 kJ/mol** for all compared terms. This validates the preservation of the liquid’s 1–4 scaling for this example; repeat the comparison for new liquid parameter sets. The isolated-molecule test does not establish bulk-liquid accuracy or equilibration of the cellulose–ethanol interface.

The example validation script currently looks for source parameters in `cellulose_ethanol/ethanol_parameters/`, while this checkout stores them directly in `cellulose_ethanol/`. Before rerunning the example, change its `SOURCE` assignment to `SOURCE = ROOT`. Adapting it to another liquid also requires updating molecule/file names and any hard-coded native scaling factors.

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

## Pure-ethanol reference density

A separate bulk-liquid calculation lives in `cellulose_ethanol/pure_ethanol/`. It contains 1,000 ethanol molecules (9,000 atoms), with the same GAFF2/AM1-BCC parameters and charge-rounding correction as the interface liquid. Its standalone topology uses native GAFF2 type-1 pairs and global scaling of 0.5 for LJ-14 and 0.83333333 for Coulomb-14. There is no cellulose and no position restraint.

```bash
cd cellulose_ethanol/pure_ethanol
sbatch run_pure_ethanol.sbatch
```

The job inserts the molecules into a 5.1 nm cubic starting box, verifies the atom count, minimizes, runs 1 ns of initial NPT equilibration, and then runs a separate 2 ns NPT sampling stage. Both NPT stages use isotropic C-rescale coupling at 1 bar and V-rescale temperature coupling at 298.15 K, with a 2 fs timestep and hydrogen-bond constraints. All three box lengths can relax together; fixed lateral dimensions are specific to the cellulose interface, not this bulk reference.

At completion, `density.xvg` records the sampling-stage density versus time. `density_summary.txt` contains the GROMACS average, block-based error estimate, RMSD, and drift for that stage. Density is reported in kg/m³; divide by 1,000 to obtain g/cm³. `equilibration_observables.xvg` and `sampling_observables.xvg`, with their corresponding summary files, provide density, temperature, pressure, and volume for the two stages. The sampling-stage time axis starts at 1,000 ps and ends at 3,000 ps. These analysis files are generated by the job and are not tracked in the repository.

These durations are initial choices, not an automatic convergence guarantee. Inspect density and volume traces and compare averages over successive time blocks. If they are still drifting, extend equilibration and exclude the transient period from the density estimate. Compare a stable bulk reference with the ethanol density plateau away from cellulose surfaces, using matching temperature, force field, and nonbonded settings. The model's density should be measured rather than forced to equal the initial packing target.

## Extending equilibration and restarting production

Keep the existing results and continue from saved NPT checkpoints:

```bash
cd cellulose_ethanol
sbatch extend_interface.sbatch
cd pure_ethanol
sbatch extend_pure_ethanol.sbatch
```

- **Interface:** `extend_interface.sbatch` extends the original restrained NPT by 2 ns, giving 2.5 ns total NPT, with X/Y fixed and Z free to relax. It then runs a fresh 2 ns restrained NVT segment, with time axis 2.5–4.5 ns, from the extended NPT coordinates and velocities. The NVT restraint reference is reset to the final NPT geometry.
- **Pure liquid:** `extend_pure_ethanol.sbatch` adds 5 ns to the original NPT calculation, giving 8 ns total including initial equilibration. It extracts the density over 7–8 ns to `extended_last1ns_density.xvg` and `extended_last1ns_density_summary.txt`, and writes density/temperature/pressure/volume for the continuation to `extended_observables.xvg`.

Both scripts use `gmx convert-tpr -extend` and checkpoint continuation to preserve the original NPT parameters. The interface continuation therefore retains the original scaled restraint reference, rather than replacing it at the start of the NPT extension. New NPT outputs use the `npt_extended` prefix and GROMACS part numbers; the original outputs remain intact. New interface NVT outputs use `nvt_production_extended`. Do not run multiple continuation jobs in the same directory at once.

The extra durations do not guarantee equilibration. Review late-time density/volume for the pure liquid and Z/Pzz for the interface, along with profiles in successive time windows. The interface script advances automatically to NVT; treat that segment as production only after these checks support equilibration.

## Larger liquid reservoir: 1,500 ethanol molecules per surface

The larger-reservoir comparison contains **3,000 ethanol molecules total**, with 1,500 initially packed on each cellulose surface. It retains the original GLYCAM cellulose slab and the same GAFF2/AM1-BCC ethanol parameters, including the per-pair 1–4 scaling representation. The generated primary and secondary molecule ITPs were checked to be byte-identical to those of the original interface.

All case files remain directly in `cellulose_ethanol/` and use `_3000` names, preserving the previous simulations. The initial box is approximately **4.1751 × 4.1520 × 19.5250 nm**, with an initial geometric liquid thickness of **8.69 nm per face**. Packing used the simulated pure-liquid density of **0.761501 g/cm³** (7–8 ns average of the extended pure-ethanol NPT run) as an initial density target, not a constraint on the final NPT density. `interface_3000_build_report.json` records the construction checks, including 31,704 total atoms and 1,500 liquid molecules per face.

Packing command from the project root (activate an environment with Open Babel and Packmol, or provide their executable paths):

```bash
python pack_liquid_cellulose.py \
  --cellulose-coordinates cellulose_ethanol/cellulose.gro \
  --smiles CCO --residue-name ETO \
  --molecules 3000 --density 0.761501 \
  --seed 20261016 --maxit 200 --nloop 20 \
  --output cellulose_ethanol/interface_3000.pdb
```

The assembled starting files are `interface_3000.gro`, `system_3000.top`, `forcefield_3000.itp`, `primary_3000.itp`, and `secondary_3000.itp`. The combined topology uses the same cellulose molecule counts as before and `ETO 3000`; its include paths refer to the case-specific ITP filenames.

```bash
cd cellulose_ethanol
sbatch run_interface_3000.sbatch
```

This job runs restrained steepest-descent minimization followed by **2.5 ns NPT** at 298.15 K and 1 bar. It uses C-rescale semi-isotropic pressure coupling with `compressibility = 0.0 4.5e-5`, keeping X/Y fixed while Z relaxes. Middle-layer heavy-atom position restraints remain enabled with `POSRES_CELLULOSE`. Output prefixes are `em_3000` and `npt_3000`. This job stops after NPT; it does not automatically launch NVT. Compare late-time liquid density profiles with the smaller reservoir and pure-liquid reference before deciding the next stage.

## Other liquids and troubleshooting

- Change SMILES, residue name, molecular charge, molecule count, and target packing density consistently for another liquid. This parameterizer specifically uses GAFF2/AM1-BCC; assess whether that model is suitable for the molecule and study.
- For missing executables, activate the relevant environment or supply the scripts' executable-path options. Use `python SCRIPT.py --help` for all supported flags.
- If ParmEd cannot be imported, run parameterization with the Python interpreter in your AmberTools/ParmEd environment.
- If the cellulose builder cannot find CelluloseForge, install the 1.1.0 release as above, or pass `--cellulose-forge-executable` with the path to `cellulose-forge`.
- If packing fails, review liquid count, target density, tolerance, and retry settings.
- If restraints are inactive, confirm that the topology contains `#ifdef POSRES_CELLULOSE` and that the processed MDP enables it.
- If GROMACS reports an atom-order/count mismatch, check the coordinate/topology pair and consistent liquid atom/residue naming.

## Licence

This project is licensed under the GNU General Public License v3.0 or later (GPL-3.0-or-later), matching CelluloseForge. See [`LICENSE`](LICENSE).
