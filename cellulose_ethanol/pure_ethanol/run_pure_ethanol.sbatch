#!/usr/bin/env bash
#SBATCH --job-name=pure_ethanol
#SBATCH --chdir=/users/ass2009/sharedscratch/Cellulose_interface/cellulose_ethanol/pure_ethanol
#SBATCH --partition=nodes
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=12:00:00
#SBATCH --mem=16G
#SBATCH --output=pure_ethanol.%j.out
#SBATCH --error=pure_ethanol.%j.err

set -eo pipefail
flight start
# Flight activation can fail transiently while generating its shell environment.
activated=false
for attempt in 1 2 3; do
    if flight env activate gridware; then
        activated=true
        break
    fi
    echo "Flight activation failed (attempt $attempt/3)" >&2
    sleep 5
done
[[ "$activated" == true ]] || { echo "Cannot activate Gridware" >&2; exit 1; }
module load apps/gromacs_cuda/2026.1
set -u

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-32}"
GMX="${GMX:-gmx_mpi}"
command -v "$GMX" >/dev/null
# Start dilute enough to insert 1000 molecules; NPT determines the final density.
if [[ ! -s initial.gro ]]; then
    "$GMX" insert-molecules -ci ethanol.gro -nmol 1000 -box 5.1 5.1 5.1 \
        -try 1000 -seed 20261016 -o initial.gro
fi
# GROMACS may return success even if fewer molecules were inserted.
python3 - <<'CHECK'
from pathlib import Path
lines = Path('initial.gro').read_text().splitlines()
if int(lines[1]) != 9000:
    raise SystemExit('Expected exactly 1000 nine-atom ethanol molecules (9000 atoms).')
CHECK

run_stage() {
    local stage="$1"
    shift
    if [[ -s "$stage.gro" ]]; then
        echo "Completed stage: $stage"
        return
    fi
    if [[ ! -s "$stage.tpr" ]]; then
        "$GMX" grompp -f "$stage.mdp" -p system.top -o "$stage.tpr" \
            -po "$stage.processed.mdp" "$@"
    fi
    local restart=()
    if [[ -s "$stage.cpt" ]]; then
        restart=(-cpi "$stage.cpt" -append)
    fi
    "$GMX" mdrun -deffnm "$stage" -ntomp "$OMP_NUM_THREADS" -pin on \
        -cpt 5 "${restart[@]}"
    [[ -s "$stage.gro" ]] || { echo "Stage $stage did not finish" >&2; exit 1; }
}

run_stage em -c initial.gro
run_stage npt_equil -c em.gro
run_stage npt_sample -c npt_equil.gro -t npt_equil.cpt
# This summary covers only the separate 2 ns sampling stage.
printf 'Density\n0\n' | "$GMX" energy -f npt_sample.edr -o density.xvg > density_summary.txt 2>&1
printf 'Density\nTemperature\nPressure\nVolume\n0\n' | "$GMX" energy \
    -f npt_equil.edr -o equilibration_observables.xvg > equilibration_summary.txt 2>&1
printf 'Density\nTemperature\nPressure\nVolume\n0\n' | "$GMX" energy \
    -f npt_sample.edr -o sampling_observables.xvg > sampling_summary.txt 2>&1
echo 'Completed pure ethanol NPT. Density (kg/m^3):'
cat density_summary.txt
