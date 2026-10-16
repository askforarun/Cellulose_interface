#!/usr/bin/env bash
#SBATCH --job-name=interface_3000
#SBATCH --chdir=/mnt/scratch/users/ass2009/Cellulose_interface/cellulose_ethanol
#SBATCH --partition=nodes
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=24:00:00
#SBATCH --mem=32G
#SBATCH --output=interface_3000.%j.out
#SBATCH --error=interface_3000.%j.err

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
# A new case: existing 1000-molecule inputs and trajectories remain unchanged.
for input in interface_3000.gro system_3000.top forcefield_3000.itp primary_3000.itp secondary_3000.itp; do
    [[ -s "$input" ]] || { echo "Missing input: $input" >&2; exit 1; }
done
python3 - <<'CHECK'
from pathlib import Path
assert int(Path('interface_3000.gro').read_text().splitlines()[1]) == 31704
CHECK

run_stage() {
    local stage="$1"
    shift
    if [[ -s "$stage.gro" ]]; then
        echo "Completed stage: $stage"
        return
    fi
    if [[ ! -s "$stage.tpr" ]]; then
        "$GMX" grompp -f "$stage.mdp" -p system_3000.top -o "$stage.tpr" \
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

run_stage em_3000 -c interface_3000.gro -r interface_3000.gro
run_stage npt_3000 -c em_3000.gro -r em_3000.gro
echo 'Finished 3000-ethanol case: restrained minimization and 2.5 ns fixed-X/Y NPT.'
