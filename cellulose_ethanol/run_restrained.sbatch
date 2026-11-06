#!/usr/bin/env bash
#SBATCH --job-name=cellulose_ethanol
#SBATCH --chdir=/users/ass2009/sharedscratch/Cellulose_interface/cellulose_ethanol
#SBATCH --partition=nodes
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=48:00:00
#SBATCH --mem=32G
#SBATCH --output=cellulose_restrained.%j.out
#SBATCH --error=cellulose_restrained.%j.err

# Submit from cellulose_ethanol/: sbatch run_restrained.sbatch
# All inputs and outputs live directly in cellulose_ethanol/.
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
for input in interface.gro system.top forcefield.itp slab.itp liquid.itp em.mdp npt.mdp nvt_production.mdp; do
    [[ -s "$input" ]] || { echo "Missing input: $PWD/$input; submit from cellulose_ethanol/" >&2; exit 1; }
done

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

# All stages activate POSRES_CELLULOSE: heavy atoms of middle layers 2 and 3.
run_stage em -c interface.gro -r interface.gro
run_stage npt -c em.gro -r em.gro
# Reset restraint reference to the pressure-relaxed cellulose geometry.
# Carry NPT velocities via the checkpoint; keep the final NPT box fixed.
run_stage nvt_production -c npt.gro -t npt.cpt -r npt.gro
echo "Finished restrained EM -> 500 ps fixed-X/Y NPT -> 2 ns restrained NVT (default)."
