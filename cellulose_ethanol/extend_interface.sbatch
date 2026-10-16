#!/usr/bin/env bash
#SBATCH --job-name=extend_interface
#SBATCH --chdir=/mnt/scratch/users/ass2009/Cellulose_interface/cellulose_ethanol
#SBATCH --partition=nodes
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=12:00:00
#SBATCH --mem=16G
#SBATCH --output=extend_interface.%j.out
#SBATCH --error=extend_interface.%j.err

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
# Convert the original TPR to retain its exact parameters and restraint reference.
if [[ ! -f npt_extended.done ]]; then
    if [[ ! -s npt_extended.tpr ]]; then
        "$GMX" convert-tpr -s npt.tpr -extend 2000 -o npt_extended.tpr
    fi
    checkpoint=npt.cpt
    [[ ! -s npt_extended.cpt ]] || checkpoint=npt_extended.cpt
    "$GMX" mdrun -s npt_extended.tpr -deffnm npt_extended \
        -cpi "$checkpoint" -noappend -ntomp "$OMP_NUM_THREADS" -pin on -cpt 5
    touch npt_extended.done
fi

# -noappend preserves old results and writes part-numbered coordinates/energies.
final_npt_gro=$(python3 - <<'COORDS'
from pathlib import Path
files = list(Path('.').glob('npt_extended*.gro'))
if not files:
    raise SystemExit('No extended NPT final coordinates found')
print(max(files, key=lambda p: p.stat().st_mtime))
COORDS
)
if [[ ! -s nvt_production_extended.gro ]]; then
    if [[ ! -s nvt_production_extended.tpr ]]; then
        "$GMX" grompp -f nvt_production_extended.mdp -c "$final_npt_gro" \
            -t npt_extended.cpt -r "$final_npt_gro" -p system.top \
            -o nvt_production_extended.tpr -po nvt_production_extended.processed.mdp
    fi
    restart=()
    [[ ! -s nvt_production_extended.cpt ]] || restart=(-cpi nvt_production_extended.cpt -append)
    "$GMX" mdrun -deffnm nvt_production_extended -ntomp "$OMP_NUM_THREADS" \
        -pin on -cpt 5 "${restart[@]}"
fi
echo 'Finished: interface NPT extended to 2.5 ns, then fresh 2 ns restrained NVT.'
