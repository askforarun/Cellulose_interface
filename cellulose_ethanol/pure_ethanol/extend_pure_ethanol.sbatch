#!/usr/bin/env bash
#SBATCH --job-name=extend_pure_ethanol
#SBATCH --chdir=/mnt/scratch/users/ass2009/Cellulose_interface/cellulose_ethanol/pure_ethanol
#SBATCH --partition=nodes
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=12:00:00
#SBATCH --mem=16G
#SBATCH --output=extend_pure_ethanol.%j.out
#SBATCH --error=extend_pure_ethanol.%j.err

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
        "$GMX" convert-tpr -s npt_sample.tpr -extend 5000 -o npt_extended.tpr
    fi
    checkpoint=npt_sample.cpt
    [[ ! -s npt_extended.cpt ]] || checkpoint=npt_extended.cpt
    "$GMX" mdrun -s npt_extended.tpr -deffnm npt_extended \
        -cpi "$checkpoint" -noappend -ntomp "$OMP_NUM_THREADS" -pin on -cpt 5
    touch npt_extended.done
fi

# Join any continuation parts before analysing the final 1 ns (7–8 ns).
parts=(npt_extended*.edr)
"$GMX" eneconv -f "${parts[@]}" -o npt_extended_combined.edr
printf 'Density\n0\n' | "$GMX" energy -f npt_extended_combined.edr \
    -b 7000 -e 8000 -o extended_last1ns_density.xvg > extended_last1ns_density_summary.txt 2>&1
printf 'Density\nTemperature\nPressure\nVolume\n0\n' | "$GMX" energy \
    -f npt_extended_combined.edr -o extended_observables.xvg > extended_summary.txt 2>&1
echo 'Finished: pure ethanol NPT extended to 8 ns total. Final 1 ns density:'
cat extended_last1ns_density_summary.txt
