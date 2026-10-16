#!/usr/bin/env bash
#SBATCH --job-name=extended_density
#SBATCH --chdir=/mnt/scratch/users/ass2009/Cellulose_interface/cellulose_ethanol
#SBATCH --partition=nodes
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:10:00
#SBATCH --mem=4G
#SBATCH --output=extended_density.%j.out
#SBATCH --error=extended_density.%j.err
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
for window in 2500:2999.99 3000:3499.99 3500:3999.99 4000:4500 3500:4500; do
    begin=${window%:*}
    end=${window#*:}
    tag=${end%%.*}
    printf 'Cellulose\nLiquid\nCellulose\n' | "$GMX" density \
        -f nvt_production_extended.xtc -s nvt_production_extended.tpr -n density_index.ndx \
        -o "extended_density_${begin}_${tag}.xvg" -b "$begin" -e "$end" \
        -d Z -dens mass -center -ng 2 -sl 200
done
