#!/usr/bin/env bash
#SBATCH --job-name=interface_density
#SBATCH --chdir=/mnt/scratch/users/ass2009/Cellulose_interface/cellulose_ethanol
#SBATCH --partition=nodes
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:10:00
#SBATCH --mem=4G
#SBATCH --output=interface_density.%j.out
#SBATCH --error=interface_density.%j.err
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
printf 'Cellulose\nLiquid\nCellulose\n' | "$GMX" density \
    -f nvt_production.xtc -s nvt_production.tpr -n density_index.ndx \
    -o density_z_last1ns.xvg -b 1000 -e 2000 -d Z -dens mass -center -ng 2 -sl 200
