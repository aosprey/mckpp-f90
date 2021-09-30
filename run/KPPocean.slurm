#!/bin/bash
#
#SBATCH --job-name=KPP
#SBATCH --time=1:00:0
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=64
#SBATCH --account=n02-cms
#SBATCH --partition=standard
#SBATCH --qos=standard

module load epcc-job-env

export OMP_NUM_THREADS=64
export OMP_PLACES=cores
export CRAY_OMP_CHECK_AFFINITY=TRUE
export OMP_WAIT_POLICY=PASSIVE

srun --hint=nomultithread --distribution=block:block ./KPP_ocean
