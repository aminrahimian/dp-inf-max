#!/bin/bash
#SBATCH --job-name="dp-alg"
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=40
#SBATCH --mem-per-cpu=3G
#SBATCH --cluster=mpi
#SBATCH--partition=mpi
#SBATCH --time=04:00:00


#SBATCH -o /ihome/arahimian/cah259/dp_inf_max/output.out
#SBATCH -e /ihome/arahimian/cah259/dp_inf_max/errors.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cah259@pitt.edu

module load gcc/8.2.0 python/anaconda3.9-2021.11
python3 /ihome/arahimian/cah259/dp_inf_max/computing_mean_spread_size.py ${SLURM_ARRAY_TASK_ID}

