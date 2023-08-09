#!/bin/bash
#SBATCH -J nest_test
#SBATCH -o nest_test_lrz_%x.%j.out
#SBATCH -D .
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --mail-type=none
#SBATCH -M serial
#SBATCH -p serial_std
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=14
#SBATCH --time=02:00:00

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
echo "NTASKS=$SLURM_NTASKS    THREADS_PER_TASK=$OMP_NUM_THREADS"
module restore nest_env

module load slurm_setup

mpiexec python3 nest_simulation.py
