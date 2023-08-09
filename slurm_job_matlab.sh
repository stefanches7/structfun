#!/bin/bash
#SBATCH -J matlab_izh
#SBATCH -o matlab_izh_%x.%j.out
#SBATCH -D .
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --mail-type=none
#SBATCH -M serial
#SBATCH -p serial_std
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=28
#SBATCH --time=08:00:00

module remove intel-mpi
module load intel-mpi/2018.4.274
module load matlab/R2022a_Update5-generic
module load slurm_setup

mpiexec matlab -nodesktop -r modIzhikevich
