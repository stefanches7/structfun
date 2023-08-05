#!/bin/bash
# Job Name and Files (also --job-name)
#SBATCH -J nest_izhik
#Output and error (also --output, --error):
#SBATCH -o ./%x.%j.out
#SBATCH -e ./%x.%j.err
#Initial working directory (also --chdir):
#SBATCH -D ./
#SBATCH --time=24:00:00
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --clusters=cm2
#SBATCH --partition=cm2_std
#SBATCH --qos=cm2_std
#SBATCH --account=ge72puf2

#SBATCH --nodes=8
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2

#SBATCH --time=01:00:00
#SBATCH --hint=nomultithread


#Important
module load slurm_setup

#Run the program:
export OMP_NUM_THREADS=2
export OMP_PROC_BIND=TRUE

mpiexec -n 8 -bootstrap=ssh python3 nest_simulation.py

