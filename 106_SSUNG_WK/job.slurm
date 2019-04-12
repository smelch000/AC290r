#!/bin/bash
#SBATCH -n 4                # Number of cores
#SBATCH -N 4                # Ensure that all cores are on one machine
#SBATCH -t 10          # Runtime in minutes
#SBATCH -p test             # Partition to submit to
#SBATCH --mem=100           # Memory pool for all cores (see also --mem-per-cpu)
###SBATCH -o job_%j.out
###SBATCH -e job_%j.err

#module load intel/15.0.0-fasrc01 openmpi/1.10.0-fasrc01
#module load MYPROGRAM

# module load gcc/8.2.0-fasrc01 openmpi/3.1.1-fasrc01  cuda/10.0.130-fasrc01
module load gcc/7.1.0-fasrc01 openmpi/3.1.1-fasrc01  cuda/10.0.130-fasrc01

srun -n $SLURM_NTASKS --mpi=pmi2 ./run2.py > out.txt 2> err.txt
