#!/bin/bash
#SBATCH -n 1                # Number of cores
#SBATCH -N 1                # Ensure that all cores are on one machine
#SBATCH -t 10          # Runtime in D-HH:MM, minimum of 10 minutes
###SBATCH -p gpu             # Partition to submit to
#SBATCH -p fas_gpu             # Partition to submit to
#SBATCH --gres=gpu:1        # gres=gpu:n where n=1-8. This 1 node has 24 cores and is equipped with 8 x Tesla K20Xm
#SBATCH --mem=8000           # Memory pool for all cores (see also --mem-per-cpu)
###SBATCH --mem=100           # Memory pool for all cores (see also --mem-per-cpu)
###SBATCH -o jobgpu_serial_%j.out
###SBATCH -e jobgpu_serial_%j.err

#module load intel/15.0.0-fasrc01 openmpi/1.10.0-fasrc01
#module load MYPROGRAM

export MOEBIUS_ROOT=~ac290ru1901/SM/MAGIC
export PYTHONPATH=$MOEBIUS_ROOT/BACKEND/SHOP:$MOEBIUS_ROOT/BACKEND/SCRIPTS

#module load gcc/8.2.0-fasrc01 openmpi/3.1.1-fasrc01  cuda/10.0.130-fasrc01
module load gcc/7.1.0-fasrc01 openmpi/3.1.1-fasrc01  cuda/10.0.130-fasrc01

srun ./run2.py -x gpu > out.txt 2> err.txt
