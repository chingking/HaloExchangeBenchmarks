#!/bin/bash
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=2
#SBATCH --constraint=k20c
#SBATCH --cpus-per-task=1
#SBATCH --partition=normal
#SBATCH --ntasks-per-core=1
#SBATCH --mem=24000
#SBATCH --time=01:00:00

export OMP_NUM_THREADS=1
export MV2_USE_CUDA=1
mpiexec.hydra -np 8 -ppn 2 ./CommTest.exe --nprocx 1 --nprocy 8 --ie 520 --je 350 --enable-GCL --ntracer-perHandler 3 --nGCLHandlers 30 --lperi_y --disable-validation-report --nbl_exchg 2
################################################################
