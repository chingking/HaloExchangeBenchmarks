#!/bin/bash
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=1
#SBATCH --partition=dev
#SBATCH --ntasks-per-core=1
#SBATCH --time=00:05:00
#SBATCH --output=CommTest.out

export OMP_NUM_THREADS=1

aprun -B ./CommTest.exe --nprocx 2 --nprocy 4 --ie 52 --je 35 --enable-GCL --ntracer-perHandler 3 --nGCLHandlers 30 --lperi_y --disable-validation-report --nbl_exchg 2
################################################################
