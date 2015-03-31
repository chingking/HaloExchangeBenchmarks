#!/bin/bash
#SBATCH --job-name="halo_exchg" 
#SBATCH --nodes=4
#SBATCH --ntasks=4
#SBATCH --time=00:05:00 
#SBATCH --account=s83
#SBATCH --partition=dev
#SBATCH --gres=gpu:1
#SBATCH --output=run_ntrac20_gclHndl10_version0_8OMP_todi_4x4_opt_b.out
#SBATCH --error=run_ntrac20_gclHndl10_version0_8OMP_todi_4x4_opt_b.out
#SBATCH --cpus-per-task=1

#======START=============================== 
echo "On which nodes it executes" 
echo $SLURM_JOB_NODELIST
#module load cudatoolkit
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/scratch/todi/cosuna/GCL/HaloExchangeTests 
export OMP_NUM_THREADS=8
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/apps/todi/binutils/2.22/gnu_462/lib
export MPICH_RDMA_ENABLED_CUDA=1
#1D
aprun -n 8 -N 1 -d 1 ./CommTest.exe --nprocx 1 --nprocy 8 --ie 520 --je 350 --enable-GCL --ntracer-perHandler 3 --nGCLHandlers 30 --lperi_y --disable-validation-report --nbl_exchg 2
#2D
aprun -n 8 -N 1 -d 1 ./CommTest.exe --nprocx 2 --nprocy 4 --ie 520 --je 350 --enable-GCL --ntracer-perHandler 3 --nGCLHandlers 30 --lperi_y --disable-validation-report --nbl_exchg 2


#======END=================================
