#!/bin/bash
#SBATCH --job-name="halo_exchg" 
#SBATCH --nodes=2
#SBATCH --ntasks=2
#SBATCH --time=00:05:00 
#SBATCH --account=h03
#SBATCH --gres=gpu:1
#SBATCH --output=run_CPU_1OMP_512x32_ntrac1_gclHndl50_version0_castor_2x1_1ts_t12_test2.out
#SBATCH --error=run_CPU_1OMP_512x32_ntrac1_gclHndl50_version0_castor_2x1_1ts_t12_test2.out
#SBATCH --partition=express
#SBATCH --constraint=fermi
#SBATCH --cpu_bind=cores
#SBATCH --threads-per-core=1
#SBATCH --ntasks-per-core=1
#SBATCH --mem=24000

#======START=============================== 
echo "On which nodes it executes" 
echo $SLURM_JOB_NODELIST 
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/users/cosuna/projects/source/GCL/INTERNODE/HaloExchangeTests
export MV2_ENABLE_AFFINITY=0
export OMP_NUM_THREADS=1

#1D
mpirun -n 8  ./CommTest.exe --nprocx 1 --nprocy 8 --ie 520 --je 350 --enable-GCL --ntracer-perHandler 3 --nGCLHandlers 30 --lperi_y --disable-validation-report --nbl_exchg 2
#2D
mpirun -n 8  ./CommTest.exe --nprocx 2 --nprocy 4 --ie 520 --je 350 --enable-GCL --ntracer-perHandler 3 --nGCLHandlers 30 --lperi_y --disable-validation-report --nbl_exchg 2

#======END=================================
