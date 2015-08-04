# Compilers
F90    =  gfortran

LD       = g++
CC       =   g++
cc       =   gcc

MPI_ROOT = /users/hamidouc/GPU/mvapich2/install/

MPI_INC = -I$(MPI_ROOT)/include

FOPT   =   -cpp  -DPOLLEN -DNETCDF -DHP2C -fno-fast-math  -ffree-form  -ffree-line-length-none  -fno-backslash  -fimplicit-none  -fopenmp -pthread  -fbacktrace  -fdump-core  -ffpe-trap=invalid,zero,overflow  -fbounds-check  -fcheck-array-temporaries  -finit-real=nan    -c -cpp  -DPOLLEN -DNETCDF -DHP2C -fno-fast-math  -ffree-form  -ffree-line-length-none  -fno-backslash  -fimplicit-none  -fopenmp -pthread 
COPT   =  -Wall
LOPT   =  -fopenmp

# includes 
BOOST_INC = -I/lus/scratch/olifu/kesch/BUILD/boost_1.49.0/include
GCL_L2_INC   =  -I../../gcl/L2/include
GCL_INC  =  $(GCL_L2_INC) -I../../gcl/L3/include -I../../gcl

#libraries
MPI_LIB  =  -L$(MPI_ROOT)/lib -lmpl -lmpichf90 -lmpich
SLIB      =  -lgfortran -ldl -lm   -lz  

CUDA_ROOT = /global/opt/nvidia/cudatoolkit/6.5.14/
CUDA_LIB = -L$(CUDA_ROOT)/lib64/ -lcudart
CUDA_INC = -I$(CUDA_ROOT)/include
