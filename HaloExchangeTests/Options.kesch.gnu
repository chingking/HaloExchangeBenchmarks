# Compilers
F90    =  gfortran

LD       = g++
CC       =   g++
cc       =   gcc


## NON GDR
MPI_ROOT = /apps/escha/easybuild/software/MVAPICH2/2.0.1-GCC-4.8.2-EB/
MPI_LIBDIR = $(MPI_ROOT)/lib
## GDR
#MPI_ROOT = /opt/mvapich2/gdr/2.1/cuda6.5/gnu/
#MPI_LIBDIR = $(MPI_ROOT)/lib64

MPI_INC = -I$(MPI_ROOT)/include

FOPT   =   -cpp  -DPOLLEN -DNETCDF -DHP2C -fno-fast-math  -ffree-form  -ffree-line-length-none  -fno-backslash  -fimplicit-none  -fopenmp -pthread  -fbacktrace  -fdump-core  -ffpe-trap=invalid,zero,overflow  -fbounds-check  -fcheck-array-temporaries  -finit-real=nan    -c -cpp  -DPOLLEN -DNETCDF -DHP2C -fno-fast-math  -ffree-form  -ffree-line-length-none  -fno-backslash  -fimplicit-none  -fopenmp -pthread 
COPT   =  -Wall
LOPT   =  -fopenmp

# includes 
BOOST_INC = -I/scratch/olifu/kesch/BUILD/boost_1.49.0/include
GCL_L2_INC   =  -I../../gcl/L2/include
GCL_INC  =  $(GCL_L2_INC) -I../../gcl/L3/include -I../../gcl

#libraries
MPI_LIB  =  -L$(MPI_LIBDIR) -lmpl -lmpichf90 -lmpich
SLIB      =  -lgfortran -ldl -lm   -lz  

CUDA_ROOT = /global/opt/nvidia/cudatoolkit/6.5.14/
CUDA_LIB = -L$(CUDA_ROOT)/lib64/ -lcudart
CUDA_INC = -I$(CUDA_ROOT)/include
