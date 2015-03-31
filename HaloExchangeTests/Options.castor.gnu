# Compilers
F90    =  gfortran

LD       = g++
CC       =   g++
cc       =   gcc

MPI_INC = -I$(MPI_ROOT)/include

# optimization flags
FOPT   =   -cpp  -DPOLLEN -DNETCDF -DHP2C -fno-fast-math  -ffree-form  -ffree-line-length-none  -fno-backslash  -fimplicit-none  -fopenmp -pthread -O0  -g  -fbacktrace  -fdump-core  -ffpe-trap=invalid,zero,overflow  -fbounds-check  -fcheck-array-temporaries  -finit-real=nan    -c -cpp  -DPOLLEN -DNETCDF -DHP2C -fno-fast-math  -ffree-form  -ffree-line-length-none  -fno-backslash  -fimplicit-none  -fopenmp -pthread -O0  -g  -fbacktrace  -fdump-core  -ffpe-trap=invalid,zero,overflow  -fbounds-check  -fcheck-array-temporaries  -finit-real=nan
COPT   =  -Wall 
LOPT   =  -fopenmp

# includes 
BOOST_INC = -I/users/mbianco/boost_1_47_0
GCL_L2_INC   =  -I../../gcl/L2/include
GCL_INC  =  $(GCL_L2_INC) -I../../gcl/L3/include -I../../gcl

#libraries
MPI_LIB  =  -L$(MPI_ROOT)/lib -lmpl -lmpichf90
SLIB      =  -lgfortran -ldl -lm  -lbfd  -lz -liberty 

CUDA_LIB = -L/apps/castor/CUDA-4.2/cuda/lib64 -lcudart
CUDA_INC = -I/apps/castor/CUDA-4.2/cuda/include/
