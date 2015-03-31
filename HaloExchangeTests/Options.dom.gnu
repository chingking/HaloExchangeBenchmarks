# Compilers
F90    =  gfortran

LD       = g++
CC       =   g++
cc       =   gcc

GCL_FLAGS = -D_GCL_MPI_

CPPFLAGS =  -fopenmp $(GCL_FLAGS)

MPI_INC = -I$(MPI_ROOT)/include

FOPT   =   -cpp  -DPOLLEN -DNETCDF -DHP2C -fno-fast-math  -ffree-form  -ffree-line-length-none  -fno-backslash  -fimplicit-none  -fopenmp -pthread  -fbacktrace  -fdump-core  -ffpe-trap=invalid,zero,overflow  -fbounds-check  -fcheck-array-temporaries  -finit-real=nan    -c -cpp  -DPOLLEN -DNETCDF -DHP2C -fno-fast-math  -ffree-form  -ffree-line-length-none  -fno-backslash  -fimplicit-none  -fopenmp -pthread 
COPT   =  -Wall -03 
LOPT   =  -O3 -fopenmp

# includes 
BOOST_INC = -I/apps/dom/boost-1.52.0/include
GCL_L2_INC   =  -I../../gcl/L2/include
GCL_INC  =  $(GCL_L2_INC) -I../../gcl/L3/include -I../../gcl

#libraries
GCL_LIB  =   -L../../gcl_build/lib -lgcl
MPI_LIB  =  -L$(MPI_ROOT)/lib -lmpl -lmpichf90 -lmpich
SLIB      =  -lgfortran -ldl -lm   -lz -liberty 

CUDA_LIB = -L$(CUDALIB) -lcudart
CUDA_INC = -I/apps/castor/CUDA-5.0/include
