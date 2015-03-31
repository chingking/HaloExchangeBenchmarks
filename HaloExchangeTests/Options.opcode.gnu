# Compilers
F90    =  gfortran

LD       = g++
CC       =   g++
cc       =   gcc

MPI_INC = -I$(MPI_ROOT)/include

# optimization flags
FOPT   =   -cpp -O3 -DPOLLEN -DNETCDF -DHP2C -fno-fast-math  -ffree-form  -ffree-line-length-none  -fno-backslash  -fimplicit-none  -fopenmp -pthread  -ffpe-trap=invalid,zero,overflow  -fbounds-check  -fcheck-array-temporaries  -finit-real=nan    -c -cpp  -DPOLLEN -DNETCDF -DHP2C -fno-fast-math  -ffree-form  -ffree-line-length-none  -fno-backslash  -fimplicit-none  -fopenmp -pthread  -ffpe-trap=invalid,zero,overflow  -fbounds-check  -fcheck-array-temporaries  -finit-real=nan
COPT   =  -Wall 
LOPT   =  -fopenmp

# includes 
BOOST_INC = -I/apps/dom/boost-1.52.0/include
GCL_L2_INC   =  -I../../gcl/L2/include
GCL_INC  =  $(GCL_L2_INC) -I../../gcl/L3/include -I../../gcl

#libraries
MPI_LIB  =  -L$(MPI_ROOT)/lib  -lmpich -lmpichcxx -lmpl 
SLIB      =  -lgfortran -ldl -lm  -lbfd  -lz -liberty 

CUDA_LIB = -L$(CUDALIB) -lcudart
CUDA_INC = -I/apps/opcode/CUDA-5.0/include
