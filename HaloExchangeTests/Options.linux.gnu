# Compilers
F90    =  gfortran
CC       =   g++
cc       =   gcc

MPI_INC = -I$(MPI_ROOT)/include

# optimization flags
FOPT   =   -cpp  -DPOLLEN -DNETCDF -DHP2C -fno-fast-math  -ffree-form  -ffree-line-length-none  -fno-backslash  -fimplicit-none  -fopenmp -pthread -O0  -g  -fbacktrace  -fdump-core  -ffpe-trap=invalid,zero,overflow  -fbounds-check  -fcheck-array-temporaries  -finit-real=nan    -c -cpp  -DPOLLEN -DNETCDF -DHP2C -fno-fast-math  -ffree-form  -ffree-line-length-none  -fno-backslash  -fimplicit-none  -fopenmp -pthread -O0  -g  -fbacktrace  -fdump-core  -ffpe-trap=invalid,zero,overflow  -fbounds-check  -fcheck-array-temporaries  -finit-real=nan
COPT   =  -Wall -O0 -g
LOPT   =  -O0 -g $(OMPFLAGS)

MPI_ROOT = /usr/lib64/mpi/gcc/openmpi/

# includes 
BOOST_INC = -I/users/mbianco/boost_1_47_0
MPI_INC  = -I$(MPI_ROOT)/include/
GCL_L2_INC   =  -I../../gcl/L2/include -I/usr/lib64/mpi/gcc/openmpi/include/
GCL_INC  =  $(GCL_L2_INC) -I../../gcl/L3/include -I../../gcl

#libraries
MPI_LIB  =  -L$(MPI_ROOT)/lib64/ -lmpi_f90 -lmpi_f77 -lmpi -lopen-rte -lopen-pal -ldl -Wl,--export-dynamic -lnsl -lutil -lm -ldl -lmpi_cxx
SLIB      =  -lgfortran -ldl -lm -L/apps/todi/binutils/2.22/gnu_462/lib  -lbfd  -lz -liberty

CUDA_LIB = -L/usr/local/cuda/lib64/ -lcudart
CUDA_INC = -I/usr/local/cuda/include/
