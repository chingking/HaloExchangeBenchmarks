#module load fftw
#module switch xt-libsci/11.1.00.11 xt-libsci/11.0.03
#module switch gcc/4.7.0 gcc/4.6.1

# Compilers
F90    =  ftn
CC       =   CC
cc       =   cc
LD = CC

# optimization flags
FOPT   =  -ffree-form -ffree-line-length-none -fno-backslash -fimplicit-none -O3 -fbounds-check
COPT   =  -Wall
LOPT   =  

# includes 
BOOST_INC = -I/apps/todi/boost/1.51.0/gnu_471/include/
GCL_L2_INC   =  -I../../gcl/L2/include
GCL_INC  =  $(GCL_L2_INC) -I../../gcl/L3/include -I../../gcl

#libraries
GCL_LIB  =   -L../../gcl/build/lib -lgcl
#MPI_LIB  =  -L/apps/todi/mpiP/3.3/gnu_462/lib -lmpiP
#SLIB      =   -ldl -lm -L/apps/todi/binutils/2.22/gnu_462/lib  -lbfd  -lz -liberty
SLIB      =   -ldl -lm -lz 

CUDA_INC = -I$(CUDATOOLKIT_HOME)/include/
CUDA_LIB = -L$(CUDATOOLKIT_HOME)/lib/ -lcudart
