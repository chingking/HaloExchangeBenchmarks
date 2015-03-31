# Compilers
F90    =  ftn
CC       =   CC
cc       =   cc

# optimization flags
FOPT   =  -ffree-form -ffree-line-length-none -fno-backslash -fimplicit-none -O0 -g -fbounds-check
COPT   =  -Wall -O0 -g
LOPT   =  -O0 -g

# includes 
BOOST_INC = -I/users/mbianco/boost_1_47_0
GCL_L2_INC   =  -I../gcl/L2/include
GCL_INC  =  $(GCL_L2_INC) -I../gcl/L3/include -I../gcl

#libraries
GCL_LIB  =   -L../gcl/lib -lgcl
MPI_LIB  =  
SLIB      =  -lgfortran -ldl -lm 
