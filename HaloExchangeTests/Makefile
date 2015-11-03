# directory definitions
SRCDIR     := src
OBJDIR     := obj
ifeq (0,${MAKELEVEL})
  ROOT      := $(shell pwd)
  VPATH     := .:$(ROOT)/$(SRCDIR)
endif

-include $(ROOT)/Options
ifndef F90
  $(error Must $(ROOT) ${MAKELEVEL} create platform specific Options file; Some default Options files are provided, link any of them: ln -s Options.flavour Options)
endif


GCL_FLAGS = -D_GCL_MPI_
GCL_LIB  =   -L../../gcl/build/lib -lgcl

LIB    := $(MPI_LIB) $(SLIB) $(GCL_LIB)
CPPFLAGS   += $(GCL_INC) $(BOOST_INC) $(MPI_INC) $(GCL_FLAGS) -O3 -DNDEBUG -I../src
#DEBUG = 1

#CPPFLAGS += -DNDEBUG 

ifdef GPU
#  CPPFLAGS = -D_GCL_MPI_ -D_GCL_GPU_ -DDEVS_PER_NODE=1
#  COPT = ""
  CPPFLAGS += -D_GCL_GPU_ -DDEVS_PER_NODE=1 
  CPPFLAGS += $(CUDA_INC) -I$(MPICH_DIR)/include
  NVCPPFLAGS = $(GCL_FLAGS) $(CPPFLAGS) -D_GCL_GPU_ -DDEVS_PER_NODE=1 -arch=sm_35 --x=cu 
  LIB      += $(CUDA_LIB)
else
## comment in the next line for gcl manual packing (instead of mpi datatypes)
#  CPPFLAGS += -DMANUAL_PACKING
   CPPFLAGS += -fopenmp
   NVCC = $(CC)
   NVCPPFLAGS = $(CPPFLAGS) 
   CPPFLAGS += $(COPT) -O3 -DNDEBUG
endif


export ROOT VPATH GPU 

.PHONY : gpu cpu

.IGNORE: clean

default: cpu

gpu:
	@$(MAKE) -C $(OBJDIR) -f $(ROOT)/Makefile GPU=1 all

cpu:
	@$(MAKE) -C $(OBJDIR) -f $(ROOT)/Makefile all 


all: $(ROOT)/CommTest.exe

$(ROOT)/CommTest.exe: fortran_comm.o CommTest.o gcl_functions.o GnuplotLogger.o GCLGenericHandler.o GCLDynamicHandler.o CommStatsCollector.o
	$(LD) ${LOPT} -O3 $^ -o $@ $(MPIF_LIB) $(MPICXX_LIB) $(LIB)


.SUFFIXES: .o .f90 .cpp

%.o: %.f90
	$(F90) -c $(FOPT) $(MPI_INC) $< -o $@
%.o: %.c
	$(cc) -c $(NVCPPFLAGS) $< -o $@
%.o: %.cpp
	$(NVCC) $(NVCPPFLAGS) -c $< -o $@

CommTest.o: CommTest.cpp
	$(CC) $(CPPFLAGS) -c $< -o $@

clean:
	rm -f $(OBJDIR)/*.o CommTest.exe 

