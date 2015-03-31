//#include <GCL.h>
#include <iostream>
#include <vector>
#include "mpi.h"
#include "GCL.h"

int pid, nprocs;
MPI_Comm CartComm;


// The following setup assumes a C cartesian communicator is provided by client
void gcl_setup_cart(MPI_Comm cCartComm) {
  GCL::GCL_Init();

  int pid;
  MPI_Comm_rank(GCL::GCL_WORLD, &pid);
  int nprocs;
  MPI_Comm_size(GCL::GCL_WORLD, &nprocs);

  CartComm = cCartComm;

  int dims[3]={0,0,0}, periods[3]={true,true,true}, coords[3]={0,0,0};
    MPI_Cart_get( CartComm, 3, dims, periods/*does not really care*/, coords);

}

void gcl_setup() {
  GCL::GCL_Init();

  int pid;
  MPI_Comm_rank(GCL::GCL_WORLD, &pid);
  int nprocs;
  MPI_Comm_size(GCL::GCL_WORLD, &nprocs);


  int dims[3] = {0,0,1};
  MPI_Dims_create(nprocs, 2, &(dims[0])); 
  int t=dims[0];
  dims[0]=dims[1];
  dims[1]=t;
  int period[3] = {1, 1, 1};
  MPI_Cart_create(GCL::GCL_WORLD, 3, dims, period, true, &CartComm);
}

void gcl_barrier() {
  MPI_Barrier(GCL::GCL_WORLD);
}

void gcl_pid(int * pid) {
  MPI_Comm_rank(GCL::GCL_WORLD, pid);
  //  *pid = GCL::PID;
}

