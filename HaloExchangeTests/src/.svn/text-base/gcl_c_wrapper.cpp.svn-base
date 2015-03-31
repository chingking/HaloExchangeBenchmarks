#include <GCL.h>
#include <descriptors.h>
#include <iostream>
#include <vector>
#include <gcl_c_wrapper.h>
//#include "halo_exchange.h"

int pid, nprocs;
MPI_Comm CartComm;

typedef GCL::gcl_utils::boollist<3> CYC;

CYC *cyc;

GCL::handler_manager_ut<
  double, 3, 
  GCL::Halo_Exchange_3D<GCL::MPI_3D_process_grid_t<CYC > > > hndl_mngr;

typedef GCL::hndlr_descriptor_ut<
  double,3,  
  GCL::Halo_Exchange_3D<GCL::MPI_3D_process_grid_t<CYC > > > handler_d;

typedef  GCL::hndlr_dynamic_ut<
  double,3, 
  GCL::Halo_Exchange_3D<GCL::MPI_3D_process_grid_t<CYC > > > dynamic_type;
//typedef GCL::halo_exchange_dynamic_ut<typename GCL::layout_map<0, 1, 2>, GCL::layout_map<0, 1, 2>, double, 3, GCL::gcl_cpu> dynamic_type;

std::vector<dynamic_type*> dynamic_list;

//void gcl_setup(MPI_Comm) {....}

// The following setup assumes a C cartesian communicator is provided by client
void gcl_setup_cart(MPI_Comm cCartComm) {
  GCL::GCL_Init();

  int pid;
  MPI_Comm_rank(GCL::GCL_WORLD, &pid);
  int nprocs;
  MPI_Comm_size(GCL::GCL_WORLD, &nprocs);

  std::cout << pid << " " << nprocs << "\n";

  CartComm = cCartComm;

  int dims[3]={0,0,0}, periods[3]={true,true,true}, coords[3]={0,0,0};
    MPI_Cart_get( CartComm, 3, dims, periods/*does not really care*/, coords);
printf("here coords for %d: %d,%d,%d\n",pid,coords[0],coords[1],coords[2]);

}

void gcl_setup() {
  GCL::GCL_Init();

  int pid;
  MPI_Comm_rank(GCL::GCL_WORLD, &pid);
  int nprocs;
  MPI_Comm_size(GCL::GCL_WORLD, &nprocs);

  std::cout << pid << " " << nprocs << "\n";

  int dims[3] = {0,0,1};
  MPI_Dims_create(nprocs, 2, &(dims[0])); 
  int t=dims[0];
  dims[0]=dims[1];
  dims[1]=t;
  int period[3] = {1, 1, 1};
  std::cout << "@" << pid << "@ MPI GRID SIZE " << dims[0] << " - " << dims[1] << " - " << dims[2] << "\n";
  std::cout << "@" << pid << "   ff " << period[0] << "  " << period[1] << "  " << period[2] << std::endl;
  MPI_Cart_create(GCL::GCL_WORLD, 3, dims, period, true, &CartComm);
}

void gcl_barrier() {
  MPI_Barrier(GCL::GCL_WORLD);
}

void gcl_pid(int * pid) {
  MPI_Comm_rank(GCL::GCL_WORLD, pid);
  std::cout << *pid << "----------------------------------------------------------------\n";
  //  *pid = GCL::PID;
}

void gcl_obtain_handler_d(void** handler) {
  cyc = new CYC(true,true,true);
  handler_d *hd = &(hndl_mngr.create_handler(*cyc,CartComm));
  *handler = reinterpret_cast<void*>(hd);
}

void gcl_register_field_d(void** handler,  double* field_pointer, int* field_index) {
  *field_index = reinterpret_cast<handler_d*>(*handler)->register_field(field_pointer);
}


void gcl_register_halo_d(void** handler, 
                     int* field_index, 
                     int* dim_index, 
                     int* minus, 
                     int* plus, 
                     int* begin, 
                     int* end, 
                     int* len) {

  int pid;
    MPI_Comm_rank(GCL::GCL_WORLD, &pid);

#ifndef NDEBUG
  std::cout << pid << " " 
            << "tuple : "
            << dim_index << " "
            << *minus << " " 
            << *plus << " " 
            << *begin << " " 
            << *end << " " 
            << *len << "\n"; 
#endif
  reinterpret_cast<handler_d*>(*handler)->register_halo(*field_index, 
                                                      *dim_index,
                                                      *minus,
                                                      *plus,
                                                      *begin,
                                                      *end,
                                                      *len);
}
 
void gcl_prepare_pattern(void **handler) {
  reinterpret_cast<handler_d*>(*handler)->allocate_buffers();
}

void gcl_perform_exchange(void **handler) {
  reinterpret_cast<handler_d*>(*handler)->pack();
  reinterpret_cast<handler_d*>(*handler)->exchange();
  reinterpret_cast<handler_d*>(*handler)->unpack();
}

void gcl_dynamic_handler_d(int* _hd,bool* lperi_x,bool* lperi_y) {
  cyc = new CYC(*lperi_x,*lperi_y,true);
  dynamic_type *hd = new dynamic_type(*cyc,CartComm);
  dynamic_list.push_back(hd);
  *_hd = dynamic_list.size()-1;
}

void gcl_dynamic_halo_d(int* hd, 
                        int* dim_index, 
                        int* minus, 
                        int* plus, 
                        int* begin, 
                        int* end, 
                        int* len) {
  dynamic_list[*hd]->halo.add_halo(*dim_index,
                                   *minus,
                                   *plus,
                                   (*begin)-1,
                                   (*end)-1,
                                   *len);
  int pid;
    MPI_Comm_rank(GCL::GCL_WORLD, &pid);

#ifndef NDEBUG
  std::cout << pid << " " 
            << "tuple : "
            << *dim_index << " "
            << *minus << " " 
            << *plus << " " 
            << (*begin)-1 << " " 
            << (*end)-1 << " " 
            << *len << "\n"; 
#endif
}

void gcl_prepare_dynamic_d(int* hd, int* howmany) {
  dynamic_list[*hd]->allocate_buffers(*howmany);
}

void gcl_dynamic_exchange_d_1(int* handler, double* v) {
  dynamic_list[*handler]->pack(v);
  dynamic_list[*handler]->exchange();
  dynamic_list[*handler]->unpack(v);
} 



// Functions to query proc grids for position information
// Am I at boundary?
// How much overlap?
