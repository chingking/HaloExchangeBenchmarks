#include "GCLGenericHandler.h"

void GCLGenericHandler::Init(int numHandlers, int numTracerPerHandler, int nbl_exchg, int istart, int jstart, int iend, int jend, int ie,int je,int ke, MPI_Comm icomm_cart, bool lperi_x, bool lperi_y)
{
  GCLHandlerBase::InitImpl(numHandlers, numTracerPerHandler, ie, je, ke, icomm_cart, lperi_x, lperi_y);

  halo_dsc_[2] = GCL::halo_descriptor(0, 0, 0, ke-1, ke);
  halo_dsc_[1] = GCL::halo_descriptor(nbl_exchg, nbl_exchg, jstart -1, jend-1, je);
  halo_dsc_[0] = GCL::halo_descriptor(nbl_exchg, nbl_exchg, istart-1, iend-1, ie);

  for(int i=0; i < numHandlers_; ++i) {
    gcl_handler_type* gclHandler = new gcl_handler_type(*pPeriodicity_, icomm_cart);
    gclHandler->setup(numTracerPerHandler, GCL::field_on_the_fly<int,GCL::layout_map<0,1,2>, gcl_handler_type::traits>(NULL,halo_dsc_), sizeof(double) );
    gclHandlerList_.push_back(gclHandler);
  }
}

void GCLGenericHandler::DoExchange(double* exchgData)
{

  double start_t = MPI_Wtime();

  for(int h=0; h < numHandlers_; ++h) {
    gcl_handler_type* gclHandler = gclHandlerList_[h];

    std::vector<GCL::field_on_the_fly<double, GCL::layout_map<0, 1, 2>, gcl_handler_type::traits> > tracerList;
    for (int l=0; l<numTracerPerHandler_; l++) {
      int offset = ie_*je_*ke_*(l+h*numTracerPerHandler_);
      GCL::field_on_the_fly<double, GCL::layout_map<0, 1, 2>, gcl_handler_type::traits> field_otf((double*) &exchgData[offset], halo_dsc_);
      tracerList.push_back( field_otf );
    }
    double start, end;
    MPI_Barrier(icomm_cart_);
    start = MPI_Wtime();
    gclHandler->pack( tracerList );
#ifdef _GCL_GPU_
    cudaDeviceSynchronize();
#endif
    end = MPI_Wtime();
    etime_gcl_[comm_pack] += end - start;

    MPI_Barrier(icomm_cart_);

    start = MPI_Wtime();
    gclHandler->start_exchange();
    end = MPI_Wtime();
    etime_gcl_[comm_send] +=  end - start;

    start  = MPI_Wtime();
    gclHandler->wait();
#ifdef _GCL_GPU_
    cudaDeviceSynchronize();
#endif
    end = MPI_Wtime();
    etime_gcl_[comm_wait] += end - start;
      
    MPI_Barrier(icomm_cart_);
    start = MPI_Wtime();
    gclHandler->unpack( tracerList );
#ifdef _GCL_GPU_
    cudaDeviceSynchronize();
#endif
    end = MPI_Wtime();
    etime_gcl_[comm_unpack] += end - start;

  }

  etime_gcl_[comm_total] += MPI_Wtime() - start_t;

}


