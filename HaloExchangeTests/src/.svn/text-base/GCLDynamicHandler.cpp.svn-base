#include "GCLDynamicHandler.h"

void GCLDynamicHandler::Init(int numHandlers, int numTracerPerHandler, int nbl_exchg, int istart, int jstart, int iend, int jend, int ie,int je,int ke, MPI_Comm icomm_cart, bool lperi_x, bool lperi_y)
{
  GCLHandlerBase::InitImpl(numHandlers, numTracerPerHandler, ie, je, ke, icomm_cart, lperi_x, lperi_y);

  for(int i=0; i < numHandlers_; ++i) {

    gcl_handler_type* gclHandler = new gcl_handler_type(*pPeriodicity_, icomm_cart);
    gclHandler->add_halo<2>( // dimension index (k running first)
                0, 0, // no boundary in k dimension
                // first position in K of inner domain plus K halo region
                0,
                // last position in K of inner domain plus K halo region
                ke-1, 
                ke // total size of a stride
                );

  
    gclHandler->add_halo<1>( // dimension index (k running first)
                nbl_exchg, nbl_exchg, // no boundary in k dimension
                // first position in K of inner domain plus K halo region
                jstart-1,
                // last position in K of inner domain plus K halo region
                jend-1, 
                je // total size of a stride
                );
    gclHandler->add_halo<0>( // dimension index (k running first)
                nbl_exchg, nbl_exchg, // no boundary in k dimension
                // first position in K of inner domain plus K halo region
                istart-1,
                // last position in K of inner domain plus K halo region
                iend-1, 
                ie // total size of a stride
                );

    gclHandler->setup(numTracerPerHandler_);
    gclHandlerList_.push_back(gclHandler);
  }


}

void GCLDynamicHandler::DoExchange(double* exchgData)
{
  double start_t = MPI_Wtime();
  for(int h=0; h < numHandlers_; ++h) {
    gcl_handler_type* gclHandler = gclHandlerList_[h];

    std::vector<double*> tracerList;
    for (int l=0; l<numTracerPerHandler_; l++) {
      int offset = ie_*je_*ke_*(l+h*numTracerPerHandler_);
      tracerList.push_back( &exchgData[offset] );
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


