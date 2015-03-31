#pragma once
#include "halo_exchange.h"
#include "Definitions.h"
#include "GCLHandlerBase.h"

class GCLGenericHandler : public GCLHandlerBase
{
public:
  typedef GCL::halo_exchange_generic<GCL::layout_map<0,1,2>, 3, arch, GCL::version_manual > gcl_handler_type;

  GCLGenericHandler(){}
  ~GCLGenericHandler(){}  
 
  void Init(int numHandlers, int numTracerPerHandler, int nbl_exchg, int istart, int jstart, int iend, int jend, int ie,int je,int ke, MPI_Comm icomm_cart, 
           bool lperi_x, bool lperi_y);

  void DoExchange(double* exchgData);
private:
  std::vector<gcl_handler_type*> gclHandlerList_;
  GCL::array<GCL::halo_descriptor,3> halo_dsc_;

};

