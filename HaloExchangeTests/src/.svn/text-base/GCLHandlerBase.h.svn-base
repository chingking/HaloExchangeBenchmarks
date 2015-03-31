#pragma once

class GCLHandlerBase
{
public:
  GCLHandlerBase()
  {
    ResetTimers();
  }
  ~GCLHandlerBase()
  {
    delete pPeriodicity_;
  }

  void ResetTimers()
  {
    for(int i=0; i < comm_nitems; ++i) {
      etime_gcl_[i] = 0;
    }
  }


  virtual void Init(int numHandlers, int numTracerPerHandler, int nbl_exchg, int istart, int jstart, int iend, int jend, int ie,int je,int ke, 
             MPI_Comm icomm_cart, bool lperi_x, bool lperi_y) = 0;

  virtual void DoExchange(double* exchgData) =0;


  void InitImpl(int numHandlers, int numTracerPerHandler,int ie,int je,int ke, MPI_Comm icomm_cart, bool lperi_x, bool lperi_y)
  {
    numHandlers_ = numHandlers;
    numTracerPerHandler_ = numTracerPerHandler;
    ie_ = ie;
    je_ = je;
    ke_ = ke;
    icomm_cart_ = icomm_cart;
    pPeriodicity_ = new GCL::gcl_utils::boollist<3>(lperi_x, lperi_y, true);
  }
  double* pEtime_gcl() {return &etime_gcl_[0];}

protected:
  GCL::gcl_utils::boollist<3>* pPeriodicity_;
  int numHandlers_;
  int numTracerPerHandler_;
  int ie_, je_, ke_;
  double etime_gcl_[comm_nitems];
  MPI_Comm icomm_cart_;
 
};

