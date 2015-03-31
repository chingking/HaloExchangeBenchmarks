#ifndef _GCL_C_WRAPPER_H_
#define _GCL_C_WRAPPER_H_

#ifdef __cplusplus
extern "C" {
#endif
  void gcl_setup_cart( MPI_Comm cCartComm );
  void gcl_setup();
  void gcl_barrier();
  void gcl_pid(int *);

  void gcl_obtain_handler_d(void**);
  void gcl_register_field_d(void**, double*, int*);
  void gcl_register_halo_d(void**, int*, int*, int*, int*, int*, int*, int*);
  void gcl_prepare_pattern(void**);
  void gcl_perform_exchange(void**);

  void gcl_dynamic_handler_d(int*,bool*,bool*);
  void gcl_dynamic_halo_d(int*, int*, int*, int*, int*, int*, int*);
  void gcl_prepare_dynamic_d(int*, int*);
  void gcl_dynamic_exchange_d_1(int*, double*);

#ifdef __cplusplus
}
#endif

#endif
