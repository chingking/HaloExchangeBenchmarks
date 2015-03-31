#include <stdbool.h>

extern "C" {

#define i_global __parallel_utilities_MOD_i_global
  int i_global(int *i_loc);

#define j_global __parallel_utilities_MOD_j_global
extern int j_global(int *j_loc);

#define init_par_utilities __parallel_utilities_MOD_init_par_utilities
extern void init_par_utilities(int *numpex, int *numpey, int ipositions[],
  int *nbounds, int *icart_id );

#define halo_exchange __parallel_utilities_MOD_halo_exchange
extern void halo_exchange(double *sendbuf, int *sendbuflen, int *imp_type, int *icomm, int *num_compute,
  int *idim, int *jdim, int *kdim, int *jstartpar, int *jendpar, int *nlines, 
  int neighbors[], bool *lperi_x, bool *lperi_y, int *ntag, 
  double* packing_time, double* sending_time, double* waiting_time, double* unpacking_time, 
  double *var01 );

#define create_communicator __parallel_utilities_MOD_create_communicator
extern void create_communicator( int* nprocx, int* nprocy, bool* lperi_x, bool* lperi_y, int* icomm_cart );

}
