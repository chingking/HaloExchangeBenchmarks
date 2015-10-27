/*------------------------------------------------------------------------------------

 CommTest.cpp  -  halo exchange standalone

 Author: Carlos Osuna
 Email : carlos.osuna@env.ethz.ch

 This program performs halo exchanges of 3D fields using GCL & COSMO exchange boundaries subroutine
 It can be used in order to validate GCL exchanges or to extract performance numbers for the communication.

 The program can be configured in to run fortran exchanges, GCL exchanges or both.

 GCL configuration:
   * packing_version (in Definitions.h) defines a mode for packing in GCL:
        GCL::version_manual -> it uses custom kernels that pack all halos in a linear buffer used for communication.
        GCL::version_datatype -> it uses mpi datatypes. 

   ntracer_perHandler sets the number of fields being exchanged per communication (gcl handler).
   nGCLHandlers sets the number of communications (each one exchanging ntracer_perHandler fields).

   random flag will generate random data and compare the fortan with the GCL exchange. If non random data is chosen,
     then fields are filled with global indices (associated to the position of the point in the grid) and exchanges are
     validated with expected results.

   For cpu compilation, 'make cpu'
   For gpu compilation, 'make gpu'   (only packing_version 0 is permitted)

-------------------------------------------------------------------------------------*/
#ifdef _GCL_GPU_
#include <cuda.h>
#include <cuda_runtime.h>
#endif

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
// Not sure if this will be available under Windows
// needed for getopt
#include <unistd.h>
#include <getopt.h>
#include <sstream>

#include "mpi.h"
#include "fortran_comm.h"
#include "gcl_functions.h"
#include "halo_exchange.h"
#include "GnuplotLogger.h"
#include "Definitions.h"
#ifdef generic_pattern
#include "GCLGenericHandler.h"
#else
#include "GCLDynamicHandler.h"
#endif
#include "CommStatsCollector.h"

#define MIN(X,Y) ((X) < (Y) ? : (X) : (Y))
#define MAX(X,Y) ((X) > (Y) ? : (X) : (Y))



// MPI environment
int my_world_id;
MPI_Group igroup_world;
MPI_Comm icomm_cart, icomm_world;

// fortran cartesian communicator
int Ficomm_cart;

int my_cart_id, my_cart_pos[2], my_cart_neigh[4];
MPI_Datatype imp_reals, imp_integers;
int Fimp_reals, Fimp_integers;

// computation grid setup
int nproc, nprocx, nprocy, nprocio;

// global grid setup
// ie_tot_plusHalos x je_tot_plusHalos correspond to a complete domain composed by all subdomains including
//     all the halos for each subdomain
int ie_tot, je_tot, ie_tot_plusHalos, je_tot_plusHalos;
int ke_tot;
// we can define number of tracer for each gcl GCL handler and the number of GCL handler (i.e. independent exchanges)
int ntracer, ntracer_perHandler, nGCLHandlers;

// local grid setup
int ie, je, ke, nboundlines, *isubpos[4], *isubposa;
int ie_max, je_max;
// sign of data of a subdomain (depends on its rank in the 2D decomposition)
int sign;
int istart, jstart, iend, jend;
int istartpar, jstartpar, iendpar, jendpar;

// flags to enable FORTRAN COSMO testing or GCL testing
int fTest,gclTest;

int disableValidationReport;

// flag to generate random data (instead of global indices).
// While using random data, verification can only by carried out by comparing fortran & GCL communication
//     (and not agains global indices).
// Therefore to be used only if fTest && gclTest
int GenerateRandoms;

// halo-exchange setup
int nbl_exchg;
bool lperi_x, lperi_y;
int iperi_x, iperi_y;

// buffers to allocate data for both types of exchanges;
double *fData, *gclData;

// Grid of subdomain size that contains 0, 1 depending on wether the grid point was validated or not
double *fvalidatedGrid;


/* Grid of full domain size that contains codes that determine wether the grid point passed validation or not.
*  0  : point of an inner subdomain with correct result
*  1  : halo point of a subdomain with correct result
*  5  : point of an inner subdomain with incorrect result
*  6  : halo point of a subdomain with incorrect result
*  -1 : Undefined (value not set)
*/
double *validationCodes_tot;

// time stepping
int nstop, nstep;

// size of the full domain matrix
size_t matrix_size;


// timings
double etime;
double etime_fortran[comm_nitems];

int data_transfered_perHandler;

// object that plots the results in a gnuplot file in case validation fails
GnuplotLogger plotter;

#ifdef generic_pattern
GCLGenericHandler gclHandler;
#else
GCLDynamicHandler gclHandler;
#endif

//std::vector<gcl_handler_type*> gclHandlerList;


// initializes the environment
void init_environment(int argc, char **argv);

// initializes the grid of PEs
void init_procgrid();

// domain decomposition
void domain_decomposition();
void grid_constants();

// produces the gnuplot summary and file in case validation fails
void ReportGnuplot();

// function that allocate a buffer with global indices, random data,
// or with a copy of a already initialized buffer if the pointer argument is not NULL
double* InitBuffers(double* );

// Do actual exchange for FORTRAN COSMO
void FortranExchange(double *data);
// Do actual exchange for GCL
void GCLExchange(double *data);

/*
* Validate the results after the exchanges
* if second argument is not null, data of the two matrices is compared for exact results,
* Otherwise, the result is compared with expected global indices after a halo exchange
*/
void Validate(double* result, double* expected=NULL);

// Gather the array of validated grid points of every subdomain
void GatherValidatedGrid();

// Gather timing values
void GatherTimings();

// update the timing values of fortran exchanges
void update_fortran_times(double fpack_time, double fsend_time, double fwait_time, double funpack_time);

/*
* combine arrays of a gathered field (from MPI_Gather) into a full domain matrix
* @param includeHalos determines wether the full matrix should contain halos of the subdomains or not
*/
void combine_subarrays( double* gathered_field, double* dest, bool includeHalos );

int main(int argc, char **argv)
{
  int id;
  int c;
  int ierror =-1;

  // ************************************
  // initializations

  // setup grid (input_lmgrid in src_setup.f90)
  //   ie_tot, je_tot, ke_tot, ntracer

  ie_tot  = 128;
  je_tot  = 112;
  nboundlines = 3;
  ke_tot  = 60;
  ntracer_perHandler = 50;
  nGCLHandlers =  1;
  // setup processor and halo-exchange configuration (input_lmgrid in src_setup.f90)
  //   nstop, nprocx, nprocy, nboundlines, nbl_exchg, lreorder,
  //   lperi_x, lperi_y

  nstop  = 1;
  nprocx = 2;
  nprocy = 2;
  nprocio = 0;
  nbl_exchg = 3;
  lperi_x = false;
  lperi_y = false;
  fTest = 1;
  gclTest = 0;
  disableValidationReport = 0;
  GenerateRandoms = 0;

  iperi_x = 0;
  iperi_y = 0;
  while (1)
  {
    static struct option long_options[] =
    {
       /* These options set a flag. */
       {"lperi_x", no_argument,       &iperi_x, 1},
       {"lperi_y", no_argument,       &iperi_y, 1},
               /* These options don't set a flag.
 *  *                   We distinguish them by their indices. */
       {"nprocx",  required_argument,       0, 'a'},
       {"nprocy",  required_argument,       0, 'b'},
       {"ie"    ,  required_argument,       0, 'c'},
       {"je"    ,  required_argument,       0, 'd'},
       {"ke"    ,  required_argument,       0, 'e'},
       {"nblines", required_argument,       0, 'f'},
       {"ntracer-perHandler", required_argument, 0, 'g'},
       {"nbl_exchg",required_argument,      0, 'h'},
       {"nGCLHandlers", required_argument, 0, 'i'},
       {"nstop",    required_argument,     0, 'j'},
       {"enable-Fortran", no_argument,      &fTest,    1 },
       {"enable-GCL",     no_argument,      &gclTest,  1 },
       {"disable-Fortran",no_argument,      &fTest,    0 },
       {"disable-GCL",    no_argument,      &gclTest,  0 },
       {"disable-validation-report", no_argument,  &disableValidationReport, 1},
       {"random",   no_argument,       &GenerateRandoms,     1 },
       {0,0,0,0}
    };
           /* getopt_long stores the option index here. */
    int option_index = 0;

    c = getopt_long (argc, argv, "a:b:c:d:e:",
                            long_options, &option_index);

    lperi_x = (bool)iperi_x;
    lperi_y = (bool)iperi_y;
    // Detect the end of the options.
    if (c == -1)
      break;

    switch (c)
    {
      case 0:
      // If this option set a flag, do nothing else now.
        if (long_options[option_index].flag != 0)
          break;
          printf ("option %s", long_options[option_index].name);
        if (optarg)
          printf (" with arg %s", optarg);
          printf ("\n");
          break;

      case 'a':
        nprocx = atoi(optarg);
        break;

      case 'b':
        nprocy = atoi(optarg);
        break;

      case 'c':
        ie_tot =  atoi(optarg);
        break;

      case 'd':
        je_tot = atoi(optarg);
        break;

      case 'e':
        ke_tot = atoi(optarg);
        break;

      case 'f':
        nboundlines = atoi(optarg);
        break;

      case 'g':
        ntracer_perHandler = atoi(optarg);
        break;

      case 'h':
        nbl_exchg = atoi(optarg);
        break;

      case 'i':
        nGCLHandlers = atoi(optarg);
        break;

      case 'j':
        nstop = atoi(optarg);
        break;

      default:
        MPI_Abort (MPI_COMM_WORLD,ierror);
    }
  }

  // If there are remaining command line arguments (not options)
  if (optind < argc)
  {
    printf ("Error specifying non-option arguments: ");
    while (optind < argc)
      printf ("%s ", argv[optind++]);
      putchar ('\n');
      MPI_Abort (MPI_COMM_WORLD,ierror);
  }


// initialize the tests
  init_environment(argc, argv);
                                  	
  if(my_world_id == 0) {
    if (nprocx*nprocy != nproc) {
      printf("ERROR *** number of MPI tasks should be set to %i!  (nprocx= %d, nprocy= %d)\n",nproc,nprocx,nprocy);
      return(1);
    }
    if( nbl_exchg > nboundlines) {
      printf("ERROR *** number of lines to exchange can not be larger than number of boundary lines\n");
      return(1);
    }
    if( !fTest && !gclTest ) {
      printf("ERROR *** Either Fortran Test or GCL Test has to be enable\n");
      MPI_Abort(MPI_COMM_WORLD,ierror);
    }
    if( GenerateRandoms &&  ( !fTest || !gclTest ) ) {
      printf("ERROR *** If random flag is switched on to generate random data, only comparison between FORTRAN COSMO and GCL makes sense. Please enable both, fortran and gcl\n");
      MPI_Abort(MPI_COMM_WORLD,ierror);
    }

  }


  // tell user what we are doing
  if (my_world_id == 0) {
    printf("HALO-EXCHANGE STANDALONE RUN\n");
    printf("  ntracer-perHandler  = %i\n",ntracer_perHandler);
    printf("  nGCLHandlers        = %i\n",nGCLHandlers);
    printf("  ie_tot              = %i\n",ie_tot);
    printf("  je_tot              = %i\n",je_tot);
    printf("  ke_tot              = %i\n",ke_tot);
    printf("  nprocx              = %i\n",nprocx);
    printf("  nprocy              = %i\n",nprocy);
    printf("  ntracer             = %i\n",ntracer);
    printf("  lperi_x             = %i\n", lperi_x);
    printf("  lperi_y             = %i\n", lperi_y);
    printf("  # bound lines       = %i\n", nboundlines);
    printf("  nbl_exchg        = %i\n", nbl_exchg);
    printf("  enable--Fortran     = %i\n", fTest);
    printf("  enable--GCL         = %i\n", gclTest);
#ifdef _GCL_GPU_
    printf("  compiled for GPU mode \n" );
#else
    printf("  compiled for CPU mode \n");
#endif
    unsigned long new_mask = 2;
    unsigned int len = sizeof(new_mask);
    cpu_set_t cur_mask;
    pid_t p = 0;
    int ret;

    ret = sched_getaffinity(p, len, &cur_mask);
    printf("  sched_getaffinity = %d, cur_mask = %08lx\n", ret, cur_mask);
  }

  for(id=0; id<nprocx*nprocy; id++) {
    MPI_Barrier(icomm_cart);
#ifdef _GCL_GPU_
    int my_dev_id;
    cudaGetDevice(&my_dev_id);
#endif
    if( id == my_cart_id) {
      printf("***************\n");
      printf("       task = %5i\n",id);
      printf("         ie = %5i  je = %5i  ke = %5i\n",ie,je,ke);
      printf("   cart-pos = %5i, %5i\n",my_cart_pos[0],my_cart_pos[1]);
      printf(" cart-neigh = %5i, %5i, %5i, %5i\n",my_cart_neigh[0],my_cart_neigh[1],my_cart_neigh[2],my_cart_neigh[3]);
      printf("      i-pos = %5i - %5i\n",isubpos[0][id],isubpos[2][id]);
      printf("      j-pos = %5i - %5i\n",isubpos[1][id],isubpos[3][id]);
      printf("      sign of data = %i\n", sign);
#ifdef _GCL_GPU_
      printf("      CUDA device = %5i\n", my_dev_id);
#endif
    }
  }

  ie_tot_plusHalos = (ie_tot+(nprocx-1)*2*nboundlines );
  je_tot_plusHalos = (je_tot+(nprocy-1)*2*nboundlines);
  ntracer = ntracer_perHandler * nGCLHandlers;

// initialize the arrays of data
  if( fTest ) {
    fData = InitBuffers(0);
  }
  if(gclTest) {
    if( GenerateRandoms )
      gclData = InitBuffers(fData);
    else
      gclData = InitBuffers(0);
  }

  if(my_cart_id==0)
    std::cout << "Setup GCL Handlers" << std::endl;
// Setup the gcl Handlers
  if(gclTest)
    gclHandler.Init(nGCLHandlers, ntracer_perHandler, nbl_exchg, istart, jstart, iend, jend, ie, je, ke, icomm_cart, lperi_x, lperi_y);

  if(my_cart_id==0)
  std::cout << "Do fortran exchange" << std::endl;
// do fortran exchange of data
  if(fTest)
    FortranExchange(fData);

  std::cout << "Do fortran exchange end" << std::endl;
  MPI_Barrier( icomm_cart);
  if( my_cart_id==0)
    std::cout << "Do GCL exchange" << std::endl;
// do gcl exchange of data
  if( gclTest) {
    GCLExchange(gclData);
  }
  std::cout << "Do GCL exchange end" << std::endl;

// validate results from the communication
  if(fTest && ! gclTest)
  {
    Validate(fData);
  }
  else if(fTest && gclTest)
  {
    if(GenerateRandoms)
      Validate(gclData, fData);
    else {
      Validate(fData);
      Validate(gclData);
    }
  }
  else if(!fTest && gclTest) {
    Validate(gclData);
  }


  if( ! disableValidationReport ) {
// gather the arrays containing the validation results
    GatherValidatedGrid();

// produce gnuplot output if validation fails
    if(my_cart_id == 0){
      ReportGnuplot();
    }
  }

// gather timing values
  GatherTimings();

#ifdef _COLLECT_STATS_
  CommStatsCollector::Report();
#endif

  MPI_Barrier(icomm_cart);

  // finalize


  if(fData) delete fData;
  if(gclData) delete gclData;
  delete isubposa;
  if(fvalidatedGrid) delete fvalidatedGrid;
  if(my_cart_id == 0 ) {
    if(validationCodes_tot) delete validationCodes_tot;
  }

  MPI_Finalize();
  return(0);

}

template<class TDataType>
void compute_moments( TDataType *data, int nelements, double* mean, double *rms)
{
  double sums = 0;
  double sums2 = 0;
  double avg_, rms_;

  for(int j=0; j < nelements; j++) {
    sums += data[j];
    sums2 += data[j]*data[j];
  }
  avg_ = sums/ ((double)nelements);
  rms_ = std::sqrt( sums2 /((double)nelements) - (avg_)*(avg_) );

  *mean = avg_;
  *rms = rms_;

}
void GatherTimings()
{

  double* etime_gcl = gclHandler.pEtime_gcl();
  double etime_gcl_tot[nproc* comm_nitems ] ;
  double etime_fortran_tot[nproc*comm_nitems];

  int data_transfered_tot[nproc];

  for(int i=0; i < comm_nitems; ++i ) {
    MPI_Gather(&etime_gcl[i], 1, MPI_DOUBLE,  &etime_gcl_tot[i*nproc], 1, MPI_DOUBLE, 0, icomm_cart);
    MPI_Gather(&etime_fortran[i], 1, MPI_DOUBLE,  &etime_fortran_tot[i*nproc], 1, MPI_DOUBLE, 0, icomm_cart);
  }
  MPI_Gather(&data_transfered_perHandler, 1, MPI_INT, data_transfered_tot, 1, MPI_INT, 0, icomm_cart);

  if( my_cart_id == 0) {

    double avg_gcl_times[comm_nitems];
    double rms_gcl_times[comm_nitems];
    double avg_fortran_times[comm_nitems];
    double rms_fortran_times[comm_nitems];

    double avg_data_transfered, rms_data_transfered;

    for(int i=0; i < comm_nitems; ++i) {

      compute_moments<double>( &etime_gcl_tot[i*nproc], nproc, &avg_gcl_times[i], &rms_gcl_times[i]);
      compute_moments<double>( &etime_fortran_tot[i*nproc], nproc, &avg_fortran_times[i], &rms_fortran_times[i]);
    }
    compute_moments<int>( data_transfered_tot, nproc, &avg_data_transfered, &rms_data_transfered);

    std::cout << std::endl << std::endl;
    std::cout << "~~~~~~~~~~~~~~~ GCL TIMES ~~~~~~~~~~~" << std::endl;
    std::cout << std::endl << "\t\t\t mean \t rms" << std::endl;
    std::cout << "pack --> \t" <<  avg_gcl_times[ comm_pack ] << " \t " << rms_gcl_times[ comm_pack ] << std::endl;
    std::cout << "send --> \t" <<  avg_gcl_times[ comm_send ] << " \t " << rms_gcl_times[ comm_send ] << std::endl;
    std::cout << "wait --> \t" <<  avg_gcl_times[ comm_wait ] << " \t " << rms_gcl_times[ comm_wait ] << std::endl;
    std::cout << "unpack --> \t" <<  avg_gcl_times[ comm_unpack ] << " \t " << rms_gcl_times[ comm_unpack ] << std::endl;
    std::cout << "Total --> \t" << avg_gcl_times[ comm_total ] << " \t " << rms_gcl_times[ comm_total ] << std::endl;

    std::cout << std::endl;
    std::cout << "~~~~~~~~~~~~~~~ FORTRAN TIMES ~~~~~~~~~~~" << std::endl;
    std::cout << std::endl << "\t\t\t mean \t rms" << std::endl;
    std::cout << "pack --> \t" <<  avg_fortran_times[ comm_pack ] << " \t " << rms_fortran_times[ comm_pack ] << std::endl;
    std::cout << "send --> \t" <<  avg_fortran_times[ comm_send ] << " \t " << rms_fortran_times[ comm_send ] << std::endl;
    std::cout << "wait --> \t" <<  avg_fortran_times[ comm_wait ] << " \t " << rms_fortran_times[ comm_wait ] << std::endl;
    std::cout << "unpack --> \t" <<  avg_fortran_times[ comm_unpack ] << " \t " << rms_fortran_times[ comm_unpack ] << std::endl;
    std::cout << "Total --> \t" <<  avg_fortran_times[ comm_total ] << " \t " << rms_fortran_times[ comm_total ] << std::endl;

    std::cout << std::endl;
    std::cout << "~~~~~~~~~~~~~~ Data Transfered ~~~~~~~~~~~~~~" << std::endl;
    std::cout <<  avg_data_transfered << " bytes (per gcl Handler) * " << nGCLHandlers << " (num of GCL Handlers) " << std::endl;
    std::cout << std::endl << std::endl;

    std::cout << "~~~~~~~~~~~~~~ Bandwidth ~~~~~~~~~~~~~~" << std::endl;
    std::cout << "pack/unpack --> \t" << avg_data_transfered * nGCLHandlers * 2 / (avg_gcl_times[ comm_pack ] + avg_gcl_times[ comm_unpack ] ) / std::pow(1024,3) << "  GB/s " << std::endl;
    std::cout << "exchange --> \t" << avg_data_transfered * nGCLHandlers * 2 / (avg_gcl_times[ comm_send ]+avg_gcl_times[ comm_wait ]) / std::pow(1024,3) << "  GB/s " << std::endl;

  }

}

void ReportGnuplot()
{
  for(int h=0; h < nGCLHandlers; ++h) {
    for (int l=0; l<ntracer_perHandler; l++) {
      int offset_plusHalos = ie_tot_plusHalos*je_tot_plusHalos*ke_tot*(l+h*ntracer_perHandler);
      std::string filename = "validatedGrid_gclHandler";

      std::stringstream ss;
      ss << h;
      filename = filename + ss.str();
      filename = filename +"_tracer";
      ss.str("");
      ss << l;
      filename = filename + ss.str();
      plotter.Init( 10, ie_tot_plusHalos, je_tot_plusHalos, ke_tot, filename, 0, ke_tot-1, 1,1);
      plotter.ConvertValidatedMatrix( &validationCodes_tot[offset_plusHalos] );
      std::cout << "GCL Handler id: " << h << " TRACER id : " << l << std::endl;
      std::cout << plotter.getResult() << std::endl;
    }
  }


}

void GatherValidatedGrid()
{
  int offset_max=0;
  int offset_tot_plusHalos = 0;

  double gathered_field[ ie_max*je_max*nproc];
  double* temp_valid = new double[ie_max*je_max*ke*ntracer];
  for( int k=0; k<ke; ++k )
  {
    for(int l=0; l < ntracer; ++l)
    {
     for(int j=0; j<je_max; ++j)
     {
       int offset = ie*je*(k+ke*l);
       int offset_tmp = ie_max*je_max*(k+ke*l);
       for(int i=0; i < ie_max; ++i)
       {
         if( i<ie && j<je ){
           temp_valid[ i + ie_max*j + offset_tmp] = fvalidatedGrid[i+ie*j + offset];
         }
// next value should never be loaded in gathered array
// just set a negative value to indicate a problem if this is collected
         else
           temp_valid[i + ie_max*j + offset_tmp] = -2;
       }
     }
   }
  }



  for(int l=0; l < ntracer; ++l) {
    for(int k=0; k < ke; ++k)
    {
      offset_max = ie_max*je_max*( k + ke*l);
      offset_tot_plusHalos = ie_tot_plusHalos*je_tot_plusHalos*(k+ke*l);

      MPI_Gather(&temp_valid[offset_max], ie_max*je_max, MPI_DOUBLE, gathered_field, ie_max*je_max, MPI_DOUBLE, 0, icomm_cart);
      if( my_cart_id == 0 )
        combine_subarrays( gathered_field, &validationCodes_tot[ offset_tot_plusHalos ], true );
    }
  }

  delete temp_valid;
}

void combine_subarrays( double* gathered_field, double* dest, bool includeHalos )
{

  int iz_is, iz_js, iz_ie, iz_je, ig_offset, jg_offset;
  int dest_offset_i = 0;
  int dest_offset_j = 0;
  int nproc_column, nproc_row;
  int iStride;

  if(includeHalos) iStride = ie_tot_plusHalos;
  else iStride = ie_tot;

  for(int np=0; np < nproc; ++np) {
    if(includeHalos) {
      nproc_column = np / nprocy;
      nproc_row = np % nprocy;
      iz_is = isubpos[0][np] - nboundlines -1;
      iz_js = isubpos[1][np] - nboundlines -1;
      iz_ie = isubpos[2][np] + nboundlines -1;
      iz_je = isubpos[3][np] + nboundlines -1;
      if(nproc_column == 0)
        dest_offset_i = 0;
      else
        dest_offset_i = nboundlines*2*nproc_column;

       if(nproc_row == 0)
        dest_offset_j = 0;
       else
        dest_offset_j = nboundlines*2*nproc_row;

    }
    else {
      if( np < nprocy || includeHalos )
        iz_is = isubpos[0][np] - nboundlines -1;
      else
        iz_is = isubpos[0][np] -1;
      if( np % nprocy == 0 || includeHalos)
        iz_js = isubpos[1][np] - nboundlines -1;
      else
        iz_js = isubpos[1][np] -1;

      if( ((np >= nproc - nprocy) && (np <= nproc-1)) || includeHalos )
        iz_ie = isubpos[2][np] + nboundlines -1;
      else
        iz_ie = isubpos[2][np] -1;

      if( ((np+1)%nprocy == 0) || includeHalos)
        iz_je = isubpos[3][np] + nboundlines -1;
      else
        iz_je = isubpos[3][np] -1;
    }

    ig_offset = -isubpos[0][np] + nboundlines+1;
    jg_offset = -isubpos[1][np] + nboundlines+1;


    for(int i=iz_is; i <=  iz_ie; ++i) {
      for(int j=iz_js; j <= iz_je; ++j) {

        if( gathered_field[ (i+ig_offset)+ie_max*((j+jg_offset)+je_max*(np))  ] == 0 ) {
          if( includeHalos  && ( (i - iz_is ) < nboundlines || (iz_ie - i) < nboundlines ||
                             (j - iz_js) < nboundlines || (iz_je -j) < nboundlines ) )
             dest[(j+dest_offset_j) * iStride + (i+dest_offset_i)] = 1;
          else
             dest[(j+dest_offset_j) * iStride + (i+dest_offset_i)] = 0;


        }
        else {
          if( includeHalos  && ( (i - iz_is ) < nboundlines || (iz_ie - i) < nboundlines ||
                             (j - iz_js) < nboundlines || (iz_je -j) < nboundlines ) )
            dest[(j+dest_offset_j) * iStride + (i+dest_offset_i)] = 6;
          else
            dest[(j+dest_offset_j) * iStride + (i+dest_offset_i)] = 5;
        }
      }
    }

  }

}


// function that checks correctness of the communication outcome
// It can check Fortran against expected, GCL agains expected or
// Fortran against GCL (depending on wether fortran and GCL were enabled).
void Validate( double* result, double* expected)
{

  int id,i,j,k,l,ii,ip1,jp1,icount;
  double b,diff;

  // ************************************
  // check correctness

  icount = 0;
  diff = -99;

  for (id = 0; id<nprocx*nprocy; id++) {
    MPI_Barrier(icomm_cart);
    if (id == my_cart_id) {

      printf("***************\n");
      printf(" task %i took %6.3f seconds\n", my_world_id, etime_fortran[comm_total]);
      icount=0;
      for (l=0; l<ntracer; l++) {
        for (k=0; k<ke; k++) {
          for (j=0; j<je; j++) {
            ii = ie*(j+je*(k+ke*l));
            for (i=0; i<ie; i++) {
              ip1 = i+1;
              jp1 = j+1;
              // b = llkkjjjiii

              b = ((double)i_global(&ip1))*1.0E0L +
                  ((double)j_global(&jp1))*1.0E3L +
                  ((double)(k+1)         )*1.0E6L +
                  ((double)(l+1)         )*1.0E8L;
              b *= sign;

// switch sign if at the right boundary
              if( my_cart_neigh[2] != -1 && !lperi_x )
                if( ( ( i > ie - nboundlines -1)  && (i < ie - nboundlines + nbl_exchg) ) &&
                    ( (j >  nboundlines -1 ) && (j < je - nboundlines ) ) )
                    b *= -1;

// switch sign if at the bottom boundary
              if( my_cart_neigh[1] != -1 && !lperi_y)
                if( (( j > je - nboundlines - 1) && (j < je - nboundlines + nbl_exchg ) ) &&
                   ( ( i > nboundlines - 1) && ( i < ie - nboundlines) ) )
                    b *= -1;

// switch sign if at the left boundary
              if( my_cart_neigh[0] != -1 && !lperi_x)
                if( ( (i < nboundlines ) && ( i > nboundlines - nbl_exchg -1 ) ) &&
                  ( (j >  nboundlines -1 ) && (j < je - nboundlines ) ) )
                    b *= -1;

// switch sign if at the top boundary
              if( my_cart_neigh[3] != -1 && !lperi_y)
                if( ( (j < nboundlines ) && ( j > nboundlines - nbl_exchg -1 ) ) &&
                  ( ( i > nboundlines - 1) && ( i < ie - nboundlines) ) )
                    b *= -1;

              if( expected != NULL) {
                diff = fabs( result[ii] - expected[ii] );
              }
              else {
                diff = fabs( result[ii] - b );
              }
              if ( GenerateRandoms ||
                  ! (  (lperi_x && ( i < (nboundlines + nbl_exchg) || i > ie -1 - ( nboundlines + nbl_exchg) )) ||
                   (lperi_y && ( j < (nboundlines + nbl_exchg) || j > je -1 - (nboundlines + nbl_exchg) ) ) ) ) {

                if (diff > 0.1L ) {

                  icount++;
                  printf("id: %d, i: %d, j: %d  sign: %d   ,  ie-nbl: %d,   je-nbl: %d,  ii: %d\n ",my_cart_id,i,j,sign,ie-nbl_exchg,je-nbl_exchg,ii);
                  if(gclTest && !fTest)
                    printf("target=%12.0f actual=%12.0f linear=%7i position=%2.2i%2.2i%3.3i%3.3i  i:%d, j: %d \n",b,gclData[ii],ii+1,l+1,k+1,j+1,i+1,i,j);
                  else if(fTest && !gclTest)
                    printf("target=%12.0f actual=%12.0f linear=%7i position=%2.2i%2.2i%3.3i%3.3i  i:%d, j: %d \n",b,fData[ii],ii+1,l+1,k+1,j+1,i+1,i,j);
                  else
                    printf("ntracer=%d Fdata=%12.0f GCLdata=%12.0f linear=%7i position=%2.2i%2.2i%3.3i%3.3i  i:%d, j: %d \n",l,fData[ii],gclData[ii],ii+1,l+1,k+1,j+1,i+1,i,j);
                  fvalidatedGrid[ii] = 1;
                }
              }
              ii++;
            }
          }
        } //end k-loop
      } // end tracer-loop
      if (icount) {
        printf("ERROR: %d failed correctness check\n",id);
      }
      else
        printf("ID: %d validated \n",id);
    }
  }


}

void BindDevice() {
#ifdef _GCL_GPU_
  int local_rank, num_local_procs;
  int dev_count, use_dev_count, my_dev_id;
  char *str;

  /* GET THE MPI RANK */
  if( ( str = getenv ( "MV2_COMM_WORLD_LOCAL_RANK") ) != NULL ) {
      local_rank = atoi ( str );
      printf( "MV2_COMM_WORLD_LOCAL_RANK %s\n", str );
  }
  if( ( str = getenv ("MPISPAWN_LOCAL_NPROCS") ) != NULL ) {
      num_local_procs = atoi( str );
      printf( "MPISPAWN_LOCAL_NPROCS %s\n", str );
  }
  /* NUMBER OF CUDA DEVICES AVAILABLE */
  cudaGetDeviceCount( &dev_count );
  // NUM_GPU_DEVICES allows to explicitly select the maximum
  // number of devices to use
  if( ( str = getenv ("SLURM_NPROCS") ) != NULL ) {
     use_dev_count = atoi( str );
     printf( "SLURM_NPROCS %s\n", str );
  } else {
     use_dev_count = dev_count;
  }
  /* ASSIGN DEVICE ID TO PROCESS */
  my_dev_id = local_rank % use_dev_count;
  printf( "local rank = %d dev id = %d\n", local_rank, my_dev_id );
  cudaSetDevice( my_dev_id );
#endif
}


void init_environment(int argc, char **argv)
{

  int ierror;

#ifdef _GCL_GPU_
  BindDevice();
#endif

  MPI_Init(&argc, &argv);


  // setup world
  icomm_world = MPI_COMM_WORLD;

  MPI_Comm_size(icomm_world, &nproc);

  if( nproc != nprocx*nprocy) {
    printf("Error nproc*nprocy does not match number of PE launched by mpi\n");
    ierror=-9;
    MPI_Abort (MPI_COMM_WORLD,ierror);
  }

  MPI_Comm_rank(icomm_world, &my_world_id);

  MPI_Comm_group(icomm_world, &igroup_world);

  // setup data types
  imp_reals   =  MPI_DOUBLE_PRECISION;
  imp_integers=  MPI_INTEGER;

  Fimp_reals = (int) MPI_Type_c2f( imp_reals );
  Fimp_integers = (int) MPI_Type_c2f( imp_integers );


  fData = NULL;
  gclData = NULL;
  fvalidatedGrid = NULL;
  validationCodes_tot = NULL;

  init_procgrid();


  if ( gclTest ) {
    gcl_setup_cart(icomm_cart);
  }


  // determine domain decomposition (domain_decomposition in src_setup.f90)
  //   isubpos, ie, je, ke, ie_max, je_max
  domain_decomposition();

  // determine local domain constants (grid_constants in src_setup.f90)
  //   istart, jstart, iend, jend, istartpar, jstartpar, iendpar, jendpar
  grid_constants();

  // call init_par_utilities in parallel_utilities.f90
  init_par_utilities(&nprocx, &nprocy, isubposa, &nboundlines, &my_cart_id);

// determine sign on this PE
  if( (my_cart_pos[0] + my_cart_pos[1])%2 ) {
    sign = -1;
  } else {
    sign =  1;
  }

  if( gclTest)
    gcl_pid(&my_cart_id);

}

void init_procgrid()
{
  int nznumdims,  ij[2];

  // create Cartesian processor grid
  nznumdims = 3;

  int ierror;
  // Instead of creating a cartesian communicator here, we call the fortran subroutine that
  // returns a int cartComm and convert it here, to test the MPI_Comm_f2c feature,
  // and make sure that we get a valid comm when passed from fortran

  create_communicator( &nprocx,&nprocy,&lperi_x, &lperi_y, &Ficomm_cart );

  icomm_cart = MPI_Comm_f2c(Ficomm_cart);

  // determine grid rank and position
  MPI_Comm_rank(icomm_cart, &my_cart_id);
  MPI_Cart_coords(icomm_cart, my_cart_id, nznumdims, my_cart_pos);

  // determine left neighbour
  if (my_cart_pos[0] == 0) {
    if (lperi_x) {
      ij[0] = nprocx-1;
      ij[1] = my_cart_pos[1];
      MPI_Cart_rank(icomm_cart, ij, &my_cart_neigh[0]);
    } else {
      my_cart_neigh[0] = -1;
    }
  } else {
    ij[0] = my_cart_pos[0]-1;
    ij[1] = my_cart_pos[1];
    MPI_Cart_rank(icomm_cart, ij, &my_cart_neigh[0]);
  }

  // determine right neighbour
  if (my_cart_pos[0] == nprocx-1) {
    if (lperi_x) {
      ij[0] = 0;
      ij[1] = my_cart_pos[1];
      MPI_Cart_rank(icomm_cart, ij, &my_cart_neigh[2]);
    } else {
      my_cart_neigh[2] = -1;
    }
  } else {
    ij[0] = my_cart_pos[0]+1;
    ij[1] = my_cart_pos[1];
    MPI_Cart_rank(icomm_cart, ij, &my_cart_neigh[2]);
  }

  // determine top neighbour
  if (my_cart_pos[1] == nprocy-1) {
    if (lperi_y) {
      ij[0] = my_cart_pos[0];
      ij[1] = 0;
      MPI_Cart_rank(icomm_cart, ij, &my_cart_neigh[1]);
    } else {
      my_cart_neigh[1] = -1;
    }
  } else {
    ij[0] = my_cart_pos[0];
    ij[1] = my_cart_pos[1]+1;
    MPI_Cart_rank(icomm_cart, ij, &my_cart_neigh[1]);
  }

  // determine bottom neighbour
  if (my_cart_pos[1] == 0) {
    if (lperi_y) {
      ij[0] = my_cart_pos[0];
      ij[1] = nprocy-1;
      MPI_Cart_rank(icomm_cart, ij, &my_cart_neigh[3]);
    } else {
      my_cart_neigh[3] = -1;
    }
  } else {
    ij[0] = my_cart_pos[0];
    ij[1] = my_cart_pos[1]-1;
    MPI_Cart_rank(icomm_cart, ij, &my_cart_neigh[3]);
  }


  // we need to pass the cart comm back to fortran again
  // here we also test MPI_Comm_c2f
  Ficomm_cart = (int) MPI_Comm_c2f(icomm_cart);
}

void domain_decomposition()
{
  int nzcompi, nzcompj, nzsubi, nzsubj, nzix2, nzix1, nzjy2, nzjy1,
      nzix2left, nzjy2lower, ix, iy, nz1d;

  // precompute sizes of local grids
  nzcompi = ie_tot - 2*nboundlines;
  nzcompj = je_tot - 2*nboundlines;

  nzsubi  = nzcompi / nprocx;
  nzsubj  = nzcompj / nprocy;

  nzix2   = nzcompi - nprocx * nzsubi;
  nzix1   = nprocx - nzix2;

  nzjy2   = nzcompj - nprocy * nzsubj;
  nzjy1   = nprocy - nzjy2;

  nzix2left  = nzix1 / 2;
  nzjy2lower = nzjy1 / 2;

  // allocate memory
  isubposa = (int *) malloc(nprocx*nprocy*4*sizeof(int));
  if (!isubposa) {
    printf("ERROR: could not allocate memory for isubposa\n");
    exit(1);
  }
  isubpos[0] = &isubposa[0*nprocx*nprocy];
  isubpos[1] = &isubposa[1*nprocx*nprocy];
  isubpos[2] = &isubposa[2*nprocx*nprocy];
  isubpos[3] = &isubposa[3*nprocx*nprocy];

  // determine position of sub-domain on global grid
  for (ix=0; ix<nprocx; ix++) {
    for (iy=0; iy<nprocy; iy++) {

      nz1d = ix*nprocy+iy;

      if ( (0 <= iy) && (iy <= nzjy2lower-1) ) {
        isubpos[1][nz1d] =  iy    *  nzsubj + nboundlines + 1;
        isubpos[3][nz1d] = (iy+1) *  nzsubj + nboundlines;
      } else if ( (nzjy2lower <= iy) && (iy <= nzjy2lower+nzjy2-1) ) {
        isubpos[1][nz1d] =  iy    * (nzsubj+1) - nzjy2lower + nboundlines + 1;
        isubpos[3][nz1d] = (iy+1) * (nzsubj+1) - nzjy2lower + nboundlines;
      } else if ( (nzjy2lower+nzjy2 <= iy) && (iy <= nprocy-1) ) {
        isubpos[1][nz1d] =  iy    *  nzsubj + nzjy2 + nboundlines + 1;
        isubpos[3][nz1d] = (iy+1) *  nzsubj + nzjy2 + nboundlines;
      }

      if ( (0 <= ix) && (ix <= nzix2left-1) ) {
        isubpos[0][nz1d] =  ix    *  nzsubi + nboundlines + 1;
        isubpos[2][nz1d] = (ix+1) *  nzsubi + nboundlines;
      } else if ( (nzix2left <= ix) && (ix <= nzix2left+nzix2-1) ) {
        isubpos[0][nz1d] =  ix    * (nzsubi+1) - nzix2left + nboundlines + 1;
        isubpos[2][nz1d] = (ix+1) * (nzsubi+1) - nzix2left + nboundlines;
      } else if ( (nzix2left+nzix2 <= ix) && (ix <= nprocx-1) ) {
        isubpos[0][nz1d] =  ix    *  nzsubi + nzix2 + nboundlines + 1;
        isubpos[2][nz1d] = (ix+1) *  nzsubi + nzix2 + nboundlines;
      }

    }
  }

  // compute sub-domain size
  ie = isubpos[2][my_cart_id] - isubpos[0][my_cart_id] + 1 + 2*nboundlines;
  je = isubpos[3][my_cart_id] - isubpos[1][my_cart_id] + 1 + 2*nboundlines;
  ke = ke_tot;

  // determine maximum local domain sizes
  ie_max = (nzsubi+1) + 2*nboundlines;
  je_max = (nzsubj+1) + 2*nboundlines;

}

void grid_constants()
{

  istart   =  1 + nboundlines;
  jstart   =  1 + nboundlines;
  iend     = ie - nboundlines;
  jend     = je - nboundlines;


// Here the logic to set istartpar & jstartpar changes differs from COSMO. Up-down exchanges
// do not exchange a full horizontal line if there is no neighbour at the left or at the right.
// right-left exchanges do not exchange a full vertical line if there is not neighbour at the top
// or at the bottom. Modifications of istartpar has no effect, because it is not considered in
// halo_exchange. For this, variables in halo_exchange were also modified.
//
//  if (my_cart_neigh[0] == -1 || (lperi_x && (my_cart_pos[0] == 0) ) ) {
//    istartpar = 1;
  if (my_cart_neigh[0] == -1 ) {
    if( lperi_x && my_cart_pos[0] == 0) {
      istartpar = 1;
    }
    else {
      istartpar = 1 + nboundlines;
    }
  } else {
    istartpar = 1 + nboundlines;
  }

//  if (my_cart_neigh[2] == -1 || (lperi_x && (my_cart_pos[0] == nprocx-1) ) ) {
//    iendpar   = ie;
  if (my_cart_neigh[2] == -1) {
    if(lperi_x && my_cart_pos[0] == nprocx-1) {
      iendpar = ie;
    }
    else {
      iendpar = ie - nboundlines;
    }
  } else {
    iendpar   = ie - nboundlines;
  }

//  if (my_cart_neigh[1] == -1 || (lperi_y && (my_cart_pos[1] == nprocy-1) ) ) {
//    jendpar   = je;

  if (my_cart_neigh[1] == -1 ) {
    if(lperi_y && my_cart_pos[1] == nprocy-1) {
      jendpar = je;
    }
    else {
      jendpar = je - nboundlines;
    }
  } else {
    jendpar   = je - nboundlines;
  }

//  if (my_cart_neigh[3] == -1 || (lperi_y && (my_cart_pos[1] == 0) ) ) {
//    jstartpar = 1;
  if(my_cart_neigh[3] == -1 ) {
    if(lperi_y && my_cart_pos[1] == 0) {
      jstartpar = 1;
    }
    else {
      jstartpar = 1 + nboundlines;
    }
  } else {
    jstartpar = 1 + nboundlines;
  }

}

void GCLExchange(double *data) {


  int ntstep;

  gcl_pid(&my_cart_id);

  double* exchgData;

#ifdef _GCL_GPU_
    cudaMalloc( &exchgData, matrix_size );
    cudaMemcpy(  exchgData, data, matrix_size, cudaMemcpyHostToDevice);
#else
    exchgData = data;
#endif


  /* barrier and start timing */
  MPI_Barrier(icomm_world);

  data_transfered_perHandler = 0;
  if( lperi_x || my_cart_neigh[0] != -1 )
    data_transfered_perHandler += nbl_exchg*(je - nboundlines*2);
  if( lperi_x || my_cart_neigh[2] != -1 )
    data_transfered_perHandler += nbl_exchg*(je - nboundlines*2);
  if( lperi_y || my_cart_neigh[1] != -1 )
    data_transfered_perHandler += nbl_exchg*(ie);
  if( lperi_y || my_cart_neigh[3] != -1 )
    data_transfered_perHandler += nbl_exchg*(ie);

  data_transfered_perHandler *= ke*ntracer_perHandler*8;

  // since mvapich2 2.0 there is a dynamic initialization of buffers in the GPU, that happen
  // the first time a user call to mpi happen (instead of in the MPI_Init)
  // This mean that the timers of the first event are expected to show bad performance numbers.
  // We run once and reset the timers to avoid measuring this effect
  gclHandler.DoExchange(exchgData);
  gclHandler.ResetTimers();

  /* loop over number of steps */
  for (ntstep=1; ntstep <= nstop; ntstep++) {

    gclHandler.DoExchange(exchgData);
  }

#ifdef _GCL_GPU_
  cudaMemcpy( data, exchgData, matrix_size, cudaMemcpyDeviceToHost);
#endif


}

double* InitBuffers(double *sourceBuff) {

  double *data;
  double b;
  int id, l,k,j,i,ii,ip1,jp1;

  // ************************************
  // allocate and initialize memory

  int asize = ie*je*ke*ntracer;
  int asize_tot_plusHalos = ie_tot_plusHalos*je_tot_plusHalos*ke_tot*ntracer;

  matrix_size = asize*sizeof(double);
  data = (double*) malloc( matrix_size );
  if (!data) {
    printf("ERROR: could not allocate memory for data\n");
    exit(1);
  }
  for( l=0; l<asize; l++) {
      data[l]=-1.0E30L;
  }
  if( !fvalidatedGrid ) {
    fvalidatedGrid = (double*) malloc( asize*sizeof(double) );
    if( !fvalidatedGrid ) {
      std::cout << "ERROR: could not allocate memory for validation grid" << std::endl;
      exit(1);
    }
    for (l=0; l<asize; l++) {
      fvalidatedGrid[l] = 0;
    }
  }

  if( my_cart_id ==0) {
    if( !validationCodes_tot && ! disableValidationReport ) {
      validationCodes_tot = (double*) malloc( asize_tot_plusHalos *sizeof(double) );
      if( !validationCodes_tot) {
        std::cout << "ERROR: could not allocate memory for gathered validation grid" << std::endl;
        exit(1);
      }
      for(int h=0; h<asize_tot_plusHalos; ++h) {
        validationCodes_tot[h] = -1.0E30L;
      }
    }
  }

  // ************************************
  // setup matrix for unit test

  for (id = 0; id<nprocx*nprocy; id++) {
    MPI_Barrier(icomm_cart);
    if (id == my_cart_id) {
      printf("***************\n");
      printf(" setup on task = %5i\n",id);
      for (l=0; l<ntracer; l++) {
        for (k=0; k<ke; k++) {
          for (j=0; j<je; j++) {
            ii = ie*(j+je*(k+ke*l));
            for (i=0; i<ie; i++) {
              ip1 = i+1;
              jp1 = j+1;
              // b = llkkjjjiii
              b = sign*( ((double)i_global(&ip1))*1.0E0L +
                  ((double)j_global(&jp1))*1.0E3L +
                  ((double)(k+1)         )*1.0E6L +
                  ((double)(l+1)         )*1.0E8L );
              if( sourceBuff != 0)
                data[ii] = sourceBuff[ii];
              else if(GenerateRandoms)
                data[ii] = (double) rand();
              else
                data[ii] = b;
              ii++;
            }
          }
        }
      }
    }
  }
  return data;

}

void update_fortran_times(double fpack_time, double fsend_time, double fwait_time, double funpack_time)
{
  etime_fortran[comm_pack] += fpack_time;
  etime_fortran[comm_send] += fsend_time;
  etime_fortran[comm_wait] += fwait_time;
  etime_fortran[comm_unpack] += funpack_time;
}

void FortranExchange( double* data)  {

  int ntstep,ii;
  double *sendbuf;
  int isendbuflen;
  int ntag;
  int ierror =-1;
  double fpack_time, fsend_time, fwait_time, funpack_time;
  double start;

  for(int i=0; i < comm_nitems; ++i) {
    etime_fortran[i] = 0;
  }

  if (ie_tot/nprocx > je_tot/nprocy) {
    isendbuflen = (ie_tot/nprocx+1+2*nboundlines)*nboundlines*(ke_tot+1)*24;
  } else {
    isendbuflen = (je_tot/nprocy+1+2*nboundlines)*nboundlines*(ke_tot+1)*24;
  }
  sendbuf = (double *) malloc(isendbuflen*sizeof(double)*8);
  if (!sendbuf) {
    printf("ERROR: could not allocate memory for sendbuf\n");
    MPI_Abort (MPI_COMM_WORLD,ierror);
  }


  // ************************************
  // init communications

  MPI_Barrier(icomm_cart);

  ntag = 3000;

  start = MPI_Wtime();

  halo_exchange(sendbuf, &isendbuflen, &Fimp_reals,  &Ficomm_cart, &nproc,
    &ie, &je, &ke, &jstartpar, &jendpar,
    &nbl_exchg, my_cart_neigh,
    &lperi_x, &lperi_y,
    &ntag,
    &fpack_time, &fsend_time, &fwait_time, &funpack_time,
    data );

  // ************************************
  // do halo-exchange

  MPI_Barrier(icomm_world);

  // barrier and start timing
  MPI_Barrier(icomm_world);
  start = MPI_Wtime();
  // loop over number of steps
  for (ntstep=1; ntstep <= nstop; ntstep++) {

    for(int h=0; h < nGCLHandlers; ++h) {
      for (int l=0; l<ntracer_perHandler; l++) {

        ntag = 3000+ntstep;
        ii = ie*je*ke*(l+h*ntracer_perHandler);
        halo_exchange(sendbuf, &isendbuflen, &Fimp_reals, &Ficomm_cart, &nproc,
          &ie, &je, &ke, &jstartpar, &jendpar,
          &nbl_exchg, my_cart_neigh,
          &lperi_x, &lperi_y,
          &ntag,
          &fpack_time, &fsend_time, &fwait_time, &funpack_time,
          &data[ii] );
        update_fortran_times( fpack_time, fsend_time, fwait_time, funpack_time );
      }
    }
  }
  etime_fortran[comm_total] += MPI_Wtime() - start;

  // end timing and barrier
  MPI_Barrier(icomm_world);

  delete sendbuf;

}
