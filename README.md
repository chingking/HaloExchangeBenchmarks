This document guides through the compilation and execution of the standalone
program for the measurement of the communication performance of GCL.


1. Checkout
===========
The sources can be found in the SVN repository:
https://github.com/cosunae/HaloExchangeBenchmarks.git

The directory contains:
* gcl: a snapshot of the Generic Communication Library
* HaloExchangeTests: the standalone program

Before continuing, make sure the environment is prepared for the compilation.
This includes switching to the desired programming environment (only
PrgEnv-gnu has been tested on Cray machines), loading the cmake module, boost and
the CUDA toolkit.


2. Building the GCL
===================
Enter the gcl directory and create a "build" directory (the exact name matters).
Locate the path to the boost headers.
$ boost=<path-to-boost-headers>

Enter the build directory and compile with

$ export CXX=CC (or any g++ compiler for non Cray machines)

$ cmake  ..  -DCMAKE_BUILD_TYPE=Release  -DBoost_INCLUDE_DIR=$boost  -DGCL_MPI=ON  -DGCL_GPU=ON (-DUSE_MPI_COMPILER=ON for non Cray machines)

$ make

As a result, the library gcl/lib/libgcl.a should exist.


3. Building the standalone executable
=====================================
Enter the HaloExchangeTest directory.  Locate the suitable Options.### file and
create the link:

$ ln  -s  Options.###  Options

Make sure the environment is consistent with the Options file, then compile:

CPU:
$ make cpu

OR

GPU:
$ make gpu

This will create the CommTest executable.


4. Running the standalone executable
====================================
The standalone program provides the following customization switches:

    Flag                    Meaning                                                 Default

    --lperi_x               Pass to enable periodicity in x                         false
    --lperi_y               Pass to enable periodicity in y                         false
    --l2dim                 Pass to do the two-dimensional test                     false

    --nprocx=n              Size of the processes grid in x direction               1
    --nprocy=n              Size of the processes grid in y direction               1

    --ie=n                  Size of the field in x direction                        128
    --je=n                  Size of the field in y directio                         112
    --ke=n                  Size of the field in z direction                        60
    --nblines=n             Number of boundary lines (halo)                         3
    --nbl_exchg=n           Number of boundary lines to exchange (halo)             3

    --enable-Fortran        Enble the fotran test                                   enabled
    --enable-GCL            Enable the GCL test                                     disabled
    --nGCLHandlers=n        Number of GCL handlers                                  1
    --ntracer-perHander=n   How many tracers are exchanged by each GCL handler      50

For the GPU version to work, the environment variable MPICH_RDMA_ENABLED_CUDA
must be exported with the value 1:

$ export MPICH_RDMA_ENABLED_CUDA=1


