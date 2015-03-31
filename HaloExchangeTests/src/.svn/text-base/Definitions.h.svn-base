#pragma once

#ifdef _GCL_GPU_
    typedef GCL::gcl_gpu arch;
#else
    typedef GCL::gcl_cpu arch;
#endif
#define packing_version GCL::version_manual
//#define generic_pattern
#define dynamic_pattern
#define _COLLECT_STATS_

// communication phases used for timing
enum comm_phases
{
  comm_total =0,
  comm_pack = 1,
  comm_send = 2,
  comm_wait = 3,
  comm_unpack = 4,
  comm_nitems = 5
};

