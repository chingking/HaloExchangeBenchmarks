#pragma once
#include <map>
#include <iostream>


// In order to collect statistics of messages
// one must add the following lines to gcl/L3/include/Halo_Exchange_3D.h
// before the MPI_Isend call 
//
// #ifdef _COLLECT_STATS_
//     CommStatsCollector::AddMessageSize(m_send_buffers.size(I,J,K));
// #endif
//  MPI_Isend(static_cast<char*>(m_send_buffers.buffer(I,J,K)),
//
class CommStatsCollector
{
public:
  CommStatsCollector(){}
  ~CommStatsCollector(){}

  static void AddMessageSize(int size)
  {
    if(messageSizes_.count(size) != 0)
      messageSizes_[size] += 1;
    else
      messageSizes_[size] = 1;
  }

  static void Report()
  {
    std::cout << "Summary of messages sent"  << std::endl;
    std::cout << "Message size (bytes)     |    # count " << std::endl;
    for(std::map<int,int>::const_iterator iter = messageSizes_.begin(); iter != messageSizes_.end(); ++iter)
    {
      std::cout << iter->first << "            " << iter->second << std::endl;
    }
  }

private:

  static std::map<int, int> messageSizes_;
};

