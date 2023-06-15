//! \file MemoryMonitoring.hpp
#pragma once

#include <unistd.h>
#include <ios>
#include <fstream>
#include <string>

//! Utility to get Virtual Memory and Resident Set usage
//! \param vm_usage Return the Virtual Memory usage
//! \param resident_set Return the Resident Set usage
void mem_usage(double& vm_usage, double& resident_set) {
  std::ifstream stat_stream("/proc/self/stat", std::ios_base::in); //get info from proc directory

  //create some variables to get info
  std::string ignore;
  std::size_t vsize, rss;

  stat_stream >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore 
              >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore 
              >> vsize >> rss;
  stat_stream.close();

  vm_usage = vsize / 1024.0;

  double page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024.0; // for x86-64 is configured to use 2MB pages
  resident_set = rss * page_size_kb;
}