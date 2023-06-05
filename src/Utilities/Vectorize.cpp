//! \file Vectorize.cpp

#include "Utilities/Vectorize.hpp"

#include <thread>

//! The number of threads used, initialized to the number of hardware threads
std::size_t Nthreads = std::thread::hardware_concurrency();
