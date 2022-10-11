#include "Utilities/Vectorize.hpp"

#include <thread>

std::size_t Nthreads( std::thread::hardware_concurrency() );
