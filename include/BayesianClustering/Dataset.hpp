#pragma once

/* ===== C++ ===== */
#include <array>
#include <vector>

/* ===== Cluster sources ===== */
#include "BayesianClustering/Data.hpp"

typedef std::array< std::array< std::vector<Data> , 512 > , 512 > Dataset;


