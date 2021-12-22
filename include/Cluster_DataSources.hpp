#pragma once

/* ===== Cluster sources ===== */
#include "Cluster_Data.hpp"

/* ===== C++ ===== */
#include <vector>


/* ===== Utility function for creating a vector of data ===== */
std::vector< Data > CreatePseudoData( const int& aBackgroundCount , const int& aClusterCount , const int& aClusterSize , const double& aClusterScale );

/* ===== Function for loading data from CSV file ===== */
std::vector< Data > LoadCSV( const std::string& aFilename , const double& c_x , const double& c_y );