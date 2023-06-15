//! \file Cluster.cxx

/* ===== Cluster sources ===== */
#include "BayesianClustering/API.hpp"
#include "BayesianClustering/Cluster.hpp"
#include "BayesianClustering/RoIproxy.hpp"
#include "BayesianClustering/Configuration.hpp"

// /* ===== C++ ===== */
#include <map>
#include <vector>
#include <iostream>

//! Callback to report clusters
//! \param aClusters A vector of clusters
void ReportClusters( const std::vector< ClusterWrapper >& aClusters )
{
  printf( "+----------------+----------------+----------------+----------------+----------------+\n" );
  printf( "| Localizations  |      Area      |    Perimeter   |   Centroid x   |   Centroid y   |\n" );
  printf( "+----------------+----------------+----------------+----------------+----------------+\n" );
  for( const auto& i : aClusters ) printf( "| %14ld | %14Le | %14Le | %+14e | %+14e |\n", i.localizations , i.area , i.perimeter , i.centroid_x , i.centroid_y );
  printf( "+----------------+----------------+----------------+----------------+----------------+\n" );
}


//! The main function
//! \param argc The number of commandline arguments
//! \param argv The commandline arguments
//! \return     The exit code
int main(int argc, char** argv)
{
  std::cout << "+------------------------------------+" << std::endl;
  std::cout << "| Cluster. Andrew W. Rose. 2022      |" << std::endl;
  std::cout << "+------------------------------------+" << std::endl;
  AuxConfiguration lMasterConfig( argc, argv );
  std::cout << "+------------------------------------+" << std::endl;

  AutoRoi_Cluster_SimpleCallback( lMasterConfig.inputFile() , lMasterConfig.ClusterR(), lMasterConfig.ClusterT() , &ReportClusters );
}
