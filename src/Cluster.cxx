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

/* ===== Local utilities ===== */
#include "Utilities/ProgressBar.hpp"


//! Callback to report clusters
//! \param aProxy The RoI to report
void ReportClusters( RoIproxy& aProxy )
{
  std::map< const Cluster*, std::vector< const Data* > > lClusters;

  for( auto& i : aProxy.mData ) {
    lClusters[ i.mCluster ? i.mCluster->GetParent() : NULL ].push_back( i.mData );
  }

  std::cout << lClusters.size() << " Clusters" << std::endl;

  for( auto& i : lClusters ) {
    if( i.first ) std::cout << " > Cluster of " << i.second.size() << " localizations" << std::endl;
    else          std::cout << " > " << i.second.size() << " background localizations" << std::endl;
  }
}


//! The main function
//! \param argc The number of commandline arguments
//! \param argv The commandline arguments
//! \return     The exit code
int main(int argc, char** argv)
{
  std::cout << "+------------------------------------+" << std::endl;
  ProgressBar2 lBar( "| Cluster. Andrew W. Rose. 2022      |", 1 );
  std::cout << "+------------------------------------+" << std::endl;
  AuxConfiguration lMasterConfig( argc, argv );
  std::cout << "+------------------------------------+" << std::endl;

  AutoRoi_Cluster_Callback( lMasterConfig.inputFile() , lMasterConfig.ClusterR(), lMasterConfig.ClusterT() , &ReportClusters );
}
