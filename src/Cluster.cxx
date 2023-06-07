//! \file Cluster.cxx

/* ===== Cluster sources ===== */
#include "BayesianClustering/LocalizationFile.hpp"
#include "BayesianClustering/Cluster.hpp"
#include "BayesianClustering/RoI.hpp"
#include "BayesianClustering/RoIproxy.hpp"
#include "BayesianClustering/Configuration.hpp"

// /* ===== C++ ===== */
#include <map>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>

/* ===== Local utilities ===== */
#include "Utilities/ProgressBar.hpp"
#include "Utilities/ListComprehension.hpp"

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

//! A callback for handling each RoI
//! \param aRoI        The current RoI
//! \param aR          The current R position of the scan
//! \param aT          The current T position of the scan
void RoIcallback( RoI& aRoI, const double& aR, const double& aT )
{
  std::cout << "Clusterizing RoI with " << aRoI.mData.size() << " localizations" << std::endl;
  aRoI.Clusterize( aR, aT, &ReportClusters );
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

  auto lDataset = LocalizationFile( lMasterConfig.inputFile() );
  lDataset.ExtractRoIs( [&]( RoI& aRoI ) { RoIcallback( aRoI, lMasterConfig.ClusterR(), lMasterConfig.ClusterT() ); } );
}
