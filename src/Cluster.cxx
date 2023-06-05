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
// \param aRoI The RoI to draw
void ReportClusters( const RoIproxy& aProxy )
{
  std::map< const Cluster* , std::vector< const Data* > > lClusters;

  for( auto& i : aProxy.mData )
  { 
    lClusters[ i.mCluster ? i.mCluster->GetParent() : NULL ].push_back( i.mData );
  }

  std::cout << lClusters.size() << " Clusters" << std::endl;

  for( auto& i : lClusters )
  { 
    if( i.first ) std::cout << " > Cluster of " << i.second.size() << " localizations" << std::endl;
    else          std::cout << " > " << i.second.size() << " background localizations" << std::endl;    
  } 

}


void RoIcallback( RoI& aRoI , const double& aR , const double& aT )
{
    std::cout << "Clusterizing RoI with " << aRoI.mData.size() << " localizations" << std::endl;
    aRoI.Clusterize( aR , aT , &ReportClusters ); 
}


/* ===== Main function ===== */
int main(int argc, char **argv)
{

  std::cout << "+------------------------------------+" << std::endl;
  ProgressBar2 lBar( "| Cluster. Andrew W. Rose. 2022 |" , 1 );
  std::cout << "+------------------------------------+" << std::endl;
  AuxConfiguration lMasterConfig;
  lMasterConfig.FromCommandline( argc , argv );
  std::cout << "+------------------------------------+" << std::endl;

  const std::string& lInputFilename = lMasterConfig.inputFile();
  if( lInputFilename.size() == 0 ) throw std::runtime_error( "No input file specified" ); 
  auto lDataset = LocalizationFile( lInputFilename );

  lDataset.ExtractRoIs( [&]( RoI& aRoI ){ RoIcallback( aRoI , lMasterConfig.ClusterR() , lMasterConfig.ClusterT() ); } );

}
