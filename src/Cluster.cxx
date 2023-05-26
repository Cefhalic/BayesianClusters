/* ===== Cluster sources ===== */
#include "BayesianClustering/API.hpp"
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


/* ===== Main function ===== */
int main(int argc, char **argv)
{

  std::cout << "+------------------------------------+" << std::endl;
  ProgressBar2 lBar( "| Cluster. Andrew W. Rose. 2022 |" , 1 );
  std::cout << "+------------------------------------+" << std::endl;
  Configuration lMasterConfig;
  lMasterConfig.FromCommandline( argc , argv );
  lMasterConfig.SetSigmaParameters( 0 , 0 , 0 , [&]( const double& ){ return 0.0; } ); 
  std::cout << "+------------------------------------+" << std::endl;

  const std::string& lInputFilename = lMasterConfig.inputFile();
  if( lInputFilename.size() == 0 ) throw std::runtime_error( "No input file specified" ); 
  auto lDataset = LoadLocalizationFile( lInputFilename );

  std::cout << "+------------------------------------+" << std::endl;
  SetCurrentConfiguration( lMasterConfig );
  for( auto& lRoI : ExtractRoIs( lDataset , Auto  ) )    
  {
    std::cout << "+------------------------------------+" << std::endl;
    SetCurrentConfiguration( lRoI.mConfiguration );
    std::cout << "Clusterizing RoI with " << lRoI.mData.size() << " localizations" << std::endl;
    lRoI.Clusterize( 
      CurrentConfiguration().ClusterR() , 
      CurrentConfiguration().ClusterT() , 
      &ReportClusters
    ); 
    lRoI.mData.clear();
  }

}
