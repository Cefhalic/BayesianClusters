/* ===== Cluster sources ===== */
#include "BayesianClustering/API.hpp"
#include "BayesianClustering/Cluster.hpp"
#include "BayesianClustering/RoI.hpp"
#include "BayesianClustering/RoIproxy.hpp"
#include "BayesianClustering/Configuration.hpp"

// /* ===== C++ ===== */
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
  Configuration::Instance.FromCommandline( argc , argv );
  Configuration::Instance.SetRBins( 0 , 0 , Configuration::Instance.ClusterR() );
  std::cout << "+------------------------------------+" << std::endl;

  const std::string& lInputFilename = Configuration::Instance.inputFile();
  if( lInputFilename.size() == 0 ) throw std::runtime_error( "No input file specified" ); 
  auto lRoI = LoadLocalizationFile( lInputFilename );

  auto lRoIs = ExtractRoIs( lRoI , Auto  );  
  lRoI.clear();

  std::sort( lRoIs.begin() , lRoIs.end() , []( RoI& a , RoI& b ){ return a.mData.size() < b.mData.size(); } );

  std::cout << "+------------------------------------+" << std::endl;
  for( auto& lRoI : lRoIs )    
  {
    std::cout << "Clusterizing RoI with " << lRoI.mData.size() << " localizations" << std::endl;
    lRoI.Clusterize( 
      Configuration::Instance.ClusterR() , 
      Configuration::Instance.ClusterT() , 
      &ReportClusters
    ); 
    std::cout << "+------------------------------------+" << std::endl;
  }

}
