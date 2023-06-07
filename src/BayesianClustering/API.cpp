//! \file API.cpp

/* ===== Cluster sources ===== */
#include "BayesianClustering/API.hpp"
#include "BayesianClustering/LocalizationFile.hpp"
#include "BayesianClustering/Configuration.hpp"
#include "BayesianClustering/RoI.hpp"
#include "BayesianClustering/RoIproxy.hpp"


// /* ===== Local utilities ===== */
// #include "Utilities/ProgressBar.hpp"
// #include "Utilities/Vectorize.hpp"
// #include "Utilities/Units.hpp"

// /* ===== C++ ===== */
#include <algorithm>
#include <mutex>

void AutoRoi_Scan_FullCallback( const std::string& aInFile , const ScanConfiguration& aScanConfig, const std::function< void( RoIproxy&, const double&, const double& ) >& aCallback )
{
  LocalizationFile( aInFile ).ExtractRoIs( [&]( RoI& aRoI ) { aRoI.ScanRT( aScanConfig, aCallback ); } );
}


void AutoRoi_Scan_SimpleCallback( const std::string& aInFile , const ScanConfiguration& aScanConfig, const std::function< void( const std::vector< ScanEntry >&  ) >& aCallback )
{
  std::mutex lMtx;
  std::vector< ScanEntry > lResults;
  AutoRoi_Scan_FullCallback( aInFile , aScanConfig , [&]( const RoIproxy& aRoI, const double& aR, const double& aT ) { lMtx.lock(); lResults.push_back( { aR, aT, aRoI.mLogP } ); lMtx.unlock(); } );
  std::sort( lResults.begin(), lResults.end() );
  aCallback( lResults );  
}


void ScanCallback_Json( const std::vector< ScanEntry >& aVector, const std::string& aOutFile )
{
  char lOutFileName[256];
  static int RoIid( 0 );  
  sprintf( lOutFileName , "RoI%d.%s", RoIid++ , aOutFile.c_str() );
  FILE *fptr = fopen( lOutFileName , "w" );
  if (fptr == NULL) throw std::runtime_error("Could not open file");
  fprintf( fptr , "[\n" );
  for( auto& lIt : aVector ) fprintf( fptr , "  { \"r\":%.5e , \"t\":%.5e , \"logP\":%.5e },\n" , lIt.r , lIt.t , lIt.score );
  fseek( fptr, -2, SEEK_CUR ); // Delete the last comma
  fprintf( fptr , "\n]\n" );
  fclose(fptr); 
}


void AutoRoi_Scan_ToJson( const std::string& aInFile , const ScanConfiguration& aScanConfig, const std::string& aOutFile )
{
  AutoRoi_Scan_SimpleCallback( aInFile , aScanConfig , [&]( const std::vector< ScanEntry >& aVector ){ ScanCallback_Json( aVector , aOutFile ); } );
}


void AutoRoi_Cluster_Callback( const std::string& aInFile , const double& aR, const double& aT, const std::function< void( RoIproxy& ) >& aCallback )
{  
  LocalizationFile( aInFile ).ExtractRoIs( [&]( RoI& aRoI ) { aRoI.Clusterize( aR, aT, aCallback ); } );
}
