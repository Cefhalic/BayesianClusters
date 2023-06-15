//! \file API.cpp

/* ===== Cluster sources ===== */
#include "BayesianClustering/API.hpp"
#include "BayesianClustering/LocalizationFile.hpp"
#include "BayesianClustering/Configuration.hpp"
#include "BayesianClustering/RoI.hpp"
#include "BayesianClustering/RoIproxy.hpp"

/* ===== C++ ===== */
#include <map>
#include <algorithm>
#include <mutex>

/* ===== BOOST C++ ===== */
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/polygon.hpp>


//! A callback to dump a scan to a JSON file
//! \param aVector   A vector of scan results
//! \param aOutFile  The name of the output JSON file
void _ScanCallback_Json_( const std::vector< ScanEntry >& aVector, const std::string& aOutFile )
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

//! A callback to neatly package the scan results for easy consumption
//! \param aRoI        The region of interest
//! \param aScanConfig The configuration for the scan
//! \param aCallback   The simple callback to be applied
void _FullScanToSimpleScan_( RoI& aRoI , const ScanConfiguration& aScanConfig , const std::function< void( const std::vector< ScanEntry >&  ) >& aCallback )
{
  std::mutex lMtx;
  std::vector< ScanEntry > lResults;
  aRoI.ScanRT( aScanConfig, [&]( const RoIproxy& aRoI, const double& aR, const double& aT ) { lMtx.lock(); lResults.push_back( { aR, aT, aRoI.mLogP } ); lMtx.unlock(); } );
  std::sort( lResults.begin(), lResults.end() );
  aCallback( lResults );  
}


//! A callback to neatly package the scan results for easy consumption
//! \param aRoIproxy   The region-proxy containing the clusters
//! \param aCallback   The simple callback to be applied
void _FullClusterToSimpleCluster_( RoIproxy& aRoIproxy , const std::function< void( const std::vector< ClusterWrapper >&  ) >& aCallback )
{
  typedef boost::geometry::model::point<double, 2, boost::geometry::cs::cartesian> geo_point;
  typedef boost::geometry::model::ring<geo_point> geo_polygon;

  std::map< Cluster* , geo_polygon > lMap;
  for( auto& i : aRoIproxy.mData )
  {
    if( ! i.mCluster ) continue;
    boost::geometry::append( lMap[ i.mCluster->GetParent() ] , geo_point( i.mData->x , i.mData->y ) );
  }


  std::vector< ClusterWrapper > lResults;
  for ( auto& i : lMap )
  {
    boost::geometry::correct( i.second );
    geo_polygon lHull;
    boost::geometry::convex_hull( i.second , lHull );   
    geo_point lCentroid ( 0 , 0 );
    boost::geometry::centroid( i.second , lCentroid );    
    lResults.emplace_back( ClusterWrapper{ i.second.size() , boost::geometry::area( lHull ) , boost::geometry::perimeter( lHull ) , boost::geometry::get<0>( lCentroid ) , boost::geometry::get<1>( lCentroid ) } );
  }

  aCallback( lResults );  
}



__attribute__((flatten))
void AutoRoi_Scan_FullCallback( const std::string& aInFile , const ScanConfiguration& aScanConfig, const std::function< void( RoIproxy&, const double&, const double& ) >& aCallback )
{
  LocalizationFile( aInFile ).ExtractRoIs( [&]( RoI& aRoI ) { aRoI.ScanRT( aScanConfig, aCallback ); } );
}

__attribute__((flatten))
void AutoRoi_Scan_SimpleCallback( const std::string& aInFile , const ScanConfiguration& aScanConfig, const std::function< void( const std::vector< ScanEntry >&  ) >& aCallback )
{
  LocalizationFile( aInFile ).ExtractRoIs( [&]( RoI& aRoI ) { _FullScanToSimpleScan_( aRoI , aScanConfig, aCallback ); } );
}

__attribute__((flatten))
void AutoRoi_Scan_ToJson( const std::string& aInFile , const ScanConfiguration& aScanConfig, const std::string& aOutFile )
{
  AutoRoi_Scan_SimpleCallback( aInFile , aScanConfig , [&]( const std::vector< ScanEntry >& aVector ){ _ScanCallback_Json_( aVector , aOutFile ); } );
}



__attribute__((flatten))
void AutoRoi_Cluster_FullCallback( const std::string& aInFile , const double& aR, const double& aT, const std::function< void( RoIproxy& ) >& aCallback )
{  
  LocalizationFile( aInFile ).ExtractRoIs( [&]( RoI& aRoI ) { aRoI.Clusterize( aR, aT, aCallback ); } );
}

__attribute__((flatten))
void AutoRoi_Cluster_SimpleCallback( const std::string& aInFile , const double& aR, const double& aT, const std::function< void( const std::vector< ClusterWrapper >& ) >& aCallback )
{
  AutoRoi_Cluster_FullCallback( aInFile , aR , aT , [&]( RoIproxy& aRoIproxy ){ _FullClusterToSimpleCluster_( aRoIproxy , aCallback ); } );
}



__attribute__((flatten))
void ManualRoi_Scan_FullCallback( const std::string& aInFile , const ManualRoI& aManualRoI , const ScanConfiguration& aScanConfig, const std::function< void( RoIproxy&, const double&, const double& ) >& aCallback )
{
  LocalizationFile( aInFile ).ExtractRoIs( aManualRoI , [&]( RoI& aRoI ) { aRoI.ScanRT( aScanConfig, aCallback ); } );
}

__attribute__((flatten))
void ManualRoi_Scan_SimpleCallback( const std::string& aInFile , const ManualRoI& aManualRoI , const ScanConfiguration& aScanConfig, const std::function< void( const std::vector< ScanEntry >&  ) >& aCallback )
{
  LocalizationFile( aInFile ).ExtractRoIs( aManualRoI , [&]( RoI& aRoI ) { _FullScanToSimpleScan_( aRoI , aScanConfig, aCallback ); } );
}

__attribute__((flatten))
void ManualRoi_Scan_ToJson( const std::string& aInFile , const ManualRoI& aManualRoI , const ScanConfiguration& aScanConfig, const std::string& aOutFile )
{
  ManualRoi_Scan_SimpleCallback( aInFile , aManualRoI , aScanConfig , [&]( const std::vector< ScanEntry >& aVector ){ _ScanCallback_Json_( aVector , aOutFile ); } );
}



__attribute__((flatten))
void ManualRoi_Cluster_FullCallback( const std::string& aInFile , const ManualRoI& aManualRoI , const double& aR, const double& aT, const std::function< void( RoIproxy& ) >& aCallback )
{  
  LocalizationFile( aInFile ).ExtractRoIs( aManualRoI , [&]( RoI& aRoI ) { aRoI.Clusterize( aR, aT, aCallback ); } );
}

__attribute__((flatten))
void ManualRoi_Cluster_SimpleCallback( const std::string& aInFile , const ManualRoI& aManualRoI , const double& aR, const double& aT, const std::function< void( const std::vector< ClusterWrapper >& ) >& aCallback )
{
  ManualRoi_Cluster_FullCallback( aInFile , aManualRoI , aR , aT , [&]( RoIproxy& aRoIproxy ){ _FullClusterToSimpleCluster_( aRoIproxy , aCallback ); } );
}
