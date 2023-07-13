//! \file API.cpp

/* ===== Cluster sources ===== */
#include "BayesianClustering/API.hpp"
#include "BayesianClustering/Configuration.hpp"
#include "BayesianClustering/RoIproxy.hpp"

/* ===== C++ ===== */
#include <map>
#include <algorithm>
#include <mutex>

/* ===== BOOST C++ ===== */
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/filesystem.hpp>

/* ===== FMT library ===== */
#include <fmt/format.h>

//! A callback to dump a scan to a JSON file
//! \param aRoiId       The RoI ID
//! \param aVector   A vector of scan results
//! \param aInFile   The name of the localization file
//! \param aOutputPattern  The name of the output JSON file
void _ScanCallback_Json_( const std::string& aRoiId , const std::vector< ScanEntry >& aVector, const std::string& aInFile , const std::string& aOutputPattern )
{
  using namespace fmt::literals;
  auto lOutFileName = boost::filesystem::path( fmt::format( aOutputPattern , "input"_a = boost::filesystem::path( aInFile ).stem().string() , "roi"_a = aRoiId ) );
  boost::filesystem::create_directories( lOutFileName.parent_path() );

  FILE *fptr = fopen( lOutFileName.c_str() , "w" );
  if (fptr == NULL) throw std::runtime_error("Could not open file");
  fprintf( fptr , "[\n" );
  for( auto& lIt : aVector ) fprintf( fptr , "  { \"r\":%.5e , \"t\":%.5e , \"logP\":%.5e },\n" , lIt.r , lIt.t , lIt.score );
  fseek( fptr, -2, SEEK_CUR ); // Delete the last comma
  fprintf( fptr , "\n]\n" );
  fclose(fptr); 
}


// ScanEntry _BestScore_( const ScanConfiguration& aScanConfig , const std::vector< ScanEntry >& aResults )
// {
//   std::size_t lMax( 0 );
//   for( std::size_t i(1) ; i!= aResults.size() ; ++i )
//     if( aResults[i].score > aResults[lMax].score ) lMax = i;

//   std::size_t R( lMax / aScanConfig.tbounds.bins ), T( lMax % aScanConfig.tbounds.bins );

//   ScanEntry lMean{ 0.0 , 0.0 , 0.0 };
//   for( int r=std::max(0,R-2) ; r!=std::min(aScanConfig.rbounds.bins,R+3) ; ++r ){
//     for( int t=std::max(0,T-2) ; t!=std::min(aScanConfig.tbounds.bins,T+3) ; ++t ){
//       std::size_t lIndex = ( aScanConfig.tbounds.bins * r ) + t;
//       auto& lBin = aResults[lIndex];
//       lMean.r += ( lBin.score * lBin.r );
//       lMean.t += ( lBin.score * lBin.t );
//       lMean.score += lBin.score;
//     }
//   } 

//   lMean.r /= lMean.score;
//   lMean.t /= lMean.score;

//   return lMean;
// }


//! A callback to neatly package the scan results for easy consumption
//! \param aRoI        The region of interest
//! \param aScanConfig The configuration for the scan
//! \param aCallback   The simple callback to be applied
void _FullScanToSimpleScan_( RoI& aRoI , const ScanConfiguration& aScanConfig , const tSimpleScanCallback& aCallback )
{
  std::mutex lMtx;
  std::vector< ScanEntry > lResults;
  aRoI.ScanRT( aScanConfig, [&]( const RoIproxy& aRoI, const double& aR, const double& aT ) { lMtx.lock(); lResults.push_back( { aR, aT, aRoI.mLogP } ); lMtx.unlock(); } );
  std::sort( lResults.begin(), lResults.end() );

  aCallback( aRoI.id() , lResults );  
}

//! A callback to dump a clustering run to a JSON file
//! \param aRoiId    The RoI ID
//! \param aVector   A vector of cluster-wrappers
//! \param aInFile   The name of the localization file
//! \param aOutputPattern  The name of the output JSON file
void _ClusterCallback_Json_( const std::string& aRoiId , const std::vector< ClusterWrapper >& aVector, const std::string& aInFile , const std::string& aOutputPattern )
{
  using namespace fmt::literals;
  auto lOutFileName = boost::filesystem::path( fmt::format( aOutputPattern , "input"_a = boost::filesystem::path( aInFile ).stem().string() , "roi"_a = aRoiId ) );
  boost::filesystem::create_directories( lOutFileName.parent_path() );

  FILE *fptr = fopen( lOutFileName.c_str() , "w" );
  if (fptr == NULL) throw std::runtime_error("Could not open file");
  fprintf( fptr , "[\n" );
  for( auto& lIt : aVector ) fprintf( fptr , "  { \"localizations\":%ld , \"area\":%.5Le , \"perimeter\":%.5Le , \"centroid_x\":%.5e , \"centroid_y\":%.5e },\n" , lIt.localizations , lIt.area , lIt.perimeter , lIt.centroid_x , lIt.centroid_y );
  fseek( fptr, -2, SEEK_CUR ); // Delete the last comma
  fprintf( fptr , "\n]\n" );
  fclose(fptr); 
}


//! A callback to neatly package the scan results for easy consumption
//! \param aRoIproxy   The region-proxy containing the clusters
//! \param aCallback   The simple callback to be applied
void _FullClusterToSimpleCluster_( RoIproxy& aRoIproxy , const tSimpleClusterCallback& aCallback )
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

  std::sort( lResults.begin(), lResults.end() );
  aCallback( aRoIproxy.mRoI.id() , lResults );  
}
