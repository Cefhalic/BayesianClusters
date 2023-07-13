//! \file API.hpp
#pragma once

/* ===== C++ ===== */
#include <string>
#include <vector>
#include <functional>

/* ===== Cluster sources ===== */
#include "BayesianClustering/LocalizationFile.hpp"
#include "BayesianClustering/RoI.hpp"

class RoIproxy;
class ScanConfiguration;


// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//! A struct for storing a result of an individual scan configuration
struct ScanEntry {  
  double r; //!< The R parameter  
  double t; //!< The T parameter
  double score; //!< The score

  //! Comparison operator for sorting
  //! \return Whether we are smaller than the other
  //! \param aOther Another ScanEntry to compare against
  bool operator< ( const ScanEntry& aOther )
  {
    if ( r != aOther.r ) return r < aOther.r;
    return ( t < aOther.t );
  }

  //! Equality operator required by boost python
  //! \return Whether we are equal to the other
  //! \param aOther Another ScanEntry to compare against
  bool operator== ( const ScanEntry& aOther )
  {
    return ( r == aOther.r ) and ( t == aOther.t ) and ( score == aOther.score );
  }

};
// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//! A struct for storing extracted parameters from a cluster
struct ClusterWrapper {  
  std::size_t localizations; //!< The number of localizations in the cluster
  long double area; //!< The area of the spanning convex hull
  long double perimeter; //!< The perimeter of the spanning convex hull
  double centroid_x; //!< The x-position of the centroid
  double centroid_y; //!< The y-position of the centroid

  //! Comparison operator for sorting
  //! \return Whether we are smaller than the other
  //! \param aOther Another ClusterWrapper to compare against
  bool operator< ( const ClusterWrapper& aOther )
  {
    if ( localizations != aOther.localizations ) return localizations < aOther.localizations;
    if ( area != aOther.area ) return area < aOther.area;
    return ( perimeter < aOther.perimeter );
  }
  
  //! Equality operator required by boost python
  //! \return Whether we are equal to the other
  //! \param aOther Another ClusterWrapper to compare against
  bool operator== ( const ClusterWrapper& aOther )
  {
    return ( centroid_x == aOther.centroid_x ) and ( centroid_y == aOther.centroid_y );
  }  
};
// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

//! Typedef the full scan callback for clarity
typedef std::function< void( RoIproxy&, const double&, const double& ) > tFullScanCallback;
//! Typedef the simplified scan callback for clarity
typedef std::function< void( const std::string& , const std::vector< ScanEntry >&  ) > tSimpleScanCallback;
//! Typedef the full clustering callback for clarity
typedef std::function< void( RoIproxy& ) > tFullClusterCallback;
//! Typedef the simplified clustering callback for clarity
typedef std::function< void( const std::string& , const std::vector< ClusterWrapper >& ) > tSimpleClusterCallback;



//! A callback to dump a scan to a JSON file
//! \param aRoiId       The RoI ID
//! \param aVector   A vector of scan results
//! \param aInFile   The name of the localization file
//! \param aOutputPattern  The name of the output JSON file
void _ScanCallback_Json_( const std::string& aRoiId , const std::vector< ScanEntry >& aVector, const std::string& aInFile , const std::string& aOutputPattern );

// ScanEntry _BestScore_( const ScanConfiguration& aScanConfig , const std::vector< ScanEntry >& aResults );

//! A callback to neatly package the scan results for easy consumption
//! \param aRoI        The region of interest
//! \param aScanConfig The configuration for the scan
//! \param aCallback   The simple callback to be applied
void _FullScanToSimpleScan_( RoI& aRoI , const ScanConfiguration& aScanConfig , const tSimpleScanCallback& aCallback );

//! A callback to dump a clustering run to a JSON file
//! \param aRoiId    The RoI ID
//! \param aVector   A vector of cluster-wrappers
//! \param aInFile   The name of the localization file
//! \param aOutputPattern  The name of the output JSON file
void _ClusterCallback_Json_( const std::string& aRoiId , const std::vector< ClusterWrapper >& aVector, const std::string& aInFile , const std::string& aOutputPattern );

//! A callback to neatly package the scan results for easy consumption
//! \param aRoIproxy   The region-proxy containing the clusters
//! \param aCallback   The simple callback to be applied
void _FullClusterToSimpleCluster_( RoIproxy& aRoIproxy , const tSimpleClusterCallback& aCallback );


// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//! Automatically extract RoI, run scan and apply a full call-back
//! \param aInFile     The name of the localization file
//! \param aScanConfig The configuration for the scan
//! \param aCallback   The full callback to be applied
template< typename RoIConfig >
void RunScan( const std::string& aInFile , const RoIConfig&& aRoIConfig , const ScanConfiguration& aScanConfig, const tFullScanCallback& aCallback )
{
  LocalizationFile( aInFile ).ExtractRoIs( std::move( aRoIConfig ) , [&]( RoI& aRoI ) { aRoI.ScanRT( aScanConfig, aCallback ); } );
}

//! Automatically extract RoI, run scan and apply a simple call-back
//! \param aInFile     The name of the localization file
//! \param aScanConfig The configuration for the scan
//! \param aCallback   The simple callback to be applied
template< typename RoIConfig >
void RunScan( const std::string& aInFile , const RoIConfig&& aRoIConfig , const ScanConfiguration& aScanConfig, const tSimpleScanCallback& aCallback )
{
  LocalizationFile( aInFile ).ExtractRoIs( std::move( aRoIConfig ) , [&]( RoI& aRoI ) { _FullScanToSimpleScan_( aRoI , aScanConfig, aCallback ); } );
}

//! Automatically extract RoI, run scan and dump to JSON file
//! \param aInFile     The name of the localization file
//! \param aScanConfig The configuration for the scan
//! \param aOutputPattern    A formattable-string specifying the name of the output JSON files. Substitutable fields are {input} (giving the stem of the input file name) and {roi} (giving the RoI id).
template< typename RoIConfig >
void RunScan( const std::string& aInFile , const RoIConfig&& aRoIConfig , const ScanConfiguration& aScanConfig, const std::string& aOutputPattern )
{
  RunScan( aInFile , std::move( aRoIConfig ) , aScanConfig , [&]( const std::string& aRoiId , const std::vector< ScanEntry >& aVector ){ _ScanCallback_Json_( aRoiId , aVector , aInFile , aOutputPattern ); } );
}
// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//! Automatically extract RoI, clusterize and apply a full call-back
//! \param aInFile     The name of the localization file
//! \param aR          The R value of the clusterizer
//! \param aT          The T value of the clusterizer
//! \param aCallback   The callback to be applied
template< typename RoIConfig >
void RunClustering( const std::string& aInFile , const RoIConfig&& aRoIConfig , const double& aR, const double& aT, const tFullClusterCallback& aCallback )
{  
  LocalizationFile( aInFile ).ExtractRoIs( std::move( aRoIConfig ) , [&]( RoI& aRoI ) { aRoI.Clusterize( aR, aT, aCallback ); } );
}

//! Automatically extract RoI, clusterize and apply a full call-back
//! \param aInFile     The name of the localization file
//! \param aR          The R value of the clusterizer
//! \param aT          The T value of the clusterizer
//! \param aCallback   The callback to be applied
template< typename RoIConfig >
void RunClustering( const std::string& aInFile , const RoIConfig&& aRoIConfig , const double& aR, const double& aT, const tSimpleClusterCallback& aCallback )
{
  RunClustering( aInFile , std::move( aRoIConfig ) , aR , aT , [&]( RoIproxy& aRoIproxy ){ _FullClusterToSimpleCluster_( aRoIproxy , aCallback ); } );
}

//! Automatically specify RoI, clusterize and apply a full call-back
//! \param aInFile     The name of the localization file
//! \param aR          The R value of the clusterizer
//! \param aT          The T value of the clusterizer
//! \param aOutputPattern    A formattable-string specifying the name of the output JSON files. Substitutable fields are {input} (giving the stem of the input file name) and {roi} (giving the RoI id).
template< typename RoIConfig >
void RunClustering( const std::string& aInFile , const RoIConfig&& aRoIConfig , const double& aR, const double& aT, const std::string& aOutputPattern )
{
  RunClustering( aInFile , std::move( aRoIConfig ) , aR , aT , [&]( const std::string& aRoiId , const std::vector< ClusterWrapper >& aVector ){ _ClusterCallback_Json_( aRoiId , aVector , aInFile , aOutputPattern ); } );
}
// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

