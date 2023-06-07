//! \file API.hpp
#pragma once

/* ===== C++ ===== */
#include <string>
#include <vector>
#include <functional>

/* ===== Cluster sources ===== */

// class Data;
class RoIproxy;
class ScanConfiguration;

// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//! A struct for storing a result of an individual scan configuration
struct ScanEntry {
  //! The R parameter
  double r;
  //! The T parameter
  double t;
  //! The score
  double score;

  inline bool operator< ( const ScanEntry& aOther )
  {
    if ( r != aOther.r ) return r < aOther.r;
    return ( t < aOther.t );
  }
};
// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


void AutoRoi_Scan_FullCallback( const std::string& aInFile , const ScanConfiguration& aScanConfig, const std::function< void( RoIproxy&, const double&, const double& ) >& aCallback );

void AutoRoi_Scan_SimpleCallback( const std::string& aInFile , const ScanConfiguration& aScanConfig, const std::function< void( const std::vector< ScanEntry >&  ) >& aCallback );

void AutoRoi_Scan_ToJson( const std::string& aInFile , const ScanConfiguration& aScanConfig, const std::string& aOutFile );

void AutoRoi_Cluster_Callback( const std::string& aInFile , const double& aR, const double& aT, const std::function< void( RoIproxy& ) >& aCallback );
