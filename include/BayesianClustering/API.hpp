//! \file API.hpp
#pragma once

/* ===== C++ ===== */
#include <string>
#include <vector>
#include <functional>

class RoIproxy;
class ScanConfiguration;
class ManualRoI;

// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//! A struct for storing a result of an individual scan configuration
struct ScanEntry {  
  double r; //!< The R parameter  
  double t; //!< The T parameter
  double score; //!< The score

  //! Comparison operator for sorting
  //! \return Whether we are smaller than the others
  //! \param aOther Another ScanEntry to compare against
  inline bool operator< ( const ScanEntry& aOther )
  {
    if ( r != aOther.r ) return r < aOther.r;
    return ( t < aOther.t );
  }
};
// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//! Automatically extract RoI, run scan and apply a full call-back
//! \param aInFile     The name of the localization file
//! \param aScanConfig The configuration for the scan
//! \param aCallback   The full callback to be applied
void AutoRoi_Scan_FullCallback( const std::string& aInFile , const ScanConfiguration& aScanConfig, const std::function< void( RoIproxy&, const double&, const double& ) >& aCallback );

//! Automatically extract RoI, run scan and apply a simple call-back
//! \param aInFile     The name of the localization file
//! \param aScanConfig The configuration for the scan
//! \param aCallback   The simple callback to be applied
void AutoRoi_Scan_SimpleCallback( const std::string& aInFile , const ScanConfiguration& aScanConfig, const std::function< void( const std::vector< ScanEntry >&  ) >& aCallback );

//! Automatically extract RoI, run scan and dump to JSON file
//! \param aInFile     The name of the localization file
//! \param aScanConfig The configuration for the scan
//! \param aOutFile    The name of the output JSON file
void AutoRoi_Scan_ToJson( const std::string& aInFile , const ScanConfiguration& aScanConfig, const std::string& aOutFile );
// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//! Automatically extract RoI, clusterize and apply a call-back
//! \param aInFile     The name of the localization file
//! \param aR          The R value of the clusterizer
//! \param aT          The T value of the clusterizer
//! \param aCallback   The callback to be applied
void AutoRoi_Cluster_Callback( const std::string& aInFile , const double& aR, const double& aT, const std::function< void( RoIproxy& ) >& aCallback );
// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//! Manually extract RoI, run scan and apply a full call-back
//! \param aInFile     The name of the localization file
//! \param aManualRoI  The manually-specified RoI window
//! \param aScanConfig The configuration for the scan
//! \param aCallback   The full callback to be applied
void ManualRoi_Scan_FullCallback( const std::string& aInFile , const ManualRoI& aManualRoI , const ScanConfiguration& aScanConfig, const std::function< void( RoIproxy&, const double&, const double& ) >& aCallback );

//! Manually extract RoI, run scan and apply a simple call-back
//! \param aInFile     The name of the localization file
//! \param aManualRoI  The manually-specified RoI window
//! \param aScanConfig The configuration for the scan
//! \param aCallback   The simple callback to be applied
void ManualRoi_Scan_SimpleCallback( const std::string& aInFile , const ManualRoI& aManualRoI , const ScanConfiguration& aScanConfig, const std::function< void( const std::vector< ScanEntry >&  ) >& aCallback );

//! Manually extract RoI, run scan and dump to JSON file
//! \param aInFile     The name of the localization file
//! \param aManualRoI  The manually-specified RoI window
//! \param aScanConfig The configuration for the scan
//! \param aOutFile    The name of the output JSON file
void ManualRoi_Scan_ToJson( const std::string& aInFile , const ManualRoI& aManualRoI , const ScanConfiguration& aScanConfig, const std::string& aOutFile );
// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//! Manually extract RoI, clusterize and apply a call-back
//! \param aInFile     The name of the localization file
//! \param aManualRoI  The manually-specified RoI window
//! \param aR          The R value of the clusterizer
//! \param aT          The T value of the clusterizer
//! \param aCallback   The callback to be applied
void ManualRoi_Cluster_Callback( const std::string& aInFile , const ManualRoI& aManualRoI , const double& aR, const double& aT, const std::function< void( RoIproxy& ) >& aCallback );
// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
