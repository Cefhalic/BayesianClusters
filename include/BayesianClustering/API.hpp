//! \file API.hpp
#pragma once

/* ===== C++ ===== */
#include <string>
#include <vector>
#include <functional>

class RoIproxy;
class ScanConfiguration;
class ManualRoI;
class LocalizationFile;

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


// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//! Automatically extract RoI, run scan and apply a full call-back
//! \param aInFile     The name of the localization file
//! \param aScanConfig The configuration for the scan
//! \param aCallback   The full callback to be applied
void AutoRoi_Scan_FullCallback( const std::string& aInFile , const ScanConfiguration& aScanConfig, const tFullScanCallback& aCallback );

//! Automatically extract RoI, run scan and apply a simple call-back
//! \param aInFile     The name of the localization file
//! \param aScanConfig The configuration for the scan
//! \param aCallback   The simple callback to be applied
void AutoRoi_Scan_SimpleCallback( const std::string& aInFile , const ScanConfiguration& aScanConfig, const tSimpleScanCallback& aCallback );

//! Automatically extract RoI, run scan and dump to JSON file
//! \param aInFile     The name of the localization file
//! \param aScanConfig The configuration for the scan
//! \param aOutputPattern    A formattable-string specifying the name of the output JSON files. Substitutable fields are {input} (giving the stem of the input file name) and {roi} (giving the RoI id).
void AutoRoi_Scan_ToJson( const std::string& aInFile , const ScanConfiguration& aScanConfig, const std::string& aOutputPattern );
// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//! Automatically extract RoI, clusterize and apply a full call-back
//! \param aInFile     The name of the localization file
//! \param aR          The R value of the clusterizer
//! \param aT          The T value of the clusterizer
//! \param aCallback   The callback to be applied
void AutoRoi_Cluster_FullCallback( const std::string& aInFile , const double& aR, const double& aT, const tFullClusterCallback& aCallback );

//! Automatically extract RoI, clusterize and apply a full call-back
//! \param aInFile     The name of the localization file
//! \param aR          The R value of the clusterizer
//! \param aT          The T value of the clusterizer
//! \param aCallback   The callback to be applied
void AutoRoi_Cluster_SimpleCallback( const std::string& aInFile , const double& aR, const double& aT, const tSimpleClusterCallback& aCallback );

//! Automatically specify RoI, clusterize and apply a full call-back
//! \param aInFile     The name of the localization file
//! \param aR          The R value of the clusterizer
//! \param aT          The T value of the clusterizer
//! \param aOutputPattern    A formattable-string specifying the name of the output JSON files. Substitutable fields are {input} (giving the stem of the input file name) and {roi} (giving the RoI id).
void AutoRoi_Cluster_ToJson( const std::string& aInFile , const double& aR, const double& aT, const std::string& aOutputPattern );
// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//! Manually specify RoI, run scan and apply a full call-back
//! \param aInFile     The name of the localization file
//! \param aManualRoI  The manually-specified RoI window
//! \param aScanConfig The configuration for the scan
//! \param aCallback   The full callback to be applied
void ManualRoi_Scan_FullCallback( const std::string& aInFile , const ManualRoI& aManualRoI , const ScanConfiguration& aScanConfig, const tFullScanCallback& aCallback );

//! Manually specify RoI, run scan and apply a simple call-back
//! \param aInFile     The name of the localization file
//! \param aManualRoI  The manually-specified RoI window
//! \param aScanConfig The configuration for the scan
//! \param aCallback   The simple callback to be applied
void ManualRoi_Scan_SimpleCallback( const std::string& aInFile , const ManualRoI& aManualRoI , const ScanConfiguration& aScanConfig, const tSimpleScanCallback& aCallback );

//! Manually specify RoI, run scan and dump to JSON file
//! \param aInFile     The name of the localization file
//! \param aManualRoI  The manually-specified RoI window
//! \param aScanConfig The configuration for the scan
//! \param aOutputPattern    A formattable-string specifying the name of the output JSON files. Substitutable fields are {input} (giving the stem of the input file name) and {roi} (giving the RoI id).
void ManualRoi_Scan_ToJson( const std::string& aInFile , const ManualRoI& aManualRoI , const ScanConfiguration& aScanConfig, const std::string& aOutputPattern );
// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//! Manually specify RoI, clusterize and apply a full call-back
//! \param aInFile     The name of the localization file
//! \param aManualRoI  The manually-specified RoI window
//! \param aR          The R value of the clusterizer
//! \param aT          The T value of the clusterizer
//! \param aCallback   The callback to be applied
void ManualRoi_Cluster_FullCallback( const std::string& aInFile , const ManualRoI& aManualRoI , const double& aR, const double& aT, const tFullClusterCallback& aCallback );

//! Manually specify RoI, clusterize and apply a full call-back
//! \param aInFile     The name of the localization file
//! \param aManualRoI  The manually-specified RoI window
//! \param aR          The R value of the clusterizer
//! \param aT          The T value of the clusterizer
//! \param aCallback   The callback to be applied
void ManualRoi_Cluster_SimpleCallback( const std::string& aInFile , const ManualRoI& aManualRoI , const double& aR, const double& aT, const tSimpleClusterCallback& aCallback );

//! Manually specify RoI, clusterize and apply a full call-back
//! \param aInFile     The name of the localization file
//! \param aManualRoI  The manually-specified RoI window
//! \param aR          The R value of the clusterizer
//! \param aT          The T value of the clusterizer
//! \param aOutputPattern    A formattable-string specifying the name of the output JSON files. Substitutable fields are {input} (giving the stem of the input file name) and {roi} (giving the RoI id).
void ManualRoi_Cluster_ToJson( const std::string& aInFile , const ManualRoI& aManualRoI , const double& aR, const double& aT, const std::string& aOutputPattern );
// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//! Extract RoI using an image-map, run scan and apply a full call-back
//! \param aInFile     The name of the localization file
//! \param aImageJfile  The name of an ImageJ RoI file file
//! \param aScale The size of the LSB in the ImageJ file
//! \param aScanConfig The configuration for the scan
//! \param aCallback   The full callback to be applied
void ImageJRoi_Scan_FullCallback( const std::string& aInFile , const std::string& aImageJfile , const double& aScale , const ScanConfiguration& aScanConfig, const tFullScanCallback& aCallback );

//! Extract RoI using an image-map, run scan and apply a simple call-back
//! \param aInFile     The name of the localization file
//! \param aImageJfile  The name of an ImageJ RoI file file
//! \param aScale The size of the LSB in the ImageJ file
//! \param aScanConfig The configuration for the scan
//! \param aCallback   The simple callback to be applied
void ImageJRoi_Scan_SimpleCallback( const std::string& aInFile , const std::string& aImageJfile , const double& aScale , const ScanConfiguration& aScanConfig, const tSimpleScanCallback& aCallback );

//! Extract RoI using an image-map, run scan and dump to JSON file
//! \param aInFile     The name of the localization file
//! \param aImageJfile  The name of an ImageJ RoI file file
//! \param aScale The size of the LSB in the ImageJ file
//! \param aScanConfig The configuration for the scan
//! \param aOutputPattern    A formattable-string specifying the name of the output JSON files. Substitutable fields are {input} (giving the stem of the input file name) and {roi} (giving the RoI id).
void ImageJRoi_Scan_ToJson( const std::string& aInFile , const std::string& aImageJfile , const double& aScale , const ScanConfiguration& aScanConfig, const std::string& aOutputPattern );
// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//! Extract RoI using an image-map, clusterize and apply a full call-back
//! \param aInFile     The name of the localization file
//! \param aImageJfile  The name of an ImageJ RoI file file
//! \param aScale The size of the LSB in the ImageJ file
//! \param aR          The R value of the clusterizer
//! \param aT          The T value of the clusterizer
//! \param aCallback   The callback to be applied
void ImageJRoi_Cluster_FullCallback( const std::string& aInFile , const std::string& aImageJfile , const double& aScale , const double& aR, const double& aT, const tFullClusterCallback& aCallback );

//! Extract RoI using an image-map, clusterize and apply a full call-back
//! \param aInFile     The name of the localization file
//! \param aImageJfile  The name of an ImageJ RoI file file
//! \param aScale The size of the LSB in the ImageJ file
//! \param aR          The R value of the clusterizer
//! \param aT          The T value of the clusterizer
//! \param aCallback   The callback to be applied
void ImageJRoi_Cluster_SimpleCallback( const std::string& aInFile , const std::string& aImageJfile , const double& aScale , const double& aR, const double& aT, const tSimpleClusterCallback& aCallback );

//! Extract RoI using an image-map, clusterize and apply a full call-back
//! \param aInFile     The name of the localization file
//! \param aImageJfile  The name of an ImageJ RoI file file
//! \param aScale The size of the LSB in the ImageJ file
//! \param aR          The R value of the clusterizer
//! \param aT          The T value of the clusterizer
//! \param aOutputPattern    A formattable-string specifying the name of the output JSON files. Substitutable fields are {input} (giving the stem of the input file name) and {roi} (giving the RoI id).
void ImageJRoi_Cluster_ToJson( const std::string& aInFile , const std::string& aImageJfile , const double& aScale , const double& aR, const double& aT, const std::string& aOutputPattern );
// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

// YA
void SegmentedImage_Cluster_FullCallback(const std::string& aInFile, const std::string& aSegmentedImagefile, const double& aScale, const double& aR, const double& aT, const tFullClusterCallback& aCallback);
void SegmentedImage_Cluster_SimpleCallback(const std::string& aInFile, const std::string& aSegmentedImagefile, const double& aScale, const double& aR, const double& aT, const tSimpleClusterCallback& aCallback);
void SegmentedImage_Cluster_ToJson(const std::string& aInFile, const std::string& aSegmentedImagefile, const double& aScale, const double& aR, const double& aT, const std::string& aOutputPattern);

