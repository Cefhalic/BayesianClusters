#pragma once

/* ===== C++ ===== */
#include <vector>
#include <functional>
#include <string>

/* ===== Cluster sources ===== */
#include "BayesianClustering/Data.hpp"
#include "BayesianClustering/Configuration.hpp"

class RoIproxy;


// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//! A class which holds the raw RoI data and global parameters
class RoI
{
public:

  //! A struct for storing a result of an individual scan configuration
  struct ScanEntry
  {
    //! The R parameter
    double r;
    //! The T parameter
    double t;
    //! The score 
    PRECISION score;
  };

  //! Default Constructor
  //! \param aData The set of data-points in the RoI
  RoI( std::vector<Data>&& aData );

  //! Deleted copy constructor
  RoI( const RoI& aOther /*!< Anonymous argument */ ) = delete;

  //! Deleted assignment operator
  //! \return Reference to this, for chaining calls
  RoI& operator= (const RoI& aOther /*!< Anonymous argument */ ) = delete;

  //! Default move constructor
  RoI( RoI&& aOther /*!< Anonymous argument */ ) = default;

  //! Default move-assignment constructor
  //! \return Reference to this, for chaining calls
  RoI& operator= ( RoI&& aOther /*!< Anonymous argument */ ) = default;

  //! All the necessary pre-processing to get the RoI ready for an RT-scan
  //! \param aMaxR The maximum radius out to which we should pre-process
  //! \param aSigmabins2 The number of sigma bins
  void Preprocess( const double& aMaxR , const std::vector< double >& aSigmabins2 );
  
  //! Run the scan
  //! \param aScanConfig      The configuration parameters for the scan
  //! \param aCallback A callback for each RT-scan result
  void ScanRT( const ScanConfiguration& aScanConfig , const std::function< void( const RoIproxy& , const double& , const double& , std::pair<int,int>  ) >& aCallback );

  //! Run the scan
  //! \param aScanConfig      The configuration parameters for the scan
  //! \param aCallback A callback for each RT-scan result
  void ScanRT( const ScanConfiguration& aScanConfig , const std::function< void( const std::vector< ScanEntry >&  ) >& aCallback  );


  //! Run clusterization for a specific choice of R and T
  //! \param R The R parameter for clusterization
  //! \param T The T parameter for clusterization
  //! \param aCallback A callback for the clusterization results
  void Clusterize( const double& R , const double& T , const std::function< void( const RoIproxy& ) >& aCallback );
  
  // //! Save an RoI to a file
  // //! \param aFilename The name of the file to which to save   
  // void WriteCSV( const std::string& aFilename );

public:
  //! Setter for the centre of the scan window
  //! \param aPhysicalCentreX The x-coordinate of the centre of the window in physical units (becomes 0 in algorithm units)
  //! \param aPhysicalCentreY The y-coordinate of the centre of the window in physical units (becomes 0 in algorithm units)
  void SetCentre( const double& aPhysicalCentreX , const double& aPhysicalCentreY );

  //! Setter for the size of the RoI window
  //! \param aWidthX The width of the window in physical units
  //! \param aWidthY The height of the window in physical units
  void SetWidth( const double& aWidthX , const double& aWidthY );
  
  //! Getter for the x-coordinate of the physical centre
  //! \return The x-coordinate of the physical centre
  double getCentreX() const {return mPhysicalCentreX;}
  //! Getter for the y-coordinate of the physical centre
  //! \return The y-coordinate of the physical centre
  double getCentreY() const {return mPhysicalCentreY;}

  //! Getter for the width of the ROI window
  //! \return The width of the ROI window
  double getWidthX() const { return mWidthX; }

  //! Getter for the height of the ROI window
  //! \return The height of the ROI window
  double getWidthY() const { return mWidthY; }

  //! Getter for the height of the ROI window
  //! \return The height of the ROI window
  double getArea() const { return mArea; }  

public:
  //! The collection of raw data points
  std::vector<Data> mData; 

private:
  //! The x-coordinate of the centre of the window in physical units
  double mPhysicalCentreX;
  //! The y-coordinate of the centre of the window in physical units 
  double mPhysicalCentreY;

  //! The width of the window in the x-direction in physical units
  double mWidthX;
  //! The width of the window in the y-direction in physical units
  double mWidthY;
  //! The area of the window in physical units 
  double mArea;  

};
// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
