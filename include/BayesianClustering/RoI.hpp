//! \file RoI.hpp
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
    //! Default Constructor
    //! \param aId The ID of the RoI
    //! \param aData The set of data-points in the RoI
    //! \param aPhysicalCentreX The x-coordinate of the centre of the window in physical units (becomes 0 in algorithm units)
    //! \param aPhysicalCentreY The y-coordinate of the centre of the window in physical units (becomes 0 in algorithm units)
    //! \param aArea The area of the RoI in physical units    
    RoI( const std::string& aId , std::vector<Data>&& aData , const double& aPhysicalCentreX, const double& aPhysicalCentreY , const double& aArea );

    //! Deleted copy constructor
    RoI( const RoI& aOther /*!< Anonymous argument */ ) = delete;

    //! Deleted assignment operator
    //! \return Reference to this, for chaining calls
    RoI& operator= (const RoI& aOther /*!< Anonymous argument */ ) = delete;

    //! Default move constructor
    RoI( RoI&& aOther /*!< Anonymous argument */ ) = default;

    //! Default destructor
    ~RoI();

    //! Default move-assignment constructor
    //! \return Reference to this, for chaining calls
    RoI& operator= ( RoI&& aOther /*!< Anonymous argument */ ) = default;

    //! All the necessary pre-processing to get the RoI ready for an RT-scan
    //! \param aMaxR The maximum radius out to which we should pre-process
    //! \param aSigmabins2 The number of sigma bins
    void Preprocess( const double& aMaxR, const std::vector< double >& aSigmabins2 );

    //! Run the scan
    //! \param aScanConfig      The configuration parameters for the scan
    //! \param aCallback A callback for each RT-scan result
    void ScanRT( const ScanConfiguration& aScanConfig, const std::function< void( RoIproxy&, const double&, const double& ) >& aCallback );

    //! Run clusterization for a specific choice of R and T
    //! \param R The R parameter for clusterization
    //! \param T The T parameter for clusterization
    //! \param aCallback A callback for the clusterization results
    void Clusterize( const double& R, const double& T, const std::function< void( RoIproxy& ) >& aCallback );

    // //! Save an RoI to a file
    // //! \param aFilename The name of the file to which to save
    // void WriteCSV( const std::string& aFilename );

  public:
    //! Getter for the x-coordinate of the physical centre
    //! \return The x-coordinate of the physical centre
    inline const double& getCentreX() const
    {
      return mPhysicalCentreX;
    }

    //! Getter for the y-coordinate of the physical centre
    //! \return The y-coordinate of the physical centre
    inline const double& getCentreY() const
    {
      return mPhysicalCentreY;
    }

    //! Getter for the height of the ROI window
    //! \return The height of the ROI window
    inline const double& getArea() const
    {
      return mArea;
    }

    //! Accessor to the raw data
    //! \return Reference to the raw data
    inline const std::vector< Data >& data() const
    {
      return mData;
    }

    //! Accessor to the RoI ID
    //! \return Reference to the RoI ID
    inline const std::string& id() const
    {
      return mId;
    }

  private:
    friend class RoIproxy;

    //! The ID of the ROI
    std::string mId;

    //! The collection of raw data points
    std::vector<Data> mData;

    //! The x-coordinate of the centre of the window in physical units
    double mPhysicalCentreX;
    //! The y-coordinate of the centre of the window in physical units
    double mPhysicalCentreY;

    //! The area of the window in physical units
    double mArea;

};

// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
