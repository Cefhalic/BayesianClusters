//! \file Data.hpp
#pragma once

/* ===== C++ ===== */
#include <math.h>
#include <vector>

/* ===== Cluster sources ===== */
#include "BayesianClustering/Precision.hpp"
#include "BayesianClustering/Configuration.hpp"

class Cluster;
class ProgressBar;

// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//! A class to store the raw data-points
class Data
{
  public:
    //! Constructor
    //! \param aX The x-position of the data-point in algorithm units
    //! \param aY The y-position of the data-point in algorithm units
    //! \param aS The sigma of the data-point in algorithm units
    Data( const PRECISION& aX, const PRECISION& aY, const PRECISION& aS );

    //! Deleted copy constructor
    Data( const Data& aOther /*!< Anonymous argument */ ) = delete;

    //! Deleted assignment operator
    //! \return Reference to this, for chaining calls
    Data& operator = (const Data& aOther /*!< Anonymous argument */ ) = delete;

    //! Default move constructor
    Data( Data&& aOther /*!< Anonymous argument */ ) = default;

    //! Default move-assignment constructor
    //! \return Reference to this, for chaining calls
    Data& operator = ( Data&& aOther /*!< Anonymous argument */ ) = default;

    //! Destructor
    virtual ~Data();

    //! Comparison operator for sorting data-points by distance from the origin
    //! \return Whether this data-point is closer to the origin than another
    //! \param aOther A data-point to compare against
    inline bool operator< ( const Data& aOther ) const
    {
      return r < aOther.r;
    }

    //! Return the squared-distance of this data-points from another
    //! \return The squared-distance of this data-points from another
    //! \param aOther A data-point to compare against
    inline PRECISION dR2( const Data& aOther ) const
    {
      PRECISION dX( x - aOther.x ), dY( y - aOther.y );
      return ( dX*dX ) + ( dY*dY );
    }

    //! Return the distance of this data-points from another
    //! \return The distance of this data-points from another
    //! \param aOther A data-point to compare against
    inline PRECISION dR( const Data& aOther ) const
    {
      return sqrt( dR2( aOther ) );
    }

    //! Return the angle between this data-points and another
    //! \return The angle between this data-points and another
    //! \param aOther A data-point to compare against
    inline PRECISION dPhi( const Data& aOther ) const
    {
      return fabs( phi - aOther.phi );
    }

    //! All the necessary pre-processing to get this data-point ready for an RT-scan
    //! \param aData   The collection of data-points
    //! \param aIndex  The index of the current data-point
    //! \param aMax2R  Twice the maximum radius out to which we will cluster
    //! \param aMax2R2 Square of twice the maximum radius out to which we will cluster
    //! \param aSigmabins2 The sigma-bins for initializing clusters
    //! \param aProgressBar     The progress bar to update
    void Preprocess( std::vector<Data>& aData, const std::size_t& aIndex, const double& aMax2R, const double& aMax2R2, const std::vector< double >& aSigmabins2, ProgressBar& aProgressBar );

    //! Calculate the localization score from the local neighbourhood
    //! \todo Remind myself how this works and what the difference is with below
    //! \param aData        ?
    //! \param aScanConfig  The configuration parameters for the scan
    //! \param aArea        The area of the window for normalizing the log score
    void PreprocessLocalizationScores( std::vector<Data>& aData, const ScanConfiguration& aScanConfig, const double& aArea );

    //! Calculate the localization score from the local neighbourhood
    //! \todo Remind myself how this works and what the difference is with above
    //! \param aData ?
    //! \param R ?
    //! \param aArea        The area of the window for normalizing the log score
    //! \return The localization score
    PRECISION CalculateLocalizationScore( const std::vector<Data>& aData, const double& R, const double& aArea ) const;

  public:
    //! The x-position of the data-point
    PRECISION x;
    //! The y-position of the data-point
    PRECISION y;
    //! The sigma of the data-point
    PRECISION s;
    //! The squared radial distance of the data-point
    PRECISION r2;
    //! The radial distance of the data-point
    PRECISION r;
    //! The phi-position of the data-point
    PRECISION phi;
    //! The locaalization scores, one per R-bin
    std::vector< PRECISION > mLocalizationScores;
    //! The list of neighbours as a pair of squared-distance and index into the list of points
    std::vector< std::pair< PRECISION, std::size_t > > mNeighbours;
    //! A cluster containing only this data-point
    Cluster* mProtoCluster;
};
// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

