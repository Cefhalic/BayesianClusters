//! \file RoIproxy.hpp
#pragma once

/* ===== C++ ===== */
#include <vector>
#include <functional>
#include <cstdint>

/* ===== Cluster sources ===== */
#include "BayesianClustering/Cluster.hpp"
#include "BayesianClustering/DataProxy.hpp"
#include "BayesianClustering/Configuration.hpp"

class RoI;
class ProgressBar;

// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//! A lightweight wrapper for the RoI to store clusters for a given scan
class RoIproxy
{
  public:
    //! Default constructor
    //! \param aRoI An RoI for which this is a lightweight proxy
    RoIproxy( RoI& aRoI );

    //! Deleted copy constructor
    RoIproxy( const RoIproxy& aOther /*!< Anonymous argument */ ) = delete;

    //! Deleted assignment operator
    //! \return Reference to this, for chaining calls
    RoIproxy& operator = (const RoIproxy& aOther /*!< Anonymous argument */ ) = delete;

    //! Default move constructor
    RoIproxy( RoIproxy&& aOther /*!< Anonymous argument */ ) = default;

    //! Default move-assignment constructor
    //! \return Reference to this, for chaining calls
    RoIproxy& operator = ( RoIproxy&& aOther /*!< Anonymous argument */ ) = default;

    //! Default destructor
    ~RoIproxy();

    //! Run validation tests on the clusters
    //! \param R The R of the last run scan
    //! \param T The T of the last run scan
    void CheckClusterization( const double& R, const double& T );

    //! Run an RT-scan
    //! \param aScanConfig      The configuration parameters for the scan
    //! \param aCallback        A callback for each RT-scan result
    //! \param aProgressBar     The progress bar to update
    //! \param aParallelization The stride with which we will iterate across RT parameters
    //! \param aOffset          The starting point for the strides as we iterate across RT parameters
    //! \param aValidate        Run validation of the score calculation
    void ScanRT( const ScanConfiguration& aScanConfig, const std::function< void( RoIproxy&, const double&, const double& ) >& aCallback, ProgressBar& aProgressBar, const uint8_t& aParallelization = 1, const uint8_t& aOffset = 0, const bool& aValidate = false );

    //! Run clusterization for a specific choice of R and T
    //! \param R The R parameter for clusterization
    //! \param T The T parameter for clusterization
    //! \param aCallback A callback for the clusterization results
    void Clusterize( const double& R, const double& T, const std::function< void( RoIproxy& ) >& aCallback );

    //! Update log-probability after a scan
    //! \param aScanConfig      The configuration parameters for the scan
    void UpdateLogScore( const ScanConfiguration& aScanConfig );

    //! Sean's validation code for testing when the running log-score fails
    //! \param aScanConfig      The configuration parameters for the scan
    void ValidateLogScore( const ScanConfiguration& aScanConfig );

    //! Get the proxy for the Nth neighbour of this data-point
    //! \return A reference to the neighbour data-proxy
    //! \param aIndex The index of the neighbour we are looking for
    inline DataProxy& GetData( const std::size_t& aIndex )
    {
      return mData.at( aIndex );
    }

  public:
    //! The collection of lightweight data-point wrappers used by this RoI wrapper
    std::vector< DataProxy > mData;

    //! The collection of clusters found by this scan
    std::vector< Cluster > mClusters;

    //! The number of clustered data-points
    std::size_t mClusteredCount;

    //! The number of background data-points
    std::size_t mBackgroundCount;

    //! The number of non-Null clusters
    std::size_t mClusterCount;

    //! The log-probability density associated with the last scan
    double mLogP;

    //! The underlying RoI this is a proxy to
    const RoI& mRoI;

};
// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
