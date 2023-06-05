//! \file DataProxy.hpp
#pragma once

/* ===== C++ ===== */
#include <cstddef>

/* ===== Cluster sources ===== */
#include "BayesianClustering/Precision.hpp"
#include "BayesianClustering/Cluster.hpp"

class RoIproxy;
class Cluster;
class Data;


// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//! A light-weight proxy for the raw data-points
class DataProxy
{
public:
  //! Default constructor
  //! \param aData The data-point for which this is the proxy
  DataProxy( Data& aData );

  //! Deleted copy constructor
  DataProxy( const DataProxy& aOther /*!< Anonymous argument */ ) = delete;

  //! Deleted assignment operator
  //! \return Reference to this, for chaining calls  
  DataProxy& operator = (const DataProxy& aOther /*!< Anonymous argument */ ) = delete;

  //! Default move constructor
  DataProxy( DataProxy&& aOther /*!< Anonymous argument */ ) = default;

  //! Default move-assignment constructor
  //! \return Reference to this, for chaining calls  
  DataProxy& operator = ( DataProxy&& aOther /*!< Anonymous argument */ ) = default;

  //! Entry point clusterization function - a new cluster will be created
  //! \param a2R2   The clusterization radius
  //! \param aRoI The RoI-proxy in which we are running
  void Clusterize( const PRECISION& a2R2 , RoIproxy& aRoI );
  
  //! Recursive clusterization function
  //! \param a2R2     The clusterization radius
  //! \param aRoI     The RoI-proxy in which we are running  
  //! \param aCluster The cluster we are building
  //! \param d        The recursion depth
  void Clusterize( const PRECISION& a2R2 , RoIproxy& aRoI , Cluster* aCluster , const std::size_t& d = 0 );
  
  //! Get a pointer to this data-proxy's ultimate parent cluster (or null if unclustered
  //! \return A pointer to this data-proxy's ultimate parent cluster  
  inline Cluster* GetCluster()
  {
    if( ! mCluster ) return NULL;
    return mCluster = mCluster->GetParent();
  }

public:
  //! The data-point for which this is the proxy
  Data* mData;
  
  //! This data-proxy's immediate parent cluster
  Cluster* mCluster;
  
  //! Whether this data-point is to be included in the clusterization
  bool mExclude;
};
// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
