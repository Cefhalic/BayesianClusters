//! \file DataProxy.cpp

/* ===== Cluster sources ===== */
#include "BayesianClustering/DataProxy.hpp"
#include "BayesianClustering/Data.hpp"
#include "BayesianClustering/Cluster.hpp"
#include "BayesianClustering/RoIproxy.hpp"
#include "BayesianClustering/RoI.hpp"

/* ===== C++ ===== */
#include <iostream>

//! The maximum depth for recursive clustering
#define RECURSION_LIMIT 75000

// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
DataProxy::DataProxy( Data& aData ) :
mData( &aData ),
mCluster( NULL ),
mExclude( true )
{}

void DataProxy::Clusterize( const PRECISION& a2R2 , RoIproxy& aRoI ) // We are at the top-level
{
  if( mCluster || mExclude ) return;
  aRoI.mClusters.emplace_back( mData->mProtoCluster->mParams.size() );
  try
  {
    Clusterize( a2R2 , aRoI , &aRoI.mClusters.back() );
  }
  catch ( const std::overflow_error& ex )
  { 
    std::cout << "Recursion limit" << std::endl; 
  }
}

void DataProxy::Clusterize( const PRECISION& a2R2 , RoIproxy& aRoI , Cluster* aCluster , const std::size_t& d )
{
  if( mCluster )
  {
    if( GetCluster() == aCluster ) return;
    *aCluster += *mCluster;
    mCluster->mParent = aCluster;
    mCluster->mClusterSize = 0;
    mCluster = aCluster;
  }
  else
  {
    if( mExclude ) return;

    *aCluster += *(mData->mProtoCluster);
    mCluster = aCluster;

    if( d > RECURSION_LIMIT ) throw std::overflow_error( "Recursion limit" );

    for( auto& i : mData->mNeighbours )
    {
      if( i.first > a2R2 ) break;
      aRoI.GetData( i.second ).Clusterize( a2R2 , aRoI , aCluster , d+1 );
    }  
  }
}
// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
