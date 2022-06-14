
/* ===== Cluster sources ===== */
#include "BayesianClustering/DataProxy.hpp"
#include "BayesianClustering/Data.hpp"
#include "BayesianClustering/Cluster.hpp"
#include "BayesianClustering/EventProxy.hpp"


// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
DataProxy::DataProxy( Data& aData ) :
mData( &aData ),
mCluster( NULL )
{}

void DataProxy::Clusterize( const PRECISION& a2R2 , const PRECISION& aT , EventProxy& aEvent ) // We are at the top-level
{
  if( mCluster || mExclude ) return;

  aEvent.mClusters.emplace_back();
  Clusterize( a2R2 , aT , aEvent , &aEvent.mClusters.back() );
}

void DataProxy::Clusterize( const PRECISION& a2R2 , const PRECISION& aT , EventProxy& aEvent , Cluster* aCluster )
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

    for( auto& i : mData->mNeighbours )
    {
      if( i.first > a2R2 ) break;
      aEvent.GetData( i.second ).Clusterize( a2R2 , aT , aEvent , aCluster );
    }  
  }
}
// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
