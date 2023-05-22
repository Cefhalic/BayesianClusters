
/* ===== Cluster sources ===== */
#include "BayesianClustering/DataProxy.hpp"
#include "BayesianClustering/Data.hpp"
#include "BayesianClustering/Cluster.hpp"
#include "BayesianClustering/RoI.hpp"


// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
DataProxy::DataProxy( Data& aData ) :
mData( &aData ),
mCluster( NULL )
{}

void DataProxy::Clusterize( const PRECISION& a2R2 , RoI& aDataset ) // We are at the top-level
{
  if( mCluster || mExclude ) return;

  aDataset.mClusters.emplace_back();
  Clusterize( a2R2 , aDataset , &aDataset.mClusters.back() );
}

void DataProxy::Clusterize( const PRECISION& a2R2 , RoI& aDataset , Cluster* aCluster )
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
      aDataset.GetData( i.second ).Clusterize( a2R2 , aDataset , aCluster );
    }  
  }
}
// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
