
/* ===== For Root ===== */
#include "Math/SpecFunc.h" 

/* ===== Cluster sources ===== */
#include "BayesianClustering/EventProxy.hpp"
#include "BayesianClustering/Event.hpp"

// /* ===== C++ ===== */
#include <iostream>

// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
EventProxy::EventProxy( Event& aEvent ) :
  mBackgroundCount( 0 )
{
  mClusters.reserve( aEvent.mData.size() );  // Reserve as much space for clusters as there are data points - prevent pointers being invalidated!
  mData.reserve( aEvent.mData.size() );
  for( auto& i : aEvent.mData ) mData.emplace_back( i );
}

void EventProxy::CheckClusterization( const double& R , const double& T )
{
  const auto lRlimit = 4.0 * R * R;

  uint32_t lClusterCount( 0 );

  for( auto& i : mClusters )
  {
    if( i.mClusterSize ) ++lClusterCount;
  }

  if( lClusterCount != mClusterCount )
  {
    std::cout << "\nR = " << R << ", T = " << T << " | #Clusters = " << mClusterCount << " | Expected #Clusters = " << lClusterCount << std::endl;
    throw std::runtime_error( "Check failed" );
  }  

  uint32_t lBackgroundCount( 0 );
  uint32_t lPointsInClusters( 0 );

  uint32_t lExpected( 0 );
  uint32_t lNotClustered( 0 );
  uint32_t lNeighbourNotClustered( 0 );
  uint32_t lWrongNeighbour( 0 );

  for( auto& i : mData )
  {
    if( i.mExclude ){ 
      lBackgroundCount++;
      continue;
    }
    
    lExpected++;
    if( ! i.GetCluster() ){ lNotClustered++ ; continue; }
    for( auto& j : i.mData->mNeighbours )
    {
      if( j.first > lRlimit ) break;
      auto& lNeighbour( GetData( j.second ) );  

      if( lNeighbour.mExclude ) continue;

      if( ! i.mCluster ){ lNeighbourNotClustered++; continue; }
      if ( lNeighbour.mCluster->GetParent() != i.mCluster )
      { 
        lWrongNeighbour++;
        continue; 
      }
    }    
  }

  for( auto& i : mClusters ) lPointsInClusters += i.mClusterSize;

  if( lBackgroundCount != mBackgroundCount )
  {
    std::cout << "\nR = " << R << ", T = " << T << " | Background = " << mBackgroundCount  << " | Expected Background = " << lBackgroundCount << std::endl;
    throw std::runtime_error( "Check failed" ); 
  }  

  if( lPointsInClusters + lBackgroundCount != mData.size() )
  {
    std::cout << "\nR = " << R << ", T = " << T << " | Points In Clusters = " << lPointsInClusters  << " | Background = " << lBackgroundCount << " | Total = " << mData.size() << std::endl;
    throw std::runtime_error( "Check failed" ); 
  }  

  if( lNotClustered or lNeighbourNotClustered or lWrongNeighbour )
  {
    std::cout << "\nR = " << R << ", T = " << T << " | Not Clustered = " << lNotClustered << "/" << lExpected << " | Neighbour Not Clustered = " << lNeighbourNotClustered << " | Wrong Neighbour = " << lWrongNeighbour << std::endl;
    throw std::runtime_error( "Check failed" );
  }
}

__attribute__((flatten))
void EventProxy::ScanRT( const std::function< void( const EventProxy& , const double& , const double& ) >& aCallback , const uint8_t& aParallelization , const uint8_t& aOffset )
{
  double dR( aParallelization * Event::mParameters.dR() );
  double R( Event::mParameters.minScanR() + ( aOffset * Event::mParameters.dR() ) ) , R2( 0 ) , twoR2( 0 ) , T( 0 );

  for( uint32_t i( aOffset ) ; i<Event::mParameters.Rbins() ; i+=aParallelization , R+=dR )
  {
    R2 = R * R;
    twoR2 = 4.0 * R2;
    T = Event::mParameters.maxScanT();

    mClusters.clear();
    for( auto& k : mData ) k.mCluster = NULL;

    for( uint32_t j(0) ; j!=Event::mParameters.Tbins() ; ++j , T-=Event::mParameters.dT() )
    {
      for( auto& k : mData ) k.mExclude = ( k.mData->mLocalizationScores[ i ] < T ) ;
      for( auto& k : mData ) k.Clusterize( twoR2 , T , *this );
      UpdateLogScore();
      if( Event::mParameters.validate() ) CheckClusterization( R , T ) ;
      aCallback( *this , R , T );
    }
  }

  mClusters.clear();
  for( auto& k : mData ) k.mCluster = NULL; // Clear cluster pointers which will be invalidated when we leave the function
}

void EventProxy::UpdateLogScore()
{
  mClusterCount = mClusteredCount = 0;
  mLogP = 0.0;

  for( auto& i: mClusters )
  {
    if( i.mClusterSize == 0 ) continue;
    i.UpdateLogScore();
    mClusterCount += 1;
    mClusteredCount += i.mClusterSize;
    mLogP += ROOT::Math::lgamma( i.mClusterSize );
  }

  mBackgroundCount = mData.size() - mClusteredCount;
  mLogP += ( mBackgroundCount * Event::mParameters.logPb() ) 
         + ( mClusteredCount * Event::mParameters.logPbDagger() )
         + ( Event::mParameters.logAlpha() * mClusterCount )
         + Event::mParameters.logGammaAlpha()
         - ROOT::Math::lgamma( Event::mParameters.alpha() + mClusteredCount );  
}
// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

