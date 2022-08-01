
/* ===== For Root ===== */
#include "Math/SpecFunc.h" 

/* ===== Cluster sources ===== */
#include "BayesianClustering/EventProxy.hpp"
#include "BayesianClustering/Event.hpp"
#include "BayesianClustering/Configuration.hpp"

// /* ===== C++ ===== */
#include <iostream>

// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
EventProxy::EventProxy( Event& aEvent ) :
  mBackgroundCount( 0 ) , mEvent( aEvent )
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
void EventProxy::ScanRT( const std::function< void( const EventProxy& , const double& , const double& , std::vector<uint32_t>  ) >& aCallback , const uint8_t& aParallelization , const uint8_t& aOffset )
{
  double dR( aParallelization * Configuration::Instance.dR() );
  double R( Configuration::Instance.minScanR() + ( aOffset * Configuration::Instance.dR() ) ) , twoR2( 0 ) , T( 0 );

  for( uint32_t i( aOffset ) ; i<Configuration::Instance.Rbins() ; i+=aParallelization , R+=dR )
  {
    twoR2 = 4.0 * R * R;
    T = Configuration::Instance.maxScanT();

    mClusters.clear();
    for( auto& k : mData ) k.mCluster = NULL;

    std::vector<uint32_t> lCurrentIJ(2, 0);

    for( uint32_t j(0) ; j!=Configuration::Instance.Tbins() ; ++j , T-=Configuration::Instance.dT() )
    {
      for( auto& k : mData ) k.mExclude = ( k.mData->mLocalizationScores[ i ] < T ) ;
      for( auto& k : mData ) k.Clusterize( twoR2 , *this );
      UpdateLogScore();
      if( Configuration::Instance.validate() ){
        CheckClusterization( R , T ) ;
        ValidateLogScore();
        }

      //place to store current ij
      lCurrentIJ[0] = i;
      lCurrentIJ[1] = j;
 
      aCallback( *this , R , T, lCurrentIJ );
    }
  }

  mClusters.clear();
  for( auto& k : mData ) k.mCluster = NULL; // Clear cluster pointers which will be invalidated when we leave the function
}


void EventProxy::Clusterize( const double& R , const double& T , const std::function< void( const EventProxy& ) >& aCallback )
{
  auto twoR2 = 4.0 * R * R;

  mClusters.clear();
  for( auto& k : mData )
  { 
    k.mCluster = NULL;
    k.mExclude = ( k.mData->CalculateLocalizationScore( mEvent.mData , R ) < T ) ;
  }

  for( auto& k : mData ) k.Clusterize( twoR2 , *this );

  UpdateLogScore();
  aCallback( *this );

}


void EventProxy::ValidateLogScore()
{
  for ( auto& i : mClusters)
  {
    if (i.mClusterSize == 0 ) continue;
    for (auto& j : i.mParams)
    {
      j.weightedCentreX = j.Bx / j.A;
      j.weightedCentreY = j.By / j.A;
    }
  }
  //iterate over dPoints here, update cluster S2
  Cluster* parent;
  Data* datapoint;
  double x, y;

  for (auto& i : mData){
    parent = i.GetCluster();
    datapoint = i.mData; //i is of type DataProxy, with member variable mData, which is of type Data*

    if (!parent) continue; //continue if no parent

    //get the coord centres
    x = datapoint->x;
    y = datapoint->y;
    auto s = datapoint->s;
    auto s2 = s * s; //bad naming! please redo
    double weightedCentre, weightedCentreX, weightedCentreY; 

    //update S2 for each sigma hypothesis
    //we need to recalculate w here i think 
    
    auto lIt(parent->mParams.begin());
    auto lSig2It( Configuration::Instance.sigmabins2().begin() );
    for ( ; lIt != parent->mParams.end() ; ++lIt, ++lSig2It){
      //we need to add on w_i here - which comes with each point in the cluster
      double w = 1.0 / (s2 + *lSig2It); //these are found in the protoclusters, inside datapoint
      // weightedCentreX = lIt -> nuBarX - x;
      // weightedCentreY = lIt -> nuBarY - y;
      weightedCentreX = (lIt -> weightedCentreX) - x;
      weightedCentreY = (lIt -> weightedCentreY) - y;
      weightedCentre = weightedCentreX*weightedCentreX + weightedCentreY*weightedCentreY;
      lIt->S2 += w*weightedCentre; 
  }
}
  //NEXT - we perform an alternate log_score 
  //and compare it with the usual log_score

  double fastLogScore, valLogScore;
  for (auto& i : mClusters)
  {
    if (i.mClusterSize == 0) continue;
    for( std::size_t j(0) ; j!=Configuration::Instance.sigmacount() ; ++j )
    {
      fastLogScore = i.mParams[j].log_score();
      valLogScore = i.mParams[j].alt_log_score();
      std::cout << "here" << std::endl;
      if (abs(fastLogScore - valLogScore) > 0.001) throw std::runtime_error("logscore check failed");
    }
  }

}


void EventProxy::UpdateLogScore()
{
  mClusterCount = mClusteredCount = 0;
  double lLogPl = 0.0;
  mLogP = 0.0;
  for( auto& i: mClusters ) // here we operate on each of the identified clusters
  {
    if( i.mClusterSize == 0 ) continue;
    
    i.UpdateLogScore();
    mClusterCount += 1;
    mClusteredCount += i.mClusterSize;
    mLogP += i.mClusterScore;
    lLogPl += ROOT::Math::lgamma( i.mClusterSize ); //this was omitted before - why?
    }
  

  mBackgroundCount = mData.size() - mClusteredCount;
  lLogPl += ( mBackgroundCount * Configuration::Instance.logPb() ) 
         + ( mClusteredCount * Configuration::Instance.logPbDagger() )
         + ( Configuration::Instance.logAlpha() * mClusterCount )
         + Configuration::Instance.logGammaAlpha()
         - ROOT::Math::lgamma( Configuration::Instance.alpha() + mClusteredCount );  

  mLogP += (-log(4.0) * mBackgroundCount) + lLogPl;
}
// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

