
/* ===== BOOST libraries ===== */
#include <boost/math/special_functions/gamma.hpp>

/* ===== Cluster sources ===== */
#include "BayesianClustering/RoIproxy.hpp"
#include "BayesianClustering/RoI.hpp"

/* ===== Local utilities ===== */
#include "Utilities/ProgressBar.hpp"

// /* ===== C++ ===== */
#include <iostream>

// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
RoIproxy::RoIproxy( RoI& aRoI ) :
  mBackgroundCount( 0 ) , mRoI( aRoI )
{
  mClusters.reserve( aRoI.mData.size() );  // Reserve as much space for clusters as there are data points - prRoI pointers being invalidated!
  mData.reserve( aRoI.mData.size() );
  for( auto& i : aRoI.mData ) mData.emplace_back( i );
}

void RoIproxy::CheckClusterization( const double& R , const double& T )
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
void RoIproxy::ScanRT( const Configuration::tBounds& R , const Configuration::tBounds& T , const std::function< void( const RoIproxy& , const double& , const double& , std::pair<int,int>  ) >& aCallback , const uint8_t& aParallelization , const uint8_t& aOffset )
{
  double dR( aParallelization * R.spacing );
  double lR( R.min + ( aOffset * R.spacing ) ) , twoR2( 0 ) , lT( 0 );

  for( uint32_t i( aOffset ) ; i<R.bins ; i+=aParallelization , lR+=dR )
  {
    twoR2 = 4.0 * lR * lR;
    lT = T.max;

    mClusters.clear();
    for( auto& k : mData ) k.mCluster = NULL;

    std::pair<int,int> lCurrentIJ;

    for( uint32_t j(0) ; j!=T.bins ; ++j , lT-=T.spacing )
    {
      for( auto& k : mData ) k.mExclude = ( k.mData->mLocalizationScores[ i ] < lT ) ;
      for( auto& k : mData ) k.Clusterize( twoR2 , *this );
      UpdateLogScore();
      if( mRoI.mConfiguration.validate() ){
        CheckClusterization( lR , lT ) ;
        ValidateLogScore();
      }

      //place to store current ij
      lCurrentIJ.first = i;
      lCurrentIJ.second = j;
 
      aCallback( *this , lR , lT, lCurrentIJ );
    }
  }

  mClusters.clear();
  for( auto& k : mData ) k.mCluster = NULL; // Clear cluster pointers which will be invalidated when we leave the function
}


void RoIproxy::Clusterize( const double& R , const double& T , const std::function< void( const RoIproxy& ) >& aCallback )
{
  {
    ProgressBar2 lProgressBar( "Clusterize"  , 0 );  
    auto twoR2 = 4.0 * R * R;

    mClusters.clear();
    for( auto& k : mData )
    { 
      k.mCluster = NULL;
      k.mExclude = ( k.mData->CalculateLocalizationScore( mRoI.mData , R , mRoI.getArea() ) < T ) ;
    }

    for( auto& k : mData ) k.Clusterize( twoR2 , *this );

    UpdateLogScore();
  }

  aCallback( *this );
}


void RoIproxy::ValidateLogScore()
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
    auto lSig2It( mRoI.mConfiguration.sigmabins2().begin() );
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
    for( std::size_t j(0) ; j!=mRoI.mConfiguration.sigmacount() ; ++j )
    {
      fastLogScore = i.mParams[j].log_score();
      valLogScore = i.mParams[j].alt_log_score();
      if (abs(fastLogScore - valLogScore) > 5) throw std::runtime_error("logscore check failed");
    }
  }

}


void RoIproxy::UpdateLogScore()
{
  if( mRoI.mConfiguration.sigmabins().size() == 0 ) return;

  mClusterCount = mClusteredCount = 0;
  double lLogPl = 0.0;
  mLogP = 0.0;
  for( auto& i: mClusters ) // here we operate on each of the identified clusters
  {
    if( i.mClusterSize == 0 ) continue;
    
    i.UpdateLogScore( mRoI.mConfiguration.sigmabins() , mRoI.mConfiguration.log_probability_sigma() );
    mClusterCount += 1;
    mClusteredCount += i.mClusterSize;
    mLogP += i.mClusterScore;
    lLogPl += boost::math::lgamma( i.mClusterSize ); //this was omitted before - why?
  }
  
  mBackgroundCount = mData.size() - mClusteredCount;
  lLogPl += ( mBackgroundCount * mRoI.mConfiguration.logPb() ) 
         + ( mClusteredCount * mRoI.mConfiguration.logPbDagger() )
         + ( mRoI.mConfiguration.logAlpha() * mClusterCount )
         + mRoI.mConfiguration.logGammaAlpha()
         - boost::math::lgamma( mRoI.mConfiguration.alpha() + mClusteredCount );  

  mLogP += (-log(4.0) * mBackgroundCount) + lLogPl;
}
// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

