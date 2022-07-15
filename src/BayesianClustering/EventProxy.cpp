
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
void EventProxy::ScanRT( const std::function< void( const EventProxy& , const double& , const double& ) >& aCallback , const uint8_t& aParallelization , const uint8_t& aOffset )
{
  double dR( aParallelization * Configuration::Instance.dR() );
  double R( Configuration::Instance.minScanR() + ( aOffset * Configuration::Instance.dR() ) ) , twoR2( 0 ) , T( 0 );

  for( uint32_t i( aOffset ) ; i<Configuration::Instance.Rbins() ; i+=aParallelization , R+=dR )
  {
    twoR2 = 4.0 * R * R;
    T = Configuration::Instance.maxScanT();

    mClusters.clear();
    for( auto& k : mData ) k.mCluster = NULL;

    for( uint32_t j(0) ; j!=Configuration::Instance.Tbins() ; ++j , T-=Configuration::Instance.dT() )
    {
      for( auto& k : mData ) k.mExclude = ( k.mData->mLocalizationScores[ i ] < T ) ;
      for( auto& k : mData ) k.Clusterize( twoR2 , *this );
      UpdateLogScore();
      if( Configuration::Instance.validate() ) CheckClusterization( R , T ) ;
      aCallback( *this , R , T );
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



// void EventProxy::LogScore(){
// //we operate on the mData objects
// // for each data point, get the parent cluster
// // update the parent cluster - we aim to construct the mean central point
//     Cluster* parent;
//     Data* datapoint;
    
//   for (auto& i : mData){
//     //i is of type data proxy
//     parent = i.GetCluster();
//     datapoint = i.mData; //will this be confusing?

//     //we want to get the x and y position of each and then operate on them
//     parent->nuBarX /= parent->nTilde //I want to divide the stored value of nuBar by nTilde
//     parent->nuBarY /= parent->nTilde 

//     //next build S2
//     //we need w_i for each point 
//   }
// }

// inline double CDF( const double& aArg )
// {
//   // Above or below ~8 are indistinguishable from 0 and 1 respectively
//   if( aArg > 8.0 ) return 1.0;
//   if( aArg < -8.0 ) return 0.0;
//   return  ROOT::Math::normal_cdf( aArg );
// }


void EventProxy::UpdateLogScore()
{
  //  iterate over all of the clusters, calculate the mean values
    
  for( auto& i: mClusters ){
    //for each of the sigma paras in the cluster
    for (auto& j : i.mParams){
      j.nuBarX /= j.nTilde;
      j.nuBarY /= j.nTilde;
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

      weightedCentreX = lIt -> nuBarX - x;
      weightedCentreY = lIt -> nuBarY - y;
      weightedCentre = weightedCentreX*weightedCentreX + weightedCentreY*weightedCentreY;
      lIt->S2 += w*weightedCentre;
    }
  }

  // from here, call Cluster::UpdateLogScore

  //we should be ready to calculate p(nu, sigma), mLogScore 

  mClusterCount = mClusteredCount = 0;
  double lLogPl = 0.0; //this is what i label p(nu, sigma)
  mLogP = 0.0;
  // std::vector< double > lMuIntegral( Configuration::Instance.sigmacount()); //should we put this in the cluster class?
  // double cdfArgumentX, cdfArgumentY;
  double lNTilde, lSqrtNTilde, lMuIntegral;
  static constexpr double pi = atan(1)*4;

  // std::size_t lClusterSize;

  for( auto& i: mClusters ) // here we operate on each of the identified clusters
  // call i.updatelogscore here and put all this code into it
  {
    if( i.mClusterSize == 0 ) continue;
    
    // for (auto& j : i.mParams){
    i.UpdateLogScore();
      
    mLogP += i.mClusterScore;
    }
  

  mBackgroundCount = mData.size() - mClusteredCount;
  lLogPl += ( mBackgroundCount * Configuration::Instance.logPb() ) 
         + ( mClusteredCount * Configuration::Instance.logPbDagger() )
         + ( Configuration::Instance.logAlpha() * mClusterCount )
         + Configuration::Instance.logGammaAlpha()
         - ROOT::Math::lgamma( Configuration::Instance.alpha() + mClusteredCount );  

   // To be implemented...

  mLogP += -log(4.0) * mBackgroundCount + lLogPl; //taking background density to be just the area of the ROI

}
// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

