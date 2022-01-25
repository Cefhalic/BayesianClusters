

/* ===== C++ ===== */
#include <stdlib.h>
#include <iostream>
#include <iomanip>

/* ===== For Root ===== */
#include "Math/ProbFunc.h" 
#include "Math/Interpolator.h" 
#include "Math/SpecFunc.h" 

/* ===== Cluster sources ===== */
#include "Cluster_GlobalVars.hpp"
#include "Cluster_Data.hpp"

/* ===== Local utilities ===== */
#include "ListComprehension.hpp"
#include "ProgressBar.hpp"
#include "Vectorize.hpp"







Cluster::Parameter::Parameter() : 
A(0.0) , Bx(0.0) , By(0.0) , C(0.0) , logF(0.0)
{}

    
Cluster::Parameter& Cluster::Parameter::operator+= ( const Cluster::Parameter& aOther )
{
  A += aOther.A;
  Bx += aOther.Bx;
  By += aOther.By;
  C += aOther.C;
  logF += aOther.logF;
  return *this;
}


inline double CDF( const double& aArg )
{
  // Above or below ~8 are indistinguishable from 0 and 1 respectively
  if( aArg > 8.0 ) return 1.0;
  if( aArg < 8.0 ) return 0.0;
  return  ROOT::Math::normal_cdf( aArg );
}

__attribute__((flatten))
double Cluster::Parameter::log_score() const
{
  auto sqrt_A( sqrt( A ) ) , inv_A( 1.0 / A );
  auto Dx( Bx * inv_A ) , Dy( By * inv_A );
  auto E( C - ( Bx * Dx ) - ( By * Dy ) );

  double log_sum = logF - double( log( A ) ) + ( 0.5 * E );

  // We place explicit bounds checks to prevent calls to expensive functions
  auto Gx = CDF( sqrt_A * (1.0-Dx) ) - CDF( sqrt_A * (-1.0-Dx) );
  if( Gx != 1.0 ) log_sum += log( Gx );
  auto Gy = CDF( sqrt_A * (1.0-Dy) ) - CDF( sqrt_A * (-1.0-Dy) );
  if( Gy != 1.0 ) log_sum += log( Gy );

  return log_sum;
}



Cluster::Cluster(): mParams( Event::mParameters.sigmacount() ),
mClusterSize( 0 ) , mLastClusterSize( 0 ) , mClusterScore( 0.0 ) , 
mParent( NULL )
{}


Cluster::Cluster( const Data& aData ): mParams( Event::mParameters.sigmacount() ),
mClusterSize( 1 ) , mLastClusterSize( 0 ) , mClusterScore( 0.0 ) , 
mParent( NULL )
{ 
  const auto s2 = aData.s * aData.s;
  auto lIt( mParams.begin() ) ;
  auto lSig2It( Event::mParameters.sigmabins2().begin() );

  for( ; lIt != mParams.end() ; ++lIt , ++lSig2It )
  {
    double w = 1.0 / ( s2 + *lSig2It );
    lIt->A = w;
    lIt->Bx = (w * aData.x);
    lIt->By = (w * aData.y);
    lIt->C = (w * aData.r2);
    lIt->logF = PRECISION( log( w ) );
  }
}


void Cluster::UpdateLogScore()
{
  static constexpr double pi = atan(1)*4;
  static constexpr double log2pi = log( 2*pi );

  if( mClusterSize <= mLastClusterSize ) return; // We were not bigger than the previous size when we were evaluated - score is still valid
  mLastClusterSize = mClusterSize;

  thread_local static std::vector< double > MuIntegral( Event::mParameters.sigmacount() , 1.0 );

  // double constant( mParams[0].log_score() + Event::mParameters.log_probability_sigma( 0 ) );
  // for( std::size_t i(1) ; i!=Event::mParameters.sigmacount() ; ++i ) MuIntegral[i] = exp( mParams[i].log_score() + Event::mParameters.log_probability_sigma( i ) - constant );
  for( std::size_t i(0) ; i!=Event::mParameters.sigmacount() ; ++i ) MuIntegral[i] = exp( mParams[i].log_score() + Event::mParameters.log_probability_sigma( i ) );

  thread_local static ROOT::Math::Interpolator lInt( Event::mParameters.sigmacount() , ROOT::Math::Interpolation::kLINEAR );
  lInt.SetData( Event::mParameters.sigmabins() , MuIntegral );

  static const double Lower( Event::mParameters.sigmabins(0) ) , Upper( Event::mParameters.sigmabins(Event::mParameters.sigmacount()-1) );
  // mClusterScore = double( log( lInt.Integ( Lower , Upper ) ) ) + constant - double( log( 4.0 ) ) + (log2pi * (1.0-mClusterSize));  
  mClusterScore = double( log( lInt.Integ( Lower , Upper ) ) ) - double( log( 4.0 ) ) + (log2pi * (1.0-mClusterSize));  
}


Cluster& Cluster::operator+= ( const Cluster& aOther )
{
  if( &aOther == this ) throw std::runtime_error( "Error #1" );
  if( aOther.mClusterSize == 0 ) throw std::runtime_error( "Error #2" );
  if( aOther.mParent ) throw std::runtime_error( "Error #3" );

  auto lIt( mParams.begin() );
  auto lIt2( aOther.mParams.begin() );

  for(  ; lIt != mParams.end() ; ++lIt , ++lIt2 ) *lIt += *lIt2;
  mClusterSize += aOther.mClusterSize;
  return *this;
}


Cluster* Cluster::GetParent()
{
  if( mParent ) return mParent = mParent->GetParent();
  return this;
}



Data::Data( const PRECISION& aX , const PRECISION& aY , const PRECISION& aS ) : 
x(aX) , y(aY) , s(aS) , r2( (aX*aX) + (aY*aY) ), r( sqrt( r2 ) ), phi( atan2( aY , aX ) ),
mCluster{ NULL } , 
mProtoCluster( NULL )
{}





__attribute__((flatten))
void Data::Preprocess( Event& aEvent )
{

  const std::size_t lIndex( this - aEvent.mData.data() );
  std::vector<Data>::iterator aPlusIt(  aEvent.mData.begin() + lIndex + 1 ) , aPlusEnd( aEvent.mData.end() );
  std::vector<Data>::reverse_iterator aMinusIt( aEvent.mData.rbegin() + aEvent.mData.size() - lIndex ) , aMinusEnd( aEvent.mData.rend() );

  static constexpr double pi = atan(1)*4;  

  // -----------------------------------------------------------------------------------------

  auto dphi = aEvent.mParameters.max2R() / ( r - aEvent.mParameters.max2R() );
  auto dphi2 = (2*pi) - dphi;

  // Iterate over other hits and populate the mNeighbour list
  for( ; aPlusIt != aPlusEnd ; aPlusIt++ )
  {
    if( ( aPlusIt->r - r ) > aEvent.mParameters.max2R() ) break; // aPlusIt is always further out than curent 
    auto lPhi = dPhi( *aPlusIt );
    if( lPhi > dphi and lPhi < dphi2 ) continue;
    PRECISION ldR2 = dR2( *aPlusIt );
    if( ldR2 < aEvent.mParameters.max2R2() ) mNeighbours.push_back( std::make_pair( ldR2 , &*aPlusIt ) );
  }

  for( ; aMinusIt != aMinusEnd ; aMinusIt++ )
  {
    if( ( r - aMinusIt->r ) > aEvent.mParameters.max2R() ) break; // curent is always further out than aMinusIn
    auto lPhi = dPhi( *aMinusIt );
    if( lPhi > dphi and lPhi < dphi2 ) continue;
    PRECISION ldR2 = dR2( *aMinusIt );    
    if( ldR2 < aEvent.mParameters.max2R2() ) mNeighbours.push_back( std::make_pair( ldR2 , &*aMinusIt ) );
  }

  std::sort( mNeighbours.begin() , mNeighbours.end() );

  // -----------------------------------------------------------------------------------------

  const double lLocalizationConstant( 4.0 / ( pi * ( aEvent.mData.size() - 1 ) ) ); 
  const PRECISION eX( 1 - fabs( x ) ) , eY( 1 - fabs( y ) );

  auto lNeighbourit( mNeighbours.begin() );
  PRECISION lLocalizationSum( 0 ) , lLastLocalizationSum( 0 ) , lLocalizationScore( 0 ) , lDist( 0 ) , lWeight( 0 );
  mLocalizationScores.reserve( aEvent.mParameters.Rbins() );

  bool lCreateProtocluster( false );

  double R( aEvent.mParameters.minScanR() ) , R2( 0 );
  for( uint32_t i(0) ; i!=aEvent.mParameters.Rbins() ; ++i , R+=aEvent.mParameters.dR() )
  {
    R2 = R * R;

    for(  ; lNeighbourit != mNeighbours.end() ; ++lNeighbourit )
    { 
      if( lNeighbourit->first > R2 ) break;
      lDist = sqrt( lNeighbourit->first );

      // Noticeably faster polynomial approximation of the edge-correction
      lWeight = 1.0;
      if( eX < lDist )  lWeight *= ( 1 + pow( acos( eX/lDist ) * (2/pi) , 4 ) );
      if( eY < lDist )  lWeight *= ( 1 + pow( acos( eY/lDist ) * (2/pi) , 4 ) );
      lLocalizationSum += lWeight;
    }

    if( lLastLocalizationSum != lLocalizationSum )
    {
      lLocalizationScore = sqrt( lLocalizationConstant * lLocalizationSum );
      if ( lLocalizationScore >= aEvent.mParameters.minScanT() ) lCreateProtocluster = true; // The scan will cluster us at some point
      lLastLocalizationSum = lLocalizationSum;
    }

    mLocalizationScores.push_back( lLocalizationScore );
  }

  if( lCreateProtocluster ) mProtoCluster = new Cluster( *this );
  // -----------------------------------------------------------------------------------------
}






// We are at the top-level
__attribute__((flatten))
void Data::Clusterize( const PRECISION& a2R2 , Event::Instance& aInstance )
{
  if( mCluster[ aInstance.mIndex ] ) return;
  if( mExclude[ aInstance.mIndex ] ) return;

  // if one of our mNeighbours is already a cluster, join that
  for( auto& j : mNeighbours )
  {
    if( j.first > a2R2 ) break;
    if( ! j.second->GetCluster( aInstance.mIndex ) ) continue;  
    return Clusterize( a2R2 , aInstance , j.second->mCluster[ aInstance.mIndex ] );
  }

  // else create a new cluster
  aInstance.mClusters.emplace_back();
  Clusterize( a2R2 , aInstance , &aInstance.mClusters.back() );
}


void Data::Clusterize( const PRECISION& a2R2 , Event::Instance& aInstance , Cluster* aCluster )
{
  if( GetCluster( aInstance.mIndex ) )
  {
    if( mCluster[ aInstance.mIndex ] == aCluster ) return;
    *aCluster += *mCluster[ aInstance.mIndex ];
    mCluster[ aInstance.mIndex ]->mParent = aCluster;
    mCluster[ aInstance.mIndex ]->mClusterSize = 0;
    mCluster[ aInstance.mIndex ] = aCluster;
  }
  else
  {
    if( mExclude[ aInstance.mIndex ] ) return;

    if( !mProtoCluster ) throw std::runtime_error( "NULL ProtoCluster!" );
    *aCluster += *mProtoCluster;
    mCluster[ aInstance.mIndex ] = aCluster;

    for( auto& i : mNeighbours )
    {
      if( i.first > a2R2 ) break;
      i.second->Clusterize( a2R2 , aInstance , aCluster );
    }  
  }
}







GlobalVars Event::mParameters;

Event::Event( const double& aPhysicalCentreX , const double& aPhysicalCentreY ) :
  mPhysicalCentreX( aPhysicalCentreX ) , mPhysicalCentreY( aPhysicalCentreY ) , 
  mInstances( PARALLELIZATION )
{}


void Event::Instance::CheckClusterization( std::vector< Data >& aData , const double& R , const double& T )
{
  if( ! mParameters.validate() ) return;

  const auto lRlimit = 4.0 * R * R;

  uint32_t lClusterCount( 0 );
  uint32_t lBackgroundCount( 0 );
  uint32_t lPointsInClusters( 0 );

  uint32_t lExpected( 0 );
  uint32_t lNotClustered( 0 );
  uint32_t lNeighbourNotClustered( 0 );
  uint32_t lWrongNeighbour( 0 );

  for( auto& i : mClusters )
  {
    if( i.mClusterSize ) ++lClusterCount;
    lPointsInClusters += i.mClusterSize;
  }



  for( auto& i : aData )
  {
    if( i.mExclude[mIndex] ){ 
      lBackgroundCount++;
      continue;
    }
    
    lExpected++;
    if( ! i.GetCluster( mIndex ) ){ lNotClustered++ ; continue; }
    for( auto& j : i.mNeighbours )
    {
      if( j.first > lRlimit ) break;
      if( j.second->mExclude[mIndex] ) continue;

      if( ! j.second->mCluster[mIndex] ){ lNeighbourNotClustered++; continue; }
      if ( j.second->mCluster[mIndex]->GetParent() != i.mCluster[mIndex] )
      { 
        lWrongNeighbour++;
        continue; 
      }
    }    
  }

  //std::cout << "R = " << R << ", T = " << T << " | Clusters = " << lClusterCount << " | Clustered Points = " << lPointsInClusters << " | Background Points = " << lBackgroundCount << std::endl;

  if( lClusterCount != mClusterCount )
  {
    std::cout << "R = " << R << ", T = " << T << " | #Clusters = " << mClusterCount << " | Expected #Clusters = " << lClusterCount << std::endl;
    throw std::runtime_error( "Check failed" );
  }  

  if( lBackgroundCount != mBackgroundCount )
  {
    std::cout << "R = " << R << ", T = " << T << " | Background = " << mBackgroundCount  << " | Expected Background = " << lBackgroundCount << std::endl;
    throw std::runtime_error( "Check failed" ); 
  }  

  if( lPointsInClusters + lBackgroundCount != aData.size() )
  {
    std::cout << "R = " << R << ", T = " << T << " | Points In Clusters = " << lPointsInClusters  << " | Background = " << lBackgroundCount << " | Total = " << aData.size() << std::endl;
    throw std::runtime_error( "Check failed" ); 
  }  

  if( lNotClustered or lNeighbourNotClustered or lWrongNeighbour )
  {
    std::cout << "R = " << R << ", T = " << T << " | Not Clustered = " << lNotClustered << "/" << lExpected << " | Neighbour Not Clustered = " << lNeighbourNotClustered << " | Wrong Neighbour = " << lWrongNeighbour << std::endl;
    throw std::runtime_error( "Check failed" );
  }

}




// __attribute__((flatten))
// void Event::Clusterize( const double& R , const double& T )
// {
//   ProgressBar2 lProgressBar( "Clustering" , 0 );

//   const std::size_t N( mData.size()-1 );
//   // Reset ahead of UpdateLocalization and Clusterize
//   []( Data& i ){ i.mCluster = NULL; i.mLocalizationSum = i.mLocalizationScore = 0.0; i.mNeighbourit = i.mNeighbours.begin(); } || mData;

//   // Update the localization score
//   [ &R , &N ]( Data& i ){ i.UpdateLocalization( R * R , N ); } || mData; // Use interleaving threading to average over systematic radial scaling

//   // And clusterize
//   auto twoR2 = 4.0*R*R;
//   for( auto& i : mData ) i.Clusterize( twoR2 , T , mClusters );
// }


__attribute__((flatten))
void Event::ScanRT( const std::function< void( const Event::Instance& , const double& , const double& ) >& aCallback )
{
  const std::size_t Nminus1( mData.size() - 1 );

  auto DR = PARALLELIZATION * mParameters.dR();

  auto lGenerateLambda = [&]( const std::size_t& offset )
  {
    return [ & , offset ](){
      // std::cout << offset << std::endl;
      auto& lInstance = mInstances.at( offset );
      double R( mParameters.minScanR() + ( offset * mParameters.dR() ) ) , twoR2( 0 ) , T( 0 );

      for( uint32_t i( offset ) ; i<mParameters.Rbins() ; i+=PARALLELIZATION , R+=DR )
      {
        lInstance.mClusters.clear();
        for( auto& k : mData ) k.mCluster[offset] = NULL;

        twoR2 = 4.0 * R * R;
        T = mParameters.maxScanT();

        for( uint32_t j(0) ; j!=mParameters.Tbins() ; ++j , T-=mParameters.dT() )
        {
          for( auto& k : mData ) k.mExclude[offset] = ( k.mLocalizationScores[ i ] < T ) ;
          for( auto& k : mData ) k.Clusterize( twoR2 , lInstance );
          lInstance.UpdateLogScore( mData );
          lInstance.CheckClusterization( mData , R , T ) ;
          //aCallback( lInstance , R , T );
        }
      }

      lInstance.mClusters.clear();
      for( auto& k : mData ) k.mCluster[offset] = NULL;
    };
  };

  auto Thread( ThreadPool.begin() );
  ProgressBar2 lProgressBar( "Scan over RT" , mParameters.Rbins() * mParameters.Tbins() );
  for( std::size_t offset(0) ; offset!=(PARALLELIZATION-1) ; ++offset , ++Thread ) (**Thread).submit( lGenerateLambda( offset ) );
  WrappedThread::run_and_wait( lGenerateLambda( PARALLELIZATION-1 ) );

}



void Event::Preprocess()
{
  {
    ProgressBar2 lProgressBar( "Populating neighbourhood and calculating localization" , mData.size() );
    [&]( Data& k ){ k.Preprocess( *this ); } || mData ;  // Interleave threading since processing time increases with radius from origin
  }

  for( std::size_t i(0) ; i != Concurrency+1 ; ++i )
  {
    auto& lInstance = mInstances[i];    
    lInstance.mClusters.reserve( mData.size() );  // Reserve as much space for clusters as there are data points - prevent pointers being invalidated!
    lInstance.mIndex = i;
  }

}



void Event::Instance::UpdateLogScore( const std::vector< Data >& aData )
{
  mClusterCount = mClusteredCount = 0;
  mLogP = 0.0;

  // Independent - could be parallelized?
  for( auto& i: mClusters )
  {
    if( i.mClusterSize == 0 ) continue;
    i.UpdateLogScore();
    mClusterCount += 1;
    mClusteredCount += i.mClusterSize;
    mLogP += ROOT::Math::lgamma( i.mClusterSize );
  }

  mBackgroundCount = aData.size() - mClusteredCount;
  mLogP += ( mBackgroundCount * mParameters.logPb() ) 
         + ( mClusteredCount * mParameters.logPbDagger() )
         + ( mParameters.logAlpha() * mClusterCount )
         + mParameters.logGammaAlpha()
         - ROOT::Math::lgamma( mParameters.alpha() + mClusteredCount );  
}