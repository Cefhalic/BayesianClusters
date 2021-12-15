

/* ===== C++ ===== */
#include <stdlib.h>
#include <iostream>

/* ===== For Root ===== */
#include "Math/ProbFunc.h" 
#include "Math/Interpolator.h" 

/* ===== Cluster sources ===== */
#include "Cluster_GlobalVars.hpp"
#include "Cluster_Data.hpp"

/* ===== Local utilities ===== */
#include "ListComprehension.hpp"
#include "ProgressBar.hpp"
#include "Vectorize.hpp"


// // Burmann approximation of the Gaussian cumulative-distribution function
// PRECISION Phi( const PRECISION& x )
// {
//   auto y = exp( -1.0*x*x );

//   constexpr PRECISION pi = atan(1)*4;
//   constexpr PRECISION rt_pi = sqrt( pi );
//   constexpr PRECISION inv_rt_pi = 1.0 / rt_pi;

//   PRECISION lRet = ( rt_pi/2.0 ) + ( y * 31.0/200.0 ) - ( y*y * 341.0/8000.0 );
//   lRet       *= inv_rt_pi * ( (x > 0) - (x < 0) ) * sqrt( 1-y );

//   return 0.5 + lRet;
// }


std::vector< Data::Cluster > Data::Clusters;
std::mutex Data::ClusterMutex;



Data::ClusterParameter::ClusterParameter() : 
A(0.0) , Bx(0.0) , By(0.0) , C(0.0) , logF(0.0)
{}


// Data::ClusterParameter::ClusterParameter( const PRECISION& aW , const PRECISION& aX , const PRECISION& aY , const PRECISION& aR2 ) : 
// A( aW ) , Bx( aW * aX ) , By( aW * aY ) , C( aW * aR2 ) , logF( log( w ) )
// {}
    
Data::ClusterParameter& Data::ClusterParameter::operator+= ( const Data::ClusterParameter& aOther )
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
double Data::ClusterParameter::log_score() const
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



Data::Cluster::Cluster(): Params( Parameters.sigmacount() ),
ClusterSize( 0 ) , LastClusterSize( 0 ) , ClusterScore( 0.0 ),
mParent( NULL )
{}


double Data::Cluster::log_score()
{
  static constexpr double pi = atan(1)*4;
  static constexpr double log2pi = log( 2*pi );

  if( ClusterSize <= LastClusterSize ) return ClusterScore; // We were not bigger than the previous size when we were evaluated - score is still valid
  LastClusterSize = ClusterSize;

  thread_local static std::vector< double > MuIntegral( Parameters.sigmacount() , 1.0 );

  // double constant( Params[0].log_score() + Parameters.log_probability_sigma( 0 ) );
  // for( std::size_t i(1) ; i!=Parameters.sigmacount() ; ++i ) MuIntegral[i] = exp( Params[i].log_score() + Parameters.log_probability_sigma( i ) - constant );
  for( std::size_t i(0) ; i!=Parameters.sigmacount() ; ++i ) MuIntegral[i] = exp( Params[i].log_score() + Parameters.log_probability_sigma( i ) );

  thread_local static ROOT::Math::Interpolator lInt( Parameters.sigmacount() , ROOT::Math::Interpolation::kLINEAR );
  lInt.SetData( Parameters.sigmabins() , MuIntegral );

  static const double Lower( Parameters.sigmabins(0) ) , Upper( Parameters.sigmabins(Parameters.sigmacount()-1) );
  // return ClusterScore = double( log( lInt.Integ( Lower , Upper ) ) ) + constant - double( log( 4.0 ) ) + (log2pi * (1.0-ClusterSize));  
  return ClusterScore = double( log( lInt.Integ( Lower , Upper ) ) ) - double( log( 4.0 ) ) + (log2pi * (1.0-ClusterSize));  
}


Data::Cluster& Data::Cluster::operator+= ( const Data& aData )
{
  auto lIt2( aData.w_i.begin() );
  for( auto lIt( Params.begin() ) ; lIt != Params.end() ; ++lIt , ++lIt2 )
  {
    auto& w = *lIt2;
    lIt->A += w;
    lIt->Bx += (w * aData.x);
    lIt->By += (w * aData.y);
    lIt->C += (w * aData.r2);
    lIt->logF += PRECISION( log( w ) );
  }
  ClusterSize += 1;

  return *this;  
}

Data::Cluster& Data::Cluster::operator+= ( Data::Cluster& aOther )
{
  Data::Cluster* lThis = GetParent();
  Data::Cluster* lOther = aOther.GetParent();

  if( lThis == lOther ) return *this;

  for( auto lIt( lThis->Params.begin() ) , lIt2( lOther->Params.begin() ) ; lIt != lThis->Params.end() ; ++lIt , ++lIt2 ) *lIt += *lIt2;
  lThis->ClusterSize += lOther->ClusterSize;
  lThis->mBoundaries.insert( lOther->mBoundaries.begin() , lOther->mBoundaries.end() );
 
  lOther->ClusterSize = 0;
  lOther->mParent = lThis;
  lOther->mBoundaries.clear();

  return *this;
}


Data::Cluster* Data::Cluster::GetParent()
{
  if( mParent ) return mParent = mParent->GetParent();
  return this;
}


Data::Data( const PRECISION& aX , const PRECISION& aY , const PRECISION& aS ) : 
x(aX) , y(aY) , r2( (aX*aX) + (aY*aY) ), r( sqrt( r2 ) ), phi( atan2( aY , aX ) ),
w_i( [ &aS ]( const double& sig2 ){ return PRECISION( 1.0 / ( (aS*aS) + sig2 ) ); } | Parameters.sigmabins2() ),
localizationsum( 0.0 ) , localizationscore( 0.0 ),
neighbourit( neighbours.end() ),
mCluster( NULL )
{}





__attribute__((flatten))
void Data::PopulateNeighbours( std::vector<Data>::iterator aPlusIt , const std::vector<Data>::iterator& aPlusEnd , std::vector<Data>::reverse_iterator aMinusIt , const std::vector<Data>::reverse_iterator& aMinusEnd )
{
  const auto dphi = Parameters.maxR() / ( r - Parameters.maxR() );
  const auto dphi2 = Parameters.max2R() / ( r - Parameters.max2R() );

  // Iterate over other hits and populate the neighbour list
  for( ; aPlusIt != aPlusEnd ; aPlusIt++ )
  {
    if( ( aPlusIt->r - r ) > Parameters.max2R() ) break; // aPlusIt is always further out than curent 
    if( fabs( aPlusIt->phi - phi ) > dphi ) continue;
    PRECISION ldR2 = dR2( *aPlusIt );
    if( ldR2 < Parameters.max2R2() ) neighbours.push_back( std::make_pair( ldR2 , &*aPlusIt ) );
  }

  for( ; aMinusIt != aMinusEnd ; aMinusIt++ )
  {
    if( ( r - aMinusIt->r ) > Parameters.max2R() ) break; // curent is always further out than aMinusIn
    if( fabs( aPlusIt->phi - phi ) > dphi ) continue;
    PRECISION ldR2 = dR2( *aMinusIt );    
    if( ldR2 < Parameters.max2R2() ) neighbours.push_back( std::make_pair( ldR2 , &*aMinusIt ) );
  }

  std::sort( neighbours.begin() , neighbours.end() );
  neighbourit = neighbours.begin();
}


void Data::UpdateLocalization( const PRECISION& aR2 , const size_t& Nminus1  )
{
  constexpr PRECISION pi = atan(1)*4;

  const PRECISION last_localizationsum( localizationsum );

  for( ; neighbourit != neighbours.end() ; ++neighbourit )
  { 
    if( neighbourit->first > aR2 ) break;
    PRECISION lDist = sqrt( neighbourit->first );

    // Approximation of the edge-correction
    PRECISION Weight( 1.0 );
    PRECISION eX( 1 - fabs( x ) ) , eY( 1 - fabs( y ) );
    if( eX < lDist )  Weight *= ( 1 + pow( acos( eX/lDist ) * (2/pi) , 4) );
    if( eY < lDist )  Weight *= ( 1 + pow( acos( eY/lDist ) * (2/pi) , 4) );

    localizationsum += Weight;
  }

  if( last_localizationsum == localizationsum ) return;

  const PRECISION LocalizationConstant( 4.0 / ( pi * Nminus1 ) ); 
  localizationscore = sqrt( LocalizationConstant * localizationsum );
}


void Data::ResetClusters()
{
  mCluster = NULL;
}


__attribute__((flatten))
void Data::Clusterize( const PRECISION& a2R2 , const PRECISION& aT , const Data* aLower , const Data* aUpper , Cluster* aCluster )
{
  if( mCluster ) return;

  if( localizationscore < aT ) return;

  // if a cluster has been specified for us
  if( aCluster )
  {
    mCluster = aCluster;
    goto breakpoint;
  }

  // if one of our neighbours is already a cluster
  for( auto& j : neighbours )
  {
    if( j.first > a2R2 ) break;
    if( ! j.second->mCluster ) continue;
    if( j.second < aLower or j.second >= aUpper ) continue; // Save us from having to mutex the clusters
  
    mCluster = j.second->mCluster;
    goto breakpoint;
  }

  // else create a new cluster
  {
    std::unique_lock< std::mutex > x( ClusterMutex );
    Clusters.emplace_back();
    mCluster = &Clusters.back();
  }

breakpoint:
  *mCluster += *this;

  for( auto& j : neighbours )
  {
    if( j.first > a2R2 ) break;
    if( j.second->localizationscore < aT ) continue;

    if( j.second < aLower or j.second >= aUpper ){ // Save us from having to mutex the clusters
      if( j.second->mCluster ) mCluster->mBoundaries.insert( j.second->mCluster );
      continue; 
    }

    Data::Cluster* lCluster = j.second->GetCluster();
    if( lCluster ) *mCluster += *lCluster; // Neighbour is already part of a cluster
    else j.second->Clusterize( a2R2 , aT , aLower , aUpper , mCluster );
    
  }
}


Data::Cluster* Data::GetCluster()
{
  if( mCluster ) return mCluster->GetParent();
  return NULL;
}


// __attribute__((flatten))
// void Data::Finalize( const PRECISION& a2R2 , const PRECISION& aT , Cluster* aCluster )
// {
//   if( !incomplete ) return;

//   for( auto& j : neighbours )
//   {
//     if( j.first > a2R2 ) break;
//     if( j.second->localizationscore < aT ) continue;
//     if( j.second->mCluster != mCluster and j.second->mCluster->ClusterSize ) *mCluster += *(j.second->mCluster);
//   }
//   incomplete = false;  
// }







void Clusterize( std::vector<Data>& aData , const double& twoR2 , const double& T )
{
  // for( auto& i : aData ) i.Clusterize( twoR2 , T , &(*aData.begin()) , &(*aData.end()) );
  // return;

  static const std::size_t lChunksize( ceil( double(aData.size()) / Concurrency ) );
  using tIt = std::vector< Data >::iterator;
  static std::vector< std::pair< tIt , tIt > > lBounds;

  if( lBounds.empty() )
  {
    auto A( aData.begin() ) , B( aData.begin() + lChunksize );
    for( ; B < aData.end() ; A = B , B+=lChunksize ) lBounds.push_back( std::make_pair( A , B ) );
    lBounds.push_back( std::make_pair( A , aData.end() ) );
  }

  // Multithreaded clustering within self-contained chunks
  auto Thread( ThreadPool.begin() );
  for( auto lChunk( lBounds.begin() ) ; lChunk != lBounds.end() ; ++Thread , ++lChunk ) (**Thread).submit( [ &twoR2 , &T , lChunk ](){ for( auto i( lChunk->first ) ; i != lChunk->second ; ++i ) i->Clusterize( twoR2 , T , &*lChunk->first , &*lChunk->second ); } );
  WrappedThread::wait();


  // std::size_t cnt(0);
  // for( auto& i : Data::Clusters )
  // {
  //   std::cout << &i << "\t" << i.ClusterSize ;
  //   for ( auto& j : i.mBoundaries ) std::cout << "\t" << j;
  //   std::cout << std::endl; 
  // }

  // // Handle boundaries between self-contained chunks
  // for( auto lChunk( lBounds.begin() + 1 ) ; lChunk != lBounds.end() ; ++Thread , ++lChunk )
  // {
  //   auto r_limit = lChunk->first->r + twoR2;
  //   for( auto i( lChunk->first ) ; i != lChunk->second ; ++i )
  //   {
  //     if( i->r > r_limit ) break;
  //     i->Clusterize2( twoR2 , T );
  //   }
  // }
}

__attribute__((flatten))
void Cluster( std::vector<Data>& aData , const double& R , const double& T )
{
  ProgressBar2 lProgressBar( "Clustering" , 0 );

  const std::size_t N( aData.size()-1 );
  // Reset ahead of UpdateLocalization and Clusterize
  []( Data& i ){ i.ResetClusters(); i.localizationsum = i.localizationscore = 0.0; i.neighbourit = i.neighbours.begin(); } || aData;

  // Update the localization score
  [ &R , &N ]( Data& i ){ i.UpdateLocalization( R * R , N ); } || aData; // Use interleaving threading to average over systematic radial scaling

  // And clusterize
  Clusterize( aData , 4.0*R*R , T );
}


__attribute__((flatten))
void ScanRT( std::vector<Data>& aData , const std::function< void( const double& , const double& ) >& aCallback )
{
  const std::size_t Nminus1( aData.size() - 1 );
  double R( Parameters.minScanR() ) , R2( 0 ) , twoR2( 0 ) , T( 0 );

  ProgressBar lProgressBar( "Scan over RT" , Parameters.Rbins() * Parameters.Tbins() );

  for( uint32_t i(0) ; i!=Parameters.Rbins() ; ++i , R+=Parameters.dR() )
  {
    R2 = R * R;
    twoR2 = 4.0 * R2;
    T = Parameters.maxScanT();

    [ &Nminus1 , &R2 ]( Data& i ){ i.UpdateLocalization( R2 , Nminus1 ); } || aData; // Use interleaving threading to average over systematic radial scaling

    Data::Clusters.clear();
    []( Data& i ){ i.ResetClusters(); } || aData;

    for( uint32_t j(0) ; j!=Parameters.Tbins() ; ++j , T-=Parameters.dT() , ++lProgressBar )
    {
      Clusterize( aData , twoR2 , T );   
      // []( Data& i ){ i.UpdateClusterScore(); } || aData; // Use interleaving threading to average over systematic radial scaling
      aCallback( R , T );
    }

  }
}



void PrepData( std::vector<Data>& aData )
{

  // Should already be sorted, but...
  // std::sort( aData.begin() , aData.end() );

  {
    ProgressBar2 lProgressBar( "Preprocessing" , aData.size() );

    // Populate neighbour lists  
    [ &aData ]( const std::size_t& i ){ aData.at( i ).PopulateNeighbours( aData.begin() + i + 1 , aData.end() , aData.rbegin() + aData.size() - i , aData.rend() ); } || range( aData.size() );  // Interleave threading since processing time increases with radius from origin
  }

  // Reserve as much space for clusters as there are data points - prevent pointers being invalidated!
  Data::Clusters.reserve( aData.size() );

  // ... any other preparation ...
}
