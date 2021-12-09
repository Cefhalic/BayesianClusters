

/* ===== C++ ===== */
#include <math.h>
#include <stdlib.h>
// #include <iostream>

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



Data::ClusterParameter::ClusterParameter( const PRECISION& aW ) : 
w( aW )
{}
    
void Data::ClusterParameter::Reset( const PRECISION& aX , const PRECISION& aY , const PRECISION& aR2 )
{
  A = w;
  Bx = w * aX;
  By = w * aY;
  C = w * aR2;
  sum_logw = log( w );
}

Data::ClusterParameter& Data::ClusterParameter::operator+= ( const Data::ClusterParameter& aOther )
{
  A += aOther.A;
  Bx += aOther.Bx;
  By += aOther.By;
  C += aOther.C;
  sum_logw += aOther.sum_logw;
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
double Data::ClusterParameter::log_score( const std::size_t& n ) const
{
  static constexpr double pi = atan(1)*4;
  static constexpr double log2pi = log( 2*pi );

  auto sqrt_A( sqrt( A ) ) , inv_A( 1.0 / A );
  auto Dx( Bx * inv_A ) , Dy( By * inv_A );
  auto E( C - ( Bx * Dx ) - ( By * Dy ) );

  double log_sum = double( -log( 4.0 ) ) // Area term
                   + (log2pi * (1.0-n))
                   + sum_logw
                   - double( log( A ) ) // A is equivalent to n_tilde
                   + (0.5 * E);

  // We place explicit bounds checks to prevent calls to expensive functions
  auto arg1 = CDF( sqrt_A * (1.0-Dx) ) - CDF( sqrt_A * (-1.0-Dx) );
  if( arg1 != 1.0 ) log_sum += log( arg1 );
  auto arg2 = CDF( sqrt_A * (1.0-Dy) ) - CDF( sqrt_A * (-1.0-Dy) );
  if( arg2 != 1.0 ) log_sum += log( arg2 );

  return log_sum;
}


Data::Data( const PRECISION& aX , const PRECISION& aY , const PRECISION& aS ) : 
x(aX) , y(aY) , r2( (aX*aX) + (aY*aY) ), r( sqrt( r2 ) ), phi( atan2( aY , aX ) ),
localizationsum( 0.0 ) , localizationscore( 0.0 ),
neighbourit( neighbours.end() ),
parent( NULL ),
ClusterSize( 0 ) , ClusterScore( 0.0 )
{
  for( auto& sig2 : Parameters.sigmabins2() ) ClusterParams.emplace_back( 1.0 / ( (aS*aS) + sig2 ) );

  neighbours.reserve( 1024 );    
}


bool Data::operator< ( const Data& aOther ) const
{ 
  return r < aOther.r; 
}

PRECISION Data::dR2( const Data& aOther ) const
{
  PRECISION dX( x - aOther.x ), dY( y - aOther.y );
  return ( dX*dX ) + ( dY*dY );
}

PRECISION Data::dR( const Data& aOther ) const
{
  return sqrt( dR2( aOther ) );
}


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
  parent = NULL;
  LastClusterSize = 0;
  ClusterSize = 1;
  ClusterScore = 0.0;
  for( auto& i : ClusterParams ) i.Reset( x , y , r2 );
}


__attribute__((flatten))
void Data::Clusterize( const PRECISION& a2R2 , const PRECISION& aT , const Data* aLower , const Data* aUpper )
{
  if( parent ) return;
  if( localizationscore < aT ) return;

  for( auto& j : neighbours )
  {
    if( j.first > a2R2 ) break;
    if( j.second < aLower or j.second >= aUpper ) continue; // Save us from having to mutex the parent
    if( j.second->localizationscore > aT ) j.second->ClusterInto( this );
  }
}


// __attribute__((flatten))
// void Data::Clusterize2( const PRECISION& a2R2 , const PRECISION& aT )
// {
//   if( localizationscore < aT ) return;

//   for( auto& j : neighbours )
//   {
//     if( j.first > a2R2 ) break;
//     if( j.second->localizationscore > aT ) j.second->ClusterInto( this , false );
//   }
// }

#include <iostream>

void Data::ClusterInto( Data* aParent )
{
  if ( parent == aParent ) return; // If the proposed parent is already our parent, return
  if ( parent ) return parent->ClusterInto( aParent ); // We already belong to another - they already have all our clusterparams - Hand them over to claiming cluster and return

  // Set the claiming cluster as parent
  parent = aParent;

  // Add our cluster params into the parent's  
  parent->ClusterSize += ClusterSize;
  ClusterSize = 0;
  for( auto lIt( ClusterParams.begin() ) , lIt2( parent->ClusterParams.begin() ) ; lIt != ClusterParams.end() ; ++lIt , ++lIt2 ) *lIt2 += *lIt;
}


Data* Data::GetParent()
{
  if( !parent ) return this;
  return parent = parent->GetParent();
}





void Data::UpdateClusterScore()
{
  if( ClusterSize <= LastClusterSize ) return; // We were not bigger than the previous size when we were evaluated - score is still valid
  LastClusterSize = ClusterSize;

  thread_local static std::vector< double > integral( Parameters.sigmacount() , 1.0 );

  double constant;

  for( std::size_t i(0) ; i!=Parameters.sigmacount() ; ++i )
  {
    double log_sum = ClusterParams[i].log_score( ClusterSize ) + Parameters.log_probability_sigma( i );
    if( i == 0 ) constant = log_sum;
    else         integral[i] = exp( log_sum - constant );
  }

  thread_local static ROOT::Math::Interpolator lInt( Parameters.sigmacount() , ROOT::Math::Interpolation::kLINEAR );
  lInt.SetData( Parameters.sigmabins() , integral );

  static const double Lower( Parameters.sigmabins(0) ) , Upper( Parameters.sigmabins(Parameters.sigmacount()-1) );
  ClusterScore = PRECISION( log( lInt.Integ( Lower , Upper ) ) ) + constant;  
}






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

    []( Data& i ){ i.ResetClusters(); } || aData;

    for( uint32_t j(0) ; j!=Parameters.Tbins() ; ++j , T-=Parameters.dT() , ++lProgressBar )
    {
      Clusterize( aData , twoR2 , T );   
      []( Data& i ){ i.UpdateClusterScore(); } || aData; // Use interleaving threading to average over systematic radial scaling
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

  // ... any other preparation ...
}
