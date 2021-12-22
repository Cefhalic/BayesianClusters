

/* ===== C++ ===== */
#include <stdlib.h>
#include <iostream>
#include <iomanip>

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



Cluster::Cluster(): mParams( Parameters.sigmacount() ),
mClusterSize( 0 ) , mLastClusterSize( 0 ) , mClusterScore( 0.0 ) , 
mParent( NULL )
{}


double Cluster::log_score()
{
  static constexpr double pi = atan(1)*4;
  static constexpr double log2pi = log( 2*pi );

  if( mClusterSize <= mLastClusterSize ) return mClusterScore; // We were not bigger than the previous size when we were evaluated - score is still valid
  mLastClusterSize = mClusterSize;

  thread_local static std::vector< double > MuIntegral( Parameters.sigmacount() , 1.0 );

  // double constant( mParams[0].log_score() + Parameters.log_probability_sigma( 0 ) );
  // for( std::size_t i(1) ; i!=Parameters.sigmacount() ; ++i ) MuIntegral[i] = exp( mParams[i].log_score() + Parameters.log_probability_sigma( i ) - constant );
  for( std::size_t i(0) ; i!=Parameters.sigmacount() ; ++i ) MuIntegral[i] = exp( mParams[i].log_score() + Parameters.log_probability_sigma( i ) );

  thread_local static ROOT::Math::Interpolator lInt( Parameters.sigmacount() , ROOT::Math::Interpolation::kLINEAR );
  lInt.SetData( Parameters.sigmabins() , MuIntegral );

  static const double Lower( Parameters.sigmabins(0) ) , Upper( Parameters.sigmabins(Parameters.sigmacount()-1) );
  //return mClusterScore = double( log( lInt.Integ( Lower , Upper ) ) ) + constant - double( log( 4.0 ) ) + (log2pi * (1.0-mClusterSize));  
  return mClusterScore = double( log( lInt.Integ( Lower , Upper ) ) ) - double( log( 4.0 ) ) + (log2pi * (1.0-mClusterSize));  
}


void Cluster::operator+= ( const Data& aData )
{
  auto w( aData.mWeights.begin() );
  for( auto lIt( mParams.begin() ) ; lIt != mParams.end() ; ++lIt , ++w )
  {
    lIt->A += *w;
    lIt->Bx += (*w * aData.x);
    lIt->By += (*w * aData.y);
    lIt->C += (*w * aData.r2);
    lIt->logF += PRECISION( log( *w ) );
  }
  mClusterSize += 1;
}

void Cluster::operator+= ( Cluster& aOther )
{
  // if( &aOther == this ) throw std::runtime_error( "Error #1" );
  // if( aOther.mClusterSize == 0 ) throw std::runtime_error( "Error #2" );
  // if( aOther.mParent ) throw std::runtime_error( "Error #3" );

  for( auto lIt( mParams.begin() ) , lIt2( aOther.mParams.begin() ) ; lIt != mParams.end() ; ++lIt , ++lIt2 ) *lIt += *lIt2;
  mClusterSize += aOther.mClusterSize;
  aOther.mClusterSize = 0;
  aOther.mParent = this;
}


Cluster* Cluster::GetParent()
{
  if( mParent ) return mParent = mParent->GetParent();
  return this;
}



Data::Data( const PRECISION& aX , const PRECISION& aY , const PRECISION& aS ) : 
x(aX) , y(aY) , r2( (aX*aX) + (aY*aY) ), r( sqrt( r2 ) ), phi( atan2( aY , aX ) ),
mWeights( [ &aS ]( const double& sig2 ){ return PRECISION( 1.0 / ( (aS*aS) + sig2 ) ); } | Parameters.sigmabins2() ),
mLocalizationSum( 0.0 ) , mLocalizationScore( 0.0 ),
mNeighbourit( mNeighbours.end() ),
mCluster( NULL )
{}





__attribute__((flatten))
void Data::PopulateNeighbours( std::vector<Data>::iterator aPlusIt , const std::vector<Data>::iterator& aPlusEnd , std::vector<Data>::reverse_iterator aMinusIt , const std::vector<Data>::reverse_iterator& aMinusEnd )
{
  const auto dphi = Parameters.maxR() / ( r - Parameters.maxR() );
  const auto dphi2 = Parameters.max2R() / ( r - Parameters.max2R() );

  // Iterate over other hits and populate the mNeighbour list
  for( ; aPlusIt != aPlusEnd ; aPlusIt++ )
  {
    if( ( aPlusIt->r - r ) > Parameters.max2R() ) break; // aPlusIt is always further out than curent 
    if( fabs( aPlusIt->phi - phi ) > dphi ) continue;
    PRECISION ldR2 = dR2( *aPlusIt );
    if( ldR2 < Parameters.max2R2() ) mNeighbours.push_back( std::make_pair( ldR2 , &*aPlusIt ) );
  }

  for( ; aMinusIt != aMinusEnd ; aMinusIt++ )
  {
    if( ( r - aMinusIt->r ) > Parameters.max2R() ) break; // curent is always further out than aMinusIn
    if( fabs( aPlusIt->phi - phi ) > dphi ) continue;
    PRECISION ldR2 = dR2( *aMinusIt );    
    if( ldR2 < Parameters.max2R2() ) mNeighbours.push_back( std::make_pair( ldR2 , &*aMinusIt ) );
  }

  std::sort( mNeighbours.begin() , mNeighbours.end() );
  mNeighbourit = mNeighbours.begin();
}


void Data::UpdateLocalization( const PRECISION& aR2 , const size_t& Nminus1  )
{
  constexpr PRECISION pi = atan(1)*4;

  const PRECISION lLastLocalizationSum( mLocalizationSum );

  for( ; mNeighbourit != mNeighbours.end() ; ++mNeighbourit )
  { 
    if( mNeighbourit->first > aR2 ) break;
    PRECISION lDist = sqrt( mNeighbourit->first );

    // Approximation of the edge-correction
    PRECISION Weight( 1.0 );
    PRECISION eX( 1 - fabs( x ) ) , eY( 1 - fabs( y ) );
    if( eX < lDist )  Weight *= ( 1 + pow( acos( eX/lDist ) * (2/pi) , 4) );
    if( eY < lDist )  Weight *= ( 1 + pow( acos( eY/lDist ) * (2/pi) , 4) );

    mLocalizationSum += Weight;
  }

  if( lLastLocalizationSum == mLocalizationSum ) return;

  const PRECISION LocalizationConstant( 4.0 / ( pi * Nminus1 ) ); 
  mLocalizationScore = sqrt( LocalizationConstant * mLocalizationSum );
}



// We are at the top-level
__attribute__((flatten))
void Data::Clusterize( const PRECISION& a2R2 , const PRECISION& aT , std::vector< Cluster >& aClusters )
{
  if( mCluster ) return;
  if( mLocalizationScore < aT ) return;

  // if one of our mNeighbours is already a cluster, join that
  for( auto& j : mNeighbours )
  {
    if( j.first > a2R2 ) break;
    if( ! j.second->mCluster ) continue;  
    return Clusterize( a2R2 , aT , j.second->mCluster->GetParent() );
  }

  // else create a new cluster
  aClusters.emplace_back();
  Clusterize( a2R2 , aT , &aClusters.back() );
}



void Data::Clusterize( const PRECISION& a2R2 , const PRECISION& aT , Cluster* aCluster )
{
  if( mCluster )
  {
    if( ( mCluster = mCluster->GetParent() ) == aCluster ) return;
    *aCluster += *mCluster;
    mCluster = aCluster;
  }
  else
  {
    if( mLocalizationScore < aT ) return;
    *aCluster += *this;
    mCluster = aCluster;

    for( auto& i : mNeighbours )
    {
      if( i.first > a2R2 ) break;
      i.second->Clusterize( a2R2 , aT , aCluster );
    }  
  }

}



bool CheckClusterization( std::vector<Data>& aData , const double& R , const double& T )
{
  const auto lRlimit = 4.0 * R * R;

  uint32_t lNotClustered( 0 );
  uint32_t lNeighbourNotClustered( 0 );
  uint32_t lWrongNeighbour( 0 );

  for( auto& i : aData )
  {
    if( i.mLocalizationScore > T )
    {
      if( ! i.mCluster ){ lNotClustered++ ; continue; }
      i.mCluster = i.mCluster->GetParent(); // Update cache
      for( auto& j : i.mNeighbours )
      {
        if( j.first > lRlimit ) break;
        if( j.second->mLocalizationScore < T ) continue;

        if( ! j.second->mCluster ){ lNeighbourNotClustered++; continue; }
        if ( j.second->mCluster->GetParent() != i.mCluster ){ lWrongNeighbour++; continue; }
      }
    }
  }

  if( lNotClustered or lNeighbourNotClustered or lWrongNeighbour )
  {
    std::cout << "Not Clustered = " << std::setw(8) << lNotClustered << " | Neighbour Not Clustered = " << std::setw(8) << lNeighbourNotClustered << " | Wrong Neighbour = " << std::setw(8) << lWrongNeighbour << std::endl;
    return false;
  }

  return true;
}




__attribute__((flatten))
void Clusterize( std::vector<Data>& aData , std::vector< Cluster >& aClusters , const double& R , const double& T )
{
  ProgressBar2 lProgressBar( "Clustering" , 0 );

  const std::size_t N( aData.size()-1 );
  // Reset ahead of UpdateLocalization and Clusterize
  []( Data& i ){ i.mCluster = NULL; i.mLocalizationSum = i.mLocalizationScore = 0.0; i.mNeighbourit = i.mNeighbours.begin(); } || aData;

  // Update the localization score
  [ &R , &N ]( Data& i ){ i.UpdateLocalization( R * R , N ); } || aData; // Use interleaving threading to average over systematic radial scaling

  // And clusterize
  auto twoR2 = 4.0*R*R;
  for( auto& i : aData ) i.Clusterize( twoR2 , T , aClusters );
}


__attribute__((flatten))
void ScanRT( std::vector<Data>& aData , const std::function< void( const std::vector< Cluster >& , const double& , const double& ) >& aCallback )
{
  const std::size_t Nminus1( aData.size() - 1 );
  double R( Parameters.minScanR() ) , R2( 0 ) , twoR2( 0 ) , T( 0 );

  std::vector< Cluster > lClusters;
  lClusters.reserve( aData.size() );  // Reserve as much space for clusters as there are data points - prevent pointers being invalidated!

  auto UpdateLocalizationExpr = [ &Nminus1 , &R2 ]( Data& i ){ i.UpdateLocalization( R2 , Nminus1 ); };
  auto ResetClustersExpr = []( Data& i ){ i.mCluster = NULL; };
  auto UpdateScoreExpr = []( Cluster& i ){ if( i.mClusterSize ) i.log_score(); };

  ProgressBar lProgressBar( "Scan over RT" , Parameters.Rbins() * Parameters.Tbins() );

  for( uint32_t i(0) ; i!=Parameters.Rbins() ; ++i , R+=Parameters.dR() )
  {
    R2 = R * R;
    twoR2 = 4.0 * R2;
    T = Parameters.maxScanT();

    UpdateLocalizationExpr || aData; // Use interleaving threading to average over systematic radial scaling

    lClusters.clear();
    ResetClustersExpr || aData;

    for( uint32_t j(0) ; j!=Parameters.Tbins() ; ++j , T-=Parameters.dT() , ++lProgressBar )
    {
      for( auto& i : aData ) i.Clusterize( twoR2 , T , lClusters );
      UpdateScoreExpr || lClusters; 
      aCallback( lClusters , R , T );
    }
  }

  lClusters.clear();
  ResetClustersExpr || aData; // Delete cluster pointers which will be invalidated when we leave the function

}



void PrepData( std::vector<Data>& aData )
{
  {
    // Populate mNeighbour lists  
    ProgressBar2 lProgressBar( "Populating neighbourhood" , aData.size() );
    [ &aData ]( const std::size_t& i ){ aData.at( i ).PopulateNeighbours( aData.begin() + i + 1 , aData.end() , aData.rbegin() + aData.size() - i , aData.rend() ); } || range( aData.size() );  // Interleave threading since processing time increases with radius from origin
  }
}
