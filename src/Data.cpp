
/* ===== C++ ===== */
#include <stdlib.h>
#include <algorithm>

/* ===== Cluster sources ===== */
#include "BayesianClustering/Data.hpp"
#include "BayesianClustering/Cluster.hpp"
#include "BayesianClustering/Event.hpp"


// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Data::Data( const PRECISION& aX , const PRECISION& aY , const PRECISION& aS ) : 
x(aX) , y(aY) , s(aS) , r2( (aX*aX) + (aY*aY) ), r( sqrt( r2 ) ), phi( atan2( aY , aX ) ),
mProtoCluster( NULL )
{}

Data::~Data()
{
  if( mProtoCluster ) delete mProtoCluster;
  mProtoCluster = NULL;
}

// Although the neighbourhood calculation is reciprocal (if I am your neighbour then you are mine) and we can, in fact, use that to halve the number of calculations,
// doing so requires arbitration between threads or a single-threaded reciprocation step, both of which take longer than brute-forcing it
__attribute__((flatten))
void Data::Preprocess( std::vector<Data>& aData , const std::size_t& aIndex )
{
  static constexpr double pi = atan(1)*4;  
  auto dphi = asin( Event::mParameters.max2R() / r ); // Event::mParameters.max2R() / ( r - aParameters.max2R() );
  auto dphi2 = (2*pi) - dphi;

  std::size_t i( aIndex + 1 );
  std::vector<Data>::const_iterator aPlusIt( aData.begin() + aIndex + 1 );
  const std::vector<Data>::const_iterator aPlusEnd( aData.end() );

  // Iterate over other hits and populate the mNeighbour list
  for( ; aPlusIt != aPlusEnd ; ++aPlusIt , ++i )
  {
    if( ( aPlusIt->r - r ) > Event::mParameters.max2R() ) break; // aPlusIt is always further out than curent 
    auto lPhi = dPhi( *aPlusIt );
    if( lPhi > dphi and lPhi < dphi2 ) continue;
    PRECISION ldR2 = dR2( *aPlusIt );
    if( ldR2 < Event::mParameters.max2R2() ) mNeighbours.push_back( std::make_pair( ldR2 , i ) );
  }

  i = aIndex - 1;
  std::vector<Data>::const_reverse_iterator aMinusIt( aData.rbegin() + aData.size() - aIndex );
  const std::vector<Data>::const_reverse_iterator aMinusEnd( aData.rend() );

  for( ; aMinusIt != aMinusEnd ; ++aMinusIt , --i )
  {
    if( ( r - aMinusIt->r ) > Event::mParameters.max2R() ) break; // curent is always further out than aMinusIn
    auto lPhi = dPhi( *aMinusIt );
    if( lPhi > dphi and lPhi < dphi2 ) continue;
    PRECISION ldR2 = dR2( *aMinusIt );    
    if( ldR2 < Event::mParameters.max2R2() ) mNeighbours.push_back( std::make_pair( ldR2 , i ) );
  }

  std::sort( mNeighbours.begin() , mNeighbours.end() );

  // -------------------------------------------------------------------------------------

  mProtoCluster = new Cluster( *this );

  // -------------------------------------------------------------------------------------

  const double lLocalizationConstant( 4.0 / ( pi * ( aData.size() - 1 ) ) ); 
  const PRECISION eX( 1 - fabs( x ) ) , eY( 1 - fabs( y ) );

  auto lNeighbourit( mNeighbours.begin() );
  PRECISION lLocalizationSum( 0 ) , lLastLocalizationSum( 0 ) , lLocalizationScore( 0 ) , lDist( 0 ) , lWeight( 0 );
  mLocalizationScores.reserve( Event::mParameters.Rbins() );

  double R( Event::mParameters.minScanR() ) , R2( 0 );
  for( uint32_t i(0) ; i!=Event::mParameters.Rbins() ; ++i , R+=Event::mParameters.dR() )
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
      lLastLocalizationSum = lLocalizationSum;
    }

    mLocalizationScores.push_back( lLocalizationScore );
  }
}
// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


