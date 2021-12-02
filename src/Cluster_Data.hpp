#pragma once

/* ===== C++ ===== */
#include <vector>
#include <list>
#include <math.h>
#include <stdlib.h>
#include <iostream>

/* ===== Cluster sources ===== */
#include "Cluster_Sigmas.hpp"

/* ===== Struct for storing data ===== */
class Data
{
public:
  Data( const double& aX , const double& aY , const double& aS ) : x(aX) , y(aY) , s(aS) , r( sqrt( (aX*aX) + (aY*aY) ) ), phi( atan2( aY , aX ) ),
  eX( 1 - fabs( aX ) ) , eY( 1 - fabs( aY ) ) , 
  w( [ aS ]( const double& sig2 ){ return 1 / ( (aS*aS) + sig2 ); } | sigmabins2 ), log_w( []( const double& w){ return log(w); } | w ),
  localizationsum( 0.0 ) , localizationscore( 0.0 ),
  neighbourit( neighbours.end() ),
  parent( NULL )
  {}

  bool operator< ( const Data& aOther ) const { return r < aOther.r; }

  inline double dR2( const Data& aOther ) const
  {
    double dX( x - aOther.x ), dY( y - aOther.y );
    return ( dX*dX ) + ( dY*dY );
  }

  inline double dR( const Data& aOther ) const
  {
    return sqrt( dR2( aOther ) );
  }

  inline void UpdateLocalization( const double& aR2 , const size_t& Nminus1  )
  {
    constexpr double pi = atan(1)*4;

    const double last_localizationsum( localizationsum );

    for( ; neighbourit != neighbours.end() ; ++neighbourit )
    { 
      if( neighbourit->first > aR2 ) break;
      double lDist = sqrt( neighbourit->first );
      double Weight( 1.0 );
      if( eX < lDist )  Weight *= ( 1 + pow( acos( eX/lDist ) * (2/pi) , 4) );
      if( eY < lDist )  Weight *= ( 1 + pow( acos( eY/lDist ) * (2/pi) , 4) );
      localizationsum += Weight;
    }

    if( last_localizationsum == localizationsum ) return;

    const double LocalizationConstant( 4.0 / ( pi * Nminus1 ) ); 
    localizationscore = sqrt( LocalizationConstant * localizationsum );
  }


  inline void Clusterize( const double& a2R2 , const double& aT )
  {
    if( parent ) return;
    if( localizationscore < aT ) return;

    Data* lParent( this );
    
    for( auto& j : neighbours )
    {
      if( j.first > a2R2 ) break;
      if( (!j.second->parent) or (j.second->localizationscore < aT) ) continue;
      lParent = j.second->GetParent();
      goto breakpoint;
    }

    for( auto& j : neighbours2 )
    {
      if( j.first > a2R2 ) break;
      if( (!j.second->parent) or (j.second->localizationscore < aT) ) continue;
      lParent = j.second->GetParent();
      goto breakpoint;
    }

    SetParent( lParent );

    breakpoint:

    for( auto& j : neighbours )
    {
      if( j.first > a2R2 ) break;
      if( j.second->localizationscore < aT) continue;
      j.second->SetParent( lParent );
    }

    for( auto& j : neighbours2 )
    {
      if( j.first > a2R2 ) break;
      if( j.second->localizationscore < aT) continue;
      j.second->SetParent( lParent );
    }

  }


  inline Data* GetParent()
  {
    if( parent == this or parent == NULL ) return parent;
    return parent->GetParent();
  }

  inline void SetParent( Data* aParent )
  {
    if( parent == aParent ) return;
    if( parent == this or parent == NULL ) parent = aParent;
    else parent->SetParent( aParent );
  }

  // void Cluster( Data& aData )
  // {
  //   SetParent( aData.GetParent() );
  // }

public:
  double x, y, s , r, phi;
  double eX , eY;
  std::vector< double > w , log_w;
  double localizationsum , localizationscore;

  std::vector< std::pair< double , Data* > > neighbours;
  std::vector< std::pair< double , Data* > >::iterator neighbourit;

  std::vector< std::pair< double , Data* > > neighbours2;
  Data* parent;
};