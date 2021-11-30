#pragma once

/* ===== C++ ===== */
#include <math.h>
#include <stdlib.h>

/* ===== Struct for storing data ===== */
class Data
{
public:
  Data( double aX , double aY ) : x(aX) , y(aY) , r( sqrt( (aX*aX) + (aY*aY) ) ), phi( atan2( aY , aX ) ), parent( NULL ){}

  bool operator< ( const Data& aOther ) const { return r < aOther.r; }

  double dR2( const Data& aOther ) const
  {
    double dX( x - aOther.x ), dY( y - aOther.y );
    return ( dX*dX ) + ( dY*dY );
  }

  double dR( const Data& aOther ) const
  {
    return sqrt( dR2( aOther ) );
  }

  Data* GetParent()
  {
    if( parent == NULL ) return NULL;    
    if( parent == this ) return this;
    return parent->GetParent();
  }

  void SetParent( Data* aParent )
  {
    if( parent == this ) parent = aParent;
    else parent->SetParent( aParent );
  }

  void Cluster( Data& aData )
  {
    SetParent( aData.GetParent() );
  }

public:
  double x, y, r, phi;
  Data* parent;
};