#pragma once

/* ===== Cluster sources ===== */
#include "Cluster_Data.hpp"

/* ===== C++ ===== */
#include <math.h>


/* ===== Approximation to correct area of circle overlapping the edge of the square ===== */
double EdgeCorrectedWeight( const Data& aPt , const double& aDist )
{
  constexpr double pi = atan(1)*4;

  double Weight( 1.0 );
  const double X( 1 - fabs( aPt.x ) ) , Y( 1 - fabs( aPt.y ) );
  if( X < aDist )  Weight *= ( 1 + pow( acos( X/aDist ) * (2/pi) , 4) );
  if( Y < aDist )  Weight *= ( 1 + pow( acos( Y/aDist ) * (2/pi) , 4) );
  return Weight;
}
