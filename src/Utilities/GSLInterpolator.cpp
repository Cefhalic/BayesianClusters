//! \file GSLInterpolator.cpp

// Andrew W. Rose, 2022
// based on an original implementation by
// Authors: L. Moneta, A. Zsenei 08/2005, Copyright (c) 2004 ROOT Foundation, CERN/PH-SFT 

#include "Utilities/GSLInterpolator.hpp"

GSLInterpolator::GSLInterpolator ( const gsl_interp_type *type , const unsigned int& size ) :
   nErrors(0),
   fAccel(nullptr),
   fSpline(nullptr),
   fInterpType(type)
{
   if (size >= fInterpType->min_size) fSpline = gsl_spline_alloc( fInterpType , size );
}


GSLInterpolator::GSLInterpolator( const gsl_interp_type *type, const std::vector<double> & x, const std::vector<double> & y ) :
   nErrors(0),
   fAccel(nullptr),
   fSpline(nullptr),
   fInterpType(type)
{
   assert( x.size() == y.size() );
   auto size = x.size();

   if (size >= fInterpType->min_size) fSpline = gsl_spline_alloc( fInterpType , size );
   SetData( size , &x.front() , &y.front() );
}

GSLInterpolator::~GSLInterpolator()
{
   if ( fSpline ) gsl_spline_free( fSpline );
   if ( fAccel ) gsl_interp_accel_free( fAccel );
}

bool  GSLInterpolator::SetData( const unsigned int& size , const double *x , const double *y )
{
   if ( !fSpline ) fSpline = gsl_spline_alloc( fInterpType , size );
   else if ( size != fSpline->interp->size )
   {
      gsl_spline_free(fSpline);
      fSpline = gsl_spline_alloc( fInterpType, size);
   }
   if ( !fSpline ) return false;

   int iret = gsl_spline_init( fSpline , x , y , size );
   if ( iret ) return false;

   if( !fAccel ) fAccel = gsl_interp_accel_alloc() ;
   else          gsl_interp_accel_reset(fAccel);

   assert ( fSpline );
   assert ( fAccel );
   nErrors = 0;   // reset counter for error messages
   return true;
}


