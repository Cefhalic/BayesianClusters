// @(#)root/mathmore:$Id$
// Authors: L. Moneta, A. Zsenei   08/2005

 /**********************************************************************
  *                                                                    *
  * Copyright (c) 2004 ROOT Foundation,  CERN/PH-SFT                   *
  *                                                                    *
  * This library is free software; you can redistribute it and/or      *
  * modify it under the terms of the GNU General Public License        *
  * as published by the Free Software Foundation; either version 2     *
  * of the License, or (at your option) any later version.             *
  *                                                                    *
  * This library is distributed in the hope that it will be useful,    *
  * but WITHOUT ANY WARRANTY; without even the implied warranty of     *
  * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU   *
  * General Public License for more details.                           *
  *                                                                    *
  * You should have received a copy of the GNU General Public License  *
  * along with this library (see file COPYING); if not, write          *
  * to the Free Software Foundation, Inc., 59 Temple Place, Suite      *
  * 330, Boston, MA 02111-1307 USA, or contact the author.             *
  *                                                                    *
  **********************************************************************/

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


bool  GSLInterpolator::SetData( const unsigned int& size , const double *x , const double *y ) {
   // initialize interpolation object with the given data
   // if given size is different a new interpolator object is created
   if ( !fSpline ) fSpline = gsl_spline_alloc( fInterpType , size );
   else if ( size != fSpline->interp->size )
   {
      gsl_spline_free(fSpline);
      fSpline = gsl_spline_alloc( fInterpType, size);
   }
   if ( !fSpline ) return false;

   int iret = gsl_spline_init( fSpline , x , y , size );
   if ( !iret ) return false;

   if( !fAccel ) fAccel = gsl_interp_accel_alloc() ;
   else          gsl_interp_accel_reset(fAccel);

   assert ( fSpline );
   assert ( fAccel );
   nErrors = 0;   // reset counter for error messages
   return true;
}

GSLInterpolator::~GSLInterpolator()
{
   // free gsl objects
   if ( fSpline ) gsl_spline_free(fSpline);
   if ( fAccel ) gsl_interp_accel_free( fAccel);
}

