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

// Header file for class GSLInterpolator
//
// Created by: moneta  at Fri Nov 26 15:31:41 2004
//
// Last update: Fri Nov 26 15:31:41 2004
//
#pragma once

#include <vector>
#include <string>
#include <cassert>
#include <iostream>
#include <functional>

#include "gsl/gsl_interp.h"
#include "gsl/gsl_spline.h"
#include "gsl/gsl_errno.h"

/**
Interpolation class based on GSL interpolation functions
 @ingroup Interpolation
 */

class GSLInterpolator {

public:

   GSLInterpolator( const gsl_interp_type *type , const unsigned int& ndata );

   GSLInterpolator( const gsl_interp_type *type, const std::vector<double> & x, const std::vector<double> & y );
   virtual ~GSLInterpolator();

   GSLInterpolator(const GSLInterpolator &) = delete;
   GSLInterpolator & operator = (const GSLInterpolator &) = delete;

public:

   inline bool SetData( const std::vector<double> & x , const std::vector<double> & y )
   {
      assert( x.size() == y.size() );
      return SetData( x.size() , &x.front() , &y.front() );
   }

   bool SetData( const unsigned int& ndata, const double *x, const double *y );


   inline double Evaluate( const std::function< int( double& ) >& aFunction , const std::string& aName )
   {
      assert(fAccel);

      double RetVal(0);
      auto ierr = aFunction( RetVal );
      if ( ierr ) {
         ++nErrors;
         if(nErrors <= 4) {
            std::cout << "GSLInterpolator::" << aName << ": " << gsl_strerror(ierr) << std::endl;
            if(nErrors == 4) std::cout << "GSLInterpolator::" << aName << ": " << "Suppressing additional warnings" << std::endl;
         }
      }
      return RetVal;
   }

   inline double Eval( const double& x )
   {
      return Evaluate( [&]( double& RetVal ){ return gsl_spline_eval_e( fSpline, x, fAccel, &RetVal ); } , "Eval" );
   }

   inline double Deriv( const double& x )
   {
      return Evaluate( [&]( double& RetVal ){ return gsl_spline_eval_deriv_e(fSpline, x, fAccel, &RetVal); } , "Deriv" );
   }

   inline double Deriv2( const double& x )
   {
      return Evaluate( [&]( double& RetVal ){ return gsl_spline_eval_deriv2_e(fSpline, x, fAccel, &RetVal); } , "Deriv2" );
   }

   inline double Integ( const double&  a, const double&  b )
   {
      if (a > b) return -Integ(b, a);
      return Evaluate( [&]( double& RetVal ){ return gsl_spline_eval_integ_e(fSpline, a, b, fAccel, &RetVal); } , "Integ" );
   }

private:
   unsigned int nErrors;
   gsl_interp_accel * fAccel;
   gsl_spline * fSpline;
   const gsl_interp_type * fInterpType;

};
