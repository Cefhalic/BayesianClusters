// Andrew W. Rose, 2022
// based on an original implementation by
// Authors: L. Moneta, A. Zsenei 08/2005, Copyright (c) 2004 ROOT Foundation, CERN/PH-SFT 

#pragma once

#include <vector>
#include <string>
#include <cassert>
#include <iostream>
#include <functional>

#include "gsl/gsl_interp.h"
#include "gsl/gsl_spline.h"
#include "gsl/gsl_errno.h"

class GSLInterpolator {

public:

   GSLInterpolator( const gsl_interp_type *type , const unsigned int& ndata );

   GSLInterpolator( const gsl_interp_type *type, const std::vector<double> & x, const std::vector<double> & y );
   virtual ~GSLInterpolator();

   GSLInterpolator(const GSLInterpolator &) = delete;
   GSLInterpolator& operator= (const GSLInterpolator &) = delete;
   GSLInterpolator& operator= ( GSLInterpolator && ) = default;

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
      return Evaluate( [&]( double& RetVal ){ return gsl_spline_eval_deriv_e( fSpline, x, fAccel, &RetVal); } , "Deriv" );
   }

   inline double Deriv2( const double& x )
   {
      return Evaluate( [&]( double& RetVal ){ return gsl_spline_eval_deriv2_e( fSpline, x, fAccel, &RetVal); } , "Deriv2" );
   }

   inline double Integ( const double&  a, const double&  b )
   {
      if (a > b) return -Integ(b, a);
      return Evaluate( [&]( double& RetVal ){ return gsl_spline_eval_integ_e( fSpline, a, b, fAccel, &RetVal); } , "Integ" );
   }

private:
   unsigned int nErrors;
   gsl_interp_accel * fAccel;
   gsl_spline * fSpline;
   const gsl_interp_type * fInterpType;

};
