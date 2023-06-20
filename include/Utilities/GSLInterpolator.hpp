//! \file GSLInterpolator.hpp
// Andrew W. Rose, 2022
// based on an original implementation by
// Authors: L. Moneta, A. Zsenei 08/2005, Copyright (c) 2004 ROOT Foundation, CERN/PH-SFT

#pragma once

#include <map>
#include <vector>
#include <string>
#include <cassert>
#include <iostream>
#include <functional>

#include "gsl/gsl_interp.h"
#include "gsl/gsl_spline.h"
#include "gsl/gsl_errno.h"

//! A utility wrapper around the GSL interpolator to give it a clean C++ interface
class GSLInterpolator
{

  public:

    //! Empty splice constructor
    //! \param type The spline type
    //! \param ndata The number of points that will be added to the spline
    GSLInterpolator( const gsl_interp_type* type, const unsigned int& ndata );

    //! Initialised splice constructor
    //! \param type The spline type
    //! \param x The points on the x-axis
    //! \param y The points on the y-axis
    GSLInterpolator( const gsl_interp_type* type, const std::vector<double>& x, const std::vector<double>& y );

    //! Initialised splice constructor
    //! \param type The spline type
    //! \param data Data points along the spline
    GSLInterpolator( const gsl_interp_type* type, const std::map<double,double>& data );

    //! Destructor
    virtual ~GSLInterpolator();

    //! Deleted copy constructor
    GSLInterpolator(const GSLInterpolator& aOther /*!< Anonymous argument */) = delete;

    //! Deleted assignment operator
    //! \return Reference to this, for chaining calls
    GSLInterpolator& operator= (const GSLInterpolator& aOther /*!< Anonymous argument */) = delete;

    //! Default move constructor
    GSLInterpolator( GSLInterpolator&& aOther /*!< Anonymous argument */ ) = default;

    //! Default move-assignment constructor
    //! \return Reference to this, for chaining calls
    GSLInterpolator& operator= ( GSLInterpolator&& aOther /*!< Anonymous argument */ ) = default;

  public:

    //! Set the spline data points
    //! \param x The x-coordinates of the datapoints
    //! \param y The y-coordinates of the datapoints
    //! \return success or fail
    inline bool SetData( const std::vector<double>& x, const std::vector<double>& y )
    {
      assert( x.size() == y.size() );
      return SetData( x.size(), &x.front(), &y.front() );
    }

    //! Set the spline data points
    //! \param ndata The number of data points
    //! \param x Pointer to the first element of an array of x-coordinates
    //! \param y Pointer to the first element of an array of y-coordinates
    //! \return success or fail
    bool SetData( const unsigned int& ndata, const double* x, const double* y );

    //! Utility function that runs the GSL function that has been wrapped in a lambda below
    //! \param aFunction A lambda that will be evaluated
    //! \param aName The operation name for the debugging messages
    //! \return The interpolated value
    inline double Evaluate( const std::function< int( double& ) >& aFunction, const std::string& aName )
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

    //! Evaluate the spline at the given x
    //! \param x The x-coordinate at which to evaluate the spline
    //! \return The value of the spline at the given x-coordinate
    inline double Eval( const double& x )
    {
      return Evaluate( [&]( double& RetVal ) {
        return gsl_spline_eval_e( fSpline, x, fAccel, &RetVal );
      }, "Eval" );
    }

    //! The first derivative of the spline at the given x
    //! \param x The x-coordinate at which to evaluate the derivative
    //! \return The first derivative of the spline at the given x-coordinate
    inline double Deriv( const double& x )
    {
      return Evaluate( [&]( double& RetVal ) {
        return gsl_spline_eval_deriv_e( fSpline, x, fAccel, &RetVal);
      }, "Deriv" );
    }

    //! The second derivative of the spline at the given x
    //! \param x The x-coordinate at which to evaluate the derivative
    //! \return The second derivative of the spline at the given x-coordinate
    inline double Deriv2( const double& x )
    {
      return Evaluate( [&]( double& RetVal ) {
        return gsl_spline_eval_deriv2_e( fSpline, x, fAccel, &RetVal);
      }, "Deriv2" );
    }

    //! The integral over the spline between two bounds
    //! \param a The lower bound of the integral
    //! \param b The upper bound of the integral
    //! \return The integral over the spline between a and b
    inline double Integ( const double&  a, const double&  b )
    {
      if (a > b) return -Integ(b, a);
      return Evaluate( [&]( double& RetVal ) {
        return gsl_spline_eval_integ_e( fSpline, a, b, fAccel, &RetVal);
      }, "Integ" );
    }

  private:
    //! An error counter to suppress excess messages
    unsigned int nErrors;
    //! Underlying GSL machinery
    gsl_interp_accel* fAccel;
    //! Underlying GSL machinery for the spline itself
    gsl_spline* fSpline;
    //! Underlying GSL machinery for the interpolation type
    const gsl_interp_type* fInterpType;

};
