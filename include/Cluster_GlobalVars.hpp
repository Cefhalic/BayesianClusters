#pragma once

/* ===== C++ ===== */
#include <functional>
#include <map>
#include <vector>
#include <sstream>
#include <string>

/* ===== BOOST libraries ===== */
#include "boost/program_options.hpp"
namespace po = boost::program_options;


//! Define a constant for converting nanometers to meters
constexpr double nanometer  = 1e-9;
//! Define a constant for converting micrometers to meters
constexpr double micrometer = 1e-6;
//! Define a constant for converting millimeters to meters
constexpr double millimeter = 1e-3;
//! Define a constant for converting meters to meters
constexpr double meter      = 1e-0;

//! A map for converting string representations of SI units to scaling factors
const std::map< std::string , double > UnitMap{ {"nm",nanometer} , {"um",micrometer} , {"mm",millimeter} , {"m",meter} };

//! User-defined literals fot nanometer quantities
//! \param aVal The specified value
//! \return The literal value
constexpr long double operator"" _nanometer( long double aVal )
{
	return aVal * nanometer;
}

//! User-defined literals fot nanometer quantities
//! \param aVal The specified value
//! \return The literal value
constexpr long double operator"" _nanometer( unsigned long long aVal )
{
	return aVal * nanometer;
}

//! User-defined literals fot micrometer quantities
//! \param aVal The specified value
//! \return The literal value
constexpr long double operator"" _micrometer( long double aVal )
{
	return aVal * micrometer;
}

//! User-defined literals fot micrometer quantities
//! \param aVal The specified value
//! \return The literal value
constexpr long double operator"" _micrometer( unsigned long long aVal )
{
	return aVal * micrometer;
}

//! Convert a string representation to a distance
//! \param aStr A string representation of a distance
//! \return The literal value
inline long double StrToDist( const std::string& aStr )
{
	std::stringstream lStr;
	lStr << aStr;
	double lVal;
	std::string lUnits;
	lStr >> lVal >> lUnits;
	return lVal * UnitMap.at( lUnits );
}




//! Class for storing the configuration parameters
class GlobalVars
{
public:
  //! Default constructor
	GlobalVars();

  //! Setter for the centre of the scan window
  //! \param aPhysicalCentreX The x-coordinate of the centre of the window in physical units (becomes 0 in algorithm units)
  //! \param aPhysicalCentreY The y-coordinate of the centre of the window in physical units (becomes 0 in algorithm units)
	void SetCentre( const double& aPhysicalCentreX , const double& aPhysicalCentreY );
  
  //! Setter for the half-width of the scan window
  //! \param aScale The scale of the window in physical units (becomes ±1 in algorithm units)
	void SetZoom( const double& aScale );

  //! Setter for the sigma-bins to be integrated over
  //! \param aSigmacount   The number of sigma bins
  //! \param aSigmaMin     The lowest sigma bin
  //! \param aSigmaMax     The highest sigma bin
  //! \param aInterpolator Function-object to generate the probability of any given sigma
	void SetSigmaParameters( const std::size_t& aSigmacount , const double& aSigmaMin , const double& aSigmaMax , const std::function< double( const double& ) >& aInterpolator );

  //! Setter for the maximum allowed R parameter
  //! \param aMaxR The maximum allowed R parameter
	void SetMaxR( const double& aMaxR );

  //! Setter for the R and T bins for the RT scan
  //! \param aRbins    The number of R bins to scan over
  //! \param aTbins    The number of T bins to scan over
  //! \param aMinScanR The lowest value of R to scan
  //! \param aMaxScanR The largest value of R to scan
  //! \param aMinScanT The lowest value of T to scan
  //! \param aMaxScanT The largest value of T to scan
	void SetBins( const std::size_t& aRbins , const std::size_t& aTbins , const double& aMinScanR = 0.0 , const double& aMaxScanR = -1  , const double& aMinScanT = 0.0 , const double& aMaxScanT = -1 );

  //! Setter for the P_b and alpha parameters
  //! \param aPB    The P_b parameter
  //! \param aAlpha The alpha parameter
	void SetPbAlpha( const double& aPB , const double& aAlpha );
  
  //! Set whether to validate clusterization
  //! \param aValidate Whether to validate clusterization
	void SetValidate( const bool& aValidate );

  //! Parse the parameters when passed in as commandline arguments
  //! \param argc The number of commandline arguments
  //! \param argv The commandline arguments
  //! \return The name of the event file
	std::string FromCommandline( int argc , char **argv );

public:
  //! Getter for the scale-parameter squared
  //! \return The scale-parameter squared
	inline const double& scale2() const { return mScale2; }

  //! Getter for the sigma count
  //! \return The sigma count
	inline const std::size_t& sigmacount() const { return mSigmacount; }

  //! Getter for the sigma spacing
  //! \return The sigma spacing  
	inline const double& sigmaspacing() const { return mSigmaspacing; }

  //! Getter for the values of sigma
  //! \return The values of sigma
	inline const std::vector< double >& sigmabins( ) const { return mSigmabins; }
  //! Getter for the values of sigma squared
  //! \return The values of sigma squared 
	inline const std::vector< double >& sigmabins2( ) const { return mSigmabins2; }
  //! Getter for the probabilities of a given sigma
  //! \return The probabilities of given sigma
	inline const std::vector< double >& probability_sigma( ) const { return mProbabilitySigma; }
  //! Getter for the log of the probabilities of a given sigma
  //! \return The log of the probabilities of given sigma
	inline const std::vector< double >& log_probability_sigma( ) const { return mLogProbabilitySigma; }

	inline const double& sigmabins( const std::size_t& i ) const { return mSigmabins[i]; }
	inline const double& sigmabins2( const std::size_t& i ) const { return mSigmabins2[i]; }

	inline const double& probability_sigma( const std::size_t& i ) const { return mProbabilitySigma[i]; }
	inline const double& log_probability_sigma( const std::size_t& i ) const { return mLogProbabilitySigma[i]; }

	inline const double& maxR() const { return mMaxR; }
	inline const double& maxR2() const { return mMaxR2; }
	inline const double& max2R() const { return mMax2R; }
	inline const double& max2R2() const { return mMax2R2; }

	inline const double& minScanR() const { return mMinScanR; }
	inline const double& maxScanR() const { return mMaxScanR; }
	inline const double& minScanT() const { return mMinScanT; }
	inline const double& maxScanT() const { return mMaxScanT; }

	inline const double& dR() const { return mDR; }
	inline const std::size_t& Rbins() const { return mRbins; }
	inline const double& dT() const { return mDT; }
	inline const std::size_t& Tbins() const { return mTbins; }

	inline const double& logPb() const { return mLogPb; }
	inline const double& logPbDagger() const { return mLogPbDagger; }
	inline const double& alpha() const { return mAlpha; }
	inline const double& logAlpha() const { return mLogAlpha; }
	inline const double& logGammaAlpha() const { return mLogGammaAlpha; }

	inline const bool& validate() const { return mValidate; }

	inline double toPhysicalUnits( const double& aAlgorithmUnits ) const 
	{
		return aAlgorithmUnits / mScale;
	}

	inline double toAlgorithmUnits( const double& aPhysicalUnits ) const
	{
	  return aPhysicalUnits * mScale;
	}


	inline double toPhysicalX( const double& aAlgorithmX ) const  
	{
		return toPhysicalUnits( aAlgorithmX ) + mPhysicalCentreX;
	}

	inline double toAlgorithmX( const double& aPhysicalX ) const
	{
		return toAlgorithmUnits( aPhysicalX - mPhysicalCentreX );
	}

	inline double toPhysicalY( const double& aAlgorithmY ) const 
	{
		return toPhysicalUnits( aAlgorithmY ) + mPhysicalCentreY;
	}

	inline double toAlgorithmY( const double& aPhysicalY ) const
	{
		return toAlgorithmUnits( aPhysicalY - mPhysicalCentreY );
	}

private:
  //! The scale parameter
	double mScale;
  
  //! The scale parameter squared
  double mScale2;

  //! The x-coordinate of the centre of the window in physical units
  double mPhysicalCentreX;
  //! The y-coordinate of the centre of the window in physical units 
  double mPhysicalCentreY;

  //! The number of sigma bins
	std::size_t mSigmacount;
  //! The spacing of sigma bins
	double mSigmaspacing;
  //! The values of sigma
	std::vector< double > mSigmabins;
  //! The values of sigma squared
  std::vector< double > mSigmabins2;
  //! The probability of a given sigma
	std::vector< double > mProbabilitySigma;
  //! The log-probability of a gievn sigma
  std::vector< double > mLogProbabilitySigma;

  //! The maximum value of R
	double mMaxR;
  //! The maximum value of R squared
  double mMaxR2;
  //! The maximum value of 2R  
  double mMax2R;
  //! The maximum value of 2R squared
  double mMax2R2;

  //! The lowest value of R to scan 
	double mMinScanR;
  //! The largest value of R to scan
  double mMaxScanR;
  //! The lowest value of T to scan
  double mMinScanT;
  //! The largest value of T to scan
  double mMaxScanT;
  //! The spacing of value of R to scan 
	double mDR;
  //! The spacing of value of T to scan 
  double mDT;
  //! The number of R values to scan 
	std::size_t mRbins;
  //! The number of T values to scan 
	std::size_t mTbins;

  //! The alpha parameter
	double mAlpha;
  //! Logarithm of the alpha parameter
  double mLogAlpha;
  //! Logarithm of the gamma function of alpha parameter  
  double mLogGammaAlpha;
  
  //! Logarithm of the P_b parameter  
  double mLogPb;
  //! Logarithm of the( 1- P_b ) parameter  
  double mLogPbDagger;

  //! Run the validation on the clustering 
	bool mValidate;
};
