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

	void SetCentre( const double& aPhysicalCentreX , const double& aPhysicalCentreY );
	void SetZoom( const double& aScale );
	void SetSigmaParameters( const std::size_t& aSigmacount , const double& aSigmaMin , const double& aSigmaMax , const std::function< double( const double& ) >& aInterpolator );
	void SetMaxR( const double& aMaxR );
	void SetBins( const std::size_t& aRbins , const std::size_t& aTbins , const double& aMinScanR = 0.0 , const double& aMaxScanR = -1  , const double& aMinScanT = 0.0 , const double& aMaxScanT = -1 );
	void SetPbAlpha( const double& aPB , const double& aAlpha );
	void SetValidate( const bool& aValidate );

	std::string FromCommandline( int argc , char **argv );

public:
	inline const double& scale2() const { return mScale2; }

	inline const std::size_t& sigmacount() const { return mSigmacount; }
	inline const double& sigmaspacing() const { return mSigmaspacing; }

	inline const std::vector< double >& sigmabins( ) const { return mSigmabins; }
	inline const std::vector< double >& sigmabins2( ) const { return mSigmabins2; }

	inline const std::vector< double >& probability_sigma( ) const { return mProbabilitySigma; }
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
	double mScale;
  double mScale2;
  double mPhysicalCentreX;
  double mPhysicalCentreY;

	std::size_t mSigmacount;
	double mSigmaspacing;

	std::vector< double > mSigmabins;
  std::vector< double > mSigmabins2;
	std::vector< double > mProbabilitySigma;
  std::vector< double > mLogProbabilitySigma;

	double mMaxR;
  double mMaxR2;
  double mMax2R;
  double mMax2R2;

	double mMinScanR;
  double mMaxScanR;
  double mMinScanT;
  double mMaxScanT;
	double mDR;
  double mDT;
	std::size_t mRbins;
	std::size_t mTbins;

	double mAlpha;
  double mLogAlpha;
  double mLogGammaAlpha;
  double mLogPb;
  double mLogPbDagger;

  //! Run the validation on the clustering 
	bool mValidate;
};
