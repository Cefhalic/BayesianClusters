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


constexpr double nanometer  = 1e-9;
constexpr double micrometer = 1e-6;
constexpr double millimeter = 1e-3;
constexpr double meter      = 1e-0;
const std::map< std::string , double > UnitMap{ {"nm",nanometer} , {"um",micrometer} , {"mm",millimeter} , {"m",meter} };


constexpr long double operator"" _nanometer( long double aVal )
{
	return aVal * nanometer;
}

constexpr long double operator"" _nanometer( unsigned long long aVal )
{
	return aVal * nanometer;
}

constexpr long double operator"" _micrometer( long double aVal )
{
	return aVal * micrometer;
}

constexpr long double operator"" _micrometer( unsigned long long aVal )
{
	return aVal * micrometer;
}

inline long double StrToDist( const std::string& aStr )
{
	std::stringstream lStr;
	lStr << aStr;
	double lVal;
	std::string lUnits;
	lStr >> lVal >> lUnits;
	return lVal * UnitMap.at( lUnits );
}





class GlobalVars
{
public:
	GlobalVars();

	void SetCentre( const double& aPhysicalCentreX , const double& aPhysicalCentreY );
	void SetZoom( const double& aScale );
	void SetSigmaParameters( const std::size_t& aSigmacount , const double& aSigmaMin , const double& aSigmaMax , const std::function< double( const double& ) >& aInterpolator );
	void SetRBins( const std::size_t& aRbins , const double& aMinScanR , const double& aMaxScanR );
	void SetTBins( const std::size_t& aTbins , const double& aMinScanT , const double& aMaxScanT );
	void SetPb( const double& aPB );
	void SetAlpha( const double& aAlpha );
	void SetValidate( const bool& aValidate );
	void SetNormalization( const double& aNormalization );

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

	inline const double& normalization() const { return mNormalization; }

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
	double mScale , mScale2;
  	double mPhysicalCentreX , mPhysicalCentreY;

	std::size_t mSigmacount;
	double mSigmaspacing;

	std::vector< double > mSigmabins, mSigmabins2;
	std::vector< double > mProbabilitySigma, mLogProbabilitySigma;

	double mNormalization;

	double mMaxR, mMaxR2, mMax2R, mMax2R2;

	double mMinScanR , mMaxScanR , mMinScanT , mMaxScanT;
	double mDR , mDT;
	std::size_t mRbins ,  mTbins;

	double mAlpha , mLogAlpha , mLogGammaAlpha , mLogPb , mLogPbDagger;

	bool mValidate;
};
