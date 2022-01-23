#pragma once

/* ===== C++ ===== */
#include <functional>
#include <vector>

constexpr double nanometer = 1e-9;
constexpr double micrometer = 1e-6;

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



class GlobalVars
{
public:
	GlobalVars();

	void SetZoom( const double& aScale );
	// void SetScale( const double& aScale );
	void SetSigmaParameters( const std::size_t& aSigmacount , const double& aSigmaMin , const double& aSigmaMax , const std::function< double( const double& ) >& aInterpolator );
	void SetMaxR( const double& aMaxR );
	void SetBins( const std::size_t& aRbins , const std::size_t& aTbins , const double& aMinScanR = 0.0 , const double& aMaxScanR = -1  , const double& aMinScanT = 0.0 , const double& aMaxScanT = -1 );
	void SetPbAlpha( const double& aPB , const double& aAlpha );
	void SetValidate( const bool& aValidate );

public:
	//inline const double& scale() const { return mScale; }
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

private:	
	double mScale , mScale2;

	std::size_t mSigmacount;
	double mSigmaspacing;

	std::vector< double > mSigmabins, mSigmabins2;
	std::vector< double > mProbabilitySigma, mLogProbabilitySigma;

	double mMaxR, mMaxR2, mMax2R, mMax2R2;

	double mMinScanR , mMaxScanR , mMinScanT , mMaxScanT;
	double mDR , mDT;
	std::size_t mRbins ,  mTbins;

	double mAlpha , mLogAlpha , mLogGammaAlpha , mLogPb , mLogPbDagger;

	bool mValidate;
};
