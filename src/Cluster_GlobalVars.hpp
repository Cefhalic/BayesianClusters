#pragma once

/* ===== Local utilities ===== */
#include "ListComprehension.hpp"


constexpr double nanometer = 1e-9;

constexpr long double operator"" _nanometer( long double aVal )
{
	return aVal * nanometer;
}

constexpr long double operator"" _nanometer( unsigned long long aVal )
{
	return aVal * nanometer;
}


class GlobalVars
{
public:
	void SetScale( const double& aScale )
	{
		mScale = aScale;
	}

	void SetSigmaParameters( const std::size_t& aSigmacount , const double& aSigmaMin , const double& aSigmaMax , const std::function< double( const double& ) >& aInterpolator )
	{
		mSigmacount = aSigmacount;
		mSigmaspacing = ( aSigmaMax - aSigmaMin ) / aSigmacount;
		auto lSigmabins = [ & ]( const int& i ){ return ( i * mSigmaspacing ) + aSigmaMin;  } | range( mSigmacount );

		mSigmabins = [ & ]( const double& i ){ return i * mScale; } | lSigmabins;
  		mSigmabins2= []( const double& i ){ return i * i; } | mSigmabins;
		mProbabilitySigma = aInterpolator | lSigmabins;
		mLogProbabilitySigma = []( const double& w){ return log(w); } | mProbabilitySigma;

		mSigmaspacing *= mScale;
	}

	void SetMaxR( const double& aMaxR )
	{
		mMaxR = aMaxR * mScale;
		mMaxR2 = mMaxR * mMaxR;
		mMax2R = 2.0 * mMaxR;
		mMax2R2 = mMax2R * mMax2R;
	}

	void SetBins( const std::size_t& aRbins , const std::size_t& aTbins , const double& aMinScanR = 0.0 , const double& aMaxScanR = -1  , const double& aMinScanT = 0.0 , const double& aMaxScanT = -1 )
	{
		mRbins = aRbins;
		mTbins = aTbins;

		mMinScanR = aMinScanR;
		if( aMaxScanR < 0 ) mMaxScanR = mMaxR;
		else mMaxScanR = aMaxScanR;

		mMinScanT = aMinScanT;
		if( aMaxScanT < 0 ) mMaxScanT = 2.5 * mMaxR;
		else mMaxScanT = aMaxScanT;

		mDR = ( mMaxScanR - mMinScanR ) / mRbins;
		mDT = ( mMaxScanT - mMinScanT ) / mTbins;

	}

	inline const double& scale() const { return mScale; }

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

private:	
	double mScale;

	std::size_t mSigmacount;
	double mSigmaspacing;

	std::vector< double > mSigmabins, mSigmabins2;
	std::vector< double > mProbabilitySigma, mLogProbabilitySigma;

	double mMaxR, mMaxR2, mMax2R, mMax2R2;

	double mMinScanR , mMaxScanR , mMinScanT , mMaxScanT;
	double mDR , mDT;
	std::size_t mRbins ,  mTbins;

};

GlobalVars Parameters;