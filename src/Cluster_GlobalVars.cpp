

/* ===== Local utilities ===== */
#include "ListComprehension.hpp"
#include "Cluster_GlobalVars.hpp"


void GlobalVars::SetZoom( const double& aScale )
{
	mScale = 2.0 / aScale;
	mScale2 = mScale * mScale;
}

// void GlobalVars::SetScale( const double& aScale )
// {
// 	mScale = aScale;
// 	mScale2 = mScale * mScale;
// }

void GlobalVars::SetSigmaParameters( const std::size_t& aSigmacount , const double& aSigmaMin , const double& aSigmaMax , const std::function< double( const double& ) >& aInterpolator )
{
	mSigmacount = aSigmacount;
	mSigmaspacing = ( aSigmaMax - aSigmaMin ) / aSigmacount;
	auto lSigmabins = [ & ]( const int& i ){ return ( i * mSigmaspacing ) + aSigmaMin;  } | range( mSigmacount );

	mSigmabins = [ & ]( const double& i ){ return i * mScale; } | lSigmabins;
	mSigmabins2 = []( const double& i ){ return i * i; } | mSigmabins;
	mProbabilitySigma = aInterpolator | lSigmabins;
	mLogProbabilitySigma = []( const double& w){ return log(w); } | mProbabilitySigma;

	mSigmaspacing *= mScale;
}

void GlobalVars::SetMaxR( const double& aMaxR )
{
	mMaxR = aMaxR * mScale;
	mMaxR2 = mMaxR * mMaxR;
	mMax2R = 2.0 * mMaxR;
	mMax2R2 = mMax2R * mMax2R;
}

void GlobalVars::SetBins( const std::size_t& aRbins , const std::size_t& aTbins , const double& aMinScanR , const double& aMaxScanR  , const double& aMinScanT , const double& aMaxScanT )
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

void GlobalVars::SetPbAlpha( const double& aPB , const double& aAlpha )
{
	mPB = aPB;
	mAlpha = aAlpha;
}

GlobalVars Parameters;
