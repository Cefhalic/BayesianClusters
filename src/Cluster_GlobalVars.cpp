

/* ===== Local utilities ===== */
#include "ListComprehension.hpp"
#include "Cluster_GlobalVars.hpp"

/* ===== For Root ===== */
#include "Math/SpecFunc.h" 

GlobalVars::GlobalVars() :
	mScale(-1) , mScale2(-1),
	mSigmacount(-1), mSigmaspacing(-1),
	mMaxR(-1), mMaxR2(-1), mMax2R(-1), mMax2R2(-1),
	mMinScanR(-1), mMaxScanR(-1), mMinScanT(-1), mMaxScanT(-1),
	mDR(-1), mDT(-1),
	mRbins(-1),  mTbins(-1),
	mLogPb(-1), mLogPbDagger(-1), 
	mAlpha(-1), mLogAlpha(-1), mLogGammaAlpha(-1),
	mValidate(false)
{}

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
	if( mScale < 0 ) throw std::runtime_error( "Scale must be set before setting Sigma parameters" );

	mSigmacount = aSigmacount;
	mSigmaspacing = ( aSigmaMax - aSigmaMin ) / aSigmacount;
	auto lSigmabins = [ & ]( const int& i ){ return ( i * mSigmaspacing ) + aSigmaMin;  } | range( mSigmacount );

	mSigmabins = [ & ]( const double& i ){ return toAlgorithmUnits( i ); } | lSigmabins;
	mSigmabins2 = []( const double& i ){ return i * i; } | mSigmabins;
	mProbabilitySigma = aInterpolator | lSigmabins;
	mLogProbabilitySigma = []( const double& w){ return log(w); } | mProbabilitySigma;

	mSigmaspacing *= mScale;
}

void GlobalVars::SetMaxR( const double& aMaxR )
{
	if( mScale < 0 ) throw std::runtime_error( "Scale must be set before setting Sigma parameters" );

	mMaxR = toAlgorithmUnits( aMaxR );
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
	mLogPb = log( aPB );
	mLogPbDagger = log( 1-aPB );
	mAlpha = aAlpha;
	mLogAlpha = log( aAlpha );
	mLogGammaAlpha = ROOT::Math::lgamma( aAlpha );
}

void GlobalVars::SetValidate( const bool& aValidate )
{
	mValidate = aValidate;
}

