//! \file Cluster.cpp

/* ===== Local utilities ===== */
#include "Utilities/GSLInterpolator.hpp"

/* ===== BOOST libraries ===== */
#include <boost/math/special_functions/erf.hpp>

/* ===== Cluster sources ===== */
#include "BayesianClustering/Cluster.hpp"
#include "BayesianClustering/Data.hpp"
#include "BayesianClustering/Configuration.hpp"


// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//! Evaluate the Gaussian normal_cdf at a given position
//! Copied from the CERN ROOT implementaion, swap ROOT::Math::erfc and ROOT::Math::erf for the boost::math version
//! \param x The position to evaluate the normal_cdf at
//! \param sigma The standard-deviation of the Gaussian
//! \param x0 The mean of the Gaussian
//! \return the value of the Gaussian normal_cdf at x
inline double normal_cdf( const double& x, const double& sigma = 1, const double& x0 = 0 )
{
  double z = ( x - x0 ) / ( sigma * sqrt(2) );
  if ( z < -1. ) return 0.5 * boost::math::erfc(-z);
  else           return 0.5 * ( 1.0 + boost::math::erf(z) );
}
// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Cluster::Parameter::Parameter() : 
A(0.0) , Bx(0.0) , By(0.0) , C(0.0) , logF(0.0),
S2(0.0)
{}
    
Cluster::Parameter& Cluster::Parameter::operator+= ( const Cluster::Parameter& aOther )
{
  A += aOther.A;
  Bx += aOther.Bx;
  By += aOther.By;
  C += aOther.C;
  logF += aOther.logF;
  return *this;
}

double Cluster::Parameter::alt_log_score() const
{
  const double pi = atan(1)*4;
  const double log2pi = log( 2*pi );

  double log_sum;
  auto inv_A( 1.0 / A ); // "A" is equivalent to nTilde
  auto sqrt_A( sqrt( A ) ); 
  double lLogMuIntegral;
  double lNubarX = Bx * inv_A; //NuBarX is Bx, nTilde is A
  double lNubarY = By * inv_A;

  lLogMuIntegral = log2pi + log(inv_A) 
                  +log( normal_cdf(sqrt_A * (1 - lNubarX)) -
                        normal_cdf(sqrt_A * (-1 - lNubarX)))
                  +log( normal_cdf(sqrt_A * (1 - lNubarY)) -
                        normal_cdf(sqrt_A * (-1 - lNubarY)));

  log_sum = logF
            -S2 / 2.0 
            +lLogMuIntegral;
  return log_sum;
}

__attribute__((flatten))
double Cluster::Parameter::log_score() const
{
  auto sqrt_A( sqrt( A ) ) , inv_A( 1.0 / A );
  auto Dx( Bx * inv_A ) , Dy( By * inv_A );
  auto E( C - ( Bx * Dx ) - ( By * Dy ) );

  double log_sum = logF - double( log( A ) ) - ( 0.5 * E );

  // We place explicit bounds checks to prRoI calls to expensive functions
  auto Gx = normal_cdf( sqrt_A * (1.0-Dx) ) - normal_cdf( sqrt_A * (-1.0-Dx) );
  if( Gx != 1.0 ) log_sum += log( Gx );
  auto Gy = normal_cdf( sqrt_A * (1.0-Dy) ) - normal_cdf( sqrt_A * (-1.0-Dy) );
  if( Gy != 1.0 ) log_sum += log( Gy );

  return log_sum;
}
// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Cluster::Cluster( const std::size_t& aParamSize ): mParams( aParamSize ),
mClusterSize( 0 ) , mLastClusterSize( 0 ) , mClusterScore( 0.0 ) , 
mParent( NULL ) ,
mData()
{}

Cluster::Cluster( const Data& aData , const std::vector< double >& aSigmabins2 ): mParams( aSigmabins2.size() ),
mClusterSize( 1 ) , mLastClusterSize( 0 ) , mClusterScore( 0.0 ) , 
mParent( NULL ) ,
mData()
{ 
  const auto s2 = aData.s * aData.s;
  auto lIt( mParams.begin() ) ;
  auto lSig2It( aSigmabins2.begin() );

  for( ; lIt != mParams.end() ; ++lIt , ++lSig2It )
  {
    double w = 1.0 / ( s2 + *lSig2It );
    lIt->A = w;
    lIt->Bx = (w * aData.x);
    lIt->By = (w * aData.y);
    lIt->C = (w * aData.r2);
    lIt->logF = PRECISION( log( w ) );
  }
}

void Cluster::UpdateLogScore( const ScanConfiguration& aScanConfig )
{
  const auto& lSigmaBins           = aScanConfig.sigmabins();
  const auto& lLogProbabilitySigma = aScanConfig.log_probability_sigma();

  static constexpr double pi = atan(1)*4;
  static constexpr double log2pi = log( 2*pi );

  if( mClusterSize <= mLastClusterSize ) return; // We were not bigger than the previous size when we were evaluated - score is still valid
  mLastClusterSize = mClusterSize;

  const std::size_t lSigmaCount( lSigmaBins.size() );
  thread_local static std::vector< double > MuIntegral( lSigmaCount , 1.0 );
  thread_local static std::vector< double > integralArguments( lSigmaCount , 1.0 );
  double lMaxValue(-9E99);
  for( std::size_t i(0) ; i!=lSigmaCount ; ++i ) {
    double lValue =  mParams[i].log_score() + lLogProbabilitySigma[i];
    if (lValue > lMaxValue) lMaxValue = lValue;
    integralArguments[i] = lValue;
  }
  //pass again to set the MuIntegral Correctly
  for( std::size_t i(0) ; i!=lSigmaCount ; ++i ) MuIntegral[i] = exp(integralArguments[i] - lMaxValue);

  thread_local static GSLInterpolator lInt( gsl_interp_linear , lSigmaCount );
  lInt.SetData( lSigmaBins , MuIntegral );

  const double Lower( lSigmaBins[0] ) , Upper( lSigmaBins[lSigmaCount-1] );
  mClusterScore = double( log( lInt.Integ( Lower , Upper ) ) ) + lMaxValue - double( log( 4.0 ) ) + (log2pi * (1.0-mClusterSize));  
  mClusterScore += log(0.25) -(mClusterSize * log2pi);
}

Cluster& Cluster::operator+= ( const Cluster& aOther )
{
  auto lIt( mParams.begin() );
  auto lIt2( aOther.mParams.begin() );

  for( ; lIt != mParams.end() ; ++lIt , ++lIt2 ) *lIt += *lIt2;
  mClusterSize += aOther.mClusterSize;
  return *this;
}

Cluster* Cluster::GetParent()
{
  if( mParent ) return mParent = mParent->GetParent();
  return this;
}
// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

