/* ===== For Root ===== */
#include "Math/ProbFunc.h" 
#include "Math/Interpolator.h" 

/* ===== Cluster sources ===== */
#include "BayesianClustering/Cluster.hpp"
#include "BayesianClustering/Data.hpp"
#include "BayesianClustering/Configuration.hpp"


// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Cluster::Parameter::Parameter() : 
A(0.0) , Bx(0.0) , By(0.0) , C(0.0) , logF(0.0)
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

inline double CDF( const double& aArg )
{
  // Above 8 or below -8 are indistinguishable from 0 and 1 respectively
  // if( aArg > 8.0 ) return 1.0;
  // if( aArg < -8.0 ) return 0.0;
  return  ROOT::Math::normal_cdf( aArg );
}

__attribute__((flatten))
double Cluster::Parameter::log_score() const
{
  auto sqrt_A( sqrt( A ) ) , inv_A( 1.0 / A );
  auto Dx( Bx * inv_A ) , Dy( By * inv_A );
  auto E( C - ( Bx * Dx ) - ( By * Dy ) );

  double log_sum = logF - double( log( A ) ) - ( 0.5 * E );

  // We place explicit bounds checks to prevent calls to expensive functions
  auto Gx = CDF( sqrt_A * (1.0-Dx) ) - CDF( sqrt_A * (-1.0-Dx) );
  if( Gx != 1.0 ) log_sum += log( Gx );
  auto Gy = CDF( sqrt_A * (1.0-Dy) ) - CDF( sqrt_A * (-1.0-Dy) );
  if( Gy != 1.0 ) log_sum += log( Gy );

  return log_sum;
}
// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Cluster::Cluster(): mParams( Configuration::Instance.sigmacount() ),
mClusterSize( 0 ) , mLastClusterSize( 0 ) , mClusterScore( 0.0 ) , 
mParent( NULL )
{}

Cluster::Cluster( const Data& aData ): mParams( Configuration::Instance.sigmacount() ),
mClusterSize( 1 ) , mLastClusterSize( 0 ) , mClusterScore( 0.0 ) , 
mParent( NULL )
{ 
  const auto s2 = aData.s * aData.s;
  auto lIt( mParams.begin() ) ;
  auto lSig2It( Configuration::Instance.sigmabins2().begin() );

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

void Cluster::UpdateLogScore()
{
  static constexpr double pi = atan(1)*4;
  static constexpr double log2pi = log( 2*pi );

  if( mClusterSize <= mLastClusterSize ) return; // We were not bigger than the previous size when we were evaluated - score is still valid
  mLastClusterSize = mClusterSize;

  thread_local static std::vector< double > MuIntegral( Configuration::Instance.sigmacount() , 1.0 );
  thread_local static std::vector< double > integralArguments( Configuration::Instance.sigmacount() , 1.0 );
  double largestArg(-9E99);
  // double constant( mParams[0].log_score() + Configuration::Instance.log_probability_sigma( 0 ) );
  // for( std::size_t i(1) ; i!=Configuration::Instance.sigmacount() ; ++i ) MuIntegral[i] = exp( mParams[i].log_score() + Configuration::Instance.log_probability_sigma( i ) - constant );
  double tempArg;
  for( std::size_t i(0) ; i!=Configuration::Instance.sigmacount() ; ++i ) {
    tempArg =  mParams[i].log_score() + Configuration::Instance.log_probability_sigma( i );
    if (tempArg > largestArg) largestArg = tempArg;

    integralArguments[i] = tempArg;
  }
  //pass again to set the MuIntegral Correctly
  for( std::size_t i(0) ; i!=Configuration::Instance.sigmacount() ; ++i ) MuIntegral[i] = exp(integralArguments[i] - largestArg);

  thread_local static ROOT::Math::Interpolator lInt( Configuration::Instance.sigmacount() , ROOT::Math::Interpolation::kLINEAR );
  lInt.SetData( Configuration::Instance.sigmabins() , MuIntegral );

  static const double Lower( Configuration::Instance.sigmabins(0) ) , Upper( Configuration::Instance.sigmabins(Configuration::Instance.sigmacount()-1) );
  // mClusterScore = double( log( lInt.Integ( Lower , Upper ) ) ) + constant - double( log( 4.0 ) ) + (log2pi * (1.0-mClusterSize));  
  mClusterScore = double( log( lInt.Integ( Lower , Upper ) ) ) + largestArg - double( log( 4.0 ) ) + (log2pi * (1.0-mClusterSize));  
  mClusterScore += log(0.25) -(mClusterSize * log2pi);
}

Cluster& Cluster::operator+= ( const Cluster& aOther )
{
  // if( &aOther == this ) throw std::runtime_error( "Error #1" );
  // if( aOther.mClusterSize == 0 ) throw std::runtime_error( "Error #2" );
  // if( aOther.mParent ) throw std::runtime_error( "Error #3" );

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

