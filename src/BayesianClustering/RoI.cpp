//! \file RoI.cpp

/* ===== Cluster sources ===== */
#include "BayesianClustering/RoI.hpp"
#include "BayesianClustering/RoIproxy.hpp"
#include "BayesianClustering/Configuration.hpp"

/* ===== Local utilities ===== */
#include "Utilities/ProgressBar.hpp"
#include "Utilities/Vectorize.hpp"

// /* ===== C++ ===== */
#include <iostream>

// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

RoI::RoI( std::vector<Data>&& aData ):
  mData( std::move( aData ) ),
  mPhysicalCentreX(0), mPhysicalCentreY(0),
  mWidthX(0), mWidthY(0),
  mArea(0)
{
  std::sort( mData.begin(), mData.end() );
  std::cout << "+------------------------------------+" << std::endl;
  std::cout << "Constructed RoI with " << mData.size() << " points" << std::endl;
}

RoI::~RoI()
{
  mData.clear();
  mData.shrink_to_fit();
}

void RoI::Preprocess( const double& aMaxR, const std::vector< double >& aSigmabins2 )
{
  const double lMax2R  = 2.0 * aMaxR;
  const double lMax2R2 = lMax2R * lMax2R;

  ProgressBar lProgressBar( "Populating neighbourhood", mData.size() );
  [ & ]( const std::size_t& i ) { mData.at( i ).Preprocess( mData, i, lMax2R, lMax2R2, aSigmabins2 , lProgressBar ); } || range( mData.size() );  // Interleave threading since processing time increases with radius from origin
}

void RoI::ScanRT( const ScanConfiguration& aScanConfig, const std::function< void( RoIproxy&, const double&, const double& ) >& aCallback )
{
  auto& R = aScanConfig.Rbounds();
  Preprocess( R.max, aScanConfig.sigmabins2() );

  [&]( const std::size_t& i ) { mData.at( i ).PreprocessLocalizationScores( mData, aScanConfig, getArea() ); } || range( mData.size() );  // Interleave threading since processing time increases with radius from origin

  std::vector< RoIproxy > lRoIproxys;
  lRoIproxys.reserve( Nthreads );
  for( std::size_t i(0) ; i!=Nthreads ; ++i ) lRoIproxys.emplace_back( *this );
  ProgressBar lProgressBar( "Scan over RT", aScanConfig.Rbounds().bins * aScanConfig.Tbounds().bins );
  [&]( const std::size_t& i ) { lRoIproxys.at(i).ScanRT( aScanConfig, aCallback , lProgressBar, Nthreads, i , false ); } || range( Nthreads );
}

void RoI::Clusterize( const double& R, const double& T, const std::function< void( RoIproxy& ) >& aCallback )
{
  if( R < 0 ) throw std::runtime_error( "R must be specified and non-negative" );
  if( T < 0 ) throw std::runtime_error( "T must be specified and non-negative" );

  Preprocess( R, std::vector<double>() );

  RoIproxy lProxy( *this );
  lProxy.Clusterize( R,  T, aCallback );
}


void RoI::SetCentre( const double& aPhysicalCentreX, const double& aPhysicalCentreY )
{
  std::cout << "Centre: x=" << aPhysicalCentreX << ", y=" << aPhysicalCentreY << std::endl;
  mPhysicalCentreX = aPhysicalCentreX;
  mPhysicalCentreY = aPhysicalCentreY;
}

void RoI::SetWidth( const double& aWidthX, const double& aWidthY )
{
  std::cout << "Width: x=" << aWidthX << ", y=" << aWidthY << std::endl;
  mWidthX = aWidthX;
  mWidthY = aWidthY;
  mArea = mWidthX * mWidthY;
}

// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
