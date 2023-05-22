
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

RoI::RoI( std::vector<Data>&& aData , const Configuration& aConfiguration ) : 
  mData( std::move( aData ) ),
  mConfiguration( aConfiguration )
{
  std::sort( mData.begin() , mData.end() );
  std::cout << "Constructed RoI with " << mData.size() << " points" << std::endl;  
}

void RoI::Preprocess()
{
  // Populate mNeighbour lists  
  ProgressBar2 lProgressBar( "Populating neighbourhood" , mData.size() );
  [&]( const std::size_t& i ){ mData.at( i ).Preprocess( mData , i ); } || range( mData.size() );  // Interleave threading since processing time increases with radius from origin
}

void RoI::ScanRT( const std::function< void( const RoIproxy& , const double& , const double& , std::pair<int,int>  ) >& aCallback ) 
{
  Preprocess();    

  {
    ProgressBar2 lProgressBar( "Populating localization scores" , mData.size() );
    [&]( const std::size_t& i ){ mData.at( i ).PreprocessLocalizationScores( mData ); } || range( mData.size() );  // Interleave threading since processing time increases with radius from origin
  }

  std::vector< RoIproxy > lRoIproxys;
  lRoIproxys.reserve( Nthreads );
  for( int i(0) ; i!=Nthreads ; ++i ) lRoIproxys.emplace_back( *this );
  ProgressBar2 lProgressBar( "Scan over RT"  , 0 );
  [&]( const std::size_t& i ){ lRoIproxys.at(i).ScanRT( aCallback , Nthreads , i ); } || range( Nthreads );
}

void RoI::Clusterize( const double& R , const double& T , const std::function< void( const RoIproxy& ) >& aCallback )
{
  if( R < 0 ) throw std::runtime_error( "R must be specified and non-negative" );
  if( T < 0 ) throw std::runtime_error( "T must be specified and non-negative" );

  Preprocess();    

  RoIproxy lProxy( *this );
  lProxy.Clusterize( R ,  T , aCallback );
}


// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
