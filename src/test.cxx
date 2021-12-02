/* ===== Cluster sources ===== */
#include "Cluster_Data.hpp"
#include "Cluster_DataSources.hpp"
#include "Cluster_PlotTools.hpp"

// /* ===== C++ ===== */
#include <thread>
#include <iostream>
#include <math.h>
#include <vector>
#include <map>

/* ===== For Root ===== */
// #include "TH2D.h"
#include "Math/ProbFunc.h" 

/* ===== Local utilities ===== */
#include "ListComprehension.hpp"
#include "ProgressBar.hpp"
#include "RootWindow.hpp"
#include "Vectorize.hpp"


// // Burmann approximation of the Gaussian cumulative-distribution function
// double Phi( const double& x )
// {
//   auto y = exp( -1.0*x*x );

//   constexpr double pi = atan(1)*4;
//   constexpr double rt_pi = sqrt( pi );
//   constexpr double inv_rt_pi = 1.0 / rt_pi;

//   double lRet = ( rt_pi/2.0 ) + ( y * 31.0/200.0 ) - ( y*y * 341.0/8000.0 );
//   lRet       *= inv_rt_pi * ( (x > 0) - (x < 0) ) * sqrt( 1-y );

//   return 0.5 + lRet;
// }


#define PRINT( x ) std::cout << ( #x ) << " = " << ( x ) << std::endl;

double Score( std::pair< const Data* const , std::vector< Data* > >& aCluster )
{
  if( aCluster.first == NULL ) return 0.0;

  constexpr double pi = atan(1)*4;
  constexpr double logA( -1.0 * log( 4.0 ) );
  const double n( aCluster.second.size() );
  const double pi_term( log( 2*pi ) * (1.0-n) ); // IF THE 1-n is integers, this gives the wrong answer!!!! WHY????????

  auto IntegrandExpr = [ &n , &aCluster , &logA , &pi_term ]( const uint32_t& i )
  {
    double n_tilde( 0.0 ) , sum_log_w( 0.0 );
    double nu_bar_x( 0.0 ) , nu_bar_y( 0.0 );
    for( auto& j : aCluster.second )
    {
      auto& w = j->w[ i ];
      n_tilde  += w;
      sum_log_w += j->log_w[ i ];      
      nu_bar_x += ( w * j->x );
      nu_bar_y += ( w * j->y );
    }

    auto sqrt_n_tilde( sqrt( n_tilde ) ) , inv_n_tilde( 1.0 / n_tilde );

    nu_bar_x *= inv_n_tilde;
    nu_bar_y *= inv_n_tilde;

    double S2( 0.0 );
    for( auto& j : aCluster.second )
    { 
      auto& w = j->w[ i ];
      auto x( j->x - nu_bar_x ) , y( j->y - nu_bar_y );
      auto t( (x*x) + (y*y) );
      S2 += w*t;
    }

    double phi_x( ROOT::Math::normal_cdf( sqrt_n_tilde * (1.0-nu_bar_x) ) - ROOT::Math::normal_cdf( sqrt_n_tilde * (-1.0-nu_bar_x) ) ) , phi_y( ROOT::Math::normal_cdf( sqrt_n_tilde * (1.0-nu_bar_y) ) - ROOT::Math::normal_cdf( sqrt_n_tilde * (-1.0-nu_bar_y) ) );
    auto sum = logA + pi_term - log( n_tilde ) + sum_log_w - ( 0.5 * S2 ) + log( phi_x ) + log( phi_y ) + Parameters.log_probability_sigma( i );

    // std::cout << "==================================" << std::endl;
    // PRINT( logA );
    // PRINT( pi_term );
    // PRINT( - log( n_tilde ) );
    // PRINT( sum_log_w );
    // PRINT( - ( 0.5 * S2 ) );
    // PRINT( log( phi_x ) );
    // PRINT( log( phi_y ) );
    // PRINT( log_p_sigma[ i ] );
    // std::cout << "---------------" << std::endl;
    // PRINT( sum );

    return sum;
  };    

  auto Integrands = IntegrandExpr | range( Parameters.sigmacount() );
  auto& Max = *std::max_element( Integrands.begin() , Integrands.end() );

  double integral( 0.0 );
  for( auto& i : Integrands ) integral += exp( i - Max );
  integral *= Parameters.sigmaspacing();  

  return log( integral ) + Max; 

}







double __Cluster__( std::vector<Data>& aData , const double& R , const double& T , std::map< const Data* , std::vector< Data* > >& aClusters )
{
  // -----------------------------------------------------------------------------------------
  // Unclear why this is quicker than maintaining clusters from a scan with a higher threshold
  for( auto& i : aData ) i.parent = NULL;
  aClusters.clear();
  // -----------------------------------------------------------------------------------------

  const double twoR2 = 4.0 * R * R;
  for( auto& i : aData ) i.Clusterize( twoR2 , T ); 

  for( auto& i : aData ) aClusters[ i.GetParent() ].push_back( &i ); // <-- Bodge - could be done internally analagous to parent calculation

  // double lScore ( 0.0 );
  // for( auto& k : aClusters ) lScore += Score( k );

  using type = decltype( *aClusters.begin() );
  auto CreateLambda2 = []( type& aData ){ return [ &aData ](){ return Score( aData ); }; }; // Create a lambda which returns a lambda. Then use list comprehension to create
  auto lRet = Vectorize( CreateLambda2 | aClusters );                                       // a vector of lambdas and pass to the threadpool. Nerd level cranked to 11 
  auto lScore = std::accumulate( lRet.begin() , lRet.end() , 0.0 );

  //std::cout << lScore << std::endl;


  return lScore;
}



void Cluster( std::vector<Data>& aData , const double& R , const double& T , std::map< const Data* , std::vector< Data* > >& aClusters )
{
  const std::size_t Nminus1( aData.size() - 1 );
  const auto R2 = R * R;

  // Reset ahead of UpdateLocalization and Clusterize
  for( auto& i : aData ){
    // i.parent = NULL;    // Uncomment this if resetting parent is removed from __Cluster__
    i.neighbourit = i.neighbours.begin();
  }

  //for( auto &i : aData ) i.UpdateLocalization( R2 , Nminus1 );
  auto CreateLambda = [ &Nminus1 , &R2 ]( Data& aData ){ return [ &aData,  &Nminus1 , &R2 ](){ aData.UpdateLocalization( R2 , Nminus1 ); }; }; // Create a lambda which returns a lambda. Then use list comprehension to create
  Vectorize( CreateLambda | aData );                                                                                                           // a vector of lambdas and pass to the threadpool. Nerd level cranked to 11 

  __Cluster__( aData , R , T , aClusters );
}



TH2D *Nclust, *ClustSize, *ClustScore;


void ScanRT( std::vector<Data>& aData )
{
  const std::size_t Nminus1( aData.size() - 1 );

  double R( Parameters.minScanR() ) , R2( 0 ) , T( 0 );

  ProgressBar lProgressBar( "Scan over RT" , Parameters.Rbins() * Parameters.Tbins() );

  for( uint32_t i(0) ; i!=Parameters.Rbins() ; ++i , R+=Parameters.dR() )
  {
    R2 = R * R;
    T = Parameters.maxScanT();

    //for( auto &i : aData ) i.UpdateLocalization( R2 , Nminus1 );
    auto CreateLambda = [ &Nminus1 , &R2 ]( Data& aData ){ return [ &aData,  &Nminus1 , &R2 ](){ aData.UpdateLocalization( R2 , Nminus1 ); }; }; // Create a lambda which returns a lambda. Then use list comprehension to create
    Vectorize( CreateLambda | aData );                                                                                                           // a vector of lambdas and pass to the threadpool. Nerd level cranked to 11 

    // Thought this would be quicker than resetting it in each call to __Cluster__, apparently not
    // for( auto& i : aData ) i.parent = NULL;               
    std::map< const Data* , std::vector< Data* > > lClusters;

    for( uint32_t j(0) ; j!=Parameters.Tbins() ; ++j , T-=Parameters.dT() , ++lProgressBar )
    {
      auto lScore = __Cluster__( aData , R , T , lClusters );

      ClustScore -> Fill( R , T , lScore );

      Nclust -> Fill( R , T , lClusters.size() );
      double mean( 0.0 );
      for ( auto& i: lClusters ) if( i.first ) mean += i.second.size();
      ClustSize -> Fill( R , T , mean / lClusters.size() );

    }

  }

}




void PrepData( std::vector<Data>& aData )
{
  ProgressBar2 lProgressBar( "Preprocessing" , aData.size() );

  // Crucial step is to sort!
  std::sort( aData.begin() , aData.end() );

  // Now populate the neighbour lists
  auto CreateLambda = [ &aData ]( const uint32_t& aIndex ){ return [ &aIndex , &aData ](){
    aData.at( aIndex ).PopulateNeighbours( aData.begin() + aIndex + 1 , aData.end() , aData.rbegin() + aData.size() - aIndex , aData.rend() );
  }; };                                              // Create a lambda which returns a lambda. Then use list comprehension to create
  Vectorize( CreateLambda | range( aData.size() ) ); // a vector of lambdas and pass to the threadpool. Nerd level cranked to 11 

  // ... any other preparation ...
}



/* ===== Main function ===== */
int main(int argc, char **argv)
{
  Parameters.SetSigmaCountAndSpacing( 10 , 1e-3 );
  Parameters.SetProbabilitySigma( { 0.03631079, 0.110302441, 0.214839819, 0.268302465, 0.214839819, 0.110302441, 0.03631079, 0.007664194, 0.001037236, 9.00054E-05 } );
  Parameters.SetMaxR( 0.02 );
  Parameters.SetBins( 40 , 40 );

//  auto lData = LoadCSV( "1_un_red.csv" , 1./64000. , -1. , -1. ); // Full file
  // auto lData = LoadCSV( "1_un_red.csv" , 1./10000. , 87000. , 32000. ); // One cluster
//  auto lData = LoadCSV( "1_un_red.csv" , 1./1000. , 87000. , 32000. ); // Very zoomed


  auto lData = CreatePseudoData( 10007 , 500 , 500 , .01 );
  //auto lData = CreatePseudoData( 70000 , 500 , 500 , .01 );

  //InteractiveDisplay( [ lData ](){ DrawPoints( lData ); } );

  PrepData( lData );


  auto Rlo = Parameters.minScanR() - ( 0.5 * Parameters.dR() );
  auto Rhi = Parameters.maxScanR() - ( 0.5 * Parameters.dR() );
  auto Tlo = Parameters.minScanT() - ( 0.5 * Parameters.dT() );
  auto Thi = Parameters.maxScanT() - ( 0.5 * Parameters.dT() );

  Nclust     = new TH2D( "Nclust" , "N_{clusters};r;T" , Parameters.Rbins() , Rlo , Rhi , Parameters.Tbins() , Tlo , Thi );
  ClustSize  = new TH2D( "ClustSize" , "<N_{points}>;r;T" , Parameters.Rbins() , Rlo , Rhi , Parameters.Tbins() , Tlo , Thi );
  ClustScore = new TH2D( "ClustScore" , "Score;r;T" , Parameters.Rbins() , Rlo , Rhi , Parameters.Tbins() , Tlo , Thi );


  ScanRT( lData );
  // int x , y , z;
  // ClustScore->GetMaximumBin( x , y , z );
  // std::cout << x*dR << " " << y*dT << " " << z << std::endl;

  // std::map< const Data* , std::vector< Data* > > lClusters;
  // Cluster( lData , 0.0075 , 0.02 , lClusters );

  InteractiveDisplay( [ &Nclust ](){ DrawHisto( Nclust ); } , [ &ClustSize ](){ DrawHisto( ClustSize ); } , [ &ClustScore ](){ DrawHisto( ClustScore ); } );

  //InteractiveDisplay( [ &Nclust ](){ DrawHisto( Nclust ); } , [ &ClustSize ](){ DrawHisto( ClustSize ); } , [ &ClustScore ](){ DrawHisto( ClustScore ); } , [ lClusters ](){ DrawPoints( lClusters ); } );


  // InteractiveDisplay( [ lData ](){ DrawPoints( lData ); } );

  //ClusterAllData( lData , 0.0035 , 0.007 );

  // auto lHist = ProfileAllData( lData , 3e-2 );
  // InteractiveDisplay( [ lData ](){ DrawPoints( lData ); } , [ lHist ](){ DrawHisto( lHist ); } );

  // InteractiveDisplay( [](){ DrawWeights(); } );
  return 0;
}
