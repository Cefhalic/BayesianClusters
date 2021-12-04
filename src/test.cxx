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


#define PRINT( x ) std::cout << ( #x ) << " = " << ( x ) << std::endl;

double __Cluster__( std::vector<Data>& aData , const double& R , const double& T )
{
  const double twoR2 = 4.0 * R * R;

  // For real multi-threading want to use blocks, not interleaved, to minimize collisions
  // [ &twoR2 , &T ]( Data& i ){ i.Clusterize( twoR2 , T ); } && aData;

  for( auto& i : aData ) i.Clusterize( twoR2 , T ); 
  // for( auto& i : aData ) i.GetParent(); // Force update on all clusters;


  double lScore ( 0.0 );
  []( Data& i ){ i.UpdateClusterScore(); } || aData; 

  //std::cout << lScore << std::endl;

  return lScore;
}



void Cluster( std::vector<Data>& aData , const double& R , const double& T )
{
  const std::size_t Nminus1( aData.size() - 1 );
  const auto R2 = R * R;

  // Reset ahead of UpdateLocalization and Clusterize
  for( auto& i : aData ){
    i.ResetClusters();  
    i.neighbourit = i.neighbours[0].begin();
  }

  [ &Nminus1 , &R2 ]( Data& i ){ i.UpdateLocalization( R2 , Nminus1 ); } || aData; 

  __Cluster__( aData , R , T );
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

    [ &Nminus1 , &R2 ]( Data& i ){ i.UpdateLocalization( R2 , Nminus1 ); } || aData; 

    []( Data& i ){ i.ResetClusters(); } || aData;

    for( uint32_t j(0) ; j!=Parameters.Tbins() ; ++j , T-=Parameters.dT() , ++lProgressBar )
    {
      auto lScore = __Cluster__( aData , R , T );
      // ClustScore -> Fill( R , T , lScore );
      // Nclust -> Fill( R , T , lClusters.size() );
      // double mean( 0.0 );
      // for ( auto& i: lClusters ) if( i.first ) mean += i.second.size();
      // ClustSize -> Fill( R , T , mean / lClusters.size() );
    }

  }

}




void PrepData( std::vector<Data>& aData )
{
  ProgressBar2 lProgressBar( "Preprocessing" , aData.size() );

  // Crucial step is to sort!
  std::sort( aData.begin() , aData.end() );

  // Populate neighbour lists
  // Interleave threading since processing time increases with radius from origin
  [ &aData ]( const std::size_t& i ){ aData.at( i ).PopulateNeighbours( aData.begin() + i + 1 , aData.end() , aData.rbegin() + aData.size() - i , aData.rend() ); } || range( aData.size() ); 

  // ... any other preparation ...
}




/* ===== Main function ===== */
int main(int argc, char **argv)
{


  Parameters.SetSigmaCountAndSpacing( 10 , 1e-3 , 3e-2 );
  Parameters.SetProbabilitySigma( { 0.0000 , 0.0050 , 0.010 , 0.015 , 0.020 , 0.025 , 0.030 , 0.035 , 0.040 , 0.045 } , 
                                  { 0.03631079, 0.110302441, 0.214839819, 0.268302465, 0.214839819, 0.110302441, 0.03631079, 0.007664194, 0.001037236, 9.00054E-05 } );
  Parameters.SetMaxR( 0.02 );
  Parameters.SetBins( 40 , 40 );

//  auto lData = LoadCSV( "1_un_red.csv" , 1./64000. , -1. , -1. ); // Full file
   // auto lData = LoadCSV( "1_un_red.csv" , 1./10000. , 87000. , 32000. ); // One cluster
  //auto lData = LoadCSV( "1_un_red.csv" , 1./1000. , 87000. , 32000. ); // Very zoomed


  auto lData = CreatePseudoData( 10000 , 500 , 500 , .015 );
  // auto lData = CreatePseudoData( 100 , 20 , 500 , .01 );
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

  //InteractiveDisplay( [ &Nclust ](){ DrawHisto( Nclust ); } , [ &ClustSize ](){ DrawHisto( ClustSize ); } , [ &ClustScore ](){ DrawHisto( ClustScore ); } );

  //InteractiveDisplay( [ &Nclust ](){ DrawHisto( Nclust ); } , [ &ClustSize ](){ DrawHisto( ClustSize ); } , [ &ClustScore ](){ DrawHisto( ClustScore ); } , [ lClusters ](){ DrawPoints( lClusters ); } );


  // InteractiveDisplay( [ lData ](){ DrawPoints( lData ); } );

  //ClusterAllData( lData , 0.0035 , 0.007 );

  // auto lHist = ProfileAllData( lData , 3e-2 );
  // InteractiveDisplay( [ lData ](){ DrawPoints( lData ); } , [ lHist ](){ DrawHisto( lHist ); } );

  // InteractiveDisplay( [](){ DrawWeights(); } );
  return 0;
}
