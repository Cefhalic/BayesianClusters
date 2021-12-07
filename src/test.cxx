/* ===== Cluster sources ===== */
#include "Cluster_Data.hpp"
#include "Cluster_DataSources.hpp"
#include "Cluster_PlotTools.hpp"
#include "Cluster_GlobalVars.hpp"

// /* ===== C++ ===== */
#include <thread>
#include <iostream>
#include <math.h>
#include <vector>
#include <map>
#include <set>

/* ===== For Root ===== */
#include "TH2D.h"
#include "Math/Interpolator.h" 

/* ===== Local utilities ===== */
#include "ListComprehension.hpp"
#include "ProgressBar.hpp"
#include "RootWindow.hpp"
#include "Vectorize.hpp"


#define PRINT( x ) std::cout << ( #x ) << " = " << ( x ) << std::endl;


void Cluster( std::vector<Data>& aData , const double& R , const double& T )
{
  // Reset ahead of UpdateLocalization and Clusterize
  const std::size_t N( aData.size()-1 );
  []( Data& i ){ i.ResetClusters(); i.localizationsum = i.localizationscore = 0.0; i.neighbourit = i.neighbours[0].begin(); } || aData;
  [ &R , &N ]( Data& i ){ i.UpdateLocalization( R * R , N ); } || aData; // Use interleaving threading to average over systematic radial scaling

  //for( auto& i: aData ) i.Clusterize( 4.0 * R * R , T );
  [ &R , &T ]( Data& i ){ i.Clusterize( 4.0 * R * R , T ); } && aData; // Use block threading to minimize mutex collisions
}



TH2D *Nclust, *ClustSize, *ClustScore;


void ScanRT( std::vector<Data>& aData )
{
  const std::size_t Nminus1( aData.size() - 1 );

  double R( Parameters.minScanR() ) , R2( 0 ) , twoR2( 0 ) , T( 0 );

  ProgressBar lProgressBar( "Scan over RT" , Parameters.Rbins() * Parameters.Tbins() );

  for( uint32_t i(0) ; i!=Parameters.Rbins() ; ++i , R+=Parameters.dR() )
  {
    R2 = R * R;
    twoR2 = 4.0 * R2;
    T = Parameters.maxScanT();

    [ &Nminus1 , &R2 ]( Data& i ){ i.UpdateLocalization( R2 , Nminus1 ); } || aData; // Use interleaving threading to average over systematic radial scaling

    // for( auto& i : aData ) i.ResetClusters();
    []( Data& i ){ i.ResetClusters(); } || aData;

    for( uint32_t j(0) ; j!=Parameters.Tbins() ; ++j , T-=Parameters.dT() , ++lProgressBar )
    {

      [ &twoR2 , &T ]( Data& i ){ i.Clusterize( twoR2 , T ); } && aData; // Use block threading to minimize mutex collisions
      []( Data& i ){ i.UpdateClusterScore(); } || aData; // Use interleaving threading to average over systematic radial scaling

      double lScore( 0.0 ) , lMean( 0.0 );
      std::size_t lCount( 0 );

      for( auto& i : aData )
      {
        if( i.children.size() > 1 )
        {
          lScore += i.ClusterScore;
          lCount += 1;
          lMean += i.children.size();
        }
      }

      ClustScore -> Fill( R , T , lScore );
      Nclust -> Fill( R , T , lCount );
      ClustSize -> Fill( R , T , lMean / lCount );
    }

  }

}




void PrepData( std::vector<Data>& aData )
{

  // Should already be sorted, but...
  // std::sort( aData.begin() , aData.end() );

  {
    ProgressBar2 lProgressBar( "Preprocessing" , aData.size() );

    // Populate neighbour lists  
    [ &aData ]( const std::size_t& i ){ aData.at( i ).PopulateNeighbours( aData.begin() + i + 1 , aData.end() , aData.rbegin() + aData.size() - i , aData.rend() ); } || range( aData.size() );  // Interleave threading since processing time increases with radius from origin
  }

  // ... any other preparation ...
}




/* ===== Main function ===== */
int main(int argc, char **argv)
{
  if (argc < 2) throw std::runtime_error( "Expecting a filename" );

  // ROOT::Math::Interpolator lInt( { 5_nanometer , 15_nanometer , 25_nanometer , 50_nanometer , 75_nanometer , 100_nanometer , 125_nanometer , 150_nanometer , 175_nanometer , 200_nanometer } , 
  //                                { 0.03631079  , 0.110302441  , 0.214839819  , 0.268302465  , 0.214839819  , 0.110302441  , 0.03631079   , 0.007664194  , 0.001037236  , 9.00054E-05 } ); // Default to cubic spline interpolation

  ROOT::Math::Interpolator lInt( { 0_nanometer , 200_nanometer , 300_nanometer , 400_nanometer , 500_nanometer , 600_nanometer , 700_nanometer , 800_nanometer , 900_nanometer , 1000_nanometer } , 
                                 { 0.03631079  , 0.110302441  , 0.214839819  , 0.268302465  , 0.214839819  , 0.110302441  , 0.03631079   , 0.007664194  , 0.001037236  , 9.00054E-05 } ); // Default to cubic spline interpolation
  
  
  Parameters.SetZoom( 20_micrometer );
  Parameters.SetSigmaParameters( 100 , 5_nanometer , 50_nanometer , [ &lInt ]( const double& aPt ){ return lInt.Eval( aPt ); } );
  Parameters.SetMaxR( 20_nanometer );
  //Parameters.SetMaxR( 30_micrometer );  
  Parameters.SetBins( 20 , 20 );

//  auto lData = LoadCSV( "1_un_red.csv" , 1./64000. , -1. , -1. ); // Full file
  auto lData = LoadCSV( argv[1] , 87_micrometer , 32_micrometer ); // One cluster
  //auto lData = LoadCSV( "1_un_red.csv" , 1./1000. , 87000. , 32000. ); // Very zoomed


  // auto lData = CreatePseudoData( 10000 , 500 , 500 , 100_nanometer );
  //auto lData = CreatePseudoData( 0 , 10 , 100 , 10_micrometer );
  //auto lData = CreatePseudoData( 70000 , 700 , 700 , .005 );

  // InteractiveDisplay( [ &lData ](){ DrawPoints( lData ); } );

  PrepData( lData );

  // for( auto& i : lData ) std::cout << (i.neighbours[0].size()?sqrt(i.neighbours[0].begin()->first)*20_micrometer:0.0) << " " << (i.neighbours[1].size()?sqrt(i.neighbours[1].begin()->first)*20_micrometer:0.0) << std::endl;


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

  // InteractiveDisplay( [ &Nclust ](){ DrawHisto( Nclust ); } , [ &ClustSize ](){ DrawHisto( ClustSize ); } , [ &ClustScore ](){ DrawHisto( ClustScore ); } );


  // const double R( 15_micrometer*Parameters.scale() );
  // const double T( 20_micrometer*Parameters.scale() );

  // const double R( 50_nanometer*Parameters.scale() );
  // const double T( 70_nanometer*Parameters.scale() );
  // std::cout << std::dec << "R=" << R << " | T=" << T << std::endl;

  // Cluster( lData , R , T );


  // std::set< Data* > lClusterable;
  // for( auto& i : lData ) 
  //   if( i.localizationscore > T )
  //     for( auto& j : i.neighbours )
  //       for( auto& k : j )
  //         if( k.first < 4.0*R*R and k.second->localizationscore > T)
  //         { 
  //           lClusterable.insert( k.second );
  //           lClusterable.insert( &i );
  //         }

  // std::cout << "Brute force "<< std::dec<< lClusterable.size() << std::endl;

  // std::map< const Data* , std::vector< Data* > > lClusters , lClusters2;
  // for( auto& i : lData )
  //   for( auto& j : i.children )
  //     lClusters[ ( i.children.size() > 1)?(&i):(NULL) ].push_back( j );

  // // for( auto& i : lData ) lClusters2[ ( lClusterable.find( &i ) != lClusterable.end() )? (Data*)(1) :(NULL) ].push_back( &i );

  // // double sum(0);
  // // for( auto& i : lClusters ) if( i.first ) sum += i.second.size();
  // // std::cout << "Recursive " << std::dec << lClusters.size()-1 << " : " << sum/(lClusters.size()-1) << std::endl;


  // // InteractiveDisplay( [ &lClusters2 ](){ DrawPoints( lClusters2 ); } , [ &lClusters ](){ DrawPoints( lClusters ); } );



  // InteractiveDisplay( [ &lData ](){ DrawPoints( lData ); } , [ &lClusters ](){ DrawPoints( lClusters ); } );

  // //InteractiveDisplay( [ &lData ](){ DrawPoints( lData ); } , [ &lClusters2 ](){ DrawPoints( lClusters2 ); } );

  // // InteractiveDisplay( [ &Nclust ](){ DrawHisto( Nclust ); } , [ &ClustSize ](){ DrawHisto( ClustSize ); } , [ &ClustScore ](){ DrawHisto( ClustScore ); } , [ lClusters ](){ DrawPoints( lClusters ); } );
  // // InteractiveDisplay( [](){ DrawWeights(); } );
  return 0;
}
