/* ===== Cluster sources ===== */
#include "Cluster_Data.hpp"
#include "Cluster_DataSources.hpp"
#include "Cluster_PlotTools.hpp"
#include "Cluster_GlobalVars.hpp"

// /* ===== C++ ===== */
#include <vector>

/* ===== For Root ===== */
#include "TH2D.h"
#include "Math/Interpolator.h" 
  
/* ===== Local utilities ===== */
#include "RootWindow.hpp"


void RTscanCallback( std::vector<Data>& aData , const double& aR , const double& aT , TH2D* ClustScore , TH2D* Nclust , TH2D* ClustSize )
{
  double lScore( 0.0 ) , lMean( 0.0 );
  std::size_t lCnt( 0 );

  for( auto& i : Data::Clusters )
  {
    if( i.mClusterSize > 1 )
    {
      lScore += i.mClusterScore;
      lCnt += 1;
      lMean += i.mClusterSize;
    }
  }

  // std::cout << std::setw(10) << aR << std::setw(10) << aT << std::setw(10) << lCount << std::setw(10) << lMean / lCount << std::endl;

  ClustScore -> Fill( aR , aT , lScore );
  Nclust -> Fill( aR , aT , lCnt );
  ClustSize -> Fill( aR , aT , lMean / lCnt );
}




/* ===== Main function ===== */
int main(int argc, char **argv)
{
  if (argc < 2) throw std::runtime_error( "Expecting a filename" );

  // ROOT::Math::Interpolator lInt( { 5_nanometer , 15_nanometer , 25_nanometer , 50_nanometer , 75_nanometer , 100_nanometer , 125_nanometer , 150_nanometer , 175_nanometer , 200_nanometer } , 
  //                                { 0.03631079  , 0.110302441  , 0.214839819  , 0.268302465  , 0.214839819  , 0.110302441  , 0.03631079   , 0.007664194  , 0.001037236  , 9.00054E-05 } ); // Default to cubic spline interpolation

  ROOT::Math::Interpolator lInt( { 0_nanometer , 20_nanometer , 30_nanometer , 40_nanometer , 50_nanometer , 60_nanometer , 70_nanometer , 80_nanometer , 90_nanometer , 100_nanometer } , 
                                 { 0.03631079  , 0.110302441  , 0.214839819  , 0.268302465  , 0.214839819  , 0.110302441  , 0.03631079   , 0.007664194  , 0.001037236  , 9.00054E-05 } ); // Default to cubic spline interpolation
  
  
  Parameters.SetZoom( 20_micrometer );
  Parameters.SetSigmaParameters( 100 , 5_nanometer , 100_nanometer , [ &lInt ]( const double& aPt ){ return lInt.Eval( aPt ); } );
  Parameters.SetMaxR( 200_nanometer );  
  Parameters.SetBins( 100 , 100 );
  Parameters.SetPbAlpha( 0.2 , 20 );

//  auto lData = LoadCSV( "1_un_red.csv" , 1./64000. , -1. , -1. ); // Full file
  auto lData = LoadCSV( argv[1] , 87_micrometer , 32_micrometer ); // One cluster
  //auto lData = LoadCSV( "1_un_red.csv" , 1./1000. , 87000. , 32000. ); // Very zoomed


  // // auto lData = CreatePseudoData( 10000 , 500 , 500 , 100_nanometer );
  // //auto lData = CreatePseudoData( 0 , 10 , 100 , 10_micrometer );
  // //auto lData = CreatePseudoData( 70000 , 700 , 700 , .005 );

  // // InteractiveDisplay( [ &lData ](){ DrawPoints( lData ); } );

  PrepData( lData );

  // // for( auto& i : lData ) std::cout << (i.neighbours[0].size()?sqrt(i.neighbours[0].begin()->first)*20_micrometer:0.0) << " " << (i.neighbours[1].size()?sqrt(i.neighbours[1].begin()->first)*20_micrometer:0.0) << std::endl;

  auto Rlo = Parameters.minScanR() - ( 0.5 * Parameters.dR() );
  auto Rhi = Parameters.maxScanR() - ( 0.5 * Parameters.dR() );
  auto Tlo = Parameters.minScanT() - ( 0.5 * Parameters.dT() );
  auto Thi = Parameters.maxScanT() - ( 0.5 * Parameters.dT() );

  auto Nclust     = new TH2D( "Nclust" , "N_{clusters};r;T" , Parameters.Rbins() , Rlo , Rhi , Parameters.Tbins() , Tlo , Thi );
  auto ClustSize  = new TH2D( "ClustSize" , "<N_{points}>;r;T" , Parameters.Rbins() , Rlo , Rhi , Parameters.Tbins() , Tlo , Thi );
  auto ClustScore = new TH2D( "ClustScore" , "Score;r;T" , Parameters.Rbins() , Rlo , Rhi , Parameters.Tbins() , Tlo , Thi );

  ScanRT( lData , [&]( const double& aR , const double& aT ){ RTscanCallback( lData , aR , aT , ClustScore , Nclust , ClustSize ); } );

  // // int x , y , z;
  // // ClustScore->GetMaximumBin( x , y , z );
  // // std::cout << x*dR << " " << y*dT << " " << z << std::endl;

  InteractiveDisplay( [&](){ DrawHisto( Nclust ); } , [&](){ DrawHisto( ClustSize ); } , [&](){ DrawHisto( ClustScore ); } );


  // // const double R( 15_micrometer*Parameters.scale() );
  // // const double T( 20_micrometer*Parameters.scale() );

  // const double R( 50_nanometer*Parameters.scale() );
  // const double T( 70_nanometer*Parameters.scale() );
  // std::cout << std::dec << "R=" << R << " | T=" << T << std::endl;

  // Cluster( lData , R , T );


  // // // std::set< Data* > lClusterable;
  // // // for( auto& i : lData ) 
  // // //   if( i.localizationscore > T )
  // // //     for( auto& j : i.neighbours )
  // // //       for( auto& k : j )
  // // //         if( k.first < 4.0*R*R and k.second->localizationscore > T)
  // // //         { 
  // // //           lClusterable.insert( k.second );
  // // //           lClusterable.insert( &i );
  // // //         }

  // // // std::cout << "Brute force "<< std::dec<< lClusterable.size() << std::endl;

  // std::map< const Data::Cluster* , std::vector< Data* > > lClusters;
  // // for( auto& i : Data::Cluster::Clusters ) std::cout << i.ClusterSize << std::endl;
  // for( auto& i : lData )
  //   lClusters[ i.mCluster ].push_back( &i );

  // // // for( auto& i : lData ) lClusters2[ ( lClusterable.find( &i ) != lClusterable.end() )? (Data*)(1) :(NULL) ].push_back( &i );

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
