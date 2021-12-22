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


void RTscanCallback( const std::vector< Cluster >& aClusters , const double& aR , const double& aT , TH2D* ClustScore , TH2D* Nclust , TH2D* ClustSize )
{
  double lScore( 0.0 ) , lMean( 0.0 );
  std::size_t lCnt( 0 );

  for( auto& i : aClusters )
  {
    if( i.mClusterSize )
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

  auto Rlo = Parameters.minScanR() - ( 0.5 * Parameters.dR() );
  auto Rhi = Parameters.maxScanR() - ( 0.5 * Parameters.dR() );
  auto Tlo = Parameters.minScanT() - ( 0.5 * Parameters.dT() );
  auto Thi = Parameters.maxScanT() - ( 0.5 * Parameters.dT() );

  auto Nclust     = new TH2D( "Nclust" , "N_{clusters};r;T" , Parameters.Rbins() , Rlo , Rhi , Parameters.Tbins() , Tlo , Thi );
  auto ClustSize  = new TH2D( "ClustSize" , "<N_{points}>;r;T" , Parameters.Rbins() , Rlo , Rhi , Parameters.Tbins() , Tlo , Thi );
  auto ClustScore = new TH2D( "ClustScore" , "Score;r;T" , Parameters.Rbins() , Rlo , Rhi , Parameters.Tbins() , Tlo , Thi );

  ScanRT( lData , [&]( const std::vector< Cluster >& aClusters , const double& aR , const double& aT ){ RTscanCallback( aClusters , aR , aT , ClustScore , Nclust , ClustSize ); } );


  // const double R( 50_nanometer*Parameters.scale() );
  // const double T( 70_nanometer*Parameters.scale() );
  // std::cout << std::dec << "R=" << R << " | T=" << T << std::endl;
  // std::vector< Cluster > lClusters;
  // lClusters.reserve( lData.size() );  // Reserve as much space for clusters as there are data points - prevent pointers being invalidated!
  // Clusterize( lData , lClusters , R , T );
  // CheckClusterization( lData , R , T );
  // InteractiveDisplay( [ &lData ](){ DrawPoints( lData ); } , [ &lData ](){ DrawClusters( lData ); } );

   return 0;
}
