/* ===== Cluster sources ===== */
#include "Cluster_Data.hpp"
#include "Cluster_DataSources.hpp"
#include "Cluster_PlotTools.hpp"

// /* ===== C++ ===== */
#include <vector>

/* ===== For Root ===== */
#include "TH2D.h"
#include "Math/Interpolator.h" 
  
/* ===== Local utilities ===== */
#include "RootWindow.hpp"


void RTscanCallback( const Event& aEvent , const double& aR , const double& aT , TH2D* ClustScore , TH2D* Nclust , TH2D* ClustSize )
{
  double lScore( 0.0 ) , lMean( 0.0 );
  std::size_t lCnt( 0 );

  for( auto& i : aEvent.mClusters )
  {
    if( i.mClusterSize )
    {
      lScore += i.mClusterScore;
      lCnt += 1;
      lMean += i.mClusterSize;
    }
  }

  ClustScore -> Fill( aR , aT , lScore );
  Nclust -> Fill( aR , aT , lCnt );
  ClustSize -> Fill( aR , aT , lMean / lCnt );
}




/* ===== Main function ===== */
int main(int argc, char **argv)
{
  if (argc < 2) throw std::runtime_error( "Expecting a filename" );


  ROOT::Math::Interpolator lInt( { 0_nanometer , 20_nanometer , 30_nanometer , 40_nanometer , 50_nanometer , 60_nanometer , 70_nanometer , 80_nanometer , 90_nanometer , 100_nanometer } , 
                                 { 0.03631079  , 0.110302441  , 0.214839819  , 0.268302465  , 0.214839819  , 0.110302441  , 0.03631079   , 0.007664194  , 0.001037236  , 9.00054E-05 } ); // Default to cubic spline interpolation

  Event::mParameters.SetZoom( 2_micrometer );
  Event::mParameters.SetMaxR( 200_nanometer );  
  Event::mParameters.SetBins( 35 , 35 );
  Event::mParameters.SetPbAlpha( 0.2 , 20 );
  Event::mParameters.SetValidate( true );

  Event::mParameters.SetSigmaParameters( 100 , 5_nanometer , 100_nanometer , [ &lInt ]( const double& aPt ){ return lInt.Eval( aPt ); } );
  // for( const auto& i :  Event::mParameters.sigmabins() ) std::cout << i << std::endl;


  Event lEvent( 87_micrometer , 32_micrometer );  
  LoadCSV( argv[1] , lEvent ); // One cluster

  // WriteCSV( "trunc_"+std::string(argv[1]) , lData );
  // return 0;

  // // auto lData = CreatePseudoData( 10000 , 500 , 500 , 100_nanometer );
  // //auto lData = CreatePseudoData( 0 , 10 , 100 , 10_micrometer );
  // //auto lData = CreatePseudoData( 70000 , 700 , 700 , .005 );

  // // InteractiveDisplay( [ &lData ](){ DrawPoints( lData ); } );

  lEvent.PrepData();

  auto Rlo = Event::mParameters.minScanR() - ( 0.5 * Event::mParameters.dR() );
  auto Rhi = Event::mParameters.maxScanR() - ( 0.5 * Event::mParameters.dR() );
  auto Tlo = Event::mParameters.minScanT() - ( 0.5 * Event::mParameters.dT() );
  auto Thi = Event::mParameters.maxScanT() - ( 0.5 * Event::mParameters.dT() );

  auto Nclust     = new TH2D( "Nclust" ,     "N_{clusters};r;T" , Event::mParameters.Rbins() , Rlo , Rhi , Event::mParameters.Tbins() , Tlo , Thi );
  auto ClustSize  = new TH2D( "ClustSize" ,  "<N_{points}>;r;T" , Event::mParameters.Rbins() , Rlo , Rhi , Event::mParameters.Tbins() , Tlo , Thi );
  auto ClustScore = new TH2D( "ClustScore" , "Score;r;T" ,        Event::mParameters.Rbins() , Rlo , Rhi , Event::mParameters.Tbins() , Tlo , Thi );

  lEvent.ScanRT( [&]( const Event& aEvent , const double& aR , const double& aT ){ RTscanCallback( aEvent , aR , aT , ClustScore , Nclust , ClustSize ); } );


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
