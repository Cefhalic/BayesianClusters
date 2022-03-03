/* ===== Cluster sources ===== */
#include "Cluster_Data.hpp"
#include "Cluster_PlotTools.hpp"

// /* ===== C++ ===== */
#include <vector>

/* ===== For Root ===== */
#include "TH2D.h"
  
/* ===== Local utilities ===== */
#include "RootWindow.hpp"
#include "ProgressBar.hpp"


void RTscanCallback( const EventProxy& aEvent , const double& aR , const double& aT , TH2D* ClustScore , TH2D* Nclust , TH2D* ClustSize )
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


  std::cout << "+------------------------------------+" << std::endl;
  ProgressBar2 lBar( "| Cluster Scan. Andrew W. Rose. 2022 |" , 1 );
  std::cout << "+------------------------------------+" << std::endl;
  std::string lInputFile = Event::mParameters.FromCommandline( argc , argv );
  std::cout << "+------------------------------------+" << std::endl;

  Event lEvent( lInputFile );  
  
  // WriteCSV( std::string("trunc_")+argv[1] , lData );
  // return 0;

  // // auto lData = CreatePseudoData( 10000 , 500 , 500 , 100_nanometer );
  // //auto lData = CreatePseudoData( 0 , 10 , 100 , 10_micrometer );
  // //auto lData = CreatePseudoData( 70000 , 700 , 700 , .005 );

  // // InteractiveDisplay( [ &lData ](){ DrawPoints( lData ); } );

  auto Rlo = Event::mParameters.minScanR() - ( 0.5 * Event::mParameters.dR() );
  auto Rhi = Event::mParameters.maxScanR() - ( 0.5 * Event::mParameters.dR() );
  auto Tlo = Event::mParameters.minScanT() - ( 0.5 * Event::mParameters.dT() );
  auto Thi = Event::mParameters.maxScanT() - ( 0.5 * Event::mParameters.dT() );

  auto Nclust     = new TH2D( "Nclust" ,     "N_{clusters};r;T" , Event::mParameters.Rbins() , Rlo , Rhi , Event::mParameters.Tbins() , Tlo , Thi );
  auto ClustSize  = new TH2D( "ClustSize" ,  "<N_{points}>;r;T" , Event::mParameters.Rbins() , Rlo , Rhi , Event::mParameters.Tbins() , Tlo , Thi );
  auto ClustScore = new TH2D( "ClustScore" , "Score;r;T" ,        Event::mParameters.Rbins() , Rlo , Rhi , Event::mParameters.Tbins() , Tlo , Thi );


  lEvent.ScanRT( [&]( const EventProxy& aEvent , const double& aR , const double& aT ){ RTscanCallback( aEvent , aR , aT , ClustScore , Nclust , ClustSize ); } );
  // const double R( 50_nanometer*Parameters.scale() );
  // const double T( 70_nanometer*Parameters.scale() );
  // std::cout << std::dec << "R=" << R << " | T=" << T << std::endl;
  // std::vector< Cluster > lClusters;
  // lClusters.reserve( lData.size() );  // Reserve as much space for clusters as there are data points - prevent pointers being invalidated!
  // Clusterize( lData , lClusters , R , T );
  // CheckClusterization( lData , R , T );
  // InteractiveDisplay( [ &lData ](){ DrawPoints( lData ); } , [ &lData ](){ DrawClusters( lData ); } );


  std::cout << "+------------------------------------+" << std::endl;

}
