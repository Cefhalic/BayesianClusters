/* ===== Cluster sources ===== */
#include "Cluster_Data.hpp"
#include "Cluster_PlotTools.hpp"

// /* ===== C++ ===== */
#include <vector>
#include <stdio.h>
#include <mutex>

/* ===== For Root ===== */
// #include "TH2D.h"
  
/* ===== Local utilities ===== */
#include "RootWindow.hpp"
#include "ProgressBar.hpp"


// void RTscanCallback( const EventProxy& aEvent , const double& aR , const double& aT , TH2D* ClustScore , TH2D* Nclust , TH2D* ClustSize )
// {
//   double lScore( 0.0 ) , lMean( 0.0 );
//   std::size_t lCnt( 0 );

//   for( auto& i : aEvent.mClusters )
//   {
//     if( i.mClusterSize )
//     {
//       lScore += i.mClusterScore;
//       lCnt += 1;
//       lMean += i.mClusterSize;
//     }
//   }

//   ClustScore -> Fill( aR , aT , lScore );
//   Nclust -> Fill( aR , aT , lCnt );
//   ClustSize -> Fill( aR , aT , lMean / lCnt );
// }



std::mutex mtx; // mutex for critical section
// void RTscanCallback2( const EventProxy& aEvent , const double& aR , const double& aT , FILE* aFile )
// {
//   uint32_t lCount(0);
//   for( auto& i : aEvent.mClusters )
//   {
//     if( i.mClusterSize ) lCount++;
//   }

//   mtx.lock();
//   fprintf( aFile , "<Scan R='%.9f' T='%.9f' Score='%.9f'>\n" , aR , aT , aEvent.mLogP );
//   fprintf( aFile , " <Localizations Clustered='%u' Background='%u' />\n" , aEvent.mClusteredCount , aEvent.mBackgroundCount );
//   fprintf( aFile , " <Clusters Total='%u' Valid='%u'>\n" , aEvent.mClusterCount , lCount );
//   for( auto& i : aEvent.mClusters )
//   {
//     if( i.mClusterSize ) fprintf( aFile , "  <Cluster Points='%u' Score='%.9f' />\n" , i.mClusterSize , i.mClusterScore );
//   }
//   fprintf( aFile , " </Clusters>\n</Scan>\n" , aR , aT , aEvent.mLogP );

//   mtx.unlock();  
// }


void RTscanCallback3( const EventProxy& aEvent , const double& aR , const double& aT , std::map< std::pair< double , double > , EventProxy >& aRes )
{
  uint32_t lCount(0);
  for( auto& i : aEvent.mClusters )
  {
    if( i.mClusterSize ) lCount++;
  }

  mtx.lock();
  aRes.emplace( std::make_pair( std::make_pair( aR , aT ) , aEvent.LightweightClone() ) );
  mtx.unlock();  
}



/* ===== Main function ===== */
int main(int argc, char **argv)
{

  // std::cout << "\xe2\x95\x94" << std::string( "\xe2\x95\x90" , 100 ) << std::endl;
  // return 1;

  std::cout << "+------------------------------------+" << std::endl;
  ProgressBar2 lBar( "| Cluster Scan. Andrew W. Rose. 2022 |" , 1 );
  std::cout << "+------------------------------------+" << std::endl;
  std::string lInputFile = Event::mParameters.FromCommandline( argc , argv );
  std::cout << "+------------------------------------+" << std::endl;

  Event lEvent( lInputFile );  
  
  // // InteractiveDisplay( [ &lData ](){ DrawPoints( lData ); } );

  // auto Rlo = Event::mParameters.minScanR() - ( 0.5 * Event::mParameters.dR() );
  // auto Rhi = Event::mParameters.maxScanR() - ( 0.5 * Event::mParameters.dR() );
  // auto Tlo = Event::mParameters.minScanT() - ( 0.5 * Event::mParameters.dT() );
  // auto Thi = Event::mParameters.maxScanT() - ( 0.5 * Event::mParameters.dT() );

  // auto Nclust     = new TH2D( "Nclust" ,     "N_{clusters};r;T" , Event::mParameters.Rbins() , Rlo , Rhi , Event::mParameters.Tbins() , Tlo , Thi );
  // auto ClustSize  = new TH2D( "ClustSize" ,  "<N_{points}>;r;T" , Event::mParameters.Rbins() , Rlo , Rhi , Event::mParameters.Tbins() , Tlo , Thi );
  // auto ClustScore = new TH2D( "ClustScore" , "Score;r;T" ,        Event::mParameters.Rbins() , Rlo , Rhi , Event::mParameters.Tbins() , Tlo , Thi );

  // lEvent.ScanRT( [&]( const EventProxy& aEvent , const double& aR , const double& aT ){ RTscanCallback( aEvent , aR , aT , ClustScore , Nclust , ClustSize ); } );



  FILE* lFile = fopen( "ScanResults.txt" , "w" );
  // lEvent.ScanRT( [&]( const EventProxy& aEvent , const double& aR , const double& aT ){ RTscanCallback2( aEvent , aR , aT , lFile ); } );

  std::map< std::pair< double , double > , EventProxy > lRes;
  lEvent.ScanRT( [&]( const EventProxy& aEvent , const double& aR , const double& aT ){ RTscanCallback3( aEvent , aR , aT , lRes ); } );

  for( auto& i : lRes )
  {
    fprintf( lFile , "<Scan R='%.9f' T='%.9f' Score='%.9f'>\n" , i.first.first , i.first.second , i.second.mLogP );
    fprintf( lFile , " <Localizations Clustered='%u' Background='%u' />\n" , i.second.mClusteredCount , i.second.mBackgroundCount );
    fprintf( lFile , " <Clusters Count='%u'>\n" , i.second.mClusters.size() );
    for( auto& j : i.second.mClusters )
    {
      if( j.mClusterSize ) fprintf( lFile , "  <Cluster Points='%u' Score='%.9f' />\n" , j.mClusterSize , j.mClusterScore );
    }
    fprintf( lFile , " </Clusters>\n</Scan>\n" );
  }
  fclose( lFile );



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
