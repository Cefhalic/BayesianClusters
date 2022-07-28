/* ===== Cluster sources ===== */
#include "BayesianClustering/Cluster.hpp"
#include "BayesianClustering/Event.hpp"
#include "BayesianClustering/EventProxy.hpp"
#include "BayesianClustering/Configuration.hpp"

// /* ===== C++ ===== */
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <mutex>
  
/* ===== Local utilities ===== */
#include "Utilities/ProgressBar.hpp"


std::mutex mtx; // mutex for critical section



void XmlCallback( const EventProxy& aEvent , const double& aR , const double& aT, std::vector<uint32_t>& aCurrentIJ , std::stringstream& aOutput, 
                  std::vector<std::vector<double>>& aRTScores, std::vector<uint32_t>& aMaxScorePosition, double& aMaxRTScore )
{
  mtx.lock();
  aOutput << "  { R:" << aR << ", T:" << aT << ", Score:" << aEvent.mLogP << ", NumClusteredPts:" << aEvent.mClusteredCount << ", NumBackgroundPts:" << aEvent.mBackgroundCount << ", Clusters:[\n";

  for( auto& i : aEvent.mClusters )
  {
    if( i.mClusterSize ) aOutput << "    { Points:" << i.mClusterSize << ",  Score:" << i.mClusterScore << " },\n";
  }
  
  // double lLogP = aEvent.mLogP;
  // //score setting stuff
  // if (lLogP > aMaxRTScore){
  //   aMaxScorePosition = aCurrentIJ;
  //   aMaxRTScore = lLogP;
  // }

  // uint32_t p = aCurrentIJ[0], q = aCurrentIJ[1];
  // aRTScores[p][q] = lLogP;

  aOutput << "  ] },\n";
  mtx.unlock();  
}


void JsonCallback( const EventProxy& aEvent , const double& aR , const double& aT, std::vector<uint32_t>& aCurrentIJ , std::stringstream& aOutput, 
                  std::vector<std::vector<double>>& aRTScores, std::vector<uint32_t>& aMaxScorePosition, double& aMaxRTScore )
{
  mtx.lock();
  aOutput << "  { R:" << aR << ", T:" << aT << ", Score:" << aEvent.mLogP << ", NumClusteredPts:" << aEvent.mClusteredCount << ", NumBackgroundPts:" << aEvent.mBackgroundCount << ", Clusters:[\n";

  for( auto& i : aEvent.mClusters )
  {
    if( i.mClusterSize ) aOutput << "    { Points:" << i.mClusterSize << ",  Score:" << i.mClusterScore << " },\n";
  }
  
  // double lLogP = aEvent.mLogP;
  // //score setting stuff
  // if (lLogP > aMaxRTScore){
  //   aMaxScorePosition = aCurrentIJ;
  //   aMaxRTScore = lLogP;
  // }

  // uint32_t p = aCurrentIJ[0], q = aCurrentIJ[1];
  // aRTScores[p][q] = lLogP;

  aOutput << "  ] },\n";
  mtx.unlock();  
}






/* ===== Main function ===== */
int main(int argc, char **argv)
{


  std::cout << "+------------------------------------+" << std::endl;
  ProgressBar2 lBar( "| Cluster Scan. Andrew W. Rose. 2022 |" , 1 );
  std::cout << "+------------------------------------+" << std::endl;
  Configuration::Instance.FromCommandline( argc , argv );
  std::cout << "+------------------------------------+" << std::endl;

  Event lEvent;  
  std::vector<std::vector<double>> lRTScores(Configuration::Instance.Rbins()
                                  ,std::vector<double>(Configuration::Instance.Tbins()/*, 1*/));
  std::vector<uint32_t> lMaxScorePosition(2,0);
  double lMaxRTScore = -9E99;
  //the above will store our scores - it needs to end up in the callback

  const std::string& lFilename = Configuration::Instance.outputFile();

  if( lFilename.size() == 0 )
  {
    std::cout << "Warning: Running scan without callback" << std::endl;
    lEvent.ScanRT( [&]( const EventProxy& aEvent , const double& aR , const double& aT, std::vector<uint32_t> aCurrentIJ){} ); // Null callback
  }
  else if( lFilename.size() > 4 and lFilename.substr(lFilename.size() - 4) == ".xml" )
  {
    std::stringstream lOutput;
    lEvent.ScanRT( [&]( const EventProxy& aEvent , const double& aR , const double& aT, std::vector<uint32_t> aCurrentIJ ){ JsonCallback( aEvent , aR , aT, aCurrentIJ  , lOutput, lRTScores, lMaxScorePosition, lMaxRTScore); } );
    std::ofstream lOutFile( lFilename );
    lOutFile << "<Results>\n" << lOutput.str() << "</Results>\n";
  }
  else if( lFilename.size() > 5 and lFilename.substr(lFilename.size() - 5) == ".json" )
  {
    std::stringstream lOutput;
    lEvent.ScanRT( [&]( const EventProxy& aEvent , const double& aR , const double& aT, std::vector<uint32_t> aCurrentIJ ){ JsonCallback( aEvent , aR , aT, aCurrentIJ  , lOutput, lRTScores, lMaxScorePosition, lMaxRTScore); } );
    std::ofstream lOutFile( lFilename );
    lOutFile << "{\nResults:[\n" << lOutput.str() << "]\n}";
  }
  else
  {
    throw std::runtime_error( "No handler for specified output-file" );
  }

  std::cout << "+------------------------------------+" << std::endl;

}
