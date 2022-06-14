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
#include "ProgressBar.hpp"


std::mutex mtx; // mutex for critical section



void XmlCallback( const EventProxy& aEvent , const double& aR , const double& aT , std::stringstream& aOutput )
{
  mtx.lock();
  aOutput << "  <Scan R='" << aR << "' T='" << aT << "' Score='" << aEvent.mLogP << " NumClusteredPts='" << aEvent.mClusteredCount << "' NumBackgroundPts='" << aEvent.mBackgroundCount << "'>\n";

  for( auto& i : aEvent.mClusters )
  {
    if( i.mClusterSize ) aOutput << "      <Cluster Points='" << i.mClusterSize << "' Score='" << i.mClusterScore << "' />\n";
  }

  aOutput << "  </Scan>\n";
  mtx.unlock();  
}


void JsonCallback( const EventProxy& aEvent , const double& aR , const double& aT , std::stringstream& aOutput )
{
  mtx.lock();
  aOutput << "  { R:" << aR << ", T:" << aT << ", Score:" << aEvent.mLogP << ", NumClusteredPts:" << aEvent.mClusteredCount << ", NumBackgroundPts:" << aEvent.mBackgroundCount << ", Clusters:[\n";

  for( auto& i : aEvent.mClusters )
  {
    if( i.mClusterSize ) aOutput << "    { Points:" << i.mClusterSize << ",  Score:" << i.mClusterScore << " },\n";
  }

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


  const std::string& lFilename = Configuration::Instance.outputFile();

  if( lFilename.size() == 0 )
  {
    std::cout << "Warning: Running scan without callback" << std::endl;
    lEvent.ScanRT( [&]( const EventProxy& aEvent , const double& aR , const double& aT ){} ); // Null callback
  }
  else if( lFilename.size() > 4 and lFilename.substr(lFilename.size() - 4) == ".xml" )
  {
    std::stringstream lOutput;
    lEvent.ScanRT( [&]( const EventProxy& aEvent , const double& aR , const double& aT ){ XmlCallback( aEvent , aR , aT , lOutput ); } );
    std::ofstream lOutFile( lFilename );
    lOutFile << "<Results>\n" << lOutput.str() << "</Results>\n";
  }
  else if( lFilename.size() > 5 and lFilename.substr(lFilename.size() - 5) == ".json" )
  {
    std::stringstream lOutput;
    lEvent.ScanRT( [&]( const EventProxy& aEvent , const double& aR , const double& aT ){ JsonCallback( aEvent , aR , aT , lOutput ); } );
    std::ofstream lOutFile( lFilename );
    lOutFile << "{\nResults:[\n" << lOutput.str() << "]\n}";
  }
  else
  {
    throw std::runtime_error( "No handler for specified output-file" );
  }

  std::cout << "+------------------------------------+" << std::endl;

}
