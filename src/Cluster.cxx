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
  
  double lLogP = aEvent.mLogP;
  //score setting stuff
  if (lLogP > aMaxRTScore){
    aMaxScorePosition = aCurrentIJ;
    aMaxRTScore = lLogP;
  }

  uint32_t p = aCurrentIJ[0], q = aCurrentIJ[1];
  aRTScores[p][q] = lLogP;

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
  
  double lLogP = aEvent.mLogP;
  //score setting stuff
  if (lLogP > aMaxRTScore){
    aMaxScorePosition = aCurrentIJ;
    aMaxRTScore = lLogP;
  }

  uint32_t p = aCurrentIJ[0], q = aCurrentIJ[1];
  aRTScores[p][q] = lLogP;

  aOutput << "  ] },\n";
  mtx.unlock();  
}

//just prints the best R, T at the end
std::vector<double> bestRT(std::vector<uint32_t>& aMaxScorePosition, std::vector<std::vector<double>>& aRTScores){
  int i = aMaxScorePosition[0], j = aMaxScorePosition[1];
  if ((i < 2) or (j < 2)) return {1.0, 1.0};
  if ((i > Configuration::Instance.Rbins() - 3 ) or (j > Configuration::Instance.Tbins() - 3)) return {1.0, 1.0};

  // double lRIndex(0), lTIndex(0);
  double lRValue(0), lTValue(0);
  double lRValueSum(0), lTValueSum(0);
  for(int index(-2); index < 3 ; ++index ){
    if(index == 0) continue;
    lRValueSum += aRTScores[i + index][j];
    lTValueSum += aRTScores[i][j + index];
    lRValue += (index + i) * aRTScores[i + index][j];
    lTValue += (index + i) * aRTScores[i][j + index];
  }

  double lRIndex = double(lRValue / lRValueSum);
  double lTIndex = double (lTValue /  lTValueSum);

  int lR = int(lRIndex), lT = int(lTIndex);

  double outputR = Configuration::Instance.minScanR() + (lR * Configuration::Instance.dR());
  double outputT = Configuration::Instance.minScanT() + (lT * Configuration::Instance.dT());

  return {Configuration::Instance.toPhysicalUnits(outputR), Configuration::Instance.toPhysicalUnits(outputT)};
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

  std::cout << "max score was: " << lMaxRTScore << std::endl;
  std::cout << "at position (" << lMaxScorePosition[0] << ", " << lMaxScorePosition[1] << ")"<< std::endl;
  std::vector<double> a(2);
  a = bestRT(lMaxScorePosition, lRTScores);
  std::cout << a[0] << " " << a[1] << std::endl;
  std::cout << "+------------------------------------+" << std::endl;

}
