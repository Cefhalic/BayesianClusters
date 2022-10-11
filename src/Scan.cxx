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



void XmlCallback( const EventProxy& aEvent , const double& aR , const double& aT, std::pair<int, int>& aCurrentIJ , std::stringstream& aOutput, 
                  std::vector<std::vector<double>>& aRTScores, std::pair<int,int>& aMaxScorePosition, double& aMaxRTScore )
{
  mtx.lock();
  // aOutput << "  { R:" << aR << ", T:" << aT << ", Score:" << aEvent.mLogP << ", NumClusteredPts:" << aEvent.mClusteredCount << ", NumBackgroundPts:" << aEvent.mBackgroundCount << ", Clusters:[\n";
  aOutput << "  { R:" << aR << ", T:" << aT << ", Score:" << aEvent.mLogP << ", NumClusteredPts:" << aEvent.mClusteredCount << ", NumBackgroundPts:" << aEvent.mBackgroundCount << "}\n";

  // for( auto& i : aEvent.mClusters )
  // {
  //   if( i.mClusterSize ) aOutput << "    { Points:" << i.mClusterSize << ",  Score:" << i.mClusterScore << " },\n";
  // }
  
  double lLogP = aEvent.mLogP;
  //score setting stuff
  if (lLogP > aMaxRTScore){
    aMaxScorePosition = aCurrentIJ;
    aMaxRTScore = lLogP;
  }

  uint32_t p = aCurrentIJ.first, q = aCurrentIJ.second;
  aRTScores[p][q] = lLogP;

  // aOutput << "   },\n";
  mtx.unlock();  
}


void JsonCallback( const EventProxy& aEvent , const double& aR , const double& aT, std::pair<int,int>& aCurrentIJ , std::stringstream& aOutput, 
                  std::vector<std::vector<double>>& aRTScores, std::pair<int,int>& aMaxScorePosition, double& aMaxRTScore )
{
  mtx.lock();
  // aOutput << "  { R:" << aR << ", T:" << aT << ", Score:" << aEvent.mLogP << ", NumClusteredPts:" << aEvent.mClusteredCount << ", NumBackgroundPts:" << aEvent.mBackgroundCount << ", Clusters:[\n";
  aOutput << "  { R:" << aR << ", T:" << aT << ", Score:" << aEvent.mLogP << ", NumClusteredPts:" << aEvent.mClusteredCount << ", NumBackgroundPts:" << aEvent.mBackgroundCount << "}\n";

  // for( auto& i : aEvent.mClusters )
  // {
  //   if( i.mClusterSize ) aOutput << "    { Points:" << i.mClusterSize << ",  Score:" << i.mClusterScore << " },\n";
  // }
  
  double lLogP = aEvent.mLogP;
  //score setting stuff
  if (lLogP > aMaxRTScore){
    aMaxScorePosition = aCurrentIJ;
    aMaxRTScore = lLogP;
  }

  uint32_t p = aCurrentIJ.first, q = aCurrentIJ.second;
  aRTScores[p][q] = lLogP;

  // aOutput << "   },\n";
  mtx.unlock();  
}

//just prints the best R, T at the end
std::pair<double,double> bestRT(std::pair<int, int>& aMaxScorePosition, std::vector<std::vector<double>>& aRTScores){
  int i = aMaxScorePosition.first, j = aMaxScorePosition.second;
  double lRValue(0), lTValue(0);
  double lValueSum(0);
  for(int I(-2); I < 3 ; ++I ){
    if ((i < 2) or (j < 2)) break;
    if ((i > Configuration::Instance.Rbins() - 3 ) or (j > Configuration::Instance.Tbins() - 3)) break;
    for(int J(-2); J < 3 ; ++J ){
      auto lVal = aRTScores[i + I][j + J];
      lValueSum += lVal;
      lRValue += (i+I) * lVal;
      lTValue += (j+J) * lVal;
  }}

  double lRIndex = lRValue / lValueSum;
  double lTIndex = lTValue /  lValueSum;

  double outputR = Configuration::Instance.minScanR() + (lRIndex * Configuration::Instance.dR());
  double outputT = Configuration::Instance.maxScanT() - (lTIndex * Configuration::Instance.dT());

  return std::make_pair(Configuration::Instance.toPhysicalUnits(outputR), Configuration::Instance.toPhysicalUnits(outputT));
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
  std::vector<std::vector<double>> lRTScores(Configuration::Instance.Rbins(),
                                            std::vector<double>(Configuration::Instance.Tbins()/*, 1*/));
  std::pair<int, int> lMaxScorePosition;
  double lMaxRTScore = -9E99;
  //the above will store our scores - it needs to end up in the callback

  const std::string& lFilename = Configuration::Instance.outputFile();

  if( lFilename.size() == 0 )
  {
    std::cout << "Warning: Running scan without callback" << std::endl;
    lEvent.ScanRT( [&]( const EventProxy& aEvent , const double& aR , const double& aT, std::pair<int, int> aCurrentIJ){} ); // Null callback
  }
  else if( lFilename.size() > 4 and lFilename.substr(lFilename.size() - 4) == ".xml" )
  {
    std::stringstream lOutput;
    lEvent.ScanRT( [&]( const EventProxy& aEvent , const double& aR , const double& aT, std::pair<int, int> aCurrentIJ ){ XmlCallback( aEvent , aR , aT, aCurrentIJ  , lOutput, lRTScores, lMaxScorePosition, lMaxRTScore); } );
    std::ofstream lOutFile( lFilename );
    lOutFile << "<Results>\n" << lOutput.str() << "</Results>\n";
  }
  else if( lFilename.size() > 5 and lFilename.substr(lFilename.size() - 5) == ".json" )
  {
    std::stringstream lOutput;
    lEvent.ScanRT( [&]( const EventProxy& aEvent , const double& aR , const double& aT, std::pair<int, int> aCurrentIJ ){ JsonCallback( aEvent , aR , aT, aCurrentIJ  , lOutput, lRTScores, lMaxScorePosition, lMaxRTScore); } );
    std::ofstream lOutFile( lFilename );
    lOutFile << "{\nResults:[\n" << lOutput.str() << "]\n}";
  }
  else
  {
    throw std::runtime_error( "No handler for specified output-file" );
  }

  std::cout << "max score was: " << lMaxRTScore << std::endl;
  std::cout << "at position (" << lMaxScorePosition.first << ", " << lMaxScorePosition.second << ")"<< std::endl;
  std::cout << "out of a possible " << Configuration::Instance.Rbins() << " R Bins"
  << " and " << Configuration::Instance.Tbins() << " T Bins" << std::endl;
  std::pair<double,double> a;
  a = bestRT(lMaxScorePosition, lRTScores);
  std::cout << "best R value is: " << a.first << " and the best T value is: " << a.second << std::endl;

}
