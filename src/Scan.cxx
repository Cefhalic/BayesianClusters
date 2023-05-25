/* ===== Cluster sources ===== */
#include "BayesianClustering/API.hpp"
#include "BayesianClustering/Cluster.hpp"
#include "BayesianClustering/RoI.hpp"
#include "BayesianClustering/RoIproxy.hpp"
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



// void XmlCallback( const RoIproxy& aRoI , const double& aR , const double& aT, std::pair<int, int>& aCurrentIJ , std::stringstream& aOutput, 
//                   std::vector<std::vector<double>>& aRTScores, std::pair<int,int>& aMaxScorePosition, double& aMaxRTScore )
// {
//   mtx.lock();
//   // aOutput << "  { R:" << aR << ", T:" << aT << ", Score:" << aRoI.mLogP << ", NumClusteredPts:" << aRoI.mClusteredCount << ", NumBackgroundPts:" << aRoI.mBackgroundCount << ", Clusters:[\n";
//   aOutput << "  { R:" << aR << ", T:" << aT << ", Score:" << aRoI.mLogP << ", NumClusteredPts:" << aRoI.mClusteredCount << ", NumBackgroundPts:" << aRoI.mBackgroundCount << "}\n";

//   // for( auto& i : aRoI.mClusters )
//   // {
//   //   if( i.mClusterSize ) aOutput << "    { Points:" << i.mClusterSize << ",  Score:" << i.mClusterScore << " },\n";
//   // }
  
//   double lLogP = aRoI.mLogP;
//   //score setting stuff
//   if (lLogP > aMaxRTScore){
//     aMaxScorePosition = aCurrentIJ;
//     aMaxRTScore = lLogP;
//   }

//   uint32_t p = aCurrentIJ.first, q = aCurrentIJ.second;
//   aRTScores[p][q] = lLogP;

//   // aOutput << "   },\n";
//   mtx.unlock();  
// }


void JsonCallback( const RoIproxy& aRoI , const double& aR , const double& aT, std::pair<int,int>& aCurrentIJ , std::stringstream& aOutput, 
                  std::vector<std::vector<double>>& aRTScores, std::pair<int,int>& aMaxScorePosition, double& aMaxRTScore )
{
  mtx.lock();


  std::cout << "R " << aCurrentIJ.first  << " " << CurrentConfiguration().minScanR() + (aCurrentIJ.first  * CurrentConfiguration().dR()) << " " << aR << std::endl;
  std::cout << "T " << aCurrentIJ.second << " " << CurrentConfiguration().maxScanT() - (aCurrentIJ.second * CurrentConfiguration().dT()) << " " << aT << std::endl;

  // aOutput << "  { R:" << aR << ", T:" << aT << ", Score:" << aRoI.mLogP << ", NumClusteredPts:" << aRoI.mClusteredCount << ", NumBackgroundPts:" << aRoI.mBackgroundCount << ", Clusters:[\n";
  aOutput << "  { R:" << aR << ", T:" << aT << ", Score:" << aRoI.mLogP << ", NumClusteredPts:" << aRoI.mClusteredCount << ", NumBackgroundPts:" << aRoI.mBackgroundCount << "}\n";

  // for( auto& i : aRoI.mClusters )
  // {
  //   if( i.mClusterSize ) aOutput << "    { Points:" << i.mClusterSize << ",  Score:" << i.mClusterScore << " },\n";
  // }
  
  double lLogP = aRoI.mLogP;
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

// //just prints the best R, T at the end
// std::pair<double,double> bestRT(std::pair<int, int>& aMaxScorePosition, std::vector<std::vector<double>>& aRTScores){
//   int i = aMaxScorePosition.first, j = aMaxScorePosition.second;
//   double lRValue(0), lTValue(0);
//   double lValueSum(0);
//   for(int I(-2); I < 3 ; ++I ){
//     if ((i < 2) or (j < 2)) break;
//     if ((i > CurrentConfiguration().Rbins() - 3 ) or (j > CurrentConfiguration().Tbins() - 3)) break;
//     for(int J(-2); J < 3 ; ++J ){
//       auto lVal = aRTScores[i + I][j + J];
//       lValueSum += lVal;
//       lRValue += (i+I) * lVal;
//       lTValue += (j+J) * lVal;
//   }}

//   double lRIndex = lRValue / lValueSum;
//   double lTIndex = lTValue /  lValueSum;

//   double outputR = CurrentConfiguration().minScanR() + (lRIndex * CurrentConfiguration().dR());
//   double outputT = CurrentConfiguration().maxScanT() - (lTIndex * CurrentConfiguration().dT());

//   return std::make_pair(CurrentConfiguration().toPhysicalUnits(outputR), CurrentConfiguration().toPhysicalUnits(outputT));
// }



/* ===== Main function ===== */
int main(int argc, char **argv)
{
  std::cout << "+------------------------------------+" << std::endl;
  ProgressBar2 lBar( "| Cluster Scan. Andrew W. Rose. 2022 |" , 1 );
  std::cout << "+------------------------------------+" << std::endl;
  Configuration lMasterConfig;  
  lMasterConfig.FromCommandline( argc , argv );
  std::cout << "+------------------------------------+" << std::endl;

  const std::string& lInputFilename = lMasterConfig.inputFile();
  if( lInputFilename.size() == 0 ) throw std::runtime_error( "No input file specified" ); 
  auto lDataset = LoadLocalizationFile( lInputFilename );

  SetCurrentConfiguration( lMasterConfig );
  for( auto& lRoI : ExtractRoIs( lDataset , FromConfigFile ) )    
  {
    SetCurrentConfiguration( lRoI.mConfiguration );

    std::vector<std::vector<double>> lRTScores( lMasterConfig.Rbins() , std::vector<double>( lMasterConfig.Tbins() /*, 1*/));
    std::pair<int, int> lMaxScorePosition;
    double lMaxRTScore = -9E99;
    //the above will store our scores - it needs to end up in the callback

    const std::string& lOutputFilename = lMasterConfig.outputFile();

    if( lOutputFilename.size() == 0 )
    {
      std::cout << "Warning: Running scan without callback" << std::endl;
      lRoI.ScanRT( lMasterConfig.Rbounds() , lMasterConfig.Tbounds() , [&]( const RoIproxy& aRoI , const double& aR , const double& aT, std::pair<int, int> aCurrentIJ){} ); // Null callback
      continue;
    }
    // else if( lOutputFilename.size() > 4 and lOutputFilename.substr(lOutputFilename.size() - 4) == ".xml" )
    // {
    //   std::stringstream lOutput;
    //   lRoI.ScanRT( lMasterConfig.Rbounds() , lMasterConfig.Tbounds() , [&]( const RoIproxy& aRoI , const double& aR , const double& aT, std::pair<int, int> aCurrentIJ ){ XmlCallback( aRoI , aR , aT, aCurrentIJ  , lOutput, lRTScores, lMaxScorePosition, lMaxRTScore); } );
    //   std::ofstream lOutFile( lOutputFilename );
    //   lOutFile << "<Results>\n" << lOutput.str() << "</Results>\n";
    // }
    else if( lOutputFilename.size() > 5 and lOutputFilename.substr(lOutputFilename.size() - 5) == ".json" )
    {
      std::stringstream lOutput;
      lRoI.ScanRT( lMasterConfig.Rbounds() , lMasterConfig.Tbounds() , [&]( const RoIproxy& aRoI , const double& aR , const double& aT, std::pair<int, int> aCurrentIJ ){ JsonCallback( aRoI , aR , aT, aCurrentIJ  , lOutput, lRTScores, lMaxScorePosition, lMaxRTScore); } );
      std::ofstream lOutFile( lOutputFilename );
      lOutFile << "{\nResults:[\n" << lOutput.str() << "]\n}";
    }
    else
    {
      throw std::runtime_error( "No handler for specified output-file" );
    }

    // std::cout << "max score was: " << lMaxRTScore << std::endl;
    // std::cout << "at position (" << lMaxScorePosition.first << ", " << lMaxScorePosition.second << ")"<< std::endl;
    // std::cout << "out of a possible " << CurrentConfiguration().Rbins() << " R Bins"
    // << " and " << CurrentConfiguration().Tbins() << " T Bins" << std::endl;
    // std::pair<double,double> a;
    // a = bestRT(lMaxScorePosition, lRTScores);
    // std::cout << "best R value is: " << a.first << " and the best T value is: " << a.second << std::endl;

  }
}
