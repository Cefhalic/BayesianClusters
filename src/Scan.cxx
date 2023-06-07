//! \file Scan.cxx

/* ===== Cluster sources ===== */
#include "BayesianClustering/API.hpp"
#include "BayesianClustering/LocalizationFile.hpp"
// #include "BayesianClustering/Cluster.hpp"
#include "BayesianClustering/RoI.hpp"
// #include "BayesianClustering/RoIproxy.hpp"
// #include "BayesianClustering/Configuration.hpp"

/* ===== C++ ===== */
#include <vector>
#include <iostream>
#include <stdio.h>

/* ===== Local utilities ===== */
#include "Utilities/ProgressBar.hpp"

// /* ===== GSL libraries ===== */
// #include <gsl/gsl_math.h>
// #include <gsl/gsl_interp2d.h>
// #include <gsl/gsl_spline2d.h>

// //! mutex for critical section
// std::mutex mtx;


// //! A callback for writing the cluster info to XML file
// //! aRoI              The current RoIproxy
// //! aR                The current R position of the scan
// //! aT                The current T position of the scan
// //! aCurrentIJ        The current position of the scan in integer units
// //! aOutput           The output stream
// //! aRTscores         A 2D map of scores to fill
// //! aMaxScorePosition The position of the maximum score
// //! aMaxRTScore       The maximum score
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


// //! A callback for writing the cluster info to JSON file
// //! \param aRoI              The current RoIproxy
// //! \param aR                The current R position of the scan
// //! \param aT                The current T position of the scan
// //! \param aCurrentIJ        The current position of the scan in integer units
// //! \param aOutput           The output stream
// //! \param aRTScores         A 2D map of scores to fill
// //! \param aMaxScorePosition The position of the maximum score
// //! \param aMaxRTScore       The maximum score
// void JsonCallback( const RoIproxy& aRoI, const double& aR, const double& aT, std::pair<int,int>& aCurrentIJ, std::stringstream& aOutput,
//                    std::vector<std::vector<double>>& aRTScores, std::pair<int,int>& aMaxScorePosition, double& aMaxRTScore )
// {
//   mtx.lock();


//   // std::cout << "R " << aCurrentIJ.first  << " " << CurrentConfiguration().minScanR() + (aCurrentIJ.first  * CurrentConfiguration().dR()) << " " << aR << std::endl;
//   // std::cout << "T " << aCurrentIJ.second << " " << CurrentConfiguration().maxScanT() - (aCurrentIJ.second * CurrentConfiguration().dT()) << " " << aT << std::endl;

//   // aOutput << "  { R:" << aR << ", T:" << aT << ", Score:" << aRoI.mLogP << ", NumClusteredPts:" << aRoI.mClusteredCount << ", NumBackgroundPts:" << aRoI.mBackgroundCount << ", Clusters:[\n";
//   aOutput << "  { R:" << aR << ", T:" << aT << ", Score:" << aRoI.mLogP << ", NumClusteredPts:" << aRoI.mClusteredCount << ", NumBackgroundPts:" << aRoI.mBackgroundCount << "}\n";

//   // for( auto& i : aRoI.mClusters )
//   // {
//   //   if( i.mClusterSize ) aOutput << "    { Points:" << i.mClusterSize << ",  Score:" << i.mClusterScore << " },\n";
//   // }

//   double lLogP = aRoI.mLogP;
//   //score setting stuff
//   if (lLogP > aMaxRTScore) {
//     aMaxScorePosition = aCurrentIJ;
//     aMaxRTScore = lLogP;
//   }

//   uint32_t p = aCurrentIJ.first, q = aCurrentIJ.second;
//   aRTScores[p][q] = lLogP;

//   // aOutput << "   },\n";
//   mtx.unlock();
// }

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


// //! A callback for handling each RoI
// //! \param aRoI        The current RoI
// //! \param aOutFile    The name of the output file
// //! \param aScanConfig The configuration parameters for the scan
// void RoIcallback( RoI& aRoI, const std::string& aOutFile, const ScanConfiguration& aScanConfig )
// {
//   std::cout << "Scanning RoI with " << aRoI.mData.size() << " localizations" << std::endl;

//   std::vector<std::vector<double>> lRTScores( aScanConfig.Rbounds().bins, std::vector<double>( aScanConfig.Tbounds().bins /*, 1*/));
//   std::pair<int, int> lMaxScorePosition;
//   double lMaxRTScore = -9E99;
//   // the above will store our scores - it needs to end up in the callback

//   if( aOutFile.size() == 0 ) {
//     std::cout << "Warning: Running scan without callback" << std::endl;
//     aRoI.ScanRT( aScanConfig, [&]( const RoIproxy& aRoI, const double& aR, const double& aT, std::pair<int, int> aCurrentIJ) {} );   // Null callback
//     return;
//   }
//   // else if( lOutputFilename.size() > 4 and lOutputFilename.substr(lOutputFilename.size() - 4) == ".xml" )
//   // {
//   //   std::stringstream lOutput;
//   //   aRoI.ScanRT( lAuxConfig.Rbounds() , lAuxConfig.Tbounds() , [&]( const RoIproxy& aRoI , const double& aR , const double& aT, std::pair<int, int> aCurrentIJ ){ XmlCallback( aRoI , aR , aT, aCurrentIJ  , lOutput, lRTScores, lMaxScorePosition, lMaxRTScore); } );
//   //   std::ofstream lOutFile( aOutFile );
//   //   lOutFile << "<Results>\n" << lOutput.str() << "</Results>\n";
//   // }
//   else if( aOutFile.size() > 5 and aOutFile.substr(aOutFile.size() - 5) == ".json" ) {
//     std::stringstream lOutput;
//     aRoI.ScanRT( aScanConfig, [&]( const RoIproxy& aRoI, const double& aR, const double& aT, std::pair<int, int> aCurrentIJ ) {
//       JsonCallback( aRoI, aR, aT, aCurrentIJ, lOutput, lRTScores, lMaxScorePosition, lMaxRTScore);
//     } );
//     std::ofstream lOutFile( aOutFile );
//     lOutFile << "{\nResults:[\n" << lOutput.str() << "]\n}";
//   } else {
//     throw std::runtime_error( "No handler for specified output-file" );
//   }

//   std::cout << "max score was: " << lMaxRTScore << std::endl;
//   std::cout << "at position (" << lMaxScorePosition.first << ", " << lMaxScorePosition.second << ")"<< std::endl;
//   std::cout << "out of a possible " << aScanConfig.Rbounds().bins << " R Bins"
//             << " and " << aScanConfig.Tbounds().bins << " T Bins" << std::endl;
//   // std::pair<double,double> a;
//   // a = bestRT(lMaxScorePosition, lRTScores);
//   // std::cout << "best R value is: " << a.first << " and the best T value is: " << a.second << std::endl;

// }


void ScanCallback_Json( const std::vector< ScanEntry >& aVector, const std::string& aOutFile )
{
  // ===========================================================================================================
  char lOutFileName[256];
  static int RoIid( 0 );  
  sprintf( lOutFileName , "RoI%d.%s", RoIid++ , aOutFile.c_str() );
  FILE *fptr = fopen( lOutFileName , "w" );
  if (fptr == NULL) throw std::runtime_error("Could not open file");
  fprintf( fptr , "[\n" );
  for( auto& lIt : aVector ) fprintf( fptr , "  { r : %.5e , t : %.5e , logP : %+.5e },\n" , lIt.r , lIt.t , lIt.score );
  fprintf( fptr , "]\n" );
  fclose(fptr); 
  // ===========================================================================================================

  // ===========================================================================================================
  // std::vector< double > x , y , z;
  // ScanEntry lMax{ 0 , 0 , -9e99 };

  // for( auto& lIt : aVector )
  // {
  //     if( lIt.score > lMax.score ) lMax = lIt;
  //     z.push_back( lIt.score );    
  // }

  // for( std::size_t i(0) ; i!=aScanConfig.Tbounds().bins ; i+=1                          ) x.push_back( aVector.at(i).t );
  // for( std::size_t i(0) ; i!=aVector.size()             ; i+=aScanConfig.Tbounds().bins ) y.push_back( aVector.at(i).r );

  // // for( auto& i : x ) std::cout << "t" << i << std::endl;
  // // for( auto& i : y ) std::cout << "r" << i << std::endl;

  // const gsl_interp2d_type *T = gsl_interp2d_bicubic;
  // gsl_spline2d *spline = gsl_spline2d_alloc( T , aScanConfig.Rbounds().bins , aScanConfig.Tbounds().bins );
  // gsl_interp_accel *xacc = gsl_interp_accel_alloc();
  // gsl_interp_accel *yacc = gsl_interp_accel_alloc();

  // gsl_spline2d_init( spline, &x.front() , &y.front() , &z.front() , aScanConfig.Rbounds().bins , aScanConfig.Tbounds().bins );

  // for( int i(0) ; i!=10 ; ++i )
  // {
  //   auto dx = gsl_spline2d_eval_deriv_x( spline , lMax.r , lMax.t , xacc, yacc );
  //   auto dy = gsl_spline2d_eval_deriv_y( spline , lMax.r , lMax.t , xacc, yacc );
  //   auto val = gsl_spline2d_eval( spline , lMax.r , lMax.t , xacc, yacc );
  //   printf( "Fit : %e %e : %e %e %e \n", lMax.r , lMax.t , dx , dy , val );
  //   lMax.r += (1e-27 * dx);
  //   lMax.t += (1e-27 * dy);
  // }
  // // printf("1.0 and -1.5 : %f\n", gsl_spline2d_eval(spline,1.0,-1.5,xacc, yacc));
  // ===========================================================================================================
}

//! The main function
//! \param argc The number of commandline arguments
//! \param argv The commandline arguments
//! \return     The exit code
int main(int argc, char** argv)
{
  std::cout << "+------------------------------------+" << std::endl;
  ProgressBar2 lBar( "| Cluster Scan. Andrew W. Rose. 2022 |", 1 );
  std::cout << "+------------------------------------+" << std::endl;
  AuxConfiguration lAuxConfig( argc, argv );
  std::cout << "+------------------------------------+" << std::endl;

  AutoRoi_Scan_SimpleCallback( lAuxConfig.inputFile() , 
                               lAuxConfig.configFile() , 
                               [&]( const std::vector< ScanEntry >& aVector ){ ScanCallback_Json( aVector , lAuxConfig.outputFile() ); }
                               );
}
