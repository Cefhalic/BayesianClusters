//! \file Configuration.cpp

/* ===== Local utilities ===== */
#include "Utilities/GSLInterpolator.hpp"
#include "Utilities/ListComprehension.hpp"
#include "Utilities/Vectorize.hpp"
#include "Utilities/Units.hpp"
#include "BayesianClustering/Configuration.hpp"

/* ===== C++ ===== */
#include <iostream>
#include <fstream>
#include <streambuf>

/* ===== BOOST libraries ===== */
#include <boost/math/special_functions/gamma.hpp>
#include "boost/algorithm/string.hpp"
#include "boost/program_options.hpp"
namespace po = boost::program_options;


// ScanConfiguration::ScanConfiguration( int argc, char** argv ) :
//   mSigmacount(-1), mSigmaspacing(-1),
//   mRbounds{-1,-1,-1,UINT_MAX}, mTbounds{-1,-1,-1,UINT_MAX},
//   mAlpha(-1), mLogAlpha(-1), mLogGammaAlpha(-1),
//   mLogPb(-1), mLogPbDagger(-1)
// {
//   std::vector< std::string > lTemp( argv+1, argv+argc );
//   FromVector( lTemp );
// }

// ScanConfiguration::ScanConfiguration( const std::vector< std::string >& aArgs ) :
//   mSigmacount(-1), mSigmaspacing(-1),
//   mRbounds{-1,-1,-1,UINT_MAX}, mTbounds{-1,-1,-1,UINT_MAX},
//   mAlpha(-1), mLogAlpha(-1), mLogGammaAlpha(-1),
//   mLogPb(-1), mLogPbDagger(-1)
// {
//   FromVector( aArgs );
// }

ScanConfiguration::ScanConfiguration( const std::string& aCfgFile ) :
  mSigmacount(-1), mSigmaspacing(-1),
  mRbounds{-1,-1,-1,UINT_MAX}, mTbounds{-1,-1,-1,UINT_MAX},
  mAlpha(-1), mLogAlpha(-1), mLogGammaAlpha(-1),
  mLogPb(-1), mLogPbDagger(-1)
{
  std::ifstream lFstr( aCfgFile.c_str() );
  std::string lStr( (std::istreambuf_iterator<char>(lFstr)), std::istreambuf_iterator<char>() );

  std::vector<std::string> lStrs( 1, "Config-file" );
  boost::split( lStrs, lStr, boost::is_any_of( " \t\r\n" ) );
  FromVector( lStrs );
}

ScanConfiguration::ScanConfiguration( 
  const std::size_t& aSigmacount, const double& aSigmaMin, const double& aSigmaMax, const std::function< double( const double& ) >& aInterpolator ,
  const std::size_t& aRbins, const double& aMinScanR, const double& aMaxScanR ,
  const std::size_t& aTbins, const double& aMinScanT, const double& aMaxScanT ,
  const double& aPB ,
  const double& aAlpha )
{
  SetSigmaParameters( aSigmacount, aSigmaMin, aSigmaMax, aInterpolator );
  SetRBins( aRbins, aMinScanR, aMaxScanR );
  SetTBins( aTbins, aMinScanT, aMaxScanT );
  SetPb( aPB );
  SetAlpha( aAlpha );
}

void ScanConfiguration::SetSigmaParameters( const std::size_t& aSigmacount, const double& aSigmaMin, const double& aSigmaMax, const std::function< double( const double& ) >& aInterpolator )
{
  std::cout << "Sigma-integral: " << aSigmaMin << " to " << aSigmaMax << " in " << aSigmacount << " steps" << std::endl;

  mSigmacount = aSigmacount;

  if( aSigmacount == 0 ) {
    mSigmabins.clear();
    mSigmabins2.clear();
    mProbabilitySigma.clear();
    mLogProbabilitySigma.clear();
    mSigmaspacing = 0;
    return;
  }

  mSigmaspacing = ( aSigmaMax - aSigmaMin ) / aSigmacount;
  mSigmabins = [ & ]( const int& i ) {
    return ( i * mSigmaspacing ) + aSigmaMin;
  } | range( mSigmacount );
  mSigmabins2 = []( const double& i ) {
    return i * i;
  } | mSigmabins;
  mProbabilitySigma = aInterpolator | mSigmabins;
  mLogProbabilitySigma = []( const double& w) {
    return log(w);
  } | mProbabilitySigma;
}

void ScanConfiguration::SetRBins( const std::size_t& aRbins, const double& aMinScanR, const double& aMaxScanR )
{
  mRbounds.bins = aRbins;
  mRbounds.min = aMinScanR ;
  mRbounds.max = aMaxScanR ;
  mRbounds.spacing = ( mRbounds.max - mRbounds.min ) / mRbounds.bins;

  std::cout << "R-bins: " << aRbins << " bins from " << aMinScanR << " to " << aMaxScanR << " in steps of " << mRbounds.spacing << std::endl;
}

void ScanConfiguration::SetTBins( const std::size_t& aTbins, const double& aMinScanT, const double& aMaxScanT )
{
  mTbounds.bins = aTbins;
  mTbounds.min = aMinScanT ;
  mTbounds.max = aMaxScanT ;
  mTbounds.spacing = ( mTbounds.max - mTbounds.min ) / mTbounds.bins;

  std::cout << "T-bins: " << aTbins << " bins from " << aMinScanT << " to " << aMaxScanT << " in steps of " << mTbounds.spacing << std::endl;
}

void ScanConfiguration::SetPb( const double& aPB )
{
  std::cout << "Pb: " << aPB << std::endl;
  mLogPb = log( aPB );
  mLogPbDagger = log( 1-aPB );
}

void ScanConfiguration::SetAlpha( const double& aAlpha )
{
  std::cout << "Alpha: " << aAlpha << std::endl;
  mAlpha = aAlpha;
  mLogAlpha = log( aAlpha );
  mLogGammaAlpha = boost::math::lgamma( aAlpha );
}







// //! Utility function for parsing a config file containing "commandline" parameters
// //! \param aDesc A boost program-options description for parsing the values
// //! \param aFilename The name of the config file
// void config_file( const po::options_description& aDesc, const std::string& aFilename )
// {
//   std::ifstream lFstr( aFilename.c_str() );
//   std::string lStr( (std::istreambuf_iterator<char>(lFstr)), std::istreambuf_iterator<char>() );

//   std::vector<std::string> lStrs( 1, "Config-file" );
//   boost::split( lStrs, lStr, boost::is_any_of( " \t\r\n" ) );

//   po::variables_map lVm;
//   po::store( po::command_line_parser( lStrs ).options( aDesc ).allow_unregistered().run(), lVm );
//   po::notify( lVm );
// }



// void ScanConfiguration::FromCommandline( int argc, char** argv )
// {
//   std::vector< std::string > lTemp( argv+1, argv+argc );
//   FromVector( lTemp );
// }

void ScanConfiguration::FromVector( const std::vector< std::string >& aArgs )
{
  typedef std::string tS;
  typedef std::vector<std::string> tVS;
  typedef unsigned tU;
  typedef double tD;
  typedef std::vector<double> tVD;

  tD sigLo, sigHi, rLo, rHi, tLo, tHi;
  tU Nsig(0), Nr(0), Nt(0);
  tVD SigKeys, SigVals;

  po::options_description lDesc("General options");
  lDesc.add_options()
  ( "help",         po::bool_switch()                          ->notifier( [&]( const bool& aArg ) { if( aArg ) { std::cout << lDesc << std::endl; exit(0); } } ), "produce help message" )  
  ( "sigma-bins",   po::value<tU>(&Nsig), "Number of sigma bins" )  
  ( "sigma-low",    po::value<tS>()                            ->notifier( [&]( const   tS& aArg ) { sigLo=StrToDist(aArg); } ), "Lower sigma integration bound" )  
  ( "sigma-high",   po::value<tS>()                            ->notifier( [&]( const   tS& aArg ) { sigHi=StrToDist(aArg); } ), "High sigma integration bound" )  
  ( "sigma-curve",  po::value<tVS>()->composing()->multitoken()->notifier( [&]( const  tVS& aArg ) {  for( auto& i : aArg ) { 
                                                                                                        std::vector<std::string> lStrs;      
                                                                                                        boost::split( lStrs, i, [](char c) { return c==':'; } );      
                                                                                                        SigKeys.push_back( StrToDist( lStrs.at(0) ) );      
                                                                                                        SigVals.push_back( std::stod( lStrs.at(1) ) );
                                                                                                      }  
                                                                                                    } ), "Parameterized sigma probability curve (list of colon-separated size-probability pairs)" )  
  ( "r-bins",       po::value<tU>(&Nr), "Number of R bins" )  
  ( "r-low",        po::value<tS>()                             ->notifier( [&]( const   tS& aArg ) { rLo=StrToDist(aArg); } ), "Lower R bound" )  
  ( "r-high",       po::value<tS>()                             ->notifier( [&]( const   tS& aArg ) { rHi=StrToDist(aArg); } ), "High R bound" )  
  ( "t-bins",       po::value<tU>(&Nt), "Number of T bins" )  
  ( "t-low",        po::value<tS>()                             ->notifier( [&]( const   tS& aArg ) { tLo=StrToDist(aArg); } ), "Lower T bound" )  
  ( "t-high",       po::value<tS>()                             ->notifier( [&]( const   tS& aArg ) { tHi=StrToDist(aArg); } ), "High T bound" )  
  ( "pb",           po::value<tD>()                             ->notifier( [&]( const   tD& aArg ) { SetPb(aArg); } ), "pb parameter" )  
  ( "alpha",        po::value<tD>()                             ->notifier( [&]( const   tD& aArg ) { SetAlpha(aArg); } ), "alpha parameter" )
  ;

  po::variables_map lVm;
  po::store( po::command_line_parser( aArgs ).options( lDesc ).allow_unregistered().run(), lVm );
  po::notify( lVm );

  if( Nr ) SetRBins( Nr, rLo, rHi );
  if( Nt ) SetTBins( Nt, tLo, tHi );

  if( Nsig ) {
    GSLInterpolator lInterpolator( gsl_interp_cspline, SigKeys, SigVals );
    SetSigmaParameters( Nsig, sigLo, sigHi, [&]( const double& aPt ) { return lInterpolator.Eval( aPt ); } );
  }

}




AuxConfiguration::AuxConfiguration( int argc, char** argv ) :
  mValidate(false),
  mInputFile(""), mOutputFile(""),
  mClusterR( -1 ), mClusterT(-1)
{
  std::vector< std::string > lTemp( argv+1, argv+argc );
  FromVector( lTemp );
}

AuxConfiguration::AuxConfiguration( const std::vector< std::string >& aArgs ) :
  mValidate(false),
  mInputFile(""), mOutputFile(""),
  mClusterR( -1 ), mClusterT(-1)
{
  FromVector( aArgs );
}

void AuxConfiguration::SetValidate( const bool& aValidate )
{
  if( aValidate ) std::cout << "Validate: TRUE" << std::endl;
  mValidate = aValidate;
}

void AuxConfiguration::SetInputFile( const std::string& aFileName )
{
  std::cout << "Input file: " << aFileName << std::endl;
  mInputFile = aFileName;
}

void AuxConfiguration::SetOutputFile( const std::string& aFileName )
{
  std::cout << "Output file: " << aFileName << std::endl;
  mOutputFile = aFileName;
}


void AuxConfiguration::SetConfigFile( const std::string& aFileName )
{
  std::cout << "Config file: " << aFileName << std::endl;
  mConfigFile = aFileName;
}

// void AuxConfiguration::FromCommandline( int argc, char** argv )
// {
//   std::vector< std::string > lTemp( argv+1, argv+argc );
//   FromVector( lTemp );
// }

void AuxConfiguration::FromVector( const std::vector< std::string >& aArgs )
{
  typedef std::string tS;
  typedef std::vector<double> tVD;
  typedef std::size_t tZ;

  tVD SigKeys, SigVals;

  po::options_description lDesc("General options");
  lDesc.add_options()  
  ( "help",         po::bool_switch()                          ->notifier( [&]( const bool& aArg ) {    if( aArg ) {      std::cout << lDesc << std::endl;      exit(0);    }  } ), "produce help message" )  
  // ( "cfg",          po::value<tS>()                            ->notifier( [&]( const   tS& aArg ) {    config_file( lDesc, aArg );  } ), "Config file" )    
  // ( "centre",       po::value<tVS>()->composing()->multitoken()->notifier( [&]( const  tVS& aArg ){ SetCentre( StrToDist(aArg.at(0)) , StrToDist(aArg.at(1)) ); } ) , "Centre of ROI as 'x y' pair" )  
  // ( "width",        po::value<tVS>()->composing()->multitoken()->notifier( [&]( const  tVS& aArg ){ SetWidth( StrToDist(aArg.at(0)) , StrToDist(aArg.at(1)) ); } ) , "Width of ROI as 'x y' pair" )  
  ( "validate,v",   po::bool_switch()                           ->notifier( [&]( const bool& aArg ) {    SetValidate( aArg );  } ), "validate clusters" )  
  ( "input-file,i", po::value<tS>()                             ->notifier( [&]( const   tS& aArg ) {    SetInputFile(aArg);  } ), "input file")  
  ( "output-file,o", po::value<tS>()                            ->notifier( [&]( const   tS& aArg ) {    SetOutputFile(aArg);  } ), "output file")  
  ( "cfg",          po::value<tS>()                             ->notifier( [&]( const   tS& aArg ) {    SetConfigFile(aArg);  } ), "Config file")  
  ( "r",            po::value<tS>()                             ->notifier( [&]( const   tS& aArg ) {    mClusterR=StrToDist(aArg);  } ), "R for clustering" )  
  ( "t",            po::value<tS>()                             ->notifier( [&]( const   tS& aArg ) {    mClusterT=StrToDist(aArg);  } ), "T for clustering" )  
  ( "threads",      po::value<tZ>( &Nthreads ), "Number of threads to use (default is value given by std::threads::hardware_concurrency())" )
  ;

  po::variables_map lVm;
  po::store( po::command_line_parser( aArgs ).options( lDesc ).allow_unregistered() /*.positional( lPositional )*/ .run(), lVm );
  po::notify( lVm );
}
