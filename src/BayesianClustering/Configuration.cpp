

/* ===== Local utilities ===== */
#include "Utilities/GSLInterpolator.hpp"
#include "Utilities/ListComprehension.hpp"
#include "Utilities/Vectorize.hpp"
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

Configuration::Configuration() :
	// mScale(1) , mScale2(1),
  mPhysicalCentreX(0), mPhysicalCentreY(0),
  mWidthX(0), mWidthY(0),
  mArea(0),
	mSigmacount(-1), mSigmaspacing(-1),
	mMaxR(-1), mMaxR2(-1), mMax2R(-1), mMax2R2(-1),
	mMinScanR(-1), mMaxScanR(-1), mMinScanT(-1), mMaxScanT(-1),
	mDR(-1), mDT(-1),
	mRbins(-1),  mTbins(-1),
	mLogPb(-1), mLogPbDagger(-1), 
	mAlpha(-1), mLogAlpha(-1), mLogGammaAlpha(-1),
	mValidate(false),
  mInputFile(""), mOutputFile(""),
  mClusterR( -1 ), mClusterT(-1)
{}


void Configuration::SetCentre( const double& aPhysicalCentreX , const double& aPhysicalCentreY )
{
	std::cout << "Centre: x=" << aPhysicalCentreX << ", y=" << aPhysicalCentreY << std::endl;
  mPhysicalCentreX = aPhysicalCentreX;
  mPhysicalCentreY = aPhysicalCentreY;
}

void Configuration::SetWidth( const double& aWidthX , const double& aWidthY )
{
  std::cout << "Centre: x=" << aWidthX << ", y=" << aWidthY << std::endl;
  mWidthX = aWidthX;
  mWidthY = aWidthY;
  mArea = mWidthX * mWidthY;
}

void Configuration::SetSigmaParameters( const std::size_t& aSigmacount , const double& aSigmaMin , const double& aSigmaMax , const std::function< double( const double& ) >& aInterpolator )
{
	std::cout << "Sigma-integral: " << aSigmaMin << " to " << aSigmaMax << " in " << aSigmacount << " steps" << std::endl;

	mSigmacount = aSigmacount;
	mSigmaspacing = ( aSigmaMax - aSigmaMin ) / aSigmacount;
	mSigmabins = [ & ]( const int& i ){ return ( i * mSigmaspacing ) + aSigmaMin;  } | range( mSigmacount );
	mSigmabins2 = []( const double& i ){ return i * i; } | mSigmabins;
	mProbabilitySigma = aInterpolator | mSigmabins;
	mLogProbabilitySigma = []( const double& w){ return log(w); } | mProbabilitySigma;
}

void Configuration::SetRBins( const std::size_t& aRbins , const double& aMinScanR , const double& aMaxScanR )
{
	mRbins = aRbins;
	mMinScanR = aMinScanR;
	mMaxScanR = aMaxScanR;
	mDR = ( mMaxScanR - mMinScanR ) / mRbins;

	std::cout << "R-bins: " << mRbins << " bins from " << mMinScanR << " to " << mMaxScanR << " in steps of " << mDR << std::endl;

	mMaxR = aMaxScanR;
	mMaxR2 = mMaxR * mMaxR;
	mMax2R = 2.0 * mMaxR;
	mMax2R2 = mMax2R * mMax2R;	
}

void Configuration::SetTBins( const std::size_t& aTbins , const double& aMinScanT , const double& aMaxScanT )
{
	mTbins = aTbins;
	mMinScanT = aMinScanT;
	mMaxScanT = aMaxScanT;
	mDT = ( mMaxScanT - mMinScanT ) / mTbins;

	std::cout << "T-bins: " << mTbins << " bins from " << mMinScanT << " to " << mMaxScanT << " in steps of " << mDT << std::endl;
}

void Configuration::SetPb( const double& aPB )
{
	std::cout << "Pb: " << aPB << std::endl;
	mLogPb = log( aPB );
	mLogPbDagger = log( 1-aPB );
}

void Configuration::SetAlpha( const double& aAlpha )
{
	std::cout << "Alpha: " << aAlpha << std::endl;
	mAlpha = aAlpha;
	mLogAlpha = log( aAlpha );
	mLogGammaAlpha = boost::math::lgamma( aAlpha );
}

void Configuration::SetValidate( const bool& aValidate )
{
	if( aValidate ) std::cout << "Validate: TRUE" << std::endl;

	mValidate = aValidate;
}


void Configuration::SetInputFile( const std::string& aFileName )
{ 
  std::cout << "Input file: " << aFileName << std::endl;

  mInputFile = aFileName;
}

void Configuration::SetOutputFile( const std::string& aFileName )
{ 
  std::cout << "Output file: " << aFileName << std::endl;

  mOutputFile = aFileName;
}



void config_file( const po::options_description& aDesc , const std::string& aFilename )
{
  std::ifstream lFstr( aFilename.c_str() );
  std::string lStr( (std::istreambuf_iterator<char>(lFstr)) , std::istreambuf_iterator<char>() );

  std::vector<std::string> lStrs( 1 , "Config-file" );
  boost::split( lStrs , lStr , boost::is_any_of( " \t\r\n" ) );

  po::variables_map lVm;
  po::store( po::command_line_parser( lStrs ).options( aDesc ).run() , lVm );
  po::notify( lVm );
}



void Configuration::FromCommandline( int argc , char **argv )
{
  std::vector< std::string > lTemp( argv+1 , argv+argc );
  FromVector( lTemp );
}

void Configuration::FromVector( const std::vector< std::string >& aArgs )
{
  typedef std::string tS;
  typedef std::vector<std::string> tVS;
  typedef unsigned tU;
  typedef std::vector<unsigned> tVU;
  typedef double tD;
  typedef std::vector<double> tVD;
  typedef std::size_t tZ;

  tD sigLo , sigHi , rLo , rHi , tLo , tHi;
  tU Nsig(0) , Nr(0) , Nt(0);
  tVD SigKeys, SigVals;

  po::positional_options_description lPositional;
  lPositional.add( "input-file" , 1 );

  po::options_description lDesc("General options");
  lDesc.add_options()
    ( "help",         po::bool_switch()                          ->notifier( [&]( const bool& aArg ){ if( aArg ) { std::cout << lDesc << std::endl; exit(0); } } ) , "produce help message" )
    ( "cfg",          po::value<tS>()                            ->notifier( [&]( const   tS& aArg ){ config_file( lDesc , aArg ); } )                        , "Config file" )
    ( "centre",       po::value<tVS>()->composing()->multitoken()->notifier( [&]( const  tVS& aArg ){ SetCentre( StrToDist(aArg.at(0)) , StrToDist(aArg.at(1)) ); } ) , "Centre of ROI as 'x y' pair" )
    ( "width",        po::value<tVS>()->composing()->multitoken()->notifier( [&]( const  tVS& aArg ){ SetWidth( StrToDist(aArg.at(0)) , StrToDist(aArg.at(1)) ); } ) , "Width of ROI as 'x y' pair" )
    ( "sigma-bins",   po::value<tU>(&Nsig)                                                                                                                    , "Number of sigma bins" )
    ( "sigma-low",    po::value<tS>()                            ->notifier( [&]( const   tS& aArg ){ sigLo=StrToDist(aArg); } )                              , "Lower sigma integration bound" )
    ( "sigma-high",   po::value<tS>()                            ->notifier( [&]( const   tS& aArg ){ sigHi=StrToDist(aArg); } )                              , "High sigma integration bound" )
    ( "sigma-curve",  po::value<tVS>()->composing()->multitoken()->notifier( [&]( const  tVS& aArg ){ for( auto& i : aArg ) { 
                                                                                                        std::vector<std::string> lStrs; 
                                                                                                        boost::split( lStrs , i , [](char c){return c==':';} ); 
                                                                                                        SigKeys.push_back( StrToDist( lStrs.at(0) ) ); 
                                                                                                        SigVals.push_back( std::stod( lStrs.at(1) ) ); 
                                                                                                      } 
                                                                                                    } )                                                       , "Parameterized sigma probability curve (list of colon-separated size-probability pairs)" )           
    ( "r-bins",       po::value<tU>(&Nr)                                                                                                                      , "Number of R bins" )
    ( "r-low",        po::value<tS>()                             ->notifier( [&]( const   tS& aArg ){ rLo=StrToDist(aArg); } )                               , "Lower R bound" )
    ( "r-high",       po::value<tS>()                             ->notifier( [&]( const   tS& aArg ){ rHi=StrToDist(aArg); } )                               , "High R bound" )
    ( "t-bins",       po::value<tU>(&Nt)                                                                                                                      , "Number of T bins" )
    ( "t-low",        po::value<tS>()                             ->notifier( [&]( const   tS& aArg ){ tLo=StrToDist(aArg); } )                               , "Lower T bound" )
    ( "t-high",       po::value<tS>()                             ->notifier( [&]( const   tS& aArg ){ tHi=StrToDist(aArg); } )                               , "High T bound" )
    ( "pb",           po::value<tD>()                             ->notifier( [&]( const   tD& aArg ){ SetPb(aArg); } )                                       , "pb parameter" )
    ( "alpha",        po::value<tD>()                             ->notifier( [&]( const   tD& aArg ){ SetAlpha(aArg); } )                                    , "alpha parameter" )
    ( "validate,v",   po::bool_switch()                           ->notifier( [&]( const bool& aArg ){ SetValidate( aArg ); } )                               , "validate clusters" )
    ( "input-file,i", po::value<tS>()                             ->notifier( [&]( const   tS& aArg ){ SetInputFile(aArg); } )                                , "input file")
    ( "output-file,o", po::value<tS>()                            ->notifier( [&]( const   tS& aArg ){ SetOutputFile(aArg); } )                               , "output file")

    ( "r",            po::value<tS>()                             ->notifier( [&]( const   tS& aArg ){ mClusterR=StrToDist(aArg); } )                         , "R for clustering" )
    ( "t",            po::value<tS>()                             ->notifier( [&]( const   tS& aArg ){ mClusterT=StrToDist(aArg); } )                         , "T for clustering" )
    ( "threads",      po::value<tZ>( &Nthreads )                                                                                                              , "Number of threads to use (default is value given by std::threads::hardware_concurrency())" )
  ;

  po::variables_map lVm;
  po::store( po::command_line_parser( aArgs ).options( lDesc ).positional( lPositional ).run() , lVm );
  po::notify( lVm );    
 
  if( Nr ) SetRBins( Nr , rLo , rHi );
  if( Nt ) SetTBins( Nt , tLo, tHi );

  if( Nsig )
  {
    GSLInterpolator lInterpolator( gsl_interp_cspline, SigKeys , SigVals );
    SetSigmaParameters( Nsig , sigLo , sigHi , [&]( const double& aPt ){ return lInterpolator.Eval( aPt ); } );  
  }

}
