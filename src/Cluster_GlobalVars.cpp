

/* ===== Local utilities ===== */
#include "ListComprehension.hpp"
#include "Cluster_GlobalVars.hpp"

/* ===== For Root ===== */
#include "Math/SpecFunc.h" 
#include "Math/Interpolator.h" 

/* ===== C++ ===== */
#include <iostream>
#include <fstream>
#include <streambuf>

/* ===== BOOST libraries ===== */
#include "boost/algorithm/string.hpp"


GlobalVars::GlobalVars() :
	mScale(-1) , mScale2(-1),
	mSigmacount(-1), mSigmaspacing(-1),
	mMaxR(-1), mMaxR2(-1), mMax2R(-1), mMax2R2(-1),
	mMinScanR(-1), mMaxScanR(-1), mMinScanT(-1), mMaxScanT(-1),
	mDR(-1), mDT(-1),
	mRbins(-1),  mTbins(-1),
	mLogPb(-1), mLogPbDagger(-1), 
	mAlpha(-1), mLogAlpha(-1), mLogGammaAlpha(-1),
	mValidate(false)
{}


void GlobalVars::SetCentre( const double& aPhysicalCentreX , const double& aPhysicalCentreY )
{
	std::cout << "Centre: x=" << aPhysicalCentreX << ", y=" << aPhysicalCentreY << std::endl;
  mPhysicalCentreX = aPhysicalCentreX;
  mPhysicalCentreY = aPhysicalCentreY;
}

void GlobalVars::SetZoom( const double& aScale )
{
	std::cout << "Zoom: " << aScale << std::endl;
	mScale = 2.0 / aScale;
	mScale2 = mScale * mScale;
}

void GlobalVars::SetSigmaParameters( const std::size_t& aSigmacount , const double& aSigmaMin , const double& aSigmaMax , const std::function< double( const double& ) >& aInterpolator )
{
	if( mScale < 0 ) throw std::runtime_error( "Scale must be set before setting Sigma parameters" );

	std::cout << "Sigma-integral: " << aSigmaMin << " to " << aSigmaMax << " in " << aSigmacount << " steps" << std::endl;

	mSigmacount = aSigmacount;
	auto lSigmaspacing = ( aSigmaMax - aSigmaMin ) / aSigmacount;
	auto lSigmabins = [ & ]( const int& i ){ return ( i * lSigmaspacing ) + aSigmaMin;  } | range( mSigmacount );

	mSigmabins = [ & ]( const double& i ){ return toAlgorithmUnits( i ); } | lSigmabins;
	mSigmabins2 = []( const double& i ){ return i * i; } | mSigmabins;
	mProbabilitySigma = aInterpolator | lSigmabins;
	mLogProbabilitySigma = []( const double& w){ return log(w); } | mProbabilitySigma;

	mSigmaspacing = toAlgorithmUnits( lSigmaspacing );
}

void GlobalVars::SetRBins( const std::size_t& aRbins , const double& aMinScanR , const double& aMaxScanR )
{
	mRbins = aRbins;
	mMinScanR = toAlgorithmUnits( aMinScanR );
	mMaxScanR = toAlgorithmUnits( aMaxScanR );
	mDR = ( mMaxScanR - mMinScanR ) / mRbins;

	std::cout << "R-bins: " << aRbins << " bins from " << aMinScanR << " to " << aMaxScanR << " in steps of " << toPhysicalUnits( mDR ) << std::endl;

	mMaxR = toAlgorithmUnits( aMaxScanR );
	mMaxR2 = mMaxR * mMaxR;
	mMax2R = 2.0 * mMaxR;
	mMax2R2 = mMax2R * mMax2R;	
}

void GlobalVars::SetTBins( const std::size_t& aTbins , const double& aMinScanT , const double& aMaxScanT )
{
	mTbins = aTbins;
	mMinScanT = toAlgorithmUnits( aMinScanT );
	mMaxScanT = toAlgorithmUnits( aMaxScanT );
	mDT = ( mMaxScanT - mMinScanT ) / mTbins;

	std::cout << "T-bins: " << aTbins << " bins from " << aMinScanT << " to " << aMaxScanT << " in steps of " << toPhysicalUnits( mDT ) << std::endl;
}

void GlobalVars::SetPb( const double& aPB )
{
	std::cout << "pb: " << aPB << std::endl;
	mLogPb = log( aPB );
	mLogPbDagger = log( 1-aPB );
}

void GlobalVars::SetAlpha( const double& aAlpha )
{
	std::cout << "alpha: " << aAlpha << std::endl;
	mAlpha = aAlpha;
	mLogAlpha = log( aAlpha );
	mLogGammaAlpha = ROOT::Math::lgamma( aAlpha );
}

void GlobalVars::SetValidate( const bool& aValidate )
{
	std::cout << "validate: " << aValidate << std::endl;

	mValidate = aValidate;
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



std::string GlobalVars::FromCommandline( int argc , char **argv )
{
  typedef std::string tS;
  typedef std::vector<std::string> tVS;
  typedef unsigned tU;
  typedef std::vector<unsigned> tVU;
  typedef double tD;
  typedef std::vector<double> tVD;

  tD Cx , Cy , Z , sigLo , sigHi , pb , alpha , normalization , rLo , rHi , tLo , tHi;
  tU Nsig , Nr , Nt;
  tVD SigKeys, SigVals;
  bool val( false );
  tS lInput;

  po::positional_options_description lPositional;
  lPositional.add( "input-file" , 1 );

  po::options_description lDesc("General options");
  lDesc.add_options()
    ( "help",         po::bool_switch()                          ->notifier( [&]( const bool& aArg ){ if( aArg ) { std::cout << lDesc << std::endl; exit(0); } } ) , "produce help message" )
    ( "cfg",          po::value<tS>()                            ->notifier( [&]( const   tS& aArg ){ config_file( lDesc , aArg ); } )                        , "Config file" )
    ( "centre",       po::value<tVS>()->composing()->multitoken()->notifier( [&]( const  tVS& aArg ){ Cx=StrToDist(aArg.at(0)); Cy=StrToDist(aArg.at(1)); } ) , "Centre of ROI as 'x y' pair" )
    ( "zoom",         po::value<tS>()                            ->notifier( [&]( const   tS& aArg ){ Z=StrToDist(aArg); } )                                  , "Dimension of ROI" )
    ( "sigma-bins",   po::value<tU>(&Nsig)                                                                                                                    , "Number of sigma bins" )
    ( "sigma-low",    po::value<tS>()                            ->notifier( [&]( const   tS& aArg ){ sigLo=StrToDist(aArg); } )                              , "Lower sigma integration bound" )
    ( "sigma-high",   po::value<tS>()                            ->notifier( [&]( const   tS& aArg ){ sigHi=StrToDist(aArg); } )                              , "High sigma integration bound" )
    ( "sigma-curve",  po::value<tVS>()->composing()->multitoken()->notifier( [&]( const  tVS& aArg ){ for( auto& i : aArg ) { 
                                                                                                        std::vector<std::string> lStrs; 
                                                                                                        boost::split( lStrs , i , [](char c){return c==':';} ); 
                                                                                                        SigKeys.push_back( StrToDist( lStrs.at(0) ) ); 
                                                                                                        SigVals.push_back( std::stoll( lStrs.at(1) ) ); 
                                                                                                      } 
                                                                                                    } )                                                       , "Parameterized sigma probability curve (list of colon-separated size-probability pairs)" )           
    ( "r-bins",       po::value<tU>(&Nr)                                                                                                                      , "Number of R bins" )
    ( "r-low",        po::value<tS>()                             ->notifier( [&]( const   tS& aArg ){ rLo=StrToDist(aArg); } )                               , "Lower R bound" )
    ( "r-high",       po::value<tS>()                             ->notifier( [&]( const   tS& aArg ){ rHi=StrToDist(aArg); } )                               , "High R bound" )
    ( "t-bins",       po::value<tU>(&Nt)                                                                                                                      , "Number of T bins" )
    ( "t-low",        po::value<tS>()                             ->notifier( [&]( const   tS& aArg ){ tLo=StrToDist(aArg); } )                               , "Lower T bound" )
    ( "t-high",       po::value<tS>()                             ->notifier( [&]( const   tS& aArg ){ tHi=StrToDist(aArg); } )                               , "High T bound" )
    ( "pb",           po::value<tD>(&pb)                                                                                                                      , "pb parameter" )
    ( "alpha",        po::value<tD>(&alpha)                                                                                                                   , "alpha parameter" )
    ( "validate,v",   po::bool_switch(&val)                                                                                                                   , "validate clusters" )
    ( "input-file,i", po::value<tS>(&lInput)                                                                                                                  , "input file")
  ;

  po::variables_map lVm;
  po::store( po::command_line_parser( argc , argv ).options( lDesc ).positional( lPositional ).run() , lVm );
  po::notify( lVm );    

  SetCentre( Cx , Cy );   
  SetZoom( Z );
  SetPb( pb );
  SetAlpha( alpha );
  SetRBins( Nr , rLo , rHi );
  SetTBins( Nt , tLo, tHi );
  SetValidate( val );

  ROOT::Math::Interpolator lInterpolator( SigKeys , SigVals ); // Default to cubic spline interpolation
  SetSigmaParameters( Nsig , sigLo , sigHi , [&]( const double& aPt ){ return lInterpolator.Eval( aPt ); } );  



  return lVm["input-file"].as<tS>();
}
