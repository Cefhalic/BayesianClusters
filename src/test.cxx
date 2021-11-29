/* ===== C++ ===== */
#include <iostream>
#include <vector>
#include <numeric>
#include <math.h> 

/* ===== For Root ===== */
#include "TRandom3.h"
#include "TAxis.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TH2F.h"

/* ===== Local utilities ===== */
#include "ListComprehension.hpp"
#include "ProgressBar.hpp"
#include "RootWindow.hpp"
#include "Vectorize.hpp"


constexpr double pi = atan(1)*4;


/* ===== Struct for storing data ===== */
class Data
{
public:
  Data( double aX , double aY ) : x(aX) , y(aY) , r( sqrt( (aX*aX) + (aY*aY) ) ), phi( atan2( aY , aX ) ){}
  bool operator< ( const Data& aOther ) const { return r < aOther.r; }

  double dR2( const Data& aOther ) const
  {
    double dX( x - aOther.x ), dY( y - aOther.y );
    return ( dX*dX ) + ( dY*dY );
  }

  double dR( const Data& aOther ) const
  {
    return sqrt( dR2( aOther ) );
  }

public:
  double x, y, r, phi;
};


/* ===== Utility function for creating a vector of data ===== */
std::vector< Data > CreatePseudoData( const int& aBackgroundCount , const int& aClusterCount , const int& aClusterSize , const double& aClusterScale )
{
  std::cout << "Generating Pseudodata" << std::endl;

  std::vector< Data > lData;
  lData.reserve( aBackgroundCount + ( aClusterCount * aClusterSize ) );

  TRandom3 lRand( 23423 );

  for( int i(0); i!= aBackgroundCount; ++i )
  {
    double x( lRand.Uniform( -1.0 , 1.0 ) ) , y( lRand.Uniform( -1.0 , 1.0 ) );
    lData.emplace_back( x , y );
  }

  for( int i(0); i!= aClusterCount; ++i )
  {
    double x( lRand.Uniform( -1.0 , 1.0 ) ) , y( lRand.Uniform( -1.0 , 1.0 ) );
    double sigma( fabs( lRand.Gaus( aClusterScale , aClusterScale/3 ) ) );
    for( int j(0) ; j!= aClusterSize ; /* */ )
    {
      double x2( lRand.Gaus( x , sigma ) ) , y2( lRand.Gaus( y , sigma ) );  
      if( x2 > 1 or x2 < -1 or y2 > 1 or y2 < -1 ) continue;    
      lData.emplace_back( x2 , y2 );
      ++j;
    }
  }

  return lData;
}


/* ===== Function for plotting data ===== */
void DrawPoints( const std::vector< Data >& aData )
{
  gPad -> SetMargin( 0.01 , 0.01 , 0.01 , 0.01 );
  TGraph* lGraph = new TGraph( aData.size() , ( &Data::x | aData ).data() , ( &Data::y | aData ).data() );
  auto FormatAxis = []( TAxis* aAxis ){ aAxis->SetRangeUser(-1,1); aAxis->SetLabelSize(0); aAxis->SetTickLength(0); };
  FormatAxis( lGraph->GetXaxis() );
  FormatAxis( lGraph->GetYaxis() );
  lGraph->Draw( "ap" );
}


/* ===== Approximation to correct area of circle overlapping the edge of the square ===== */
double EdgeCorrectedWeight( const Data& aPt , const double& aDist )
{
  double Weight( 1.0 );
  const double X( 1 - fabs( aPt.x ) ) , Y( 1 - fabs( aPt.y ) );
  if( X < aDist )  Weight *= ( 1 + pow( acos( X/aDist ) * (2/pi) , 4) );
  if( Y < aDist )  Weight *= ( 1 + pow( acos( Y/aDist ) * (2/pi) , 4) );
  return Weight;
}


/* ===== Function for plotting weights ===== */
void DrawWeights()
{
  int Counter(0);
  TGraph2D* lGraph = new TGraph2D( 201 * 201 );
  for( int i(-100) ; i!=101 ; ++i )
    for( int j(-100) ; j!=101 ; ++j )
      lGraph-> SetPoint( Counter++ , i/100.0 , j/100.0 , EdgeCorrectedWeight( Data( i/100.0 , j/100.0 ) , 0.2 ) );

  auto FormatAxis = []( TAxis* aAxis ){ aAxis->SetRangeUser(-1,1); aAxis->SetLabelSize(0); aAxis->SetTickLength(0); };
  FormatAxis( lGraph->GetXaxis() );
  FormatAxis( lGraph->GetYaxis() );
  lGraph->GetZaxis()->SetLimits(0,4);
  lGraph->Draw( "surf1" );  
}


/* ===== Function for plotting output hist ===== */
void DrawHisto( TH2F* aHist )
{
  gPad -> SetLeftMargin( 0.15);
  gPad -> SetRightMargin( 0.15);
  gPad->SetLogz();
  aHist->Draw("colz");  
}




/* ===== Perform clustering ===== */
TH2F* ClusterOne( const int& aIndex , const std::vector<Data>& aData , const double& maxdR )
{
  const double maxdR2( maxdR * maxdR ); // Could be calculated outside and passed in, but cost is minimal and code is cleaner this way

  auto lPlus( aData.begin() + aIndex );
  auto lMinus( aData.rbegin() + aData.size() - 1 - aIndex );

   // Create a container for storing distances to neighbours
  std::vector< double > NeighbourR;
  NeighbourR.reserve( aData.size() - 1 );

  // Iterate over other hits and populate the neighbour list
  for( auto lPlusIt( lPlus + 1 ); lPlusIt != aData.end() ; lPlusIt++ )
  {
    if( ( lPlusIt->r - lPlus->r ) > maxdR ) break; // lPlusIt is always further out than lPlus
    double dR2( lPlus->dR2( *lPlusIt ) );
    if( dR2 < maxdR2 ) NeighbourR.push_back( sqrt( dR2 ) );
  }

  for( auto lMinusIt( lMinus + 1 ); lMinusIt != aData.rend() ; lMinusIt++ )
  {
    if( ( lMinus->r - lMinusIt->r ) > maxdR ) break; // lMinus is always further out than lMinusIn
    double dR2( lMinus->dR2( *lMinusIt ) );
    if( dR2 < maxdR2 ) NeighbourR.push_back( sqrt( dR2 ) );
  }

  // Sort the list
  std::sort( NeighbourR.begin() , NeighbourR.end() );

  // Iterate over sorted distances, updating the localization score for each distance
  double lCumWeight( 0.0 );
  const double LocalizationConstant( 4.0 / ( pi * (aData.size()-1) ) ); 


  static thread_local const std::string lTitle( "Histo"+std::to_string( std::hash<std::thread::id>{}(std::this_thread::get_id()) ) );
  static thread_local TH2F* lRet = new TH2F( lTitle.c_str() , lTitle.c_str() , 101 , 0.0 , maxdR , 101 , 0.0 , 2.5 * maxdR );

  const double lRstep( maxdR / 100 );
  double lR( lRstep );
  auto lIt( NeighbourR.begin() );

  for( int i(0) ; i!=100 ; ++i , lR += lRstep )
  { 
    while( lIt != NeighbourR.end() )
    {
      if ( *lIt > lR ) break;
      lCumWeight += EdgeCorrectedWeight( *lPlus , *lIt++ );
    }
    if( lCumWeight > 0 ) lRet->Fill( lR , sqrt( LocalizationConstant * lCumWeight ) );
  }

  return lRet;
}


/* ===== Perform clustering ===== */
TH2F* ClusterAll( std::vector<Data>& aData , const double& maxdR )
{
  std::cout << "Sorting" << std::endl;
  std::sort( aData.begin() , aData.end() );

  std::vector< TH2F*> lRet;

  {
    ProgressBar2 lProgressBar( "Clustering" , aData.size() );
    auto CreateLambda = [ &aData, &maxdR ]( const uint32_t& i ){ return [ &i, &aData, &maxdR ](){ return ClusterOne( i, aData, maxdR ); }; }; // Create a lambda which returns a lambda
    lRet = Vectorize( CreateLambda | range( aData.size() ) );                                                                            // Then use list comprehension to create a vector of lambdas and pass to the threadpool
  }                                                                                                                                        // Nerd level cranked to 11

  std::cout << "Aggregating output" << std::endl;
  std::unordered_set<TH2F*> s;
  for ( auto i : lRet ) s.insert(i); // Apparently much faster than using the constructor!!

  TH2F* lHist = new TH2F( "Hist" , "Hist;r;L(r)" , 101 , 0.0 , maxdR , 101 , 0.0 , 2.5 * maxdR );
  for( auto i: s ) *lHist = *lHist + *i;

  return lHist;
}


/* ===== Main function ===== */
int main(int argc, char **argv)
{
  auto lData = CreatePseudoData( 70000 , 500 , 500 , .01 );
  auto lHist = ClusterAll( lData , 0.04 );


  // InteractiveDisplay( [](){ DrawWeights(); } );
  // InteractiveDisplay(  );
  InteractiveDisplay( [ lData ](){ DrawPoints( lData ); } , [ lHist ](){ DrawHisto( lHist ); } );

  return 0;
}
