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


/* ===== Perform clustering ===== */
double ClusterOne( const int& aIndex , const std::vector<Data>& aData , const double& maxdR , const double& maxdR2 )
{
  auto lPlus( aData.begin() + aIndex );
  auto lMinus( aData.rbegin() + aData.size() - 1 - aIndex );

   // Clear the list of neighbours
  std::vector< double > NeighbourR;
  NeighbourR.reserve( aData.size() - 1 );

  double c(0);

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
  uint32_t lCounter(0);
  for( auto lDist : NeighbourR )
  { 
    lCumWeight += EdgeCorrectedWeight( *lPlus , lDist );
    double LocalizationScore = sqrt( (4.0 / pi) * lCumWeight / ++lCounter );
    c += LocalizationScore;
  }

  return c;
}




/* ===== Perform clustering ===== */
void ClusterAll( std::vector<Data>& aData , const double& maxdR )
{
  std::cout << "Sorting" << std::endl;
  std::sort( aData.begin() , aData.end() );

  const double maxdR2( maxdR * maxdR );

// Single-threaded version
// double c(0);
// ProgressBar lProgressBar( "Clustering" , aData.size() );
// for( uint32_t i(0) ; i!=aData.size() ; ++i , ++lProgressBar ) c += ClusterOne( i , aData , maxdR , maxdR2 );
// std::cout << c << std::endl;

  ProgressBar2 lProgressBar( "Clustering" , aData.size() );
  auto CreateLambda = [ &aData, &maxdR, &maxdR2 ]( const uint32_t& i ){ return [ &i, &aData, &maxdR, &maxdR2 ](){ return ClusterOne( i, aData, maxdR, maxdR2 ); }; }; // Create a lambda which returns a lambda
  auto lRet = Vectorize( CreateLambda | range( aData.size() ) );                                                                                                      // Then use list comprehension to create a vector of lambdas and pass to the threadpool
                                                                                                                                                                      // Nerd level cranked to 11

  std::cout << std::accumulate( lRet.begin() , lRet.end() , 0.0 ) << std::endl;
}


/* ===== Main function ===== */
int main(int argc, char **argv)
{
  auto lData = CreatePseudoData( 70000 , 500 , 500 , .01 );


  ClusterAll( lData , 0.03 );

  // InteractiveDisplay( [](){ DrawWeights(); } );
  // InteractiveDisplay( [ lData ](){ DrawPoints( lData ); } );

  return 0;
}
