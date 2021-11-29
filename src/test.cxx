/* ===== C++ ===== */
#include <iostream>
#include <vector>
#include <list>
#include <numeric>
#include <math.h> 

/* ===== For Root ===== */
#include "TRandom3.h"
#include "TAxis.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TH2F.h"
#include "TEllipse.h"

/* ===== Local utilities ===== */
#include "ListComprehension.hpp"
#include "ProgressBar.hpp"
#include "RootWindow.hpp"
#include "Vectorize.hpp"


constexpr double pi = atan(1)*4;

// class Data;

// struct Cluster
// {
//   std::list< Data* > members;
//   void merge( Cluster* aOther );
// };


/* ===== Struct for storing data ===== */
class Data
{
public:
  Data( double aX , double aY ) : x(aX) , y(aY) , r( sqrt( (aX*aX) + (aY*aY) ) ), phi( atan2( aY , aX ) ), parent( NULL ){}

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

  Data* GetParent()
  {
    if( parent == NULL ) return NULL;    
    if( parent == this ) return this;
    return parent->GetParent();
  }

  void SetParent( Data* aParent )
  {
    if( parent == this ) parent = aParent;
    else parent->SetParent( aParent );
  }

  void Cluster( Data& aData )
  {
    SetParent( aData.GetParent() );
    // std::cout << parent << " | " << this << " | " << aData.parent << " | " << &aData << std::endl;
    // if( aData.parent == &aData ) aData.parent = this; // The other cluster is it's own cluster - take ownership
    // else if( parent == this ) parent = aData.parent;        // The other cluster is in a cluster and we are not - we join theirs
    //else std::cout << "2" << std::endl;
  }

public:
  double x, y, r, phi;
  Data* parent;
};


// void Cluster::merge( Cluster* aOther )
// {
//   for( auto i : aOther->members ) i->cluster = this;
//   members.splice( members.end() , aOther->members , aOther->members.begin() , aOther->members.end() );
// }



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


/* ===== Function for plotting data ===== */
void DrawPoints( const std::map< Data* , std::vector< Data* > >& aData )
{
  gPad -> SetMargin( 0.01 , 0.01 , 0.01 , 0.01 );

  auto GetX = []( const Data* i ){ return i->x; };
  auto GetY = []( const Data* i ){ return i->y; };

  auto i( aData.begin() );
  TGraph* lGraph = new TGraph( i->second.size() , ( GetX | i->second ).data() , ( GetY | i->second ).data() );

  auto FormatAxis = []( TAxis* aAxis ){ aAxis->SetRangeUser(-1,1); aAxis->SetLabelSize(0); aAxis->SetTickLength(0); };
  FormatAxis( lGraph->GetXaxis() );
  FormatAxis( lGraph->GetYaxis() );
  lGraph->SetMarkerColor( 1 );
  lGraph->Draw( "ap" );

  i++;
  int cnt(0);
  for( ; i != aData.end() ; ++i , ++cnt )
  {
    double x(0.0) , y(0.0) , x2(0.0) , y2(0.0);
    for( auto& j: i->second )
    {
      x += (j->x);
      x2 += (j->x)*(j->x);
      y += (j->y);      
      y2 += (j->y)*(j->y);
    }

    x /= i->second.size();
    x2 /= i->second.size();
    y /= i->second.size();
    y2 /= i->second.size();

    // std::cout << sqrt( x2 - (x*x) ) << " " << sqrt( y2 - (y*y) ) << std::endl;

    TEllipse *el1 = new TEllipse( x , y , 4 * sqrt( x2 - (x*x) ) , 4 * sqrt( y2 - (y*y) ) );
    el1->SetFillStyle(0);
    el1->SetLineColor( (cnt%8)+2 );    
    el1->Draw();

    TGraph* lGraph = new TGraph( i->second.size() , ( GetX | i->second ).data() , ( GetY | i->second ).data() );
    lGraph->SetMarkerColor( (cnt%8)+2 );
    lGraph->Draw( "samep" );

  }
}



/* ===== Perform cluster profiling ===== */
TH2F* ProfileOneDatum( const int& aIndex , const std::vector<Data>& aData , const double& maxdR )
{
  const double maxdR2( maxdR * maxdR ); // Could be calculated outside and passed in, but cost is minimal and code is cleaner this way

  const auto lPlus( aData.begin() + aIndex );
  const auto lMinus( aData.rbegin() + aData.size() - 1 - aIndex );

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


/* ===== Perform cluster profiling ===== */
TH2F* ProfileAllData( std::vector<Data>& aData , const double& maxdR )
{
  std::cout << "Sorting" << std::endl;
  std::sort( aData.begin() , aData.end() );

  std::vector< TH2F*> lRet;

  {
    ProgressBar2 lProgressBar( "Profiling Data" , aData.size() );
    auto CreateLambda = [ &aData, &maxdR ]( const uint32_t& i ){ return [ &i, &aData, &maxdR ](){ return ProfileOneDatum( i, aData, maxdR ); }; }; // Create a lambda which returns a lambda
    lRet = Vectorize( CreateLambda | range( aData.size() ) );                                                                                      // Then use list comprehension to create a vector of lambdas and pass to the threadpool
  }                                                                                                                                                // Nerd level cranked to 11

  std::cout << "Aggregating output" << std::endl;
  std::unordered_set<TH2F*> s;
  for ( auto& i : lRet ) s.insert(i); // Apparently much faster than using the constructor!!

  TH2F* lHist = new TH2F( "Hist" , "Hist;r;L(r)" , 101 , 0.0 , maxdR , 101 , 0.0 , 2.5 * maxdR );
  for( auto& i: s ) *lHist = *lHist + *i;

  return lHist;
}





/* ===== Perform clustering ===== */
bool ClusterOneDatum( const int& aIndex , std::vector<Data>& aData, const double& maxdR , const double& T )
{
  const double maxdR2( maxdR * maxdR ); // Could be calculated outside and passed in, but cost is minimal and code is cleaner this way

  const auto lPlus( aData.begin() + aIndex );
  const auto lMinus( aData.rbegin() + aData.size() - 1 - aIndex );

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

  // Iterate over sorted distances, updating the localization score for each distance
  double lCumWeight( 0.0 );
  for( auto& lDist : NeighbourR ) lCumWeight += EdgeCorrectedWeight( *lPlus , lDist );

  const double LocalizationConstant( 4.0 / ( pi * (aData.size()-1) ) ); 
  
  if( sqrt( LocalizationConstant * lCumWeight ) > T )
  {
    lPlus->parent = &*lPlus;
    // lPlus->cluster = new Cluster();
    // lPlus->cluster->members.push_back( &*lPlus );
    return true;
  }
  else
  {
    // lPlus->parent = NULL;
    return false;
  }
}


/* ===== Perform clustering ===== */
void ClusterOneDatum2( const int& aIndex , std::vector<Data>& aData, const double& maxdR )
{
  const auto lPlus( aData.begin() + aIndex );
  if( lPlus->parent == NULL ) return; // We have been assigned to the background

  const auto lMinus( aData.rbegin() + aData.size() - 1 - aIndex );

  const double maxdR2( maxdR * maxdR ); // Could be calculated outside and passed in, but cost is minimal and code is cleaner this way

  // Iterate over other hits and populate the neighbour list
  for( auto lPlusIt( lPlus + 1 ); lPlusIt != aData.end() ; lPlusIt++ )
  {
    if( lPlusIt->parent == NULL ) continue; // We have been assigned to the background
    if( ( lPlusIt->r - lPlus->r ) > maxdR ) break; // lPlusIt is always further out than lPlus
    if( lPlus->dR2( *lPlusIt ) < maxdR2 ) lPlusIt->Cluster( *lPlus );
  }

  for( auto lMinusIt( lMinus + 1 ); lMinusIt != aData.rend() ; lMinusIt++ )
  {
    if( lMinusIt->parent == NULL ) continue; // We have been assigned to the background
    if( ( lMinus->r - lMinusIt->r ) > maxdR ) break; // lMinus is always further out than lMinusIn
    if( lMinus->dR2( *lMinusIt ) < maxdR2 ) lMinusIt->Cluster( *lMinus );
  }
}






/* ===== Perform clustering ===== */
void ClusterAllData( std::vector<Data>& aData , const double& maxdR , const double& T )
{
  std::cout << "Sorting (unnecessary if profiling already performed)" << std::endl;
  std::sort( aData.begin() , aData.end() );

  std::vector< bool > lRet;

  {
    ProgressBar2 lProgressBar( "Clustering" , aData.size() );
    auto CreateLambda = [ &aData, &maxdR, &T ]( const uint32_t& i ){ return [ &i, &aData, &maxdR, &T ](){ return ClusterOneDatum( i, aData, maxdR, T ); }; }; // Create a lambda which returns a lambda
    lRet = Vectorize( CreateLambda | range( aData.size() ) );                                                                                      // Then use list comprehension to create a vector of lambdas and pass to the threadpool
  }                                                                                                                                                // Nerd level cranked to 11

  // uint32_t lCount( 0 );
  // for( auto i : aData ) lCount += bool( i.cluster );
  // std::cout << lCount << std::endl;

  {
    ProgressBar lProgressBar( "Clustering 2" , aData.size() );
    for( uint32_t i(0) ; i!=aData.size() ; ++i , ++lProgressBar ) ClusterOneDatum2( i , aData , 2*maxdR );
  }

  {
    std::map< Data* , std::vector< Data* > > s;
    for( auto& i : aData )
    { 
      s[ i.GetParent() ].push_back( &i );
    }
    std::cout << s.size() << std::endl;

    InteractiveDisplay( [ s ](){ DrawPoints( s ); } );
  }

}


/* ===== Main function ===== */
int main(int argc, char **argv)
{
   //auto lData = CreatePseudoData( 10000 , 500 , 500 , .01 );
   auto lData = CreatePseudoData( 70000 , 500 , 500 , .01 );

  ClusterAllData( lData , 0.005 , 0.0075 );

  //auto lHist = ProfileAllData( lData , 0.04 );

  // InteractiveDisplay( [](){ DrawWeights(); } );
  //InteractiveDisplay( [ lData ](){ DrawPoints( lData ); } );
  // InteractiveDisplay( [ lData ](){ DrawPoints( lData ); } , [ lHist ](){ DrawHisto( lHist ); } );

  return 0;
}
