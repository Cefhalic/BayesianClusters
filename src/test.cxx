/* ===== Cluster sources ===== */
#include "Cluster_Data.hpp"
#include "Cluster_DataSources.hpp"
#include "Cluster_PlotTools.hpp"
#include "Cluster_EdgeCorrection.hpp"

// /* ===== C++ ===== */
#include <thread>
#include <iostream>
#include <math.h>
#include <vector>
#include <unordered_set>
// #include <list>
// #include <numeric>

/* ===== For Root ===== */
#include "TH2F.h"

/* ===== Local utilities ===== */
#include "ListComprehension.hpp"
#include "ProgressBar.hpp"
#include "RootWindow.hpp"
#include "Vectorize.hpp"










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
  constexpr double pi = atan(1)*4;  
  const double LocalizationConstant( 4.0 / ( pi * (aData.size()-1) ) ); 


  static thread_local const std::string lTitle( "Histo"+std::to_string( std::hash<std::thread::id>{}(std::this_thread::get_id()) ) );
  static thread_local TH2F* lRet = new TH2F( lTitle.c_str() , lTitle.c_str() , 101 , 0.0 , maxdR , 101 , 0.0 , 2.5 * maxdR );

  const double lRstep( maxdR / 100 );
  const double lTstep( 2.5 * lRstep );

  double lR( lRstep );
  auto lIt( NeighbourR.begin() );

  for( int i(0) ; i!=100 ; ++i , lR += lRstep )
  { 
    while( lIt != NeighbourR.end() )
    {
      if ( *lIt > lR ) break;
      lCumWeight += EdgeCorrectedWeight( *lPlus , *lIt++ );
    }
    if( lCumWeight > 0 )
    {
      auto L_r =  sqrt( LocalizationConstant * lCumWeight );
      auto lT( lTstep );
      for( int j(0) ; j!=100 ; ++j , lT += lTstep )
        if( L_r > lT ) lRet->Fill( lR , lT );
    }
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

  TH2F* lHist = new TH2F( "Hist" , "Hist;r;T" , 101 , 0.0 , maxdR , 101 , 0.0 , 2.5 * maxdR );
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

  constexpr double pi = atan(1)*4;
  const double LocalizationConstant( 4.0 / ( pi * (aData.size()-1) ) ); 
  
  if( sqrt( LocalizationConstant * lCumWeight ) > T )
  {
    lPlus->parent = &*lPlus;
    return true;
  }
  else
  {
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
    ProgressBar2 lProgressBar( "Clustering 1" , aData.size() );
    auto CreateLambda = [ &aData, &maxdR, &T ]( const uint32_t& i ){ return [ &i, &aData, &maxdR, &T ](){ return ClusterOneDatum( i, aData, maxdR, T ); }; }; // Create a lambda which returns a lambda
    lRet = Vectorize( CreateLambda | range( aData.size() ) );                                                                                      // Then use list comprehension to create a vector of lambdas and pass to the threadpool
  }                                                                                                                                                // Nerd level cranked to 11

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
//  auto lData = LoadCSV( "1_un_red.csv" , 1./64000. , -1. , -1. ); // Full file
  auto lData = LoadCSV( "1_un_red.csv" , 1./10000. , 87000. , 32000. ); // One cluster

   //auto lData = CreatePseudoData( 10000 , 500 , 500 , .01 );
   //auto lData = CreatePseudoData( 70000 , 500 , 500 , .01 );

  InteractiveDisplay( [ lData ](){ DrawPoints( lData ); } );

  //ClusterAllData( lData , 0.0035 , 0.007 );

  // auto lHist = ProfileAllData( lData , 3e-2 );
  // InteractiveDisplay( [ lData ](){ DrawPoints( lData ); } , [ lHist ](){ DrawHisto( lHist ); } );

  // InteractiveDisplay( [](){ DrawWeights(); } );
  return 0;
}
