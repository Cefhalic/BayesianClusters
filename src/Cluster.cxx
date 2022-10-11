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
#include <algorithm>
  
/* ===== Local utilities ===== */
#include "Utilities/ProgressBar.hpp"
#include "Utilities/ListComprehension.hpp"
#include "Utilities/Display.hpp"

/* ===== For Root ===== */
#include "TAxis.h"
#include "TGraph.h"
#include "TEllipse.h"


TGraph* DrawGraphs( const std::vector< const Data* >& aPoints , const int& aColour )
{
  auto x( [](const Data* j){return j->x;} | aPoints ) , y( [](const Data* j){return j->y;} | aPoints );
  TGraph* lGraph = new TGraph( x.size() , x.data() , y.data() );
  lGraph->SetMarkerColor( aColour );
  return lGraph;
}


TEllipse* DrawEllipse( const std::vector< const Data* >& aPoints , const int& aColour )
{
  auto x( [](const Data* j){return j->x;} | aPoints ) , y( [](const Data* j){return j->y;} | aPoints );

  double x1(0.0) , y1(0.0) , x2(0.0) , y2(0.0);
  for( auto& i: aPoints )
  {
    x1 += i->x;
    x2 += i->x * i->x;
    y1 += i->y;
    y2 += i->y * i->y;
  }

  x1 /= x.size();
  x2 /= x.size();
  y1 /= y.size();
  y2 /= y.size();

  TEllipse* lEll = new TEllipse( x1 , y1 , 4 * sqrt( x2 - (x1*x1) ) , 4 * sqrt( y2 - (y1*y1) ) );
  lEll->SetFillStyle(0);
  lEll->SetLineColor( aColour );    
  lEll->Draw();

  return lEll;
}


//! Function for an event-display
// \param aEvent The event to draw
void DrawPoints( const EventProxy& aProxy )
{
  std::map< const Cluster* , std::vector< const Data* > > lClusters;
  PRECISION x0( 9e99 ) , x1( -9e99 ) , y0( 9e99 ) , y1( -9e99 );

  for( auto& i : aProxy.mData )
  { 
    Data* lData( i.mData );
    lClusters[ i.mCluster ? i.mCluster->GetParent() : NULL ].push_back( lData );
//    lClusters[ NULL ].push_back( lData );

    x0 = std::min( x0 , lData->x );
    x1 = std::max( x1 , lData->x );
    y0 = std::min( y0 , lData->y );
    y1 = std::max( y1 , lData->y );
  }

  std::cout << "Clusters = " << lClusters.size() << std::endl;

  auto lGraph = DrawGraphs( lClusters[NULL] , 1 );
  lGraph->SetTitle("Data points");
  auto lAxis = lGraph->GetXaxis();
  lAxis->SetLimits(x0,x1); lAxis->SetRangeUser(x0,x1); lAxis->SetLabelSize(0.025); lAxis->SetTickLength(0);
  lAxis = lGraph->GetYaxis();
  lAxis->SetLimits(y0,y1); lAxis->SetRangeUser(y0,y1); lAxis->SetLabelSize(0.025); lAxis->SetTickLength(0);
  lGraph->Draw( "a p" );

  int counter( 0 );
  for( auto& i : lClusters )
  { 
    if( i.first ) DrawGraphs( i.second , counter + 2 ) -> Draw( "p same" );
    counter = ( counter + 1 ) % 8;
  } 

  counter =  0;
  for( auto& i : lClusters )
  { 
    if( i.first ) DrawEllipse( i.second , counter + 2 ) -> Draw( "p same" );
    counter = ( counter + 1 ) % 8;
  } 


}




/* ===== Main function ===== */
int main(int argc, char **argv)
{

  std::cout << "+------------------------------------+" << std::endl;
  ProgressBar2 lBar( "| Cluster. Andrew W. Rose. 2022 |" , 1 );
  std::cout << "+------------------------------------+" << std::endl;
  Configuration::Instance.FromCommandline( argc , argv );
  std::cout << "+------------------------------------+" << std::endl;

  Event lEvent;  

  if( ( Configuration::Instance.ClusterR() < 0 ) or ( Configuration::Instance.ClusterT() < 0 ) ) throw std::runtime_error( "Must specify r and t" );

  const std::string& lFilename = Configuration::Instance.outputFile();

  if( lFilename.size() == 0 )
  {
    lEvent.Clusterize( 
      Configuration::Instance.ClusterR() , 
      Configuration::Instance.ClusterT() , 
      [&]( const EventProxy& aEvent ){ Display( [&](){ DrawPoints( aEvent ); } ); } 
    ); 
  }
  else
  {
    lEvent.Clusterize( 
      Configuration::Instance.ClusterR() , 
      Configuration::Instance.ClusterT() , 
      [&]( const EventProxy& aEvent ){ ToFile( lFilename , [&](){ DrawPoints( aEvent ); } ); } 
    ); 
  }

  std::cout << "+------------------------------------+" << std::endl;

}
