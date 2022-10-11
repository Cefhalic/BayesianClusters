/* ===== Cluster sources ===== */
#include "BayesianClustering/Event.hpp"
#include "BayesianClustering/Configuration.hpp"

// /* ===== C++ ===== */
#include <iostream>
#include <algorithm>
  
/* ===== Local utilities ===== */
#include "Utilities/ProgressBar.hpp"
#include "Utilities/ListComprehension.hpp"
#include "Utilities/Display.hpp"

/* ===== For Root ===== */
#include "TAxis.h"
#include "TGraph.h"

//! Function for an event-display
// \param aEvent The event to draw
void DrawPoints( const Event& aEvent )
{

  auto x( &Data::x | aEvent.mData ) , y( &Data::y | aEvent.mData );


  PRECISION x0( 9e99 ) , x1( -9e99 ) , y0( 9e99 ) , y1( -9e99 );
  for( auto& i : x )
  {
    x0 = std::min( x0 , i );
    x1 = std::max( x1 , i );
  }

  for( auto& i : y )
  {
    y0 = std::min( y0 , i );
    y1 = std::max( y1 , i );
  }

  // gPad -> SetMargin( 0.01 , 0.15 , 0.01 , 0.01 );
  TGraph* lGraph = new TGraph( x.size() , x.data() , y.data() );
  lGraph -> SetTitle("Data points");
  auto lAxis = lGraph->GetXaxis();
  lAxis->SetLimits(x0,x1); lAxis->SetRangeUser(x0,x1); lAxis->SetLabelSize(0.025); lAxis->SetTickLength(0);
  lAxis = lGraph->GetYaxis();
  lAxis->SetLimits(y0,y1); lAxis->SetRangeUser(y0,y1); lAxis->SetLabelSize(0.025); lAxis->SetTickLength(0);
  lGraph->Draw( "a p" );
}


/* ===== Main function ===== */
int main(int argc, char **argv)
{


  std::cout << "+------------------------------------+" << std::endl;
  ProgressBar2 lBar( "| Data viewer. Andrew W. Rose. 2022  |" , 1 );
  std::cout << "+------------------------------------+" << std::endl;
  Configuration::Instance.FromCommandline( argc , argv );
  std::cout << "+------------------------------------+" << std::endl;

  Event lEvent;  

  const std::string& lFilename = Configuration::Instance.outputFile();

  if( lFilename.size() == 0 )
  {
    Display( [&](){ DrawPoints( lEvent ); } ); 
  }
  else
  {
    ToFile( lFilename , [&](){ DrawPoints( lEvent ); } );
  }


  std::cout << "+------------------------------------+" << std::endl;

}
