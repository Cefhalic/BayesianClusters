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


  float x0( 9e99 ) , x1( -9e99 ) , y0( 9e99 ) , y1( -9e99 );
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
  TGraph* lGraph0 = new TGraph( x.size() , x.data() , y.data() );
  lGraph0 -> SetTitle("Data points");
  auto FormatAxis = [&]( TAxis* aAxis ){ aAxis->SetLimits(x0,x1); aAxis->SetRangeUser(y0,y1); aAxis->SetLabelSize(0.025); aAxis->SetTickLength(0); };
  FormatAxis( lGraph0->GetXaxis() );
  FormatAxis( lGraph0->GetYaxis() );
  lGraph0->Draw( "a p" );
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

  Display( [ & ](){ DrawPoints( lEvent ); } );


  std::cout << "+------------------------------------+" << std::endl;

}
