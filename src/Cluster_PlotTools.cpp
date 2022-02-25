

/* ===== Cluster sources ===== */
#include "Cluster_Data.hpp"
#include "Cluster_PlotTools.hpp"

/* ===== C++ ===== */
#include <vector>

/* ===== For Root ===== */
#include "TCanvas.h"
#include "TROOT.h"
#include "TAxis.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TH2D.h"
#include "TEllipse.h"

/* ===== Local utilities ===== */
#include "ListComprehension.hpp"


/* ===== Function for plotting data ===== */
void DrawPoints( const std::vector< Data >& aData )
{
  // auto x( &Data::x | aData ) , y( &Data::y | aData ) , z( &Data::mLocalizationScore | aData );
  // gPad -> SetMargin( 0.01 , 0.15 , 0.01 , 0.01 );
  // TGraph* lGraph0 = new TGraph( aData.size() , x.data() , y.data() );
  // TGraph2D* lGraph1 = new TGraph2D( aData.size() , x.data() , y.data() , z.data() );
  // auto FormatAxis = []( TAxis* aAxis ){ aAxis->SetLimits(-1,1); aAxis->SetRangeUser(-1,1); aAxis->SetLabelSize(0); aAxis->SetTickLength(0); };
  // FormatAxis( lGraph0->GetXaxis() );
  // FormatAxis( lGraph0->GetYaxis() );
  // lGraph0->Draw( "a p" );
  // lGraph1->Draw( "cont1z same" );
}


/* ===== Function for plotting weights ===== */
// void DrawWeights()
// {
//   int Counter(0);
//   TGraph2D* lGraph = new TGraph2D( 201 * 201 );
//   for( int i(-100) ; i!=101 ; ++i )
//     for( int j(-100) ; j!=101 ; ++j )
//       lGraph-> SetPoint( Counter++ , i/100.0 , j/100.0 , EdgeCorrectedWeight( Data( i/100.0 , j/100.0 ) , 0.2 ) );

//   auto FormatAxis = []( TAxis* aAxis ){ aAxis->SetRangeUser(-1,1); aAxis->SetLabelSize(0); aAxis->SetTickLength(0); };
//   FormatAxis( lGraph->GetXaxis() );
//   FormatAxis( lGraph->GetYaxis() );
//   lGraph->GetZaxis()->SetLimits(0,4);
//   lGraph->Draw( "surf1" );  
// }


/* ===== Function for plotting output hist ===== */
void DrawHisto( TH2D* aHist )
{
  gPad -> SetLeftMargin( 0.15 );
  gPad -> SetRightMargin( 0.15 );
  // gPad->SetLogz();
  aHist->SetContour(1e6);
  aHist->Draw("colz");  
}


/* ===== Function for plotting data ===== */
void DrawClusters( const std::vector< Data >& aData )
{
  // std::map< const Cluster* , std::vector< const Data* > > lClusters;
  // for( auto& i : aData ) lClusters[ i.mCluster ? i.mCluster->GetParent() : NULL ].push_back( &i );

  // gPad -> SetMargin( 0.01 , 0.15 , 0.01 , 0.01 );

  // auto GetX = []( const Data* i ){ return i->x; };
  // auto GetY = []( const Data* i ){ return i->y; };

  // auto i( lClusters.begin() );
  // TGraph* lGraph = new TGraph( i->second.size() , ( GetX | i->second ).data() , ( GetY | i->second ).data() );

  // auto FormatAxis = []( TAxis* aAxis ){ aAxis->SetLimits(-1,1); aAxis->SetRangeUser(-1,1); aAxis->SetLabelSize(0); aAxis->SetTickLength(0); };
  // FormatAxis( lGraph->GetXaxis() );
  // FormatAxis( lGraph->GetYaxis() );
  // lGraph->SetMarkerColor( 1 );
  // lGraph->Draw( "ap" );

  // i++;
  // int cnt(0);
  // for( ; i != lClusters.end() ; ++i , ++cnt )
  // {
  //   double x(0.0) , y(0.0) , x2(0.0) , y2(0.0);
  //   for( auto& j: i->second )
  //   {
  //     x += (j->x);
  //     x2 += (j->x)*(j->x);
  //     y += (j->y);      
  //     y2 += (j->y)*(j->y);
  //   }

  //   x /= i->second.size();
  //   x2 /= i->second.size();
  //   y /= i->second.size();
  //   y2 /= i->second.size();

  //   TEllipse *el1 = new TEllipse( x , y , 4 * sqrt( x2 - (x*x) ) , 4 * sqrt( y2 - (y*y) ) );
  //   el1->SetFillStyle(0);
  //   el1->SetLineColor( (cnt%8)+2 );    
  //   el1->Draw();

  //   TGraph* lGraph = new TGraph( i->second.size() , ( GetX | i->second ).data() , ( GetY | i->second ).data() );
  //   lGraph->SetMarkerColor( (cnt%8)+2 );
  //   lGraph->Draw( "samep" );

  // }
}
