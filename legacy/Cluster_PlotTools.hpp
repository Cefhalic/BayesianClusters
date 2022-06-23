// #pragma once

// /* ===== C++ ===== */
// #include <vector>
// #include <map>

// class Data;
// class TH2D;


// /* ===== Function for plotting data ===== */
// void DrawPoints( const std::vector< Data >& aData );

// /* ===== Function for plotting weights ===== */
// // void DrawWeights()
// // {
// //   int Counter(0);
// //   TGraph2D* lGraph = new TGraph2D( 201 * 201 );
// //   for( int i(-100) ; i!=101 ; ++i )
// //     for( int j(-100) ; j!=101 ; ++j )
// //       lGraph-> SetPoint( Counter++ , i/100.0 , j/100.0 , EdgeCorrectedWeight( Data( i/100.0 , j/100.0 ) , 0.2 ) );

// //   auto FormatAxis = []( TAxis* aAxis ){ aAxis->SetRangeUser(-1,1); aAxis->SetLabelSize(0); aAxis->SetTickLength(0); };
// //   FormatAxis( lGraph->GetXaxis() );
// //   FormatAxis( lGraph->GetYaxis() );
// //   lGraph->GetZaxis()->SetLimits(0,4);
// //   lGraph->Draw( "surf1" );  
// // }

// /* ===== Function for plotting output hist ===== */
// void DrawHisto( TH2D* aHist );

// /* ===== Function for plotting data ===== */
// void DrawClusters( const std::vector< Data >& aData );
