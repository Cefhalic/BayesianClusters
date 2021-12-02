#pragma once

/* ===== For Root ===== */
#include "TCanvas.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TApplication.h"
#include "TRootCanvas.h"


void DisplayHelper( TCanvas* c , const int& i ){}

template< typename First , typename... Rest >
void DisplayHelper( TCanvas* c , const int& i , const First& aFirst , Rest&&... aRest )
{
  c->cd( i );
  aFirst();
  DisplayHelper(  c , i+1 , std::forward<Rest>( aRest )... );
}


/* ===== Display the data on a ROOT canvas ===== */
template< typename ... Functors >
void InteractiveDisplay( Functors&&... aFunctors )
{
  std::cout << "Plotting" << std::endl;

  int a(0);
  TApplication app("app", &a , NULL );

  gROOT->Reset ( ) ; // re−initialize ROOT
  gROOT->SetStyle ( "Plain" ) ; // set empty TStyle ( nicer o npaper )
  gStyle->SetOptStat ( 0 ) ; // print statistics on plots , ( 0 ) for no output
  gStyle->SetPalette ( 1 ) ; // set nicer colors than default
  // gStyle->SetOptTitle ( 0 ) ; // suppress title box

  constexpr std::size_t cnt( sizeof...(Functors) );

  TCanvas* c = new TCanvas( "c" , "" , 0 , 0 , 1600 , 800 );

  std::size_t i = ceil( sqrt( cnt ) );
  std::size_t j = ceil( double(cnt) / i );

  c->Divide( i , j );
  DisplayHelper( c , 1 , std::forward<Functors>( aFunctors)... );

  TRootCanvas *rc = (TRootCanvas *)c->GetCanvasImp();
  rc->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");
  app.Run();
}