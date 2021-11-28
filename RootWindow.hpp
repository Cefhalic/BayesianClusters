#pragma once

/* ===== For Root ===== */
#include "TCanvas.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TApplication.h"
#include "TRootCanvas.h"

/* ===== Display the data on a ROOT canvas ===== */
template< typename Functor >
void InteractiveDisplay( const Functor& aFunctor )
{
  std::cout << "Plotting" << std::endl;

  int a(0);
  TApplication app("app", &a , NULL );

  gROOT->Reset ( ) ; // re−initialize ROOT
  gROOT->SetStyle ( "Plain" ) ; // set empty TStyle ( nicer o npaper )
  gStyle->SetOptStat ( 11111 ) ; // print statistics on plots , ( 0 ) for no output
  gStyle->SetPalette ( 1 ) ; // set nicer colors than default
  gStyle->SetOptTitle ( 0 ) ; // suppress title box

  TCanvas* c = new TCanvas( "c" , "" , 0 , 0 , 800 , 800 );
  c->SetMargin( 0.01 , 0.01 , 0.01 , 0.01 );

  aFunctor();

  TRootCanvas *rc = (TRootCanvas *)c->GetCanvasImp();
  rc->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");
  app.Run();
}