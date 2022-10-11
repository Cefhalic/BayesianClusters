#pragma once

/* ===== C++ ===== */
#include <iostream>

/* ===== For Root ===== */
#include "TCanvas.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TApplication.h"
#include "TRootCanvas.h"

//! Helper function to terminate recursion on variadic template calls
//! \param c A ROOT canvas to which we are drawing
//! \param i The index of the subcanvas
inline void DisplayHelper( TCanvas* c , const int& i )
{}

//! Helper function to fill subcanvases via a variadic series of callbacks
//! \tparam First The type of the callback we will issue on this iteration
//! \tparam Rest A variadic tail
//! \param c A ROOT canvas to which we are drawing
//! \param i The index of the subcanvas
//! \param aFirst The callback we will issue on this iteration
//! \param aRest A variadic tail
template< typename First , typename... Rest >
void DisplayHelper( TCanvas* c , const int& i , const First& aFirst , Rest&&... aRest )
{
  c->cd( i );
  aFirst();
  DisplayHelper(  c , i+1 , std::forward<Rest>( aRest )... );
}


//! Create a display window using CERN ROOT and fill it via a variadic series of callbacks 
//! \tparam Functors A parameter pack for a variadic series of callbacks 
//! \param A variadic series of callbacks
template< typename ... Functors >
void Display( Functors&&... aFunctors )
{
  std::cout << "Plotting" << std::endl;

  int a(0);
  TApplication app("app", &a , NULL );

  gROOT->Reset ( ) ; // re−initialize ROOT
  gROOT->SetStyle ( "Plain" ) ; // set empty TStyle ( nicer on paper )
  gStyle->SetOptStat ( 0 ) ; // print statistics on plots , ( 0 ) for no output
  gStyle->SetPalette ( 1 ) ; // set nicer colors than default
  // gStyle->SetOptTitle ( 0 ) ; // suppress title box
  gStyle->SetTitleSize( 0.025 , "t" );
  gStyle->SetTitleSize( 0.025 , "xyz" );

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



template< typename ... Functors >
void ToFile( const std::string& aFilename , Functors&&... aFunctors )
{
  gROOT->Reset ( ) ; // re−initialize ROOT
  gROOT->SetBatch( kTRUE );
  gROOT->SetStyle ( "Plain" ) ; // set empty TStyle ( nicer on paper )
  gStyle->SetOptStat ( 0 ) ; // print statistics on plots , ( 0 ) for no output
  gStyle->SetPalette ( 1 ) ; // set nicer colors than default
  // gStyle->SetOptTitle ( 0 ) ; // suppress title box
  gStyle->SetTitleSize( 0.025 , "t" );
  gStyle->SetTitleSize( 0.025 , "xyz" );

  constexpr std::size_t cnt( sizeof...(Functors) );

  TCanvas* c = new TCanvas( "c" , "" , 0 , 0 , 3200 , 1600 );

  std::size_t i = ceil( sqrt( cnt ) );
  std::size_t j = ceil( double(cnt) / i );

  c->Divide( i , j );
  DisplayHelper( c , 1 , std::forward<Functors>( aFunctors)... );

  c->SaveAs( aFilename.c_str() );
}