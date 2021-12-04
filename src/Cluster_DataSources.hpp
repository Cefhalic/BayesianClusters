#pragma once

/* ===== Cluster sources ===== */
#include "Cluster_Data.hpp"

/* ===== C++ ===== */
#include <vector>
#include <iostream>

/* ===== For Root ===== */
#include "TRandom3.h"

/* ===== Local utilities ===== */
#include "ProgressBar.hpp"



/* ===== Utility function for creating a vector of data ===== */
std::vector< Data > CreatePseudoData( const int& aBackgroundCount , const int& aClusterCount , const int& aClusterSize , const double& aClusterScale )
{
  std::cout << "Generating Pseudodata" << std::endl;

  std::vector< Data > lData;
  lData.reserve( aBackgroundCount + ( aClusterCount * aClusterSize ) );

  TRandom3 lRand( 23423 );

  for( int i(0); i!= aBackgroundCount; ++i )
  {
    double x( lRand.Uniform( -1.0 , 1.0 ) ) , y( lRand.Uniform( -1.0 , 1.0 ) ) , s( lRand.Gaus( aClusterScale/10 , aClusterScale/30 ) );
    lData.emplace_back( x , y , s );
  }

  for( int i(0); i!= aClusterCount; ++i )
  {
    double x( lRand.Uniform( -1.0 , 1.0 ) ) , y( lRand.Uniform( -1.0 , 1.0 ) );
    double sigma( fabs( lRand.Gaus( aClusterScale , aClusterScale/3 ) ) );
    for( int j(0) ; j!= aClusterSize ; /* */ )
    {
      double x2( lRand.Gaus( x , sigma ) ) , y2( lRand.Gaus( y , sigma ) ) , s( lRand.Gaus( aClusterScale/10 , aClusterScale/30 ) );  
      if( x2 > 1 or x2 < -1 or y2 > 1 or y2 < -1 ) continue;    
      lData.emplace_back( x2 , y2 , s );
      ++j;
    }
  }

  return lData;
}


/* ===== Function for loading data from CSV file ===== */
std::vector< Data > LoadCSV( const std::string& aFilename , const double& m , const double& c_x , const double& c_y )
{
  auto f = fopen( aFilename.c_str() , "r");
  if ( f == NULL ) throw std::runtime_error( "File is not available" );

  fseek(f, 0, SEEK_END); // seek to end of file
  auto lSize = ftell(f); // get current file pointer
  fseek(f, 0, SEEK_SET); // seek back to beginning of file

  std::vector< Data > lData;
  lData.reserve( 3e6 );

  {
    char ch[256];
    char* lPtr( ch );
    ProgressBar lProgressBar( "Reading File" , lSize );

    auto ReadUntil = [ &ch , &f , &lPtr , &lProgressBar ]( const char& aChar ){
      lPtr = ch;
      while ( ( *lPtr = fgetc(f)) != EOF )
      {
        lProgressBar++;
        if( *lPtr == aChar ) return;
        lPtr++;
      }
    };

    ReadUntil( '\n' ); // Throw away first line
    while( true )
    {
      ReadUntil( ',' ); //"id"
      if( *lPtr == EOF ) break;
      ReadUntil( ',' ); //"frame"
      ReadUntil( ',' ); //"x [nm]"
      double x = m*(strtod( ch , &lPtr ) - c_x);
      ReadUntil( ',' ); //"y [nm]"
      double y = m*(strtod( ch , &lPtr ) - c_y);      
      ReadUntil( ',' ); //"sigma [nm]"      
      ReadUntil( ',' ); //"intensity [photon]"
      ReadUntil( ',' ); //"offset [photon]"
      ReadUntil( ',' ); //"bkgstd [photon]"
      ReadUntil( ',' ); //"chi2"
      ReadUntil( '\n' ); //"uncertainty_xy [nm]"
      double s = m * strtod( ch , &lPtr );      

      if( fabs(x) < 1 and fabs(y) < 1 ) lData.emplace_back( x , y , s );
    }
  }
  fclose(f);

  std::cout << "Read " << lData.size() << " points" << std::endl;

  return lData;
}
