#pragma once

/* ===== Cluster sources ===== */
#include "Cluster_Data.hpp"

/* ===== C++ ===== */
#include <vector>


/* ===== Utility function for creating a vector of data ===== */
std::vector< Data > CreatePseudoData( const int& aBackgroundCount , const int& aClusterCount , const int& aClusterSize , const double& aClusterScale );




// /* ===== Function for loading data from CSV file ===== */
// void __LoadCSV__( const std::string& aFilename , const double& c_x , const double& c_y , std::vector< Data >& aData , const std::size_t& aOffset , int aCount )
// {
//   auto f = fopen( aFilename.c_str() , "r");
//   fseek(f, aOffset, SEEK_SET); // seek to start_point

//   char ch[256];
//   char* lPtr( ch );

//   auto ReadUntil = [ & ]( const char& aChar ){
//     lPtr = ch;
//     while ( ( *lPtr = fgetc(f)) != EOF )
//     {
//       if( *lPtr == aChar ) return;
//       lPtr++;
//       aCount--;
//     }
//   };

//   ReadUntil( '\n' ); // Throw away first line, or any partial lines (other thread will handle it)
//   while( aCount > 0 )
//   {
//     ReadUntil( ',' ); //"id"
//     if( *lPtr == EOF ) break;
//     ReadUntil( ',' ); //"frame"
//     ReadUntil( ',' ); //"x [nm]"
//     double x = Parameters.scale() * ( (strtod( ch , &lPtr ) * nanometer ) - c_x);
//     ReadUntil( ',' ); //"y [nm]"
//     double y = Parameters.scale() * ( (strtod( ch , &lPtr ) * nanometer ) - c_y);      
//     ReadUntil( ',' ); //"sigma [nm]"      
//     ReadUntil( ',' ); //"intensity [photon]"
//     ReadUntil( ',' ); //"offset [photon]"
//     ReadUntil( ',' ); //"bkgstd [photon]"
//     ReadUntil( ',' ); //"chi2"
//     ReadUntil( '\n' ); //"uncertainty_xy [nm]"
//     double s = Parameters.scale() * ( strtod( ch , &lPtr ) * nanometer );      

//     if( fabs(x) < 1 and fabs(y) < 1 ) aData.emplace_back( x , y , s );
//   }
  
//   fclose(f);

//   // std::sort( aData );
// }


// /* ===== Function for loading data from CSV file ===== */
// std::vector< Data > LoadCSV( const std::string& aFilename , const double& c_x , const double& c_y )
// {
//   auto f = fopen( aFilename.c_str() , "r");
//   if ( f == NULL ) throw std::runtime_error( "File is not available" );
//   fseek(f, 0, SEEK_END); // seek to end of file
//   auto lSize = ftell(f); // get current file pointer
//   fclose(f);


//   std::size_t lChunkSize = ceil( double(lSize) / Concurrency );
//   std::vector< std::vector< Data > > lData( Concurrency );

// {
//   ProgressBar2 lProgressBar( "Reading File" , lSize );
//   for( std::size_t i(0); i!=Concurrency ; ++i ) ThreadPool.at(i)->submit( [ i , &lChunkSize , &lData , &aFilename , &c_x , &c_y ](){ __LoadCSV__( aFilename , c_x , c_y , lData[i] , i*lChunkSize , lChunkSize ); } );
// }

//   for( auto& i : lData ) std::cout << i.size() << std::endl;
//   //std::cout << "Read " << lData.size() << " points" << std::endl;

//   std::vector< Data > lData2;
//   return lData2;
// }




/* ===== Function for loading data from CSV file ===== */
std::vector< Data > LoadCSV( const std::string& aFilename , const double& c_x , const double& c_y );