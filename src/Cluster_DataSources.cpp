/* ===== Cluster sources ===== */
#include "Cluster_Data.hpp"
#include "Cluster_GlobalVars.hpp"

/* ===== C++ ===== */
#include <iostream>
#include <algorithm>

/* ===== For Root ===== */
#include "TRandom3.h"

/* ===== Local utilities ===== */
#include "ProgressBar.hpp"
#include "Vectorize.hpp"




// enum sort_direction { up, down };

// template < sort_direction aDir , typename T >
// void BitonicSort( const typename std::vector<T>::iterator& aDataStart, const typename std::vector<T>::iterator& aDataEnd );

// template < sort_direction aDir , typename T >
// void BitonicMerge( const typename std::vector<T>::iterator& aDataStart, const typename std::vector<T>::iterator& aDataEnd);


// template < sort_direction aDir , typename T >
// void BitonicSort( const typename std::vector<T>::iterator& aDataStart, const typename std::vector<T>::iterator& aDataEnd) {
//   uint32_t lSize(aDataEnd - aDataStart);
//   if (lSize > 1) {
//     typename std::vector<T>::iterator lMidpoint(aDataStart + (lSize >> 1));
//     if (aDir == down) {
//       BitonicSort<up,T>(aDataStart, lMidpoint);
//       BitonicSort<down,T>( lMidpoint, aDataEnd);
//     } else {
//       BitonicSort<down,T>(aDataStart, lMidpoint);
//       BitonicSort<up,T>( lMidpoint, aDataEnd);
//     }
//     BitonicMerge<aDir,T>(aDataStart, aDataEnd);
//   }
// }

// template < sort_direction aDir , typename T >
// void BitonicMerge( const typename std::vector<T>::iterator& aDataStart, const typename std::vector<T>::iterator& aDataEnd) {
//   uint32_t lSize(aDataEnd - aDataStart);
//   if (lSize > 1) {
//     uint32_t lPower2(1);
//     while (lPower2 < lSize) lPower2 <<= 1;

//     typename std::vector<T>::iterator lMidpoint(aDataStart + (lPower2 >> 1)) , lFirst(aDataStart) , lSecond(lMidpoint);

//     for (; lSecond != aDataEnd; ++lFirst, ++lSecond) {
//       if (((*lSecond) < (*lFirst)) == (aDir == up)) {
//         std::swap(*lFirst, *lSecond);
//       }
//     }

//     BitonicMerge<aDir,T>( aDataStart, lMidpoint);
//     BitonicMerge<aDir,T>( lMidpoint, aDataEnd);
//   }
// }





/* ===== Utility function for creating a vector of data ===== */
std::vector< Data > CreatePseudoData( const int& aBackgroundCount , const int& aClusterCount , const int& aClusterSize , const double& aClusterScale )
{
  std::cout << "Generating Pseudodata" << std::endl;

  const auto lClusterScale( aClusterScale * Parameters.scale() );

  std::vector< Data > lData;
  lData.reserve( aBackgroundCount + ( aClusterCount * aClusterSize ) );

  TRandom3 lRand( 2345234534 );

  for( int i(0); i!= aBackgroundCount; ++i )
  {
    double x( lRand.Uniform( -1.0 , 1.0 ) ) , y( lRand.Uniform( -1.0 , 1.0 ) ) , s( lRand.Gaus( lClusterScale/10 , lClusterScale/30 ) );
    lData.emplace_back( x , y , s );
  }

  for( int i(0); i!= aClusterCount; ++i )
  {
    double x( lRand.Uniform( -1.0 , 1.0 ) ) , y( lRand.Uniform( -1.0 , 1.0 ) );
    double sigma( fabs( lRand.Gaus( lClusterScale , lClusterScale/3 ) ) );
    for( int j(0) ; j!= aClusterSize ; /* */ )
    {
      double x2( lRand.Gaus( x , sigma ) ) , y2( lRand.Gaus( y , sigma ) ) , s( lRand.Gaus( lClusterScale/10 , lClusterScale/30 ) );  
      if( x2 > 1 or x2 < -1 or y2 > 1 or y2 < -1 ) continue;    
      lData.emplace_back( x2 , y2 , s );
      ++j;
    }
  }

  std::sort( lData.begin() , lData.end() );
  return lData;
}



// #include <sstream>

/* ===== Function for loading data from CSV file ===== */
void __LoadCSV__( const std::string& aFilename , const double& c_x , const double& c_y , std::vector< Data >& aData , const std::size_t& aOffset , int aCount )
{
  auto f = fopen( aFilename.c_str() , "rb");
  if (fseek(f, aOffset, SEEK_SET)) throw std::runtime_error( "Fseek failed" ); // seek to offset from start_point

  char ch[256];
  char* lPtr( ch );

  // std::stringstream lStr;
  // lStr << std::this_thread::get_id() << " =>\n";

  auto ReadUntil = [ & ]( const char& aChar ){
    lPtr = ch;
    while ( ( *lPtr = fgetc(f)) != EOF )
    {
      // lStr << *lPtr;
      aCount--;
      if( *lPtr == aChar ) return;
      lPtr++;
    }
  };

  // lStr << "Head : ";
  ReadUntil( '\n' ); // Throw away first line, or any partial lines (other thread will handle it)
  while( aCount > 0 )
  {
    // lStr << "Body : ";
    ReadUntil( ',' ); //"id"
    if( *lPtr == EOF ) break;
    ReadUntil( ',' ); //"frame"
    ReadUntil( ',' ); //"x [nm]"
    double x = Parameters.scale() * ( (strtod( ch , &lPtr ) * nanometer ) - c_x);
    ReadUntil( ',' ); //"y [nm]"
    double y = Parameters.scale() * ( (strtod( ch , &lPtr ) * nanometer ) - c_y);      
    ReadUntil( ',' ); //"sigma [nm]"      
    ReadUntil( ',' ); //"intensity [photon]"
    ReadUntil( ',' ); //"offset [photon]"
    ReadUntil( ',' ); //"bkgstd [photon]"
    ReadUntil( ',' ); //"chi2"
    ReadUntil( '\n' ); //"uncertainty_xy [nm]"
    double s = Parameters.scale() * ( strtod( ch , &lPtr ) * nanometer );      

    if( fabs(x) < 1 and fabs(y) < 1 ) aData.emplace_back( x , y , s );
  }
  
  // lStr << " \n\n";
  // std ::cout << lStr.str() << std::flush;
  fclose(f);

  // BitonicSort< up , Data >( aData.begin() , aData.end() );
  std::sort( aData.begin() , aData.end() );
}


/* ===== Function for loading data from CSV file ===== */
std::vector< Data > LoadCSV( const std::string& aFilename , const double& c_x , const double& c_y )
{
  auto f = fopen( aFilename.c_str() , "rb");
  if ( f == NULL ) throw std::runtime_error( "File is not available" );
  fseek(f, 0, SEEK_END); // seek to end of file
  auto lSize = ftell(f); // get current file pointer
  fclose(f);

  int lChunkSize = ceil( double(lSize) / Concurrency );
  std::vector< std::vector< Data > > lData( Concurrency );

  ProgressBar2 lProgressBar( "Reading File" , lSize );
  for( std::size_t i(0); i!=Concurrency ; ++i ) ThreadPool.at(i)->submit( [ i , &lChunkSize , &lData , &aFilename , &c_x , &c_y ](){ __LoadCSV__( aFilename , c_x , c_y , lData[i] , i*lChunkSize , lChunkSize ); } );
  WrappedThread::wait();

  std::size_t lSize2( 0 );
  for( auto& i : lData ) lSize2 += i.size();

  std::vector< Data > lData2;
  lData2.reserve( lSize2 );

  for( auto& i : lData )
  {
    lSize2 = lData2.size();
    lData2.insert( lData2.end() , std::make_move_iterator( i.begin() ) , std::make_move_iterator( i.end() ) );
    i.erase( i.begin() , i.end() );
    std::inplace_merge ( lData2.begin() , lData2.begin()+lSize2 , lData2.end() );  
  }

  std::cout << "Read " << lData2.size() << " points" << std::endl;

  return lData2;
}




// /* ===== Function for loading data from CSV file ===== */
// std::vector< Data > LoadCSV( const std::string& aFilename , const double& c_x , const double& c_y )
// {
//   auto f = fopen( aFilename.c_str() , "r");
//   if ( f == NULL ) throw std::runtime_error( "File is not available" );

//   fseek(f, 0, SEEK_END); // seek to end of file
//   auto lSize = ftell(f); // get current file pointer
//   fseek(f, 0, SEEK_SET); // seek back to beginning of file

//   std::vector< Data > lData;
//   lData.reserve( 3e6 );

//   {
//     char ch[256];
//     char* lPtr( ch );
//     ProgressBar lProgressBar( "Reading File" , lSize );

//     auto ReadUntil = [ &ch , &f , &lPtr , &lProgressBar ]( const char& aChar ){
//       lPtr = ch;
//       while ( ( *lPtr = fgetc(f)) != EOF )
//       {
//         lProgressBar++;
//         if( *lPtr == aChar ) return;
//         lPtr++;
//       }
//     };

//     ReadUntil( '\n' ); // Throw away first line
//     while( true )
//     {
//       ReadUntil( ',' ); //"id"
//       if( *lPtr == EOF ) break;
//       std::size_t i = strtoul( ch , &lPtr , 10 );
//       ReadUntil( ',' ); //"frame"
//       ReadUntil( ',' ); //"x [nm]"
//       double x = Parameters.scale() * ( (strtod( ch , &lPtr ) * nanometer ) - c_x);
//       ReadUntil( ',' ); //"y [nm]"
//       double y = Parameters.scale() * ( (strtod( ch , &lPtr ) * nanometer ) - c_y);      
//       ReadUntil( ',' ); //"sigma [nm]"      
//       ReadUntil( ',' ); //"intensity [photon]"
//       ReadUntil( ',' ); //"offset [photon]"
//       ReadUntil( ',' ); //"bkgstd [photon]"
//       ReadUntil( ',' ); //"chi2"
//       ReadUntil( '\n' ); //"uncertainty_xy [nm]"
//       double s = Parameters.scale() * ( strtod( ch , &lPtr ) * nanometer );      

//       if( fabs(x) < 1 and fabs(y) < 1 ) lData.emplace_back( i , x , y , s );
//     }
//   }
//   fclose(f);

//   std::sort( lData.begin() , lData.end() );
  
//   std::cout << "Read " << lData.size() << " points" << std::endl;

//   return lData;
// }
