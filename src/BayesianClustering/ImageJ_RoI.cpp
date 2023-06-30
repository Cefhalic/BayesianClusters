//! \file ImageJ_RoI.cpp

/* ===== C++ ===== */
#include <iostream>
#include <stdlib.h>
#include <cstring>
#include <zip.h>
#include <stdexcept>

/* ===== Cluster sources ===== */
#include "BayesianClustering/ImageJ_RoI.hpp"


roi_polygon DecodeBinaryRoI( const uint8_t* const aData )
{
  const auto GetShort = []( const uint8_t* const aData , const std::size_t& aIndex ){ return ( uint16_t( aData[aIndex] ) << 8 ) | uint16_t( aData[aIndex+1] ); };

  if( strncmp( (char*) aData , "Iout" , 4 ) ) throw std::runtime_error( "Not an ImageJ RoI file" );
  if( GetShort( aData , 6 ) != 0x00 ) throw std::runtime_error( "Only polygon data-type supported" ); 
  if( GetShort( aData , 4 ) != 0xE4 ) throw std::runtime_error( "Only tested on version 0xE4" ); 
  roi_polygon lRing;
  uint16_t Top( GetShort( aData , 8 ) ), Left( GetShort( aData , 10 ) ), lCount( GetShort( aData , 16 ) );
  for( std::size_t i(0) , i0( 64 ) , i1( 64 + (2*lCount) ) ; i!=lCount ; ++i , i0+=2 , i1+=2 ) boost::geometry::append( lRing , roi_point( GetShort( aData , i0 ) + Left , GetShort( aData , i1 ) + Top ) );    
  return lRing;
}


std::map< std::string , roi_polygon > OpenRoiZipfile( const std::string& aZipFileName )
{   
  std::map< std::string , roi_polygon > lRet;

  // -------------------------------------------------------------------------------
  int lError(0);
  zip_t* lArchive = zip_open( aZipFileName.c_str() , 0 , &lError ); // initializes a pointer to a zip archive , sets that pointer to the zip file
  struct zip_stat* lMemberInfo = (struct zip_stat*) malloc( 1024 ); // the zip_stat structure contains information such as file name, size, compressed size
  zip_stat_init( lMemberInfo );                                     // create and "initialize" the structure according to documentation

  uint8_t* lBuffer( NULL );
  for ( std::size_t lIndex( 0 ); !zip_stat_index( lArchive , lIndex , 0 , lMemberInfo ) ; ++lIndex ) { // we open the file at the lIndex'th index inside the archive we loop and print every file and its contents, stopping when zip_stat_index did not return 0
    lBuffer = (uint8_t*) realloc( lBuffer , lMemberInfo->size + 1 );   // allocate room for the entire file contents
    zip_file_t* lMemberPtr = zip_fopen_index( lArchive , lIndex , 0 ); // opens file at lIndex index 
    zip_fread( lMemberPtr , lBuffer , lMemberInfo->size );             // reads lMemberInfo->size bytes from lMemberPtr into lBuffer buffer 
    lRet[ lMemberInfo->name ] = DecodeBinaryRoI( lBuffer );  
  } 
  free( lBuffer ); // free allocated buffer
  // -------------------------------------------------------------------------------

  return lRet;
}