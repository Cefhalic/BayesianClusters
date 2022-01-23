/* ===== Cluster sources ===== */
#include "Cluster_Data.hpp"
#include "Cluster_GlobalVars.hpp"

/* ===== C++ ===== */
#include <iostream>
#include <algorithm>

/* ===== Local utilities ===== */
#include "ProgressBar.hpp"
#include "Vectorize.hpp"

/* ===== Function for loading data from CSV file ===== */
void __LoadCSV__( const std::string& aFilename , const Event& aEvent , std::vector< Data >& aData , const std::size_t& aOffset , int aCount )
{
  auto f = fopen( aFilename.c_str() , "rb");
  if (fseek(f, aOffset, SEEK_SET)) throw std::runtime_error( "Fseek failed" ); // seek to offset from start_point

  char ch[256];
  char* lPtr( ch );

  auto ReadUntil = [ & ]( const char& aChar ){
    lPtr = ch;
    while ( ( *lPtr = fgetc(f)) != EOF )
    {
      aCount--;
      if( *lPtr == aChar ) return;
      lPtr++;
    }
  };

  ReadUntil( '\n' ); // Throw away first line, or any partial lines (other thread will handle it)
  while( aCount > 0 )
  {
    ReadUntil( ',' ); //"id"
    if( *lPtr == EOF ) break;
    ReadUntil( ',' ); //"frame"
    ReadUntil( ',' ); //"x [nm]"
    double x = aEvent.toAlgorithmX( strtod( ch , &lPtr ) * nanometer );
    ReadUntil( ',' ); //"y [nm]"
    double y = aEvent.toAlgorithmY( strtod( ch , &lPtr ) * nanometer );      
    ReadUntil( ',' ); //"sigma [nm]"      
    ReadUntil( ',' ); //"intensity [photon]"
    ReadUntil( ',' ); //"offset [photon]"
    ReadUntil( ',' ); //"bkgstd [photon]"
    ReadUntil( ',' ); //"chi2"
    ReadUntil( '\n' ); //"uncertainty_xy [nm]"
    double s = Event::mParameters.toAlgorithmUnits( strtod( ch , &lPtr ) * nanometer );      

    if( fabs(x) < 1 and fabs(y) < 1 ) aData.emplace_back( x , y , s );
  }
  
  fclose(f);

  // BitonicSort< up , Data >( aData.begin() , aData.end() );
  std::sort( aData.begin() , aData.end() );
}


/* ===== Function for loading data from CSV file ===== */
void LoadCSV( const std::string& aFilename , Event& aEvent )
{
  auto f = fopen( aFilename.c_str() , "rb");
  if ( f == NULL ) throw std::runtime_error( "File is not available" );
  fseek(f, 0, SEEK_END); // seek to end of file
  auto lSize = ftell(f); // get current file pointer
  fclose(f);

  int lChunkSize = ceil( double(lSize) / (Concurrency+1) );
  std::vector< std::vector< Data > > lData( Concurrency+1 );

  ProgressBar2 lProgressBar( "Reading File" , lSize );
  for( std::size_t i(0); i!=Concurrency ; ++i ) ThreadPool.at(i)->submit( [ i , &lChunkSize , &lData , &aFilename , &aEvent ](){ __LoadCSV__( aFilename , aEvent , lData[i] , i*lChunkSize , lChunkSize ); } );
  WrappedThread::run_and_wait( [ &lChunkSize , &lData , &aFilename , &aEvent ](){ __LoadCSV__( aFilename , aEvent , lData[Concurrency] , Concurrency*lChunkSize , lChunkSize ); } );
  // WrappedThread::wait();

  std::size_t lSize2( 0 );
  for( auto& i : lData ) lSize2 += i.size();

  aEvent.mData.clear();
  aEvent.mData.reserve( lSize2 );

  for( auto& i : lData )
  {
    lSize2 = aEvent.mData.size();
    aEvent.mData.insert( aEvent.mData.end() , std::make_move_iterator( i.begin() ) , std::make_move_iterator( i.end() ) );
    i.erase( i.begin() , i.end() );
    std::inplace_merge ( aEvent.mData.begin() , aEvent.mData.begin()+lSize2 , aEvent.mData.end() );  
  }

  std::cout << "Read " << aEvent.mData.size() << " points" << std::endl;
}



void WriteCSV( const std::string& aFilename , const Event& aEvent )
{
  auto f = fopen( aFilename.c_str() , "w");
  if ( f == NULL ) throw std::runtime_error( "File is not available" );

  fprintf( f , "id,frame,x [nm],y [nm],sigma [nm],intensity [photon],offset [photon],bkgstd [photon],chi2,uncertainty_xy [nm]\n" );

  ProgressBar lProgressBar( "Writing File" , aEvent.mData.size() );
  for( auto& i : aEvent.mData ){
    fprintf( f , ",,%f,%f,,,,,,%f\n" , aEvent.toPhysicalX(i.x)/nanometer , aEvent.toPhysicalY(i.y)/nanometer , aEvent.mParameters.toPhysicalUnits(i.s)/nanometer );
    lProgressBar++;
  }

  fclose(f);
}

