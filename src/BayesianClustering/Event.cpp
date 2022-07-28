
/* ===== Cluster sources ===== */
#include "BayesianClustering/Event.hpp"
#include "BayesianClustering/EventProxy.hpp"
#include "BayesianClustering/Configuration.hpp"

/* ===== Local utilities ===== */
#include "Utilities/ProgressBar.hpp"
#include "Utilities/Vectorize.hpp"

// /* ===== C++ ===== */
#include <iostream>


// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Configuration Configuration::Instance;

Event::Event()
{
  const std::string& lFilename = Configuration::Instance.inputFile();
  if( lFilename.size() == 0 ) throw std::runtime_error( "No input file specified" ); 

  LoadCSV( lFilename );
}

void Event::Preprocess()
{
  // Populate mNeighbour lists  
  ProgressBar2 lProgressBar( "Populating neighbourhood" , mData.size() );
  [&]( const std::size_t& i ){ mData.at( i ).Preprocess( mData , i ); } || range( mData.size() );  // Interleave threading since processing time increases with radius from origin
}

void Event::ScanRT( const std::function< void( const EventProxy& , const double& , const double& , std::vector<uint32_t>  ) >& aCallback ) 
{
  Preprocess();    

  {
    ProgressBar2 lProgressBar( "Populating localization scores" , mData.size() );
    [&]( const std::size_t& i ){ mData.at( i ).PreprocessLocalizationScores( mData ); } || range( mData.size() );  // Interleave threading since processing time increases with radius from origin
  }

  std::vector< EventProxy > lEventProxys;
  lEventProxys.reserve( Nthreads );
  for( int i(0) ; i!=Nthreads ; ++i ) lEventProxys.emplace_back( *this );
  ProgressBar2 lProgressBar( "Scan over RT"  , 0 );
  [&]( const std::size_t& i ){ lEventProxys.at(i).ScanRT( aCallback , Nthreads , i ); } || range( Nthreads );
}

void Event::Clusterize( const double& R , const double& T , const std::function< void( const EventProxy& ) >& aCallback )
{
  Preprocess();    

  EventProxy lProxy( *this );
  ProgressBar2 lProgressBar( "Clusterize"  , 0 );
  lProxy.Clusterize( R ,  T , aCallback );
}

/* ===== Function for loading a chunk of data from CSV file ===== */
void __LoadCSV__( const std::string& aFilename , Event& aEvent , std::vector< Data >& aData , const std::size_t& aOffset , int aCount )
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
    double x = Configuration::Instance.toAlgorithmX( strtod( ch , &lPtr ) * nanometer );
    ReadUntil( ',' ); //"y [nm]"
    double y = Configuration::Instance.toAlgorithmY( strtod( ch , &lPtr ) * nanometer );      
    ReadUntil( ',' ); //"sigma [nm]"      
    ReadUntil( ',' ); //"intensity [photon]"
    ReadUntil( ',' ); //"offset [photon]"
    ReadUntil( ',' ); //"bkgstd [photon]"
    ReadUntil( ',' ); //"chi2"
    ReadUntil( '\n' ); //"uncertainty_xy [nm]"
    double s = Configuration::Instance.toAlgorithmUnits( strtod( ch , &lPtr ) * nanometer );      

    if( fabs(x) < 1 and fabs(y) < 1 ) aData.emplace_back( x , y , s );
  }
  
  fclose(f);

  std::sort( aData.begin() , aData.end() );
}

void Event::LoadCSV( const std::string& aFilename )
{
  auto f = fopen( aFilename.c_str() , "rb");
  if ( f == NULL ) throw std::runtime_error( "File is not available" );
  fseek(f, 0, SEEK_END); // seek to end of file
  auto lSize = ftell(f); // get current file pointer
  fclose(f);

  int lChunkSize = ceil( double(lSize) / Nthreads );
  std::vector< std::vector< Data > > lData( Nthreads );

  ProgressBar2 lProgressBar( "Reading File" , lSize );
  [ & ]( const std::size_t& i ){ __LoadCSV__( aFilename , *this , lData[i] , i*lChunkSize , lChunkSize ); } && range( Nthreads );

  std::size_t lSize2( 0 );
  for( auto& i : lData ) lSize2 += i.size();

  mData.reserve( lSize2 );

  for( auto& i : lData )
  {
    lSize2 = mData.size();
    mData.insert( mData.end() , std::make_move_iterator( i.begin() ) , std::make_move_iterator( i.end() ) );
    i.erase( i.begin() , i.end() );
    std::inplace_merge ( mData.begin() , mData.begin()+lSize2 , mData.end() );  
  }

  // // Validate the in-place merge
  // PRECISION last = -1;
  // for( auto& i : mData )
  // {
  //   if( i.r < last ) throw std::runtime_error( "Out of order!" );
  //   last = i.r;
  // }

  std::cout << "Read " << mData.size() << " points" << std::endl;
}

void Event::WriteCSV( const std::string& aFilename )
{
  auto f = fopen( aFilename.c_str() , "w");
  if ( f == NULL ) throw std::runtime_error( "File is not available" );

  fprintf( f , "id,frame,x [nm],y [nm],sigma [nm],intensity [photon],offset [photon],bkgstd [photon],chi2,uncertainty_xy [nm]\n" );

  ProgressBar lProgressBar( "Writing File" , mData.size() );
  for( auto& i : mData ){
    fprintf( f , ",,%f,%f,,,,,,%f\n" , Configuration::Instance.toPhysicalX(i.x)/nanometer , Configuration::Instance.toPhysicalY(i.y)/nanometer , Configuration::Instance.toPhysicalUnits(i.s)/nanometer );
    lProgressBar++;
  }

  fclose(f);
}
// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
