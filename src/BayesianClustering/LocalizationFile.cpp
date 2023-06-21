//! \file LocalizationFile.cpp

/* ===== C++ ===== */
#include <iostream>

/* ===== Cluster sources ===== */
#include "BayesianClustering/LocalizationFile.hpp"
#include "BayesianClustering/Data.hpp"
#include "BayesianClustering/RoI.hpp"

/* ===== Local utilities ===== */
#include "Utilities/ProgressBar.hpp"
#include "Utilities/Vectorize.hpp"
#include "Utilities/Units.hpp"

// #include "Utilities/MemoryMonitoring.hpp"

//! Multithreading handler for loading a chunk of data from CSV file
//! \param aFilename The name of the file to open
//! \param aData     A vector into which to fill data
//! \param aOffset   The offset into the file
//! \param aCount    The (approximate) number of bytes to be handled by this handler
void __LoadCSV__( const std::string& aFilename, std::vector< Data >& aData, const std::size_t& aOffset, int aCount )
{
  aData.reserve( aCount / 50.0 ); // The minimum line-length appears to be mid-50's bytes long, so by reserving this many entries, we should never have to reallocate

  auto f = fopen( aFilename.c_str(), "rb");
  if (fseek(f, aOffset, SEEK_SET)) throw std::runtime_error( "Fseek failed" ); // seek to offset from start_point

  char ch[256];
  char* lPtr( ch );

  auto ReadUntil = [ & ]( const char& aChar ) {
    lPtr = ch;
    while ( ( *lPtr = fgetc(f)) != EOF ) {
      aCount--;
      if( *lPtr == aChar ) return;
      lPtr++;
    }
  };

  ReadUntil( '\n' ); // Throw away first line, or any partial lines (other thread will handle it)
  while( aCount > 0 ) {
    ReadUntil( ',' ); //"id"
    if( *lPtr == EOF ) break;
    ReadUntil( ',' ); //"frame"
    ReadUntil( ',' ); //"x [nm]"
    double x = strtod( ch, &lPtr ) * nanometer;
    ReadUntil( ',' ); //"y [nm]"
    double y = strtod( ch, &lPtr ) * nanometer;
    ReadUntil( ',' ); //"sigma [nm]"
    double sigma = strtod( ch, &lPtr );
    ReadUntil( ',' ); //"intensity [photon]"
    ReadUntil( ',' ); //"offset [photon]"
    ReadUntil( ',' ); //"bkgstd [photon]"
    ReadUntil( ',' ); //"chi2"
    ReadUntil( '\n' ); //"uncertainty_xy [nm]"
    double s = strtod( ch, &lPtr ) * nanometer;
    if ( ( sigma < 100 ) or ( sigma  > 300) ) continue;
    aData.emplace_back( x, y, s );
  }

  fclose(f);
}


// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
LocalizationFile::LocalizationFile( const std::string& aFilename )
{
  if( aFilename.size() == 0 ) throw std::runtime_error( "No input file specified" );

  auto f = fopen( aFilename.c_str(), "rb");
  if ( f == NULL ) throw std::runtime_error( "File is not available" );
  fseek(f, 0, SEEK_END); // seek to end of file
  auto lSize = ftell(f); // get current file pointer
  fclose(f);

  int lChunkSize = ceil( double(lSize) / Nthreads );
  std::vector< std::vector< Data > > lData( Nthreads );

  {
    ProgressTimer lProgressTimer( "Reading File" );
    [ & ]( const std::size_t& i ) { __LoadCSV__( aFilename, lData[i], i*lChunkSize, lChunkSize ); } && range( Nthreads );

    std::size_t lSize2( 0 );
    for( auto& i : lData ) lSize2 += i.size();
    mData.reserve( lSize2 );

    for( auto& i : lData ) {
      mData.insert( mData.end(), std::make_move_iterator( i.begin() ), std::make_move_iterator( i.end() ) );
    }
  }

  std::cout << "Read " << mData.size() << " points" << std::endl;

}
// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//! Typedef an array for histogramming a Localization File
typedef std::array< std::array< int, 512 >, 512 > tArray;

//! Recursively search histogram for continuously connected regions over threshold
//! \param aHist The histogram being searched
//! \param aRoIid The id of the region being allocated
//! \param i The horizontal index of the current cell in the histogram
//! \param j The vertical index of the current cell in the histogram
void __RecursiveSearch__( tArray& aHist, const int& aRoIid, const int& i, const int& j )
{
  aHist[i][j] = aRoIid;

  if( i != 0   and aHist[i-1][j] < 0 ) __RecursiveSearch__( aHist, aRoIid, i-1, j );
  if( i != 511 and aHist[i+1][j] < 0 ) __RecursiveSearch__( aHist, aRoIid, i+1, j );
  if( j != 0   and aHist[i][j-1] < 0 ) __RecursiveSearch__( aHist, aRoIid, i, j-1 );
  if( j != 511 and aHist[i][j+1] < 0 ) __RecursiveSearch__( aHist, aRoIid, i, j+1 );
}

void LocalizationFile::ExtractRoIs( const std::function< void( RoI& ) >& aCallback ) const
{
  // Calculate our scaling factor
  std::pair< double, double> lXbound( std::make_pair( 9e99, -9e99 ) ), lYbound( std::make_pair( 9e99, -9e99 ) );

  for( auto& k : mData ) {
    if ( k.x < lXbound.first  ) lXbound.first  = k.x;
    if ( k.x > lXbound.second ) lXbound.second = k.x;
    if ( k.y < lYbound.first  ) lYbound.first  = k.y;
    if ( k.y > lYbound.second ) lYbound.second = k.y;
  }

  auto lXScale( 512.0 / ( lXbound.second - lXbound.first ) );
  auto lYScale( 512.0 / ( lYbound.second - lYbound.first ) );

  // Fill a histogram
  tArray lHist;
  for( int i(0) ; i!=512 ; ++i ) for( int j(0) ; j!=512 ; ++j ) lHist[i][j] = 0; // Ensure array is zero'ed to start with

  for( auto& k : mData ) {
    std::size_t x = ( k.x - lXbound.first ) * lXScale;
    std::size_t y = ( k.y - lYbound.first ) * lYScale;
    lHist[x][y] += 1;
  }

  //Apply a Gaussian Blur
  constexpr int cGaussianSize( 7 );
  constexpr int cMaskSize( (2*cGaussianSize)+1 );
  constexpr double cScaling( -1.0 / (cMaskSize*cMaskSize) );
  std::array< std::array< float, cMaskSize >, cMaskSize > lMask;
  for( int i(-cGaussianSize) ; i<=cGaussianSize ; ++i ) for( int j(-cGaussianSize) ; j<=cGaussianSize ; ++j ) {
      lMask[i+cGaussianSize][j+cGaussianSize] = exp( cScaling * ( (i*i)+(j*j) ) );
    }

  tArray lHist2;
  int Max( INT_MIN );
  for( int i(0) ; i!=512 ; ++i ) for( int j(0) ; j!=512 ; ++j ) {
      lHist2[i][j] = 0;
      for( int k(-cGaussianSize) ; k<=cGaussianSize ; ++k ) for( int l(-cGaussianSize) ; l<=cGaussianSize ; ++l ) {
          auto i2( i+k ), j2( j+l );
          if( i2 < 0 or i2 > 511 or j2 < 0 or j2 > 511 ) continue;
          lHist2[i][j] += ( lHist[i2][j2] * lMask[k+cGaussianSize][l+cGaussianSize] );
        }
      if ( Max < lHist2[i][j] ) Max = lHist2[i][j];
    }

  // Threshold the histogram
  int lThreshold = 0.2 * Max;
  for( int i(0) ; i!=512 ; ++i ) for( int j(0) ; j!=512 ; ++j ) {
      lHist2[i][j] = - int( lHist2[i][j] > lThreshold );
    }

  // Depth-first search for continuously connected regions
  int lRoIid( 0 );
  for( int i(0) ; i!=512 ; ++i ) for( int j(0) ; j!=512 ; ++j ) {
      if( lHist2[i][j] < 0 ) __RecursiveSearch__( lHist2, ++lRoIid, i, j );
    }


  // Local record to store the size, the bounds and the datapoints
  struct tRecord {
    std::size_t Id;
    std::size_t Size;
    std::pair< double, double > X, Y;
    std::vector< const Data* > Ptrs;
  };

  std::vector < tRecord > lRecords( lRoIid+1, { 0 , 0, std::make_pair( 9e99, -9e99 ), std::make_pair( 9e99, -9e99 ), std::vector< const Data* >() } );

  // Fill in holes and calculate the area
  for( int i(0) ; i!=512 ; ++i ) for( int j(0) ; j!=512 ; ++j ) {
    if( lHist2[i][j] ) {
      lRecords[ lHist2[i][j] ].Size++;
      continue;
    }

    int L(0), R(0), U(0), D(0);

    for( int k(i-1) ; k>=0 ; --k )
      if( (L = lHist2[k][j]) ) break;

    for( int k(i+1) ; k!=512 ; ++k )
      if( (R = lHist2[k][j]) ) break;

    for( int k(j-1) ; k>=0 ; --k )
      if( (D = lHist2[i][k]) ) break;

    for( int k(j+1) ; k!=512 ; ++k )
      if( (U = lHist2[i][k]) ) break;

    if( U == L and R == L and D == L ) lRecords[ lHist2[i][j] = L ].Size++;
  }

  // Find the bounds and the datapoints
  for( auto& k : mData ) {
    std::size_t x = ( k.x - lXbound.first ) * lXScale;
    std::size_t y = ( k.y - lYbound.first ) * lYScale;

    auto& Id = lHist2[x][y];
    if( !Id ) continue;

    auto& lRecord = lRecords[Id];

    if( lRecord.Size < 500 ) continue; // Cull micro RoIs

    lRecord.Id = Id;
    lRecord.Ptrs.push_back( &k );
    if ( k.x < lRecord.X.first  ) lRecord.X.first  = k.x;
    if ( k.x > lRecord.X.second ) lRecord.X.second = k.x;
    if ( k.y < lRecord.Y.first  ) lRecord.Y.first  = k.y;
    if ( k.y > lRecord.Y.second ) lRecord.Y.second = k.y;
  }

  std::sort( lRecords.begin() , lRecords.end() , []( const tRecord& a , const tRecord& b ){ return a.Ptrs.size() < b.Ptrs.size(); } );

  for( auto& lRecord : lRecords ) {
    if( !lRecord.Ptrs.size() ) continue;

    // double vm_usage, resident_set;
    // mem_usage( vm_usage , resident_set );    
    // std::cout << vm_usage << " " << resident_set << std::endl;

    double lCentreX( ( lRecord.X.second + lRecord.X.first ) / 2.0 ), lCentreY( ( lRecord.Y.second + lRecord.Y.first ) / 2.0 );
    double lWidthX( lRecord.X.second - lRecord.X.first ), lWidthY( lRecord.Y.second - lRecord.Y.first );

    std::vector< Data > lData;

    for( auto& k : lRecord.Ptrs ) {
      double x = k->x - lCentreX;
      double y = k->y - lCentreY;
      lData.emplace_back( x, y, k->s );
    }

    RoI lRoI( std::move( lData ) );
    lRoI.SetCentre( lCentreX, lCentreY );
    lRoI.SetWidth( lWidthX, lWidthY );

    aCallback( lRoI );

    // int fd = open("/proc/sys/vm/drop_caches", O_WRONLY);
    // write(fd, "3", 1 );
    // fsync(fd);
    // close(fd);

  }

}


void LocalizationFile::ExtractRoIs( const ManualRoI& aRoI , const std::function< void( RoI& ) >& aCallback ) const
{
  auto lMaxX = aRoI.width / 2.0;
  auto lMaxY = aRoI.height / 2.0;

  std::vector< Data > lData;  
  for( auto& k : mData ) {
      double x = k.x - aRoI.x;
      double y = k.y - aRoI.y;
      if( fabs(x) < lMaxX and fabs(y) < lMaxY ) lData.emplace_back( x, y, k.s );    
  }

  RoI lRoI( std::move( lData ) );
  lRoI.SetCentre( aRoI.x , aRoI.y );
  lRoI.SetWidth( aRoI.width , aRoI.height );

  aCallback( lRoI );        
}
// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
