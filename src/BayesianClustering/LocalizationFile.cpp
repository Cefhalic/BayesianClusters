//! \file LocalizationFile.cpp

/* ===== C++ ===== */
#include <iostream>

/* ===== Cluster sources ===== */
#include "BayesianClustering/LocalizationFile.hpp"
#include "BayesianClustering/Data.hpp"
#include "BayesianClustering/RoI.hpp"
#include "BayesianClustering/ImageJ_RoI.hpp"

/* ===== Local utilities ===== */
#include "Utilities/ProgressBar.hpp"
#include "Utilities/Vectorize.hpp"
#include "Utilities/Units.hpp"
// #include "Utilities/MemoryMonitoring.hpp"

/* ===== BOOST C++ ===== */
#include <boost/filesystem.hpp>
#include <boost/geometry.hpp>
#include <boost/gil/extension/io/bmp.hpp>
#include <boost/gil/image.hpp>

enum class LocalizationTableType { ThunderSTORM, X_Y_Index, Unknown };

// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
LocalizationTableType getLocalizationTableType(const std::string& aFilename)
{
    auto f = fopen(aFilename.c_str(), "rb");
    if (fseek(f, 0, SEEK_SET))
        throw std::runtime_error("Fseek failed");

    char cbuf[256];
    fgets(cbuf, sizeof(cbuf), f);

    std::string header(cbuf), 
                delimiters = ",\t\r\n";
    
    std::vector<std::string> caption(1, "Caption");
    boost::split(caption, header, boost::is_any_of(delimiters));
    
    if (caption.size() >= 9) 
    {   // ThunderSTORM
        if (boost::algorithm::contains(caption[0], "id") &&
            boost::algorithm::contains(caption[1], "frame") &&
            boost::algorithm::contains(caption[2], "x [nm]") &&
            boost::algorithm::contains(caption[3], "y [nm]") &&
            boost::algorithm::contains(caption[4], "sigma [nm]") &&
            boost::algorithm::contains(caption[5], "intensity [photon]") &&
            boost::algorithm::contains(caption[6], "offset [photon]") &&
            boost::algorithm::contains(caption[7], "bkgstd [photon]") &&
            boost::algorithm::contains(caption[8], "uncertainty_xy [nm]"))
                return LocalizationTableType::ThunderSTORM;
    }    
    else if (caption.size() >= 3) 
    {   // X_Y_Index
        if (boost::algorithm::contains(caption[0], "x") &&
            boost::algorithm::contains(caption[1], "y") &&
            boost::algorithm::contains(caption[2], "index"))
                return LocalizationTableType::X_Y_Index;
    }

    return LocalizationTableType::Unknown;
}

//! Multithreading handler for loading a chunk of data from CSV file
//! \param aFilename The name of the file to open
//! \param aData     A vector into which to fill data
//! \param aOffset   The offset into the file
//! \param aCount    The (approximate) number of bytes to be handled by this handler
void __LoadCSV__(const std::string& aFilename, LocalizationTableType aFileType, std::vector< Data >& aData, const std::size_t& aOffset, int aCount)
{
    aData.reserve(aCount / 50.0); // The minimum line-length appears to be mid-50's bytes long, so by reserving this many entries, we should never have to reallocate

    auto f = fopen(aFilename.c_str(), "rb");
    if (fseek(f, aOffset, SEEK_SET)) throw std::runtime_error("Fseek failed"); // seek to offset from start_point

    char ch[256];
    char* lPtr(ch);

    auto ReadUntil = [&](const char& aChar) {
        lPtr = ch;
        while ((*lPtr = fgetc(f)) != EOF) {
            aCount--;
            if (*lPtr == aChar) return;
            lPtr++;
        }
    };

    ReadUntil('\n'); // Throw away first line, or any partial lines (other thread will handle it)

    switch (aFileType)
    {
        case LocalizationTableType::ThunderSTORM :
            while (aCount > 0) {
                ReadUntil(','); //"id"
                if (*lPtr == EOF) break;
                ReadUntil(','); //"frame"
                ReadUntil(','); //"x [nm]"
                double x = strtod(ch, &lPtr) * nanometer;
                ReadUntil(','); //"y [nm]"
                double y = strtod(ch, &lPtr) * nanometer;
                ReadUntil(','); //"sigma [nm]"
                double sigma = strtod(ch, &lPtr);
                ReadUntil(','); //"intensity [photon]"
                ReadUntil(','); //"offset [photon]"
                ReadUntil(','); //"bkgstd [photon]"
                ReadUntil(','); //"chi2"
                ReadUntil('\n'); //"uncertainty_xy [nm]"
                double s = strtod(ch, &lPtr) * nanometer;
                if ((sigma < 100) or (sigma > 300)) continue;
                aData.emplace_back(x, y, s);
            }
        case LocalizationTableType::X_Y_Index :
            while (aCount > 0) {
                ReadUntil(','); //"x"
                if (*lPtr == EOF) break;
                double x = strtod(ch, &lPtr) * nanometer;
                ReadUntil(','); //"y"
                double y = strtod(ch, &lPtr) * nanometer;
                ReadUntil('\n'); // ignore last value
                aData.emplace_back(x, y, 1.);
            }
        case LocalizationTableType::Unknown: ; // cannot happen
    }

    fclose(f);
}

// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
LocalizationFile::LocalizationFile( const std::string& aFilename ) : mFilename( aFilename )
{
  if( aFilename.size() == 0 ) throw std::runtime_error( "No input file specified" );

  auto f = fopen( aFilename.c_str(), "rb");
  if ( f == NULL ) throw std::runtime_error( "File is not available" );
  fseek(f, 0, SEEK_END); // seek to end of file
  auto lSize = ftell(f); // get current file pointer
  fclose(f);

  int lChunkSize = ceil( double(lSize) / Nthreads );
  std::vector< std::vector< Data > > lData( Nthreads );

  LocalizationTableType aFileType = getLocalizationTableType(aFilename);

  {
    ProgressTimer lProgressTimer( "Reading File" );
    [ & ]( const std::size_t& i ) { __LoadCSV__( aFilename, aFileType, lData[i], i*lChunkSize, lChunkSize ); } && range( Nthreads );

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
void LocalizationFile::ExtractRoIs( const ManualRoI& aRoI , const std::function< void( RoI& ) >& aCallback ) const
{
  auto lMaxX( aRoI.width / 2.0 ) , lMaxY( aRoI.height / 2.0 );

  std::vector< Data > lData;  
  for( auto& k : mData ) {
      double x = k.x - aRoI.x;
      double y = k.y - aRoI.y;
      if( fabs(x) < lMaxX and fabs(y) < lMaxY ) lData.emplace_back( x, y, k.s );    
  }

  RoI lRoI( "" , std::move( lData ) , aRoI.x , aRoI.y , aRoI.width * aRoI.height );
  aCallback( lRoI );        
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
  auto lBinArea( 1.0 / ( lXScale * lYScale ) );

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
  for( int i(-cGaussianSize) ; i<=cGaussianSize ; ++i ) for( int j(-cGaussianSize) ; j<=cGaussianSize ; ++j )
  {
    lMask[i+cGaussianSize][j+cGaussianSize] = exp( cScaling * ( (i*i)+(j*j) ) );
  }

  tArray lHist2;
  int Max( INT_MIN );
  for( int i(0) ; i!=512 ; ++i ) for( int j(0) ; j!=512 ; ++j ) {
      lHist2[i][j] = 0;
      for( int k(-cGaussianSize) ; k<=cGaussianSize ; ++k ) for( int l(-cGaussianSize) ; l<=cGaussianSize ; ++l )
      {
        auto i2( i+k ), j2( j+l );
        if( i2 < 0 or i2 > 511 or j2 < 0 or j2 > 511 ) continue;
        lHist2[i][j] += ( lHist[i2][j2] * lMask[k+cGaussianSize][l+cGaussianSize] );
      }
      if ( Max < lHist2[i][j] ) Max = lHist2[i][j];
    }

  // Threshold the histogram
  int lThreshold = 0.2 * Max;
  for( int i(0) ; i!=512 ; ++i ) for( int j(0) ; j!=512 ; ++j )
  {
    lHist2[i][j] = - int( lHist2[i][j] > lThreshold );
  }

  // Depth-first search for continuously connected regions
  int lRoIid( 0 );
  for( int i(0) ; i!=512 ; ++i ) for( int j(0) ; j!=512 ; ++j )
  {
    if( lHist2[i][j] < 0 ) __RecursiveSearch__( lHist2, ++lRoIid, i, j );
  }


  //! Local record to store the size, the bounds and the datapoints
  struct tRecord {
    std::size_t Size; //!< The number of histogram cells in the RoI
    double CentreX; //!< The mean X of the RoI
    double CentreY; //!< The mean Y of the RoI
    std::vector< const Data* > Ptrs; //!< The data points in the RoI
  };

  std::vector < tRecord > lRecords( lRoIid+1, { 0, 0.0 , 0.0, std::vector< const Data* >() } );

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

    if ( ( ( U == L ) or ( D == L ) ) and ( R == L ) )
    {
      lRecords[ lHist2[i][j] = L ].Size++;
      continue;
    } 

    if ( ( ( L == U ) or ( R == U ) ) and ( D == U ) ) 
    {
      lRecords[ lHist2[i][j] = U ].Size++;
      continue;
    } 
  }

  // // =============================================================================
  // // Export RoI's to BMP
  // {  
  //   namespace bg = boost::gil;
  //   auto img = bg::rgb8_image_t { 512 , 512 , bg::rgb8_pixel_t{255,255,255} };
  //   auto view = bg::view(img);

  //   auto ToRGB = []( const int& aVal )
  //   {
  //     uint32_t lVal(0) , lMask( 1 ) , lOutShift( 7 );
  //     for( int i(0) ; i!=24 ; ++i )
  //     {
  //       if( aVal & lMask ) lVal |= ( 1 << lOutShift );
  //       lMask <<= 1;
  //       lOutShift = ( lOutShift + 8 ) % 25;  
  //     }
  //     return bg::rgb8_pixel_t{ (uint8_t)(lVal>>0) , (uint8_t)(lVal>>8) , (uint8_t)(lVal>>16) };
  //   };

  //   std::map< int , bg::rgb8_pixel_t > lColours;

  //   for( int i(0) ; i!=512 ; ++i ) for( int j(0) ; j!=512 ; ++j ) {
  //     auto& Id = lHist2[i][j];
  //     if( !Id ) continue;
  //     auto& lRecord = lRecords[Id];
  //     if( lRecord.Size < 500 ) continue; // Cull micro RoIs

  //     auto lIt = lColours.find( Id );
  //     if( lIt == lColours.end() ) lIt = lColours.emplace( std::make_pair( Id , ToRGB( lColours.size() ) ) ).first;
  //     view( i , j ) = lIt->second;    
  //   }
  //   auto lOutFileName = boost::filesystem::path( mFilename ).stem().string() + ".AutoRoI.bmp";
  //   bg::write_view( lOutFileName , bg::const_view(img), bg::bmp_tag());
  // }
  // // =============================================================================


  // Find the centre and the datapoints
  for( auto& k : mData ) {
    std::size_t x( ( k.x - lXbound.first ) * lXScale ) , y( ( k.y - lYbound.first ) * lYScale );

    auto& Id = lHist2[x][y];
    if( !Id ) continue;

    auto& lRecord = lRecords[Id];
    if( lRecord.Size < 500 ) continue; // Cull micro RoIs
    lRecord.Ptrs.push_back( &k );
    lRecord.CentreX += k.x;
    lRecord.CentreY += k.y;    
  }

  std::sort( lRecords.begin() , lRecords.end() , []( const tRecord& a , const tRecord& b ){ return a.Ptrs.size() < b.Ptrs.size(); } );

  for( std::size_t i(1) ; i!=lRecords.size() ; ++i )
  {
    auto& lRecord = lRecords[i];    
    if( !lRecord.Ptrs.size() ) continue;

    lRecord.CentreX /= lRecord.Ptrs.size();
    lRecord.CentreY /= lRecord.Ptrs.size();

    std::vector< Data > lData;
    for( auto& k : lRecord.Ptrs ) lData.emplace_back( k->x - lRecord.CentreX , k->y - lRecord.CentreY , k->s );

    RoI lRoI( "ROI"+std::to_string( i ) , std::move( lData ) , lRecord.CentreX , lRecord.CentreY , lRecord.Size * lBinArea );
    aCallback( lRoI );
  }

}


void LocalizationFile::ExtractRoIs( const std::string& aImageJfile , const double& aScale , const std::function< void( RoI& ) >& aCallback ) const
{
  std::map< std::string , roi_polygon > lRoIs = OpenRoiZipfile( aImageJfile );

  typedef boost::geometry::model::point<double, 2, boost::geometry::cs::cartesian> geo_point;
  typedef boost::geometry::model::ring<geo_point> geo_polygon;

  for( const auto& j : lRoIs )
  {
    geo_polygon lPoly;
    for( const auto& point : j.second ) boost::geometry::append( lPoly , geo_point( boost::geometry::get<0>(point) * aScale , boost::geometry::get<1>(point) * aScale ) );
    boost::geometry::correct( lPoly ); // Add closing points, etc.

    geo_point lCentroid( 0 , 0 );
    boost::geometry::centroid( lPoly , lCentroid );
    double lCentreX( boost::geometry::get<0>(lCentroid) ), lCentreY( boost::geometry::get<1>(lCentroid) );

    std::vector< Data > lData;

    for( const auto& i : mData )
    {
      geo_point lPoint( i.x , i.y );    
      if( boost::geometry::within( lPoint , lPoly ) ) lData.emplace_back( i.x - lCentreX , i.y - lCentreY , i.s );
    }

    RoI lRoI( j.first , std::move( lData ) , lCentreX , lCentreY , boost::geometry::area( lPoly ) );
    aCallback( lRoI );
  }


}


void LocalizationFile::ExtractRoIsFromSegmentedImage(const std::string& aSegmentedImagefile, const double& aScale, const std::function< void(RoI&) >& aCallback) const
{  
       namespace bg = boost::gil;
       bg::rgb8_image_t img;
       bg::read_image(aSegmentedImagefile, img , bg::bmp_tag() );
       auto view = bg::view(img);

       std::size_t lwidth = view.width();
       std::size_t lheight = view.height();

       // find number of ROIs in the image
       int roi_number = 0;
       for (int y = 0; y < view.height(); ++y)
       {
           for (int x = 0; x < view.width(); ++x)
           {
               int label = get_color(view(x,y),bg::red_t());
               if (label > roi_number) roi_number = label;
           }
       }
       //
       std::vector<double>  Area(roi_number, 0), 
                            CX(roi_number, 0), 
                            CY(roi_number, 0);
       // calculate areas and centroids of these ROIs - in segmented image pixels
        for (int y = 0; y < view.height(); ++y)
        {
            for (int x = 0; x < view.width(); ++x)
            {
                int label = get_color(view(x, y), bg::red_t());
                if (0 != label)
                {
                    Area[label-1]++;
                    CX[label-1] += x;
                    CY[label-1] += y;
                }
            }
        }
        // normalize centroids accumulators
        for (int k = 0; k < roi_number; k++)
        {
            CX[k] = CX[k] / Area[k];
            CY[k] = CY[k] / Area[k];
        }        
        //
        //do actual job
        for (int k = 0; k < roi_number; k++)
        {
            std::vector< Data > lData;

            double  lCentreX = CX[k] * aScale,
                    lCentreY = CY[k] * aScale;

            for (const auto& i : mData)
            {
                int x = (int)(i.x / aScale), // coordinates at segmented image
                    y = (int)(i.y / aScale);
                //
                if (x < 0 || y < 0 || x >= (int)lwidth || y >= (int)lheight) continue;
                //
                if ((k+1)==(int)get_color(view(x, y), bg::red_t()))
                    lData.emplace_back(i.x - lCentreX, i.y - lCentreY, i.s);
            }

            RoI lRoI(std::to_string(k+1), std::move(lData), lCentreX, lCentreY, Area[k]*aScale*aScale);
            aCallback(lRoI);
        } 
}

// void LocalizationFile::ExtractRoIs( const std::string& aImageMap , const std::function< void( RoI& ) >& aCallback ) const
// {
//   namespace bg = boost::gil;
//   bg::rgb8_image_t img;
//   bg::read_image( aImageMap , img , bg::bmp_tag() );
//   auto view = bg::view(img);

//   // Calculate our scaling factor
//   std::pair< double, double> lXbound( std::make_pair( 9e99, -9e99 ) ), lYbound( std::make_pair( 9e99, -9e99 ) );

//   for( auto& k : mData ) {
//     if ( k.x < lXbound.first  ) lXbound.first  = k.x;
//     if ( k.x > lXbound.second ) lXbound.second = k.x;
//     if ( k.y < lYbound.first  ) lYbound.first  = k.y;
//     if ( k.y > lYbound.second ) lYbound.second = k.y;
//   }

//   std::size_t lWidth = view.width();
//   std::size_t lHeight = view.height();

//   auto lXScale( (lWidth - 1e-12) / ( lXbound.second - lXbound.first ) );
//   auto lYScale( (lHeight - 1e-12) / ( lYbound.second - lYbound.first ) );

//   std::map < int , tRecord > lRecords;

//   // Find the bounds and the datapoints
//   for( auto& k : mData ) {
//     std::size_t x = ( k.x - lXbound.first ) * lXScale;
//     std::size_t y = ( k.y - lYbound.first ) * lYScale;

//     // std::cout << lWidth << " " << lHeight << " " << x << " " << y << std::endl;

//     auto& pixel = view( x , y );
//     auto Id = bg::semantic_at_c<0>( pixel ) + ( 256 * bg::semantic_at_c<1>( pixel ) ) + ( 65536 * bg::semantic_at_c<2>( pixel ) );
//     if( Id == 0xFFFFFF ) continue;

//     auto lIt = lRecords.find( Id );
//     if( lIt == lRecords.end() ) lIt = lRecords.emplace( std::make_pair( Id , tRecord{ 0, std::make_pair( 9e99, -9e99 ), std::make_pair( 9e99, -9e99 ), std::vector< const Data* >() } ) ).first;
//     auto& lRecord = lIt->second;

//     lRecord.Ptrs.push_back( &k );
//     if ( k.x < lRecord.X.first  ) lRecord.X.first  = k.x;
//     if ( k.x > lRecord.X.second ) lRecord.X.second = k.x;
//     if ( k.y < lRecord.Y.first  ) lRecord.Y.first  = k.y;
//     if ( k.y > lRecord.Y.second ) lRecord.Y.second = k.y;
//   }

//   std::vector < tRecord* > lRecords2;
//   for( auto& i : lRecords ) lRecords2.push_back( &i.second );
//   std::sort( lRecords2.begin() , lRecords2.end() , []( const tRecord* a , const tRecord* b ){ return a->Ptrs.size() < b->Ptrs.size(); } );

//   for( auto& lIt : lRecords2 ) {
//     auto& lRecord = *lIt;
//     double lCentreX( ( lRecord.X.second + lRecord.X.first ) / 2.0 ), lCentreY( ( lRecord.Y.second + lRecord.Y.first ) / 2.0 );
//     double lWidthX( lRecord.X.second - lRecord.X.first ), lWidthY( lRecord.Y.second - lRecord.Y.first );

//     std::vector< Data > lData;

//     for( auto& k : lRecord.Ptrs ) {
//       double x = k->x - lCentreX;
//       double y = k->y - lCentreY;
//       lData.emplace_back( x, y, k->s );
//     }

//     RoI lRoI( std::move( lData ) );
//     lRoI.SetCentre( lCentreX, lCentreY );
//     lRoI.SetWidth( lWidthX, lWidthY );

//     aCallback( lRoI );
//   }

// }
// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
