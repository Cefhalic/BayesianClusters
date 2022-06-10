

/* ===== C++ ===== */
#include <stdlib.h>
#include <iostream>
#include <iomanip>

/* ===== For Root ===== */
#include "Math/ProbFunc.h" 
#include "Math/Interpolator.h" 
#include "Math/SpecFunc.h" 

/* ===== Cluster sources ===== */
#include "Cluster_GlobalVars.hpp"
#include "Cluster_Data.hpp"

/* ===== Local utilities ===== */
#include "ListComprehension.hpp"
#include "ProgressBar.hpp"
#include "Vectorize.hpp"



// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
GlobalVars Event::mParameters;

Event::Event()
{
  const std::string& lFilename = Event::mParameters.inputFile();
  if( lFilename.size() == 0 ) throw std::runtime_error( "No input file specified" ); 

  LoadCSV( lFilename );
  Preprocess();  
}

void Event::Preprocess()
{
  // Populate mNeighbour lists  
  ProgressBar2 lProgressBar( "Populating neighbourhood" , mData.size() );
  [&]( const std::size_t& i ){ mData.at( i ).Preprocess( mData , i ); } || range( mData.size() );  // Interleave threading since processing time increases with radius from origin
}

void Event::ScanRT( const std::function< void( const EventProxy& , const double& , const double& ) >& aCallback )
{
  const auto N = Concurrency + 1;
  std::vector< EventProxy > lEventProxys;
  lEventProxys.reserve( N );
  for( int i(0) ; i!=N ; ++i ) lEventProxys.emplace_back( *this );
  ProgressBar2 lProgressBar( "Scan over RT"  , 0 );
  [&]( const std::size_t& i ){ lEventProxys.at(i).ScanRT( aCallback , N , i ); } || range( N );
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
    double x = Event::mParameters.toAlgorithmX( strtod( ch , &lPtr ) * nanometer );
    ReadUntil( ',' ); //"y [nm]"
    double y = Event::mParameters.toAlgorithmY( strtod( ch , &lPtr ) * nanometer );      
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

  std::sort( aData.begin() , aData.end() );
}

void Event::LoadCSV( const std::string& aFilename )
{
  auto f = fopen( aFilename.c_str() , "rb");
  if ( f == NULL ) throw std::runtime_error( "File is not available" );
  fseek(f, 0, SEEK_END); // seek to end of file
  auto lSize = ftell(f); // get current file pointer
  fclose(f);

  int lChunkSize = ceil( double(lSize) / (Concurrency+1) );
  std::vector< std::vector< Data > > lData( Concurrency+1 );

  ProgressBar2 lProgressBar( "Reading File" , lSize );
  for( std::size_t i(0); i!=Concurrency ; ++i ) ThreadPool.at(i)->submit( [ & , i ](){ __LoadCSV__( aFilename , *this , lData[i] , i*lChunkSize , lChunkSize ); } );
  WrappedThread::run_and_wait( [ & ](){ __LoadCSV__( aFilename , *this , lData[Concurrency] , Concurrency*lChunkSize , lChunkSize ); } );

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

  std::cout << "Read " << mData.size() << " points" << std::endl;
}

void Event::WriteCSV( const std::string& aFilename )
{
  auto f = fopen( aFilename.c_str() , "w");
  if ( f == NULL ) throw std::runtime_error( "File is not available" );

  fprintf( f , "id,frame,x [nm],y [nm],sigma [nm],intensity [photon],offset [photon],bkgstd [photon],chi2,uncertainty_xy [nm]\n" );

  ProgressBar lProgressBar( "Writing File" , mData.size() );
  for( auto& i : mData ){
    fprintf( f , ",,%f,%f,,,,,,%f\n" , Event::mParameters.toPhysicalX(i.x)/nanometer , Event::mParameters.toPhysicalY(i.y)/nanometer , Event::mParameters.toPhysicalUnits(i.s)/nanometer );
    lProgressBar++;
  }

  fclose(f);
}
// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
EventProxy::EventProxy( Event& aEvent ) :
  mBackgroundCount( 0 )
{
  mClusters.reserve( aEvent.mData.size() );  // Reserve as much space for clusters as there are data points - prevent pointers being invalidated!
  mData.reserve( aEvent.mData.size() );
  for( auto& i : aEvent.mData ) mData.emplace_back( i );
}

void EventProxy::CheckClusterization( const double& R , const double& T )
{
  const auto lRlimit = 4.0 * R * R;

  uint32_t lClusterCount( 0 );

  for( auto& i : mClusters )
  {
    if( i.mClusterSize ) ++lClusterCount;
  }

  if( lClusterCount != mClusterCount )
  {
    std::cout << "\nR = " << R << ", T = " << T << " | #Clusters = " << mClusterCount << " | Expected #Clusters = " << lClusterCount << std::endl;
    throw std::runtime_error( "Check failed" );
  }  

  uint32_t lBackgroundCount( 0 );
  uint32_t lPointsInClusters( 0 );

  uint32_t lExpected( 0 );
  uint32_t lNotClustered( 0 );
  uint32_t lNeighbourNotClustered( 0 );
  uint32_t lWrongNeighbour( 0 );

  for( auto& i : mData )
  {
    if( i.mExclude ){ 
      lBackgroundCount++;
      continue;
    }
    
    lExpected++;
    if( ! i.GetCluster() ){ lNotClustered++ ; continue; }
    for( auto& j : i.mData->mNeighbours )
    {
      if( j.first > lRlimit ) break;
      auto& lNeighbour( i.GetNeighbour( *this , j.second ) );  

      if( lNeighbour.mExclude ) continue;

      if( ! i.mCluster ){ lNeighbourNotClustered++; continue; }
      if ( lNeighbour.mCluster->GetParent() != i.mCluster )
      { 
        lWrongNeighbour++;
        continue; 
      }
    }    
  }

  for( auto& i : mClusters ) lPointsInClusters += i.mClusterSize;

  if( lBackgroundCount != mBackgroundCount )
  {
    std::cout << "\nR = " << R << ", T = " << T << " | Background = " << mBackgroundCount  << " | Expected Background = " << lBackgroundCount << std::endl;
    throw std::runtime_error( "Check failed" ); 
  }  

  if( lPointsInClusters + lBackgroundCount != mData.size() )
  {
    std::cout << "\nR = " << R << ", T = " << T << " | Points In Clusters = " << lPointsInClusters  << " | Background = " << lBackgroundCount << " | Total = " << mData.size() << std::endl;
    throw std::runtime_error( "Check failed" ); 
  }  

  if( lNotClustered or lNeighbourNotClustered or lWrongNeighbour )
  {
    std::cout << "\nR = " << R << ", T = " << T << " | Not Clustered = " << lNotClustered << "/" << lExpected << " | Neighbour Not Clustered = " << lNeighbourNotClustered << " | Wrong Neighbour = " << lWrongNeighbour << std::endl;
    throw std::runtime_error( "Check failed" );
  }
}

__attribute__((flatten))
void EventProxy::ScanRT( const std::function< void( const EventProxy& , const double& , const double& ) >& aCallback , const uint8_t& aParallelization , const uint8_t& aOffset )
{
  double dR( aParallelization * Event::mParameters.dR() );
  double R( Event::mParameters.minScanR() + ( aOffset * Event::mParameters.dR() ) ) , R2( 0 ) , twoR2( 0 ) , T( 0 );

  for( uint32_t i( aOffset ) ; i<Event::mParameters.Rbins() ; i+=aParallelization , R+=dR )
  {
    R2 = R * R;
    twoR2 = 4.0 * R2;
    T = Event::mParameters.maxScanT();

    mClusters.clear();
    for( auto& k : mData ) k.mCluster = NULL;

    for( uint32_t j(0) ; j!=Event::mParameters.Tbins() ; ++j , T-=Event::mParameters.dT() )
    {
      for( auto& k : mData ) k.mExclude = ( k.mData->mLocalizationScores[ i ] < T ) ;
      for( auto& k : mData ) k.Clusterize( twoR2 , T , *this );
      UpdateLogScore();
      if( Event::mParameters.validate() ) CheckClusterization( R , T ) ;
      aCallback( *this , R , T );
    }
  }

  mClusters.clear();
  for( auto& k : mData ) k.mCluster = NULL; // Clear cluster pointers which will be invalidated when we leave the function
}

void EventProxy::UpdateLogScore()
{
  mClusterCount = mClusteredCount = 0;
  mLogP = 0.0;

  for( auto& i: mClusters )
  {
    if( i.mClusterSize == 0 ) continue;
    i.UpdateLogScore();
    mClusterCount += 1;
    mClusteredCount += i.mClusterSize;
    mLogP += ROOT::Math::lgamma( i.mClusterSize );
  }

  mBackgroundCount = mData.size() - mClusteredCount;
  mLogP += ( mBackgroundCount * Event::mParameters.logPb() ) 
         + ( mClusteredCount * Event::mParameters.logPbDagger() )
         + ( Event::mParameters.logAlpha() * mClusterCount )
         + Event::mParameters.logGammaAlpha()
         - ROOT::Math::lgamma( Event::mParameters.alpha() + mClusteredCount );  
}
// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Cluster::Parameter::Parameter() : 
A(0.0) , Bx(0.0) , By(0.0) , C(0.0) , logF(0.0)
{}
    
Cluster::Parameter& Cluster::Parameter::operator+= ( const Cluster::Parameter& aOther )
{
  A += aOther.A;
  Bx += aOther.Bx;
  By += aOther.By;
  C += aOther.C;
  logF += aOther.logF;
  return *this;
}

inline double CDF( const double& aArg )
{
  // Above or below ~8 are indistinguishable from 0 and 1 respectively
  if( aArg > 8.0 ) return 1.0;
  if( aArg < 8.0 ) return 0.0;
  return  ROOT::Math::normal_cdf( aArg );
}

__attribute__((flatten))
double Cluster::Parameter::log_score() const
{
  auto sqrt_A( sqrt( A ) ) , inv_A( 1.0 / A );
  auto Dx( Bx * inv_A ) , Dy( By * inv_A );
  auto E( C - ( Bx * Dx ) - ( By * Dy ) );

  double log_sum = logF - double( log( A ) ) + ( 0.5 * E );

  // We place explicit bounds checks to prevent calls to expensive functions
  auto Gx = CDF( sqrt_A * (1.0-Dx) ) - CDF( sqrt_A * (-1.0-Dx) );
  if( Gx != 1.0 ) log_sum += log( Gx );
  auto Gy = CDF( sqrt_A * (1.0-Dy) ) - CDF( sqrt_A * (-1.0-Dy) );
  if( Gy != 1.0 ) log_sum += log( Gy );

  return log_sum;
}
// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Cluster::Cluster(): mParams( Event::mParameters.sigmacount() ),
mClusterSize( 0 ) , mLastClusterSize( 0 ) , mClusterScore( 0.0 ) , 
mParent( NULL )
{}

Cluster::Cluster( const Data& aData ): mParams( Event::mParameters.sigmacount() ),
mClusterSize( 1 ) , mLastClusterSize( 0 ) , mClusterScore( 0.0 ) , 
mParent( NULL )
{ 
  const auto s2 = aData.s * aData.s;
  auto lIt( mParams.begin() ) ;
  auto lSig2It( Event::mParameters.sigmabins2().begin() );

  for( ; lIt != mParams.end() ; ++lIt , ++lSig2It )
  {
    double w = 1.0 / ( s2 + *lSig2It );
    lIt->A = w;
    lIt->Bx = (w * aData.x);
    lIt->By = (w * aData.y);
    lIt->C = (w * aData.r2);
    lIt->logF = PRECISION( log( w ) );
  }
}

void Cluster::UpdateLogScore()
{
  static constexpr double pi = atan(1)*4;
  static constexpr double log2pi = log( 2*pi );

  if( mClusterSize <= mLastClusterSize ) return; // We were not bigger than the previous size when we were evaluated - score is still valid
  mLastClusterSize = mClusterSize;

  thread_local static std::vector< double > MuIntegral( Event::mParameters.sigmacount() , 1.0 );

  // double constant( mParams[0].log_score() + Event::mParameters.log_probability_sigma( 0 ) );
  // for( std::size_t i(1) ; i!=Event::mParameters.sigmacount() ; ++i ) MuIntegral[i] = exp( mParams[i].log_score() + Event::mParameters.log_probability_sigma( i ) - constant );
  for( std::size_t i(0) ; i!=Event::mParameters.sigmacount() ; ++i ) MuIntegral[i] = exp( mParams[i].log_score() + Event::mParameters.log_probability_sigma( i ) );

  thread_local static ROOT::Math::Interpolator lInt( Event::mParameters.sigmacount() , ROOT::Math::Interpolation::kLINEAR );
  lInt.SetData( Event::mParameters.sigmabins() , MuIntegral );

  static const double Lower( Event::mParameters.sigmabins(0) ) , Upper( Event::mParameters.sigmabins(Event::mParameters.sigmacount()-1) );
  // mClusterScore = double( log( lInt.Integ( Lower , Upper ) ) ) + constant - double( log( 4.0 ) ) + (log2pi * (1.0-mClusterSize));  
  mClusterScore = double( log( lInt.Integ( Lower , Upper ) ) ) - double( log( 4.0 ) ) + (log2pi * (1.0-mClusterSize));  
}

Cluster& Cluster::operator+= ( const Cluster& aOther )
{
  // if( &aOther == this ) throw std::runtime_error( "Error #1" );
  // if( aOther.mClusterSize == 0 ) throw std::runtime_error( "Error #2" );
  // if( aOther.mParent ) throw std::runtime_error( "Error #3" );

  auto lIt( mParams.begin() );
  auto lIt2( aOther.mParams.begin() );

  for( ; lIt != mParams.end() ; ++lIt , ++lIt2 ) *lIt += *lIt2;
  mClusterSize += aOther.mClusterSize;
  return *this;
}

Cluster* Cluster::GetParent()
{
  if( mParent ) return mParent = mParent->GetParent();
  return this;
}
// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Data::Data( const PRECISION& aX , const PRECISION& aY , const PRECISION& aS ) : 
x(aX) , y(aY) , s(aS) , r2( (aX*aX) + (aY*aY) ), r( sqrt( r2 ) ), phi( atan2( aY , aX ) ),
mProtoCluster( NULL )
{}

Data::~Data()
{
  if( mProtoCluster ) delete mProtoCluster;
  mProtoCluster = NULL;
}

// Although the neighbourhood calculation is reciprocal (if I am your neighbour then you are mine) and we can, in fact, use that to halve the number of calculations,
// doing so requires arbitration between threads or a single-threaded reciprocation step, both of which take longer than brute-forcing it
__attribute__((flatten))
void Data::Preprocess( std::vector<Data>& aData , const std::size_t& aIndex )
{
  static constexpr double pi = atan(1)*4;  
  auto dphi = asin( Event::mParameters.max2R() / r ); // Event::mParameters.max2R() / ( r - aParameters.max2R() );
  auto dphi2 = (2*pi) - dphi;

  std::size_t i( aIndex + 1 );
  std::vector<Data>::const_iterator aPlusIt( aData.begin() + aIndex + 1 );
  const std::vector<Data>::const_iterator aPlusEnd( aData.end() );

  // Iterate over other hits and populate the mNeighbour list
  for( ; aPlusIt != aPlusEnd ; ++aPlusIt , ++i )
  {
    if( ( aPlusIt->r - r ) > Event::mParameters.max2R() ) break; // aPlusIt is always further out than curent 
    auto lPhi = dPhi( *aPlusIt );
    if( lPhi > dphi and lPhi < dphi2 ) continue;
    PRECISION ldR2 = dR2( *aPlusIt );
    if( ldR2 < Event::mParameters.max2R2() ) mNeighbours.push_back( std::make_pair( ldR2 , i ) );
  }

  i = aIndex - 1;
  std::vector<Data>::const_reverse_iterator aMinusIt( aData.rbegin() + aData.size() - aIndex );
  const std::vector<Data>::const_reverse_iterator aMinusEnd( aData.rend() );

  for( ; aMinusIt != aMinusEnd ; ++aMinusIt , --i )
  {
    if( ( r - aMinusIt->r ) > Event::mParameters.max2R() ) break; // curent is always further out than aMinusIn
    auto lPhi = dPhi( *aMinusIt );
    if( lPhi > dphi and lPhi < dphi2 ) continue;
    PRECISION ldR2 = dR2( *aMinusIt );    
    if( ldR2 < Event::mParameters.max2R2() ) mNeighbours.push_back( std::make_pair( ldR2 , i ) );
  }

  std::sort( mNeighbours.begin() , mNeighbours.end() );

  // -------------------------------------------------------------------------------------

  mProtoCluster = new Cluster( *this );

  // -------------------------------------------------------------------------------------

  const double lLocalizationConstant( 4.0 / ( pi * ( aData.size() - 1 ) ) ); 
  const PRECISION eX( 1 - fabs( x ) ) , eY( 1 - fabs( y ) );

  auto lNeighbourit( mNeighbours.begin() );
  PRECISION lLocalizationSum( 0 ) , lLastLocalizationSum( 0 ) , lLocalizationScore( 0 ) , lDist( 0 ) , lWeight( 0 );
  mLocalizationScores.reserve( Event::mParameters.Rbins() );

  double R( Event::mParameters.minScanR() ) , R2( 0 );
  for( uint32_t i(0) ; i!=Event::mParameters.Rbins() ; ++i , R+=Event::mParameters.dR() )
  {
    R2 = R * R;

    for(  ; lNeighbourit != mNeighbours.end() ; ++lNeighbourit )
    { 
      if( lNeighbourit->first > R2 ) break;
      lDist = sqrt( lNeighbourit->first );

      // Noticeably faster polynomial approximation of the edge-correction
      lWeight = 1.0;
      if( eX < lDist )  lWeight *= ( 1 + pow( acos( eX/lDist ) * (2/pi) , 4 ) );
      if( eY < lDist )  lWeight *= ( 1 + pow( acos( eY/lDist ) * (2/pi) , 4 ) );
      lLocalizationSum += lWeight;
    }

    if( lLastLocalizationSum != lLocalizationSum )
    {
      lLocalizationScore = sqrt( lLocalizationConstant * lLocalizationSum );
      lLastLocalizationSum = lLocalizationSum;
    }

    mLocalizationScores.push_back( lLocalizationScore );
  }
}
// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
DataProxy::DataProxy( Data& aData ) :
mData( &aData ),
mCluster( NULL )
{}

void DataProxy::Clusterize( const PRECISION& a2R2 , const PRECISION& aT , EventProxy& aEvent ) // We are at the top-level
{
  if( mCluster || mExclude ) return;

  aEvent.mClusters.emplace_back();
  Clusterize( a2R2 , aT , aEvent , &aEvent.mClusters.back() );
}

void DataProxy::Clusterize( const PRECISION& a2R2 , const PRECISION& aT , EventProxy& aEvent , Cluster* aCluster )
{
  if( mCluster )
  {
    if( GetCluster() == aCluster ) return;
    *aCluster += *mCluster;
    mCluster->mParent = aCluster;
    mCluster->mClusterSize = 0;
    mCluster = aCluster;
  }
  else
  {
    if( mExclude ) return;

    *aCluster += *(mData->mProtoCluster);
    mCluster = aCluster;

    for( auto& i : mData->mNeighbours )
    {
      if( i.first > a2R2 ) break;
      GetNeighbour( aEvent , i.second ).Clusterize( a2R2 , aT , aEvent , aCluster );
    }  
  }
}
// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
