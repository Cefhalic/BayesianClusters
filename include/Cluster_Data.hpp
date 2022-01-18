#pragma once

/* ===== C++ ===== */
#include <math.h>
#include <vector>
#include <functional>
#include <mutex>
#include <atomic>

#define PRECISION float

// class Data;

class Cluster
{
public:
  struct Parameter
  {
    Parameter();
    Parameter& operator+= ( const Parameter& aOther );
    double log_score() const;
    PRECISION A , Bx, By, C, logF;
  }; 

  Cluster();

  Cluster( Cluster&& aOther );
  Cluster& operator=( Cluster&& aOther );

  std::vector< Parameter > mParams;
  std::size_t mClusterSize , mLastClusterSize;
  PRECISION mClusterScore;
  Cluster* mParent;
  std::mutex mMutex;
  // std::vector< Cluster* > mMergeList;

  // void AddToMergeList( Cluster* aCluster );

  // Cluster& operator+= ( const Data& aData );
  Cluster& operator+= ( Cluster& aData );

  Cluster* GetParent();

  double log_score();
};



/* ===== Struct for storing data ===== */
class Data
{
public:
  Data( const PRECISION& aX , const PRECISION& aY , const PRECISION& aS );
  virtual ~Data();

  Data( const Data& ) = delete;
  Data& operator = (const Data& ) = delete;

  Data( Data&& ) = default;
  Data& operator = ( Data&& ) = default;

  inline bool operator< ( const Data& aOther ) const
  { 
    return r < aOther.r; 
  }

  inline PRECISION dR2( const Data& aOther ) const
  {
    PRECISION dX( x - aOther.x ), dY( y - aOther.y );
    return ( dX*dX ) + ( dY*dY );
  }

  inline PRECISION dR( const Data& aOther ) const
  {
    return sqrt( dR2( aOther ) );
  }

  inline PRECISION dPhi( const Data& aOther ) const
  {
   return fabs( phi - aOther.phi );
  }

  inline Cluster* cluster()
  {
    if( mCluster ) return mCluster = mCluster->GetParent();
    return NULL;
  }

  void PopulateNeighbours( std::vector<Data>::iterator aPlusIt , const std::vector<Data>::iterator& aPlusEnd , std::vector<Data>::reverse_iterator aMinusIt , const std::vector<Data>::reverse_iterator& aMinusEnd );
  void UpdateLocalization( const PRECISION& aR2 , const size_t& Nminus1  );

  void Clusterize( const PRECISION& a2R2 , const PRECISION& aT , std::vector< Cluster >& aClusters );
  
public:
  PRECISION x, y, r2 , r, phi;
  std::vector< PRECISION > mWeights;

  PRECISION mLocalizationSum , mLocalizationScore;

  std::vector< std::pair< PRECISION , Data* > > mNeighbours;
  std::vector< std::pair< PRECISION , Data* > >::iterator mNeighbourit;
  Cluster* mCluster;
  Cluster* mProtoCluster;
};



// void HandleMergeList( std::vector< Cluster >& aClusters );
bool CheckClusterization( std::vector<Data>& aData , std::vector< Cluster >& aClusters , const double& R , const double& T );
void Clusterize( std::vector<Data>& aData , std::vector< Cluster >& aClusters , const double& R , const double& T );
void ScanRT( std::vector<Data>& aData , const std::function< void( const std::vector< Cluster >& , const double& , const double& ) >& aCallback );
void PrepData( std::vector<Data>& aData );


