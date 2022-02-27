#pragma once

/* ===== C++ ===== */
#include <math.h>
#include <vector>
#include <functional>

/* ===== Local utilities ===== */
#include "Cluster_GlobalVars.hpp"


#define PRECISION float

class Event;
class EventProxy;
class Cluster;
class Data;
class DataProxy;


// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
class Event
{
public:
  Event();  
  Event( const std::string& aFilename );

  Event( const Event& ) = delete;
  Event& operator = (const Event& ) = delete;

  Event( Event&& ) = default;
  Event& operator = ( Event&& ) = default;

  void Preprocess();
  void ScanRT( const std::function< void( const EventProxy& , const double& , const double& ) >& aCallback );
  void LoadCSV( const std::string& aFilename );
  void WriteCSV( const std::string& aFilename );

public:
  std::vector<Data> mData;
  static GlobalVars mParameters;
};
// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
class EventProxy
{
public:
  EventProxy( Event& aEvent );

  EventProxy( const EventProxy& ) = delete;
  EventProxy& operator = (const EventProxy& ) = delete;

  EventProxy( EventProxy&& ) = default;
  EventProxy& operator = ( EventProxy&& ) = default;

  void CheckClusterization( const double& R , const double& T );
  void ScanRT( const std::function< void( const EventProxy& , const double& , const double& ) >& aCallback , const uint8_t& aParallelization = 1 , const uint8_t& aOffset = 0 );
  void UpdateLogScore();

public:
  std::vector<DataProxy> mData;
  std::vector< Cluster > mClusters;

  std::size_t mClusteredCount , mBackgroundCount , mClusterCount;
  double mLogP;
};
// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
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
  Cluster( const Data& aData );

  Cluster& operator+= ( const Cluster& aOther );

  Cluster* GetParent();

  void UpdateLogScore();

public:
  std::vector< Parameter > mParams;
  std::size_t mClusterSize , mLastClusterSize;
  PRECISION mClusterScore;
  Cluster* mParent;
};


// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
class Data
{
public:
  Data( const PRECISION& aX , const PRECISION& aY , const PRECISION& aS );

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

  void Preprocess( std::vector<Data>& aData , const std::size_t& aIndex );

public:
  PRECISION x, y, s , r2 , r, phi;
  std::vector< PRECISION > mLocalizationScores;
  std::vector< std::pair< PRECISION , std::size_t > > mNeighbours;
  Cluster* mProtoCluster;  
};
// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
class DataProxy
{
public:
  DataProxy( Data& aData );

  DataProxy( const DataProxy& ) = delete;
  DataProxy& operator = (const DataProxy& ) = delete;

  DataProxy( DataProxy&& ) = default;
  DataProxy& operator = ( DataProxy&& ) = default;

  inline DataProxy& GetNeighbour( EventProxy& aEvent , const std::size_t& aIndex )
  {
    return aEvent.mData.at( aIndex );
  }

  void Clusterize( const PRECISION& a2R2 , const PRECISION& aT , EventProxy& aEvent );
  void Clusterize( const PRECISION& a2R2 , const PRECISION& aT , EventProxy& aEvent , Cluster* aCluster );
  
  inline Cluster* GetCluster()
  {
    if( ! mCluster ) return NULL;
    return mCluster = mCluster->GetParent();
  }

public:
  Data* mData;
  Cluster* mCluster;
  bool mExclude;
};
// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
