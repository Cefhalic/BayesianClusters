#pragma once

/* ===== C++ ===== */
#include <math.h>
#include <vector>
#include <functional>

/* ===== Local utilities ===== */
#include "Cluster_GlobalVars.hpp"


#define PRECISION float

class Data;
class DataWrapper;
class Cluster;


class Event
{
public:
  Event( const double& aPhysicalCentreX , const double& aPhysicalCentreY );

  std::vector<Data> mData;
  
  static GlobalVars mParameters;
  double mPhysicalCentreX , mPhysicalCentreY;


  inline double toPhysicalX( const double& aAlgorithmX ) const 
  {
    return mParameters.toPhysicalUnits( aAlgorithmX ) + mPhysicalCentreX;
  }

  inline double toAlgorithmX( const double& aPhysicalX ) const
  {
    return mParameters.toAlgorithmUnits( aPhysicalX - mPhysicalCentreX );
  }

  inline double toPhysicalY( const double& aAlgorithmY ) const 
  {
    return mParameters.toPhysicalUnits( aAlgorithmY ) + mPhysicalCentreY;
  }

  inline double toAlgorithmY( const double& aPhysicalY ) const
  {
    return mParameters.toAlgorithmUnits( aPhysicalY - mPhysicalCentreY );
  }

  void PrepData();
};




class EventWrapper
{
public:
  EventWrapper( Event& aEvent );
  
  std::vector<DataWrapper> mData;
  std::vector< Cluster > mClusters;

  std::size_t mClusteredCount , mBackgroundCount , mClusterCount;
  double mLogP;

  void CheckClusterization( const double& R , const double& T );
  // void Clusterize( const double& R , const double& T );
  void ScanRT( const std::function< void( const EventWrapper& , const double& , const double& ) >& aCallback , const uint8_t& aParallelization = 1 , const uint8_t& aOffset = 0 );
  void UpdateLogScore();
};



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

  std::vector< Parameter > mParams;
  std::size_t mClusterSize , mLastClusterSize;
  PRECISION mClusterScore;
  Cluster* mParent;

  Cluster& operator+= ( const Cluster& aOther );

  Cluster* GetParent();

  void UpdateLogScore();
};



/* ===== Struct for storing data ===== */
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
  // std::vector< PRECISION > mWeights;
  std::vector< std::pair< PRECISION , std::size_t > > mNeighbours;
  Cluster* mProtoCluster;  
};



class DataWrapper
{
public:
  DataWrapper( Data& aData );

  DataWrapper( const DataWrapper& ) = delete;
  DataWrapper& operator = (const DataWrapper& ) = delete;

  DataWrapper( DataWrapper&& ) = default;
  DataWrapper& operator = ( DataWrapper&& ) = default;

  // void UpdateLocalization( const PRECISION& aR2 , const size_t& Nminus1  );

  DataWrapper& GetNeighbour( EventWrapper& aEvent , const std::size_t& aIndex );

  void Clusterize( const PRECISION& a2R2 , const PRECISION& aT , EventWrapper& aEvent );
  void Clusterize( const PRECISION& a2R2 , const PRECISION& aT , EventWrapper& aEvent , Cluster* aCluster );
  
  inline Cluster* GetCluster()
  {
    if( ! mCluster ) return NULL;
    return mCluster = mCluster->GetParent();
  }

public:
  Data* mData;

  // PRECISION mLocalizationSum , mLocalizationScore;
  // std::vector< std::pair< PRECISION , std::size_t > >::iterator mNeighbourit;
  bool mExclude;
  Cluster* mCluster;
};







