#pragma once

/* ===== C++ ===== */
#include <math.h>
#include <vector>
#include <functional>

/* ===== Local utilities ===== */
#include "Cluster_GlobalVars.hpp"


#define PRECISION float

#define PARALLELIZATION 8


class Data;
class DataProxy;
class Cluster;


class Event
{
public:
  Event( const double& aPhysicalCentreX , const double& aPhysicalCentreY );

  std::vector< Data > mData;

  struct Instance
  {
    std::size_t mIndex;
    std::vector< Cluster > mClusters;
    std::size_t mClusteredCount , mBackgroundCount , mClusterCount;
    double mLogP;

    void UpdateLogScore( const std::vector< Data >& aData );   
    void CheckClusterization( std::vector< Data >& aData , const double& R , const double& T );
  };

  std::vector< Instance > mInstances;

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


  // void Clusterize( const double& R , const double& T );
  void ScanRT( const std::function< void( const Event::Instance& , const double& , const double& ) >& aCallback );
  void Preprocess();

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

  void Preprocess( Event& aEvent );

  void Clusterize( const PRECISION& a2R2 , Event::Instance& aInstance );
  void Clusterize( const PRECISION& a2R2 , Event::Instance& aInstance , Cluster* aCluster );
  
  inline Cluster* GetCluster( const std::size_t& aIndex )
  {
    Cluster*& lCluster = mCluster[ aIndex ];
    if( ! lCluster ) return NULL;
    return lCluster = lCluster->GetParent();
  }

  // inline Cluster* GetCluster( const std::size_t& aIndex ) const
  // {
  //   Cluster* const& lCluster = mCluster[ aIndex ];
  //   if( ! lCluster ) return NULL;
  //   return lCluster->GetParent();
  // }

public:
  PRECISION x, y, s , r2 , r, phi;
  std::vector< PRECISION > mLocalizationScores;
  std::vector< std::pair< PRECISION , Data* > > mNeighbours;

  std::vector< bool > mExclude;
  std::vector< Cluster* > mCluster;
  Cluster* mProtoCluster;
};
