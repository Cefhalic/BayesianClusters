#pragma once

/* ===== C++ ===== */
#include <math.h>
#include <vector>
#include <functional>

#define PRECISION float

/* ===== Struct for storing data ===== */
class Data
{
public:
  struct ClusterParameter
  {
    ClusterParameter();

    ClusterParameter& operator+= ( const ClusterParameter& aOther );

    double log_score() const;

    PRECISION A , Bx, By, C, logF;
  };  


  struct Cluster
  {
    Cluster();
    std::vector< ClusterParameter > Params;
    std::size_t ClusterSize , LastClusterSize;
    PRECISION ClusterScore;
    Cluster* mParent;

    void operator+= ( const Data& aData );
    void operator+= ( Cluster& aData );

    Cluster* GetParent();

    double log_score();
  };



public:
  Data( const PRECISION& aX , const PRECISION& aY , const PRECISION& aS );

  Data( const Data& ) = delete;
  Data& operator = (const Data&) = delete;

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

  void PopulateNeighbours( std::vector<Data>::iterator aPlusIt , const std::vector<Data>::iterator& aPlusEnd , std::vector<Data>::reverse_iterator aMinusIt , const std::vector<Data>::reverse_iterator& aMinusEnd );
  void UpdateLocalization( const PRECISION& aR2 , const size_t& Nminus1  );

  void Clusterize( const PRECISION& a2R2 , const PRECISION& aT );
  void ClusterInto( const PRECISION& aT , Cluster* aCluster );
  
public:
  PRECISION x, y, r2 , r, phi;
  std::vector< PRECISION > w_i;

  PRECISION localizationsum , localizationscore;

  std::vector< std::pair< PRECISION , Data* > > neighbours;
  std::vector< std::pair< PRECISION , Data* > >::iterator neighbourit;
  Cluster* mCluster;

  static std::vector< Cluster > Clusters;
};




void Clusterize( std::vector<Data>& aData , const double& twoR2 , const double& T );
void Cluster( std::vector<Data>& aData , const double& R , const double& T );
void ScanRT( std::vector<Data>& aData , const std::function< void( const double& , const double& ) >& aCallback );
void PrepData( std::vector<Data>& aData );