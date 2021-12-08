#pragma once

/* ===== C++ ===== */
#include <vector>
#include <array>

#define PRECISION float


/* ===== Struct for storing data ===== */
class Data
{
public:
  struct ClusterParameter
  {
    ClusterParameter( const PRECISION& aW );

    void Reset( const PRECISION& aX , const PRECISION& aY);

    ClusterParameter& operator+= ( const ClusterParameter& aOther );

    double log_score( const std::size_t& n ) const;

    const PRECISION w , logw;
    PRECISION A , Bx, Cx, By, Cy, sum_logw;
  };  

public:
  Data( const PRECISION& aX , const PRECISION& aY , const PRECISION& aS );

  Data( const Data& ) = delete;
  Data& operator = (const Data&) = delete;

  Data( Data&& ) = default;
  Data& operator = ( Data&& ) = default;

  bool operator< ( const Data& aOther ) const;

  PRECISION dR2( const Data& aOther ) const;
  PRECISION dR( const Data& aOther ) const;

  void PopulateNeighbours( std::vector<Data>::iterator aPlusIt , const std::vector<Data>::iterator& aPlusEnd , std::vector<Data>::reverse_iterator aMinusIt , const std::vector<Data>::reverse_iterator& aMinusEnd );
  void UpdateLocalization( const PRECISION& aR2 , const size_t& Nminus1  );
  void ResetClusters();

  void Clusterize( const PRECISION& a2R2 , const PRECISION& aT , const Data* aLower , const Data* aUpper );
  void Clusterize2( const PRECISION& a2R2 , const PRECISION& aT );

  Data* GetParent();

  void ClusterInto( Data* aParent );

  void UpdateClusterScore();


public:
  PRECISION x, y, r, phi;
  PRECISION eX , eY;
  PRECISION localizationsum , localizationscore;

  std::vector< std::pair< PRECISION , Data* > > neighbours;
  std::vector< std::pair< PRECISION , Data* > >::iterator neighbourit;

  Data* parent;
  std::size_t ClusterSize , LastClusterSize;
  PRECISION ClusterScore;

  std::vector< ClusterParameter > ClusterParams;

};



