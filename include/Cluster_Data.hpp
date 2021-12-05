#pragma once

/* ===== C++ ===== */
#include <vector>
#include <list>
#include <array>



// // Burmann approximation of the Gaussian cumulative-distribution function
// double Phi( const double& x )
// {
//   auto y = exp( -1.0*x*x );

//   constexpr double pi = atan(1)*4;
//   constexpr double rt_pi = sqrt( pi );
//   constexpr double inv_rt_pi = 1.0 / rt_pi;

//   double lRet = ( rt_pi/2.0 ) + ( y * 31.0/200.0 ) - ( y*y * 341.0/8000.0 );
//   lRet       *= inv_rt_pi * ( (x > 0) - (x < 0) ) * sqrt( 1-y );

//   return 0.5 + lRet;
// }



/* ===== Struct for storing data ===== */
class Data
{
public:
  struct ClusterParameter
  {
    ClusterParameter( const double& aW );
    
    void Reset( const double& aX , const double& aY);

    ClusterParameter& operator+= ( const ClusterParameter& aOther );

    const double w , logw;
    double n_tilde , sum_logw , nu_bar_x , nu_bar_y;
  };  

public:
  Data( const double& aX , const double& aY , const double& aS );

  bool operator< ( const Data& aOther ) const;

  double dR2( const Data& aOther ) const;
  double dR( const Data& aOther ) const;

  void PopulateNeighbours( std::vector<Data>::iterator aPlusIt , const std::vector<Data>::iterator& aPlusEnd , std::vector<Data>::reverse_iterator aMinusIt , const std::vector<Data>::reverse_iterator& aMinusEnd );
  void UpdateLocalization( const double& aR2 , const size_t& Nminus1  );
  void ResetClusters();

  void Clusterize( const double& a2R2 , const double& aT );
  void ClusterInto( Data* aParent );

  double S2( const std::size_t& index , const double& nu_bar_x , const double& nu_bar_y ) const;

  void UpdateClusterScore();


public:
  // std::unique_ptr< std::mutex > mutex;
  double x, y, s , r, phi;
  double eX , eY;
  double localizationsum , localizationscore;

  std::array< std::vector< std::pair< double , Data* > > , 2 > neighbours;
  std::vector< std::pair< double , Data* > >::iterator neighbourit;

  Data* parent;
  std::list< Data* > children;
  std::size_t ClusterSize;
  double ClusterScore;

  std::vector< ClusterParameter > ClusterParams;

};



