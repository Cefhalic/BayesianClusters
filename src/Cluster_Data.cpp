#pragma once

/* ===== C++ ===== */
#include <vector>
#include <list>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <mutex>
#include <memory>

/* ===== For Root ===== */
#include "Math/ProbFunc.h" 

/* ===== Cluster sources ===== */
#include "Cluster_GlobalVars.hpp"

/* ===== Local utilities ===== */
#include "ListComprehension.hpp"



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
    ClusterParameter( const double& aW ) : w(aW) , logw( log(aW) )
    {}
    
    inline void Reset( const double& aX , const double& aY)
    {
      n_tilde = w;
      sum_logw = logw;
      nu_bar_x = w * aX;
      nu_bar_y = w * aY;
    }

    inline ClusterParameter& operator+= ( const ClusterParameter& aOther )
    {
      n_tilde += aOther.n_tilde;
      sum_logw += aOther.sum_logw;
      nu_bar_x += aOther.nu_bar_x;
      nu_bar_y += aOther.nu_bar_y;
      return *this;
    }

    const double w , logw;
    double n_tilde , sum_logw , nu_bar_x , nu_bar_y;
  };  

public:
  Data( const double& aX , const double& aY , const double& aS ) : 
  // mutex( new std::mutex() ),
  x(aX) , y(aY) , s(aS) , r( sqrt( (aX*aX) + (aY*aY) ) ), phi( atan2( aY , aX ) ),
  eX( 1 - fabs( aX ) ) , eY( 1 - fabs( aY ) ) , 
  localizationsum( 0.0 ) , localizationscore( 0.0 ),
  neighbourit( neighbours[0].end() ),
  parent( NULL ), children( { this } ),
  ClusterSize( 0 ) , ClusterScore( 0.0 )
  {
    for( auto& sig2 : Parameters.sigmabins2() ) ClusterParams.emplace_back( 1 / ( (aS*aS) + sig2 ) );
  }

  bool operator< ( const Data& aOther ) const { return r < aOther.r; }

  inline double dR2( const Data& aOther ) const
  {
    double dX( x - aOther.x ), dY( y - aOther.y );
    return ( dX*dX ) + ( dY*dY );
  }

  inline double dR( const Data& aOther ) const
  {
    return sqrt( dR2( aOther ) );
  }


  void PopulateNeighbours( std::vector<Data>::iterator aPlusIt , const std::vector<Data>::iterator& aPlusEnd , std::vector<Data>::reverse_iterator aMinusIt , const std::vector<Data>::reverse_iterator& aMinusEnd )
  {

    // Iterate over other hits and populate the neighbour list
    for( ; aPlusIt != aPlusEnd ; aPlusIt++ )
    {
      if( ( aPlusIt->r - r ) > Parameters.maxR() ) break; // aPlusIt is always further out than curent 
      double ldR2 = dR2( *aPlusIt );
      if( ldR2 < Parameters.maxR2() ) neighbours[0].push_back( std::make_pair( ldR2 , &*aPlusIt ) );
    }

    for( ; aMinusIt != aMinusEnd ; aMinusIt++ )
    {
      if( ( r - aMinusIt->r ) > Parameters.maxR() ) break; // curent is always further out than aMinusIn
      double ldR2 = dR2( *aMinusIt );    
      if( ldR2 < Parameters.maxR2() ) neighbours[0].push_back( std::make_pair( ldR2 , &*aMinusIt ) );
    }

    std::sort( neighbours[0].begin() , neighbours[0].end() );
    neighbourit = neighbours[0].begin();

    // Iterate over other hits and populate the neighbour2 list
    for( ; aPlusIt != aPlusEnd ; aPlusIt++ )
    {
      if( ( aPlusIt->r - r ) > Parameters.max2R() ) break; // aPlusIt is always further out than curent  
      double ldR2 = dR2( *aPlusIt );
      if( ldR2 < Parameters.max2R2() ) neighbours[1].push_back( std::make_pair( ldR2 , &*aPlusIt ) );
    }

    for( ; aMinusIt != aMinusEnd ; aMinusIt++ )
    {
      if( ( r - aMinusIt->r ) > Parameters.max2R() ) break; // curent is always further out than aMinusIn
      double ldR2 = dR2( *aMinusIt );    
      if( ldR2 < Parameters.max2R2() ) neighbours[1].push_back( std::make_pair( ldR2 , &*aMinusIt ) );
    }

    std::sort( neighbours[1].begin() , neighbours[1].end() );
  }


  inline void UpdateLocalization( const double& aR2 , const size_t& Nminus1  )
  {
    constexpr double pi = atan(1)*4;

    const double last_localizationsum( localizationsum );

    for( ; neighbourit != neighbours[0].end() ; ++neighbourit )
    { 
      if( neighbourit->first > aR2 ) break;
      double lDist = sqrt( neighbourit->first );

      // Approximation of the edge-correction
      double Weight( 1.0 );
      if( eX < lDist )  Weight *= ( 1 + pow( acos( eX/lDist ) * (2/pi) , 4) );
      if( eY < lDist )  Weight *= ( 1 + pow( acos( eY/lDist ) * (2/pi) , 4) );

      localizationsum += Weight;
    }

    if( last_localizationsum == localizationsum ) return;

    const double LocalizationConstant( 4.0 / ( pi * Nminus1 ) ); 
    localizationscore = sqrt( LocalizationConstant * localizationsum );
  }


  inline void ResetClusters()
  {
    parent = NULL;

    for( auto lIt( children.begin() ) ; lIt != children.end() ;  ){
      auto lIt_copy = lIt++;
      if( *lIt_copy == this ) continue;
      // std::unique_lock< std::mutex > lLock( *(**lIt_copy).mutex );
      (**lIt_copy).children.splice( (**lIt_copy).children.end() , children , lIt_copy );
    } 

    ClusterSize = 1;
    ClusterScore = 0.0;
    for( auto& i : ClusterParams ) i.Reset( x , y );
  }


  inline void Clusterize( const double& a2R2 , const double& aT )
  {
    if( parent ) return;
    if( localizationscore < aT ) return;
  
    parent = this;

    for( auto& i : neighbours )
    {
      for( auto& j : i )
      {
        if( j.first > a2R2 ) break;
        if( j.second->localizationscore > aT ) j.second->ClusterInto( this );
      }
    }
  }


  inline void ClusterInto( Data* aParent )
  {
    if ( parent == this ) return;

    // We already belong to another - they already have all our children and clusterparams - Hand them over to claiming cluster and return
    if ( parent ) return parent-> ClusterInto( aParent );

    // If the existing parent is smaller than us, merge it into us
    if( aParent->children.size() < children.size() ) return aParent->ClusterInto( this );

    // std::unique_lock< std::mutex > lLock1( *mutex );
    // std::unique_lock< std::mutex > lLock2( *(aParent->mutex) );

    // Set the claiming cluster as parent
    parent = aParent;

    // Hand over all our children
    for( auto& child : children ) {
      // std::unique_lock< std::mutex > lLock1( *(child->mutex) );
      child->parent = aParent;
    }
    parent->children.splice( parent->children.end() , children );

    // Add our cluster params into the parent's
    for( auto lIt( ClusterParams.begin() ) , lIt2( parent->ClusterParams.begin() ) ; lIt != ClusterParams.end() ; ++lIt , ++lIt2 ) *lIt2 += *lIt;
  }

  inline double S2( const std::size_t& index , const double& nu_bar_x , const double& nu_bar_y ) const
  {
    auto lX( x - nu_bar_x ) , lY( y - nu_bar_y );
    auto d2( (lX*lX) + (lY*lY) );
    return ClusterParams[index].w * d2;
  }


  inline void UpdateClusterScore()
  {
    if( ! (children.size() > ClusterSize) ) return; // We were not bigger than the previous size when we were evaluated - score is still valid

    ClusterSize = children.size();
    ClusterScore = 0.0;

    static constexpr double pi = atan(1)*4;
    static constexpr double logA( -1.0 * log( 4.0 ) );
    static constexpr double log2pi( log( 2*pi ) );

    const double n( children.size() );
    const double pi_term( log2pi * (1.0-n) ); // IF THE 1-n is integers, this gives the wrong answer!!!! WHY????????


    for( std::size_t i(0) ; i!=Parameters.sigmacount() ; ++i )
    {
      auto sqrt_n_tilde( sqrt( ClusterParams[i].n_tilde ) ) , inv_n_tilde( 1.0 / ClusterParams[i].n_tilde );
      auto nu_bar_x ( ClusterParams[i].nu_bar_x * inv_n_tilde ) , nu_bar_y ( ClusterParams[i].nu_bar_y * inv_n_tilde );

      double S2( 0.0 );
      for( auto& j : children ) S2 += j->S2( i , nu_bar_x , nu_bar_y );

      double phi_x( ROOT::Math::normal_cdf( sqrt_n_tilde * (1.0-nu_bar_x) ) - ROOT::Math::normal_cdf( sqrt_n_tilde * (-1.0-nu_bar_x) ) ) , phi_y( ROOT::Math::normal_cdf( sqrt_n_tilde * (1.0-nu_bar_y) ) - ROOT::Math::normal_cdf( sqrt_n_tilde * (-1.0-nu_bar_y) ) );
      auto sum = logA + pi_term - log( ClusterParams[i].n_tilde ) + ClusterParams[i].sum_logw - ( 0.5 * S2 ) + log( phi_x ) + log( phi_y ) + Parameters.log_probability_sigma( i );
    }

      // double sum( 0 );

    //   // std::cout << "==================================" << std::endl;
    //   // PRINT( logA );
    //   // PRINT( pi_term );
    //   // PRINT( - log( n_tilde ) );
    //   // PRINT( sum_log_w );
    //   // PRINT( - ( 0.5 * S2 ) );
    //   // PRINT( log( phi_x ) );
    //   // PRINT( log( phi_y ) );
    //   // PRINT( log_p_sigma[ i ] );
    //   // std::cout << "---------------" << std::endl;
    //   // PRINT( sum );

    //   return sum;
    // };    

    // auto Integrands = IntegrandExpr | range( Parameters.sigmacount() );
    // auto& Max = *std::max_element( Integrands.begin() , Integrands.end() );

    // double integral( 0.0 );
    // for( auto& i : Integrands ) integral += exp( i - Max );
    // integral *= Parameters.sigmaspacing();  
    //return log( integral ) + Max; 

  }


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



