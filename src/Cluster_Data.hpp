#pragma once

/* ===== C++ ===== */
#include <vector>
#include <list>
#include <math.h>
#include <stdlib.h>
#include <iostream>

/* ===== Cluster sources ===== */
#include "Cluster_GlobalVars.hpp"

/* ===== Local utilities ===== */
#include "ListComprehension.hpp"


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
  Data( const double& aX , const double& aY , const double& aS ) : x(aX) , y(aY) , s(aS) , r( sqrt( (aX*aX) + (aY*aY) ) ), phi( atan2( aY , aX ) ),
  eX( 1 - fabs( aX ) ) , eY( 1 - fabs( aY ) ) , 
  //w( [ aS ]( const double& sig2 ){ return  } | Parameters.sigmabins2() ), log_w( []( const double& w){ return log(w); } | w ),
  localizationsum( 0.0 ) , localizationscore( 0.0 ),
  neighbourit( neighbours.end() ),
  parent( NULL ), children( { this } )
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
      if( ( aPlusIt->r - r ) > Parameters.maxR() ) break; // aPlusIt is always further out than lPlus  
      double ldR2 = dR2( *aPlusIt );
      if( ldR2 < Parameters.maxR2() ) neighbours.push_back( std::make_pair( ldR2 , &*aPlusIt ) );
    }

    for( ; aMinusIt != aMinusEnd ; aMinusIt++ )
    {
      if( ( r - aMinusIt->r ) > Parameters.maxR() ) break; // lMinus is always further out than lMinusIn
      double ldR2 = dR2( *aMinusIt );    
      if( ldR2 < Parameters.maxR2() ) neighbours.push_back( std::make_pair( ldR2 , &*aMinusIt ) );
    }

    std::sort( neighbours.begin() , neighbours.end() );
    neighbourit = neighbours.begin();

    // Iterate over other hits and populate the neighbour2 list
    for( ; aPlusIt != aPlusEnd ; aPlusIt++ )
    {
      if( ( aPlusIt->r - r ) > Parameters.max2R() ) break; // aPlusIt is always further out than lPlus  
      double ldR2 = dR2( *aPlusIt );
      if( ldR2 < Parameters.max2R2() ) neighbours2.push_back( std::make_pair( ldR2 , &*aPlusIt ) );
    }

    for( ; aMinusIt != aMinusEnd ; aMinusIt++ )
    {
      if( ( r - aMinusIt->r ) > Parameters.max2R() ) break; // lMinus is always further out than lMinusIn
      double ldR2 = dR2( *aMinusIt );    
      if( ldR2 < Parameters.max2R2() ) neighbours2.push_back( std::make_pair( ldR2 , &*aMinusIt ) );
    }

    std::sort( neighbours2.begin() , neighbours2.end() );
  }


  inline void UpdateLocalization( const double& aR2 , const size_t& Nminus1  )
  {
    constexpr double pi = atan(1)*4;

    const double last_localizationsum( localizationsum );

    for( ; neighbourit != neighbours.end() ; ++neighbourit )
    { 
      if( neighbourit->first > aR2 ) break;
      double lDist = sqrt( neighbourit->first );
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
      (**lIt_copy).children.splice( (**lIt_copy).children.end() , children , lIt_copy );
    } 

    for( auto& i : ClusterParams ) i.Reset( x , y );
  }


  inline void Clusterize( const double& a2R2 , const double& aT )
  {
    if( parent ) return;
    if( localizationscore < aT ) return;

    Data* lParent( this );
    
    for( auto& j : neighbours )
    {
      if( j.first > a2R2 ) break;
      if( (!j.second->parent) or (j.second->localizationscore < aT) ) continue;
      lParent = j.second->GetParent();
      goto breakpoint;
    }

    for( auto& j : neighbours2 )
    {
      if( j.first > a2R2 ) break;
      if( (!j.second->parent) or (j.second->localizationscore < aT) ) continue;
      lParent = j.second->GetParent();
      goto breakpoint;
    }

    SetParent( lParent );

    breakpoint:

    for( auto& j : neighbours )
    {
      if( j.first > a2R2 ) break;
      if( j.second->localizationscore < aT) continue;
      j.second->SetParent( lParent );
    }

    for( auto& j : neighbours2 )
    {
      if( j.first > a2R2 ) break;
      if( j.second->localizationscore < aT) continue;
      j.second->SetParent( lParent );
    }

  }


  inline Data* GetParent()
  {
    if( parent == this or parent == NULL ) return parent;
    return parent = parent->GetParent();
  }

  inline void SetParent( Data* aParent )
  {
    if( parent == aParent ) return;
    if( parent != this and parent != NULL ) parent->SetParent( aParent );
    else parent = aParent;
    
    if( children.size() )
    {
      parent->children.splice( parent->children.end() , children );
      for( auto lIt( ClusterParams.begin() ) , lIt2( parent->ClusterParams.begin() ) ; lIt != ClusterParams.end() ; ++lIt , ++lIt2 ) *lIt2 += *lIt;
    }
  }


  inline double S2( const std::size_t& index , const double& nu_bar_x , const double& nu_bar_y ) const
  {
    auto lX( x - nu_bar_x ) , lY( y - nu_bar_y );
    auto d2( (lX*lX) + (lY*lY) );
    return ClusterParams[index].w * d2;
  }


  inline double ClusterScore() const
  {
    if( children.empty() ) return std::nan( "" );

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


    //   // double phi_x( ROOT::Math::normal_cdf( sqrt_n_tilde * (1.0-nu_bar_x) ) - ROOT::Math::normal_cdf( sqrt_n_tilde * (-1.0-nu_bar_x) ) ) , phi_y( ROOT::Math::normal_cdf( sqrt_n_tilde * (1.0-nu_bar_y) ) - ROOT::Math::normal_cdf( sqrt_n_tilde * (-1.0-nu_bar_y) ) );
    //   // auto sum = logA + pi_term - log( n_tilde ) + sum_log_w - ( 0.5 * S2 ) + log( phi_x ) + log( phi_y ) + Parameters.log_probability_sigma( i );
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

    return 0.0;
  }


public:
  double x, y, s , r, phi;
  double eX , eY;
  double localizationsum , localizationscore;

  std::vector< std::pair< double , Data* > > neighbours , neighbours2;
  std::vector< std::pair< double , Data* > >::iterator neighbourit;

  Data* parent;
  std::list< Data* > children;

  std::vector< ClusterParameter > ClusterParams;

};



