

/* ===== C++ ===== */
#include <math.h>
#include <stdlib.h>
// #include <iostream>

/* ===== For Root ===== */
#include "Math/ProbFunc.h" 
#include "Math/Interpolator.h" 

/* ===== Cluster sources ===== */
#include "Cluster_GlobalVars.hpp"
#include "Cluster_Data.hpp"

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



Data::ClusterParameter::ClusterParameter( const double& aW ) : w(aW) , logw( log(aW) )
{}
    
void Data::ClusterParameter::Reset( const double& aX , const double& aY)
{
  A = w;
  Bx = w * aX;
  By = w * aY;
  Cx = Bx * aX;
  Cy = By * aY;
  // n_tilde = w;
  sum_logw = logw;
  // nu_bar_x = w * aX;
  // nu_bar_y = w * aY;
}

Data::ClusterParameter& Data::ClusterParameter::operator+= ( const Data::ClusterParameter& aOther )
{
  A += aOther.A;
  Bx += aOther.Bx;
  By += aOther.By;
  Cx += aOther.Cx;
  Cy += aOther.Cy;
  // n_tilde += aOther.n_tilde;
  sum_logw += aOther.sum_logw;
  // nu_bar_x += aOther.nu_bar_x;
  // nu_bar_y += aOther.nu_bar_y;
  return *this;
}


Data::Data( const std::size_t& aI , const double& aX , const double& aY , const double& aS ) : 
mutex( new std::mutex() ),
i(aI) ,
x(aX) , y(aY) , s(aS) , r( sqrt( (aX*aX) + (aY*aY) ) ), phi( atan2( aY , aX ) ),
eX( 1 - fabs( aX ) ) , eY( 1 - fabs( aY ) ) , 
localizationsum( 0.0 ) , localizationscore( 0.0 ),
neighbourit( neighbours[0].end() ),
parent( NULL ), children(),
ClusterSize( 0 ) , ClusterScore( 0.0 )
{
  for( auto& sig2 : Parameters.sigmabins2() ) ClusterParams.emplace_back( 1.0 / ( (aS*aS) + sig2 ) );

  neighbours[0].reserve( 1024 );
  neighbours[1].reserve( 1024 );    
}


bool Data::operator< ( const Data& aOther ) const
{ 
  return r < aOther.r; 
}

double Data::dR2( const Data& aOther ) const
{
  double dX( x - aOther.x ), dY( y - aOther.y );
  return ( dX*dX ) + ( dY*dY );
}

double Data::dR( const Data& aOther ) const
{
  return sqrt( dR2( aOther ) );
}


__attribute__((flatten))
void Data::PopulateNeighbours( std::vector<Data>::iterator aPlusIt , const std::vector<Data>::iterator& aPlusEnd , std::vector<Data>::reverse_iterator aMinusIt , const std::vector<Data>::reverse_iterator& aMinusEnd )
{
  children = { this };

  auto& neighbours0 = neighbours[0];
  auto& neighbours1 = neighbours[1];

  auto dphi = Parameters.maxR() / ( r - Parameters.maxR() );
  auto dphi2 = Parameters.max2R() / ( r - Parameters.max2R() );


  // Iterate over other hits and populate the neighbour list
  for( ; aPlusIt != aPlusEnd ; aPlusIt++ )
  {
    if( ( aPlusIt->r - r ) > Parameters.maxR() ) break; // aPlusIt is always further out than curent 
    if( fabs( aPlusIt->phi - phi ) > dphi ) continue;
    double ldR2 = dR2( *aPlusIt );
    if( ldR2 < Parameters.maxR2() ) neighbours0.push_back( std::make_pair( ldR2 , &*aPlusIt ) );
  }

  for( ; aMinusIt != aMinusEnd ; aMinusIt++ )
  {
    if( ( r - aMinusIt->r ) > Parameters.maxR() ) break; // curent is always further out than aMinusIn
    if( fabs( aPlusIt->phi - phi ) > dphi ) continue;
    double ldR2 = dR2( *aMinusIt );    
    if( ldR2 < Parameters.maxR2() ) neighbours0.push_back( std::make_pair( ldR2 , &*aMinusIt ) );
  }

  std::sort( neighbours0.begin() , neighbours0.end() );
  neighbourit = neighbours0.begin();

  // Iterate over other hits and populate the neighbour2 list
  for( ; aPlusIt != aPlusEnd ; aPlusIt++ )
  {
    if( ( aPlusIt->r - r ) > Parameters.max2R() ) break; // aPlusIt is always further out than curent  
    if( fabs( aPlusIt->phi - phi ) > dphi2 ) continue;
    double ldR2 = dR2( *aPlusIt );
    if( ldR2 < Parameters.max2R2() ) neighbours1.push_back( std::make_pair( ldR2 , &*aPlusIt ) );
  }

  for( ; aMinusIt != aMinusEnd ; aMinusIt++ )
  {
    if( ( r - aMinusIt->r ) > Parameters.max2R() ) break; // curent is always further out than aMinusIn
    if( fabs( aPlusIt->phi - phi ) > dphi2 ) continue;
    double ldR2 = dR2( *aMinusIt );    
    if( ldR2 < Parameters.max2R2() ) neighbours1.push_back( std::make_pair( ldR2 , &*aMinusIt ) );
  }

  std::sort( neighbours1.begin() , neighbours1.end() );
}


void Data::UpdateLocalization( const double& aR2 , const size_t& Nminus1  )
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


void Data::ResetClusters()
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

// #include <iostream>

__attribute__((flatten))
void Data::Clusterize( const double& a2R2 , const double& aT )
{
  // std::cout << "---------------- " << r << " ----------------" << std::endl;

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


void Data::ClusterInto( Data* aParent )
{
  if ( parent == aParent ) return;

  // std::cout << std::string( recursion , ' ' ) << r << " " << this << " " << parent << " " << aParent << std::endl;

  // We already belong to another - they already have all our children and clusterparams - Hand them over to claiming cluster and return
  if ( parent and parent != this ) return parent-> ClusterInto( aParent );

  // If the existing parent is smaller than us, merge it into us
  if( aParent->children.size() < children.size() ) return aParent->ClusterInto( this );

  std::lock( *mutex , *(aParent->mutex) );

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

  mutex->unlock();
  aParent->mutex->unlock();
}

// double Data::S2( const std::size_t& index , const double& nu_bar_x , const double& nu_bar_y ) const
// {
//   auto lX( x - nu_bar_x ) , lY( y - nu_bar_y );
//   auto d2( (lX*lX) + (lY*lY) );
//   return ClusterParams[index].w * d2;
// }

// #include <sstream>
// #include <iostream>

inline double CDF( const double& aArg )
{
  // Above or below ~8 are indistinguishable from 0 and 1 respectively
  if( aArg > 8.0 ) return 1.0;
  if( aArg < 8.0 ) return 0.0;
  return  ROOT::Math::normal_cdf( aArg );
}


__attribute__((flatten))
void Data::UpdateClusterScore()
{
  if( ! (children.size() > ClusterSize) ) return; // We were not bigger than the previous size when we were evaluated - score is still valid

  ClusterSize = children.size();

  static constexpr double pi = atan(1)*4;
  static constexpr double log2pi = log( 2*pi );

  thread_local static std::vector< double > integral( Parameters.sigmacount() );

  double constant;

  for( std::size_t i(0) ; i!=Parameters.sigmacount() ; ++i )
  {
    double log_sum = 0.0;
    log_sum += double( -log( 4.0 ) ); // Area term
    log_sum += (log2pi * (1.0-children.size()));
    log_sum += ClusterParams[i].sum_logw;
    log_sum -= double( log( ClusterParams[i].A ) ); // A is equivalent to n_tilde

    double sqrt_A( sqrt( ClusterParams[i].A ) ) , inv_A( 1.0 / ClusterParams[i].A );
    double Dx( ClusterParams[i].Bx * inv_A ) , Dy( ClusterParams[i].By * inv_A );
    double Ex( ClusterParams[i].Cx - (ClusterParams[i].Bx * Dx ) ) , Ey( ClusterParams[i].Cy - (ClusterParams[i].By * Dy ) );

    log_sum += 0.5 * (Ex + Ey);

    // We place explicit bounds checks to prevent calls to expensive functions
    double arg1 = CDF( sqrt_A * (1.0-Dx) ) - CDF( sqrt_A * (-1.0-Dx) );
    if( arg1 != 1.0 ) log_sum += log( arg1 );
    double arg2 = CDF( sqrt_A * (1.0-Dy) ) - CDF( sqrt_A * (-1.0-Dy) );
    if( arg2 != 1.0 ) log_sum += log( arg2 );

    log_sum += Parameters.log_probability_sigma( i );

    if( i == 0 ) constant = log_sum;
    integral[i] = exp( log_sum - constant );
  }

  thread_local static ROOT::Math::Interpolator lInt( Parameters.sigmacount() );
  lInt.SetData( Parameters.sigmabins() , integral );

  static const double Lower( Parameters.sigmabins(0) ) , Upper( Parameters.sigmabins(Parameters.sigmacount()-1) );
  ClusterScore = double( log( lInt.Integ( Lower , Upper ) ) ) + constant;  

}


