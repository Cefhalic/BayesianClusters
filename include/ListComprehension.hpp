#pragma once

#include <functional>
#include <algorithm>
#include <vector>

/* ===== Super nerd template magic emulating list comprehension ===== */
template< typename tContainer , typename tExpr , typename T = typename std::remove_reference<tContainer>::type::value_type , typename U = decltype( std::declval<tExpr>().operator()( std::declval<T>() ) ) >
inline 
typename std::enable_if< not std::is_same<U, void>::value, std::vector< U > >::type
operator| ( tExpr&& aExpr , tContainer&& aContainer )
{
  std::vector< U > lRet;
  std::transform( aContainer.begin() , aContainer.end() , std::back_inserter(lRet) , aExpr );
  return lRet;
}


/* ===== Super nerd template magic emulating list comprehension ===== */
template< typename tContainer , typename tExpr , typename T = typename std::remove_reference<tContainer>::type::value_type , typename U = decltype( std::declval<tExpr>().operator()( std::declval<T>() ) ) >
inline 
typename std::enable_if< std::is_same<U, void>::value, void >::type
operator| ( tExpr&& aExpr , tContainer&& aContainer )
{
  std::transform( aContainer.begin() , aContainer.end() , aExpr );
}

// Handle pointer to member variable
template<typename tContainer, typename tType , typename tContainerType = typename std::remove_reference<tContainer>::type::value_type >
inline std::vector< tType > operator| ( tType tContainerType::* aPtr , tContainer&& aContainer )
{
  return [ aPtr ]( const tContainerType& i ){ return i.*aPtr; } | std::forward< tContainer >( aContainer );
}

// /* ===== Handle Function Pointer ===== */
// template< typename tContainer , typename tRet , typename tContainerType = typename std::remove_reference<tContainer>::type::value_type >
// inline auto operator| ( tRet (*aFnPtr)( const tContainerType& ) , tContainer&& aContainer ) -> std::vector< decltype( aExpr( *aContainer.begin() ) ) >
// {
//   return [ aFnPtr ]( const tContainerType& i ){ return aFnPtr.*( i ); } | std::forward< tContainer >( aContainer );
// }




inline std::vector< uint32_t > range( const uint32_t& N )
{
  std::vector< uint32_t > lVec( N );
  std::iota( lVec.begin(), lVec.end(), 0 );
  return lVec;
}


template< typename T >
struct Construct
{
  template< typename ... Args > Construct( Args... args ) : mBuilder( std::bind( &Construct::constructor<Args...> , args... ) ) {} // Lambda capture of parameter packs doesn't work in C++11 - use std::bind instead
  template< typename ... Args > static T constructor( Args... args ) { return T( args... ); } 
  template< typename U > inline T operator() ( const U& ){ return mBuilder(); }
  std::function< T() > mBuilder;
};

