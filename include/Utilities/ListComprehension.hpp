#pragma once

#include <numeric>
#include <functional>
#include <algorithm>
#include <vector>

//! Super nerd template magic emulating list comprehension for function with return type
//! \tparam tContainer A container type
//! \tparam tExpr      A function-call type
//! \tparam T          Template magic to determine the type of the data in the container
//! \tparam U          Template magic to determine the return-type of the function, given the type of the data in the container 
//! \param  aExpr      A function-call to be applied to each element of the container
//! \param  aContainer A container holding the arguments to be fed to the expression
template< typename tContainer , typename tExpr , typename T = typename std::remove_reference<tContainer>::type::value_type , typename U = decltype( std::declval<tExpr>().operator()( std::declval<T>() ) ) >
inline 
typename std::enable_if< not std::is_same<U, void>::value, std::vector< U > >::type
operator| ( tExpr&& aExpr , tContainer&& aContainer )
{
  std::vector< U > lRet;
  lRet.reserve( aContainer.size() );
  std::transform( aContainer.begin() , aContainer.end() , std::back_inserter(lRet) , aExpr );
  return lRet;
}


//! Super nerd template magic emulating list comprehension for function with void return type
//! \tparam tContainer A container type
//! \tparam tExpr      A function-call type
//! \tparam T          Template magic to determine the type of the data in the container
//! \tparam U          Template magic to determine the return-type of the function, given the type of the data in the container 
//! \param  aExpr      A function-call to be applied to each element of the container
//! \param  aContainer A container holding the arguments to be fed to the expression
template< typename tContainer , typename tExpr , typename T = typename std::remove_reference<tContainer>::type::value_type , typename U = decltype( std::declval<tExpr>().operator()( std::declval<T>() ) ) >
inline 
typename std::enable_if< std::is_same<U, void>::value, void >::type
operator| ( tExpr&& aExpr , tContainer&& aContainer )
{
  std::transform( aContainer.begin() , aContainer.end() , aExpr );
}

//! Return a container holding copies of a member-variable from each object in a container
//! \tparam tType          A container type
//! \tparam tContainerType Template magic to determine the type of the data in the container
//! \param  aPtr           A pointer-to-member-variable to be applied to each element of the container
//! \param  aContainer     A container holding the objects whose member variable is to be extracted
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



//! Emulate the python range function to generate a vector of ints
//! \param N The number of elements
//! \return A vector of ints
inline std::vector< uint32_t > range( const uint32_t& N )
{
  std::vector< uint32_t > lVec( N );
  std::iota( lVec.begin(), lVec.end(), 0 );
  return lVec;
}


// template< typename T >
// struct Construct
// {
  // template< typename ... Args > Construct( Args... args ) : mBuilder( std::bind( &Construct::constructor<Args...> , args... ) ) {} // Lambda capture of parameter packs doesn't work in C++11 - use std::bind instead
  // template< typename ... Args > static T constructor( Args... args ) { return T( args... ); } 
  // template< typename U > inline T operator() ( const U& ){ return mBuilder(); }
  // std::function< T() > mBuilder;
// };

