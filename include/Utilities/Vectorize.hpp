//! \file Vectorize.hpp
#pragma once

#include <boost/asio/post.hpp>
#include <boost/asio/thread_pool.hpp>

#include <functional>
#include <cmath>

#include "ListComprehension.hpp"

//! Utility variable for the concurrency
extern std::size_t Nthreads;

//! Syntactic sugar to allow you to interleave parallelize via operator
//! \tparam tContainer A container type
//! \tparam tExpr      A function-call type
//! \tparam tContainerType A SFINAE hack to ensure that the container is a container
//! \param  aExpr      A function-call to be applied to each element of the container
//! \param  aContainer A container holding the arguments to be distributed to the parallelized function calls
template< typename tContainer, typename tExpr, typename tContainerType = typename std::remove_reference<tContainer>::type::value_type >
inline void operator|| ( tExpr&& aExpr, tContainer&& aContainer )
{
  boost::asio::thread_pool ThreadPool( Nthreads );
  for( std::size_t offset(0) ; offset!=Nthreads ; ++offset ) boost::asio::post( ThreadPool, [ aExpr, &aContainer, offset ]() {
    for( auto i( aContainer.begin() + offset ) ; i<aContainer.end() ; i+=Nthreads ) aExpr( *i );
  }  );
  ThreadPool.join();
}

//! Syntactic sugar to allow you to block parallelize via operator
//! \tparam tContainer A container type
//! \tparam tExpr      A function-call type
//! \tparam tContainerType A SFINAE hack to ensure that the container is a container
//! \param  aExpr      A function-call to be applied to each element of the container
//! \param  aContainer A container holding the arguments to be distributed to the parallelized function calls
template< typename tContainer, typename tExpr, typename tContainerType = typename std::remove_reference<tContainer>::type::value_type >
inline void operator&& ( tExpr&& aExpr, tContainer&& aContainer )
{
  boost::asio::thread_pool ThreadPool( Nthreads );

  const std::size_t lChunksize( ceil( double( aContainer.size() ) / Nthreads ) );
  auto A( aContainer.begin() ), B( aContainer.begin() + lChunksize );
  for( ; B < aContainer.end() ; A = B, B+=lChunksize ) boost::asio::post( ThreadPool,  [ aExpr, A, B ]() {
    for( auto i( A ) ; i != B ; ++i ) aExpr( *i );
  } );
  boost::asio::post( ThreadPool, [ aExpr, A, &aContainer ]() {
    for( auto i( A ) ; i != aContainer.end() ; ++i ) aExpr( *i );
  } );
  ThreadPool.join();
}
