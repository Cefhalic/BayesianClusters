#pragma once

#include <functional>
#include <mutex>
#include <vector>
#include <thread>
#include <atomic>
#include <condition_variable>

#include "ListComprehension.hpp"

class WrappedThread
{
public:
  WrappedThread();

  WrappedThread( const WrappedThread& ) = delete;
  WrappedThread& operator = (const WrappedThread&) = delete;

  WrappedThread(WrappedThread&&) = default;
  WrappedThread& operator = (WrappedThread&&) = default;

  virtual ~WrappedThread();

  void submit( const std::function< void() >& aFunc );
 
  static void run_and_wait( const std::function< void() >& aFunc );
  static void wait();

private:
  void Runner();

  static std::atomic< std::uint64_t > sBusy;
  static std::uint64_t sInstanceCounter;
  const std::uint64_t mMask;

  std::condition_variable mConditionVariable;
  std::function< void() > mFunc;
  std::atomic< bool > mTerminate;
  std::mutex mMutex;
  std::thread mThread; // MUST BE LAST!
   
};


const std::size_t Concurrency( std::thread::hardware_concurrency() - 1 );
extern std::vector< std::unique_ptr< WrappedThread > > ThreadPool;


// Syntactic sugar
template< typename tContainer , typename tExpr, typename tContainerType = typename std::remove_reference<tContainer>::type::value_type >
inline void operator|| ( tExpr&& aExpr , tContainer&& aContainer )
{
  // Previously the main thread did nothing but loop
  // Now we increase the step-size by one (changing Concurrency to N in the inner loop), so each thread handles less
  // And handle the remaining entries in the master thread, before entering a loop to check that the children have finished
  auto Thread( ThreadPool.begin() );
  const auto N = Concurrency + 1;
  for( std::size_t offset(0) ; offset!=Concurrency ; ++offset , ++Thread ) (**Thread).submit( [ aExpr , &aContainer , offset , &N ](){ for( auto i( aContainer.begin() + offset) ; i<aContainer.end() ; i+=N ) aExpr( *i ); } );
  WrappedThread::run_and_wait( [ aExpr , &aContainer , &N ](){ for( auto i( aContainer.begin() + Concurrency ) ; i<aContainer.end() ; i+=N ) aExpr( *i ); } );
}

// Syntactic sugar
template< typename tContainer , typename tExpr, typename tContainerType = typename std::remove_reference<tContainer>::type::value_type >
inline void operator&& ( tExpr&& aExpr , tContainer&& aContainer )
{
  const std::size_t lChunksize( ceil( double( aContainer.size() ) / Concurrency ) );
  auto Thread( ThreadPool.begin() );

  auto A( aContainer.begin() ) , B( aContainer.begin() + lChunksize );
  for( ; B < aContainer.end() ; ++Thread , A = B , B+=lChunksize ){ (**Thread).submit( [ aExpr , A , B ](){ for( auto i( A ) ; i != B ; ++i ) aExpr( *i ); } ); } 
  WrappedThread::run_and_wait( [ aExpr , A , &aContainer ](){ for( auto i( A ) ; i != aContainer.end() ; ++i ) aExpr( *i ); } );
}
