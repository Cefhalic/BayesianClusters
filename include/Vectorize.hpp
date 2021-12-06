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

  // void submit( const std::function< void( const std::size_t& aIndex ) >& aFunc );

  static void wait();

private:
  void Runner();

  static std::atomic< std::uint64_t > lBusy;
  static std::uint64_t lInstanceCtr;
  const std::uint64_t lMask;

  std::condition_variable m_data_condition;
  std::function< void() > lFunc;
  std::atomic< bool > lTerminate;
  std::mutex lMutex;
  std::thread lThread; // MUST BE LAST!
   
};


const std::size_t Concurrency( std::thread::hardware_concurrency() - 1 );
extern std::vector< std::unique_ptr< WrappedThread > > ThreadPool;


// Syntactic sugar
template< typename tContainer , typename tExpr, typename tContainerType = typename std::remove_reference<tContainer>::type::value_type >
inline void operator|| ( tExpr&& aExpr , tContainer&& aContainer )
{
  auto Thread( ThreadPool.begin() );
  for( std::size_t offset(0) ; offset!=Concurrency ; ++offset , ++Thread ) (**Thread).submit( [ aExpr , &aContainer , offset ](){ for( auto i( aContainer.begin() + offset) ; i<aContainer.end() ; i+=Concurrency ) aExpr( *i ); } );
  WrappedThread::wait();
}

// Syntactic sugar
template< typename tContainer , typename tExpr, typename tContainerType = typename std::remove_reference<tContainer>::type::value_type >
inline void operator&& ( tExpr&& aExpr , tContainer&& aContainer )
{
  const std::size_t lChunksize( ceil( double(aContainer.size()) / Concurrency ) );
  auto Thread( ThreadPool.begin() );
  auto A( aContainer.begin() ) , B( aContainer.begin() + lChunksize );

  for( ; A != aContainer.end() ; ++Thread , A = B , B+=lChunksize ){
    if ( B > aContainer.end() ) B = aContainer.end();
    (**Thread).submit( [ aExpr , A , B ](){ for( auto i( A ) ; i != B ; ++i ) aExpr( *i ); } );
  } 
  WrappedThread::wait();
}
