#pragma once

#include <functional>
#include <mutex>
#include <vector>
#include <thread>
#include <atomic>

#include "ListComprehension.hpp"

class WrappedThread
{
public:
  WrappedThread(): lMask( uint64_t(0x1) << lInstanceCtr++ ) , //lInstance( lInstanceCtr++ ) , 
  lThread( &WrappedThread::Runner , this )
  {}

  WrappedThread( const WrappedThread& ) = delete;
  WrappedThread& operator = (const WrappedThread&) = delete;

  WrappedThread(WrappedThread&&) = default;
  WrappedThread& operator = (WrappedThread&&) = default;

  virtual ~WrappedThread()
  {
    lTerminate = true;
    lThread.join();
  }

  inline void submit( const std::function< void() >& aFunc )
  {
    std::unique_lock<std::mutex> lLock(lMutex);
    lFunc = aFunc;
    lBusy |= lMask;
  }

  // inline void submit( const std::function< void( const std::size_t& aIndex ) >& aFunc )
  // {
  //   std::unique_lock<std::mutex> lLock(lMutex);
  //   lFunc = [ & , aFunc ](){ aFunc( lInstance ); };
  //   lBusy |= lMask;    
  // }

  static inline void wait()
  {
    while( WrappedThread::lBusy ){}
  }

private:
  void Runner()
  {
    while (true)
    {
      std::unique_lock<std::mutex> lLock(lMutex);
      if( lTerminate ) return;
      if( !( lBusy & lMask ) ) continue;
      (lFunc)();
      lBusy &= ~lMask;
    }
  }

  static std::atomic< std::uint64_t > lBusy;
  static std::uint64_t lInstanceCtr;
  const std::uint64_t lMask;
  // const std::uint64_t lInstance;

  std::function< void() > lFunc;
  std::atomic< bool > lTerminate;
  std::mutex lMutex;
  std::thread lThread; // MUST BE LAST!
   
};

std::uint64_t WrappedThread::lInstanceCtr( 0x0 );
std::atomic< std::uint64_t > WrappedThread::lBusy( 0x0 );


const std::size_t Concurrency( std::thread::hardware_concurrency() - 1 );
std::vector< std::unique_ptr< WrappedThread > > ThreadPool( []( const int& ){ return std::unique_ptr< WrappedThread >( new WrappedThread() ); } | range( Concurrency ) );


// Syntactic sugar
template< typename tContainer , typename tExpr >
inline void operator|| ( tExpr&& aExpr , tContainer&& aContainer )
{
  auto Thread( ThreadPool.begin() );
  for( std::size_t offset(0) ; offset!=Concurrency ; ++offset , ++Thread ) (**Thread).submit( [ &aExpr , &aContainer , offset ](){ for( auto i( aContainer.begin() + offset) ; i<aContainer.end() ; i+=Concurrency ) aExpr( *i ); } );
  WrappedThread::wait();
}
