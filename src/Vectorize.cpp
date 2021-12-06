
#include "ListComprehension.hpp"

#include "Vectorize.hpp"


WrappedThread::WrappedThread(): lMask( uint64_t(0x1) << lInstanceCtr++ ) , //lInstance( lInstanceCtr++ ) , 
lThread( &WrappedThread::Runner , this )
{}

WrappedThread::~WrappedThread()
{
  lTerminate = true;
  lThread.join();
}

void WrappedThread::submit( const std::function< void() >& aFunc )
{
  std::unique_lock<std::mutex> lLock(lMutex);
  lFunc = aFunc;
  lBusy |= lMask;
}

// void WrappedThread::submit( const std::function< void( const std::size_t& aIndex ) >& aFunc )
// {
//   std::unique_lock<std::mutex> lLock(lMutex);
//   lFunc = [ & , aFunc ](){ aFunc( lInstance ); };
//   lBusy |= lMask;    
// }

void WrappedThread::wait()
{
  while( WrappedThread::lBusy ){}
}

void WrappedThread::Runner()
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

std::uint64_t WrappedThread::lInstanceCtr( 0x0 );
std::atomic< std::uint64_t > WrappedThread::lBusy( 0x0 );

std::vector< std::unique_ptr< WrappedThread > > ThreadPool( []( const int& ){ return std::unique_ptr< WrappedThread >( new WrappedThread() ); } | range( Concurrency ) );
