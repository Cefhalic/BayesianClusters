
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

