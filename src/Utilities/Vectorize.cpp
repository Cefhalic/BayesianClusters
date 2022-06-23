
#include "Utilities/ListComprehension.hpp"

#include "Utilities/Vectorize.hpp"


WrappedThread::WrappedThread(): mMask( uint64_t(0x1) << sInstanceCounter++ ) ,
mThread( &WrappedThread::Runner , this )
{}

WrappedThread::~WrappedThread()
{
  mTerminate = true;
  mConditionVariable.notify_one();  
  mThread.join();
}

void WrappedThread::submit( const std::function< void() >& aFunc )
{
  std::unique_lock<std::mutex> lLock(mMutex);
  mFunc = aFunc;
  sBusy |= mMask;
  mConditionVariable.notify_one();
}

void WrappedThread::wait()
{
  while( WrappedThread::sBusy ){}
}

void WrappedThread::run_and_wait( const std::function< void() >& aFunc )
{
  (aFunc)();
  wait();
}


void WrappedThread::Runner()
{
  while (true)
  {
    std::unique_lock<std::mutex> lLock(mMutex);
    mConditionVariable.wait( lLock, [this]() { return ( sBusy & mMask ) or mTerminate; });
    if( mTerminate ) return;
    // if( !( sBusy & mMask ) ) continue;
    (mFunc)();
    sBusy &= ~mMask;
  }
}

std::uint64_t WrappedThread::sInstanceCounter( 0x0 );
std::atomic< std::uint64_t > WrappedThread::sBusy( 0x0 );

std::vector< std::unique_ptr< WrappedThread > > ThreadPool( []( const int& ){ return std::unique_ptr< WrappedThread >( new WrappedThread() ); } | range( Concurrency ) );
