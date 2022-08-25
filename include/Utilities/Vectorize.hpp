#pragma once

#include <functional>
#include <mutex>
#include <vector>
#include <thread>
#include <atomic>
#include <condition_variable>

#include "ListComprehension.hpp"

//! A class to wrap a worker-thread
class WrappedThread
{
public:
  //! Default constructor
  WrappedThread();

  //! Deleted copy constructor
  WrappedThread( const WrappedThread& aOther /*!< Anonymous argument */ ) = delete;

  //! Deleted assignment operator
  //! \return Reference to this, for chaining calls
  WrappedThread& operator = (const WrappedThread& aOther /*!< Anonymous argument */ ) = delete;

  //! Default move constructor
  WrappedThread(WrappedThread&& aOther /*!< Anonymous argument */ ) = default;

  //! Default move-assignment constructor
  //! \return Reference to this, for chaining calls
  WrappedThread& operator = (WrappedThread&& aOther /*!< Anonymous argument */ ) = default;

  //! Destructor
  virtual ~WrappedThread();

  //! Submit a job to this thread
  //! \param aFunc The job to run
  void submit( const std::function< void() >& aFunc );
 
  //! Submit a job to the current thread and then wait for all other threads to finish
  //! \param aFunc The job to run
  static void run_and_wait( const std::function< void() >& aFunc );
  
  //! Wait for all other threads to finish  
  static void wait();

private:
  //! The function run by the raw thread 
  void Runner();

  //! An atomic static register keeping track of which threads are busy
  static std::atomic< std::uint64_t > sBusy;
  
  //! A static counter to tell us how many threads are available
  static std::uint64_t sInstanceCounter;
  
  //! A mask indicating which thread we are on
  const std::uint64_t mMask;

  //! A condition variable for talking across threads
  std::condition_variable mConditionVariable;
  
  //! The function call to be run by the thread
  std::function< void() > mFunc;
  
  //! An atomic flag indicating termination
  std::atomic< bool > mTerminate;

  //! The access mutex
  std::mutex mMutex;
  
  //! The raw thread
  std::thread mThread; // MUST BE LAST!
   
};

//! Utility variable for the concurrency
extern std::size_t Nthreads;

//! The pool of all wrapped threads
extern std::vector< std::unique_ptr< WrappedThread > > ThreadPool;


//! Syntactic sugar to allow you to interleave parallelize via operator
//! \tparam tContainer A container type
//! \tparam tExpr      A function-call type
//! \tparam tContainerType A SFINAE hack to ensure that the container is a container
//! \param  aExpr      A function-call to be applied to each element of the container
//! \param  aContainer A container holding the arguments to be distributed to the parallelized function calls
template< typename tContainer , typename tExpr, typename tContainerType = typename std::remove_reference<tContainer>::type::value_type >
inline void operator|| ( tExpr&& aExpr , tContainer&& aContainer )
{
  // Previously the main thread did nothing but loop
  // Now we increase the step-size by one (changing Nthreads - 1 to N in the inner loop), so each thread handles less
  // And handle the remaining entries in the master thread, before entering a loop to check that the children have finished
  auto Thread( ThreadPool.begin() );
  for( std::size_t offset(0) ; offset!=Nthreads - 1 ; ++offset , ++Thread ) (**Thread).submit( [ aExpr , &aContainer , offset ](){ for( auto i( aContainer.begin() + offset) ; i<aContainer.end() ; i+=Nthreads ) aExpr( *i ); } );
  WrappedThread::run_and_wait( [ aExpr , &aContainer ](){ for( auto i( aContainer.begin() + Nthreads - 1 ) ; i<aContainer.end() ; i+=Nthreads ) aExpr( *i ); } );
}

//! Syntactic sugar to allow you to block parallelize via operator
//! \tparam tContainer A container type
//! \tparam tExpr      A function-call type
//! \tparam tContainerType A SFINAE hack to ensure that the container is a container
//! \param  aExpr      A function-call to be applied to each element of the container
//! \param  aContainer A container holding the arguments to be distributed to the parallelized function calls
template< typename tContainer , typename tExpr, typename tContainerType = typename std::remove_reference<tContainer>::type::value_type >
inline void operator&& ( tExpr&& aExpr , tContainer&& aContainer )
{
  const std::size_t lChunksize( ceil( double( aContainer.size() ) / ( Nthreads ) ) );
  auto Thread( ThreadPool.begin() );

  auto A( aContainer.begin() ) , B( aContainer.begin() + lChunksize );
  for( ; B < aContainer.end() ; ++Thread , A = B , B+=lChunksize ){ (**Thread).submit( [ aExpr , A , B ](){ for( auto i( A ) ; i != B ; ++i ) aExpr( *i ); } ); } 
  WrappedThread::run_and_wait( [ aExpr , A , &aContainer ](){ for( auto i( A ) ; i != aContainer.end() ; ++i ) aExpr( *i ); } );
}
