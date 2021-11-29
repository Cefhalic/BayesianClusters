#pragma once

#include <functional>
#include <mutex>
#include <vector>
#include <thread>

#include "ListComprehension.hpp"


template< typename tContainer , typename T = typename std::remove_reference<tContainer>::type::value_type , typename U = decltype( std::declval<T>().operator()() ) >
typename std::enable_if< not std::is_same<U, void>::value, std::vector< U > >::type
Vectorize( tContainer&& aFunctionQueue )
{
  std::vector< U > lRet( aFunctionQueue.size() );

  std::mutex lQueueMutex , lRetMutex;
  auto lIt( aFunctionQueue.begin() );

  auto Runner = [ &aFunctionQueue , &lIt , &lQueueMutex , &lRet , &lRetMutex ]()
  {
    while (true)
    {
      std::unique_lock<std::mutex> lLock(lQueueMutex);
      if (lIt == aFunctionQueue.end()) return;
      std::size_t lIndex( lIt - aFunctionQueue.begin() );
      auto lFunc = *lIt++;
      lLock.unlock();
      auto lRetVal = lFunc();
      std::unique_lock<std::mutex> lLock2(lRetMutex);
      lRet[lIndex] = lRetVal;
      lLock2.unlock();
    }
  };

  std::vector<std::thread> lThreadPool( Construct< std::thread >( Runner ) | range( std::thread::hardware_concurrency() ) );
  for ( auto& i : lThreadPool ) i.join();    

  return lRet;  
}


template< typename tContainer , typename T = typename std::remove_reference<tContainer>::type::value_type , typename U = decltype( std::declval<T>().operator()() ) >
typename std::enable_if< std::is_same<U, void>::value, void >::type
Vectorize ( tContainer&& aFunctionQueue )
{
  std::mutex lQueueMutex;
  auto lIt( aFunctionQueue.begin() );

  auto Runner = [ &aFunctionQueue , &lIt , &lQueueMutex ]()
  {
    while (true)
    {
      std::unique_lock<std::mutex> lLock(lQueueMutex);
      if (lIt == aFunctionQueue.end()) return;
      auto lFunc = *lIt++;
      lLock.unlock();
      lFunc();
    }
  };

  std::vector<std::thread> lThreadPool( Construct< std::thread >( Runner ) | range( std::thread::hardware_concurrency() ) );
  for ( auto& i : lThreadPool ) i.join();    
}