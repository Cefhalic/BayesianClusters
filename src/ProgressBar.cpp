#pragma once

#include <iostream>
#include <chrono>

/* ===== Utility progress-bar ===== */
struct ProgressBar
{
  ProgressBar( const std::string& aLabel , const uint32_t& aMax ) : mBlockSize( aMax / 100.0 ) , mNextThreshold( mBlockSize ) , mCount( 0 ) ,
  mStart( std::chrono::high_resolution_clock::now() )
  {
    std::cout << std::string( 102 + aLabel.size() , '=' ) << "]\r" << aLabel << " [" << std::flush;
  }

  virtual ~ProgressBar()
  {
    std::cout << "#\n  Completed in " << (std::chrono::duration_cast< std::chrono::milliseconds>( std::chrono::high_resolution_clock::now() - mStart ).count()/1000.0) << " seconds" << std::endl;
  }

  inline void operator++ ()
  {
    if( mNextThreshold < mCount++ )
    {
      std::cout << '#' << std::flush;
      mNextThreshold += mBlockSize;
    }
  }

  inline void operator++ ( int ) { operator++ (); }

  float mBlockSize , mNextThreshold ;
  std::size_t mCount;
  std::chrono::high_resolution_clock::time_point mStart;  
};


/* ===== Utility progress-bar ===== */
struct ProgressBar2
{
  ProgressBar2( const std::string& aLabel , const uint32_t& aMax ) : 
  mStart( std::chrono::high_resolution_clock::now() )
  {
    std::cout << aLabel << std::endl;  
  }

  virtual ~ProgressBar2()
  {
    std::cout << "  Completed in " << (std::chrono::duration_cast< std::chrono::milliseconds>( std::chrono::high_resolution_clock::now() - mStart ).count()/1000.0) << " seconds" << std::endl;
  }

  inline void operator++ (){}
  inline void operator++ ( int ) { operator++ (); }

  std::chrono::high_resolution_clock::time_point mStart;  
};