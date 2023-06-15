//! \file ProgressBar.cpp

#include <iostream>
#include <iomanip>

#include "Utilities/ProgressBar.hpp"

ProgressBar::ProgressBar( const std::string& aLabel, const uint32_t& aMax ) : 
  mBlockSize( aMax / 100.0 ), 
  mNextThreshold( mBlockSize ), 
  mCount( 0 ),
  mStart( std::chrono::high_resolution_clock::now() ),
  mPercent( 0 )
{
  if( aLabel.size() > 30 ) mLabel = std::string( aLabel.begin() , aLabel.begin() + 27 ) + "...";
  else                     mLabel = aLabel + std::string( 30 - aLabel.size() , ' ' );

  print();
}

ProgressBar::~ProgressBar()
{
  mPercent = 100;
  print();
  std::cout << std::endl;
}

void ProgressBar::operator++ ()
{
  mMutex.lock();
  if( mNextThreshold < ++mCount ) {
    mPercent++;
    print();
    mNextThreshold += mBlockSize;
  }
  mMutex.unlock();
}

void ProgressBar::operator++ ( int )
{
  operator++ ();
}


void ProgressBar::print ()
{ 
  std::ios_base::fmtflags f( std::cout.flags() );
  std::cout << std::fixed << std::setprecision(3);
  std::cout << "\r" << mLabel << " [" << std::string( mPercent , '#' ) << std::string( 100 - mPercent , '=' ) << "] " << std::setw(10) << (std::chrono::duration_cast< std::chrono::milliseconds>( std::chrono::high_resolution_clock::now() - mStart ).count()/1000.0) << "s" << std::flush;
  std::cout.flags( f );
}



ProgressTimer::ProgressTimer( const std::string& aLabel ) :
  mStart( std::chrono::high_resolution_clock::now() )
{
  std::cout << aLabel << std::flush;
}

ProgressTimer::~ProgressTimer()
{
  std::cout << ": Completed in " << (std::chrono::duration_cast< std::chrono::milliseconds>( std::chrono::high_resolution_clock::now() - mStart ).count()/1000.0) << " seconds" << std::endl;
}
