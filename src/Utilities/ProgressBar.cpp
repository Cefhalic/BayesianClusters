

#include <iostream>

#include "Utilities/ProgressBar.hpp"

ProgressBar::ProgressBar( const std::string& aLabel , const uint32_t& aMax ) : mBlockSize( aMax / 100.0 ) , mNextThreshold( mBlockSize ) , mCount( 0 ) ,
mStart( std::chrono::high_resolution_clock::now() )
{
  std::cout << std::string( 102 + aLabel.size() , '=' ) << "]\r" << aLabel << " [" << std::flush;
}

ProgressBar::~ProgressBar()
{
  std::cout << "#\n  Completed in " << (std::chrono::duration_cast< std::chrono::milliseconds>( std::chrono::high_resolution_clock::now() - mStart ).count()/1000.0) << " seconds" << std::endl;
}

void ProgressBar::operator++ ()
{
  if( mNextThreshold < mCount++ )
  {
    std::cout << '#' << std::flush;
    mNextThreshold += mBlockSize;
  }
}

void ProgressBar::operator++ ( int ) { operator++ (); }




ProgressBar2::ProgressBar2( const std::string& aLabel , const uint32_t& aMax ) : 
mStart( std::chrono::high_resolution_clock::now() )
{
  std::cout << aLabel << std::endl;  
}

ProgressBar2::~ProgressBar2()
{
  std::cout << "  Completed in " << (std::chrono::duration_cast< std::chrono::milliseconds>( std::chrono::high_resolution_clock::now() - mStart ).count()/1000.0) << " seconds" << std::endl;
}

void ProgressBar2::operator++ (){}
void ProgressBar2::operator++ ( int ) { operator++ (); }
