#pragma once

#include <chrono>

//! A utility progress-bar
struct ProgressBar
{
  //! Constructor
  //! \param aLabel A description of the task being timed
  //! \param aMax   The number of calls equalling 100%
  ProgressBar( const std::string& aLabel , const uint32_t& aMax );
  
  //! Destructor
  virtual ~ProgressBar();

  //! Postfix increment
  void operator++ ();
  
  //! Prefix increment
  void operator++ ( int aDummy /*!< Anonymous argument */ );

  //! The size of each increment
  float mBlockSize;
  //! The next threshold at which we will write a block to stdout
  float mNextThreshold ;
  //! The number of times we have incremented
  std::size_t mCount;
  //! A timer for end-of-task stats
  std::chrono::high_resolution_clock::time_point mStart;  
};


//! A utility code timer
struct ProgressBar2
{
  //! Constructor
  //! \param aLabel A description of the task being timed
  //! \param aMax   The number of calls equalling 100%  
  ProgressBar2( const std::string& aLabel , const uint32_t& aMax );

  //! Destructor
  virtual ~ProgressBar2();

  //! Postfix increment
  void operator++ ();
  
  //! Prefix increment
  void operator++ ( int aDummy /*!< Anonymous argument */ );

  //! A timer for end-of-task stats
  std::chrono::high_resolution_clock::time_point mStart;  
};