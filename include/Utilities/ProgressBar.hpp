//! \file ProgressBar.hpp
#pragma once

#include <chrono>
#include <mutex>

//! A utility progress-bar
class ProgressBar
{
public:
  //! Constructor
  //! \param aLabel A description of the task being timed
  //! \param aMax   The number of calls equalling 100%
  ProgressBar( const std::string& aLabel, const uint32_t& aMax );

  //! Destructor
  virtual ~ProgressBar();

  //! Postfix increment
  void operator++ ();

  //! Prefix increment
  void operator++ ( int aDummy /*!< Anonymous argument */ );

private:
  //! Update the screen
  void print();

private:
  //! The size of each increment
  float mBlockSize;
  //! The next threshold at which we will write a block to stdout
  float mNextThreshold ;
  //! The number of times we have incremented
  std::size_t mCount;
  //! A timer for end-of-task stats
  std::chrono::high_resolution_clock::time_point mStart;
  //! A mutex for multi-threaded updates
  std::mutex mMutex;
  //! The label for the start of the line
  std::string mLabel;
  //! The current progress
  std::size_t mPercent;
};



//! A utility code timer
struct ProgressTimer {
  //! Constructor
  //! \param aLabel A description of the task being timed
  //! \param aMax   The number of calls equalling 100%
  ProgressTimer( const std::string& aLabel );

  //! Destructor
  virtual ~ProgressTimer();

  //! A timer for end-of-task stats
  std::chrono::high_resolution_clock::time_point mStart;
};