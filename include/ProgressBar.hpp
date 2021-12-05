#pragma once

#include <chrono>

/* ===== Utility progress-bar ===== */
struct ProgressBar
{
  ProgressBar( const std::string& aLabel , const uint32_t& aMax );
  virtual ~ProgressBar();

  void operator++ ();
  void operator++ ( int );

  float mBlockSize , mNextThreshold ;
  std::size_t mCount;
  std::chrono::high_resolution_clock::time_point mStart;  
};


/* ===== Utility progress-bar ===== */
struct ProgressBar2
{
  ProgressBar2( const std::string& aLabel , const uint32_t& aMax );
  virtual ~ProgressBar2();

  void operator++ ();
  void operator++ ( int );

  std::chrono::high_resolution_clock::time_point mStart;  
};