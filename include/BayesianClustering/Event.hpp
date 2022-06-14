#pragma once

/* ===== C++ ===== */
#include <vector>
#include <functional>

/* ===== Cluster sources ===== */
// #include "BayesianClustering/Configuration.hpp"
#include "BayesianClustering/Data.hpp"

class EventProxy;


// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//! A class which holds the raw event data and global parameters
class Event
{
public:
  //! Default Constructor
  Event();  

  //! Deleted copy constructor
  Event( const Event& aOther /*!< Anonymous argument */ ) = delete;

  //! Deleted assignment operator
  //! \return Reference to this, for chaining calls
  Event& operator = (const Event& aOther /*!< Anonymous argument */ ) = delete;

  //! Default move constructor
  Event( Event&& aOther /*!< Anonymous argument */ ) = default;

  //! Default move-assignment constructor
  //! \return Reference to this, for chaining calls
  Event& operator = ( Event&& aOther /*!< Anonymous argument */ ) = default;

  //! All the necessary pre-processing to get the event ready for an RT-scan
  void Preprocess();
  
  //! Run the scan
  //! \param aCallback A callback for each RT-scan result
  void ScanRT( const std::function< void( const EventProxy& , const double& , const double& ) >& aCallback );

  //! Load an event from given file
  //! \param aFilename The name of the file to load 
  void LoadCSV( const std::string& aFilename );
  
  //! Save an event to a file
  //! \param aFilename The name of the file to which to save   
  void WriteCSV( const std::string& aFilename );

public:
  //! The collection of raw data points
  std::vector<Data> mData; 
};
// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
