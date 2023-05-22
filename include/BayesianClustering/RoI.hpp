#pragma once

/* ===== C++ ===== */
#include <vector>
#include <functional>
#include <string>

/* ===== Cluster sources ===== */
#include "BayesianClustering/Dataset.hpp"
#include "BayesianClustering/Data.hpp"

class RoIproxy;


// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//! A class which holds the raw RoI data and global parameters
class RoI
{
public:
  //! Default Constructor
  //! \param aFilename The name of the file to load   
  RoI( const Dataset& aDataset );  

  //! Deleted copy constructor
  RoI( const RoI& aOther /*!< Anonymous argument */ ) = delete;

  //! Deleted assignment operator
  //! \return Reference to this, for chaining calls
  RoI& operator= (const RoI& aOther /*!< Anonymous argument */ ) = delete;

  //! Default move constructor
  RoI( RoI&& aOther /*!< Anonymous argument */ ) = default;

  //! Default move-assignment constructor
  //! \return Reference to this, for chaining calls
  RoI& operator= ( RoI&& aOther /*!< Anonymous argument */ ) = default;

  //! All the necessary pre-processing to get the RoI ready for an RT-scan
  void Preprocess();
  
  //! Run the scan
  //! \param aCallback A callback for each RT-scan result
  void ScanRT( const std::function< void( const RoIproxy& , const double& , const double& , std::pair<int,int>  ) >& aCallback  );

  //! Run clusterization for a specific choice of R and T
  //! \param R The R parameter for clusterization
  //! \param T The T parameter for clusterization
  //! \param aCallback A callback for the clusterization results
  void Clusterize( const double& R , const double& T , const std::function< void( const RoIproxy& ) >& aCallback );
  
  // //! Save an RoI to a file
  // //! \param aFilename The name of the file to which to save   
  // void WriteCSV( const std::string& aFilename );

public:
  //! The collection of raw data points
  std::vector<Data> mData; 
};
// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
