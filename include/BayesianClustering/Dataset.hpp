#pragma once

/* ===== C++ ===== */
#include <vector>
#include <functional>
#include <string>

/* ===== Cluster sources ===== */
#include "BayesianClustering/Data.hpp"

class RoI;


// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//! A class which holds the raw Dataset data and global parameters
class Dataset
{
public:
  //! Default Constructor
  Dataset();  

  //! Deleted copy constructor
  Dataset( const Dataset& aOther /*!< Anonymous argument */ ) = delete;

  //! Deleted assignment operator
  //! \return Reference to this, for chaining calls
  Dataset& operator= (const Dataset& aOther /*!< Anonymous argument */ ) = delete;

  //! Default move constructor
  Dataset( Dataset&& aOther /*!< Anonymous argument */ ) = default;

  //! Default move-assignment constructor
  //! \return Reference to this, for chaining calls
  Dataset& operator= ( Dataset&& aOther /*!< Anonymous argument */ ) = default;

  //! All the necessary pre-processing to get the Dataset ready for an RT-scan
  void Preprocess();
  
  //! Run the scan
  //! \param aCallback A callback for each RT-scan result
  void ScanRT( const std::function< void( const RoI& , const double& , const double& , std::pair<int,int>  ) >& aCallback  );

  //! Run clusterization for a specific choice of R and T
  //! \param R The R parameter for clusterization
  //! \param T The T parameter for clusterization
  //! \param aCallback A callback for the clusterization results
  void Clusterize( const double& R , const double& T , const std::function< void( const RoI& ) >& aCallback );
  
  //! Load an Dataset from given file
  //! \param aFilename The name of the file to load 
  void LoadCSV( const std::string& aFilename );
  
  //! Save an Dataset to a file
  //! \param aFilename The name of the file to which to save   
  void WriteCSV( const std::string& aFilename );

public:
  //! The collection of raw data points
  std::vector<Data> mData; 
};
// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
