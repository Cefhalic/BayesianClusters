#pragma once

/* ===== C++ ===== */
#include <vector>

/* ===== Cluster sources ===== */
#include "BayesianClustering/Precision.hpp"

class Data;


// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//! A class representing a cluster
class Cluster
{
public:

  //! A struct representing the cluster parameters
  struct Parameter
  {
    //! Default constructor
    Parameter();
    
    //! Add another set of parameters to this set
    //! \param aOther Another set of parameters to add to this set
    //! \return Reference to this, for chaining calls
    Parameter& operator+= ( const Parameter& aOther );
    
    //! Convert the parameters to a log-probability
    //! \return the log-probability of this set of cluster parameters
    // TODO: rewrie this
    double log_score() const;
    
    // //! Parameter A defined in the math
    // PRECISION A;
    // //! Parameter Bx defined in the math
    // PRECISION Bx;
    // //! Parameter By defined in the math
    // PRECISION By;
    // //! Parameter C defined in the math
    // PRECISION C;
    // //! Parameter logF defined in the math
    // PRECISION logF;

    //nTilde as defined in the paper
    PRECISION nTilde;

    PRECISION nuBarX;
    PRECISION nuBarY;
  }; 

  //! Default constructor
  Cluster();
  
  //! Construct a cluster from a single data-point
  //! \param aData A data-point with which to initialize the cluster
  Cluster( const Data& aData );

  //! Add another cluster to this one
  //! \param aOther Another cluster of parameters to add to this one
  //! \return Reference to this, for chaining calls
  Cluster& operator+= ( const Cluster& aOther );

  //! Get a pointer to this cluster's ultimate parent
  //! \return A pointer to this cluster's ultimate parent
  Cluster* GetParent();

  //! Update log-probability after a scan
  void UpdateLogScore();

  // rewriting the above function
  void EvalLogScore();

public:
  //! The collection of parameters, each corresponding to a different sigma hypothesis 
  std::vector< Parameter > mParams;
  
  //! The number of points in the current cluster
  std::size_t mClusterSize;
  
  //! The number of points in the cluster on the previous scan iteration
  std::size_t mLastClusterSize;
  
  //! The log-probability of the current cluster
  PRECISION mClusterScore;
  
  //! A pointer to the immediate parent of the current cluster
  Cluster* mParent;
};


// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


