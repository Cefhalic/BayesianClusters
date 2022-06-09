#pragma once

/* ===== C++ ===== */
#include <math.h>
#include <vector>
#include <functional>

/* ===== Local utilities ===== */
#include "Cluster_GlobalVars.hpp"


#define PRECISION float

class Event;
class EventProxy;
class Cluster;
class Data;
class DataProxy;


// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//! A class which holds the raw event data and global parameters
class Event
{
public:
  //! Default Constructor
  Event();  

  //! Destructor
  Event( const std::string& aFilename );

  //! Deleted copy constructor
  Event( const Event& ) = delete;

  //! Deleted assignment operator
  Event& operator = (const Event& ) = delete;

  //! Default move constructor
  Event( Event&& ) = default;

  //! Default move-assignment constructor
  Event& operator = ( Event&& ) = default;

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
  
  //! A single global copy of the global variables
  static GlobalVars mParameters;
};
// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//! A lightweight wrapper for the event to store clusters for a given scan
class EventProxy
{
public:
  //! Default constructor
  EventProxy( Event& aEvent );

  //! Deleted copy constructor
  EventProxy( const EventProxy& ) = delete;

  //! Deleted assignment operator
  EventProxy& operator = (const EventProxy& ) = delete;

  //! Default move constructor
  EventProxy( EventProxy&& ) = default;

  //! Default move-assignment constructor
  EventProxy& operator = ( EventProxy&& ) = default;

  //! Run validation tests on the clusters
  //! \param R The R of the last run scan
  //! \param T The T of the last run scan  
  void CheckClusterization( const double& R , const double& T );
  
  //! Run an RT-scan
  //! \param aCallback        A callback for each RT-scan result
  //! \param aParallelization The stride with which we will iterate across RT parameters
  //! \param aOffset          The starting point for the strides as we iterate across RT parameters
  void ScanRT( const std::function< void( const EventProxy& , const double& , const double& ) >& aCallback , const uint8_t& aParallelization = 1 , const uint8_t& aOffset = 0 );
  
  //! Update log-probability after a scan
  void UpdateLogScore();

public:
  //! The collection of lightweight data-point wrappers used by this event wrapper
  std::vector< DataProxy > mData;
  
  //! The collection of clusters found by this scan
  std::vector< Cluster > mClusters;

  //! The number of clustered data-points
  std::size_t mClusteredCount;
  
  //! The number of background data-points
  std::size_t mBackgroundCount;
  
  //! The number of non-Null clusters
  std::size_t mClusterCount;
  
  //! The log-probability density associated with the last scan
  double mLogP;
};
// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



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
    double log_score() const;
    
    //! Parameter A defined in the math
    PRECISION A;
    //! Parameter Bx defined in the math
    PRECISION Bx;
    //! Parameter By defined in the math
    PRECISION By;
    //! Parameter C defined in the math
    PRECISION C;
    //! Parameter logF defined in the math
    PRECISION logF;
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



// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//! A class to store the raw data-points
class Data
{
public:
  //! Default constructor
  Data( const PRECISION& aX , const PRECISION& aY , const PRECISION& aS );

  //! Deleted copy constructor
  Data( const Data& ) = delete;

  //! Deleted assignment operator
  Data& operator = (const Data& ) = delete;

  //! Default move constructor
  Data( Data&& ) = default;

  //! Default move-assignment constructor
  Data& operator = ( Data&& ) = default;

  //! Comparison operator for sorting data-points by distance from the origin
  //! \return Whether this data-point is closer to the origin than another
  //! \param aOther A data-point to compare against 
  inline bool operator< ( const Data& aOther ) const
  { 
    return r < aOther.r; 
  }

  //! Return the squared-distance of this data-points from another
  //! \return The squared-distance of this data-points from another
  //! \param aOther A data-point to compare against 
  inline PRECISION dR2( const Data& aOther ) const
  {
    PRECISION dX( x - aOther.x ), dY( y - aOther.y );
    return ( dX*dX ) + ( dY*dY );
  }

  //! Return the distance of this data-points from another
  //! \return The distance of this data-points from another
  //! \param aOther A data-point to compare against 
  inline PRECISION dR( const Data& aOther ) const
  {
    return sqrt( dR2( aOther ) );
  }

  //! Return the angle between this data-points and another
  //! \return The angle between this data-points and another
  //! \param aOther A data-point to compare against 
  inline PRECISION dPhi( const Data& aOther ) const
  {
   return fabs( phi - aOther.phi );
  }

  //! All the necessary pre-processing to get this data-point ready for an RT-scan
  //! \param aData  The collection of data-points 
  //! \param aIndex The index of the current data-point
  void Preprocess( std::vector<Data>& aData , const std::size_t& aIndex );

public:
  //! The x-position of the data-point
  PRECISION x;
  //! The y-position of the data-point
  PRECISION y;
  //! The sigma of the data-point  
  PRECISION s;
  //! The squared radial distance of the data-point
  PRECISION r2;
  //! The radial distance of the data-point
  PRECISION r;
  //! The phi-position of the data-point
  PRECISION phi;
  //! The locaalization scores, one per R-bin
  std::vector< PRECISION > mLocalizationScores;
  //! The list of neighbours as a pair of squared-distance and index into the list of points
  std::vector< std::pair< PRECISION , std::size_t > > mNeighbours;
  //! A cluster containing only this data-point 
  Cluster* mProtoCluster;  
};
// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//! A light-weight proxy for the raw data-points
class DataProxy
{
public:
  //! Default constructor
  //! \param aData The data-point for which this is the proxy
  DataProxy( Data& aData );

  //! Deleted copy constructor
  DataProxy( const DataProxy& ) = delete;

  //! Deleted assignment operator
  DataProxy& operator = (const DataProxy& ) = delete;

  //! Default move constructor
  DataProxy( DataProxy&& ) = default;

  //! Default move-assignment constructor
  DataProxy& operator = ( DataProxy&& ) = default;

  //! Get the proxy for the Nth neighbour of this data-point
  //! \return A reference to the neighbour data-proxy
  //! \param aEvent The event-proxy in which we are running
  //! \param aIndex The index of the neighbour we are looking for 
  inline DataProxy& GetNeighbour( EventProxy& aEvent , const std::size_t& aIndex )
  {
    return aEvent.mData.at( aIndex );
  }

  //! Entry point clusterization function - a new cluster will be created
  //! \param a2R2   The clusterization radius
  //! \param aT     The clusterization threshold
  //! \param aEvent The event-proxy in which we are running
  void Clusterize( const PRECISION& a2R2 , const PRECISION& aT , EventProxy& aEvent );
  
  //! Recursive clusterization function
  //! \param a2R2     The clusterization radius
  //! \param aT       The clusterization threshold
  //! \param aEvent   The event-proxy in which we are running  
  //! \param aCluster The cluster we are building
  void Clusterize( const PRECISION& a2R2 , const PRECISION& aT , EventProxy& aEvent , Cluster* aCluster );
  
  //! Get a pointer to this data-proxy's ultimate parent cluster (or null if unclustered
  //! \return A pointer to this data-proxy's ultimate parent cluster  
  inline Cluster* GetCluster()
  {
    if( ! mCluster ) return NULL;
    return mCluster = mCluster->GetParent();
  }

public:
  //! The data-point for which this is the proxy
  Data* mData;
  
  //! This data-proxy's immediate parent cluster
  Cluster* mCluster;
  
  //! Whether this data-point is to be included in the clusterization
  bool mExclude;
};
// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
