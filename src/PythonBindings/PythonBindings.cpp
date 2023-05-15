//! \file PythonBindings.cpp
//! Self-contained sourcefile for producing python-bindings

#include <boost/python.hpp>
using namespace boost::python;

#include "BayesianClustering/Configuration.hpp"
#include "BayesianClustering/Event.hpp"
#include "BayesianClustering/EventProxy.hpp"
#include "BayesianClustering/Cluster.hpp"

#include "BayesianClustering/Data.hpp"
#include "BayesianClustering/DataProxy.hpp"


#include <iostream>

//! Utility function to convert a python list to STL vector
//! \tparam T The object type in the STL container
//! \param aIterable The python list to convert
//! \return The python data in an STL vector
template<typename T>
inline std::vector< T > py_list_to_std_vector( const boost::python::object& aIterable )
{ 
  return std::vector< T >( boost::python::stl_input_iterator< T >( aIterable ) , boost::python::stl_input_iterator< T >() ); 
}

//! Utility function to convert an STL vector to python list
//! \tparam T The object type in the STL container
//! \param aVector The STL vector to convert
//! \return The STL data in a python list
template <class T>
inline boost::python::list std_vector_to_py_list( const std::vector<T>& aVector )
{
  boost::python::list lList;
  for ( const auto& i : aVector ) lList.append( boost::ref( i ) );
  return lList;
}

// ---------------------------------------------------------------------------------------------
//! A python iterator over a C++ container
//! \todo There must be an out-of-the box way, but I can't find it
//! \tparam U The type of the data in the C++ container
template< typename U >
struct PyIterator
{
  //! The current location of the iterator
  typename std::vector< U >::const_iterator mIt;
  //! The end of the underlying container
  const typename std::vector< U >::const_iterator mEnd;

  //! Constructor
  //! \param aData The underlying data to be iterated over
  PyIterator( const std::vector< U >& aData ) : mIt( aData.begin() ) , mEnd( aData.end() ) {}

  //! Return the current value and advance the iterator
  //! \return The current value
  const U& next()
  {    
    if ( mIt == mEnd ) {
      PyErr_SetString(PyExc_StopIteration, "No more data.");
      throw_error_already_set();
    }

    return *mIt++;
  }

};

//! A partial specialization creating a python iterator over a C++ container of pointers
//! \todo There must be an out-of-the box way, but I can't find it
//! \tparam U The type of the pointers in the C++ container
template< typename U >
struct PyIterator< U* >
{
  //! The current location of the iterator
  typename std::vector< U* >::const_iterator mIt;
  //! The end of the underlying container
  const typename std::vector< U* >::const_iterator mEnd;

  //! Constructor
  //! \param aData The underlying data to be iterated over
  PyIterator( const std::vector< U* >& aData ) : mIt( aData.begin() ) , mEnd( aData.end() ) {}

  //! Return the current value and advance the iterator
  //! \return The current value
  const U& next()
  {    
    if ( mIt == mEnd ) {
      PyErr_SetString(PyExc_StopIteration, "No more data.");
      throw_error_already_set();
    }

    return **mIt++;
  }

};
// ---------------------------------------------------------------------------------------------

//! Set the Bayesian-clustering configuration from a vector of arguments
//! \param aList A list of strings to parse as config arguments
void ConfigFromVector( const boost::python::object& aList )
{
  Configuration::Instance.FromVector( py_list_to_std_vector< std::string >( aList ) );
}

//! Run a 1-pass clustering for a specified R & T and pass the results to a callback function
//! \param aCallback A callback to which results are passed
void OneStopGetClusters( const boost::python::object& aCallback )
{ 
  Event lEvent;

  lEvent.Clusterize( Configuration::Instance.ClusterR() , Configuration::Instance.ClusterT() , 
    [&]( const EventProxy& aEventProxy ){

      boost::python::list lClusters;
      boost::python::list lBackground;

      for ( const auto& i : aEventProxy.mClusters  )
      {
        if( !i.mParent ) lClusters.append( boost::ref( i ) );
      }

      for( auto& i : aEventProxy.mData )
      { 
        if( i.mCluster ) i.mCluster->GetParent()->mData.push_back( i.mData );
        else             lBackground.append( boost::ref( i.mData ) );
      }

      aCallback( lClusters , lBackground );
    }
  ); 
}

//! Utility function to get a python iterator over all the data points in a clusters
//! \param aCluster The cluster over which we are iterating
//! \return An iterator object pointing to a member of the cluster
PyIterator<Data*> Cluster_GetIterator( const Cluster& aCluster ) { return PyIterator<Data*>( aCluster.mData ); }

//! Utility function to get the number of data points in a clusters
//! \param aCluster The cluster we are inspecting
//! \return The number of data points in the clusters
std::size_t Cluster_GetSize( Cluster& aCluster )
{ 
  return aCluster.mData.size();
}

//! Utility function to get a python iterator over all the data points in an event
//! \param aEvent The event over which we are iterating
//! \return An iterator object pointing to a member of the event
PyIterator<Data> Event_GetIterator( const Event& aEvent ) { return PyIterator<Data>( aEvent.mData ); }

//! Utility function to get the number of data points in an event
//! \param aEvent The event we are inspecting
//! \return The number of data points in the event
std::size_t Event_GetSize( Event& aEvent )
{ 
  return aEvent.mData.size();
}

// boost::python::object Data_GetNearestNeighbour( const Data& aData , const Event& aEvent )
// {
//   for ( const auto& i : aData.mNeighbours )
//   {
//     if( i.first > 0 ) return boost::python::make_tuple( i.first , boost::ref( aEvent.mData[ i.second ] ) );
//   }

//   return boost::python::object();
// }

//! Boost Python Wrapper providing bindings for our C++ functions
BOOST_PYTHON_MODULE( BayesianClustering )
{

  def( "OneStopGetClusters", &OneStopGetClusters );

	class_< Configuration >( "Configuration" )
    .def( "FromVector" , &ConfigFromVector ).staticmethod("FromVector")
    ;

	class_< Event, boost::noncopyable >( "Event" )
    .def( "__iter__" , &Event_GetIterator )
    .def( "__len__" , &Event_GetSize ) 
    .def( "Preprocess" , &Event::Preprocess )       
    ;

  class_< PyIterator<Data> >( "DataIterator", no_init )
    .def("__next__" , &PyIterator<Data>::next , return_value_policy<reference_existing_object>() )
    ;        

  class_< EventProxy, boost::noncopyable >( "EventProxy", init< Event& >() )
    ;   

  class_< Cluster >( "Cluster" )
    .def( "__iter__" , &Cluster_GetIterator )
    .def( "__len__" , &Cluster_GetSize )    
    ;        

  class_< PyIterator<Data*> >( "DataPtrIterator", no_init )
    .def("__next__" , &PyIterator<Data*>::next , return_value_policy<reference_existing_object>() )
    ;        

  class_< Data, boost::noncopyable >( "Data" , no_init )
    .def_readonly("x", &Data::x)
    .def_readonly("y", &Data::y)
    // .def( "NearestNeighbour" , &Data_GetNearestNeighbour )      
    ;   

  class_< DataProxy, boost::noncopyable >( "DataProxy", init< Data& >() )
    ;   
}

