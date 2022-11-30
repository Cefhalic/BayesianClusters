#include <boost/python.hpp>
using namespace boost::python;

#include "BayesianClustering/Configuration.hpp"
#include "BayesianClustering/Event.hpp"
#include "BayesianClustering/EventProxy.hpp"
#include "BayesianClustering/Cluster.hpp"

#include "BayesianClustering/Data.hpp"
#include "BayesianClustering/DataProxy.hpp"


#include <iostream>


template<typename T>
inline std::vector< T > py_list_to_std_vector( const boost::python::object& aIterable )
{ 
  return std::vector< T >( boost::python::stl_input_iterator< T >( aIterable ) , boost::python::stl_input_iterator< T >() ); 
}

template <class T>
inline boost::python::list std_vector_to_py_list( const std::vector<T>& aVector )
{
  boost::python::list lList;
  for ( const auto& i : aVector ) lList.append( boost::ref( i ) );
  return lList;
}

// ---------------------------------------------------------------------------------------------
// There must be an out-of-the box way, but I can't find it
template< typename U >
struct PyIterator
{
  typename std::vector< U >::const_iterator mIt;
  const typename std::vector< U >::const_iterator mEnd;

  PyIterator( const std::vector< U >& aData ) : mIt( aData.begin() ) , mEnd( aData.end() ) {}

  const U& next()
  {    
    if ( mIt == mEnd ) {
      PyErr_SetString(PyExc_StopIteration, "No more data.");
      throw_error_already_set();
    }

    return *mIt++;
  }

};

template< typename U >
struct PyIterator< U* >
{
  typename std::vector< U* >::const_iterator mIt;
  const typename std::vector< U* >::const_iterator mEnd;

  PyIterator( const std::vector< U* >& aData ) : mIt( aData.begin() ) , mEnd( aData.end() ) {}

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



void ConfigFromVector( const boost::python::object& aList )
{
  Configuration::Instance.FromVector( py_list_to_std_vector< std::string >( aList ) );
}

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



PyIterator<Data*> Cluster_GetIterator( const Cluster& aCluster ) { return PyIterator<Data*>( aCluster.mData ); }

std::size_t Cluster_GetSize( Cluster& aCluster )
{ 
  return aCluster.mData.size();
}


PyIterator<Data> Event_GetIterator( const Event& aEvent ) { return PyIterator<Data>( aEvent.mData ); }

std::size_t Event_GetSize( Event& aEvent )
{ 
  return aEvent.mData.size();
}


boost::python::object Data_GetNearestNeighbour( const Data& aData , const Event& aEvent )
{
  for ( const auto& i : aData.mNeighbours )
  {
    if( i.first > 0 ) return boost::python::make_tuple( i.first , boost::ref( aEvent.mData[ i.second ] ) );
  }

  return boost::python::object();
}


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
    .def( "NearestNeighbour" , &Data_GetNearestNeighbour )      
    ;   

  class_< DataProxy, boost::noncopyable >( "DataProxy", init< Data& >() )
    ;   
}

