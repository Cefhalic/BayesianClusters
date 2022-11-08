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



struct ClusterIterator
{
  const Cluster& mCluster;
  std::vector< Data* >::const_iterator mIt;

  ClusterIterator( const Cluster& aCluster ) : mCluster( aCluster ) , mIt( mCluster.mData.begin() ) {}

  static ClusterIterator Create( const Cluster& aCluster ) { return ClusterIterator( aCluster ); }

  Data& next()
  {    
    if ( mIt == mCluster.mData.end() ) {
      PyErr_SetString(PyExc_StopIteration, "No more data.");
      throw_error_already_set();
    }

    return **mIt++;
  }

};


std::size_t Cluster_GetSize( Cluster& aCluster )
{ 
  return aCluster.mData.size();
}



BOOST_PYTHON_MODULE( BayesianClustering )
{

  def( "OneStopGetClusters", &OneStopGetClusters );

	class_< Configuration >( "Configuration" )
        .def( "FromVector" , &ConfigFromVector ).staticmethod("FromVector")
    ;

	class_< Event, boost::noncopyable >( "Event" )
    ;

  class_< EventProxy, boost::noncopyable >( "EventProxy", init< Event& >() )
    ;   

  class_< Cluster >( "Cluster" )
    .def("__iter__" , &ClusterIterator::Create )
    .def("__len__" , &Cluster_GetSize )    
    ;        

  class_< ClusterIterator >( "ClusterIterator", no_init )
    .def("__next__" , &ClusterIterator::next , return_value_policy<reference_existing_object>() )
    ;        


  class_< Data, boost::noncopyable >( "Data" , no_init )
    .def_readonly("x", &Data::x)
    .def_readonly("y", &Data::y)
    ;   

  class_< DataProxy, boost::noncopyable >( "DataProxy", init< Data& >() )
    ;   
}

