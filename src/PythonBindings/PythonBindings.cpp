//! \file PythonBindings.cpp
//! Self-contained sourcefile for producing python-bindings

/* ===== C++ libraries ===== */
#include <iostream>

/* ===== BOOST libraries ===== */
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
using namespace boost::python;

/* ===== Cluster sources ===== */
#include "BayesianClustering/API.hpp"
#include "BayesianClustering/LocalizationFile.hpp"
#include "BayesianClustering/Configuration.hpp"
// #include "BayesianClustering/RoI.hpp"
#include "BayesianClustering/RoIproxy.hpp"
// #include "BayesianClustering/Cluster.hpp"
// #include "BayesianClustering/RoI.hpp"
#include "BayesianClustering/Data.hpp"
// #include "BayesianClustering/DataProxy.hpp"


//! Wrapper to automatically extract RoI, run scan and apply a full python callback
//! \param aInFile     The name of the localization file
//! \param aScanConfig The configuration for the scan
//! \param aCallback   The full python callback to be applied 
void _AutoRoi_Scan_FullCallback_( const std::string& aInFile , const ScanConfiguration& aScanConfig , const boost::python::object& aCallback )
{
  AutoRoi_Scan_FullCallback( aInFile , aScanConfig, [&]( RoIproxy& aRoIproxy, const double& aR , const double& aT ){ aCallback( aRoIproxy , aR , aT ); } );
}

//! Wrapper to automatically extract RoI, run scan and apply a simple python callback
//! \param aInFile     The name of the localization file
//! \param aScanConfig The configuration for the scan
//! \param aCallback   The simple python callback to be applied
void _AutoRoi_Scan_SimpleCallback_( const std::string& aInFile , const ScanConfiguration& aScanConfig , const boost::python::object& aCallback )
{
  AutoRoi_Scan_SimpleCallback( aInFile , aScanConfig, [&]( const std::vector< ScanEntry >& aScanResults ){ aCallback( aScanResults ); } );
}

//! Wrapper to manually extract RoI, run scan and apply a full python callback
//! \param aInFile     The name of the localization file
//! \param aManualRoI  The manually-specified RoI window
//! \param aScanConfig The configuration for the scan
//! \param aCallback   The full python callback to be applied
void _ManualRoi_Scan_FullCallback_( const std::string& aInFile , const ManualRoI& aManualRoI , const ScanConfiguration& aScanConfig , const boost::python::object& aCallback )
{
  ManualRoi_Scan_FullCallback( aInFile , aManualRoI , aScanConfig, [&]( RoIproxy& aRoIproxy, const double& aR , const double& aT ){ aCallback( aRoIproxy , aR , aT ); } );
}

//! Wrapper to manually extract RoI, run scan and apply a simple python callback
//! \param aInFile     The name of the localization file
//! \param aManualRoI  The manually-specified RoI window
//! \param aScanConfig The configuration for the scan
//! \param aCallback   The simple python callback to be applied
void _ManualRoi_Scan_SimpleCallback_( const std::string& aInFile , const ManualRoI& aManualRoI , const ScanConfiguration& aScanConfig , const boost::python::object& aCallback )
{
  ManualRoi_Scan_SimpleCallback( aInFile , aManualRoI , aScanConfig , [&]( const std::vector< ScanEntry >& aScanResults ){ aCallback( aScanResults ); } );
}


//! Helper Macro to simplify defining functions
//! \param X The function being defined
#define          FN( X , ... ) def( #X , &X    , args( __VA_ARGS__ ) )

//! Helper Macro to simplify defining functions with python callbacks
//! \param X The function being defined
#define CALLBACK_FN( X , ... ) def( #X , &_##X##_ , args( __VA_ARGS__ ) )

//! Boost Python Wrapper providing bindings for our C++ functions
BOOST_PYTHON_MODULE( BayesianClustering )
{

  class_< ManualRoI >( "ManualRoI" )
    .def_readwrite( "x"      , &ManualRoI::x )
    .def_readwrite( "y"      , &ManualRoI::y )
    .def_readwrite( "width"  , &ManualRoI::width )
    .def_readwrite( "height" , &ManualRoI::height );

  class_< ScanConfiguration, boost::noncopyable >( "ScanConfiguration" , init< const std::string& >( arg( "aCfgFile" ) ) )
    .def( init< const std::size_t& , const double& , const double& , const std::function< double( const double& ) >&  ,
                const std::size_t& , const double& , const double&  ,
                const std::size_t& , const double& , const double&  ,
                const double&  , const double&  >
                ( args( "aSigmacount", "aSigmaMin", "aSigmaMax", "aInterpolator", "aRbins", "aMinScanR", "aMaxScanR", "aTbins", "aMinScanT", "aMaxScanT", "aPB", "aAlpha" ) ) ); 

  class_< ScanEntry >( "ScanEntry" )  
      .def_readwrite( "r" , &ScanEntry::r )
      .def_readwrite( "t" , &ScanEntry::t )
      .def_readwrite( "score" , &ScanEntry::score );

  class_< std::vector<ScanEntry> >( "VectorScanEntry" )
      .def( vector_indexing_suite< std::vector<ScanEntry> >() );

  CALLBACK_FN( AutoRoi_Scan_FullCallback , "aInFile" , "aScanConfig" , "aCallback" );
  CALLBACK_FN( AutoRoi_Scan_SimpleCallback , "aInFile" , "aScanConfig" , "aCallback" );
  FN( AutoRoi_Scan_ToJson , "aInFile" , "aScanConfig" , "aOutFile" );

  // CALLBACK_FN( AutoRoi_Cluster_Callback );

  CALLBACK_FN( ManualRoi_Scan_FullCallback , "aInFile" , "aManualRoI" , "aScanConfig" , "aCallback" );
  CALLBACK_FN( ManualRoi_Scan_SimpleCallback , "aInFile" , "aManualRoI" , "aScanConfig" , "aCallback" );
  FN( ManualRoi_Scan_ToJson , "aInFile" , "aManualRoI" , "aScanConfig" , "aOutFile" );

  // CALLBACK_FN( ManualRoi_Cluster_Callback );

}

#undef CALLBACK_FN
#undef FN


// // //! Utility function to convert a python list to STL vector
// // //! \tparam T The object type in the STL container
// // //! \param aIterable The python list to convert
// // //! \return The python data in an STL vector
// // template<typename T>
// // inline std::vector< T > py_list_to_std_vector( const boost::python::object& aIterable )
// // {
// //   return std::vector< T >( boost::python::stl_input_iterator< T >( aIterable ) , boost::python::stl_input_iterator< T >() );
// // }

// // //! Utility function to convert an STL vector to python list
// // //! \tparam T The object type in the STL container
// // //! \param aVector The STL vector to convert
// // //! \return The STL data in a python list
// // template <class T>
// // inline boost::python::list std_vector_to_py_list( const std::vector<T>& aVector )
// // {
// //   boost::python::list lList;
// //   for ( const auto& i : aVector ) lList.append( boost::ref( i ) );
// //   return lList;
// // }

// // // ---------------------------------------------------------------------------------------------
// // //! A python iterator over a C++ container
// // //! \todo There must be an out-of-the box way, but I can't find it
// // //! \tparam U The type of the data in the C++ container
// // template< typename U >
// // struct PyIterator
// // {
// //   //! The current location of the iterator
// //   typename std::vector< U >::const_iterator mIt;
// //   //! The end of the underlying container
// //   const typename std::vector< U >::const_iterator mEnd;

// //   //! Constructor
// //   //! \param aData The underlying data to be iterated over
// //   PyIterator( const std::vector< U >& aData ) : mIt( aData.begin() ) , mEnd( aData.end() ) {}

// //   //! Return the current value and advance the iterator
// //   //! \return The current value
// //   const U& next()
// //   {
// //     if ( mIt == mEnd ) {
// //       PyErr_SetString(PyExc_StopIteration, "No more data.");
// //       throw_error_already_set();
// //     }

// //     return *mIt++;
// //   }

// // };

// // //! A partial specialization creating a python iterator over a C++ container of pointers
// // //! \todo There must be an out-of-the box way, but I can't find it
// // //! \tparam U The type of the pointers in the C++ container
// // template< typename U >
// // struct PyIterator< U* >
// // {
// //   //! The current location of the iterator
// //   typename std::vector< U* >::const_iterator mIt;
// //   //! The end of the underlying container
// //   const typename std::vector< U* >::const_iterator mEnd;

// //   //! Constructor
// //   //! \param aData The underlying data to be iterated over
// //   PyIterator( const std::vector< U* >& aData ) : mIt( aData.begin() ) , mEnd( aData.end() ) {}

// //   //! Return the current value and advance the iterator
// //   //! \return The current value
// //   const U& next()
// //   {
// //     if ( mIt == mEnd ) {
// //       PyErr_SetString(PyExc_StopIteration, "No more data.");
// //       throw_error_already_set();
// //     }

// //     return **mIt++;
// //   }

// // };
// // // ---------------------------------------------------------------------------------------------

// // //! Set the Bayesian-clustering configuration from a vector of arguments
// // //! \param aList A list of strings to parse as config arguments
// // void ConfigFromVector( const boost::python::object& aList )
// // {
// //   CurrentConfiguration().FromVector( py_list_to_std_vector< std::string >( aList ) );
// // }

// // //! Run a 1-pass clustering for a specified R & T and pass the results to a callback function
// // //! \param aCallback A callback to which results are passed
// // void OneStopGetClusters( const std::string& aFilename , const boost::python::object& aCallback )
// // {
// //   RoI lRoI = LoadLocalizationFile( aFilename );
// //   RoI lRoI( lRoI );

// //   Configuration::Instance.SetRBins( 0 , 0 , Configuration::Instance.ClusterR() );

// //   lRoI.Clusterize( Configuration::Instance.ClusterR() , Configuration::Instance.ClusterT() ,
// //     [&]( const RoIproxy& aRoIproxy ){

// //       boost::python::list lClusters;
// //       boost::python::list lBackground;

// //       for ( const auto& i : aRoIproxy.mClusters  )
// //       {
// //         if( !i.mParent ) lClusters.append( boost::ref( i ) );
// //       }

// //       for( auto& i : aRoIproxy.mData )
// //       {
// //         if( i.mCluster ) i.mCluster->GetParent()->mData.push_back( i.mData );
// //         else             lBackground.append( boost::ref( i.mData ) );
// //       }

// //       aCallback( lClusters , lBackground );
// //     }
// //   );
// // }

// // //! Utility function to get a python iterator over all the data points in a clusters
// // //! \param aCluster The cluster over which we are iterating
// // //! \return An iterator object pointing to a member of the cluster
// // PyIterator<Data*> Cluster_GetIterator( const Cluster& aCluster ) { return PyIterator<Data*>( aCluster.mData ); }

// // //! Utility function to get the number of data points in a clusters
// // //! \param aCluster The cluster we are inspecting
// // //! \return The number of data points in the clusters
// // std::size_t Cluster_GetSize( Cluster& aCluster )
// // {
// //   return aCluster.mData.size();
// // }

// // //! Utility function to get a python iterator over all the data points in an RoI
// // //! \param aRoI The RoI over which we are iterating
// // //! \return An iterator object pointing to a member of the RoI
// // PyIterator<Data> RoI_GetIterator( const RoI& aRoI ) { return PyIterator<Data>( aRoI.mData ); }

// // //! Utility function to get the number of data points in an RoI
// // //! \param aRoI The RoI we are inspecting
// // //! \return The number of data points in the RoI
// // std::size_t RoI_GetSize( RoI& aRoI )
// // {
// //   return aRoI.mData.size();
// // }

// // // boost::python::object Data_GetNearestNeighbour( const Data& aData , const RoI& aRoI )
// // // {
// // //   for ( const auto& i : aData.mNeighbours )
// // //   {
// // //     if( i.first > 0 ) return boost::python::make_tuple( i.first , boost::ref( aRoI.mData[ i.second ] ) );
// // //   }

// // //   return boost::python::object();
// // // }



// //   def( "OneStopGetClusters", &OneStopGetClusters );

// //  class_< Configuration >( "Configuration" )
// //     .def( "FromVector" , &ConfigFromVector ).staticmethod("FromVector")
// //     ;

// //  class_< RoI, boost::noncopyable >( "RoI", init<const RoI&>() )
// //     .def( "__iter__" , &RoI_GetIterator )
// //     .def( "__len__" , &RoI_GetSize )
// //     .def( "Preprocess" , &RoI::Preprocess )
// //     ;

// //   class_< PyIterator<Data> >( "DataIterator", no_init )
// //     .def("__next__" , &PyIterator<Data>::next , return_value_policy<reference_existing_object>() )
// //     ;

// //   class_< RoIproxy, boost::noncopyable >( "RoIproxy", init< RoI& >() )
// //     ;

// //   class_< Cluster >( "Cluster" )
// //     .def( "__iter__" , &Cluster_GetIterator )
// //     .def( "__len__" , &Cluster_GetSize )
// //     ;

// //   class_< PyIterator<Data*> >( "DataPtrIterator", no_init )
// //     .def("__next__" , &PyIterator<Data*>::next , return_value_policy<reference_existing_object>() )
// //     ;

// //   class_< Data, boost::noncopyable >( "Data" , no_init )
// //     .def_readonly("x", &Data::x)
// //     .def_readonly("y", &Data::y)
// //     // .def( "NearestNeighbour" , &Data_GetNearestNeighbour )
// //     ;

// //   class_< DataProxy, boost::noncopyable >( "DataProxy", init< Data& >() )
// //     ;

// //   class_< RoI, boost::noncopyable >( "RoI" )
// //     ;
// // }

