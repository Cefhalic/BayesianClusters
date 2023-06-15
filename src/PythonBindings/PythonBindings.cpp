//! \file PythonBindings.cpp
//! Self-contained sourcefile for producing python-bindings

/* ===== C++ libraries ===== */
#include <iostream>
#include <functional>

/* ===== BOOST libraries ===== */
#pragma GCC diagnostic push 
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#pragma GCC diagnostic pop

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
__attribute__((flatten))
void _AutoRoi_Scan_FullCallback_( const std::string& aInFile , const ScanConfiguration& aScanConfig , const boost::python::object& aCallback )
{
  AutoRoi_Scan_FullCallback( aInFile , aScanConfig, [&]( RoIproxy& aRoIproxy, const double& aR , const double& aT ){ aCallback( aRoIproxy , aR , aT ); } );
}

//! Wrapper to automatically extract RoI, run scan and apply a simple python callback
//! \param aInFile     The name of the localization file
//! \param aScanConfig The configuration for the scan
//! \param aCallback   The simple python callback to be applied
__attribute__((flatten))
void _AutoRoi_Scan_SimpleCallback_( const std::string& aInFile , const ScanConfiguration& aScanConfig , const boost::python::object& aCallback )
{
  AutoRoi_Scan_SimpleCallback( aInFile , aScanConfig, [&]( const std::vector< ScanEntry >& aScanResults ){ aCallback( aScanResults ); } );
}

//! Wrapper to automatically extract RoI, clusterize and apply a full python callback
//! \param aInFile     The name of the localization file
//! \param aR          The R value of the clusterizer
//! \param aT          The T value of the clusterizer
//! \param aCallback   The full python callback to be applied
__attribute__((flatten))
void _AutoRoi_Cluster_FullCallback_( const std::string& aInFile , const double& aR, const double& aT, const boost::python::object& aCallback )
{
  throw std::runtime_error( "RoIproxy data-type not yet adapted for python" );
  AutoRoi_Cluster_FullCallback( aInFile , aR , aT , [&]( RoIproxy& aRoIproxy ){ aCallback( aRoIproxy ); } );
}

//! Wrapper to automatically extract RoI, clusterize and apply a simple python callback
//! \param aInFile     The name of the localization file
//! \param aR          The R value of the clusterizer
//! \param aT          The T value of the clusterizer
//! \param aCallback   The simple python callback to be applied
__attribute__((flatten))
void _AutoRoi_Cluster_SimpleCallback_( const std::string& aInFile , const double& aR, const double& aT, const boost::python::object& aCallback )
{
  AutoRoi_Cluster_SimpleCallback( aInFile , aR , aT , [&]( const std::vector< ClusterWrapper >& aClusters ){ aCallback( aClusters ); } );
}



//! Wrapper to manually extract RoI, run scan and apply a full python callback
//! \param aInFile     The name of the localization file
//! \param aManualRoI  The manually-specified RoI window
//! \param aScanConfig The configuration for the scan
//! \param aCallback   The full python callback to be applied
__attribute__((flatten))
void _ManualRoi_Scan_FullCallback_( const std::string& aInFile , const ManualRoI& aManualRoI , const ScanConfiguration& aScanConfig , const boost::python::object& aCallback )
{
  ManualRoi_Scan_FullCallback( aInFile , aManualRoI , aScanConfig, [&]( RoIproxy& aRoIproxy, const double& aR , const double& aT ){ aCallback( aRoIproxy , aR , aT ); } );
}

//! Wrapper to manually extract RoI, run scan and apply a simple python callback
//! \param aInFile     The name of the localization file
//! \param aManualRoI  The manually-specified RoI window
//! \param aScanConfig The configuration for the scan
//! \param aCallback   The simple python callback to be applied
__attribute__((flatten))
void _ManualRoi_Scan_SimpleCallback_( const std::string& aInFile , const ManualRoI& aManualRoI , const ScanConfiguration& aScanConfig , const boost::python::object& aCallback )
{
  ManualRoi_Scan_SimpleCallback( aInFile , aManualRoI , aScanConfig , [&]( const std::vector< ScanEntry >& aScanResults ){ aCallback( aScanResults ); } );
}



//! Wrapper to manually extract RoI, clusterize and apply a full python callback
//! \param aInFile     The name of the localization file
//! \param aManualRoI  The manually-specified RoI window
//! \param aR          The R value of the clusterizer
//! \param aT          The T value of the clusterizer
//! \param aCallback   The full python callback to be applied
__attribute__((flatten))
void _ManualRoi_Cluster_FullCallback_( const std::string& aInFile , const ManualRoI& aManualRoI , const double& aR, const double& aT, const boost::python::object& aCallback )
{
  throw std::runtime_error( "RoIproxy data-type not yet adapted for python" );
  ManualRoi_Cluster_FullCallback( aInFile , aManualRoI , aR , aT , [&]( RoIproxy& aRoIproxy ){ aCallback( aRoIproxy ); } );
}

//! Wrapper to manually extract RoI, clusterize and apply a simple python callback
//! \param aInFile     The name of the localization file
//! \param aManualRoI  The manually-specified RoI window
//! \param aR          The R value of the clusterizer
//! \param aT          The T value of the clusterizer
//! \param aCallback   The simple python callback to be applied
__attribute__((flatten))
void _ManualRoi_Cluster_SimpleCallback_( const std::string& aInFile , const ManualRoI& aManualRoI , const double& aR, const double& aT, const boost::python::object& aCallback )
{
  ManualRoi_Cluster_SimpleCallback( aInFile , aManualRoI , aR , aT , [&]( const std::vector< ClusterWrapper >& aClusters ){ aCallback( aClusters ); } );
}


// namespace std
// {
//   template<> // Specializing class
//   template<> // Specializing constructor
//   function< void( const std::vector< ClusterWrapper >& ) >::function< boost::python::object& >( boost::python::object& aCallback ) :
//     function( [&]( const std::vector< ClusterWrapper >& aClusters ){ aCallback( aClusters ); } ) // deferred constructor
//   {}
// }


//! Helper Macro to simplify defining functions
//! \param X The function being defined
//! \param ... List of strings giving argument names
#define          FN( X , ... ) def( #X , &X    , args( __VA_ARGS__ ) )

//! Helper Macro to simplify defining functions with python callbacks
//! \param X The function being defined
//! \param ... List of strings giving argument names
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

  class_< ClusterWrapper >( "ClusterWrapper" )  
    .def_readwrite( "localizations" , &ClusterWrapper::localizations )
    .def_readwrite( "area"          , &ClusterWrapper::area )
    .def_readwrite( "perimeter"     , &ClusterWrapper::perimeter )
    .def_readwrite( "centroid_x"    , &ClusterWrapper::centroid_x )
    .def_readwrite( "centroid_y"    , &ClusterWrapper::centroid_y );

  class_< std::vector<ClusterWrapper> >( "VectorClusterWrapper" )
    .def( vector_indexing_suite< std::vector<ClusterWrapper> >() );


  CALLBACK_FN( AutoRoi_Scan_FullCallback   , "aInFile" , "aScanConfig" , "aCallback" );
  CALLBACK_FN( AutoRoi_Scan_SimpleCallback , "aInFile" , "aScanConfig" , "aCallback" );
           FN( AutoRoi_Scan_ToJson         , "aInFile" , "aScanConfig" , "aOutFile" );

  // CALLBACK_FN( AutoRoi_Cluster_FullCallback   , "aInFile" , "aR" , "aT" , "aCallback" );
  CALLBACK_FN( AutoRoi_Cluster_SimpleCallback , "aInFile" , "aR" , "aT" , "aCallback" );

  CALLBACK_FN( ManualRoi_Scan_FullCallback   , "aInFile" , "aManualRoI" , "aScanConfig" , "aCallback" );
  CALLBACK_FN( ManualRoi_Scan_SimpleCallback , "aInFile" , "aManualRoI" , "aScanConfig" , "aCallback" );
           FN( ManualRoi_Scan_ToJson         , "aInFile" , "aManualRoI" , "aScanConfig" , "aOutFile" );

  // CALLBACK_FN( ManualRoi_Cluster_FullCallback   , "aInFile" , "aManualRoI" , "aR" , "aT" , "aCallback" );
  CALLBACK_FN( ManualRoi_Cluster_SimpleCallback , "aInFile" , "aManualRoI" , "aR" , "aT" , "aCallback" );

}

#undef CALLBACK_FN
#undef FN


// //  class_< RoI, boost::noncopyable >( "RoI", init<const RoI&>() )
// //     .def( "__iter__" , &RoI_GetIterator )
// //     .def( "__len__" , &RoI_GetSize )
// //     .def( "Preprocess" , &RoI::Preprocess )
// //     ;

// //   class_< RoIproxy, boost::noncopyable >( "RoIproxy", init< RoI& >() )
// //     ;

// //   class_< Cluster >( "Cluster" )
// //     .def( "__iter__" , &Cluster_GetIterator )
// //     .def( "__len__" , &Cluster_GetSize )
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

