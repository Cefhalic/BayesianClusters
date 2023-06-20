//! \file PythonBindings.cpp
//! Self-contained sourcefile for producing python-bindings

/* ===== C++ libraries ===== */
#include <iostream>
// #include <functional>

/* ===== BOOST libraries ===== */
#pragma GCC diagnostic push 
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#pragma GCC diagnostic pop

using namespace boost::python;

#include <boost/preprocessor/stringize.hpp>
#include <boost/preprocessor/seq/for_each.hpp>


/* ===== Cluster sources ===== */
#include "Utilities/Units.hpp"
#include "BayesianClustering/API.hpp"
#include "BayesianClustering/LocalizationFile.hpp"
#include "BayesianClustering/Configuration.hpp"
// #include "BayesianClustering/RoI.hpp"
#include "BayesianClustering/RoIproxy.hpp"
// #include "BayesianClustering/Cluster.hpp"
// #include "BayesianClustering/RoI.hpp"
#include "BayesianClustering/Data.hpp"
// #include "BayesianClustering/DataProxy.hpp"


// =====================================================================================================================
//! Helper macro to define a constructor for a std::function of a given signature from a boost::python callback object
//! Using explicit specialization of the class to the specified type and explicit specialization of the constructor to take a boost::python::object (hence two template<>'s at the start),
//! create a generic variadic lambda that captures the boost::python callback object and pass the lambda to a deferred constructor. Gnarly!
//! \param SIGNATURE The signature of the std::function we are creating a constructor for
#define ADAPT_CALLBACK_CONSTRUCTOR( SIGNATURE ) template<> template<> SIGNATURE::function< boost::python::object& >( boost::python::object& aCallback ) : function( [&]( auto&&... aArgs ){ aCallback( aArgs... ); } ) {}

ADAPT_CALLBACK_CONSTRUCTOR( std::function< void( RoIproxy&, const double&, const double& ) > ); //!< Define a std::function constructor for full scan callback
ADAPT_CALLBACK_CONSTRUCTOR( std::function< void( const std::vector< ScanEntry >& ) > );         //!< Define a std::function constructor for simple scan callback
ADAPT_CALLBACK_CONSTRUCTOR( std::function< void( RoIproxy& ) > );                               //!< Define a std::function constructor for full clusterizer callback
ADAPT_CALLBACK_CONSTRUCTOR( std::function< void( const std::vector< ClusterWrapper >& ) > );    //!< Define a std::function constructor for simple clusterizer callback

#undef ADAPT_CALLBACK_CONSTRUCTOR

std::shared_ptr< ScanConfiguration > ScanConfigurationConstructor( 
      const std::size_t& aSigmaBins, const double& aSigmaMin, const double& aSigmaMax, const boost::python::object& aInterpolator ,
      const std::size_t& aRbins, const double& aMinScanR, const double& aMaxScanR ,
      const std::size_t& aTbins, const double& aMinScanT, const double& aMaxScanT ,
      const double& aPB , const double& aAlpha )
{
  extract< dict > lInterpolator( aInterpolator );
  if( lInterpolator.check() )
  { 
    dict lDict = lInterpolator;
    std::map< double , double > lMap;
    for (auto it = stl_input_iterator<tuple>( lDict.items() ); it != stl_input_iterator<tuple>(); ++it) lMap[ extract<double>((*it)[0]) ] = extract<double>((*it)[1]);
    return std::shared_ptr<ScanConfiguration>( new ScanConfiguration( aSigmaBins, aSigmaMin, aSigmaMax, lMap, aRbins, aMinScanR, aMaxScanR , aTbins, aMinScanT, aMaxScanT , aPB , aAlpha ) );
  }
  else
  {
    auto lFunc = [&]( const double& aVal ){ return boost::python::call<double>( aInterpolator.ptr() , aVal ); };
    return std::shared_ptr<ScanConfiguration>( new ScanConfiguration( aSigmaBins, aSigmaMin, aSigmaMax, lFunc, aRbins, aMinScanR, aMaxScanR , aTbins, aMinScanT, aMaxScanT , aPB , aAlpha ) );
  }
}
// =====================================================================================================================


// =====================================================================================================================
namespace Adapted
{
  // auto AutoRoi_Scan_FullCallback        = []( const std::string& aInFile , const ScanConfiguration& aScanConfig , const boost::python::object& aCallback ){ ::AutoRoi_Scan_FullCallback(      aInFile , aScanConfig , aCallback ); }; //!< Lambda to automatically extract RoI, run scan and apply a full python callback
  auto AutoRoi_Scan_SimpleCallback      = []( const std::string& aInFile , const ScanConfiguration& aScanConfig , const boost::python::object& aCallback ){ ::AutoRoi_Scan_SimpleCallback(    aInFile , aScanConfig , aCallback ); }; //!< Lambda to automatically extract RoI, run scan and apply a simple python callback
  // auto AutoRoi_Cluster_FullCallback     = []( const std::string& aInFile , const double& aR , const double& aT ,  const boost::python::object& aCallback ){ ::AutoRoi_Cluster_FullCallback(   aInFile , aR , aT ,     aCallback ); }; //!< Lambda to automatically extract RoI, clusterize and apply a full python callback
  auto AutoRoi_Cluster_SimpleCallback   = []( const std::string& aInFile , const double& aR , const double& aT ,  const boost::python::object& aCallback ){ ::AutoRoi_Cluster_SimpleCallback( aInFile , aR , aT ,     aCallback ); }; //!< Lambda to automatically extract RoI, clusterize and apply a simple python callback
  // auto ManualRoi_Scan_FullCallback      = []( const std::string& aInFile , const ManualRoI& aManualRoI , const ScanConfiguration& aScanConfig , const boost::python::object& aCallback ){ ::ManualRoi_Scan_FullCallback(      aInFile , aManualRoI , aScanConfig , aCallback ); }; //!< Lambda to manually extract RoI, run scan and apply a full python callback
  auto ManualRoi_Scan_SimpleCallback    = []( const std::string& aInFile , const ManualRoI& aManualRoI , const ScanConfiguration& aScanConfig , const boost::python::object& aCallback ){ ::ManualRoi_Scan_SimpleCallback(    aInFile , aManualRoI , aScanConfig , aCallback ); }; //!< Lambda to manually extract RoI, run scan and apply a simple python callback
  // auto ManualRoi_Cluster_FullCallback   = []( const std::string& aInFile , const ManualRoI& aManualRoI , const double& aR , const double& aT ,  const boost::python::object& aCallback ){ ::ManualRoi_Cluster_FullCallback(   aInFile , aManualRoI , aR , aT ,     aCallback ); }; //!< Lambda to manually extract RoI, clusterize and apply a full python callback
  auto ManualRoi_Cluster_SimpleCallback = []( const std::string& aInFile , const ManualRoI& aManualRoI , const double& aR , const double& aT ,  const boost::python::object& aCallback ){ ::ManualRoi_Cluster_SimpleCallback( aInFile , aManualRoI , aR , aT ,     aCallback ); }; //!< Lambda to manually extract RoI, clusterize and apply a simple python callback
} 
// =====================================================================================================================


// =====================================================================================================================
//! Helper Macro to simplify defining functions
//! \param X The function being defined
//! \param DOC A string to be used as python documentation
//! \param ... List of strings giving argument names
#define FN( X , DOC , ... ) def( #X , &X    , args( __VA_ARGS__ ) , DOC )

//! Helper Macro to simplify defining functions with python callbacks
//! \param X The function being defined
//! \param DOC A string to be used as python documentation
//! \param ... List of strings giving argument names
#define ADAPTED_FN( X , DOC , ... ) def( #X , +Adapted::X , args( __VA_ARGS__ ) , DOC )

//! Helper Macro to deal with the boilerplate when dealing with structs
//! \param r BOOST PP internal
//! \param CLASS The Class name
//! \param ARG One of the arguments
#define STRUCT_ARG( r , CLASS , ARG ) .def_readwrite( BOOST_PP_STRINGIZE( ARG ) , &CLASS::ARG )

//! Helper Macro to deal with the boilerplate when dealing with structs
//! \param CLASS The Class name
//! \param DOC A string to be used as python documentation
//! \param ARGS One of the arguments
#define EXPOSE_STRUCT( CLASS , DOC , ARGS ) class_< CLASS >( #CLASS , DOC ) BOOST_PP_SEQ_FOR_EACH( STRUCT_ARG , CLASS , ARGS )

//! Helper Macro to deal with the boilerplate when dealing with vectors of objects
//! \param CLASS The Class name
#define EXPOSE_VECTOR( CLASS ) class_< std::vector< CLASS > >( BOOST_PP_STRINGIZE( Vector##CLASS ) , BOOST_PP_STRINGIZE( An STL vector of CLASS ) ).def( vector_indexing_suite< std::vector< CLASS > >() );


//! Boost Python Wrapper providing bindings for our C++ functions
BOOST_PYTHON_MODULE( BayesianClustering )
{
  scope().attr( "nanometer" ) = nanometer;
  scope().attr( "micrometer" ) = micrometer;

  EXPOSE_STRUCT( ManualRoI , "A struct for storing the parameters of a manual RoI" , (x)(y)(width)(height) );
  EXPOSE_STRUCT( ScanEntry , "A struct for storing a result of an individual scan configuration" , (r)(t)(score) );
  EXPOSE_STRUCT( ClusterWrapper , "A struct for storing extracted parameters from a cluster" , (localizations)(area)(perimeter)(centroid_x)(centroid_y) );

  EXPOSE_VECTOR( ScanEntry );
  EXPOSE_VECTOR( ClusterWrapper );

  class_< ScanConfiguration, std::shared_ptr<ScanConfiguration>, boost::noncopyable >( "ScanConfiguration" , "A class for storing the scan configuration parameters" , no_init )
    .def( init< const std::string& >( arg( "aCfgFile" ) ) )
    .def( "__init__" , make_constructor( &ScanConfigurationConstructor, default_call_policies(), args( "aSigmaBins", "aSigmaMin", "aSigmaMax", "aInterpolator", "aRbins", "aMinScanR", "aMaxScanR", "aTbins", "aMinScanT", "aMaxScanT", "aPB", "aAlpha" ) ) ); 

  // ADAPTED_FN( AutoRoi_Scan_FullCallback , "Automatically extract RoI, run scan and apply a full call-back"   , "aInFile" , "aScanConfig" , "aCallback" );
  ADAPTED_FN( AutoRoi_Scan_SimpleCallback  , "Automatically extract RoI, run scan and apply a simple call-back" , "aInFile" , "aScanConfig" , "aCallback" );
          FN( AutoRoi_Scan_ToJson          , "Automatically extract RoI, run scan and dump to JSON file"        , "aInFile" , "aScanConfig" , "aOutFile" );

  // ADAPTED_FN( AutoRoi_Cluster_FullCallback , "Automatically extract RoI, clusterize and apply a full call-back"   , "aInFile" , "aR" , "aT" , "aCallback" );
  ADAPTED_FN( AutoRoi_Cluster_SimpleCallback  , "Automatically extract RoI, clusterize and apply a simple call-back" , "aInFile" , "aR" , "aT" , "aCallback" );

  // ADAPTED_FN( ManualRoi_Scan_FullCallback , "Manually specify RoI, run scan and apply a full call-back"   , "aInFile" , "aManualRoI" , "aScanConfig" , "aCallback" );
  ADAPTED_FN( ManualRoi_Scan_SimpleCallback  , "Manually specify RoI, run scan and apply a simple call-back" , "aInFile" , "aManualRoI" , "aScanConfig" , "aCallback" );
          FN( ManualRoi_Scan_ToJson          , "Manually specify RoI, run scan and dump to JSON file"        , "aInFile" , "aManualRoI" , "aScanConfig" , "aOutFile" );

  // ADAPTED_FN( ManualRoi_Cluster_FullCallback , "Manually specify RoI, clusterize and apply a full call-back"   , "aInFile" , "aManualRoI" , "aR" , "aT" , "aCallback" );
  ADAPTED_FN( ManualRoi_Cluster_SimpleCallback  , "Manually specify RoI, clusterize and apply a simple call-back" , "aInFile" , "aManualRoI" , "aR" , "aT" , "aCallback" );
}

#undef ADAPTED_FN
#undef FN
// =====================================================================================================================


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

