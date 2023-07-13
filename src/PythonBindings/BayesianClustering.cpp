//! \file BayesianClustering.cpp
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
#include <boost/preprocessor/seq/enum.hpp>
#include <boost/preprocessor/seq/for_each_product.hpp>

/* ===== Cluster sources ===== */
#include "Utilities/Units.hpp"
#include "BayesianClustering/API.hpp"
#include "BayesianClustering/LocalizationFile.hpp"
#include "BayesianClustering/Configuration.hpp"
#include "BayesianClustering/RoI.hpp"
#include "BayesianClustering/RoIproxy.hpp"
// #include "BayesianClustering/Cluster.hpp"
#include "BayesianClustering/Data.hpp"
// #include "BayesianClustering/DataProxy.hpp"


// =====================================================================================================================
//! Factory function to construct a ScanConfiguration which take the parameters directly in python
//! \param aSigmaBins   The number of sigma bins
//! \param aSigmaMin     The lowest sigma bin
//! \param aSigmaMax     The highest sigma bin
//! \param aInterpolator A python function call or python dictionary containing a set of points from which to create an interpolator
//! \param aRbins    The number of R bins to scan over
//! \param aMinScanR The lowest value of R to scan
//! \param aMaxScanR The largest value of R to scan
//! \param aTbins    The number of T bins to scan over
//! \param aMinScanT The lowest value of T to scan
//! \param aMaxScanT The largest value of T to scan
//! \param aPB    The P_b parameter
//! \param aAlpha The alpha parameter
//! \return a shared pointer to the new ScanConfiguration
std::shared_ptr< ScanConfiguration > ScanConfigurationConstructor( 
      const std::size_t& aSigmaBins, const double& aSigmaMin, const double& aSigmaMax, const object& aInterpolator ,
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
    auto lFunc = [&]( const double& aVal ){ return call<double>( aInterpolator.ptr() , aVal ); };
    return std::shared_ptr<ScanConfiguration>( new ScanConfiguration( aSigmaBins, aSigmaMin, aSigmaMax, lFunc, aRbins, aMinScanR, aMaxScanR , aTbins, aMinScanT, aMaxScanT , aPB , aAlpha ) );
  }
}
// =====================================================================================================================



// =====================================================================================================================
//! Helper Macro to deal with the boilerplate when dealing with structs
//! \param r BOOST PP internal
//! \param CLASS The Class name
//! \param ARG One of the arguments
#define STRUCT_ARG( r , CLASS , ARG ) .def_readwrite( BOOST_PP_STRINGIZE( ARG ) , &CLASS::ARG )

//! Helper Macro to deal with the boilerplate when dealing with structs
//! \param CLASS The Class name
//! \param DOC A string to be used as python documentation
//! \param ARGS Sequence of arguments
#define EXPOSE_STRUCT_NO_CONSTRUCTOR( CLASS , DOC , ARGS ) class_< CLASS >( #CLASS , DOC ) BOOST_PP_SEQ_FOR_EACH( STRUCT_ARG , CLASS , ARGS )

//! Helper Macro to deal with the boilerplate when dealing with structs
//! \param CLASS The Class name
//! \param DOC A string to be used as python documentation
//! \param CONSTRUCTOR_ARGS Constructor argument types
//! \param ARGS Sequence of arguments
#define EXPOSE_STRUCT( CLASS , DOC , CONSTRUCTOR_ARGS , ARGS ) class_< CLASS >( #CLASS , DOC , init< BOOST_PP_SEQ_ENUM( CONSTRUCTOR_ARGS ) >() ) BOOST_PP_SEQ_FOR_EACH( STRUCT_ARG , CLASS , ARGS )

//! Helper Macro to deal with the boilerplate when dealing with vectors of objects
//! \param CLASS The Class name
#define EXPOSE_VECTOR( CLASS ) class_< std::vector< CLASS > >( BOOST_PP_STRINGIZE( Vector##CLASS ) , BOOST_PP_STRINGIZE( An STL vector of CLASS ) ).def( vector_indexing_suite< std::vector< CLASS > >() );

//! The ROI configs available 
#define ROICONFIGS       (ImageJRoI)(AutoRoI)(ManualRoI) 

//! The scan callbacks available 
#define SCANCALLBACKS    (std::string)(tSimpleScanCallback)(tFullScanCallback)

//! The clustering callbacks available 
#define CLUSTERCALLBACKS (std::string)(tSimpleClusterCallback)(tFullClusterCallback)

//! Macro to produce all permutations of RunScan
//! \param r UNUSED
//! \param product Sequence of each possible products
#define RUNSCAN( r , product ) def( "RunScan" , static_cast<void (*)( const std::string& , const BOOST_PP_SEQ_ELEM(0,product)& , const ScanConfiguration& , const BOOST_PP_SEQ_ELEM(1,product)& aCallback )>(&RunScan) , args( "aInFile" , "aRoIConfig" , "aScanConfig" , "aHandler" ) ); 

//! Macro to produce all permutations of RunScan
//! \param r UNUSED
//! \param product Sequence of each possible products
#define RUNCLUSTER( r , product ) def( "RunClustering" , static_cast<void (*)( const std::string& , const BOOST_PP_SEQ_ELEM(0,product)& , const double& , const double& , const BOOST_PP_SEQ_ELEM(1,product)& aCallback )>(&RunClustering) , args( "aInFile" , "aRoIConfig" , "aR" , "aT" , "aHandler" ) ); 



//! Boost Python Wrapper providing bindings for our C++ functions
BOOST_PYTHON_MODULE( BayesianClustering )
{
  scope().attr( "nanometer" ) = nanometer;
  scope().attr( "micrometer" ) = micrometer;

  EXPOSE_STRUCT( ManualRoI , "A struct for storing the parameters of a manual RoI" , (const double&)(const double&)(const double&)(const double&) , (x)(y)(width)(height) );
  EXPOSE_STRUCT_NO_CONSTRUCTOR( AutoRoI , "A struct for storing the parameters for Auto RoI extraction" , );
  EXPOSE_STRUCT( ImageJRoI , "A struct for storing the parameters of an ImageJ RoI file" , (const std::string&)(const double&) , (filename)(scale) );

  EXPOSE_STRUCT_NO_CONSTRUCTOR( ScanEntry , "A struct for storing a result of an individual scan configuration" , (r)(t)(score) );
  EXPOSE_STRUCT_NO_CONSTRUCTOR( ClusterWrapper , "A struct for storing extracted parameters from a cluster" , (localizations)(area)(perimeter)(centroid_x)(centroid_y) );

  EXPOSE_VECTOR( ScanEntry );
  EXPOSE_VECTOR( ClusterWrapper );

  class_< ScanConfiguration, std::shared_ptr<ScanConfiguration>, boost::noncopyable >( "ScanConfiguration" , "A class for storing the scan configuration parameters" , no_init )
    .def( init< const std::string& >( arg( "aCfgFile" ) ) )
    .def( "__init__" , make_constructor( &ScanConfigurationConstructor, default_call_policies(), args( "aSigmaBins", "aSigmaMin", "aSigmaMax", "aInterpolator", "aRbins", "aMinScanR", "aMaxScanR", "aTbins", "aMinScanT", "aMaxScanT", "aPB", "aAlpha" ) ) ); 

  BOOST_PP_SEQ_FOR_EACH_PRODUCT( RUNSCAN , (ROICONFIGS)(SCANCALLBACKS) );
  BOOST_PP_SEQ_FOR_EACH_PRODUCT( RUNCLUSTER , (ROICONFIGS)(CLUSTERCALLBACKS) );
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

