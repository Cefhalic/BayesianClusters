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
#include "BayesianClustering/RoI.hpp"
#include "BayesianClustering/RoIproxy.hpp"
// #include "BayesianClustering/Cluster.hpp"
#include "BayesianClustering/Data.hpp"
// #include "BayesianClustering/DataProxy.hpp"
#include "BayesianClustering/ImageJ_RoI.hpp"


// =====================================================================================================================
//! Helper macro to define a constructor for a std::function of a given signature from a boost::python callback object
//! Using explicit specialization of the class to the specified type and explicit specialization of the constructor to take a object (hence two template<>'s at the start),
//! create a generic variadic lambda that captures the boost::python callback object and pass the lambda to a deferred constructor. Gnarly!
//! \param SIGNATURE The signature of the std::function we are creating a constructor for
#define ADAPT_CALLBACK_CONSTRUCTOR( SIGNATURE ) template<> template<> SIGNATURE::function< const object& >( const object& aCallback ) : function( [&]( auto&&... aArgs ){ aCallback( aArgs... ); } ) {}

ADAPT_CALLBACK_CONSTRUCTOR( tFullScanCallback ); //!< Define a std::function constructor for full scan callback
ADAPT_CALLBACK_CONSTRUCTOR( tSimpleScanCallback );         //!< Define a std::function constructor for simple scan callback
ADAPT_CALLBACK_CONSTRUCTOR( tFullClusterCallback );                               //!< Define a std::function constructor for full clusterizer callback
ADAPT_CALLBACK_CONSTRUCTOR( tSimpleClusterCallback );    //!< Define a std::function constructor for simple clusterizer callback

#undef ADAPT_CALLBACK_CONSTRUCTOR


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
namespace Adapted
{
  // auto AutoRoi_Scan_FullCallback      = []( const std::string& aInFile , const ScanConfiguration& aScanConfig , const object& aCallback ){ ::AutoRoi_Scan_FullCallback(      aInFile , aScanConfig , aCallback ); }; //!< Lambda to automatically extract RoI, run scan and apply a full python callback
  auto AutoRoi_Scan_SimpleCallback       = []( const std::string& aInFile , const ScanConfiguration& aScanConfig , const object& aCallback ){ ::AutoRoi_Scan_SimpleCallback(    aInFile , aScanConfig , aCallback ); }; //!< Lambda to automatically extract RoI, run scan and apply a simple python callback
  // auto AutoRoi_Cluster_FullCallback   = []( const std::string& aInFile , const double& aR , const double& aT ,  const object& aCallback ){ ::AutoRoi_Cluster_FullCallback(   aInFile , aR , aT ,     aCallback ); }; //!< Lambda to automatically extract RoI, clusterize and apply a full python callback
  auto AutoRoi_Cluster_SimpleCallback    = []( const std::string& aInFile , const double& aR , const double& aT ,  const object& aCallback ){ ::AutoRoi_Cluster_SimpleCallback( aInFile , aR , aT ,     aCallback ); }; //!< Lambda to automatically extract RoI, clusterize and apply a simple python callback

  // auto ManualRoi_Scan_FullCallback    = []( const std::string& aInFile , const ManualRoI& aManualRoI , const ScanConfiguration& aScanConfig , const object& aCallback ){ ::ManualRoi_Scan_FullCallback(      aInFile , aManualRoI , aScanConfig , aCallback ); }; //!< Lambda to manually specify RoI, run scan and apply a full python callback
  auto ManualRoi_Scan_SimpleCallback     = []( const std::string& aInFile , const ManualRoI& aManualRoI , const ScanConfiguration& aScanConfig , const object& aCallback ){ ::ManualRoi_Scan_SimpleCallback(    aInFile , aManualRoI , aScanConfig , aCallback ); }; //!< Lambda to manually specify RoI, run scan and apply a simple python callback
  // auto ManualRoi_Cluster_FullCallback = []( const std::string& aInFile , const ManualRoI& aManualRoI , const double& aR , const double& aT ,  const object& aCallback ){ ::ManualRoi_Cluster_FullCallback(   aInFile , aManualRoI , aR , aT ,     aCallback ); }; //!< Lambda to manually specify RoI, clusterize and apply a full python callback
  auto ManualRoi_Cluster_SimpleCallback  = []( const std::string& aInFile , const ManualRoI& aManualRoI , const double& aR , const double& aT ,  const object& aCallback ){ ::ManualRoi_Cluster_SimpleCallback( aInFile , aManualRoI , aR , aT ,     aCallback ); }; //!< Lambda to manually specify RoI, clusterize and apply a simple python callback

  // auto ImageJRoi_Scan_FullCallback    = []( const std::string& aInFile , const std::string& aImageJ , const double& aScale , const ScanConfiguration& aScanConfig , const object& aCallback ){ ::ImageJRoi_Scan_FullCallback(      aInFile , aImageJ , aScale , aScanConfig , aCallback ); }; //!< Lambda to extract RoI via an ImagJ RoI file, run scan and apply a full python callback
  auto ImageJRoi_Scan_SimpleCallback     = []( const std::string& aInFile , const std::string& aImageJ , const double& aScale , const ScanConfiguration& aScanConfig , const object& aCallback ){ ::ImageJRoi_Scan_SimpleCallback(    aInFile , aImageJ , aScale , aScanConfig , aCallback ); }; //!< Lambda to extract RoI via an ImagJ RoI file, run scan and apply a simple python callback
  // auto ImageJRoi_Cluster_FullCallback = []( const std::string& aInFile , const std::string& aImageJ , const double& aScale , const double& aR , const double& aT ,  const object& aCallback ){ ::ImageJRoi_Cluster_FullCallback(   aInFile , aImageJ , aScale , aR , aT ,     aCallback ); }; //!< Lambda to extract RoI via an ImagJ RoI file, clusterize and apply a full python callback
  auto ImageJRoi_Cluster_SimpleCallback  = []( const std::string& aInFile , const std::string& aImageJ , const double& aScale , const double& aR , const double& aT ,  const object& aCallback ){ ::ImageJRoi_Cluster_SimpleCallback( aInFile , aImageJ , aScale , aR , aT ,     aCallback ); }; //!< Lambda to extract RoI via an ImagJ RoI file, clusterize and apply a simple python callback
} 
// =====================================================================================================================


// =====================================================================================================================
//! Debugging tool to get the raw x-y coordinates
//! \param aFile The name of the localizations file
//! \return The x-coordinates and the y-coordinates of the raw points as a python tuple (optimised for displaying in MatPlotLib)
boost::python::tuple GetLocalizations( const std::string& aFile )
{
  LocalizationFile lFile( aFile );
  boost::python::list x , y;
  for( auto& i : lFile.data() )
  {
    x.append( i.x );
    y.append( i.y );
  }
  return boost::python::make_tuple( x , y );
}

//! Debugging tool to get the raw coordinates from an ImageJ RoI file  
//! \param aFile The name of an ImageJ RoI file
//! \return A list of python tuples, each containing the x-coordinates and the y-coordinates of the polygon points (optimised for displaying in MatPlotLib)
boost::python::list GetRoIs( const std::string& aFile )
{
  boost::python::list lRet;
  std::map< std::string , roi_polygon > lRoIs = OpenRoiZipfile( aFile );
  for( const auto& i : lRoIs )
  {
    boost::python::list x , y;
    for( const auto& point : i.second )
    {
      x.append( boost::geometry::get<0>(point) );
      y.append( boost::geometry::get<1>(point) );
    }
    lRet.append( boost::python::make_tuple( x , y ) );
  }

  return lRet;
}

//! Debugging tool to get the raw x-y coordinates and which RoI they are included in
//! \param aFile The name of the localizations file
//! \param aRoIFile  The name of an ImageJ RoI file file
//! \param aScale The size of the LSB in the ImageJ file
//! \return A python tuple of the raw localizations and a list of tuples containing the x-coordinates and the y-coordinates of the localizations in each RoI (both optimised for displaying in MatPlotLib) 
boost::python::tuple CheckRoIs( const std::string& aFile , const std::string& aRoIFile , const double& aScale )
{
  LocalizationFile lFile( aFile );
  boost::python::list x , y;
  for( auto& i : lFile.data() )
  {
    x.append( i.x );
    y.append( i.y );
  }
  auto lLocs = boost::python::make_tuple( x , y );

  boost::python::list lRoIs;
  lFile.ExtractRoIs( aRoIFile , aScale , [&]( RoI& aRoI ) {
                                                            boost::python::list x , y;
                                                            for( auto& i : aRoI.data() )
                                                            {
                                                              x.append( i.x + aRoI.getCentreX() );
                                                              y.append( i.y + aRoI.getCentreY() );
                                                            }
                                                            lRoIs.append( boost::python::make_tuple( x , y ) );
                                                          } );
  
  return boost::python::make_tuple( lLocs , lRoIs );
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
          FN( AutoRoi_Scan_ToJson          , "Automatically extract RoI, run scan and dump to JSON file"        , "aInFile" , "aScanConfig" , "aOutputPattern" );

  // ADAPTED_FN( AutoRoi_Cluster_FullCallback , "Automatically extract RoI, clusterize and apply a full call-back"   , "aInFile" , "aR" , "aT" , "aCallback" );
  ADAPTED_FN( AutoRoi_Cluster_SimpleCallback  , "Automatically extract RoI, clusterize and apply a simple call-back" , "aInFile" , "aR" , "aT" , "aCallback" );
          FN( AutoRoi_Cluster_ToJson          , "Automatically extract RoI, clusterize and dump to JSON file"        , "aInFile" , "aR" , "aT" , "aOutputPattern" );

  // ADAPTED_FN( ManualRoi_Scan_FullCallback , "Manually specify RoI, run scan and apply a full call-back"   , "aInFile" , "aManualRoI" , "aScanConfig" , "aCallback" );
  ADAPTED_FN( ManualRoi_Scan_SimpleCallback  , "Manually specify RoI, run scan and apply a simple call-back" , "aInFile" , "aManualRoI" , "aScanConfig" , "aCallback" );
          FN( ManualRoi_Scan_ToJson          , "Manually specify RoI, run scan and dump to JSON file"        , "aInFile" , "aManualRoI" , "aScanConfig" , "aOutputPattern" );

  // ADAPTED_FN( ManualRoi_Cluster_FullCallback , "Manually specify RoI, clusterize and apply a full call-back"   , "aInFile" , "aManualRoI" , "aR" , "aT" , "aCallback" );
  ADAPTED_FN( ManualRoi_Cluster_SimpleCallback  , "Manually specify RoI, clusterize and apply a simple call-back" , "aInFile" , "aManualRoI" , "aR" , "aT" , "aCallback" );
          FN( ManualRoi_Cluster_ToJson          , "Manually specify RoI, clusterize and dump to JSON file"        , "aInFile" , "aManualRoI" , "aR" , "aT" , "aOutputPattern" );

  // ADAPTED_FN( ImageJRoi_Scan_FullCallback , "Extract RoI via an ImagJ RoI file, run scan and apply a full call-back"   , "aInFile" , "aImageJ" , "aScale" , "aScanConfig" , "aCallback" );
  ADAPTED_FN( ImageJRoi_Scan_SimpleCallback  , "Extract RoI via an ImagJ RoI file, run scan and apply a simple call-back" , "aInFile" , "aImageJ" , "aScale" , "aScanConfig" , "aCallback" );
          FN( ImageJRoi_Scan_ToJson          , "Extract RoI via an ImagJ RoI file, run scan and dump to JSON file"        , "aInFile" , "aImageJ" , "aScale" , "aScanConfig" , "aOutputPattern" );

  // ADAPTED_FN( ImageJRoi_Cluster_FullCallback , "Extract RoI via an ImagJ RoI file, clusterize and apply a full call-back"   , "aInFile" , "aImageJ" , "aScale" , "aR" , "aT" , "aCallback" );
  ADAPTED_FN( ImageJRoi_Cluster_SimpleCallback  , "Extract RoI via an ImagJ RoI file, clusterize and apply a simple call-back" , "aInFile" , "aImageJ" , "aScale" , "aR" , "aT" , "aCallback" );
          FN( ImageJRoi_Cluster_ToJson          , "Extract RoI via an ImagJ RoI file, clusterize and dump to JSON file"        , "aInFile" , "aImageJ" , "aScale" , "aR" , "aT" , "aOutputPattern" );

  // ------------------------------------------


  def( "GetLocalizations" , &GetLocalizations );
  def( "GetRoIs" , &GetRoIs );
  def( "CheckRoIs" , &CheckRoIs );


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

