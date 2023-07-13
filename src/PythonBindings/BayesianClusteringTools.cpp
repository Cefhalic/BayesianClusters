//! \file BayesianClusteringTools.cpp
//! Self-contained sourcefile for producing python-bindings

/* ===== BOOST libraries ===== */
#pragma GCC diagnostic push 
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#include <boost/python.hpp>
#pragma GCC diagnostic pop

using namespace boost::python;

/* ===== Cluster sources ===== */
#include "Utilities/Units.hpp"
#include "BayesianClustering/LocalizationFile.hpp"
#include "BayesianClustering/Data.hpp"
#include "BayesianClustering/RoI.hpp"
#include "BayesianClustering/ImageJ_RoI.hpp"


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
  lFile.ExtractRoIs( ImageJRoI { aRoIFile , aScale } , [&]( RoI& aRoI ) {
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
//! Boost Python Wrapper providing bindings for our C++ functions
BOOST_PYTHON_MODULE( BayesianClusteringTools )
{
  scope().attr( "nanometer" ) = nanometer;
  scope().attr( "micrometer" ) = micrometer;

  def( "GetLocalizations" , &GetLocalizations , args( "aFile" ) );
  def( "GetRoIs" , &GetRoIs , args( "aFile" ) );
  def( "CheckRoIs" , &CheckRoIs , args( "aFile" , "aRoIFile" , "aScale" ) );
}
// =====================================================================================================================
