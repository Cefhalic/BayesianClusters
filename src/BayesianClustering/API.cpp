//! \file API.cpp

/* ===== Cluster sources ===== */
#include "BayesianClustering/API.hpp"
#include "BayesianClustering/LocalizationFile.hpp"
#include "BayesianClustering/Configuration.hpp"
#include "BayesianClustering/RoI.hpp"
#include "BayesianClustering/RoIproxy.hpp"
#include "Utilities/Units.hpp"


/* ===== C++ ===== */
#include <map>
#include <algorithm>
#include <mutex>

/* ===== BOOST C++ ===== */
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/filesystem.hpp>

/* ===== FMT library ===== */
#include <fmt/format.h>

// 
std::vector<std::string> strsplit(const std::string& str, const std::string& delim)
{
    std::vector<std::string> result;
    std::size_t start = 0;

    for (std::size_t found = str.find(delim); found != std::string::npos; found = str.find(delim, start))
    {
        result.emplace_back(str.begin() + start, str.begin() + found);
        start = found + delim.size();
    }
    if (start != str.size())
        result.emplace_back(str.begin() + start, str.end());
    return result;
}

//! A callback to dump a scan to a JSON file
//! \param aRoiId       The RoI ID
//! \param aVector   A vector of scan results
//! \param aInFile   The name of the localization file
//! \param aOutputPattern  The name of the output JSON file
void _ScanCallback_Json_( const std::string& aRoiId , const std::vector< ScanEntry >& aVector, const std::string& aInFile , const std::string& aOutputPattern )
{
  using namespace fmt::literals;
  auto lOutFileName = boost::filesystem::path( fmt::format( aOutputPattern , "input"_a = boost::filesystem::path( aInFile ).stem().string() , "roi"_a = aRoiId ) );
  boost::filesystem::create_directories( lOutFileName.parent_path() );

  FILE *fptr = fopen( lOutFileName.c_str() , "w" );
  if (fptr == NULL) throw std::runtime_error("Could not open file");
  fprintf( fptr , "[\n" );
  for( auto& lIt : aVector ) fprintf( fptr , "  { \"r\":%.5e , \"t\":%.5e , \"logP\":%.5e },\n" , lIt.r , lIt.t , lIt.score );
  fseek( fptr, -2, SEEK_CUR ); // Delete the last comma
  fprintf( fptr , "\n]\n" );
  fclose(fptr); 
}


// ScanEntry _BestScore_( const ScanConfiguration& aScanConfig , const std::vector< ScanEntry >& aResults )
// {
//   std::size_t lMax( 0 );
//   for( std::size_t i(1) ; i!= aResults.size() ; ++i )
//     if( aResults[i].score > aResults[lMax].score ) lMax = i;

//   std::size_t R( lMax / aScanConfig.tbounds.bins ), T( lMax % aScanConfig.tbounds.bins );

//   ScanEntry lMean{ 0.0 , 0.0 , 0.0 };
//   for( int r=std::max(0,R-2) ; r!=std::min(aScanConfig.rbounds.bins,R+3) ; ++r ){
//     for( int t=std::max(0,T-2) ; t!=std::min(aScanConfig.tbounds.bins,T+3) ; ++t ){
//       std::size_t lIndex = ( aScanConfig.tbounds.bins * r ) + t;
//       auto& lBin = aResults[lIndex];
//       lMean.r += ( lBin.score * lBin.r );
//       lMean.t += ( lBin.score * lBin.t );
//       lMean.score += lBin.score;
//     }
//   } 

//   lMean.r /= lMean.score;
//   lMean.t /= lMean.score;

//   return lMean;
// }


//! A callback to neatly package the scan results for easy consumption
//! \param aRoI        The region of interest
//! \param aScanConfig The configuration for the scan
//! \param aCallback   The simple callback to be applied
void _FullScanToSimpleScan_( RoI& aRoI , const ScanConfiguration& aScanConfig , const tSimpleScanCallback& aCallback )
{
  std::mutex lMtx;
  std::vector< ScanEntry > lResults;
  aRoI.ScanRT( aScanConfig, [&]( const RoIproxy& aRoI, const double& aR, const double& aT ) { lMtx.lock(); lResults.push_back( { aR, aT, aRoI.mLogP } ); lMtx.unlock(); } );
  std::sort( lResults.begin(), lResults.end() );

  aCallback( aRoI.id() , lResults );  
}

//! A callback to dump a clustering run to a JSON file
//! \param aRoiId    The RoI ID
//! \param aVector   A vector of cluster-wrappers
//! \param aInFile   The name of the localization file
//! \param aOutputPattern  The name of the output JSON file
void _ClusterCallback_Json_( const std::string& aRoiId , const std::vector< ClusterWrapper >& aVector, const std::string& aInFile , const std::string& aOutputPattern )
{
  using namespace fmt::literals;
  auto lOutFileName = boost::filesystem::path( fmt::format( aOutputPattern , "input"_a = boost::filesystem::path( aInFile ).stem().string() , "roi"_a = aRoiId ) );
  boost::filesystem::create_directories( lOutFileName.parent_path() );

  FILE *fptr = fopen( lOutFileName.c_str() , "w" );
  if (fptr == NULL) throw std::runtime_error("Could not open file");
  fprintf( fptr , "[\n" );
  for( auto& lIt : aVector ) fprintf( fptr , "  { \"localizations\":%ld , \"area\":%.5Le , \"perimeter\":%.5Le , \"centroid_x\":%.5e , \"centroid_y\":%.5e },\n" , lIt.localizations , lIt.area , lIt.perimeter , lIt.centroid_x , lIt.centroid_y );
  fseek( fptr, -2, SEEK_CUR ); // Delete the last comma
  fprintf( fptr , "\n]\n" );
  fclose(fptr); 
}


//! A callback to neatly package the scan results for easy consumption
//! \param aRoIproxy   The region-proxy containing the clusters
//! \param aCallback   The simple callback to be applied
void _FullClusterToSimpleCluster_( RoIproxy& aRoIproxy , const tSimpleClusterCallback& aCallback )
{
  typedef boost::geometry::model::point<double, 2, boost::geometry::cs::cartesian> geo_point;
  typedef boost::geometry::model::ring<geo_point> geo_polygon;

  std::map< Cluster* , geo_polygon > lMap;
  for( auto& i : aRoIproxy.mData )
  {
    if( ! i.mCluster ) continue;
    boost::geometry::append( lMap[ i.mCluster->GetParent() ] , geo_point( i.mData->x , i.mData->y ) );
  }

  std::vector< ClusterWrapper > lResults;
  for ( auto& i : lMap )
  {
    boost::geometry::correct( i.second );
    geo_polygon lHull;
    boost::geometry::convex_hull( i.second , lHull );   
    geo_point lCentroid ( 0 , 0 );
    boost::geometry::centroid( i.second , lCentroid );    
    lResults.emplace_back( ClusterWrapper{ i.second.size() , boost::geometry::area( lHull ) , boost::geometry::perimeter( lHull ) , boost::geometry::get<0>( lCentroid ) , boost::geometry::get<1>( lCentroid ) } );
  }

  std::sort( lResults.begin(), lResults.end() );
  aCallback( aRoIproxy.mRoI.id() , lResults );  
}


void _RoI_Info_Json_(RoI& aRoI, const std::string& aInFile, const std::string& aOutputPattern)
{
    auto aRoIid = aRoI.id();
    using namespace fmt::literals;
    auto lOutFileName = boost::filesystem::path(fmt::format(aOutputPattern, "input"_a = boost::filesystem::path(aInFile).stem().string(), "roi"_a = aRoIid));
    boost::filesystem::create_directories(lOutFileName.parent_path());

    FILE* fptr = fopen(lOutFileName.c_str(), "w");
    if (fptr == NULL) throw std::runtime_error("Could not open file");
    fprintf(fptr, "[\n");

    long int  localizations = aRoI.data().size();
    long double        area = aRoI.getArea();

    double centroid_x      = aRoI.getCentreX(),
           centroid_y      = aRoI.getCentreY();

    fprintf(fptr, "  { \"localizations\":%ld , \"area\":%.5Le , \"centroid_x\":%.5e , \"centroid_y\":%.5e },\n", localizations, area, centroid_x, centroid_y);
    fseek(fptr, -2, SEEK_CUR); // Delete the last comma
    fprintf(fptr, "\n]\n");
    fclose(fptr);
}

void _FullAnalysis_(RoI& aRoI, const ScanConfiguration& aScanConfig, const std::string& aOutputPattern_Scan, const std::string& aOutputPattern_Cluster, const std::string& aOutputPattern_Info)
{    
    auto RoIid = aRoI.id();    
    // Scan
    std::mutex lMtx;
    std::vector< ScanEntry > lResults;
    aRoI.ScanRT(aScanConfig, [&](const RoIproxy& aRoI, const double& aR, const double& aT) { lMtx.lock(); lResults.push_back({ aR, aT, aRoI.mLogP }); lMtx.unlock(); });
    std::sort(lResults.begin(), lResults.end());    
    // find optimal R,T
    long double max_logP = - INFINITY;
    for( std::size_t i(1) ; i!= lResults.size() ; ++i )
        if( lResults[i].score > max_logP) max_logP = lResults[i].score;
    //
    double sum_t = 0, sum_r = 0, cnt = 0;
    for (std::size_t i(1); i != lResults.size(); ++i)
        if (lResults[i].score == max_logP)
        {
            sum_t += lResults[i].t;
            sum_r += lResults[i].r;
            cnt += 1;
        }
    double  t_star = sum_t/cnt, 
            r_star = sum_r/cnt;
    // Clusterize with optimal R,T and save clusters' properties
    aRoI.Clusterize(r_star, t_star, [RoIid, aOutputPattern_Cluster](RoIproxy& aRoIproxy) { _FullClusterToSimpleCluster_(aRoIproxy, [RoIid,aOutputPattern_Cluster](const std::string& aRoiId, const std::vector< ClusterWrapper >& aVector) { _ClusterCallback_Json_(aRoiId, aVector,"", aOutputPattern_Cluster); }); });
    // save Scan results
    _ScanCallback_Json_(RoIid, lResults, "", aOutputPattern_Scan);
    // save ROI info
    _RoI_Info_Json_(aRoI, "", aOutputPattern_Info);
}

void _SaveClusteredPartitionCSV_(RoIproxy& aRoIproxy, const std::string& aOutputPattern_Cluster)
{
    //USING BOOST STRUCTURES - START
    //const std::string aRoiId = aRoIproxy.mRoI.id();
    //// ROI centre
    //double  CX = aRoIproxy.mRoI.getCentreX(),
    //        CY = aRoIproxy.mRoI.getCentreY();
    //// cluster map 
    //typedef boost::geometry::model::point<double, 2, boost::geometry::cs::cartesian> geo_point;
    //typedef boost::geometry::model::ring<geo_point> geo_polygon;

    //std::map< Cluster*, geo_polygon > cluster_map;
    //for (auto& i : aRoIproxy.mData)
    //{
    //    if (!i.mCluster) continue;
    //    boost::geometry::append(cluster_map[i.mCluster->GetParent()], geo_point(i.mData->x, i.mData->y));
    //}
    //// create csv filename
    //using namespace fmt::literals;
    //auto lOutFileName = boost::filesystem::path(fmt::format(aOutputPattern_Cluster, "roi"_a = aRoiId));
    //boost::filesystem::create_directories(lOutFileName.parent_path());
    //// write clusters with indices to CSV file
    //FILE* fptr = fopen(lOutFileName.c_str(), "w");
    //if (fptr == NULL) throw std::runtime_error("Could not open file");
    //fprintf(fptr, "x,y,index\n");
    ////
    //long int cluster_index = 0;
    //for (auto& i : cluster_map)
    //{
    //    cluster_index++;
    //    boost::geometry::correct(i.second);
    //    for (auto& p : i.second)
    //    {
    //        fprintf(fptr, "%f, %f, %ld\n", (CX + p.get<0>()) / nanometer, (CY + p.get<1>()) / nanometer, cluster_index);
    //    }
    //}
    //fclose(fptr);
    //USING BOOST STRUCTURES - ENDS

    const std::string aRoiId = aRoIproxy.mRoI.id();
    // ROI centre
    double  CX = aRoIproxy.mRoI.getCentreX(), 
            CY = aRoIproxy.mRoI.getCentreY();

    typedef std::vector<std::pair<double,double> >  points;

    std::set<Cluster*> set_of_clusters;    
    // define unique set_of_clusters
    for (auto& i : aRoIproxy.mData)
    {
        if (!i.mCluster) continue;
        //set_of_clusters.insert(i.mCluster->GetParent()); // ??
        set_of_clusters.insert(i.mCluster);
    }
    //
    std::map<Cluster*, points > cluster_map;    
    // initialize "cluster_map" with "set_of_clusters"
    for (auto& i : set_of_clusters)
    {
        points i_th_points = points(0); 
        cluster_map.insert(std::pair<Cluster*, points>(i, i_th_points));
    }    
    // fill cluster_map entries with localizations
    for (auto& i : aRoIproxy.mData)
    {
        cluster_map[i.mCluster].push_back(std::pair<double, double>(CX + i.mData->x, CY + i.mData->y));
    }
    // create csv filename
    using namespace fmt::literals; 
    auto lOutFileName = boost::filesystem::path(fmt::format(aOutputPattern_Cluster, "roi"_a = aRoiId));
    boost::filesystem::create_directories(lOutFileName.parent_path());
    //
    //write clusters with indices to CSV file - KLUDGE
    FILE* fptr = fopen(lOutFileName.c_str(), "w");

    if (fptr == NULL) throw std::runtime_error("Could not open file");
    fprintf(fptr, "x,y,index\n");
    //
    long int cluster_index = 0;
    for (auto& i : cluster_map)
    {
        cluster_index++;
        for (auto& p : i.second)
        {
            fprintf(fptr, "%f, %f, %ld\n", p.first / nanometer, p.second / nanometer, cluster_index);
        }
    }
    fclose(fptr);
}

void _Analyse_FOV_(RoI& aRoI,
    const std::string& aTaskDescription,
    const ScanConfiguration& aScanConfig,
    const double& aR, const double& aT,
    const std::string& aOutputPattern_Scan,
    const std::string& aOutputPattern_Cluster,
    const std::string& aOutputPattern_ClusteredLocalizations,
    const std::string& aOutputPattern_Info)
{
    if (!(0 == aTaskDescription.compare("SCAN_ONLY") || 0 == aTaskDescription.compare("SCAN_THEN_CLUSTER") || 0 == aTaskDescription.compare("CLUSTER_ONLY")))
        return;

    auto RoIid = aRoI.id();
    // will be used for "Cluster"
    double  r_star = aR, t_star = aT;

    if (0 == aTaskDescription.compare("SCAN_ONLY") || 0 == aTaskDescription.compare("SCAN_THEN_CLUSTER"))
    {
        std::mutex lMtx;
        std::vector< ScanEntry > lResults;
        aRoI.ScanRT(aScanConfig, [&](const RoIproxy& aRoI, const double& aR, const double& aT) { lMtx.lock(); lResults.push_back({ aR, aT, aRoI.mLogP }); lMtx.unlock(); });
        std::sort(lResults.begin(), lResults.end());
        // find optimal r,t
        long double max_logP = -INFINITY;
        for (std::size_t i(1); i != lResults.size(); ++i)
            if (lResults[i].score > max_logP) max_logP = lResults[i].score;
        //
        double sum_t = 0, sum_r = 0, cnt = 0;
        for (std::size_t i(1); i != lResults.size(); ++i)
            if (lResults[i].score == max_logP)
            {
                sum_t += lResults[i].t;
                sum_r += lResults[i].r;
                cnt += 1;
            }
        r_star = sum_r / cnt, t_star = sum_t / cnt;
        _ScanCallback_Json_(RoIid, lResults, "", aOutputPattern_Scan);
    };
    if (0 == aTaskDescription.compare("CLUSTER_ONLY") || 0 == aTaskDescription.compare("SCAN_THEN_CLUSTER"))
    {
        const std::string PARAM_STR = aOutputPattern_Cluster + "*" + aOutputPattern_ClusteredLocalizations;
        aRoI.Clusterize(r_star, t_star, [RoIid, PARAM_STR](RoIproxy& aRoIproxy){  
            _FullClusterToSimpleCluster_(aRoIproxy, [RoIid, PARAM_STR](const std::string& aRoiId, const std::vector< ClusterWrapper >& aVector) { _ClusterCallback_Json_(aRoiId, aVector, "", strsplit(PARAM_STR, "*")[0]); }); 
            _SaveClusteredPartitionCSV_(aRoIproxy, strsplit(PARAM_STR, "*")[1]); });
        // save ROI info
        _RoI_Info_Json_(aRoI, "", aOutputPattern_Info);
    };
}


__attribute__((flatten))
void AutoRoi_Scan_FullCallback( const std::string& aInFile , const ScanConfiguration& aScanConfig, const tFullScanCallback& aCallback )
{
  LocalizationFile( aInFile ).ExtractRoIs( [&]( RoI& aRoI ) { aRoI.ScanRT( aScanConfig, aCallback ); } );
}

__attribute__((flatten))
void AutoRoi_Scan_SimpleCallback( const std::string& aInFile , const ScanConfiguration& aScanConfig, const tSimpleScanCallback& aCallback )
{
  LocalizationFile( aInFile ).ExtractRoIs( [&]( RoI& aRoI ) { _FullScanToSimpleScan_( aRoI , aScanConfig, aCallback ); } );
}

__attribute__((flatten))
void AutoRoi_Scan_ToJson( const std::string& aInFile , const ScanConfiguration& aScanConfig, const std::string& aOutputPattern )
{
  AutoRoi_Scan_SimpleCallback( aInFile , aScanConfig , [&]( const std::string& aRoiId , const std::vector< ScanEntry >& aVector ){ _ScanCallback_Json_( aRoiId , aVector , aInFile , aOutputPattern ); } );
}



__attribute__((flatten))
void AutoRoi_Cluster_FullCallback( const std::string& aInFile , const double& aR, const double& aT, const tFullClusterCallback& aCallback )
{  
  LocalizationFile( aInFile ).ExtractRoIs( [&]( RoI& aRoI ) { aRoI.Clusterize( aR, aT, aCallback ); } );
}

__attribute__((flatten))
void AutoRoi_Cluster_SimpleCallback( const std::string& aInFile , const double& aR, const double& aT, const tSimpleClusterCallback& aCallback )
{
  AutoRoi_Cluster_FullCallback( aInFile , aR , aT , [&]( RoIproxy& aRoIproxy ){ _FullClusterToSimpleCluster_( aRoIproxy , aCallback ); } );
}

__attribute__((flatten))
void AutoRoi_Cluster_ToJson( const std::string& aInFile , const double& aR, const double& aT, const std::string& aOutputPattern )
{
  AutoRoi_Cluster_SimpleCallback( aInFile , aR , aT , [&]( const std::string& aRoiId , const std::vector< ClusterWrapper >& aVector ){ _ClusterCallback_Json_( aRoiId , aVector , aInFile , aOutputPattern ); } );
}



__attribute__((flatten))
void ManualRoi_Scan_FullCallback( const std::string& aInFile , const ManualRoI& aManualRoI , const ScanConfiguration& aScanConfig, const tFullScanCallback& aCallback )
{
  LocalizationFile( aInFile ).ExtractRoIs( aManualRoI , [&]( RoI& aRoI ) { aRoI.ScanRT( aScanConfig, aCallback ); } );
}

__attribute__((flatten))
void ManualRoi_Scan_SimpleCallback( const std::string& aInFile , const ManualRoI& aManualRoI , const ScanConfiguration& aScanConfig, const tSimpleScanCallback& aCallback )
{
  LocalizationFile( aInFile ).ExtractRoIs( aManualRoI , [&]( RoI& aRoI ) { _FullScanToSimpleScan_( aRoI , aScanConfig, aCallback ); } );
}

__attribute__((flatten))
void ManualRoi_Scan_ToJson( const std::string& aInFile , const ManualRoI& aManualRoI , const ScanConfiguration& aScanConfig, const std::string& aOutputPattern )
{
  ManualRoi_Scan_SimpleCallback( aInFile , aManualRoI , aScanConfig , [&]( const std::string& aRoiId , const std::vector< ScanEntry >& aVector ){ _ScanCallback_Json_( aRoiId , aVector , aInFile , aOutputPattern ); } );
}



__attribute__((flatten))
void ManualRoi_Cluster_FullCallback( const std::string& aInFile , const ManualRoI& aManualRoI , const double& aR, const double& aT, const tFullClusterCallback& aCallback )
{  
  LocalizationFile( aInFile ).ExtractRoIs( aManualRoI , [&]( RoI& aRoI ) { aRoI.Clusterize( aR, aT, aCallback ); } );
}

__attribute__((flatten))
void ManualRoi_Cluster_SimpleCallback( const std::string& aInFile , const ManualRoI& aManualRoI , const double& aR, const double& aT, const tSimpleClusterCallback& aCallback )
{
  ManualRoi_Cluster_FullCallback( aInFile , aManualRoI , aR , aT , [&]( RoIproxy& aRoIproxy ){ _FullClusterToSimpleCluster_( aRoIproxy , aCallback ); } );
}

__attribute__((flatten))
void ManualRoi_Cluster_ToJson( const std::string& aInFile , const ManualRoI& aManualRoI , const double& aR, const double& aT, const std::string& aOutputPattern )
{
  ManualRoi_Cluster_SimpleCallback( aInFile , aManualRoI , aR , aT , [&]( const std::string& aRoiId , const std::vector< ClusterWrapper >& aVector ){ _ClusterCallback_Json_( aRoiId , aVector , aInFile , aOutputPattern ); } );
}



__attribute__((flatten))
void ImageJRoi_Scan_FullCallback( const std::string& aInFile , const std::string& aImageJfile , const double& aScale , const ScanConfiguration& aScanConfig, const tFullScanCallback& aCallback )
{
  LocalizationFile( aInFile ).ExtractRoIs( aImageJfile , aScale , [&]( RoI& aRoI ) { aRoI.ScanRT( aScanConfig, aCallback ); } );
}

__attribute__((flatten))
void ImageJRoi_Scan_SimpleCallback( const std::string& aInFile , const std::string& aImageJfile , const double& aScale , const ScanConfiguration& aScanConfig, const tSimpleScanCallback& aCallback )
{
  LocalizationFile( aInFile ).ExtractRoIs( aImageJfile , aScale , [&]( RoI& aRoI ) { _FullScanToSimpleScan_( aRoI , aScanConfig, aCallback ); } );
}

__attribute__((flatten))
void ImageJRoi_Scan_ToJson( const std::string& aInFile , const std::string& aImageJfile , const double& aScale , const ScanConfiguration& aScanConfig, const std::string& aOutputPattern )
{
  ImageJRoi_Scan_SimpleCallback( aInFile , aImageJfile , aScale , aScanConfig , [&]( const std::string& aRoiId , const std::vector< ScanEntry >& aVector ){ _ScanCallback_Json_( aRoiId , aVector , aInFile , aOutputPattern ); } );
}



__attribute__((flatten))
void ImageJRoi_Cluster_FullCallback( const std::string& aInFile , const std::string& aImageJfile , const double& aScale , const double& aR, const double& aT, const tFullClusterCallback& aCallback )
{  
  LocalizationFile( aInFile ).ExtractRoIs( aImageJfile , aScale , [&]( RoI& aRoI ) { aRoI.Clusterize( aR, aT, aCallback ); } );
}

__attribute__((flatten))
void ImageJRoi_Cluster_SimpleCallback( const std::string& aInFile , const std::string& aImageJfile , const double& aScale , const double& aR, const double& aT, const tSimpleClusterCallback& aCallback )
{
  ImageJRoi_Cluster_FullCallback( aInFile , aImageJfile , aScale , aR , aT , [&]( RoIproxy& aRoIproxy ){ _FullClusterToSimpleCluster_( aRoIproxy , aCallback ); } );
}

__attribute__((flatten))
void ImageJRoi_Cluster_ToJson( const std::string& aInFile , const std::string& aImageJfile , const double& aScale , const double& aR, const double& aT, const std::string& aOutputPattern )
{
  ImageJRoi_Cluster_SimpleCallback( aInFile , aImageJfile , aScale , aR , aT , [&]( const std::string& aRoiId , const std::vector< ClusterWrapper >& aVector ){ _ClusterCallback_Json_( aRoiId , aVector , aInFile , aOutputPattern ); } );
}


__attribute__((flatten))
void SegmentedImage_Cluster_FullCallback(const std::string& aInFile, const std::string& aSegmentedImagefile, const double& aScale, const double& aR, const double& aT, const tFullClusterCallback& aCallback)
{
    LocalizationFile(aInFile).ExtractRoIsFromSegmentedImage(aSegmentedImagefile, aScale, [&](RoI& aRoI) { aRoI.Clusterize(aR, aT, aCallback); });
}

__attribute__((flatten))
void SegmentedImage_Cluster_SimpleCallback(const std::string& aInFile, const std::string& aSegmentedImagefile, const double& aScale, const double& aR, const double& aT, const tSimpleClusterCallback& aCallback)
{
    SegmentedImage_Cluster_FullCallback(aInFile, aSegmentedImagefile, aScale, aR, aT, [&](RoIproxy& aRoIproxy) { _FullClusterToSimpleCluster_(aRoIproxy, aCallback); });
}

__attribute__((flatten))
void SegmentedImage_Cluster_ToJson(const std::string& aInFile, const std::string& aSegmentedImagefile, const double& aScale, const double& aR, const double& aT, const std::string& aOutputPattern)
{
    SegmentedImage_Cluster_SimpleCallback(aInFile, aSegmentedImagefile, aScale, aR, aT, [&](const std::string& aRoiId, const std::vector< ClusterWrapper >& aVector) { _ClusterCallback_Json_(aRoiId, aVector, aInFile, aOutputPattern); });
}


__attribute__((flatten))
void SegmentedImage_Scan_FullCallback(const std::string& aInFile, const std::string& aSegmentedImagefile, const double& aScale, const ScanConfiguration& aScanConfig, const tFullScanCallback& aCallback)
{
    LocalizationFile(aInFile).ExtractRoIsFromSegmentedImage(aSegmentedImagefile, aScale, [&](RoI& aRoI) { aRoI.ScanRT(aScanConfig, aCallback); });
}

__attribute__((flatten))
void SegmentedImage_Scan_SimpleCallback(const std::string& aInFile, const std::string& aSegmentedImagefile, const double& aScale, const ScanConfiguration& aScanConfig, const tSimpleScanCallback& aCallback)
{
    LocalizationFile(aInFile).ExtractRoIsFromSegmentedImage(aSegmentedImagefile, aScale, [&](RoI& aRoI) { _FullScanToSimpleScan_(aRoI, aScanConfig, aCallback); });
}

__attribute__((flatten))
void SegmentedImage_Scan_ToJson(const std::string& aInFile, const std::string& aSegmentedImagefile, const double& aScale, const ScanConfiguration& aScanConfig, const std::string& aOutputPattern)
{
    SegmentedImage_Scan_SimpleCallback(aInFile, aSegmentedImagefile, aScale, aScanConfig, [&](const std::string& aRoiId, const std::vector< ScanEntry >& aVector) { _ScanCallback_Json_(aRoiId, aVector, aInFile, aOutputPattern); });
}


__attribute__((flatten))
void SegmentedImage_FullAnalysis_ToJson(const std::string& aInFile, 
                                        const std::string& aSegmentedImagefile, 
                                        const double& aScale, 
                                        const ScanConfiguration& aScanConfig, 
                                        const std::string& aOutputPattern_Scan, const std::string& aOutputPattern_Cluster,const std::string& aOutputPattern_Info)
{
    LocalizationFile(aInFile).ExtractRoIsFromSegmentedImage(aSegmentedImagefile, aScale, [&](RoI& aRoI) { _FullAnalysis_(aRoI, aScanConfig, aOutputPattern_Scan, aOutputPattern_Cluster, aOutputPattern_Info); });
}

__attribute__((flatten))
void ImageJRoi_FullAnalysis_ToJson(const std::string& aInFile,
                                        const std::string& aImageJfile,
                                        const double& aScale,
                                        const ScanConfiguration& aScanConfig,
                                        const std::string& aOutputPattern_Scan, const std::string& aOutputPattern_Cluster, const std::string& aOutputPattern_Info)
{
    LocalizationFile(aInFile).ExtractRoIs(aImageJfile, aScale, [&](RoI& aRoI) { _FullAnalysis_(aRoI, aScanConfig, aOutputPattern_Scan, aOutputPattern_Cluster, aOutputPattern_Info); });
}

__attribute__((flatten))
void Analyse_FOV_ImageJ(const std::string& aInFile,
    const std::string& aImageJfile,
    const double& aScale,
    const std::string& aTaskDescription,
    const ScanConfiguration& aScanConfig,
    const double& aR, const double& aT,
    const std::string& aOutputPattern_Scan, const std::string& aOutputPattern_Cluster, const std::string& aOutputPattern_ClusteredLocalizations, const std::string& aOutputPattern_Info)
{
    LocalizationFile(aInFile).ExtractRoIs(aImageJfile, aScale, [&](RoI& aRoI) { _Analyse_FOV_(aRoI,
        aTaskDescription,
        aScanConfig, aR, aT,
        aOutputPattern_Scan,
        aOutputPattern_Cluster,
        aOutputPattern_ClusteredLocalizations,
        aOutputPattern_Info); });
}

__attribute__((flatten))
void Analyse_FOV_SegmentedImage(const std::string& aInFile,
    const std::string& aSegmentedImagefile,
    const double& aScale,
    const std::string& aTaskDescription,
    const ScanConfiguration& aScanConfig,
    const double& aR, const double& aT,
    const std::string& aOutputPattern_Scan, const std::string& aOutputPattern_Cluster, const std::string& aOutputPattern_ClusteredLocalizations, const std::string& aOutputPattern_Info)
{
    LocalizationFile(aInFile).ExtractRoIsFromSegmentedImage(aSegmentedImagefile, aScale, [&](RoI& aRoI) { _Analyse_FOV_(aRoI,
        aTaskDescription,
        aScanConfig, aR, aT,
        aOutputPattern_Scan,
        aOutputPattern_Cluster,
        aOutputPattern_ClusteredLocalizations,
        aOutputPattern_Info); });
}
