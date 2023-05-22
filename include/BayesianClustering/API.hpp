#pragma once

/* ===== C++ ===== */
#include <string>
#include <vector>

/* ===== Cluster sources ===== */

class Data;
class RoI;
class Configuration;

//! API to load a datafile
//! \param aFilename The name of the file to load   
//! \return The vector of datapoints
std::vector< Data > LoadLocalizationFile( const std::string& aFilename ); 


class tFromConfigFile{}; static const tFromConfigFile FromConfigFile;
class tAuto{}; static const tAuto Auto;


std::vector< RoI > ExtractRoIs( const std::vector< Data >& aDataset , const tFromConfigFile& aDummy ); 
std::vector< RoI > ExtractRoIs( const std::vector< Data >& aDataset , const tAuto& aDummy ); 


