//! \file API.cpp

/* ===== Cluster sources ===== */
#include "BayesianClustering/API.hpp"
#include "BayesianClustering/Configuration.hpp"

#include "BayesianClustering/RoI.hpp"


// /* ===== Local utilities ===== */
// #include "Utilities/ProgressBar.hpp"
// #include "Utilities/Vectorize.hpp"
// #include "Utilities/Units.hpp"

// /* ===== C++ ===== */
#include <iostream>
#include <functional>
#include <array>
#include <limits>



// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

// std::vector< Data > LoadLocalizationFile( const std::string& aFilename )
// {
// }




// std::vector< RoI > ExtractRoIs( const std::vector< Data >& aDataset , const tFromConfigFile& aDummy )
// {
//   std::vector< RoI > lRet;

//   auto Xmax( CurrentConfiguration().getWidthX() / 2.0 ) , Ymax( CurrentConfiguration().getWidthY() / 2.0 )

//   std::vector< Data > lData;
//   for( auto& k : aDataset )
//   {
//     double x = k.x - CurrentConfiguration().getCentreX();
//     double y = k.x - CurrentConfiguration().getCentreY();     
//     if( fabs(x) < Xmax and fabs(y) < Ymax ) lData.emplace_back( x , y , k.s );
//   }

//   lRet.emplace_back( std::move( lData ) , CurrentConfiguration() );
//   return lRet;
// }




