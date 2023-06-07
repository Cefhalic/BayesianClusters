//! \file API.hpp
#pragma once

/* ===== C++ ===== */
#include <string>
#include <vector>
#include <functional>

/* ===== Cluster sources ===== */

// class Data;
class RoI;
class ScanConfiguration;
struct ScanEntry;

void AutoRoi_Scan_SimpleCallback( const std::string& aFilename , const ScanConfiguration& aScanConfig, const std::function< void( const std::vector< ScanEntry >&  ) >& aCallback );

// //! API to load a datafile
// //! \param aFilename The name of the file to load
// //! \return The vector of datapoints
// std::vector< Data > LoadLocalizationFile( const std::string& aFilename );


// class tFromConfigFile{}; static const tFromConfigFile FromConfigFile;
// class tAuto{}; static const tAuto Auto;


// // std::vector< RoI > ExtractRoIs( const std::vector< Data >& aDataset , const tFromConfigFile& aDummy );
// std::vector< RoI > ExtractRoIs( const std::vector< Data >& aDataset , const tAuto& aDummy );



// template < typename T > class Handler;

// class BaseHandler
// {
// public:
//  BaseHandler() : mPrevious( NULL ) , mNext( NULL )
//  {}

//  template < typename T >
//  BaseHandler& operator>> ( T&& aNext )
//  {
//    mNext = new T( std::move( aNext ) );
//    return *mNext;
//  }

//  virtual ~BaseHandler()
//  {
//    if( mNext ) delete mNext;
//    mNext = NULL;
//  }

//  virtual void run()
//  {
//    if( mPrevious ) mPrevious -> run();
//    else throw std::runtime( "No previous defined" );
//  }

//  Handler<T>& Do

// protected:
//  BaseHandler* mPrevious;

// private:
//  virtual void Dummy() = 0;
//  BaseHandler* mNext;
// };


// template < typename T >
// class Handler : public BaseHandler
// {
// public:
//  Handler() = default;

//  virtual ~Handler() = default;

//  virtual void handle( T& aArg ) = 0;

// private:
//  void Dummy(){};
// };


