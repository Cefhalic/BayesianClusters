//! \file LocalizationFile.hpp
#pragma once

/* ===== C++ ===== */
#include <vector>
#include <string>
#include <functional>


class Data;
class RoI;

// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//! A class to store the raw data-points
class LocalizationFile
{
  public:
    //! Constructor
    //! \param aFilename The name of the localizations file
    LocalizationFile( const std::string& aFilename );

    //! Deleted copy constructor
    LocalizationFile( const LocalizationFile& aOther /*!< Anonymous argument */ ) = delete;

    //! Deleted assignment operator
    //! \return Reference to this, for chaining calls
    LocalizationFile& operator = (const LocalizationFile& aOther /*!< Anonymous argument */ ) = delete;

    //! Default move constructor
    LocalizationFile( LocalizationFile&& aOther /*!< Anonymous argument */ ) = default;

    //! Default move-assignment constructor
    //! \return Reference to this, for chaining calls
    LocalizationFile& operator = ( LocalizationFile&& aOther /*!< Anonymous argument */ ) = default;

    //! Default destructor
    ~LocalizationFile() = default;

  public:
    //! Automatically extract the RoIs
    //! \param aCallback A handler for each RoI found
    void ExtractRoIs( const std::function< void( RoI& ) >& aCallback );

  private:
    //! The localizations in the file
    std::vector< Data > mData;
};
// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

