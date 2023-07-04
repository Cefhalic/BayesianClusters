//! \file LocalizationFile.hpp
#pragma once

/* ===== C++ ===== */
#include <vector>
#include <string>
#include <functional>


class Data;
class RoI;

// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//! A struct for storing the parameters of a manual RoI
struct ManualRoI {  
  double x; //!< The x-centre of the RoI
  double y; //!< The y-centre of the RoI
  double width; //!< The width of the RoI
  double height; //!< The height of the RoI
};
// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

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
    //! Manually extract an RoI 
    //! \param aRoI The manual RoI window
    //! \param aCallback A handler for each RoI found
    void ExtractRoIs( const ManualRoI& aRoI , const std::function< void( RoI& ) >& aCallback ) const;

    //! Automatically extract the RoIs
    //! \param aCallback A handler for each RoI found
    void ExtractRoIs( const std::function< void( RoI& ) >& aCallback ) const;

    //! Manually extract an RoI 
    //! \param aImageJfile The name of an ImageJ RoI file
    //! \param aScale The size of the LSB in the ImageJ file
    //! \param aCallback A handler for each RoI found
    void ExtractRoIs( const std::string& aImageJfile , const double& aScale , const std::function< void( RoI& ) >& aCallback ) const;

    //! Accessor to the raw data
    //! \return Reference to the raw data
    inline const std::vector< Data >& data()
    {
      return mData;
    }

  private:
    //! The localization file name
    std::string mFilename;

    //! The localizations in the file
    std::vector< Data > mData;
};
// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

