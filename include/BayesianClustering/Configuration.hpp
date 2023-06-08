//! \file Configuration.hpp
#pragma once

/* ===== C++ ===== */
#include <vector>
#include <string>
#include <functional>
#include <stdexcept>

//! Class for storing the scan configuration parameters
class ScanConfiguration
{
  public:

    //! A struct to store the bounds of a scan in either R or T
    struct tBounds {
      //! The lowest value of R to scan
      double min;
      //! The largest value of R to scan
      double max;
      //! The spacing of value of R to scan
      double spacing;
      //! The number of R values to scan
      std::size_t bins;
    };

    // //! Constructor which parses the parameters when passed in as commandline arguments
    // //! \param argc The number of commandline arguments
    // //! \param argv The commandline arguments
    // ScanConfiguration( int argc, char** argv );

    // //! Constructor which parses the parameters when passed in as commandline arguments
    // //! \param aArgs The commandline arguments
    // ScanConfiguration( const std::vector< std::string >& aArgs );

    //! Constructor which parses the parameters when passed in as commandline arguments
    //! \param aCfgFile A Scan-parameter config file name
    ScanConfiguration( const std::string& aCfgFile );

    //! Constructor which take the parameters directly
    //! \param aSigmacount   The number of sigma bins
    //! \param aSigmaMin     The lowest sigma bin
    //! \param aSigmaMax     The highest sigma bin
    //! \param aInterpolator Function-object to generate the probability of any given sigma
    //! \param aRbins    The number of R bins to scan over
    //! \param aMinScanR The lowest value of R to scan
    //! \param aMaxScanR The largest value of R to scan
    //! \param aTbins    The number of T bins to scan over
    //! \param aMinScanT The lowest value of T to scan
    //! \param aMaxScanT The largest value of T to scan
    //! \param aPB    The P_b parameter
    //! \param aAlpha The alpha parameter
    ScanConfiguration( 
      const std::size_t& aSigmacount, const double& aSigmaMin, const double& aSigmaMax, const std::function< double( const double& ) >& aInterpolator ,
      const std::size_t& aRbins, const double& aMinScanR, const double& aMaxScanR ,
      const std::size_t& aTbins, const double& aMinScanT, const double& aMaxScanT ,
      const double& aPB , const double& aAlpha );

    //! Deleted copy constructor
    ScanConfiguration( const ScanConfiguration& aOther /*!< Anonymous argument */ ) = delete;

    //! Deleted assignment operator
    //! \return Reference to this, for chaining calls
    ScanConfiguration& operator = (const ScanConfiguration& aOther /*!< Anonymous argument */ ) = delete;

    //! Default move constructor
    ScanConfiguration( ScanConfiguration&& aOther /*!< Anonymous argument */ ) = default;

    //! Default move-assignment constructor
    //! \return Reference to this, for chaining calls
    ScanConfiguration& operator = ( ScanConfiguration&& aOther /*!< Anonymous argument */ ) = default;

    //! Default destructor
    ~ScanConfiguration() = default;


    //! Setter for the sigma-bins to be integrated over
    //! \param aSigmacount   The number of sigma bins
    //! \param aSigmaMin     The lowest sigma bin
    //! \param aSigmaMax     The highest sigma bin
    //! \param aInterpolator Function-object to generate the probability of any given sigma
    void SetSigmaParameters( const std::size_t& aSigmacount, const double& aSigmaMin, const double& aSigmaMax, const std::function< double( const double& ) >& aInterpolator );

    //! Setter for the R bins for the RT scan
    //! \param aRbins    The number of R bins to scan over
    //! \param aMinScanR The lowest value of R to scan
    //! \param aMaxScanR The largest value of R to scan
    void SetRBins( const std::size_t& aRbins, const double& aMinScanR, const double& aMaxScanR );

    //! \param aTbins    The number of T bins to scan over
    //! \param aMinScanT The lowest value of T to scan
    //! \param aMaxScanT The largest value of T to scan
    void SetTBins( const std::size_t& aTbins, const double& aMinScanT, const double& aMaxScanT );

    //! Setter for the P_b parameter
    //! \param aPB    The P_b parameter
    void SetPb( const double& aPB );

    //! Setter for the alpha parameter
    //! \param aAlpha The alpha parameter
    void SetAlpha( const double& aAlpha );

    // //! Parse the parameters when passed in as commandline arguments
    // //! \param argc The number of commandline arguments
    // //! \param argv The commandline arguments
    // void FromCommandline( int argc, char** argv );

  private:
    //! Parse the parameters when passed in as commandline arguments
    //! \param aArgs The commandline arguments
    void FromVector( const std::vector< std::string >& aArgs );


  public:
    //! Getter for the sigma count
    //! \return The sigma count
    inline const std::size_t& sigmacount() const
    {
      return mSigmacount;
    }

    //! Getter for the sigma spacing
    //! \return The sigma spacing
    inline const double& sigmaspacing() const
    {
      return mSigmaspacing;
    }

    //! Getter for the values of sigma
    //! \return The values of sigma
    inline const std::vector< double >& sigmabins( ) const
    {
      return mSigmabins;
    }
    //! Getter for the values of sigma squared
    //! \return The values of sigma squared
    inline const std::vector< double >& sigmabins2( ) const
    {
      return mSigmabins2;
    }
    //! Getter for the probabilities of a given sigma
    //! \return The probabilities of given sigma
    inline const std::vector< double >& probability_sigma( ) const
    {
      return mProbabilitySigma;
    }
    //! Getter for the log of the probabilities of a given sigma
    //! \return The log of the probabilities of given sigma
    inline const std::vector< double >& log_probability_sigma( ) const
    {
      return mLogProbabilitySigma;
    }

    //! Getter for the i'th value of sigma
    //! \param i The index of the value of sigma to get
    //! \return The value of sigma_i
    inline const double& sigmabins( const std::size_t& i ) const
    {
      return mSigmabins[i];
    }
    //! Getter for the i'th value of sigma squared
    //! \param i The index of the value of sigma squared to get
    //! \return The value of sigma_i squared
    inline const double& sigmabins2( const std::size_t& i ) const
    {
      return mSigmabins2[i];
    }
    //! Getter for the probability of the i'th value of sigma
    //! \param i The index of the value of sigma to get the probability for
    //! \return The probability of sigma_i
    inline const double& probability_sigma( const std::size_t& i ) const
    {
      return mProbabilitySigma[i];
    }
    //! Getter for the log-probability of the i'th value of sigma
    //! \param i The index of the value of sigma to get the log-probability for
    //! \return The log-probability of sigma_i
    inline const double& log_probability_sigma( const std::size_t& i ) const
    {
      return mLogProbabilitySigma[i];
    }

    //! Getter for the bounds of R to scan
    //! \return The lbounds of R to scan
    inline const tBounds& Rbounds() const
    {
      return mRbounds;
    }
    //! Getter for the bounds of T to scan
    //! \return The lbounds of T to scan
    inline const tBounds& Tbounds() const
    {
      return mTbounds;
    }

    //! Logarithm of the P_b parameter
    //! \return Logarithm of the P_b parameter
    inline const double& logPb() const
    {
      return mLogPb;
    }
    //! Logarithm of the ( 1 - P_b ) parameter
    //! \return Logarithm of the ( 1 - P_b ) parameter
    inline const double& logPbDagger() const
    {
      return mLogPbDagger;
    }

    //! Getter for the alpha parameter
    //! \return The alpha parameter
    inline const double& alpha() const
    {
      return mAlpha;
    }
    //! Getter for the logarithm of the alpha parameter
    //! \return The logarithm of the alpha parameter
    inline const double& logAlpha() const
    {
      return mLogAlpha;
    }
    //! Getter for the logarithm of the gamma function of alpha parameter
    //! \return The logarithm of the gamma function of alpha parameter
    inline const double& logGammaAlpha() const
    {
      return mLogGammaAlpha;
    }

  private:
    //! The number of sigma bins
    std::size_t mSigmacount;
    //! The spacing of sigma bins
    double mSigmaspacing;
    //! The values of sigma
    std::vector< double > mSigmabins;
    //! The values of sigma squared
    std::vector< double > mSigmabins2;
    //! The probability of a given sigma
    std::vector< double > mProbabilitySigma;
    //! The log-probability of a gievn sigma
    std::vector< double > mLogProbabilitySigma;

    //! The bounds of R to scan
    tBounds mRbounds;
    //! The bounds of T to scan
    tBounds mTbounds;

    //! The alpha parameter
    double mAlpha;
    //! Logarithm of the alpha parameter
    double mLogAlpha;
    //! Logarithm of the gamma function of alpha parameter
    double mLogGammaAlpha;

    //! Logarithm of the P_b parameter
    double mLogPb;
    //! Logarithm of the( 1- P_b ) parameter
    double mLogPbDagger;

};


//! Class for storing the auxilliary configuration parameters
class AuxConfiguration
{
  public:
    //! Default constructor
    // AuxConfiguration();

    //! Constructor which parses the parameters when passed in as commandline arguments
    //! \param argc The number of commandline arguments
    //! \param argv The commandline arguments
    AuxConfiguration( int argc, char** argv );

    //! Constructor which parses the parameters when passed in as commandline arguments
    //! \param aArgs The commandline arguments
    AuxConfiguration( const std::vector< std::string >& aArgs );


    //! Set whether to validate clusterization
    //! \param aValidate Whether to validate clusterization
    void SetValidate( const bool& aValidate );
    //! Setter for the input file
    //! \param aFileName The name of the file
    void SetInputFile( const std::string& aFileName );
    //! Setter for the output file
    //! \param aFileName The name of the file
    void SetOutputFile( const std::string& aFileName );
    //! Setter for the config file
    //! \param aFileName The name of the file
    void SetConfigFile( const std::string& aFileName );


    //! Getter for whether or not to run the validation on the clustering
    //! \return Whether or not to run the validation on the clustering
    inline const bool& validate() const
    {
      return mValidate;
    }

    //! Getter for the input file
    //! \return The name of the input RoI file
    inline const std::string& inputFile() const
    {
      return mInputFile;
    }

    //! Getter for the output file
    //! \return The name of the output file
    inline const std::string& outputFile() const
    {
      return mOutputFile;
    }

    //! Getter for the config file
    //! \return The name of the config file
    inline const std::string& configFile() const
    {
      return mConfigFile;
    }

    //! Getter for the R value for a clusterization pass
    //! \return The R value for a clusterization pass
    inline const double& ClusterR() const
    {
      return mClusterR;
    }

    //! Getter for the T value for a clusterization pass
    //! \return The T value for a clusterization pass
    inline const double& ClusterT() const
    {
      return mClusterT;
    }

    // //! Parse the parameters when passed in as commandline arguments
    // //! \param argc The number of commandline arguments
    // //! \param argv The commandline arguments
    // void FromCommandline( int argc, char** argv );

  private:
    //! Parse the parameters when passed in as commandline arguments
    //! \param aArgs The commandline arguments
    void FromVector( const std::vector< std::string >& aArgs );

  public:
    //! Whether or not to run the validation on the clustering
    bool mValidate;

    //! The input RoI file
    std::string mInputFile;

    //! The output file
    std::string mOutputFile;

    //! The config file
    std::string mConfigFile;

    //! The value of R for clustering
    double mClusterR;
    //! The value of T for clustering
    double mClusterT;
};
