#pragma once

/* ===== C++ ===== */
#include <map>
#include <vector>
#include <sstream>
#include <string>
#include <functional>

//! Define a constant for converting nanometers to meters
constexpr double nanometer  = 1e-9;
//! Define a constant for converting micrometers to meters
constexpr double micrometer = 1e-6;
//! Define a constant for converting millimeters to meters
constexpr double millimeter = 1e-3;
//! Define a constant for converting meters to meters
constexpr double meter      = 1e-0;

//! A map for converting string representations of SI units to scaling factors
const std::map< std::string , double > UnitMap{ {"nm",nanometer} , {"um",micrometer} , {"mm",millimeter} , {"m",meter} };

//! User-defined literals fot nanometer quantities
//! \param aVal The specified value
//! \return The literal value
constexpr long double operator"" _nanometer( long double aVal )
{
	return aVal * nanometer;
}

//! User-defined literals fot nanometer quantities
//! \param aVal The specified value
//! \return The literal value
constexpr long double operator"" _nanometer( unsigned long long aVal )
{
	return aVal * nanometer;
}

//! User-defined literals fot micrometer quantities
//! \param aVal The specified value
//! \return The literal value
constexpr long double operator"" _micrometer( long double aVal )
{
	return aVal * micrometer;
}

//! User-defined literals fot micrometer quantities
//! \param aVal The specified value
//! \return The literal value
constexpr long double operator"" _micrometer( unsigned long long aVal )
{
	return aVal * micrometer;
}

//! Convert a string representation to a distance
//! \param aStr A string representation of a distance
//! \return The literal value
inline long double StrToDist( const std::string& aStr )
{
	std::stringstream lStr;
	lStr << aStr;
	double lVal;
	std::string lUnits;
	lStr >> lVal >> lUnits;
	return lVal * UnitMap.at( lUnits );
}




//! Class for storing the configuration parameters
class Configuration
{
public:
  //! Default constructor
	Configuration();

  //! Setter for the centre of the scan window
  //! \param aPhysicalCentreX The x-coordinate of the centre of the window in physical units (becomes 0 in algorithm units)
  //! \param aPhysicalCentreY The y-coordinate of the centre of the window in physical units (becomes 0 in algorithm units)
	void SetCentre( const double& aPhysicalCentreX , const double& aPhysicalCentreY );

  //! Setter for the size of the RoI window
  //! \param aWidthX The width of the window in physical units
  //! \param aWidthY The height of the window in physical units
  void SetWidth( const double& aWidthX , const double& aWidthY );

  //! Setter for the sigma-bins to be integrated over
  //! \param aSigmacount   The number of sigma bins
  //! \param aSigmaMin     The lowest sigma bin
  //! \param aSigmaMax     The highest sigma bin
  //! \param aInterpolator Function-object to generate the probability of any given sigma
	void SetSigmaParameters( const std::size_t& aSigmacount , const double& aSigmaMin , const double& aSigmaMax , const std::function< double( const double& ) >& aInterpolator );

  //! Setter for the R bins for the RT scan
  //! \param aRbins    The number of R bins to scan over
  //! \param aMinScanR The lowest value of R to scan
  //! \param aMaxScanR The largest value of R to scan
	void SetRBins( const std::size_t& aRbins , const double& aMinScanR = 0.0 , const double& aMaxScanR = -1 );

  //! \param aTbins    The number of T bins to scan over
  //! \param aMinScanT The lowest value of T to scan
  //! \param aMaxScanT The largest value of T to scan
  void SetTBins( const std::size_t& aTbins , const double& aMinScanT = 0.0 , const double& aMaxScanT = -1 );

  //! Setter for the P_b parameter
  //! \param aPB    The P_b parameter
	void SetPb( const double& aPB );

  //! Setter for the alpha parameter
  //! \param aAlpha The alpha parameter
  void SetAlpha( const double& aAlpha );
  
  //! Set whether to validate clusterization
  //! \param aValidate Whether to validate clusterization
	void SetValidate( const bool& aValidate );

  //! Setter for the input file 
  //! \param aFileName The name of the file 
  void SetInputFile( const std::string& aFileName );
  //! Setter for the output file 
  //! \param aFileName The name of the file
  void SetOutputFile( const std::string& aFileName );

  //! Parse the parameters when passed in as commandline arguments
  //! \param argc The number of commandline arguments
  //! \param argv The commandline arguments
	void FromCommandline( int argc , char **argv );

  //! Parse the parameters when passed in as commandline arguments
  //! \param aArgs The commandline arguments
  void FromVector( const std::vector< std::string >& aArgs );


public:
  //! Getter for the sigma count
  //! \return The sigma count
	inline const std::size_t& sigmacount() const { return mSigmacount; }

  //! Getter for the sigma spacing
  //! \return The sigma spacing  
	inline const double& sigmaspacing() const { return mSigmaspacing; }

  //! Getter for the values of sigma
  //! \return The values of sigma
	inline const std::vector< double >& sigmabins( ) const { return mSigmabins; }
  //! Getter for the values of sigma squared
  //! \return The values of sigma squared 
	inline const std::vector< double >& sigmabins2( ) const { return mSigmabins2; }
  //! Getter for the probabilities of a given sigma
  //! \return The probabilities of given sigma
	inline const std::vector< double >& probability_sigma( ) const { return mProbabilitySigma; }
  //! Getter for the log of the probabilities of a given sigma
  //! \return The log of the probabilities of given sigma
	inline const std::vector< double >& log_probability_sigma( ) const { return mLogProbabilitySigma; }

  //! Getter for the i'th value of sigma
  //! \param i The index of the value of sigma to get
  //! \return The value of sigma_i
	inline const double& sigmabins( const std::size_t& i ) const { return mSigmabins[i]; }
  //! Getter for the i'th value of sigma squared
  //! \param i The index of the value of sigma squared to get
  //! \return The value of sigma_i squared
	inline const double& sigmabins2( const std::size_t& i ) const { return mSigmabins2[i]; }
  //! Getter for the probability of the i'th value of sigma
  //! \param i The index of the value of sigma to get the probability for
  //! \return The probability of sigma_i
	inline const double& probability_sigma( const std::size_t& i ) const { return mProbabilitySigma[i]; }
  //! Getter for the log-probability of the i'th value of sigma
  //! \param i The index of the value of sigma to get the log-probability for
  //! \return The log-probability of sigma_i
	inline const double& log_probability_sigma( const std::size_t& i ) const { return mLogProbabilitySigma[i]; }

  //! Getter for the maximum value of R
  //! \return The maximum value of R
	inline const double& maxR() const { return mMaxR; }
  //! Getter for the maximum value of R squared
  //! \return The maximum value of R squared
	inline const double& maxR2() const { return mMaxR2; }
  //! Getter for the maximum value of 2R
  //! \return The maximum value of 2R
	inline const double& max2R() const { return mMax2R; }
  //! Getter for the maximum value of 2R squared
  //! \return The maximum value of 2R squared
	inline const double& max2R2() const { return mMax2R2; }

  //! Getter for the lowest value of R to scan 
  //! \return The lowest value of R to scan 
	inline const double& minScanR() const { return mMinScanR; }
  //! Getter for the highest value of R to scan 
  //! \return The highest value of R to scan 
	inline const double& maxScanR() const { return mMaxScanR; }
  //! Getter for the lowest value of T to scan 
  //! \return The lowest value of T to scan 
	inline const double& minScanT() const { return mMinScanT; }
  //! Getter for the highest value of T to scan 
  //! \return The highest value of T to scan 
	inline const double& maxScanT() const { return mMaxScanT; }

  //! Getter for the spacing of value of R to scan 
  //! \return The spacing of value of R to scan 
	inline const double& dR() const { return mDR; }
  //! Getter for the number of R values to scan
  //! \return The number of R values to scan
	inline const std::size_t& Rbins() const { return mRbins; }
  //! Getter for the spacing of value of T to scan 
  //! \return The spacing of value of T to scan 
	inline const double& dT() const { return mDT; }
  //! Getter for the number of T values to scan
  //! \return The number of T values to scan
	inline const std::size_t& Tbins() const { return mTbins; }

  //! Logarithm of the P_b parameter  
  //! \return Logarithm of the P_b parameter 
	inline const double& logPb() const { return mLogPb; }
  //! Logarithm of the ( 1 - P_b ) parameter  
  //! \return Logarithm of the ( 1 - P_b ) parameter  
	inline const double& logPbDagger() const { return mLogPbDagger; }

  //! Getter for the alpha parameter  
  //! \return The alpha parameter  
	inline const double& alpha() const { return mAlpha; }
  //! Getter for the logarithm of the alpha parameter  
  //! \return The logarithm of the alpha parameter  
	inline const double& logAlpha() const { return mLogAlpha; }
  //! Getter for the logarithm of the gamma function of alpha parameter  
  //! \return The logarithm of the gamma function of alpha parameter  
	inline const double& logGammaAlpha() const { return mLogGammaAlpha; }

  //! Getter for whether or not to run the validation on the clustering 
  //! \return Whether or not to run the validation on the clustering 
	inline const bool& validate() const { return mValidate; }


  //! Getter for the input file 
  //! \return The name of the input RoI file
  inline const std::string& inputFile() const { return mInputFile; }
  //! Getter for the output file 
  //! \return The name of the output file
  inline const std::string& outputFile() const { return mOutputFile; }


  //! Getter for the R value for a clusterization pass
  //! \return The R value for a clusterization pass
  inline const double& ClusterR() const { return mClusterR; }
  //! Getter for the T value for a clusterization pass
  //! \return The T value for a clusterization pass
  inline const double& ClusterT() const { return mClusterT; }

  //! Getter for the x-coordinate of the physical centre
  //! \return The x-coordinate of the physical centre
  double getCentreX() const {return mPhysicalCentreX;}
  //! Getter for the y-coordinate of the physical centre
  //! \return The y-coordinate of the physical centre
  double getCentreY() const {return mPhysicalCentreY;}

  //! Getter for the width of the ROI window
  //! \return The width of the ROI window
  double getWidthX() const { return mWidthX; }

  //! Getter for the height of the ROI window
  //! \return The height of the ROI window
  double getWidthY() const { return mWidthY; }

  //! Getter for the height of the ROI window
  //! \return The height of the ROI window
  double getArea() const { return mArea; }



  void Rezoom( const double& aScale );



public:
  //! A single global copy of the global variables
  static Configuration* Current;

//   //! Getter for the singleton instance
//   //! \return The singleton instance
//   inline static Configuration& getInstance()
//   {
//     return Instance;
//   }


private:
  //! The x-coordinate of the centre of the window in physical units
  double mPhysicalCentreX;
  //! The y-coordinate of the centre of the window in physical units 
  double mPhysicalCentreY;

  //! The width of the window in the x-direction in physical units
  double mWidthX;
  //! The width of the window in the y-direction in physical units
  double mWidthY;
  //! The area of the window in physical units 
  double mArea;

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

  //! The maximum value of R
	double mMaxR;
  //! The maximum value of R squared
  double mMaxR2;
  //! The maximum value of 2R  
  double mMax2R;
  //! The maximum value of 2R squared
  double mMax2R2;

  //! The lowest value of R to scan 
	double mMinScanR;
  //! The largest value of R to scan
  double mMaxScanR;
  //! The lowest value of T to scan
  double mMinScanT;
  //! The largest value of T to scan
  double mMaxScanT;
  //! The spacing of value of R to scan 
	double mDR;
  //! The spacing of value of T to scan 
  double mDT;
  //! The number of R values to scan 
	std::size_t mRbins;
  //! The number of T values to scan 
	std::size_t mTbins;

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

  //! Whether or not to run the validation on the clustering 
	bool mValidate;

  //! The input RoI file
  std::string mInputFile;

  //! The output file 
  std::string mOutputFile;

  //! The value of R for clustering
  double mClusterR;
  //! The value of T for clustering
  double mClusterT;
};


inline Configuration& CurrentConfiguration()
{
  if( ! Configuration::Current ) throw std::runtime_error( "Current configuration not set" );
  return *Configuration::Current;
}

inline void SetCurrentConfiguration( Configuration& aConfiguration )
{
  Configuration::Current = & aConfiguration;
}
