#pragma once

/* ===== C++ ===== */
#include <map>
#include <string>

//! Define a constant for converting nanometers to meters
constexpr double nanometer  = 1e-9;
//! Define a constant for converting micrometers to meters
constexpr double micrometer = 1e-6;
//! Define a constant for converting millimeters to meters
constexpr double millimeter = 1e-3;
//! Define a constant for converting meters to meters
constexpr double meter      = 1e-0;

//! A map for converting string representations of SI units to scaling factors
extern const std::map< std::string , double > UnitMap;

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
long double StrToDist( const std::string& aStr );
