//! \file ImageJ_RoI.hpp
#pragma once

/* ===== C++ ===== */
#include <string>
#include <map>

/* ===== BOOST C++ ===== */
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/polygon.hpp>

//! Typedef a boost::geometry type representing an ImageJ Roi point for simplicity
typedef boost::geometry::model::point<uint16_t, 2, boost::geometry::cs::cartesian> roi_point;
//! Typedef a boost::geometry type representing an ImageJ Roi polygon for simplicity
typedef boost::geometry::model::ring<roi_point> roi_polygon;

//! Decode ImageJ binary RoI data to boost::geometry polygon
//! Reverse engineered from https://github.com/imagej/ImageJ/blob/master/ij/io/RoiDecoder.java
//! \param aData A C-array containing the binary RoI data
//! \return A boost::geometry polygon containing the RoI information
roi_polygon DecodeBinaryRoI( const uint8_t* const aData );

//! Decode ImageJ zipped binary RoI data to a map of named boost::geometry polygons
//! \param aZipFileName The name of the zip-file to be opened
//! \return A map of named boost::geometry polygons containing the RoI information
std::map< std::string , roi_polygon > OpenRoiZipfile( const std::string& aZipFileName );

//! Todo Add OpenRoiFile for completeness?
