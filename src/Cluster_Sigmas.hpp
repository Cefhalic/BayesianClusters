#pragma once

/* ===== Local utilities ===== */
#include "ListComprehension.hpp"

const std::size_t sigmacount = 10;
const double sigmaspacing = 1e-3;

const auto sigmabins  = []( const int& i ){ return i * sigmaspacing; } | range( sigmacount );
const auto sigmabins2 = []( const double& i ){ return i * i; }         | sigmabins;

const std::vector< double > p_sigma = { 0.03631079, 0.110302441, 0.214839819, 0.268302465, 0.214839819, 0.110302441, 0.03631079, 0.007664194, 0.001037236, 9.00054E-05 };
const auto log_p_sigma = []( const double& w){ return log(w); } | p_sigma;