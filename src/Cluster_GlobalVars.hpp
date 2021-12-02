#pragma once

namespace GlobalVars
{
	std::size_t sigmacount = 0;
	double sigmaspacing = 0.0;

	std::vector< double > sigmabins;
	std::vector< double > sigmabins2;

	std::vector< double > probability_sigma;
	std::vector< double > log_probability_sigma;

	double maxR;
	double maxR2;
	double max2R;
	double max2R2;
}