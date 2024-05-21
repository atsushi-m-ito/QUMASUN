#pragma once
#include <cmath>

/*
* coef is 1.0 if spin is on, 2.0 if spin is off.
*/
inline double SetOccupancy(double coef, double kbT, int num_solution, const double* eigen_values, double init_mu, double* occupancy) {
	constexpr double PMAX = 80.0;
	const double beta = 1.0 / kbT;
	double mu = init_mu;
	double sum = 0.0;
	for (int s = 0; s < num_solution; ++s) {
		const double p = std::min(beta * (eigen_values[s] - mu), PMAX);
		occupancy[s] = coef / (exp(p) + 1.0);
		sum += occupancy[s];
	}
	return sum;
}
