#pragma once
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include "qumasun_single.h"
#include "SetOccupancy.h"

inline
void QUMASUN_SINGLE::mSetOccupancy() {	
	if (is_spin_on) {
		mSetOccupancySpinOn();	
	} else {
		mSetOccupancySpinOff();	
	}
}

inline
void QUMASUN_SINGLE::mSetOccupancySpinOff() {

	const double d_num_electron = (double)num_electrons;
	//const double init_mu = m_eigen_values[num_electrons / 2];
	double mu_min = m_eigen_values[0];
	double mu_max = m_eigen_values[num_solution - 1];
	mu_max = (mu_max > 0.0) ? mu_max * 16.0 : fabs(mu_max) + 16.0;
	
	//必要最低限の軌道しか解いていない場合//
	if ((num_electrons + 1) / 2 >= num_solution) {
		mu_max *= 10.0;
	}

	int step = 0;
	double mu = (mu_min + mu_max) / 2.0;
	constexpr double THRESHOLD_MU = 1.0e-9;
	while (fabs(mu_max - mu_min) > THRESHOLD_MU){
		const double num_e = SetOccupancy(2.0, m_kbT, num_solution, m_eigen_values, mu, m_occupancy);


#ifdef DEBUG_PRINT1
		printf("occupancy_num_e(%d) = %f, mu = %f, %f, %f\n", step, num_e, mu, mu_min, mu_max);
#endif

		if (num_e < d_num_electron) {
			mu_min = mu;
		} else if(num_e > d_num_electron) {
			mu_max = mu;
		} else {
			break;
		}
		mu = (mu_min + mu_max) / 2.0;
	}

//#define DEBUG_PRINT2
#ifdef DEBUG_PRINT2

	for (int s = 0; s < num_solution; ++s) {		
		printf("occupancy[%d] = %f\n", s, m_occupancy[s]);		
	}
#endif


}

inline
void QUMASUN_SINGLE::mSetOccupancySpinOn() {

	const double d_num_electron = (double)num_electrons;
	//const double init_mu = m_eigen_values[num_electrons / 2];
	double mu_min = std::min( m_eigen_values[0], m_eigen_values_down[0] );
	double mu_max = std::max(m_eigen_values[num_solution - 1], m_eigen_values_down[num_solution - 1] );
	mu_max = (mu_max > 0.0) ? mu_max * 16.0 : fabs(mu_max) + 16.0;

	//必要最低限の軌道しか解いていない場合//
	if ((num_electrons + 1) / 2 >= num_solution) {
		mu_max *= 10.0;
	}

	int step = 0;
	double mu = (mu_min + mu_max) / 2.0;
	constexpr double THRESHOLD_MU = 1.0e-9;
	while (fabs(mu_max - mu_min) > THRESHOLD_MU) {
		const double num_e_up = SetOccupancy(1.0, m_kbT, num_solution, m_eigen_values, mu, m_occupancy);
		const double num_e_dn = SetOccupancy(1.0, m_kbT, num_solution, m_eigen_values_down, mu, m_occupancy_down);
		const double num_e = num_e_up + num_e_dn;


#ifdef DEBUG_PRINT1
		printf("occupancy_num_e(%d) = %f,%f, mu = %f, %f, %f\n", step, num_e_up, num_e_dn, mu, mu_min, mu_max);
#endif


		if (num_e < d_num_electron) {
			mu_min = mu;
		} else if (num_e > d_num_electron) {
			mu_max = mu;
		} else {
			break;
		}
		mu = (mu_min + mu_max) / 2.0;
	}

//#define DEBUG_PRINT2
#ifdef DEBUG_PRINT2
	for (int s = 0; s < num_solution; ++s) {
		printf("occupancy[%d] = %f, %f\n", s, m_occupancy[s], m_occupancy_down[s]);
	}
#endif



}


