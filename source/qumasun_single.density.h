#pragma once
//#define _USE_MATH_DEFINES
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include "qumasun_single.h"
#include "vecmath.h"

inline
void QUMASUN_SINGLE::mSetDensityOne(size_t size_3d, double* rho, RspaceFunc<double> const* psi_set, const double* occupancy, int num_solution){

	
	for (int s = 0; s < num_solution; ++s) {
		const RspaceFunc<double>& psi_r = psi_set[s];

		if (occupancy[s] <= 0.0) continue;		
		
		const double factor = occupancy[s];
		// rho += conj(pr) * pr * factor, where factor is to normalize //
		for (size_t i = 0; i < size_3d; ++i) {			
			rho[i] += (psi_r[i]*psi_r[i]) * factor;
		}
				
	}

}


/*
* 更新前の電子密度との差を返す
* 
*/
inline
double QUMASUN_SINGLE::mSetDensity(bool is_mixing) {

	const double mixing_ratio = 0.5;
	vecmath::Copy<double>(m_rho_prev, m_rho, m_size_3d);

	if (!is_spin_on) {

		vecmath::SetZero<double>(m_rho, m_size_3d);
		mSetDensityOne(m_size_3d, m_rho, m_psi_set, m_occupancy, num_solution);
		
		//check/////
		{
			double total_rho = 0.0;
			for (size_t i = 0; i < m_size_3d; ++i) {
				total_rho += m_rho[i];
			}
			total_rho *= m_dx * m_dy * m_dz;
			printf("total rho = %f\n", total_rho);
		}
		
	}else{
		vecmath::SetZero<double>(m_rho, m_size_3d);
		vecmath::SetZero<double>(m_rho_diff, m_size_3d);
		mSetDensityOne(m_size_3d, m_rho, m_psi_set, m_occupancy, num_solution);
		mSetDensityOne(m_size_3d, m_rho_diff, &m_psi_set[num_solution], m_occupancy_down, num_solution);
		
		double total_rho_up = 0.0;
		double total_rho_dn = 0.0;
		for (size_t i = 0; i < m_size_3d; ++i) {
			total_rho_up += m_rho[i];
			total_rho_dn += m_rho_diff[i];

			const double plus_rho = m_rho[i] + m_rho_diff[i];
			const double diff_rho = m_rho[i] - m_rho_diff[i];
			m_rho[i] = plus_rho;
			m_rho_diff[i] = diff_rho;
		}
		total_rho_up *= m_dx * m_dy * m_dz;
		total_rho_dn *= m_dx * m_dy * m_dz;
		printf("total rho = %f + %f = %f\n", total_rho_up, total_rho_dn, total_rho_up + total_rho_dn);

	}

	double diff_abs_rho = 0.0;
	if (is_mixing) {
		for (size_t i = 0; i < m_size_3d; ++i) {
			const double diff = m_rho[i] - m_rho_prev[i];
			m_rho[i] = diff * mixing_ratio + m_rho_prev[i];
			diff_abs_rho += fabs(diff);
		}
		diff_abs_rho *= m_dx * m_dy * m_dz;
	}
	return diff_abs_rho;
}

