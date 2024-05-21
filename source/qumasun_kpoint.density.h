#pragma once
//#define _USE_MATH_DEFINES
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include "qumasun_kpoint.h"
#include "vecmath.h"
#include "symmetrize_density.h"

inline
void QUMASUN_KPOINT::mSetDensityOne(GridRange& l_grid, double* l_rho, const SoAComplex* l_psi, const double* occupancy, int num_solution) {

	const size_t size_3d = l_grid.Size3D();


	for (int s = 0; s < num_solution; ++s) {

		auto psi_l_re = l_psi[s].re;
		auto psi_l_im = l_psi[s].im;

		if (occupancy[s] > 0.0) {

			const double factor = occupancy[s] ;

			// rho += conj(pr) * pr * factor, where factor is to normalize //
			for (size_t i = 0; i < size_3d; ++i) {
				l_rho[i] += (psi_l_re[i] * psi_l_re[i] + psi_l_im[i] * psi_l_im[i]) * factor;
			}

		}
	}

}


//mpi supported//
inline
double QUMASUN_KPOINT::mSetDensity(bool is_mixing) {

	
	const bool is_root_global = IsRoot(m_mpi_comm);
		
	double diff_abs_rho = 0.0;
	const size_t size_3d = ml_grid.Size3D();


	const double mixing_ratio = m_mixing_ratio;
	if (is_root_global) {
		vecmath::Copy<double>(m_rho_prev, m_rho, m_size_3d);
	}

	vecmath::SetZero<double>(ml_rho, size_3d);
	if (is_spin_on) {
		vecmath::SetZero<double>(ml_rho_diff, size_3d);
	}

	for (int sk = 0; sk < m_num_having_spin_kpoint; ++sk) {
		if (ml_wave_set[sk].spin == SPIN::UP) {
			mSetDensityOne(ml_grid, ml_rho, ml_wave_set[sk].l_psi_set, ml_wave_set[sk].occupancy, num_solution);
		} else {
			mSetDensityOne(ml_grid, ml_rho_diff, ml_wave_set[sk].l_psi_set, ml_wave_set[sk].occupancy, num_solution);
		}
	}

	mHierarchyGatherField(m_rho.Pointer(), ml_rho);
	if (is_spin_on) {
		mHierarchyGatherField(m_rho_diff.Pointer(), ml_rho_diff);
	}
		
	if (is_root_global) {

		if (kpoint_symmetry != QUMASUN::KPOINT_SYMMETRY::NONE) {
			SymmetrizeDensity(m_rho.Pointer(), m_global_grid.SizeX(), m_global_grid.SizeY(), m_global_grid.SizeZ(), kpoint_symmetry);
			if (is_spin_on) {
				SymmetrizeDensity(m_rho_diff.Pointer(), m_global_grid.SizeX(), m_global_grid.SizeY(), m_global_grid.SizeZ(), kpoint_symmetry);
			}
		}


		if (is_spin_on) { //spin polarization//

			
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
		} else {
			//check/////		
			double total_rho = 0.0;
			for (size_t i = 0; i < m_size_3d; ++i) {
				total_rho += m_rho[i];
			}
			total_rho *= m_dx * m_dy * m_dz;
			printf("total rho = %f\n", total_rho);
		}


		if (is_mixing) {
			for (size_t i = 0; i < m_size_3d; ++i) {
				const double diff = m_rho[i] - m_rho_prev[i];
				m_rho[i] = diff * mixing_ratio + m_rho_prev[i];
				diff_abs_rho += fabs(diff);
			}
			diff_abs_rho *= m_dx * m_dy * m_dz;

#if 1
			printf("diff_abs_rho = %f\n", diff_abs_rho); fflush(stdout);
#endif

		}

	}

	mHierarchyScatterField(ml_rho, m_rho.Pointer());
	if (is_spin_on) {
		mHierarchyScatterField(ml_rho_diff, m_rho_diff.Pointer());
	}
	
	return diff_abs_rho;
}

