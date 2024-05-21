#pragma once
//#define _USE_MATH_DEFINES
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include "qumasun_mpi.h"
#include "vecmath.h"

inline
void QUMASUN_MPI::mSetDensityOne(GridRange& l_grid, double* l_rho, const double* const* l_psi, const double* occupancy, int num_solution) {

	
	//mpi//////////////////////////
	const size_t size_3d = l_grid.Size3D();



	for (int s = 0; s < num_solution; ++s) {

		auto psi_l = l_psi[s];

		if (occupancy[s] > 0.0) {

			const double factor = occupancy[s];

			// rho += conj(pr) * pr * factor, where factor is to normalize //
			for (size_t i = 0; i < size_3d; ++i) {
				l_rho[i] += (psi_l[i] * psi_l[i]) * factor;
			}

		}
	}

}


//mpi supported//
inline
double QUMASUN_MPI::mSetDensity(bool is_mixing) {
	const int proc_id = GetProcessID(m_mpi_comm);
	const bool is_root_spin = IsRoot(m_mpi_comm);
	double diff_abs_rho = 0.0;

	const double mixing_ratio = m_mixing_ratio;
	if (is_root_spin) {
		vecmath::Copy<double>(m_rho_prev, m_rho, m_size_3d);
	}

	const int local_size = ml_grid.Size3D();

	if (!is_spin_on) {
		vecmath::SetZero<double>(ml_rho, local_size);
		mSetDensityOne(ml_grid, ml_rho, ml_psi, m_occupancy, num_solution);
		mGatherField(m_rho.Pointer(), ml_rho);

		if (is_root_spin) {
			//check/////		
			double total_rho = 0.0;
			for (size_t i = 0; i < m_size_3d; ++i) {
				total_rho += m_rho[i];
			}
			total_rho *= m_dx * m_dy * m_dz;
			printf("total rho = %f\n", total_rho);// fflush(stdout);
		

			if (is_mixing) {
				for (size_t i = 0; i < m_size_3d; ++i) {
					const double diff = m_rho[i] - m_rho_prev[i];
					m_rho[i] = diff * mixing_ratio + m_rho_prev[i];
					diff_abs_rho += fabs(diff);
				}
				diff_abs_rho *= m_dx * m_dy * m_dz;
			}
#if 1
			printf("diff_abs_rho = %f\n", diff_abs_rho); fflush(stdout);
#endif
		}
        mHierarchyScatterField(ml_rho, m_rho.Pointer());

	} else { //spin polarization//

		//printf("[%d]path 1B \n", proc_id); fflush(stdout);
		const int proc_id = GetProcessID(m_mpi_comm);
		const int num_procs = GetNumProcess(m_mpi_comm);
		if (m_num_spin == 1) {
			//printf("[%d]path 1C \n", proc_id); fflush(stdout);
			vecmath::SetZero<double>(ml_rho, local_size);
			mSetDensityOne(ml_grid, ml_rho, ml_psi, m_occupancy, num_solution);
			mGatherField(m_rho.Pointer(), ml_rho);
			//spin is parallelized//
			const int TAG = 0x907;
			if (proc_id == 0) {
				MPI_Status status;
				MPI_Recv(m_rho_diff.Pointer(), m_size_3d, MPI_DOUBLE, num_procs / 2, TAG, m_mpi_comm, &status);
			} else if (proc_id == num_procs / 2) {
				MPI_Send(m_rho.Pointer(), m_size_3d, MPI_DOUBLE, 0, TAG, m_mpi_comm);
			}
		} else {
			//(m_num_spin==2)//
			//printf("[%d]path 1D \n", proc_id); fflush(stdout);
			vecmath::SetZero<double>(ml_rho, local_size);
			mSetDensityOne(ml_grid, ml_rho, ml_psi, m_occupancy, num_solution);
			mGatherField(m_rho.Pointer(), ml_rho);
			vecmath::SetZero<double>(ml_rho, local_size);
			mSetDensityOne(ml_grid, ml_rho, &ml_psi[num_solution], m_occupancy_down, num_solution);
			mGatherField(m_rho_diff.Pointer(), ml_rho);
		//	printf("[%d]m_rho_diff = 0x%zx, 0x%zx, 0x%zx\n", proc_id, m_rho.Pointer(), m_rho_diff.Pointer(), ml_rho); fflush(stdout);
		//	printf("[%d]occupancy = 0x%zx, 0x%zx, %f,%f\n", proc_id, m_occupancy, m_occupancy_down, m_occupancy[0], m_occupancy_down[0]); fflush(stdout);
		}

		const bool is_root_spin = IsRoot(m_mpi_comm);
		if (is_root_spin) {
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
			printf("[%d]total rho = %f + %f = %f\n", proc_id, total_rho_up, total_rho_dn, total_rho_up + total_rho_dn); fflush(stdout);



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

        /*
		if (m_num_spin == 1) {
			const int TAG = 2200;
			if (proc_id == 0) {
				MPI_Send(m_rho.Pointer(), m_size_3d, MPI_DOUBLE, num_procs / 2, TAG, m_mpi_comm);
			} else if (proc_id == num_procs / 2) {
				MPI_Status status;
				MPI_Recv(m_rho.Pointer(), m_size_3d, MPI_DOUBLE, 0, TAG, m_mpi_comm, &status);
			}
		}
		mScatterField(ml_rho, m_rho.Pointer());
		mScatterField(ml_rho_diff, m_rho_diff.Pointer());
        */
        mHierarchyScatterField(ml_rho, m_rho.Pointer());
        mHierarchyScatterField(ml_rho_diff, m_rho_diff.Pointer());
        
	}
	
	return diff_abs_rho;
}

