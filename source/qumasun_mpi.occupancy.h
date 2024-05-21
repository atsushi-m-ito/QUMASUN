#pragma once
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include "qumasun_mpi.h"
#include "SetOccupancy.h"

//mpi supported//
inline
void QUMASUN_MPI::mSetOccupancy() {
	const int proc_id = GetProcessID(m_mpi_comm);
	const bool is_root_spin = IsRoot(m_mpi_comm);


	if (is_spin_on) {
		if (is_root_spin) {
			mSetOccupancySpinOn();
		}


		if (m_num_spin == 1) {
			//printf("[%d] mSetOccupancy 1\n", proc_id); fflush(stdout);
			//mpi parallel for spin//
			const int TAG = 0x908; 
			const int proc_id = GetProcessID(m_mpi_comm);
			const int num_procs = GetNumProcess(m_mpi_comm);
			if (proc_id == 0) {
				MPI_Send(m_occupancy_down, num_solution, MPI_DOUBLE, num_procs/2, TAG, m_mpi_comm);
			} else if(proc_id==num_procs/2){
				MPI_Status status;
				MPI_Recv(m_occupancy, num_solution, MPI_DOUBLE, 0, TAG, m_mpi_comm, &status);
			}

			//printf("[%d] mSetOccupancy 2: 0x%zx\n", proc_id, m_occupancy); fflush(stdout);
			MPI_Bcast(m_occupancy, num_solution, MPI_DOUBLE, 0, m_ddm_comm);
			//printf("[%d] mSetOccupancy 3\n", proc_id); fflush(stdout);
			
		} else {
			//spin is not parallelized//

			//printf("[%d] mSetOccupancy 4\n", proc_id); fflush(stdout);
			//forDDM//
			MPI_Bcast(m_occupancy, num_solution, MPI_DOUBLE, 0, m_ddm_comm);
			MPI_Bcast(m_occupancy_down, num_solution, MPI_DOUBLE, 0, m_ddm_comm);
			//printf("[%d: res=]occupancy = 0x%zx, 0x%zx, %f,%f\n", proc_id, m_occupancy, m_occupancy_down, m_occupancy[num_solution - 1], m_occupancy_down[num_solution - 1]); fflush(stdout);
		}
	} else {
		if (is_root_spin) {
			mSetOccupancySpinOff();
		}
		//forDDM//
		MPI_Bcast(m_occupancy, num_solution, MPI_DOUBLE, 0, m_ddm_comm);
	}
	
}



inline
void QUMASUN_MPI::mSetOccupancySpinOff() {


	const double d_num_electron = (double)num_electrons;
	//const double init_mu = m_eigen_values[num_electrons / 2];
	double mu_min = m_eigen_values[0];
	double mu_max = m_eigen_values[num_solution - 1];
	mu_max = (mu_max > 0.0) ? mu_max * 16.0 : fabs(mu_max) + 16.0;
	
	//必要最低限の軌道しか説いていない場合//
	if ((num_electrons + 1) / 2 >= num_solution) {
		mu_max = mu_max * 10.0;
	}

	int step = 0;
	double mu = (mu_min + mu_max) / 2.0;
	constexpr double THRESHOLD_MU = 1.0e-9;


	printf("mu = %f, %f, %f\n", mu, mu_min, mu_max); fflush(stdout);

	while (fabs(mu_max - mu_min) > THRESHOLD_MU){
		const double num_e = SetOccupancy(2.0, m_kbT, num_solution, m_eigen_values, mu, m_occupancy);

//#define DEBUG_PRINT1
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
		printf("occupancy[%d] = %f\n", s, m_occupancy[s]);	fflush(stdout);
	}
#endif


}

inline
void QUMASUN_MPI::mSetOccupancySpinOn() {


	const double d_num_electron = (double)num_electrons;
	//const double init_mu = m_eigen_values[num_electrons / 2];
	double mu_min = std::min(m_eigen_values[0], m_eigen_values_down[0]);
	double mu_max = std::max(m_eigen_values[num_solution - 1], m_eigen_values_down[num_solution - 1]);
	mu_max = (mu_max > 0.0) ? mu_max * 16.0 : fabs(mu_max) + 16.0;

	//必要最低限の軌道しか説いていない場合//
	if ((num_electrons + 1) / 2 >= num_solution) {
		mu_max = mu_max * 10.0;
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


