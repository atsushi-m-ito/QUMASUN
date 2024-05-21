#pragma once
//#define _USE_MATH_DEFINES
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <complex>
#include "qumasun_mpi.h"
#include "lobpcg_d_multi_mpi.h"

//mpi supported//
inline
void QUMASUN_MPI::mSolveLOBPCG(int steps, int SCF_current_step) {
	const int proc_id = GetProcessID(m_mpi_comm);
	//printf("[%d]mSolveLOBPCG begin\n", proc_id); fflush(stdout);

#ifdef USE_SCALAPACK
	BlacsGridInfo blacs_grid = BeginBLACS(m_ddm_comm, m_mpi_split_color);
#endif

	{		

		LOBPCG::Eigen_d_multi_mpi(ml_grid, m_dx * m_dy * m_dz, num_solution, steps,
			m_eigen_values, ml_psi[0], ml_lobpcg_keep_p, m_work, 
			m_lobpcg_keep_S_matrix, SCF_current_step,
#ifdef USE_SCALAPACK
			blacs_grid,
#endif
#ifdef OPERATION_BUNDLE
			[this](double* Ax, double* x, double* l_tmp) {
				this->mHamiltonianMatrix_bundle_ddm(num_solution, Ax, this->ml_Vtot, x, l_tmp);
			},
#else
			[this](double* Ax, double* x) {
				this->mHamiltonianMatrix_ddm(Ax, this->ml_Vtot, x);
			},
#endif
			watch);
	}
	if(m_num_spin==2){
		const size_t local_size = ml_grid.Size3D();
		LOBPCG::Eigen_d_multi_mpi(ml_grid, m_dx * m_dy * m_dz, num_solution, steps,
			m_eigen_values_down, ml_psi[num_solution], ml_lobpcg_keep_p + local_size* num_solution, m_work,
			m_lobpcg_keep_S_matrix + num_solution* num_solution*3, SCF_current_step,
#ifdef USE_SCALAPACK
			blacs_grid,
#endif
#ifdef OPERATION_BUNDLE
			[this](double* Ax, double* x, double* l_tmp) {
				this->mHamiltonianMatrix_bundle_ddm(num_solution, Ax, this->ml_Vtot_down, x, l_tmp);
			},
#else
			[this](double* Ax, double* x) {
				this->mHamiltonianMatrix_ddm(Ax, this->ml_Vtot_down, x);
			}, 
#endif
			watch);
	}
	
	if(is_spin_on && (m_num_spin==1)) { //spin parallel//
		//printf("[%d]mSolveLOBPCG spin-parallel\n", proc_id); fflush(stdout);
		//const int proc_id = GetProcessID(m_mpi_comm);
		const int num_procs = GetNumProcess(m_mpi_comm);
		const int TAG = 0x904;
		if (proc_id == 0) {
			MPI_Status status;
			MPI_Recv(m_eigen_values_down, num_solution, MPI_DOUBLE, num_procs / 2, TAG, m_mpi_comm, &status);
		} else if (proc_id == num_procs / 2) {
			MPI_Send(m_eigen_values, num_solution, MPI_DOUBLE, 0, TAG, m_mpi_comm);
		}
	}
	//printf("[%d]mSolveLOBPCG fin\n", proc_id); fflush(stdout);


#ifdef USE_SCALAPACK
	EndBLACS(blacs_grid);
#endif

    MPI_Barrier(m_mpi_comm);
    watch.Record(18);
}

