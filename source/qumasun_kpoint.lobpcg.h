#pragma once
//#define _USE_MATH_DEFINES
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <complex>
#include "qumasun_kpoint.h"
#include "lobpcg_z_multi_mpi.h"


//mpi supported//
inline
void QUMASUN_KPOINT::mSolveLOBPCG(int steps, int SCF_current_step) {
	const int proc_id = GetProcessID(m_mpi_comm);
	const bool is_ddm_root = IsRoot(m_ddm_comm);
	const size_t local_size = ml_grid.Size3D();

#ifdef USE_SCALAPACK
	BlacsGridInfo blacs_grid = BeginBLACS(m_ddm_comm, m_mpi_split_color);
#endif

	for(int sk = 0; sk < m_num_having_spin_kpoint;++sk){
		const int kpoint_x = ml_wave_set[sk].kpoint_x;
		const int kpoint_y = ml_wave_set[sk].kpoint_y;
		const int kpoint_z = ml_wave_set[sk].kpoint_z;
		auto& V_tot = (ml_wave_set[sk].spin == SPIN::UP) ? ml_Vtot : ml_Vtot_down;

		if (is_ddm_root) {
			printf("[%d] sk=%d, kpoint=%d,%d,%d, spin=%s\n", proc_id, sk, ml_wave_set[sk].kpoint_x, ml_wave_set[sk].kpoint_y, ml_wave_set[sk].kpoint_z, ml_wave_set[sk].spin == SPIN::UP ? "up" : "down");
		}

		LOBPCG::Eigen_z_multi_mpi(ml_grid, m_dx * m_dy * m_dz, num_solution, steps,
			ml_wave_set[sk].eigen_values, ml_wave_set[sk].l_psi_set, ml_keep_lobpcg[sk].keep_p, m_work,
			ml_keep_lobpcg[sk].keep_S_matrix, SCF_current_step, sk, 
#ifdef USE_SCALAPACK
			blacs_grid,
#endif
			[&V_tot, &kpoint_x, &kpoint_y, &kpoint_z, &sk, this](SoAComplex& Ax, const SoAComplex& x) {
				this->mHamiltonianMatrix_ddm(Ax, V_tot, x, kpoint_x, kpoint_y, kpoint_z, sk);
			},
			watch);
	}
	
#ifdef USE_SCALAPACK
	EndBLACS(blacs_grid);
#endif

    MPI_Barrier(m_mpi_comm);
    watch.Record(18);
}

