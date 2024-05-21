#pragma once
//#define _USE_MATH_DEFINES
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include "qumasun_td.h"
#include "vecmath.h"
#include "GridDifference2nd.h"
#include "GridDifference4th.h"
#include "GridDifference6th.h"
#include "GridDifference8th.h"
#include "GridDifferenceKpoint2nd.h"
#include "GridDifferenceKpoint6th.h"
#include "GridDifferenceKpoint8th.h"
#include "actual_kpoint.h"

#define KINETIC_8TH



/*
calculate qk = H' pk for K point
where H' = -(1/2) (\Delta +2ik \Nabla - k^2) + V
*/
inline
void QUMASUN_TD::mHamiltonianMatrix_ddm(SoAComplex& l_Hp, const double* l_Vtot, const SoAComplex& l_phi, int kpoint_x, int kpoint_y, int kpoint_z, int id_spin_kpoint) {
#define ZERO_CLEAR
#ifdef ZERO_CLEAR
	SoAC::SetZero(l_Hp, ml_grid.Size3D());
#endif
	
	//calculate q = V p , where V is effective potential////////////////
	if (m_hamiltonian_type == HAMILTONIAN::KohnSham_PP) {
		mPotentialMatrix_ddm(l_Hp.re, l_Vtot, l_phi.re);
		mPotentialMatrix_ddm(l_Hp.im, l_Vtot, l_phi.im);
		watch.Record( 12);
				
		m_pp_integrator.ProjectionPP(l_Hp, l_phi, ml_grid, m_nuclei, m_num_nuclei, id_spin_kpoint);
        MPI_Barrier(ml_grid.mpi_comm);
		watch.Record(13);
	} else {
		mPotentialMatrix_ddm(l_Hp.re, l_Vtot, l_phi.re);
		mPotentialMatrix_ddm(l_Hp.im, l_Vtot, l_phi.im);
		watch.Record( 12);
	}
	
	//calculate qk += K pk , where K is kinetic energy operator. //////////
	mKineticMatrixAdd_ddm_kpint(l_Hp, l_phi, kpoint_x, kpoint_y, kpoint_z);
	watch.Record( 11);
	///////complete to create qk = Hpk //
}


inline
void QUMASUN_TD::mKineticMatrixAdd_ddm(double* Kp, const double* p) {

#ifdef KINETIC_8TH
	Laplasian8th_ddm<double>(ml_grid, Kp, p, -1.0 / 2.0, m_dx, m_dy, m_dz);
#elif defined( KINETIC_6TH)
	Laplasian6th_ddm<double>(ml_grid, Kp, p, -1.0 / 2.0, m_dx, m_dy, m_dz);
#elif defined( KINETIC_4TH)
	Laplasian4th_ddm<double>(ml_grid, Kp, p, -1.0 / 2.0, m_dx, m_dy, m_dz);
#else
	Laplasian2nd_ddm<double>(ml_grid, Kp, p, -1.0 / 2.0, m_dx, m_dy, m_dz);
#endif
}

inline
void QUMASUN_TD::mKineticMatrixAdd_ddm_kpint(SoAComplex& Kp, const SoAComplex& p, int kpoint_x, int kpoint_y, int kpoint_z) {
	
	
	const double kpx = (double)kpoint_x;
	const double kpy = (double)kpoint_y;
	const double kpz = (double)kpoint_z;

	//int proc_id = GetProcessID(m_mpi_comm);
	//printf("[%d] call mKineticMatrixAdd_ddm_kpint: %d, %f, %f\n", proc_id, m_kpoint_sampling[0], m_dkx, m_box_x); fflush(stdout);

#ifdef KINETIC_8TH
	LaplasianKpoint8th_ddm(ml_grid, Kp, p, -1.0 / 2.0, m_dx, m_dy, m_dz, kpx * m_dkx, kpy * m_dky, kpz * m_dkz);
#elif defined( KINETIC_6TH)
	LaplasianKpoint6th_ddm(ml_grid, Kp, p, -1.0 / 2.0, m_dx, m_dy, m_dz, kpx * m_dkx, kpy * m_dky, kpz * m_dkz);
#elif defined( KINETIC_4TH)
	LaplasianKpoint4th_ddm(ml_grid, Kp, p, -1.0 / 2.0, m_dx, m_dy, m_dz, kpx * m_dkx, kpy * m_dky, kpz * m_dkz);
#else
	LaplasianKpoint2nd_ddm(ml_grid, Kp, p, -1.0 / 2.0, m_dx, m_dy, m_dz, kpx * m_dkx, kpy * m_dky, kpz * m_dkz);
#endif
}

// Vpr = V(r) * pr, operate V operator in real space//
inline
void QUMASUN_TD::mPotentialMatrix_ddm(double* Vp, const double* V, const double* p) {
	const int64_t size_3d = ml_grid.Size3D();
	for (int64_t i = 0; i < size_3d; ++i) {
		Vp[i] = V[i] * p[i];
	}
}
