#pragma once
//#define _USE_MATH_DEFINES
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include "qumasun_mpi.h"
#include "vecmath.h"
#include "GridDifference2nd.h"
#include "GridDifference4th.h"
#include "GridDifference6th.h"
#include "GridDifference8th.h"
#include "w_dgemm.h"

#define KINETIC_8TH
#define ZERO_CLEAR

/*
calculate qk = H pk
*/
inline
void QUMASUN_MPI::mHamiltonianMatrix_ddm(double* l_Hp, const double* l_Vtot, const double* l_phi) {

#ifdef ZERO_CLEAR
	vecmath::SetZero<double>(l_Hp, ml_grid.Size3D());
#endif
	
	//calculate q = V p , where V is effective potential////////////////
	if (m_hamiltonian_type == HAMILTONIAN::KohnSham_PP) {
		mPotentialMatrix_ddm(l_Hp, l_Vtot, l_phi);
		watch.Record( 12);

		m_pp_integrator.ProjectionPP(l_Hp, l_phi, ml_grid, m_nuclei, m_num_nuclei);
		watch.Record(13);
	} else {
		mPotentialMatrix_ddm(l_Hp, l_Vtot, l_phi);
		watch.Record(12);
	}
	
	//calculate qk += K pk , where K is kinetic energy operator. //////////
	mKineticMatrixAdd_ddm(l_Hp, l_phi);
	watch.Record(11);
	///////complete to create qk = Hpk //
}


inline
void QUMASUN_MPI::mHamiltonianMatrix_bundle_ddm(int num_states, double* l_Hp, const double* l_Vtot, const double* l_phi, double* l_tmp) {
	const int64_t local_size = ml_grid.Size3D();

#ifdef ZERO_CLEAR
	for (int n = 0; n < num_states; ++n) {
		vecmath::SetZero<double>(l_Hp + local_size * n, local_size);
	}
#endif

	//calculate q = V p , where V is effective potential////////////////
	
	if (m_hamiltonian_type == HAMILTONIAN::KohnSham_PP) {		
#if 1
		transpose<double>(local_size, num_states, l_phi, local_size, l_tmp, num_states);
		watch.Record(10);

		m_pp_integrator.ProjectionPP_bundle(num_states, l_Hp, l_tmp, ml_grid, m_nuclei, m_num_nuclei);
		watch.Record(13);

		transpose<double>(num_states, local_size, l_Hp, num_states, l_tmp, local_size);
		watch.Record(10);

		for (int n = 0; n < num_states; ++n) {
			mPotentialMatrix_ddm(l_Hp + local_size * n, l_Vtot, l_phi + local_size * n);
		}
		watch.Record(12);
		for (int64_t i = 0; i < (int64_t)num_states* (int64_t)local_size; ++i) {
			l_Hp[i] += l_tmp[i];
		}
		watch.Record(10);

#else
		for (int n = 0; n < num_states; ++n) {
			mPotentialMatrix_ddm(l_Hp + local_size * n, l_Vtot, l_phi + local_size * n);
		}
		watch.Record( 12);
		for (int n = 0; n < num_states; ++n) {
			m_pp_integrator.ProjectionPP(l_Hp + local_size * n, l_phi + local_size * n, ml_grid, m_nuclei, m_num_nuclei);
		}
		watch.Record( 13);
#endif


	} else {
		for (int n = 0; n < num_states; ++n) {
			mPotentialMatrix_ddm(l_Hp + local_size * n, l_Vtot, l_phi + local_size * n);
		}
		watch.Record( 12);

	}

	//calculate qk += K pk , where K is kinetic energy operator. //////////
	for (int n = 0; n < num_states; ++n) {
		mKineticMatrixAdd_ddm(l_Hp + local_size * n, l_phi + local_size * n);
	}
	watch.Record( 11);
	///////complete to create qk = Hpk //
}


#if 0
void QUMASUN_MPI::mKineticMatrixAdd(RspaceFunc<double>& Kp, const RspaceFunc<double>& p) {

#ifdef KINETIC_6TH
	Laplasian6th<double>(GridRange{ 0,0,0,m_size_x, m_size_y, m_size_z },
		Kp.Pointer(), p.Pointer(), -1.0 / 2.0, m_dx, m_dy, m_dz);

#elif defined( KINETIC_4TH)
	Laplasian4th<double>(GridRange{ 0,0,0,m_size_x, m_size_y, m_size_z },
			Kp.Pointer(), p.Pointer(), -1.0 / 2.0, m_dx, m_dy, m_dz);
#else
	Laplasian2nd<double>({ 0,0,0,m_size_x, m_size_y, m_size_z },
		Kp.Pointer(), p.Pointer(), -1.0 / 2.0, m_dx, m_dy, m_dz);

#endif
}
#endif

inline
void QUMASUN_MPI::mKineticMatrixAdd_ddm(double* Kp, const double* p) {

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

// Vpr = V(r) * pr, operate V operator in real space//
inline
void QUMASUN_MPI::mPotentialMatrix(RspaceFunc<double>& Vp, const RspaceFunc<double>& V, const RspaceFunc<double>& p) {
	
	
	for (int iz = 0; iz < m_size_z; ++iz) {
		for (int iy = 0; iy < m_size_y; ++iy) {
			for (int ix = 0; ix < m_size_x; ++ix) {
				const size_t i = ix + m_size_x * (iy + (m_size_y * iz));

				Vp[i] = V[i] * p[i];

			}
		}
	}

}

// Vpr = V(r) * pr, operate V operator in real space//
inline
void QUMASUN_MPI::mPotentialMatrix_ddm(double* Vp, const double* V, const double* p) {
	const int64_t size_3d = ml_grid.Size3D();
	for (int64_t i = 0; i < size_3d; ++i) {
		Vp[i] = V[i] * p[i];
	}
}
