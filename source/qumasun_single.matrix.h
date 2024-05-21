#pragma once
//#define _USE_MATH_DEFINES
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include "qumasun_single.h"
#include "wrap_fft.h"
#include "vecmath.h"
#include "GridDifference2nd.h"
#include "GridDifference4th.h"
#include "GridDifference6th.h"
#include "GridDifference8th.h"

#define KINETIC_8TH

/*
calculate qk = H pk
*/
inline
void QUMASUN_SINGLE::mHamiltonianMatrix(RspaceFunc<double>& Hp, const RspaceFunc<double>& V, const RspaceFunc<double>& p) {
	
#define ZERO_CLEAR
#ifdef ZERO_CLEAR
	vecmath::SetZero<double>(Hp, m_size_3d);
#endif

	//calculate q = V p , where V is effective potential////////////////
	if (m_hamiltonian_type == HAMILTONIAN::KohnSham_PP) {
		mPotentialMatrix(Hp, V, p);
		watch.Record(12);
		m_pp_integrator.ProjectionPP(Hp, p, m_nuclei, m_num_nuclei);
		watch.Record(13);
	} else {
		mPotentialMatrix(Hp, V, p);
	}
	//calculate qk += K pk , where K is kinetic energy operator. //////////
	mKineticMatrixAdd(Hp, p);
	watch.Record( 11);
	///////complete to create qk = Hpk //
}

inline
void QUMASUN_SINGLE::mKineticMatrixAdd(RspaceFunc<double>& Kp, const RspaceFunc<double>& p) {

#ifdef KINETIC_6TH
	Laplasian6th<double>({ 0,0,0,m_size_x, m_size_y, m_size_z },
		Kp.Pointer(), p.Pointer(), -1.0 / 2.0, m_dx, m_dy, m_dz);

#elif defined(KINETIC_4TH)
	Laplasian4th<double>({ 0,0,0,m_size_x, m_size_y, m_size_z },
			Kp.Pointer(), p.Pointer(), -1.0 / 2.0, m_dx, m_dy, m_dz);
	
#elif defined( KINETIC_4TH)

	for (int iz = 0; iz < m_size_z; ++iz) {
		const int idz_m3 = ((iz <= 2) ? m_size_z - 3 : -3) * m_size_x * m_size_y;
		const int idz_m2 = ((iz <= 1) ? m_size_z - 2 : -2) * m_size_x * m_size_y;
		const int idz_m = ((iz == 0) ? m_size_z - 1 : -1) * m_size_x * m_size_y;
		const int idz_p = ((iz == m_size_z - 1) ? 1 - m_size_z : 1) * m_size_x * m_size_y;
		const int idz_p2 = ((iz >= m_size_z - 2) ? 2 - m_size_z : 2) * m_size_x * m_size_y;
		const int idz_p3 = ((iz >= m_size_z - 3) ? 3 - m_size_z : 3) * m_size_x * m_size_y;

		for (int iy = 0; iy < m_size_y; ++iy) {
			const int idy_m3 = ((iy <= 2) ? m_size_y - 3 : -3) * m_size_x;
			const int idy_m2 = ((iy <= 1) ? m_size_y - 2 : -2) * m_size_x;
			const int idy_m = ((iy == 0) ? m_size_y - 1 : -1) * m_size_x;
			const int idy_p = ((iy == m_size_y - 1) ? 1 - m_size_y : 1) * m_size_x;
			const int idy_p2 = ((iy >= m_size_y - 2) ? 2 - m_size_y : 2) * m_size_x;
			const int idy_p3 = ((iy >= m_size_y - 3) ? 3 - m_size_y : 3) * m_size_x;
			for (int ix = 0; ix < m_size_x; ++ix) {
				const int idx_m3 = ((ix <= 2) ? m_size_x - 3 : -3);
				const int idx_m2 = ((ix <= 1) ? m_size_x - 2 : -2);
				const int idx_m = ((ix == 0) ? m_size_x - 1 : -1);
				const int idx_p = ((ix == m_size_x - 1) ? 1 - m_size_x : 1);
				const int idx_p2 = ((ix >= m_size_x - 2) ? 2 - m_size_x : 2);
				const int idx_p3 = ((ix >= m_size_x - 3) ? 3 - m_size_x : 3);
				const size_t i = ix + m_size_x * (iy + (m_size_y * iz));

				double d2psidx2 = ((p[i + idx_p3] + p[i + idx_m3]) -54.0 * (p[i + idx_p2] + p[i + idx_m2]) +783.0*(p[i + idx_p] + p[i + idx_m]) - 1460.0 * p[i]) / (m_dx * m_dx);
				d2psidx2 += ((p[i + idy_p3] + p[i + idy_m3]) - 54.0*(p[i + idy_p2] + p[i + idy_m2]) + 783.0 * (p[i + idy_p] + p[i + idy_m]) - 1460.0 * p[i]) / (m_dy * m_dy);
				d2psidx2 += ((p[i + idz_p3] + p[i + idz_m3]) -54.0*(p[i + idz_p2] + p[i + idz_m2]) +783.0*(p[i + idz_p] + p[i + idz_m]) - 1460.0 * p[i]) / (m_dz * m_dz);
				Kp[i] += -d2psidx2 / (2.0 * 24.0* 24.0);

			}
		}
	}

#elif 1
	Laplasian2nd<double>({ 0,0,0,m_size_x, m_size_y, m_size_z },
		Kp.Pointer(), p.Pointer(), -1.0 / 2.0, m_dx, m_dy, m_dz);
#else
	for (int iz = 0; iz < m_size_z; ++iz) {
		const int idz_m = ((iz == 0) ? m_size_z - 1 : -1) * m_size_x * m_size_y;
		const int idz_p = ((iz == m_size_z - 1) ? 1 - m_size_z : 1) * m_size_x * m_size_y;
		for (int iy = 0; iy < m_size_y; ++iy) {
			const int idy_m = ((iy == 0) ? m_size_y - 1 : -1) * m_size_x;
			const int idy_p = ((iy == m_size_y - 1) ? 1 - m_size_y : 1) * m_size_x;
			for (int ix = 0; ix < m_size_x; ++ix) {
				const int idx_m = ((ix == 0) ? m_size_x - 1 : -1);
				const int idx_p = ((ix == m_size_x - 1) ? 1 - m_size_x : 1);
				const size_t i = ix + m_size_x * (iy + (m_size_y * iz));

				const double d2psidx2 = (p[i + idx_p] + p[i + idx_m] - 2.0 * p[i]) / (m_dx*m_dx)
					+ (p[i + idy_p] + p[i + idy_m] - 2.0 * p[i]) / (m_dy*m_dy)
					+ (p[i + idz_p] + p[i + idz_m] - 2.0 * p[i]) / (m_dz*m_dz);
				Kp[i] += -d2psidx2 / 2.0;

			}
		}
	}
#endif
}

/*
Vp = V(r) * p
*/
inline
void QUMASUN_SINGLE::mPotentialMatrix(RspaceFunc<double>& Vp, const RspaceFunc<double>& V, const RspaceFunc<double>& p) {
	
	// Vpr = V(r) * pr, operate V operator in real space//
	for (int iz = 0; iz < m_size_z; ++iz) {
		for (int iy = 0; iy < m_size_y; ++iy) {
			for (int ix = 0; ix < m_size_x; ++ix) {
				const size_t i = ix + m_size_x * (iy + (m_size_y * iz));

				Vp[i] = V[i] * p[i];

			}
		}
	}

}
