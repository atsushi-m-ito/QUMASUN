#pragma once
//#define _USE_MATH_DEFINES
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <complex>
#include "qumasun_single.h"
#include "wrap_fft.h"
#include "lobpcg_d_multi.h"

inline
void QUMASUN_SINGLE::mSolveLOBPCG(int steps, bool is_first_step) {

	{
		RspaceFunc<double>& V = m_Vtot;

		EigenLOBPCG_d_multi(m_size_3d, m_dx * m_dy * m_dz, num_solution, steps,
			m_eigen_values, m_psi_set[0].Pointer(), m_lobpcg_keep_p, m_work, is_first_step,
			[this, &V](double* Ax, double* x) {
				RspaceFunc<double> refAx(Ax);
				RspaceFunc<double> refx(x);
				this->mHamiltonianMatrix(refAx, V, refx);
			},watch);
	}

	if (is_spin_on) {
		RspaceFunc<double>& V = m_Vtot_down;
		EigenLOBPCG_d_multi(m_size_3d, m_dx * m_dy * m_dz, num_solution, steps,
			m_eigen_values_down, m_psi_set[num_solution].Pointer(), 
			m_lobpcg_keep_p + m_size_3d * num_solution, m_work, is_first_step,
			[this, &V](double* Ax, double* x) {
				RspaceFunc<double> refAx(Ax);
				RspaceFunc<double> refx(x);
				this->mHamiltonianMatrix(refAx, V, refx);
			},watch);
	}
	//mCheckEigenVector();
}

