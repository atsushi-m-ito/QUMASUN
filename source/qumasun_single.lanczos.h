#pragma once
//#define _USE_MATH_DEFINES
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include "qumasun_single.h"
#include "lanczos_z.h"
#include "msz.complex.h"
#include "vecmath.h"

//Solve the Kohn-Sham equation by Lanczos method////////////////

inline
void QUMASUN_SINGLE::mSolveLanczos() {

	RspaceFunc<double>& V = m_Vtot;
	

	//const int num_solusion_Lanczos = std::max<int>(num_solution * 2, (int)sqrt(m_size_3d));
	const int num_solusion_Lanczos = std::max<int>(num_solution * 2, (int)cbrt(m_size_3d));
	

#if 1
	printf("Lanczos solver is not yet supported in Real field(not complex field)\n");
#else
	EigenLanczos_z(m_size_3d, num_solution, num_solusion_Lanczos,
		m_eigen_values, m_psi_buffer.Pointer(),
		[this, &V](dcomplex* Ax, dcomplex* x) {
			RspaceFunc<dcomplex> refAx(Ax);
			RspaceFunc<dcomplex> refx(x);
			this->mHamiltonianMatrix(refAx, V, refx);
		});
#endif

	mCheckEigenVector();

}

inline
void QUMASUN_SINGLE::mCheckEigenVector() {
	using namespace msz::complex;
	RspaceFunc<double>& V = m_Vtot;

	for (int s = 0; s < num_solution; ++s) {
		RspaceFunc<double> & p = m_psi_set[s];
		RspaceFunc<double> Hp(m_work);

		mHamiltonianMatrix(Hp, V, p);

		{
			const double norm_Hp = vecmath::Norm<double>(Hp, m_size_3d);
			const double norm_p = vecmath::Norm<double>(p, m_size_3d);
			const double in = vecmath::InnerProd<double>(p, Hp, m_size_3d);

			printf("(%d) norm(p), norm(Hp), cos = %f, %f, (%f)\n", s, norm_p, norm_Hp, in / sqrt(norm_Hp*norm_p));

		}

	}

	
}
