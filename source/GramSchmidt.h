#pragma once

#include <cmath>
#include "wave_function.h"

inline void GramSchmidt(int N, int num_solution, RspaceFunc<double>* m_psi_set, double dVol) {

	for (int n = 0; n < num_solution; ++n){
		double* psi_n = m_psi_set[n];

		//diagonalize//
		for (int k = 0; k < n; ++k) {
			const double* psi_k = m_psi_set[k];

			double in = 0.0;
			for (int i = 0; i < N; ++i) {
				in += psi_n[i] * psi_k[i];
			}
			in *= dVol;

			for (int i = 0; i < N; ++i) {
				psi_n[i] -= in * psi_k[i];
			}
		}

		//normalize//
		double norm = 0.0;
		for (int i = 0; i < N; ++i) {
			norm += psi_n[i] * psi_n[i];
		}
		norm *= dVol;

		const double inorm = 1.0 / std::sqrt(norm);
		for (int i = 0; i < N; ++i) {
			psi_n[i] *= inorm;
		}

	}

}
