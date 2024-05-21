#pragma once
#ifdef USE_MPI
#include <mpi.h>
#include <cmath>
#include "mpi_helper.h"
#include "GridRange.h"
#include "soacomplex.h"

inline void GramSchmidt_ddm(GridRangeMPI& l_grid, int num_solution, double** psi_set, double dVol) {
	double* l_inner = new double[num_solution*2];
	double* sum_inner = l_inner + num_solution;

	const int local_size = l_grid.Size3D();
	for (int n = 0; n < num_solution; ++n){
		double* psi_n = psi_set[n];

		//diagonalize//
		for (int k = 0; k < n; ++k) {
			const double* psi_k = psi_set[k];

			double in = 0.0;
			for (int i = 0; i < local_size; ++i) {
				in += psi_n[i] * psi_k[i];
			}
			l_inner[k] = in * dVol;
		}
		MPI_Allreduce(l_inner, sum_inner, n, MPI_DOUBLE, MPI_SUM, l_grid.mpi_comm);

		for (int k = 0; k < n; ++k) {
			const double in = sum_inner[k];
			const double* psi_k = psi_set[k];
			for (int i = 0; i < local_size; ++i) {
				psi_n[i] -= in * psi_k[i];
			}
		}

		//normalize//
		double norm = 0.0;
		for (int i = 0; i < local_size; ++i) {
			norm += psi_n[i] * psi_n[i];
		}
		norm *= dVol;
		double sum_norm = 0.0;
		MPI_Allreduce(&norm, &sum_norm, 1, MPI_DOUBLE, MPI_SUM, l_grid.mpi_comm);

		const double inorm = 1.0 / std::sqrt(sum_norm);
		for (int i = 0; i < local_size; ++i) {
			psi_n[i] *= inorm;
		}

	}

	delete[] l_inner;
}

inline void GramSchmidt_ddm(GridRangeMPI& l_grid, int num_solution, SoAComplex* psi_set, double dVol) {
	OneComplex* l_inner = new OneComplex[num_solution * 2];
	OneComplex* sum_inner = l_inner + num_solution;

	const int local_size = l_grid.Size3D();
	for (int n = 0; n < num_solution; ++n) {
		auto& psi_n = psi_set[n];

		//diagonalize//
		for (int k = 0; k < n; ++k) {
			const auto& psi_k = psi_set[k];

			OneComplex in = SoAC::InnerProd(psi_k, psi_n, local_size);
			l_inner[k].r = in.r * dVol;
			l_inner[k].i = in.i * dVol;
		}
		MPI_Allreduce(l_inner, sum_inner, n, MPI_DOUBLE_COMPLEX, MPI_SUM, l_grid.mpi_comm);

		for (int k = 0; k < n; ++k) {
			const OneComplex in = sum_inner[k];
			const auto& psi_k = psi_set[k];
			for (int i = 0; i < local_size; ++i) {
				psi_n.re[i] -= in.r * psi_k.re[i] - in.i * psi_k.im[i];
				psi_n.im[i] -= in.r * psi_k.im[i] + in.i * psi_k.re[i];
			}
		}

		//normalize//
		double norm = SoAC::Norm(psi_n, local_size);
		norm *= dVol;
		double sum_norm = 0.0;
		MPI_Allreduce(&norm, &sum_norm, 1, MPI_DOUBLE, MPI_SUM, l_grid.mpi_comm);

		const double inorm = 1.0 / std::sqrt(sum_norm);
		SoAC::MulC(psi_n, inorm, local_size);

	}

	delete[] l_inner;
}
#endif
