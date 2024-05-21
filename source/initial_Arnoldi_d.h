#pragma once
#include <random>
#include "w_dstevd.h"
#include "GridRange.h"
#include "vecmath.h"


/*
* Lanczos method for Eigen value problem
* 
* The size_trimatrix can be set a smaller value than size of matrix N,
* and while for accuracy the size_trimatrix should be larger than the num solution(num of required solution).
* 
*/
template <typename OperationA>
inline
int InitializeByArnoldy_d(const GridRangeMPI& l_grid, const double dVol, const int num_solution, const int num_trial,
	double* state_vectors,
	OperationA OpeA) {

	using namespace vecmath;

	constexpr double SMALL_LIMIT = 1.0e-6;

	const int proc_id = GetProcessID(l_grid.mpi_comm);
	const bool is_root = (proc_id == 0);
	const size_t local_size = l_grid.Size3D();

	auto IntervalPointers = [](double* a, size_t N, size_t num_solution) {
		std::vector<double*> pointers(num_solution);
		for (size_t i = 0; i < num_solution; ++i) {
			pointers[i] = a + i * N;
		}
		return pointers;
		};


	auto x = IntervalPointers(state_vectors, local_size, num_solution);
		
	





	int gen_count = 1;
	double beta = 0.0;
	
	//Loop of Lanczos to create symmetric tridiagonal matrix T// 
	for (int t = 0; t < num_trial - 1; ++t) {

		//p is state vector in wave space////////////////////////
		const int t_next = (t + 1) % num_solution;
		double* p = x[t % num_solution];
		double* Hp = x[t_next];

		//create H p//
		OpeA(Hp, p);
		++gen_count;
		if (t > 1) {
			int t_prev = (t - 1 + num_solution) % num_solution;
			double* p_prev = x[t_prev];
			for (size_t i = 0; i < local_size; ++i) {
				Hp[i] -= beta * p_prev[i];
			}
		}

		double value[3]{ 0.0,0.0,0.0 };
		double& l_norm = value[0];
		double& l_alpha = value[1];
		double& test = value[2];
		for (size_t i = 0; i < local_size; ++i) {
			l_norm += Hp[i] * Hp[i];
			l_alpha += p[i] * Hp[i];
			test += p[i] * p[i];
		}

		double sum_value[3];
		MPI_Allreduce(value, sum_value, 3, MPI_DOUBLE, MPI_SUM, l_grid.mpi_comm);
		double alpha = sum_value[1] * dVol;
		double &l_beta = value[0] = 0.0;
		double& l_rp = value[1] = 0.0;
		sum_value[0] = 1.0 / std::sqrt(sum_value[0] * dVol);
		for (size_t i = 0; i < local_size; ++i) {

			auto r = Hp[i] - alpha * p[i];
			l_beta += r * r;
			//Hp[i] *= sum_value[0];
			Hp[i] = r;
			l_rp += r * p[i];
		}
		
		MPI_Allreduce(value, sum_value, 2, MPI_DOUBLE, MPI_SUM, l_grid.mpi_comm);
		beta = std::sqrt(sum_value[0] * dVol);
		double coef = 1.0 / beta ;
		for (size_t i = 0; i < local_size; ++i) {
			Hp[i] *= coef;
		}

		if (beta < SMALL_LIMIT) {

			break;
		}

	}
	////////////////////////////End of loop of Lanczos to create Symmetric Tridiagonal T// 

	return gen_count;
}

