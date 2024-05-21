#pragma once
//#include <complex>
#include <random>
//#include "wrap_lapack.h"
#include "w_dstevd.h"
#include "dSFMT_wrapper.h"
#include "msz.complex.h"



/*
* Lanczos method for Eigen value problem
* 
* The size_trimatrix can be set a smaller value than size of matrix N,
* and while for accuracy the size_trimatrix should be larger than the num solution(num of required solution).
* 
*/
template <typename OperationA>
inline
int EigenLanczos_d(const int N, const int num_solution, const int size_trimatrix,
	double* eigen_values, double* eigen_vectors,
	OperationA OpeA) {

//	using dcomplex = std::complex<double>;
//	using namespace msz::complex;

	dcomplex* p_all = new dcomplex[N * size_trimatrix];
	dcomplex* Hp = new dcomplex[N];
	double* diagonal = new double[size_trimatrix];
	double* off_diagonal = new double[size_trimatrix];
	double* eigen_values_tri = new double[size_trimatrix];
	double* eigen_vectors_tri = new double[size_trimatrix*size_trimatrix];

	int valid_size_matrix = size_trimatrix;
	
	const uint32_t seed = 123456789;
	//WrapdSFMT dsfmt(seed, std::min<int>(N, 10000));	//状態保持.
	std::mt19937 mt(seed);
	std::normal_distribution<> distribution(0.0, 1.0);
	
	//initialize//
	{

		dcomplex* p = p_all;
		for (size_t i = 0; i < N; ++i) {
			p[i] = distribution(mt);
		}

		//create H p//
		OpeA(Hp, p);

		for (size_t i = 0; i < N; ++i) {
			p[i] = Hp[i];
		}
		double cc = 0.0;
		for (size_t i = 0; i < N; ++i) {
			cc += norm(p[i]);
		}
		cc = 1.0 / sqrt(cc);
		for (size_t i = 0; i < N; ++i) {
			p[i] *= cc;
		}
	}


	double beta_prev = 0.0;
	//Loop of Lanczos to create symmetric tridiagonal matrix T// 
	for (int t = 0; t < size_trimatrix; ++t) {

		//p is state vector in wave space////////////////////////
		dcomplex* p = p_all + N * t;


		//create H p//
		OpeA(Hp, p);

		if (t > 0) {
			dcomplex* p_prev = p_all + N * (t - 1);
			for (size_t i = 0; i < N; ++i) {
				Hp[i] -= beta_prev * p_prev[i];
			}
		}

		double alpha = 0.0;
		for (size_t i = 0; i < N; ++i) {
			alpha += inner(p[i], Hp[i]).real();
		}
		diagonal[t] = alpha;

		if (t == size_trimatrix - 1) break;

		// off diagonal term beta //
		dcomplex* p_next = p_all + N * (t + 1);
		double beta = 0.0;
		//if (t == 0) {
		// off diagonal term beta = fabs(Ap - alpha * p) //
		for (size_t i = 0; i < N; ++i) {
			p_next[i] = Hp[i] - alpha * p[i];
			beta += norm(p_next[i]);
		}
		/*
		} else {
		// off diagonal term beta = fabs(Ap - alpha * p - beta[t-1] * p[t-1]) //
		double* p_prev = p_all + N * (t - 1);
		for (size_t i = 0; i < N; ++i) {
		p_next[i] = Hp[i] - alpha * p[i] - beta_prev * p_prev[i];
		beta += p_next[i] * p_next[i];
		}
		}
		*/
		bool checking = false;
		constexpr double SMALL_LIMIT = 1.0e-10;
		if (beta > SMALL_LIMIT) {
			beta = sqrt(beta);

			// next vector p[t+1] = (Ap - alpha * p)/beta //
			const double ibeta = 1.0 / beta;
			for (size_t i = 0; i < N; ++i) {
				p_next[i] *= ibeta;
			}
		} else {

			valid_size_matrix = t + 1;
			break;

			beta = 0.0;
			//過去の全ての解と直交するベクトルを再生成//
			for (size_t i = 0; i < N; ++i) {
				p_next[i] = distribution(mt());
			}
			checking = true;
		}

		/*
		if (t % 10 == 9) {
		checking = true;
		}
		*/
		checking = true;

		/*
		if(	checking ){
		for (int j = 0; j <= t; ++j) {
		double* p_old = p_all + N* j;
		double c = 0.0;
		for (size_t i = 0; i < N; ++i) {
		c += p_next[i] * p_old[i];
		}
		for (size_t i = 0; i < N; ++i) {
		p_next[i] -= c * p_old[i];
		}
		}
		double cc = 0.0;
		for (size_t i = 0; i < N; ++i) {
		cc += p_next[i] * p_next[i];
		}
		beta = sqrt(cc);
		cc = 1.0 / beta;
		for (size_t i = 0; i < N; ++i) {
		p_next[i] *= cc;
		}
		}
		*/

		off_diagonal[t] = beta;
		beta_prev = beta;

	}
	////////////////////////////End of loop of Lanczos to create Symmetric Tridiagonal T// 

	//Solve the eigen value problem for T ///////////////////////////////////////
	off_diagonal[size_trimatrix - 1] = 0.0;
	//対象三重対角行列の固有値と固有ベクトルを求める//
	lapack_DSTEVD(diagonal, off_diagonal, eigen_values_tri, eigen_vectors_tri, valid_size_matrix);

	//Create Eigen vector of A from that of T//

	for (int i = 0; i < num_solution; ++i) {
		eigen_values[i] = eigen_values_tri[i];
	}

	for (int i = 0; i < num_solution; ++i) {
		dcomplex* ev = eigen_vectors + N * i;
		double* ev_tri = eigen_vectors_tri + valid_size_matrix * i;

		for (int k = 0; k < N; ++k) {
			ev[k] = { 0.0, 0.0};
		}
		for (int t = 0; t < valid_size_matrix; ++t) {
			dcomplex* p = p_all + N * t;
			for (int k = 0; k < N; ++k) {
				ev[k] += p[k] * ev_tri[t];
			}
		}
		double cc = 0.0;
		for (int k = 0; k < N; ++k) {
			cc += norm(ev[k]);
		}
		cc = 1.0 / sqrt(cc);
		for (int k = 0; k < N; ++k) {
			ev[k] *= cc;
		}

	}

	//finish///////////////
	delete[] p_all;
	delete[] Hp;
	delete[] diagonal;
	delete[] off_diagonal;
	delete[] eigen_values_tri;
	delete[] eigen_vectors_tri;

	return std::min(valid_size_matrix, num_solution);
}
