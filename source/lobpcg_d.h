#pragma once
#include "wrap_lapack.h"
#include "w_dsygvx.h"

#include "vecmath.h"


//対角化を手動でせずに、小行列の一般化固有値解法に任せた方が良い//
//なぜなら、(1)matrix Aの演算が1回/iterationで済む。//

//※非一般化の固有値問題を解き、かつ、r,pのnormalizeを忘れると収束しなくなる//



/*
固有値問題を解くLOBPCG

*/
template <typename OperationA>
inline
void EigenLOBPCG_d(const int N, const double dVol, const int iter_max,
	double* eigen_values, double* eigen_vectors, double* work, double* keep, bool is_first_time,
	OperationA OpeA, int log_print_step = 0 ) {

	using namespace vecmath;

	

	auto Norm = [](const double& a) {return a * a; };

#if 1
	double* x = eigen_vectors;
	double* r = work;
	double* p = keep;
	double* Ax = work + N ;
	double* Ar = work + N * 2;
	double* Ap = work + N * 3;
#else
	double* x = new double[N * 3];
	double* r = x + N;
	double* p = x + N * 2;
	double* Ax = new double[N * 3];
	double* Ar = Ax + N;
	double* Ap = Ax + N * 2;
#endif

	double ei;

	double Sa[9];
	double Sb[9];
	double triangle_factors[9];
	double S_e_value[3];
	double S_e_vector[3 * 3];

	const double limit = 1.0e-10;



	{
		//Copy(x, eigen_vectors, N);


	}

	//calculte initial redidual r//
	{
		/*
		for (size_t i = 0; i < N; ++i) {
		double sum = 0.0;
		for (size_t j = 0; j < N; ++j) {
		sum += A[i * N + j] * x[j];
		}
		Ax[i] = sum;
		}*/
		Normalize(x, N, dVol);

		//Because Hamiltonian operator A was changed by changing electronic density rho, //
		// Ax and Ap should be re-calculated//
		OpeA(Ax, x); // calculate Ax from x 
		OpeA(Ap, p);

		ei = InnerProd(x, Ax, N) * dVol;
		//x is already normalize//
		
		double norm_r = 0.0;
		for (size_t i = 0; i < N; ++i) {
			r[i] = (Ax[i] - ei * x[i]);
			norm_r += Norm(r[i]);	//inner(r[i] , r[i]).real();
		}
		norm_r *= dVol;
		const double cc = 1.0 / sqrt(norm_r);		
		MulC(r, cc, N);

		if (is_first_time) {
			SetZero(p, N);
			SetZero(Ap, N);
		}
		
		printf("(0) ei, norm_r: %f, %f\n", ei, norm_r);
	}



//Loop of CG
	int iter;
	for (iter = 0; iter < iter_max; ++iter) {
			
		OpeA(Ar, r);
		//OpeA(Ax, x);

		for (int i = 0; i < 9; ++i) {
			Sa[i] = 0.0;
			Sb[i] = 0.0;
		}

#if 1
		Sa[0] = InnerProd(x, Ax, N) * dVol;
		Sa[1] = InnerProd(x, Ar, N) * dVol;
		Sa[2] = InnerProd(x, Ap, N) * dVol;

		Sa[3] = InnerProd(r, Ax, N) * dVol;
		Sa[4] = InnerProd(r, Ar, N) * dVol;
		Sa[5] = InnerProd(r, Ap, N) * dVol;

		Sa[6] = InnerProd(p, Ax, N) * dVol;
		Sa[7] = InnerProd(p, Ar, N) * dVol;
		Sa[8] = InnerProd(p, Ap, N) * dVol;
#if 0
		//対称化(数値誤差により非対称になっているため)
		Sa[1] = (Sa[0 * 3 + 1] + Sa[1 * 3 + 0]) / 2.0;
		Sa[3] = Sa[1];

		Sa[2] = (Sa[0 * 3 + 2] + Sa[2 * 3 + 0]) / 2.0;
		Sa[6] = Sa[2];

		Sa[5] = (Sa[1 * 3 + 2] + Sa[2 * 3 + 1]) / 2.0;
		Sa[7] = Sa[5];


#endif


		Sb[0] = InnerProd(x, x, N) * dVol;
		Sb[1] = InnerProd(x, r, N) * dVol;
		Sb[2] = InnerProd(x, p, N) * dVol;

		Sb[3] = InnerProd(r, x, N) * dVol;
		Sb[4] = InnerProd(r, r, N) * dVol;
		Sb[5] = InnerProd(r, p, N) * dVol;

		Sb[6] = InnerProd(p, x, N) * dVol;
		Sb[7] = InnerProd(p, r, N) * dVol;
		Sb[8] = InnerProd(p, p, N) * dVol;

#else
		for (size_t i = 0; i < N; ++i) {
			Sa[0] += inner(x[i] , Ax[i]);
			Sa[1] += inner(x[i] , Ar[i]);
			Sa[2] += inner(x[i] , Ap[i]);

			Sa[3] += inner(r[i] , Ax[i]);
			Sa[4] += inner(r[i] , Ar[i]);
			Sa[5] += inner(r[i] , Ap[i]);

			Sa[6] += inner(p[i] , Ax[i]);
			Sa[7] += inner(p[i] , Ar[i]);
			Sa[8] += inner(p[i] , Ap[i]);

			Sb[0] += inner(x[i] , x[i]);
			Sb[1] += inner(x[i] , r[i]);
			Sb[2] += inner(x[i] , p[i]);

			Sb[3] += inner(r[i] , x[i]);
			Sb[4] += inner(r[i] , r[i]);
			Sb[5] += inner(r[i] , p[i]);

			Sb[6] += inner(p[i] , x[i]);
			Sb[7] += inner(p[i] , r[i]);
			Sb[8] += inner(p[i] , p[i]);
		}
#endif

		double cx;
		double cr;
		double cp=0.0;
		double ei_predict;
		if (is_first_time && (iter==0)) {
			Sa[2] = Sa[3];
			Sa[3] = Sa[4];
			Sb[2] = Sb[3];
			Sb[3] = Sb[4];
			//int info = lapack_ZHEEVX(&(Sa[0]), &(S_e_value[0]), &(S_e_vector[0]), &(triangle_factors[0]), 2);
			int info = lapack_DSYGVX(&(Sa[0]), &(Sb[0]), &(S_e_value[0]), &(S_e_vector[0]), &(triangle_factors[0]), 2);
			const int min_id = (S_e_value[0] < S_e_value[1]) ? 0 : 1;
			cx = S_e_vector[min_id * 2];
			cr = S_e_vector[min_id * 2 + 1];


			//LAPACK後のSa,Sbは書き換わっているので注意//
			//以下の計算はSaが書き換わっているので間違い//
			ei_predict = 0.0;
			double sb = 0.0;
			for (int j = 0; j < 2; ++j) {
				for (int i = 0; i < 2; ++i) {
					ei_predict += S_e_vector[min_id * 2 + j] * Sa[i * 2 + j] * S_e_vector[min_id * 2 + i];
					sb += S_e_vector[min_id * 2 + j] * Sb[i * 2 + j] * S_e_vector[min_id * 2 + i];
				}
			}
			ei_predict /= sb;

		} else {
			//int info = lapack_ZHEEVX(&(Sa[0]), &(S_e_value[0]), &(S_e_vector[0]), &(triangle_factors[0]), 3);
			int info = lapack_DSYGVX(&(Sa[0]), &(Sb[0]), &(S_e_value[0]), &(S_e_vector[0]), &(triangle_factors[0]), 3);
			const int min_id = (S_e_value[0] < S_e_value[1]) ? (S_e_value[0] < S_e_value[2]) ? 0 : 2 : (S_e_value[1] < S_e_value[2]) ? 1 : 2;
			cx = S_e_vector[min_id * 3];
			cr = S_e_vector[min_id * 3 + 1];
			cp = S_e_vector[min_id * 3 + 2];


			//LAPACK後のSa,Sbは書き換わっているので注意//
			//以下の計算はSaが書き換わっているので間違い//
			ei_predict = 0.0;
			double sb = 0.0;
			for (int j = 0; j < 3; ++j) {
				for (int i = 0; i < 3; ++i) {
					ei_predict += S_e_vector[min_id * 3 + j] * Sa[i * 3 + j] * S_e_vector[min_id * 3 + i];
					sb += S_e_vector[min_id * 3 + j] * Sb[i * 3 + j] * S_e_vector[min_id * 3 + i];
				}
			}
			ei_predict /= sb;

		}

		double tetAx = vecmath::Norm(Ax, N) * dVol;
		double tetAr = vecmath::Norm(Ar, N) * dVol;

		double norm_x = 0.0;
		double norm_p = 0.0;
		for (size_t i = 0; i < N; ++i) {
			x[i] = cx * x[i] + cr * r[i] + cp * p[i];
			p[i] = cr * r[i] + cp * p[i];
			Ax[i] = cx * Ax[i] + cr * Ar[i] + cp * Ap[i];
			Ap[i] = cr * Ar[i] + cp * Ap[i];

			norm_x += Norm(x[i]);	// inner(x[i], x[i]).real();
			norm_p += Norm(p[i]);	// inner(p[i], p[i]).real();
		}
		/* for test
		{
			double a_xAx = InnerProd(x, Ax, N);
			double a_pAx = InnerProd(p, Ax, N);
			double a_xAp = InnerProd(x, Ap, N);
			double a_pAp = InnerProd(p, Ap, N);
			double ei3 = a_xAx * cx * cx;
			ei3 += a_pAx * cp * cx;
		}

		double ei2 = InnerProd(x, Ax, N);
		ei2 /= InnerProd(x, x, N);
		*/
		double norm_ax = vecmath::Norm(Ax, N) * dVol;
		double norm_ap = vecmath::Norm(Ap, N) * dVol;

		norm_x = 1.0 / sqrt(norm_x * dVol);
		norm_p = 1.0 / sqrt(norm_p * dVol);
		MulC(x, norm_x, N);
		MulC(p, norm_p, N);
		MulC(Ax, norm_x, N);
		MulC(Ap, norm_p, N);

		ei = InnerProd(x, Ax, N) * dVol;
		double norm_r = 0.0;
		for (size_t i = 0; i < N; ++i) {
			r[i] = (Ax[i] - ei * x[i]);
			norm_r += Norm(r[i]);	// inner(r[i] , r[i]).real();
		}
		norm_r *= dVol;
		//収束判定////////////////////////////////////////////

		if (norm_r < limit * dVol) {
			//		printf("CG-fin: %d\n", count+1);
			break;
		}
			
		if (log_print_step > 0) {
			if ((iter + 1) % log_print_step == 0) {
				printf("(%d) ei, ei_predict, norm_r: %f, %f, %f\n", iter + 1, ei, ei_predict, norm_r);

				auto x_r = InnerProd(x, r, N);
				auto x_p = InnerProd(x, p, N);
				auto r_p = InnerProd(r, p, N);
				printf("   <x|r> = %f\n", x_r);
				printf("   <x|p> = %f\n", x_p);
				printf("   <r|p> = %f\n", r_p);
			}
		}
	}

	eigen_values[0] = ei;

#if 0
	{
		printf("check eigen vector\n");
		const int i_end = std::min<int>(100, N);
		for (int i = 0; i < i_end; i++) {
			auto diff = Ax[i] - ei * x[i];
			printf(" diff %f\n", diff);
		}

	}
#endif

	printf("LOBPCG iteration = %d:\n", iter);
	
	////////////////////////////End of loop of Lanczos to create Symmetric Tridiagonal T// 
	
	/*
	for (int k = 0; k < N; ++k) {
		eigen_vectors[k] = x[k];
	}
	*/

	//delete[] Ax;
	//delete[] x;
}


class LOBPCG_d_Solver {
private:
	double* work = nullptr;   //for Ax, r, and Ar
	double* keep = nullptr;   //for p and Ap
	int num_solution;
	size_t one_size;
	bool is_first_step = false;
public:
	~LOBPCG_d_Solver() {
		delete[] work;
		delete[] keep;
	}

	size_t Initialize(size_t N, const int num_solution_) {
		one_size = N;
		num_solution = num_solution_;
		work = new double[N* num_solution * 4];
		keep = new double[N* num_solution ];  //

		vecmath::SetZero(keep, N);

		is_first_step = true;

		return 5 * N* num_solution;
	}

	template <class OperationA>
	void Run(const double dVol, const int num_solution, const int iter_max,
			double* eigen_values, double* eigen_vectors, 
			OperationA OpeA, int log_print_step = 0) {

		EigenLOBPCG_d(one_size, dVol, iter_max,
			eigen_values, eigen_vectors, work, keep, is_first_step,
			OpeA, log_print_step);

		is_first_step = false;
	}

};


