#pragma once
#include "wrap_lapack.h"
#include "w_dsygvx.h"

#include "vecmath.h"
#include "inverse_m.h"

//LOBPCG法で複数の固有ベクトルを求める//

//対角化を手動でせずに、小行列の一般化固有値解法に任せた方が良い//
//なぜなら、(1)matrix Aの演算が1回/iterationで済む。//

//※非一般化の固有値問題を解き、かつ、r,pのnormalizeを忘れると収束しなくなる//



/*
固有値問題を解くLOBPCG

問題点、iter_max==1でも前回の履歴からr,p,Ar,Apを継続させるべき

*/
template <class OperationA, class WATCH>
inline
void EigenLOBPCG_d_multi(const int N, const double dVol, const int num_solution, const int iter_max,
	double* eigen_values, double* eigen_vectors, double* keep, double* work, bool is_first_time,
	OperationA OpeA, WATCH& watch) {

	using namespace vecmath;

	

	auto IntervalPointers = [](double* a, size_t N, size_t num_solution) {
		std::vector<double*> pointers(num_solution);
		for (size_t i = 0; i < num_solution; ++i) {
			pointers[i] = a + i * N;
		}
		return pointers;
		};


	auto x = IntervalPointers(eigen_vectors, N, num_solution);
	auto r = IntervalPointers(work, N, num_solution);
	auto p = IntervalPointers(keep, N, num_solution); 
	auto Ax = IntervalPointers(work + N * num_solution, N, num_solution);
	auto Ar = IntervalPointers(work + 2 * N * num_solution, N, num_solution);
	auto Ap = IntervalPointers(work + 3 * N * num_solution, N, num_solution);
	auto tmpM = IntervalPointers(work + 4 * N * num_solution, N, num_solution);


	const int n3 = 3 * num_solution;
	const int nn9 = (3 * num_solution) * (3 * num_solution);
	double* Sa = new double[nn9];
	double* Sb = new double[nn9];
	double* triangle_factors = new double[nn9];
	double* S_e_value = new double[(3 * num_solution)];
	double* S_e_vector = new double[nn9];

	const double limit = 1.0e-10;



	//calculte initial redidual r//
	{
		for (int k = 0; k < num_solution; ++k) {
			Normalize(x[k], N, dVol);
			watch.Record( 10);
			//Because Hamiltonian operator A was changed by changing electronic density rho, //
			// Ax and Ap should be re-calculated//
			OpeA(Ax[k], x[k]); // calculate Ax from x 
			OpeA(Ap[k], p[k]);

			double ei = InnerProd(x[k], Ax[k], N) * dVol;
			

			for (size_t i = 0; i < N; ++i) {
				r[k][i] = (Ax[k][i] - ei * x[k][i]);
			}
			double norm_r = Normalize(r[k], N, dVol);	//inner(r[i] , r[i]).real();


			if (is_first_time) {
				SetZero(p[k], N);
				SetZero(Ap[k], N);
			}

			printf("(0, %d) ei, norm_r: %f, %f\n",k, ei, norm_r);

			eigen_values[k] = ei;
		}
	}


//Loop of CG
	int iter;
	for (iter = 0; iter < iter_max; ++iter) {
		watch.Record( 10);
		for (int k = 0; k < num_solution; ++k) {
			OpeA(Ar[k], r[k]);
		}


		for (int j = 0; j < num_solution; ++j) {
			for (int k = 0; k < num_solution; ++k) {

				Sa[j * n3 + k] = InnerProd(x[j], Ax[k], N) * dVol;
				Sa[j * n3 + k + num_solution] = InnerProd(x[j], Ar[k], N) * dVol;
				Sa[j * n3 + k + 2 * num_solution] = InnerProd(x[j], Ap[k], N) * dVol;

				Sa[(j + num_solution) * n3 + k] = InnerProd(r[j], Ax[k], N) * dVol;
				Sa[(j + num_solution) * n3 + k + num_solution] = InnerProd(r[j], Ar[k], N) * dVol;
				Sa[(j + num_solution) * n3 + k + 2 * num_solution] = InnerProd(r[j], Ap[k], N) * dVol;

				Sa[(j + 2 * num_solution) * n3 + k] = InnerProd(p[j], Ax[k], N) * dVol;
				Sa[(j + 2 * num_solution) * n3 + k + num_solution] = InnerProd(p[j], Ar[k], N) * dVol;
				Sa[(j + 2 * num_solution) * n3 + k + 2 * num_solution] = InnerProd(p[j], Ap[k], N) * dVol;


				Sb[j * n3 + k] = InnerProd(x[j], x[k], N) * dVol;
				Sb[j * n3 + k + num_solution] = InnerProd(x[j], r[k], N) * dVol;
				Sb[j * n3 + k + 2 * num_solution] = InnerProd(x[j], p[k], N) * dVol;

				Sb[(j + num_solution) * n3 + k] = InnerProd(r[j], x[k], N) * dVol;
				Sb[(j + num_solution) * n3 + k + num_solution] = InnerProd(r[j], r[k], N) * dVol;
				Sb[(j + num_solution) * n3 + k + 2 * num_solution] = InnerProd(r[j], p[k], N) * dVol;

				Sb[(j + 2 * num_solution) * n3 + k] = InnerProd(p[j], x[k], N) * dVol;
				Sb[(j + 2 * num_solution) * n3 + k + num_solution] = InnerProd(p[j], r[k], N) * dVol;
				Sb[(j + 2 * num_solution) * n3 + k + 2 * num_solution] = InnerProd(p[j], p[k], N) * dVol;


			}
		}
		watch.Record( 17);

//#define DEBUG_PRINT1
#ifdef DEBUG_PRINT1
		if(iter==0){
			FILE* fp = fopen("Sa.txt", "w");
			for (int j = 0; j < n3; ++j) {
				for (int k = 0; k < n3 - 1; ++k) {
					fprintf(fp, "%f\t", Sa[k + n3 * j]);
				}
				fprintf(fp, "%f\n", Sa[n3-1 + n3*j]);
			}
			fclose(fp);
		}
#endif



#if 0
		//対称化(数値誤差により非対称になっているため)
		Sa[1] = (Sa[0 * 3 + 1] + Sa[1 * 3 + 0]) / 2.0;
		Sa[3] = Sa[1];

		Sa[2] = (Sa[0 * 3 + 2] + Sa[2 * 3 + 0]) / 2.0;
		Sa[6] = Sa[2];

		Sa[5] = (Sa[1 * 3 + 2] + Sa[2 * 3 + 1]) / 2.0;
		Sa[7] = Sa[5];


#endif




		//double cx;
		//double cr;
		//double cp=0.0;
		//double ei_predict;

		if (is_first_time && (iter==0)) {
			const int n2 = 2 * num_solution;
			double* Sa2 = new double[n2 * n2];
			double* Sb2 = new double[n2 * n2];

			//初回はp=0のため、3n行列を2n行列に縮小//
			for (int j = 0; j < num_solution; ++j) {
				for (int k = 0; k < num_solution; ++k) {
					Sa2[j * n2 + k] = Sa[j * n3 + k];
					Sa2[j * n2 + k + num_solution] = Sa[j * n3 + k + num_solution];
			
					Sa2[(j + num_solution) * n2 + k] = Sa[(j + num_solution) * n3 + k];
					Sa2[(j + num_solution) * n2 + k + num_solution] = Sa[(j + num_solution) * n3 + k + num_solution];
					
					Sb2[j * n2 + k] = Sb[j * n3 + k];
					Sb2[j * n2 + k + num_solution] = Sb[j * n3 + k + num_solution];

					Sb2[(j + num_solution) * n2 + k] = Sb[(j + num_solution) * n3 + k];
					Sb2[(j + num_solution) * n2 + k + num_solution] = Sb[(j + num_solution) * n3 + k + num_solution];
				}
			}
			watch.Record( 10);
#ifdef DEBUG_PRINT1
			if (iter == 0) {
				FILE* fp = fopen("Sa2.txt", "w");
				for (int j = 0; j < n2; ++j) {
					for (int k = 0; k < n2 - 1; ++k) {
						fprintf(fp, "%f\t", Sa2[k + n2 * j]);
					}
					fprintf(fp, "%f\n", Sa2[n2 - 1 + n2 * j]);
				}
				fclose(fp);
			}
#endif

			int info = lapack_DSYGVX(&(Sa2[0]), &(Sb2[0]), &(S_e_value[0]), &(S_e_vector[0]), &(triangle_factors[0]), 2 * num_solution);
			//const int min_id = (S_e_value[0] < S_e_value[1]) ? 0 : 1;
			//固有値の小さい順に並んでいるもの考えてよい//
			//cx = S_e_vector[min_id * 2];
			//cr = S_e_vector[min_id * 2 + 1];
			watch.Record( 15);

			delete[] Sa2;
			delete[] Sb2;

			/*
			X = {x[0], x[1], ..., x[num_solution-1]}
			X = X * Cx + R * Cr
			P = R * Cr
			AX = AX * Cx + AR * Cr
			AP = AR * Cr
			*/

			{
				//P_next = Cr * R
				for (int k = 0; k < num_solution; ++k) {
					SetZero(p[k], N);
					const double* cr = S_e_vector + k * n2 + num_solution;
					for (int j = 0; j < num_solution; ++j) {
						AddV(p[k], r[j], cr[j], N);
					}
				}

				auto& tmp = r;//rはこの時点で使わないので作業領域にする//
				//TMP = X * Cx + P_next;
				for (int k = 0; k < num_solution; ++k) {
					Copy(tmp[k], p[k], N);
					const double* cx = S_e_vector + k * n2;
					for (int j = 0; j < num_solution; ++j) {
						AddV(tmp[k], x[j], cx[j], N);
					}
				}

				//X_next = TMP				
				for (int k = 0; k < num_solution; ++k) {
					Copy(x[k], tmp[k], N);
				}
			}


			{
				//AP_next = AR * Cr
				for (int k = 0; k < num_solution; ++k) {
					SetZero(Ap[k], N);
					const double* cr = S_e_vector + k * n2 + num_solution;
					for (int j = 0; j < num_solution; ++j) {
						AddV(Ap[k], Ar[j], cr[j], N);
					}
				}

				auto& tmp = r;//rはこの時点で使わないので作業領域にする//
				//TMP = AX * Cx + AP_next;
				for (int k = 0; k < num_solution; ++k) {
					Copy(tmp[k], Ap[k], N);
					const double* cx = S_e_vector + k * n2;
					for (int j = 0; j < num_solution; ++j) {
						AddV(tmp[k], Ax[j], cx[j], N);
					}
				}

				//AX_next = TMP				
				for (int k = 0; k < num_solution; ++k) {
					Copy(Ax[k], tmp[k], N);
				}
			}

			watch.Record( 16);

		} else {
			//int info = lapack_ZHEEVX(&(Sa[0]), &(S_e_value[0]), &(S_e_vector[0]), &(triangle_factors[0]), 3);
			int info = lapack_DSYGVX(&(Sa[0]), &(Sb[0]), &(S_e_value[0]), &(S_e_vector[0]), &(triangle_factors[0]), 3 * num_solution);
			//固有値の小さい順に並んでいるもの考えてよい//
			watch.Record( 15);
			//const int min_id = (S_e_value[0] < S_e_value[1]) ? (S_e_value[0] < S_e_value[2]) ? 0 : 2 : (S_e_value[1] < S_e_value[2]) ? 1 : 2;
			//cx = S_e_vector[min_id * 3];
			//cr = S_e_vector[min_id * 3 + 1];
			//cp = S_e_vector[min_id * 3 + 2];
#if 1
/*
X = {x[0], x[1], ..., x[num_solution-1]}
X = X * Cx + R * Cr + P * Cp
P = R * Cr + P * Cp
AX = AX * Cx + AR * Cr + AP * Cp
AP = AR * Cr + AP * Cp
*/

			auto& tmp = tmpM;
			//TMP = R * Cr + P * Cp;
			for (int k = 0; k < num_solution; ++k) {
				const double* cx = S_e_vector + k * n3;
				const double* cr = S_e_vector + k * n3 + num_solution;
				const double* cp = S_e_vector + k * n3 + num_solution * 2;
				SetZero(tmp[k], N);
				for (int j = 0; j < num_solution; ++j) {
					AddV(tmp[k], r[j], cr[j], N);
					AddV(tmp[k], p[j], cp[j], N);
				}
			}

			//P_next = TMP;
			for (int k = 0; k < num_solution; ++k) {
				Copy(p[k], tmp[k], N);
			}

			//TMP = TMP + X * Cx;
			for (int k = 0; k < num_solution; ++k) {
				const double* cx = S_e_vector + k * n3;
				const double* cr = S_e_vector + k * n3 + num_solution;
				const double* cp = S_e_vector + k * n3 + num_solution * 2;
				for (int j = 0; j < num_solution; ++j) {
					AddV(tmp[k], x[j], cx[j], N);
				}
			}

			//X_next = TMP;
			for (int k = 0; k < num_solution; ++k) {
				Copy(x[k], tmp[k], N);
			}


			//TMP = AR * Cr + AP * Cp;
			for (int k = 0; k < num_solution; ++k) {
				const double* cx = S_e_vector + k * n3;
				const double* cr = S_e_vector + k * n3 + num_solution;
				const double* cp = S_e_vector + k * n3 + num_solution * 2;
				SetZero(tmp[k], N);
				for (int j = 0; j < num_solution; ++j) {
					AddV(tmp[k], Ar[j], cr[j], N);
					AddV(tmp[k], Ap[j], cp[j], N);
				}
			}

			//AP_next = TMP;
			for (int k = 0; k < num_solution; ++k) {
				Copy(Ap[k], tmp[k], N);
			}

			//TMP = TMP + AX * Cx;
			for (int k = 0; k < num_solution; ++k) {
				const double* cx = S_e_vector + k * n3;
				const double* cr = S_e_vector + k * n3 + num_solution;
				const double* cp = S_e_vector + k * n3 + num_solution * 2;
				for (int j = 0; j < num_solution; ++j) {
					AddV(tmp[k], Ax[j], cx[j], N);
				}
			}

			//AX_next = TMP;
			for (int k = 0; k < num_solution; ++k) {
				Copy(Ax[k], tmp[k], N);
			}


#else
			/*
			X = {x[0], x[1], ..., x[num_solution-1]}
			X = X * Cx + R * Cr + P * Cp
			P = R * Cr + P * Cp
			AX = AX * Cx + AR * Cr + AP * Cp
			AP = AR * Cr + AP * Cp

			//小メモリ化の為, 以下のように変形
			変形前の上記の場合は行列行列積の格納先としてn*Nの一時バッファが必要
			変形後はn*nのCp行列の逆行列を作ることで、行列積の格納先をPやRに指定できる
			P = ( R * (Cr * Cp^{-1}) + P ) * Cp
			X = X * Cx + P
			AP = ( AR * (Cr * Cp^{-1}) + AP) * Cp
			AX = AX * Cx + AP
			*/

			auto& tmp = Ar;
			double* Cp = new double[num_solution * num_solution];
			for (int k = 0; k < num_solution; ++k) {
				const double* cp = S_e_vector + k * n3 + num_solution * 2;
				for (int j = 0; j < num_solution; ++j) {
					Cp[k * num_solution + j] = cp[j];
				}
			}
			double* invCp = new double[num_solution * num_solution];
			InverseM(Cp, num_solution, invCp);//注意:元の行列も破壊される..

			double* Cr_invCp = new double[num_solution * num_solution];
			for (int k = 0; k < num_solution; ++k) {
				const double* cr = S_e_vector + num_solution;
				for (int j = 0; j < num_solution; ++j) {
					double sum = 0.0;
					for (int i = 0; i < num_solution; ++i) {
						sum +=  cr[i * n3 + j] * invCp[num_solution * k + i];
					}
					Cr_invCp[k * num_solution + j] = sum;
				}
			}
			
			{
				//P' = ((Cp^ { -1 }*Cr) * R + P)
				for (int k = 0; k < num_solution; ++k) {
					for (int j = 0; j < num_solution; ++j) {
						AddV(p[k], r[j], Cr_invCp[k * num_solution + j], N);
					}
				}

				auto& tmp = r;
				//TMP = Cp * P
				for (int k = 0; k < num_solution; ++k) {
					SetZero(tmp[k],N);
					const double* cp = S_e_vector + k * n3 + num_solution * 2;
					for (int j = 0; j < num_solution; ++j) {
						AddV(tmp[k], p[j], cp[j], N);
					}
				}

				//P_next = TMP				
				for (int k = 0; k < num_solution; ++k) {
					Copy(p[k], tmp[k], N);
				}

				//TMP = Cx * X + P_next;
				for (int k = 0; k < num_solution; ++k) {
					Copy(tmp[k], p[k], N);
					const double* cx = S_e_vector + k * n3;
					for (int j = 0; j < num_solution; ++j) {
						AddV(tmp[k], x[j], cx[j], N);
					}
				}

				//X_next = TMP				
				for (int k = 0; k < num_solution; ++k) {
					Copy(x[k], tmp[k], N);
				}
			}

			{

				//AP' = ((Cp^ { -1 }*Cr) * AR + AP)
				for (int k = 0; k < num_solution; ++k) {
					for (int j = 0; j < num_solution; ++j) {
						AddV(Ap[k], Ar[j], Cr_invCp[k * num_solution + j], N);
					}
				}

				auto& tmp = Ar;
				//TMP = Cp * AP
				for (int k = 0; k < num_solution; ++k) {
					SetZero(tmp[k],N);
					const double* cp = S_e_vector + k * n3 + num_solution * 2;
					for (int j = 0; j < num_solution; ++j) {
						AddV(tmp[k], Ap[j], cp[j], N);
					}
				}

				//AP_next = TMP				
				for (int k = 0; k < num_solution; ++k) {
					Copy(Ap[k], tmp[k], N);
				}

				//TMP = Cx * AX + AP_next;
				for (int k = 0; k < num_solution; ++k) {
					Copy(tmp[k], Ap[k], N);
					const double* cx = S_e_vector + k * n3;
					for (int j = 0; j < num_solution; ++j) {
						AddV(tmp[k], Ax[j], cx[j], N);
					}
				}

				//X_next = TMP				
				for (int k = 0; k < num_solution; ++k) {
					Copy(Ax[k], tmp[k], N);
				}
			}

			delete[] Cp;
			delete[] invCp;
			delete[] Cr_invCp;
#endif
			watch.Record( 16);
		}

		double norm_r_max = 0.0;

		for (int j = 0; j < num_solution; ++j) {

			double norm_ax = vecmath::Norm(Ax[j], N) * dVol;
			double norm_ap = vecmath::Norm(Ap[j], N) * dVol;

			double norm_x = vecmath::Norm(x[j], N);
			double norm_p = vecmath::Norm(p[j], N);
			norm_x = 1.0 / sqrt(norm_x * dVol);
			norm_p = 1.0 / sqrt(norm_p * dVol);
			MulC(x[j], norm_x, N);
			MulC(p[j], norm_p, N);
			MulC(Ax[j], norm_x, N);
			MulC(Ap[j], norm_p, N);

			double ei = InnerProd(x[j], Ax[j], N) * dVol;
			double norm_r = 0.0;
			for (size_t i = 0; i < N; ++i) {
				r[j][i] = (Ax[j][i] - ei * x[j][i]);
				norm_r += r[j][i] * r[j][i];	// inner(r[i] , r[i]).real();
			}
			norm_r *= dVol;

			if (norm_r_max < norm_r) norm_r_max = norm_r;

			printf("(%d, %d) ei, norm_r: %f, %f, %f\n", iter + 1, j, ei, norm_r, S_e_value[j]);
			
			eigen_values[j] = ei;
		}
		
		
		//収束判定////////////////////////////////////////////

		if (norm_r_max < limit * dVol) {
			//		printf("CG-fin: %d\n", count+1);
			break;
		}
			
		
	}

	

	printf("LOBPCG iteration = %d:\n", iter);
	
	////////////////////////////End of loop of Lanczos to create Symmetric Tridiagonal T// 
	
	delete[] Sa;
	delete[] Sb;
	delete[] triangle_factors;
	delete[] S_e_value;
	delete[] S_e_vector;
	watch.Record( 10);
}


inline
size_t WorkSizeLOBPCG_d_multi(int line_size, int num_solution) {
	return 5 * line_size * num_solution;
}
/*
class LOBPCG_d_multi_Solver{
private:
	double* keep = nullptr;   //for p and Ap
	int num_solution;
	size_t one_vector_size;
	bool is_first_step = false;
public:
	~LOBPCG_d_multi_Solver() {
		delete[] keep;
	}

	size_t Initialize(size_t N, const int num_solution_) {
		one_vector_size = N;
		num_solution = num_solution_;
		keep = new double[N* num_solution];

		vecmath::SetZero(keep, N);

		is_first_step = true;

		return 5 * N* num_solution;
	}

	template <class OperationA>
	void Run(const double dVol, const int num_solution, const int iter_max,
			double* eigen_values, double* eigen_vectors, double* work,
			OperationA OpeA) {

		EigenLOBPCG_d_multi(one_vector_size, dVol, num_solution, iter_max,
			eigen_values, eigen_vectors, keep, work, is_first_step,
			OpeA);

		is_first_step = false;
	}

};
*/

