#ifdef USE_MPI
#pragma once
#include <mpi.h>
#include "mpi_helper.h"
#include "wrap_lapack.h"
#include "w_dgemm.h"
#include "w_zhegvx.h"
#include "GridRange.h"
#include "GridScatterGather.h"

#include "vecmath.h"
#include "inverse_m.h"
#include "print_matrix.h"
#include "soacomplex.h"

#ifdef USE_SCALAPACK
#include "w_pzhegvx.h"
#endif

//LOBPCG法で複数の固有ベクトルを求める//
//複素数version//

//対角化を手動でせずに、小行列の一般化固有値解法に任せた方が良い//
//なぜなら、(1)matrix Aの演算が1回/iterationで済む。//

//dVolを乗じずに内積を取って1になるように関数の入り口でxをscaleする//
//出口でxを逆scaleして大きさを戻す//
//pはスケール後のものが関数から出力される//


#define USE_BLAS
//#ifdef _NEC
//#define SELF_TRANSPOSE
//#endif

//#define DEBUG_PRINT_MATRIX
//#define LOBPCG_DEBUG_PRINT

//XR, RP of S-matrix is directly calculated from Cx,Cr,Cp. But accuracy is low., and then convergence has problem.
#define DIRECT_S_MATRIX_XR_RP

//S行列を次のステップでも使いまわす場合//
//安定しないのでしばらくpending//
//#define S_MATRIX_RELOAD


//#define EVERY_NORMALIZE_P
// 
//数値誤差が発生時にrollbackとescapeをする実装(must)
#define CHECK_NEXT_VECTOR


//#define LOBPCG_Z_PRINT_EIGEN
//#define LOBPCG_Z_PRINT_NORM




namespace LOBPCG{


    inline constexpr double LIMIT_Z_NORM_X_FOR_RESET_MATRIX = 1.0e-3;
    inline constexpr double LIMIT_Z_NORM_R_FOR_ESCAPE = 1.0e-10;

inline
int MakeMatrix_z(int local_size, int num_solution, OneComplex* Sa, OneComplex* Sb,
	const SoAComplex* x, const SoAComplex* r, const SoAComplex* p,
	const SoAComplex* Ax, const SoAComplex* Ar, const SoAComplex* Ap) {

	using namespace vecmath;

	const int n3 = num_solution * 3;

	for (int j = 0; j < num_solution; ++j) {
		for (int k = 0; k < num_solution; ++k) {

			Sa[j * n3 + k] = SoAC::InnerProd(x[k], Ax[j], local_size);
			Sa[j * n3 + k + num_solution] = SoAC::InnerProd(r[k], Ax[j], local_size);
			Sa[j * n3 + k + 2 * num_solution] = SoAC::InnerProd(p[k], Ax[j], local_size);

			Sa[(j + num_solution) * n3 + k] = SoAC::InnerProd(x[k], Ar[j], local_size);
			Sa[(j + num_solution) * n3 + k + num_solution] = SoAC::InnerProd(r[k], Ar[j], local_size);
			Sa[(j + num_solution) * n3 + k + 2 * num_solution] = SoAC::InnerProd(p[k], Ar[j], local_size);

			Sa[(j + 2 * num_solution) * n3 + k] = SoAC::InnerProd(x[k], Ap[j], local_size);
			Sa[(j + 2 * num_solution) * n3 + k + num_solution] = SoAC::InnerProd(r[k], Ap[j], local_size);
			Sa[(j + 2 * num_solution) * n3 + k + 2 * num_solution] = SoAC::InnerProd(p[k], Ap[j], local_size);


			Sb[j * n3 + k] = SoAC::InnerProd(x[k], x[j], local_size);
			Sb[j * n3 + k + num_solution] = SoAC::InnerProd(r[k], x[j], local_size);
			Sb[j * n3 + k + 2 * num_solution] = SoAC::InnerProd(p[k], x[j], local_size);

			Sb[(j + num_solution) * n3 + k] = SoAC::InnerProd(x[k], r[j], local_size);
			Sb[(j + num_solution) * n3 + k + num_solution] = SoAC::InnerProd(r[k], r[j], local_size);
			Sb[(j + num_solution) * n3 + k + 2 * num_solution] = SoAC::InnerProd(p[k], r[j], local_size);

			Sb[(j + 2 * num_solution) * n3 + k] = SoAC::InnerProd(x[k], p[j], local_size);
			Sb[(j + 2 * num_solution) * n3 + k + num_solution] = SoAC::InnerProd(r[k], p[j], local_size);
			Sb[(j + 2 * num_solution) * n3 + k + 2 * num_solution] = SoAC::InnerProd(p[k], p[j], local_size);


		}
	}

	int num_gemm_call = 18;
	return num_gemm_call;
}

#ifdef USE_BLAS
#ifdef SELF_TRANSPOSE

inline
int MakeMatrix_z_blas(int local_size, int num_solution, OneComplex* Sa, OneComplex* Sb, double* temp_mat6n6n,
	const SoAComplex* x, const SoAComplex* r, const SoAComplex* p,
	const SoAComplex* Ax, const SoAComplex* Ar, const SoAComplex* Ap,
	SoAComplex* tmpM,
	bool is_skip_x_p, bool is_skip_Ax_Ap, bool is_skip_xr_rp) {

	const int N = num_solution;
	const int n2 = num_solution * 2;
	const int n3 = num_solution * 3;
	const int n6 = n3 * 2;


	auto AoSComplexFrom2by2Real = [](int N, OneComplex* Sa, int lda, double* temp_mat6n6n, int ldt) {
#pragma ivdep
		for (int m = 0; m < N; ++m) {
			for (int n = 0; n < N; ++n) {
				Sa[n + lda * m].r = temp_mat6n6n[(n * 2) + ldt * (m * 2)] + temp_mat6n6n[(n * 2 + 1) + ldt * (m * 2 + 1)];
				Sa[n + lda * m].i = -temp_mat6n6n[(n * 2 + 1) + ldt * (m * 2)] + temp_mat6n6n[(n * 2) + ldt * (m * 2 + 1)];
			}
		}
		};


	double* X_t = tmpM[0].re;
	transpose(local_size, n2, x[0].re, local_size, X_t, n2);


	if (!is_skip_Ax_Ap) {
		blas_DGEMM_n(n2, n2, local_size, X_t, Ax[0].re, temp_mat6n6n, n6, 1.0, 0.0);
		AoSComplexFrom2by2Real(num_solution, Sa, n3, temp_mat6n6n, n6);
	}
	blas_DGEMM_n(n2, n2, local_size, X_t, Ar[0].re, temp_mat6n6n + n2 * n6, n6, 1.0, 0.0);
	AoSComplexFrom2by2Real(num_solution, Sa + N * n3, n3, temp_mat6n6n + n2 * n6, n6);
	if (!is_skip_Ax_Ap) {
		blas_DGEMM_n(n2, n2, local_size, X_t, Ap[0].re, temp_mat6n6n + n2 * (2 * n6), n6, 1.0, 0.0);
		AoSComplexFrom2by2Real(num_solution, Sa + N * (2 * n3), n3, temp_mat6n6n + n2 * (2 * n6), n6);
	}

	if (!is_skip_x_p) {
		blas_DGEMM_n(n2, n2, local_size, X_t, x[0].re, temp_mat6n6n, n6, 1.0, 0.0);
		AoSComplexFrom2by2Real(num_solution, Sb, n3, temp_mat6n6n, n6);
	}
#ifdef DIRECT_S_MATRIX_XR_RP
	if (!is_skip_xr_rp)
#endif
	{
		blas_DGEMM_n(n2, n2, local_size, X_t, r[0].re, temp_mat6n6n + n2 * n6, n6, 1.0, 0.0);
		AoSComplexFrom2by2Real(num_solution, Sb + N * n3, n3, temp_mat6n6n + n2 * n6, n6);
	}
	if (!is_skip_x_p) {
		blas_DGEMM_n(n2, n2, local_size, X_t, p[0].re, temp_mat6n6n + n2 * (2 * n6), n6, 1.0, 0.0);
		AoSComplexFrom2by2Real(num_solution, Sb + N * (2 * n3), n3, temp_mat6n6n + n2 * (2 * n6), n6);
	}


	double* R_t = tmpM[0].re;
	transpose(local_size, n2, r[0].re, local_size, R_t, n2);

	blas_DGEMM_n(n2, n2, local_size, R_t, Ar[0].re, temp_mat6n6n + n2 * (n6 + 1), n6, 1.0, 0.0);
	AoSComplexFrom2by2Real(num_solution, Sa + N * (n3 + 1), n3, temp_mat6n6n + n2 * (n6 + 1), n6);

	blas_DGEMM_n(n2, n2, local_size, R_t, Ap[0].re, temp_mat6n6n + n2 * (2 * n6 + 1), n6, 1.0, 0.0);
	AoSComplexFrom2by2Real(num_solution, Sa + N * (2 * n3 + 1), n3, temp_mat6n6n + n2 * (2 * n6 + 1), n6);

	blas_DGEMM_n(n2, n2, local_size, R_t, r[0].re, temp_mat6n6n + n2 * (n6 + 1), n6, 1.0, 0.0);
	AoSComplexFrom2by2Real(num_solution, Sb + N * (n3 + 1), n3, temp_mat6n6n + n2 * (n6 + 1), n6);

#ifdef DIRECT_S_MATRIX_XR_RP
	if (!is_skip_xr_rp)
#endif
	{
		blas_DGEMM_n(n2, n2, local_size, R_t, p[0].re, temp_mat6n6n + n2 * (2 * n6 + 1), n6, 1.0, 0.0);
		AoSComplexFrom2by2Real(num_solution, Sb + N * (2 * n3 + 1), n3, temp_mat6n6n + n2 * (2 * n6 + 1), n6);
	}


	double* P_t = tmpM[0].re;
	transpose(local_size, n2, p[0].re, local_size, P_t, n2);

	if (!is_skip_Ax_Ap) {
		blas_DGEMM_n(n2, n2, local_size, P_t, Ap[0].re, temp_mat6n6n + n2 * (2 * n6 + 2), n6, 1.0, 0.0);
		AoSComplexFrom2by2Real(num_solution, Sa + N * (2 * n3 + 2), n3, temp_mat6n6n + n2 * (2 * n6 + 2), n6);
	}
	if (!is_skip_x_p) {
		blas_DGEMM_n(n2, n2, local_size, P_t, p[0].re, temp_mat6n6n + n2 * (2 * n6 + 2), n6, 1.0, 0.0);
		AoSComplexFrom2by2Real(num_solution, Sb + N * (2 * n3 + 2), n3, temp_mat6n6n + n2 * (2 * n6 + 2), n6);
	}


	//number of call of dgemm//
	int num_gemm_call = 18 - 6; //6 means lower triangle
	if (is_skip_x_p) num_gemm_call -= 3;
	if (is_skip_Ax_Ap) num_gemm_call -= 3;
#ifdef DIRECT_S_MATRIX_XR_RP
	if (is_skip_xr_rp) num_gemm_call -= 2;
#endif
	return num_gemm_call;
}

#else
inline
int MakeMatrix_z_blas(int local_size, int num_solution, OneComplex* Sa, OneComplex* Sb, double* temp_mat6n6n,
	const SoAComplex* x, const SoAComplex* r, const SoAComplex* p,
	const SoAComplex* Ax, const SoAComplex* Ar, const SoAComplex* Ap,
	bool is_skip_x_p, bool is_skip_Ax_Ap, bool is_skip_xr_rp) {

	const int N = num_solution;
	const int n2 = num_solution*2;
	const int n3 = num_solution * 3;
	const int n6 = n3 * 2;

	auto AoSComplexFrom2by2Real = [](int N, OneComplex* Sa, int lda, double* temp_mat6n6n, int ldt) {
		#pragma ivdep
		for (int m = 0; m < N; ++m) {
			for (int n = 0; n < N; ++n) {
				Sa[n + lda * m].r = temp_mat6n6n[(n * 2) + ldt * (m * 2)] + temp_mat6n6n[(n * 2 + 1) + ldt * (m * 2 + 1)];
				Sa[n + lda * m].i = -temp_mat6n6n[(n * 2 + 1) + ldt * (m * 2)] + temp_mat6n6n[(n * 2) + ldt * (m * 2 + 1)];
			}
		}
	};


	if (!is_skip_Ax_Ap) {
		blas_DGEMM_t(n2, n2, local_size, x[0].re, Ax[0].re, temp_mat6n6n, n6, 1.0, 0.0);
		AoSComplexFrom2by2Real(num_solution, Sa, n3, temp_mat6n6n, n6);
	}
	blas_DGEMM_t(n2, n2, local_size, x[0].re, Ar[0].re, temp_mat6n6n + n2 * n6, n6, 1.0, 0.0);
	AoSComplexFrom2by2Real(num_solution, Sa + N * n3, n3, temp_mat6n6n + n2 * n6, n6);

	blas_DGEMM_t(n2, n2, local_size, r[0].re, Ar[0].re, temp_mat6n6n + n2 * (n6 + 1), n6, 1.0, 0.0);
	AoSComplexFrom2by2Real(num_solution, Sa + N * (n3 + 1), n3, temp_mat6n6n + n2 * (n6 + 1), n6);

	if (!is_skip_Ax_Ap) {
		blas_DGEMM_t(n2, n2, local_size, x[0].re, Ap[0].re, temp_mat6n6n + n2 * (2 * n6), n6, 1.0, 0.0);
		AoSComplexFrom2by2Real(num_solution, Sa + N * (2 * n3), n3, temp_mat6n6n + n2 * (2 * n6), n6);
	}
	blas_DGEMM_t(n2, n2, local_size, r[0].re, Ap[0].re, temp_mat6n6n + n2 * (2 * n6 + 1), n6, 1.0, 0.0);
	AoSComplexFrom2by2Real(num_solution, Sa + N * (2 * n3 + 1), n3, temp_mat6n6n + n2 * (2 * n6 + 1), n6);
	if (!is_skip_Ax_Ap) {
		blas_DGEMM_t(n2, n2, local_size, p[0].re, Ap[0].re, temp_mat6n6n + n2 * (2 * n6 + 2), n6, 1.0, 0.0);
		AoSComplexFrom2by2Real(num_solution, Sa + N * (2 * n3 + 2), n3, temp_mat6n6n + n2 * (2 * n6 + 2), n6);
	}




	if (!is_skip_x_p) {
		blas_DGEMM_t(n2, n2, local_size, x[0].re, x[0].re, temp_mat6n6n, n6, 1.0, 0.0);
		AoSComplexFrom2by2Real(num_solution, Sb, n3, temp_mat6n6n, n6);
	}
#ifdef DIRECT_S_MATRIX_XR_RP
	if (!is_skip_xr_rp) 
#endif
	{
		blas_DGEMM_t(n2, n2, local_size, x[0].re, r[0].re, temp_mat6n6n + n2 * n6, n6, 1.0, 0.0);
		AoSComplexFrom2by2Real(num_solution, Sb + N * n3, n3, temp_mat6n6n + n2 * n6, n6);
	}
	blas_DGEMM_t(n2, n2, local_size, r[0].re, r[0].re, temp_mat6n6n + n2 * (n6 + 1), n6, 1.0, 0.0);
	AoSComplexFrom2by2Real(num_solution, Sb + N * (n3 + 1), n3, temp_mat6n6n + n2 * (n6 + 1), n6);

	if (!is_skip_x_p) {
		blas_DGEMM_t(n2, n2, local_size, x[0].re, p[0].re, temp_mat6n6n + n2 * (2 * n6), n6, 1.0, 0.0);
		AoSComplexFrom2by2Real(num_solution, Sb + N * (2 * n3), n3, temp_mat6n6n + n2 * (2 * n6), n6);
	}
#ifdef DIRECT_S_MATRIX_XR_RP
	if (!is_skip_xr_rp) 
#endif
	{
		blas_DGEMM_t(n2, n2, local_size, r[0].re, p[0].re, temp_mat6n6n + n2 * (2 * n6 + 1), n6, 1.0, 0.0);
		AoSComplexFrom2by2Real(num_solution, Sb + N * (2 * n3 + 1), n3, temp_mat6n6n + n2 * (2 * n6 + 1), n6);
	}
	if (!is_skip_x_p) {
		blas_DGEMM_t(n2, n2, local_size, p[0].re, p[0].re, temp_mat6n6n + n2 * (2 * n6 + 2), n6, 1.0, 0.0);
		AoSComplexFrom2by2Real(num_solution, Sb + N * (2 * n3 + 2), n3, temp_mat6n6n + n2 * (2 * n6 + 2), n6);
	}


	//number of call of dgemm//
	int num_gemm_call = 18 - 6; //6 means lower triangle
	if (is_skip_x_p) num_gemm_call -= 3;
	if (is_skip_Ax_Ap) num_gemm_call -= 3;
#ifdef DIRECT_S_MATRIX_XR_RP
	if (is_skip_xr_rp) num_gemm_call -= 2;
#endif
	return num_gemm_call;
}
#endif

inline
void CopySubMat(int Nx, int Ny, OneComplex* dest, int ldd, const OneComplex* src, int lds) {
	//column major//
	#pragma ivdep
	for (int ix = 0; ix < Nx; ++ix) {
		for (int iy = 0; iy < Ny; ++iy) {
			dest[iy + ldd * ix] = src[iy + lds * ix];
		}
	}
}

inline
int NextMatrix_z(int num_solution, OneComplex* Sa, OneComplex* Sb,
	const OneComplex* CxCrCp, const double* eigen_values, bool is_use_p, bool need_H_matrix) {

	using namespace vecmath;
	DEBUG_PRINTF("Use BLAS in MakeMatrix_z_mat\n");


	const int N = num_solution;
	const int n3 = N * 3;
	const int n2 = N * 2;
	const int nn = N * N;
	const int next = n3 * n3 * 2;
	const OneComplex* Cx = CxCrCp;
	const OneComplex* Cr = CxCrCp + N;
	const OneComplex* Cp = CxCrCp + N * 2;
	const int stride_C = (is_use_p) ? num_solution * 3 : num_solution * 2;
	const int stride_S = num_solution * 3;

	const OneComplex ONE{ 1.0,0.0 };
	const OneComplex ZERO{ 0.0,0.0 };



	{
		///for direct update (P^t P), (X^t P), (X^t X) /////////////////////////////////
		OneComplex* XtX = Sb;
		OneComplex* XtR = Sb + N * n3;
		OneComplex* XtP = Sb + N * 2 * n3;
		OneComplex* RtR = Sb + N * (n3 + 1);
		OneComplex* RtP = Sb + N * (2 * n3 + 1);
		OneComplex* PtP = Sb + N * (2 * n3 + 2);

		OneComplex* tmp1 = Sb + next;
		OneComplex* tmpPtP = Sb + next + nn;
		OneComplex* tmpXtP = Sb + next + nn * 2;
		OneComplex* tmpXtX = Sb + next + nn * 3;

		/*
		P'= R * Cr + P * Cp
		(P'^t P') = Cr^t [(R^t R) Cr + (R^t P) Cp] + Cp^t [(P^t R) Cr + (P^t P) Cp]
		*/
		//tmp1 = (R ^ t R) Cr
		cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, N, N, N, &ONE, RtR, stride_S, Cr, stride_C, &ZERO, tmp1, N);
		if (is_use_p) {
			//tmp1 += (R ^ t P) Cp
			cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, N, N, N, &ONE, RtP, stride_S, Cp, stride_C, &ONE, tmp1, N);
		}
		//tmpPtP = Cr^t tmp1
		cblas_zgemm(CblasColMajor, CblasConjTrans, CblasNoTrans, N, N, N, &ONE, Cr, stride_C, tmp1, N, &ZERO, tmpPtP, N);

		if (is_use_p) {
			//tmp1 = (P ^ t R) Cr = (R ^ t P)^t Cr
			//tmp1 += (P ^ t P) Cp
			cblas_zgemm(CblasColMajor, CblasConjTrans, CblasNoTrans, N, N, N, &ONE, RtP, stride_S, Cr, stride_C, &ZERO, tmp1, N);
			cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, N, N, N, &ONE, PtP, stride_S, Cp, stride_C, &ONE, tmp1, N);
			//tmpPtP += Cp^t tmp1
			cblas_zgemm(CblasColMajor, CblasConjTrans, CblasNoTrans, N, N, N, &ONE, Cp, stride_C, tmp1, N, &ONE, tmpPtP, N);
		}


		//(X'^t P') = (P'^t P') + Cx^t [(X^t R) Cr + (X^t P) Cp]
		//tmp1 = (X^t R) Cr
		cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, N, N, N, &ONE, XtR, stride_S, Cr, stride_C, &ZERO, tmp1, N);
		if (is_use_p) {
			//tmp1 += (X^t P) Cp
			cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, N, N, N, &ONE, XtP, stride_S, Cp, stride_C, &ONE, tmp1, N);
		}

		//tmpXtP = tmpPtP
		CopySubMat(N, N, tmpXtP, N, tmpPtP, N);
		//tmpXtP += Cx^t tmp1
		cblas_zgemm(CblasColMajor, CblasConjTrans, CblasNoTrans, N, N, N, &ONE, Cx, stride_C, tmp1, N, &ONE, tmpXtP, N);

		//(X'^t X') = (X'^t P') + [(X^t X) Cx + (X^t R) Cr + (X^t P)Cp]^t Cx
		//here, tmp1 == (X^t R) Cr + (X^t P)Cp
		//tmp1 += (X^t X) Cx
		cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, N, N, N, &ONE, XtX, stride_S, Cx, stride_C, &ONE, tmp1, N);
		//tmpXtX = tmpXtP
		CopySubMat(N, N, tmpXtX, N, tmpXtP, N);
		//tmpXtX += tmp1^t Cx
		cblas_zgemm(CblasColMajor, CblasConjTrans, CblasNoTrans, N, N, N, &ONE, tmp1, N, Cx, stride_C, &ONE, tmpXtX, N);

		//P'^t P' = tmpPtP;
		CopySubMat(N, N, PtP, n3, tmpPtP, N);

		//X'^t P' = tmpXtP;
		CopySubMat(N, N, XtP, n3, tmpXtP, N);

		//X'^t X' = tmpXtX;
		CopySubMat(N, N, XtX, n3, tmpXtX, N);

	}


	int num_gemm_call = is_use_p ? 11 : 6;
	if (!need_H_matrix) return num_gemm_call;
	
	{
		///for direct update (P^t HP), (X^t HP), (X^t HX) /////////////////////////////////
		OneComplex* XtHX = Sa;
		OneComplex* XtHR = Sa + N * n3;
		OneComplex* XtHP = Sa + N * 2 * n3;
		OneComplex* RtHR = Sa + N * (n3 + 1);
		OneComplex* RtHP = Sa + N * (2 * n3 + 1);
		OneComplex* PtHP = Sa + N * (2 * n3 + 2);

		OneComplex* tmp1 = Sa + next;
		OneComplex* tmpPtHP = Sa + next + nn;
		OneComplex* tmpXtHP = Sa + next + nn * 2;
		OneComplex* tmpXtHX = Sa + next + nn * 3;

		/*
		P'= R * Cr + P * Cp
		(P'^t HP') = Cr^t [(R^t HR) Cr + (R^t HP) Cp] + Cp^t [(P^t HR) Cr + (P^t HP) Cp]
		*/
		//tmp1 = (R ^ t HR) Cr
		cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, N, N, N, &ONE, RtHR, stride_S, Cr, stride_C, &ZERO, tmp1, N);
		if (is_use_p) {
			//tmp1 += (R ^ t HP) Cp
			cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, N, N, N, &ONE, RtHP, stride_S, Cp, stride_C, &ONE, tmp1, N);
		}
		//tmpPtHP = Cr^t tmp1
		cblas_zgemm(CblasColMajor, CblasConjTrans, CblasNoTrans, N, N, N, &ONE, Cr, stride_C, tmp1, N, &ZERO, tmpPtHP, N);

		if (is_use_p) {
			//tmp1 = (P ^ t HR) Cr = (R ^ t HP)^t Cr
			//tmp1 += (P ^ t HP) Cp
			cblas_zgemm(CblasColMajor, CblasConjTrans, CblasNoTrans, N, N, N, &ONE, RtHP, stride_S, Cr, stride_C, &ZERO, tmp1, N);
			cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, N, N, N, &ONE, PtHP, stride_S, Cp, stride_C, &ONE, tmp1, N);
			//tmpPtHP += Cp^t tmp1
			cblas_zgemm(CblasColMajor, CblasConjTrans, CblasNoTrans, N, N, N, &ONE, Cp, stride_C, tmp1, N, &ONE, tmpPtHP, N);
		}


		//(X'^t HP') = (P'^t HP') + Cx^t [(X^t HR) Cr + (X^t HP) Cp]
		//tmp1 = (X^t HR) Cr
		cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, N, N, N, &ONE, XtHR, stride_S, Cr, stride_C, &ZERO, tmp1, N);
		if (is_use_p) {
			//tmp1 += (X^t HP) Cp
			cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, N, N, N, &ONE, XtHP, stride_S, Cp, stride_C, &ONE, tmp1, N);
		}

		//tmpXtHP = tmpPtHP
		CopySubMat(N, N, tmpXtHP, N, tmpPtHP, N);
		//tmpXtHP += Cx^t tmp1
		cblas_zgemm(CblasColMajor, CblasConjTrans, CblasNoTrans, N, N, N, &ONE, Cx, stride_C, tmp1, N, &ONE, tmpXtHP, N);

		//(X'^t HX') = (X'^t HP') + [(X^t HX) Cx + (X^t HR) Cr + (X^t HP)Cp]^t Cx
		//here, tmp1 == (X^t HR) Cr + (X^t HP)Cp
		//tmp1 += (X^t HX) Cx
		cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, N, N, N, &ONE, XtHX, stride_S, Cx, stride_C, &ONE, tmp1, N);
		//tmpXtHX = tmpXtP
		CopySubMat(N, N, tmpXtHX, N, tmpXtHP, N);
		//tmpXtHX += tmp1^t Cx
		cblas_zgemm(CblasColMajor, CblasConjTrans, CblasNoTrans, N, N, N, &ONE, tmp1, N, Cx, stride_C, &ONE, tmpXtHX, N);

		//P'^t HP' = tmpPtHP;
		CopySubMat(N, N, PtHP, n3, tmpPtHP, N);

		//X'^t HP' = tmpXtHP;
		CopySubMat(N, N, XtHP, n3, tmpXtHP, N);

		//X'^t HX' = tmpXtHX;
		CopySubMat(N, N, XtHX, n3, tmpXtHX, N);
	}


#ifdef DIRECT_S_MATRIX_XR_RP
	{//direct update  (X^t R) and (R^t P)
		{
			//(X'^t R') = (X'^t H X') - (X'^t X')E'//
			OneComplex* XtR = Sb + N * n3;
			int ldm = n3;
			const OneComplex* tmpXtHX = Sa + next + nn * 3;
			const OneComplex* tmpXtX = Sb + next + nn * 3;
			for (int j = 0; j < N; ++j) {
				for (int i = 0; i < N; ++i) {
					XtR[i + ldm * j].r = tmpXtHX[i + N * j].r - tmpXtX[i + N * j].r * eigen_values[j];
					XtR[i + ldm * j].i = tmpXtHX[i + N * j].i - tmpXtX[i + N * j].i * eigen_values[j];
				}
			}
		}

		{
			//(R'^t P') = (X'^t H P') - E'(X'^t P')//
			OneComplex* RtP = Sb + N * (2 * n3 + 1);
			int ldm = n3;
			const OneComplex* tmpXtHP = Sa + next + nn * 2;
			const OneComplex* tmpXtP = Sb + next + nn * 2;
			for (int j = 0; j < N; ++j) {
				for (int i = 0; i < N; ++i) {
					RtP[i + ldm * j].r = tmpXtHP[i + N * j].r - eigen_values[i] * tmpXtP[i + N * j].r;
					RtP[i + ldm * j].i = tmpXtHP[i + N * j].i - eigen_values[i] * tmpXtP[i + N * j].i;
				}
			}
		}
	}
#endif

	num_gemm_call += is_use_p ? 11 : 6;
	return num_gemm_call;
}

void ScaleHSMatrix_z(int N, OneComplex* Sa, OneComplex* Sb, const double* sqrt_norms_x, const double* sqrt_norms_p, const double* sqrt_norms_r, bool need_H_matrix) {
	const int n3 = N * 3;
	const int ldd = n3;

	auto ScaleMatrix = [](int N, OneComplex* m, int ldm, const double* scale_i, const double* scale_j) {
		if (scale_j == nullptr) {
			//nothing to do//
		} else if (scale_i == nullptr) {
			for (int j = 0; j < N; ++j) {
				for (int i = 0; i < N; ++i) {
					m[i + ldm * j].r *= scale_j[j];
					m[i + ldm * j].i *= scale_j[j];
				}
			}
		} else {
			for (int j = 0; j < N; ++j) {
				for (int i = 0; i < N; ++i) {
					m[i + ldm * j].r *= scale_i[i] * scale_j[j];
					m[i + ldm * j].i *= scale_i[i] * scale_j[j];
				}
			}
		}
		};


	//for
	OneComplex* XtX = Sb;
	OneComplex* XtP = Sb + N * 2 * n3;
	OneComplex* PtP = Sb + N * (2 * n3 + 2);
	ScaleMatrix(N, PtP, ldd, sqrt_norms_p, sqrt_norms_p);
	ScaleMatrix(N, XtP, ldd, sqrt_norms_x, sqrt_norms_p);
	ScaleMatrix(N, XtX, ldd, sqrt_norms_x, sqrt_norms_x);

	if (!need_H_matrix) return;

	OneComplex* XtHX = Sa;
	OneComplex* XtHP = Sa + N * 2 * n3;
	OneComplex* PtHP = Sa + N * (2 * n3 + 2);
	ScaleMatrix(N, PtHP, ldd, sqrt_norms_p, sqrt_norms_p);
	ScaleMatrix(N, XtHP, ldd, sqrt_norms_x, sqrt_norms_p);
	ScaleMatrix(N, XtHX, ldd, sqrt_norms_x, sqrt_norms_x);

#ifdef DIRECT_S_MATRIX_XR_RP
	OneComplex* XtR = Sb + N * n3;
	OneComplex* RtP = Sb + N * (2 * n3 + 1);
	ScaleMatrix(N, XtR, ldd, sqrt_norms_x, sqrt_norms_r);
	ScaleMatrix(N, RtP, ldd, sqrt_norms_r, sqrt_norms_p);
#endif
}




/*
X = {x[0], x[1], ..., x[num_solution-1]}
X = X * Cx + R * Cr + P * Cp
P = R * Cr + P * Cp
AX = AX * Cx + AR * Cr + AP * Cp
AP = AR * Cr + AP * Cp
* where when the first update (is_use_p==false),
X = X * Cx + R * Cr
P = R * Cr
AX = AX * Cx + AR * Cr
AP = AR * Cr
* because CG vector P == 0.
* 
* When need_Ax is true, Ap and Ap are calculated
*/
inline
int NextVector_z(int local_size, int num_solution, const OneComplex* CxCrCp,
	bool is_use_p, bool need_Ax,
	SoAComplex* x, SoAComplex* r, SoAComplex* p,
	SoAComplex* Ax, SoAComplex* Ar, SoAComplex* Ap,
	SoAComplex* tmp, OneComplex* rotC) {
	//tmp : required size is (num_solution*local_size) //
	//rotC: required size is (num_solution*num_solution*3*2) //

	const int n3 = num_solution * 3;
	//double* Cx = (double*)CxCrCp;
	//double* Cr = (double*)(CxCrCp + num_solution);
	//double* Cp = (double*)(CxCrCp + num_solution * 2);
	double* rotCx = (double*)rotC;
	double* rotCr = (double*)(rotC + num_solution);
	double* rotCp = (double*)(rotC + num_solution * 2);
	OneComplex* cnjC = rotC + num_solution * n3;
	double* cnjCx = (double*)(cnjC);
	double* cnjCr = (double*)(cnjC + num_solution);
	double* cnjCp = (double*)(cnjC + num_solution * 2);

	const int stride_C = (is_use_p) ? num_solution * 6 : num_solution * 4;

	if (is_use_p) {
		
		for (int m = 0; m < num_solution; ++m) {
			for (int n = 0; n < n3; ++n) {
				rotC[n + n3 * m].r = CxCrCp[n + n3 * m].i;
				rotC[n + n3 * m].i = CxCrCp[n + n3 * m].r;
			}
		}
		for (int m = 0; m < num_solution; ++m) {
			for (int n = 0; n < n3; ++n) {
				cnjC[n + n3 * m].r = CxCrCp[n + n3 * m].r;
				cnjC[n + n3 * m].i = -CxCrCp[n + n3 * m].i;
			}
		}
	} else {
		const int n2 = num_solution * 2;
		for (int m = 0; m < num_solution; ++m) {
			for (int n = 0; n < n2; ++n) {
				rotC[n + n2 * m].r = CxCrCp[n + n2 * m].i;
				rotC[n + n2 * m].i = CxCrCp[n + n2 * m].r;
			}
		}

		for (int m = 0; m < num_solution; ++m) {
			for (int n = 0; n < n2; ++n) {
				cnjC[n + n2 * m].r = CxCrCp[n + n2 * m].r;
				cnjC[n + n2 * m].i = -CxCrCp[n + n2 * m].i;
			}
		}
	}

	//TMP = R * Cr + P * Cp;
	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, local_size, num_solution, num_solution * 2, 1.0, r[0].re, local_size, cnjCr, stride_C, 0.0, tmp[0].re, local_size * 2);
	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, local_size, num_solution, num_solution * 2, 1.0, r[0].re, local_size, rotCr, stride_C, 0.0, tmp[0].im, local_size * 2);

	if (is_use_p) {
		cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, local_size, num_solution, num_solution * 2, 1.0, p[0].re, local_size, cnjCp, stride_C, 1.0, tmp[0].re, local_size * 2);
		cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, local_size, num_solution, num_solution * 2, 1.0, p[0].re, local_size, rotCp, stride_C, 1.0, tmp[0].im, local_size * 2);
	}

	//P_next = TMP;
	for (int k = 0; k < num_solution; ++k) {
		SoAC::Copy(p[k], tmp[k], local_size);
	}

	//TMP = TMP + X * Cx;
	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, local_size, num_solution, num_solution * 2, 1.0, x[0].re, local_size, cnjCx, stride_C, 1.0, tmp[0].re, local_size * 2);
	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, local_size, num_solution, num_solution * 2, 1.0, x[0].re, local_size, rotCx, stride_C, 1.0, tmp[0].im, local_size * 2);

	//X_next = TMP;
	for (int k = 0; k < num_solution; ++k) {
		SoAC::Copy(x[k], tmp[k], local_size);
	}

	if (need_Ax) {
		//TMP = AR * Cr + AP * Cp;
		cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, local_size, num_solution, num_solution * 2, 1.0, Ar[0].re, local_size, cnjCr, stride_C, 0.0, tmp[0].re, local_size * 2);
		cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, local_size, num_solution, num_solution * 2, 1.0, Ar[0].re, local_size, rotCr, stride_C, 0.0, tmp[0].im, local_size * 2);

		if (is_use_p) {
			cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, local_size, num_solution, num_solution * 2, 1.0, Ap[0].re, local_size, cnjCp, stride_C, 1.0, tmp[0].re, local_size * 2);
			cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, local_size, num_solution, num_solution * 2, 1.0, Ap[0].re, local_size, rotCp, stride_C, 1.0, tmp[0].im, local_size * 2);
		}

		//AP_next = TMP;
		for (int k = 0; k < num_solution; ++k) {
			SoAC::Copy(Ap[k], tmp[k], local_size);
		}

		//TMP = TMP + AX * Cx;
		cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, local_size, num_solution, num_solution * 2, 1.0, Ax[0].re, local_size, cnjCx, stride_C, 1.0, tmp[0].re, local_size * 2);
		cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, local_size, num_solution, num_solution * 2, 1.0, Ax[0].re, local_size, rotCx, stride_C, 1.0, tmp[0].im, local_size * 2);

		//AX_next = TMP;
		for (int k = 0; k < num_solution; ++k) {
			SoAC::Copy(Ax[k], tmp[k], local_size);
		}
	}

	int num_gemm_call = is_use_p ? 6 : 4;
	if(need_Ax) num_gemm_call += is_use_p ? 6 : 4;
	return num_gemm_call;
}



/*
X = {x[0], x[1], ..., x[num_solution-1]}
X' = X * Cx + R * Cr + P * Cp
P' = R * Cr + P * Cp
* where when the first update (is_use_p==false),
X' = X * Cx + R * Cr
P' = R * Cr
because CG vector P == 0.
Pay attention. P' is written in buffer of P, while X' is witten in tmp. 
Then user should be copy tmp to buffer of X after this function call.

And, this function can be used for update of AX and AP when argument is 
AX = AX * Cx + AR * Cr + AP * Cp
AP = AR * Cr + AP * Cp
AX = AX * Cx + AR * Cr
AP = AR * Cr
*
* When need_Ax is true, Ap and Ap are calculated
*/
inline
int NextVector_z2(int local_size, int num_solution, const OneComplex* CxCrCp,
    bool is_use_p, 
    SoAComplex* x, SoAComplex* r, SoAComplex* p,
    SoAComplex* tmp, OneComplex* rotC) {
    //tmp : required size is (num_solution*local_size) //
    //rotC: required size is (num_solution*num_solution*3*2) //

    const int n3 = num_solution * 3;
    //double* Cx = (double*)CxCrCp;
    //double* Cr = (double*)(CxCrCp + num_solution);
    //double* Cp = (double*)(CxCrCp + num_solution * 2);
    double* rotCx = (double*)rotC;
    double* rotCr = (double*)(rotC + num_solution);
    double* rotCp = (double*)(rotC + num_solution * 2);
    OneComplex* cnjC = rotC + num_solution * n3;
    double* cnjCx = (double*)(cnjC);
    double* cnjCr = (double*)(cnjC + num_solution);
    double* cnjCp = (double*)(cnjC + num_solution * 2);

    const int stride_C = (is_use_p) ? num_solution * 6 : num_solution * 4;

    if (is_use_p) {

        for (int m = 0; m < num_solution; ++m) {
            for (int n = 0; n < n3; ++n) {
                rotC[n + n3 * m].r = CxCrCp[n + n3 * m].i;
                rotC[n + n3 * m].i = CxCrCp[n + n3 * m].r;
            }
        }
        for (int m = 0; m < num_solution; ++m) {
            for (int n = 0; n < n3; ++n) {
                cnjC[n + n3 * m].r = CxCrCp[n + n3 * m].r;
                cnjC[n + n3 * m].i = -CxCrCp[n + n3 * m].i;
            }
        }
    }
    else {
        const int n2 = num_solution * 2;
        for (int m = 0; m < num_solution; ++m) {
            for (int n = 0; n < n2; ++n) {
                rotC[n + n2 * m].r = CxCrCp[n + n2 * m].i;
                rotC[n + n2 * m].i = CxCrCp[n + n2 * m].r;
            }
        }

        for (int m = 0; m < num_solution; ++m) {
            for (int n = 0; n < n2; ++n) {
                cnjC[n + n2 * m].r = CxCrCp[n + n2 * m].r;
                cnjC[n + n2 * m].i = -CxCrCp[n + n2 * m].i;
            }
        }
    }

    //TMP = R * Cr + P * Cp;
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, local_size, num_solution, num_solution * 2, 1.0, r[0].re, local_size, cnjCr, stride_C, 0.0, tmp[0].re, local_size * 2);
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, local_size, num_solution, num_solution * 2, 1.0, r[0].re, local_size, rotCr, stride_C, 0.0, tmp[0].im, local_size * 2);

    if (is_use_p) {
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, local_size, num_solution, num_solution * 2, 1.0, p[0].re, local_size, cnjCp, stride_C, 1.0, tmp[0].re, local_size * 2);
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, local_size, num_solution, num_solution * 2, 1.0, p[0].re, local_size, rotCp, stride_C, 1.0, tmp[0].im, local_size * 2);
    }

    //P_next = TMP;
    for (int k = 0; k < num_solution; ++k) {
        SoAC::Copy(p[k], tmp[k], local_size);
    }

    //TMP = TMP + X * Cx;
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, local_size, num_solution, num_solution * 2, 1.0, x[0].re, local_size, cnjCx, stride_C, 1.0, tmp[0].re, local_size * 2);
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, local_size, num_solution, num_solution * 2, 1.0, x[0].re, local_size, rotCx, stride_C, 1.0, tmp[0].im, local_size * 2);
    /*
    //X_next = TMP;
    for (int k = 0; k < num_solution; ++k) {
        SoAC::Copy(x[k], tmp[k], local_size);
    }
    */

    int num_gemm_call = is_use_p ? 6 : 4;
    return num_gemm_call;
}

inline
void LoadSMatrix_d(int N, OneComplex* Sb, OneComplex* src) {
	const int n3 = N * 3;
	const int ldd = n3;

	OneComplex* XtX = Sb;
	OneComplex* XtP = Sb + N * 2 * n3;
	OneComplex* PtP = Sb + N * (2 * n3 + 2);
	CopySubMat(N, N, XtX, n3, src, N);
	CopySubMat(N, N, XtP, n3, src + N * N, N);
	CopySubMat(N, N, PtP, n3, src + N * N * 2, N);

}

inline
void SaveSMatrix_d(int N, OneComplex* dest, const OneComplex* Sb) {
	const int n3 = N * 3;
	const int ldd = n3;

	const OneComplex* XtX = Sb;
	const OneComplex* XtP = Sb + N * 2 * n3;
	const OneComplex* PtP = Sb + N * (2 * n3 + 2);
	CopySubMat(N, N, dest, N, XtX, n3);
	CopySubMat(N, N, dest + N * N, N, XtP, n3);
	CopySubMat(N, N, dest + N * N * 2, N, PtP, n3);

}


#endif


/*
X = {x[0], x[1], ..., x[num_solution-1]}
X = X * Cx + R * Cr
P = R * Cr
AX = AX * Cx + AR * Cr
AP = AR * Cr
*/

inline
int NextState_first_z(int local_size, int num_solution, const OneComplex* CxCrCp,
	SoAComplex* x, SoAComplex* r, SoAComplex* p,
	SoAComplex* Ax, SoAComplex* Ar, SoAComplex* Ap,
	SoAComplex* tmp) {


	const int n2 = num_solution * 2;
	{
		//P_next = Cr * R
		for (int k = 0; k < num_solution; ++k) {
			SoAC::SetZero(p[k], local_size);
			const OneComplex* cr = CxCrCp + k * n2 + num_solution;
			for (int j = 0; j < num_solution; ++j) {
				SoAC::AddV(p[k], r[j], cr[j], local_size);
			}
		}

		auto& tmp = r;//rはこの時点で使わないので作業領域にする//
		//TMP = X * Cx + P_next;
		for (int k = 0; k < num_solution; ++k) {
			SoAC::Copy(tmp[k], p[k], local_size);
			const OneComplex* cx = CxCrCp + k * n2;
			for (int j = 0; j < num_solution; ++j) {
				SoAC::AddV(tmp[k], x[j], cx[j], local_size);
			}
		}

		//X_next = TMP				
		for (int k = 0; k < num_solution; ++k) {
			SoAC::Copy(x[k], tmp[k], local_size);
		}
	}


	{
		//AP_next = AR * Cr
		for (int k = 0; k < num_solution; ++k) {
			SoAC::SetZero(Ap[k], local_size);
			const OneComplex* cr = CxCrCp + k * n2 + num_solution;
			for (int j = 0; j < num_solution; ++j) {
				SoAC::AddV(Ap[k], Ar[j], cr[j], local_size);
			}
		}

		auto& tmp = r;//rはこの時点で使わないので作業領域にする//
		//TMP = AX * Cx + AP_next;
		for (int k = 0; k < num_solution; ++k) {
			SoAC::Copy(tmp[k], Ap[k], local_size);
			const OneComplex* cx = CxCrCp + k * n2;
			for (int j = 0; j < num_solution; ++j) {
				SoAC::AddV(tmp[k], Ax[j], cx[j], local_size);
			}
		}

		//AX_next = TMP				
		for (int k = 0; k < num_solution; ++k) {
			SoAC::Copy(Ax[k], tmp[k], local_size);
		}
	}

	return 8;
}

/*
X = {x[0], x[1], ..., x[num_solution-1]}
X = X * Cx + R * Cr + P * Cp
P = R * Cr + P * Cp
AX = AX * Cx + AR * Cr + AP * Cp
AP = AR * Cr + AP * Cp
*/
inline
int NextState_z(int local_size, int num_solution, const OneComplex* CxCrCp,
	SoAComplex* x, SoAComplex* r, SoAComplex* p,
	SoAComplex* Ax, SoAComplex* Ar, SoAComplex* Ap,
	SoAComplex* tmp) {

	const int n3 = num_solution * 3;
	//TMP = R * Cr + P * Cp;
	for (int k = 0; k < num_solution; ++k) {
		const OneComplex* cx = CxCrCp + k * n3;
		const OneComplex* cr = CxCrCp + k * n3 + num_solution;
		const OneComplex* cp = CxCrCp + k * n3 + num_solution * 2;
		SoAC::SetZero(tmp[k], local_size);
		for (int j = 0; j < num_solution; ++j) {
			SoAC::AddV(tmp[k], r[j], cr[j], local_size);
			SoAC::AddV(tmp[k], p[j], cp[j], local_size);
		}
	}

	//P_next = TMP;
	for (int k = 0; k < num_solution; ++k) {
		SoAC::Copy(p[k], tmp[k], local_size);
	}

	//TMP = TMP + X * Cx;
	for (int k = 0; k < num_solution; ++k) {
		const OneComplex* cx = CxCrCp + k * n3;
		const OneComplex* cr = CxCrCp + k * n3 + num_solution;
		const OneComplex* cp = CxCrCp + k * n3 + num_solution * 2;
		for (int j = 0; j < num_solution; ++j) {
			SoAC::AddV(tmp[k], x[j], cx[j], local_size);
		}
	}

	//X_next = TMP;
	for (int k = 0; k < num_solution; ++k) {
		SoAC::Copy(x[k], tmp[k], local_size);
	}


	//TMP = AR * Cr + AP * Cp;
	for (int k = 0; k < num_solution; ++k) {
		const OneComplex* cx = CxCrCp + k * n3;
		const OneComplex* cr = CxCrCp + k * n3 + num_solution;
		const OneComplex* cp = CxCrCp + k * n3 + num_solution * 2;
		SoAC::SetZero(tmp[k], local_size);
		for (int j = 0; j < num_solution; ++j) {
			SoAC::AddV(tmp[k], Ar[j], cr[j], local_size);
			SoAC::AddV(tmp[k], Ap[j], cp[j], local_size);
		}
	}

	//AP_next = TMP;
	for (int k = 0; k < num_solution; ++k) {
		SoAC::Copy(Ap[k], tmp[k], local_size);
	}

	//TMP = TMP + AX * Cx;
	for (int k = 0; k < num_solution; ++k) {
		const OneComplex* cx = CxCrCp + k * n3;
		const OneComplex* cr = CxCrCp + k * n3 + num_solution;
		const OneComplex* cp = CxCrCp + k * n3 + num_solution * 2;
		for (int j = 0; j < num_solution; ++j) {
			SoAC::AddV(tmp[k], Ax[j], cx[j], local_size);
		}
	}

	//AX_next = TMP;
	for (int k = 0; k < num_solution; ++k) {
		SoAC::Copy(Ax[k], tmp[k], local_size);
	}

	return 12;
}


/*
固有値問題を解くLOBPCG

問題点、iter_max==1でも前回の履歴からr,p,Ar,Apを継続させるべき

*/
template <class OperationA, class WATCH>
inline
void Eigen_z_multi_mpi(const GridRangeMPI& l_grid, const double dVol, const int num_solution, const int iter_max,
	double* eigen_values, SoAComplex* eigen_vectors, double* keep, double* work,
	OneComplex* keep_S_matrix, int SCF_step, int sk,
#ifdef USE_SCALAPACK
	const BlacsGridInfo& blacs_grid,
#endif
	OperationA OpeA, WATCH& watch) {

	using namespace vecmath;

	const int proc_id = GetProcessID(l_grid.mpi_comm);
	const bool is_root = (proc_id == 0);
	const size_t local_size = l_grid.Size3D();

	auto IntervalPointers = [](double* a, size_t N, size_t num_solution) {
		std::vector<SoAComplex> pointers(num_solution);
		for (size_t i = 0; i < num_solution; ++i) {
			pointers[i].re = a + (i * 2) * N;
			pointers[i].im = a + (i * 2 + 1) * N;
		}
		return pointers;
		};


	auto x = eigen_vectors;
	auto r = IntervalPointers(work, local_size, num_solution);
	auto p = IntervalPointers(keep, local_size, num_solution);
	auto Ax = IntervalPointers(work + local_size * num_solution * 2, local_size, num_solution);
	auto Ar = IntervalPointers(work + 2 * local_size * num_solution * 2, local_size, num_solution);
	auto Ap = IntervalPointers(work + 3 * local_size * num_solution * 2, local_size, num_solution);
	auto tmpM = IntervalPointers(work + 4 * local_size * num_solution * 2, local_size, num_solution);


	const int n3 = 3 * num_solution;
	const int nn9 = (3 * num_solution) * (3 * num_solution);
	OneComplex* Sa = new OneComplex[nn9 * 4];
	OneComplex* Sb = Sa + nn9;
	OneComplex* red_Sa = Sa + nn9*2;
	OneComplex* red_Sb = Sa + nn9*3;
	SetZero<double>((double*)Sa, nn9 * 4 * 2);



	const bool is_first_time = (SCF_step == 0);
	bool is_loaded_S_matrix = false;
	constexpr int S_MAT_REST_LIMIT = 10;//S行列はロードだけだと数値誤差が積もるので定期的に更新//
#ifdef S_MATRIX_RELOAD
	//数値誤差が積もって発散する//
	if ((keep_S_matrix != nullptr) && (!is_first_time) && (SCF_step % S_MAT_REST_LIMIT != 0)) {
		LoadSMatrix_d(num_solution, Sb, keep_S_matrix);
		is_loaded_S_matrix = true;
	}
#endif


	
	double* S_e_value = new double[(3 * num_solution)];
	OneComplex* S_e_vector = new OneComplex[nn9];

	const double limit = 1.0e-10;

	double* l_norms = new double[num_solution * 6];
	double* l_norms_p = l_norms + num_solution;
	double* l_norms_r = l_norms + num_solution * 2;
	double* sum_norms = l_norms + num_solution * 3;
	double* sum_norms_p = l_norms + num_solution * 4;
	double* sum_norms_r = l_norms + num_solution * 5;

#ifdef USE_BLAS
	double* temp_mat33 = new double[nn9*4];
#endif

	if (is_root) {
		printf("Eigen solver with LOBPCG_z_multi_mpi\n");
	}


	for (int k = 0; k < num_solution; ++k) {
		SoAC::MulC(x[k], sqrt(dVol), local_size);
	}

	//calculte initial redidual r//
	//Because Hamiltonian operator A was changed by changing electronic density rho, //
	// Ax and Ap should be re-calculated//
	for (int k = 0; k < num_solution; ++k) {
		OpeA(Ax[k], x[k]); // calculate Ax from x 	
		OpeA(Ap[k], p[k]);
	}

	bool is_large_difference = false;
	double*& l_ei = l_norms;
	double*& ei = eigen_values;
	for (int k = 0; k < num_solution; ++k) {
		auto xAx = SoAC::InnerProd(x[k], Ax[k], local_size);
		l_ei[k] = xAx.r;
	}
	watch.Record(10);
	MPI_Allreduce(l_ei, ei, num_solution, MPI_DOUBLE, MPI_SUM, l_grid.mpi_comm);
	watch.Record(14);
	
	for (int k = 0; k < num_solution; ++k) {
		double norm = 0.0;
		for (size_t i = 0; i < local_size; ++i) {
			r[k].re[i] = (Ax[k].re[i] - ei[k] * x[k].re[i]);
			r[k].im[i] = (Ax[k].im[i] - ei[k] * x[k].im[i]);
			norm += r[k].re[i] * r[k].re[i] + r[k].im[i] * r[k].im[i];
		}

		l_norms[k] = norm;
	}

	for (int k = 0; k < num_solution; ++k) {
		if (is_first_time) {
			SoAC::SetZero(p[k], local_size);
			SoAC::SetZero(Ap[k], local_size);
		}
#ifdef  LOBPCG_DEBUG_PRINT
		if (is_root) {
			printf("(0, %d) ei, norm_r: %f, %f\n", k, ei[k], sum_norms[k] );
		}
#endif
	}

	if (iter_max < 1) {
		//iterationしないのでnormalizeを元に戻して終了//
		for (int k = 0; k < num_solution; ++k) {
			SoAC::MulC(x[k], 1.0 / sqrt(dVol), local_size);
		}
		return;
	}


	watch.Record( 10);
//main loop of LOBPCG
	

	int iter;
	for (iter = 0; iter < iter_max; ++iter) {


		if (is_root) {
			printf("  LOBPCG iteration: %d\n", iter + 1); fflush(stdout);
		}

		

		for (int k = 0; k < num_solution; ++k) {
			OpeA(Ar[k], r[k]);
		}



		bool required_reduce_H_S = true;

#ifdef USE_BLAS
#ifdef SELF_TRANSPOSE
		int num_gemm1_call = 0;
		if (iter == 0) {
			num_gemm1_call = MakeMatrix_z_blas(local_size, num_solution, Sa, Sb, temp_mat33, &x[0], &r[0], &p[0], &Ax[0], &Ar[0], &Ap[0], &tmpM[0], is_loaded_S_matrix, false, false);
		} else if ((iter % S_MAT_REST_LIMIT == 0) || is_large_difference) {
			num_gemm1_call = MakeMatrix_z_blas(local_size, num_solution, Sa, Sb, temp_mat33, &x[0], &r[0], &p[0], &Ax[0], &Ar[0], &Ap[0], &tmpM[0], false, false, false);
			is_large_difference = false;
		} else {
			num_gemm1_call = MakeMatrix_z_blas(local_size, num_solution, Sa, Sb, temp_mat33, &x[0], &r[0], &p[0], &Ax[0], &Ar[0], &Ap[0], &tmpM[0], true, true, true);

			required_reduce_H_S = false;

		}
#else
		int num_gemm1_call = 0;
		if (iter == 0) {
			num_gemm1_call = MakeMatrix_z_blas(local_size, num_solution, Sa, Sb, temp_mat33, &x[0], &r[0], &p[0], &Ax[0], &Ar[0], &Ap[0], is_loaded_S_matrix, false, false);
		} else if ((iter % S_MAT_REST_LIMIT == 0) || is_large_difference) {
			num_gemm1_call = MakeMatrix_z_blas(local_size, num_solution, Sa, Sb, temp_mat33, &x[0], &r[0], &p[0], &Ax[0], &Ar[0], &Ap[0], false, false, false);
			is_large_difference = false;
		} else {
			num_gemm1_call = MakeMatrix_z_blas(local_size, num_solution, Sa, Sb, temp_mat33, &x[0], &r[0], &p[0], &Ax[0], &Ar[0], &Ap[0], true, true, true);			

			required_reduce_H_S = false;

		}
#endif
#else
		num_gemm1_call=MakeMatrix_z(local_size, num_solution, Sa, Sb, &x[0], &r[0], &p[0], &Ax[0], &Ar[0], &Ap[0]);
#endif
		watch.Record( 17, num_gemm1_call);



		//functions to compress H matrix into the empty region of S matrix and to shurink the size of MPI_Allreduce//
		const int N = num_solution;
		const int n2 = num_solution * 2;
		auto Compress_S_H = [&]() {
			const auto* XtHX = Sa;
			const auto* XtHR = Sa + N * n3;
			auto* dstPtR = Sb + N * (n3 + 2);
			auto* dstRtX = Sb + N;
			CopySubMat(N, N, dstPtR, n3, XtHX, n3);
			CopySubMat(N, n2, dstRtX, n3, XtHR, n3);

			};
		auto Expand_S_H = [&]() {
			//repair S matrix//
			auto* dstXtHX = red_Sa;
			auto* dstXtHR = red_Sa + N * n3;
			const auto* PtR = red_Sb + N * (n3 + 2);
			const auto* RtX = red_Sb + N;
			CopySubMat(N, N, dstXtHX, n3, PtR, n3);
			CopySubMat(N, n2, dstXtHR, n3, RtX, n3);

			};


#ifdef CHECK_NEXT_VECTOR
        const bool is_use_p = !( (is_first_time && (iter == 0)) || is_large_difference);
#else
		const bool is_use_p = !(is_first_time && (iter == 0));
#endif

		if (!is_use_p) {

			Compress_S_H();
			watch.Record(10);
#ifdef USE_SCALAPACK
			MPI_Allreduce(Sb, red_Sb, (N* N * 2) * 6, MPI_DOUBLE, MPI_SUM, l_grid.mpi_comm);
#else
			MPI_Reduce(Sb, red_Sb, (N * N * 2) * 6, MPI_DOUBLE, MPI_SUM, 0, l_grid.mpi_comm);
#endif
			//Expand_S_H();
			watch.Record(14);

			
			//初回はp=0のため、3n x 3n領域を2n x 2n行列に縮小//
			// 通信ではred_Sbのみを使ったので、red_Sa(9nn領域中の8nn領域)
			//3x3領域を2x2領域に見せるためにとして,初回はPtRからPtPまでの4blockを利用する
			OneComplex* Sa2 = red_Sa;
			OneComplex* Sb2 = Sa2 + n2 * n2;

			CopySubMat(N, n2, Sa2 + N*n2, n2, red_Sb + N, n3);
			CopySubMat(N, N, Sa2, n2, red_Sb + N*(n3+2), n3);
			CopySubMat(N, n2, Sb2 + N * n2, n2, red_Sb + N*n3, n3);
			CopySubMat(N, N, Sb2, n2, red_Sb, n3);


			
			{//store reduced H and S matrix
				CopySubMat(N, N, Sa, n3, red_Sb + N * (n3 + 2), n3);
				CopySubMat(N, n2, Sa + N*n3, n3, red_Sb + N, n3);
				CopySubMat(N, N, Sb, n3, red_Sb, n3);
				CopySubMat(N, n2, Sb + N * n3, n3, red_Sb + N*n3, n3);
			}


			watch.Record(10);

#ifdef USE_SCALAPACK
			int info = lapack_PZHEGVX_split(Sa2, Sb2, &(S_e_value[0]), &(S_e_vector[0]), 2 * num_solution, l_grid.mpi_comm, blacs_grid);
#else
			if (is_root) {
				int info = lapack_ZHEGVX(Sa2, Sb2, &(S_e_value[0]), &(S_e_vector[0]), 2 * num_solution);
			}
#endif
			watch.Record(15);


			MPI_Bcast(S_e_value, num_solution, MPI_DOUBLE, 0, l_grid.mpi_comm);
			MPI_Bcast(S_e_vector, num_solution * n2 * 2, MPI_DOUBLE, 0, l_grid.mpi_comm);
			watch.Record(14);
			//const int min_id = (S_e_value[0] < S_e_value[1]) ? 0 : 1;
			//固有値の小さい順に並んでいるもの考えてよい//


		} else {//(is_use_p==true)//


			auto StoreReducedMatrix33 = [](int N, OneComplex* Sb, const OneComplex* red_Sb) {
				auto n2 = N * 2;
				auto n3 = N * 3;

				OneComplex* XtX = Sb;
				OneComplex* XtR = Sb + N * n3;
				OneComplex* XtP = Sb + N * 2 * n3;
				const OneComplex* red_XtX = red_Sb;
				const OneComplex* red_XtR = red_Sb + N * n3;
				const OneComplex* red_XtP = red_Sb + N * 2 * n3;

				CopySubMat(N, N, XtX, n3, red_XtX, n3);
				CopySubMat(N, n2, XtR, n3, red_XtR, n3);
				CopySubMat(N, n3, XtP, n3, red_XtP, n3);;
				};
			if (required_reduce_H_S) {
				Compress_S_H();
				watch.Record(10);
#ifdef USE_SCALAPACK
				MPI_Allreduce(Sa + N * n3 * 2, red_Sa + N * n3 * 2, (N* N * 2) * 12, MPI_DOUBLE, MPI_SUM, l_grid.mpi_comm);
#else
				MPI_Reduce(Sa + N * n3 * 2, red_Sa + N * n3 * 2, (N* N * 2) * 12, MPI_DOUBLE, MPI_SUM, 0, l_grid.mpi_comm);
#endif
				Expand_S_H();
				watch.Record(14);

				StoreReducedMatrix33(num_solution, Sa, red_Sa);
				StoreReducedMatrix33(num_solution, Sb, red_Sb);
				

			} else {
				//この場合に通信するのは,XtHR,RtHR,RtHP,RtRの4つ//
				//ただしDIRECT_S_MATRIX_XR_RPが未定義の場合はXtR,RtPも通信が必要//

				const size_t offset_XtR = N * n3;
				const size_t offset_RtP = N * (2 * n3 + 1);
				const size_t offset_PtR = N * (n3 + 2);

				//reduce XtHR,RtHR,RtHP//
				CopySubMat(N, N, Sa + offset_PtR, n3, Sa + offset_RtP, n3);
				watch.Record(10);
				MPI_Allreduce(Sa + offset_XtR, red_Sa + offset_XtR, (N* N * 2) * 3, MPI_DOUBLE, MPI_SUM, l_grid.mpi_comm);
				watch.Record(14);
				auto* RtHP2 = red_Sa + N * (2 * n3 + 1);
				auto* PtHR2 = red_Sa + N * (n3 + 2);
				CopySubMat(N, N, red_Sa + offset_RtP, n3, red_Sa + offset_PtR, n3);

				//store for direct calculation of matrix//
				CopySubMat(N, n2, Sa + offset_XtR, n3, red_Sa + offset_XtR, n3);
				CopySubMat(N, N, Sa + offset_RtP, n3, red_Sa + offset_RtP, n3);

				//load from direct calculated matrix//
				const size_t offset_XtP = N * 2 * n3;
				const size_t offset_PtP = N * (2 * n3 + 2);
				CopySubMat(N, N, red_Sa, n3, Sa, n3);
				CopySubMat(N, N, red_Sa + offset_XtP, n3, Sa + offset_XtP, n3);
				CopySubMat(N, N, red_Sa + offset_PtP, n3, Sa + offset_PtP, n3);

#ifdef DIRECT_S_MATRIX_XR_RP
				//transfer only RtR, XtR, RtP //
				const size_t offset_RtR = N * (n3 + 1);
				//auto* RtR = Sb + N * (n3 + 1);
				CopySubMat(N, N, red_Sb, N, Sb + offset_RtR, n3);
				watch.Record(10); 
				MPI_Allreduce(red_Sb, red_Sb + N * N, (N* N * 2) * 1, MPI_DOUBLE, MPI_SUM, l_grid.mpi_comm);
				watch.Record(14);
				//auto* RtR2 = red_Sb + N * (n3 + 1);
				CopySubMat(N, N, red_Sb + offset_RtR, n3, red_Sb + N * N, N);

				//store for direct calculation of matrix//
				CopySubMat(N, N, Sb + offset_RtR, n3, red_Sb + N * N, N);

				//load from direct calculated matrix//
				CopySubMat(N, N, red_Sb + offset_XtR, n3, Sb + offset_XtR, n3);
				CopySubMat(N, N, red_Sb + offset_RtP, n3, Sb + offset_RtP, n3);

#else
				//transfer RtR //
				auto* RtP = Sb + N * (2 * n3 + 1);
				auto* PtR = Sb + N * (n3 + 2);
				CopySubMat(N, N, PtR, n3, RtP, n3);
				watch.Record(10);
				MPI_Allreduce(Sb + N * n33n, red_Sb + N * n3, (N* N * 2) * 3, MPI_DOUBLE, MPI_SUM, l_grid.mpi_comm);
				watch.Record(14);
				auto* RtP2 = red_Sb + N * (2 * n3 + 1);
				auto* PtR2 = red_Sb + N * (n3 + 2);
				CopySubMat(N, N, RtP2, n3, PtR2, n3);

				//store for direct calculation of matrix//
				CopySubMat(N, N * 2, Sb + N * n3, n3, red_Sb + N * n3, n3);
				CopySubMat(N, N, RtP, n3, RtP2, n3);
#endif

				//load from direct calculated matrix//
				CopySubMat(N, N, red_Sb, n3, Sb, n3);
				CopySubMat(N, N, red_Sb + offset_XtP, n3, Sb + offset_XtP, n3);
				CopySubMat(N, N, red_Sb + offset_PtP, n3, Sb + offset_PtP, n3);

				
			}

			watch.Record(10);
#ifdef USE_SCALAPACK
			int info = lapack_PZHEGVX_split(&(red_Sa[0]), &(red_Sb[0]), &(S_e_value[0]), &(S_e_vector[0]), 3 * num_solution, l_grid.mpi_comm, blacs_grid);
#else
			if (is_root) {
				int info = lapack_ZHEGVX(&(red_Sa[0]), &(red_Sb[0]), &(S_e_value[0]), &(S_e_vector[0]), 3 * num_solution);
			}
#endif
			watch.Record(15);



			
			MPI_Bcast(S_e_value, num_solution, MPI_DOUBLE, 0, l_grid.mpi_comm);
			MPI_Bcast(S_e_vector, num_solution *n3 * 2, MPI_DOUBLE, 0, l_grid.mpi_comm);
			watch.Record( 14);
			//固有値の小さい順に並んでいるもの考えてよい//
			//const int min_id = (S_e_value[0] < S_e_value[1]) ? (S_e_value[0] < S_e_value[2]) ? 0 : 2 : (S_e_value[1] < S_e_value[2]) ? 1 : 2;
		
		}


		const bool is_final_step = (iter + 1 == iter_max);

		/*
		X = {x[0], x[1], ..., x[num_solution-1]}
		X = X * Cx + R * Cr + P * Cp
		P = R * Cr + P * Cp
		AX = AX * Cx + AR * Cr + AP * Cp
		AP = AR * Cr + AP * Cp
		where Cx,Cr,CP which are solution S_e_vector of eigen value problem
		*/
		
#ifdef CHECK_NEXT_VECTOR
        int num_gemm2_call = NextVector_z2(local_size, num_solution, S_e_vector, is_use_p, &x[0], &r[0], &p[0], &tmpM[0], (OneComplex*)temp_mat33);
        watch.Record(16, num_gemm2_call);
        //check norm, here.
        
        for (int j = 0; j < num_solution; ++j) {
            l_norms[j] = SoAC::Norm(tmpM[j], local_size);
            l_norms_p[j] = SoAC::Norm(p[j], local_size);
        }
        watch.Record(10);

        //const int norm_size = is_final_step ? num_solution * 2 : num_solution * 3;
        const int norm_size = num_solution * 2;
        MPI_Allreduce(l_norms, sum_norms, norm_size, MPI_DOUBLE, MPI_SUM, l_grid.mpi_comm);
        watch.Record(14);

#ifdef LOBPCG_Z_PRINT_NORM
        //show norm
        {
            if (is_root) {
                for (int j = 0; j < norm_size; ++j) {
                    printf("[%d]sk=%d, norm[%d] = %g, %g\n", proc_id, sk, j, sum_norms[j], 1.0 / sqrt(sum_norms[j]));
                }
                fflush(stdout);
            }
        }
#endif

        //正しく解けた場合はxがnormalizeされている//
        //これがずれた場合は数値誤差が積もっているのでXをrollbackし、//
        //Escapeするべし(Hを更新して次のSCF-loopで再計算)//
        double max_diff_norm_x = 0.0;
        for (int k = 0; k < num_solution; ++k) {
            double diff = fabs(1.0 - sum_norms[k]);
            if (max_diff_norm_x < diff)max_diff_norm_x = diff;
        }
        uint8_t check_flag = 0;
        if (max_diff_norm_x > LIMIT_Z_NORM_X_FOR_RESET_MATRIX) {
            check_flag = 1;
        }

        uint8_t check_flag_sum = 0;
        MPI_Allreduce(&check_flag, &check_flag_sum, 1, MPI_UINT8_T, MPI_BOR, l_grid.mpi_comm);
        is_large_difference = (check_flag_sum == 1);

        if(is_large_difference){
            if (is_root) {
                printf("LOBPCG numerical stability breaking: | 1 - |psi|^2| = %g > %g\n", max_diff_norm_x, LIMIT_Z_NORM_X_FOR_RESET_MATRIX);
                printf("Rollback psi to the previous state, and escape from LOBPCG loop\n");
                fflush(stdout);
            }

            //////////////////////////////////////
            //escape from this loop
            //////////////////////////////////////
            //keep previous X and R, and H and S matrix is completely calculated in next step and clear P//

            
            //idea copy R to P for next steps//
            for (int k = 0; k < num_solution; ++k) {
                SoAC::Copy(p[k], r[k], local_size);
            }
            
            //scaling for dVol, where previous X is alearedy normalized to 1 in grid scale//
            const double coef_x = 1.0 / sqrt(dVol);
            for (int j = 0; j < num_solution; ++j) {
                SoAC::MulC(x[j], coef_x, local_size);
            }


            break;            

        }else{//(!is_large_difference)
 
            //X_next = TMP;
            for (int k = 0; k < num_solution; ++k) {
                SoAC::Copy(x[k], tmpM[k], local_size);
            }


            for (int i = 0; i < num_solution; ++i) {
                ei[i] = S_e_value[i];
            }


#ifdef EVERY_NORMALIZE_P
            //normalize p to supress nan//
            {
                for (int j = 0; j < num_solution; ++j) {
                    const double coef_p = 1.0 / sqrt(sum_norms_p[j]);
                    sum_norms_p[j] = coef_p;
                    SoAC::MulC(p[j], coef_p, local_size);
                }

                for (int j = 0; j < num_solution; ++j) {
                    SoAC::MulC(Ap[j], sum_norms_p[j], local_size);
                }

            }
#endif


            if (!is_final_step) {

                watch.Record(10);
                int num_gemm2_call = NextVector_z2(local_size, num_solution, S_e_vector, is_use_p, &Ax[0], &Ar[0], &Ap[0], &tmpM[0], (OneComplex*)temp_mat33);

                //AX_next = TMP;
                for (int k = 0; k < num_solution; ++k) {
                    SoAC::Copy(Ax[k], tmpM[k], local_size);
                }
                watch.Record(16, num_gemm2_call);



                for (int j = 0; j < num_solution; ++j) {
                    //double norm = 0.0;
                    for (size_t i = 0; i < local_size; ++i) {
                        r[j].re[i] = (Ax[j].re[i] - ei[j] * x[j].re[i]);
                        r[j].im[i] = (Ax[j].im[i] - ei[j] * x[j].im[i]);
                        //norm += r[j].re[i] * r[j].re[i] + r[j].im[i] * r[j].im[i];
                    }

                    //l_norms_r[j] = norm;
                }


                watch.Record(10);
                int num_gemm3_call = 0;
                num_gemm3_call = NextMatrix_z(num_solution, Sa, Sb, S_e_vector, S_e_value, is_use_p, !is_final_step);
                watch.Record(19, num_gemm3_call);

#ifdef EVERY_NORMALIZE_P
                //already not required//
                //normalize XtX ... in H- and S-matrix//
                ScaleHSMatrix_z(num_solution, Sa, Sb, nullptr, sum_norms_p, nullptr, !is_final_step);
#endif
                
            } else {//(is_final_step) 

                for (int j = 0; j < num_solution; ++j) {
                    const double coef_x = 1.0 / sqrt(sum_norms[j] * dVol);
                    SoAC::MulC(x[j], coef_x, local_size);
                }


            }


        

#ifdef LOBPCG_Z_PRINT_EIGEN
            if (is_root) {
                printf("LOBPCG iter = %d\n", iter + 1);
                if (is_final_step) {
                    for (int i = 0; i < std::min(20, num_solution); ++i) {
                        printf("eigen[%d] = %f\n", i, ei[i]);
                    }
                }
                else {
                    for (int i = 0; i < std::min(20, num_solution); ++i) {
                        printf("eigen[%d] = %f, %f\n", i, ei[i], S_e_value[i]);
                    }
                    /*
                    for (int i = std::max(20, num_solution / 2); i < std::min(num_solution / 2 + 20, num_solution); ++i) {
                        printf("eigen[%d] = %f, %f\n", i, ei[i], S_e_value[i]);
                    }
                    for (int i = std::max(40, num_solution - 20); i < num_solution; ++i) {
                        printf("eigen[%d] = %f, %f\n", i, ei[i], S_e_value[i]);
                    }
                    */
                }
                fflush(stdout);
            }

#endif

            watch.Record(10);
        }

#else      //!CHECK_NEXT_VECTOR

        int num_gemm2_call = 0;
#ifdef USE_BLAS
        num_gemm2_call = NextVector_z(local_size, num_solution, S_e_vector, is_use_p, !is_final_step, &x[0], &r[0], &p[0], &Ax[0], &Ar[0], &Ap[0], &tmpM[0], (OneComplex*)temp_mat33);
#else
		if (is_use_p) {
			num_gemm2_call = NextState_z(local_size, num_solution, S_e_vector, &x[0], &r[0], &p[0], &Ax[0], &Ar[0], &Ap[0], &tmpM[0]);
		} else {
			num_gemm2_call = NextState_first_z(local_size, num_solution, S_e_vector, &x[0], &r[0], &p[0], &Ax[0], &Ar[0], &Ap[0], &tmpM[0]);
		}
#endif

        watch.Record(16, num_gemm2_call);




		for (int i = 0; i < num_solution; ++i) {
			ei[i] = S_e_value[i];
		}


		if (!is_final_step) {
			for (int j = 0; j < num_solution; ++j) {
				double norm = 0.0;
				for (size_t i = 0; i < local_size; ++i) {
					r[j].re[i] = (Ax[j].re[i] - ei[j] * x[j].re[i]);
					r[j].im[i] = (Ax[j].im[i] - ei[j] * x[j].im[i]);
					norm += r[j].re[i] * r[j].re[i] + r[j].im[i] * r[j].im[i];
				}

				l_norms_r[j] = norm;
			}
		}

		for (int j = 0; j < num_solution; ++j) {
			l_norms[j] = SoAC::Norm(x[j], local_size);
			l_norms_p[j] = SoAC::Norm(p[j], local_size);
		}
		watch.Record(10);

		//const int norm_size = is_final_step ? num_solution * 2 : num_solution * 3;
        const int norm_size = num_solution * 3;
		MPI_Allreduce(l_norms, sum_norms, norm_size, MPI_DOUBLE, MPI_SUM, l_grid.mpi_comm);
		watch.Record(14);
		

#ifdef LOBPCG_Z_PRINT_NORM
		//show norm
		{
			if (is_root) {
				for (int j = 0; j < norm_size; ++j) {
					printf("[%d]sk=%d, norm[%d] = %g, %g\n", proc_id, sk, j, sum_norms[j], 1.0 / sqrt(sum_norms[j]));
				}
				fflush(stdout);
			}
		}
#endif

		//正しく解けた場合はxがnormalizeされている//
		//これがずれた場合は数値誤差が積もっているのでH,S行列を再計算するべし//
		double max_diff_norm_x = 0.0;
		for (int k = 0; k < num_solution; ++k) {
			double diff = fabs(1.0 - sum_norms[k]);
			if (max_diff_norm_x < diff)max_diff_norm_x= diff;
		}
		if (max_diff_norm_x > LIMIT_Z_NORM_X_FOR_RESET_MATRIX) {
			is_large_difference = true;

            if (is_root) {
                printf("LOBPCG | 1 - |psi|^2| = %g > %g\n", max_diff_norm_x, LIMIT_Z_NORM_X_FOR_RESET_MATRIX);
                fflush(stdout);
            }
		}
        //test////////////////////////////////////
        //is_large_difference = true;
        ////////////////////////////////////test//

        bool is_escape = false;
        double max_diff_norm_r = 0.0;
        for (int k = 0; k < num_solution; ++k) {            
            if (max_diff_norm_r < sum_norms_r[k])max_diff_norm_r = sum_norms_r[k];
        }
        if (max_diff_norm_r < LIMIT_Z_NORM_R_FOR_ESCAPE) {
            is_escape = true;
        }

#ifdef EVERY_NORMALIZE_P
        //normalize p to supress nan//
        {
            for (int j = 0; j < num_solution; ++j) {
                const double coef_p = 1.0 / sqrt(sum_norms_p[j]);
                sum_norms_p[j] = coef_p;
                SoAC::MulC(p[j], coef_p, local_size);
            }
            
            for (int j = 0; j < num_solution; ++j) {
                SoAC::MulC(Ap[j], sum_norms_p[j], local_size);
            }

        }
#endif

		if (is_final_step || is_escape){

			for (int j = 0; j < num_solution; ++j) {
				const double coef_x = 1.0 / sqrt(sum_norms[j] * dVol);
				SoAC::MulC(x[j], coef_x, local_size);
			}


		}else{// !is_final_step //
            if (is_large_difference) {
                //数値誤差によるnanの回避
                for (int j = 0; j < num_solution; ++j) {
                    const double coef_x = 1.0 / sqrt(sum_norms[j]);
                    SoAC::MulC(x[j], coef_x, local_size);
                    SoAC::MulC(Ax[j], coef_x, local_size);
                    SoAC::MulC(r[j], coef_x, local_size);
                }
#ifndef EVERY_NORMALIZE_P
                for (int j = 0; j < num_solution; ++j) {
                    const double coef_p = 1.0 / sqrt(sum_norms_p[j]);
                    SoAC::MulC(p[j], coef_p, local_size);
                    SoAC::MulC(Ap[j], coef_p, local_size);
                }
#endif

            }else{//(!is_large_difference) 
				watch.Record(10);
				int num_gemm3_call = 0;
				num_gemm3_call = NextMatrix_z(num_solution, Sa, Sb, S_e_vector, S_e_value, is_use_p, !is_final_step);
				watch.Record(19, num_gemm3_call);

#ifdef EVERY_NORMALIZE_P
                //already not required//
                //normalize XtX ... in H- and S-matrix//
                ScaleHSMatrix_z(num_solution, Sa, Sb, nullptr, sum_norms_p, nullptr, !is_final_step);
#endif
/*				
#if 0
                auto VanishDiagonalIm = [](int N, OneComplex* A, int lda) {
                    for (int i = 0; i < N; ++i) {
                        A[i + lda * i].i = 0.0;
                    }
                };
#else                
                auto VanishDiagonalIm = [](int N, OneComplex* A, int lda) {
                    for (int i = 0; i < N; ++i) {
                        for (int j = 0; j < N; ++j) {
                            const double ave_r = (A[j + lda * i].r + A[i + lda * j].r) * 0.5;
                            const double ave_i = (A[j + lda * i].i - A[i + lda * j].i) * 0.5;
                            A[j + lda * i].r = ave_r;
                            A[i + lda * j].r = ave_r;
                            A[j + lda * i].i = ave_i;
                            A[i + lda * j].i = -ave_i;
                        }
                    }
                    };
#endif                
                VanishDiagonalIm(num_solution, Sa, n3);
                VanishDiagonalIm(num_solution, Sa + N * (2*n3+2), n3);
                VanishDiagonalIm(num_solution, Sb, n3);
                VanishDiagonalIm(num_solution, Sb + N * (2 * n3 + 2), n3);
*/                
			}
		}



#ifdef LOBPCG_Z_PRINT_EIGEN
		if (is_root) {
			printf("LOBPCG iter = %d\n", iter + 1);
			if (is_final_step) {
				for (int i = 0; i < std::min(20, num_solution); ++i) {
					printf("eigen[%d] = %f\n", i, ei[i]);
				}
			} else {
				for (int i = 0; i < std::min(20, num_solution); ++i) {
					printf("eigen[%d] = %f, %f\n", i, ei[i], S_e_value[i]);
				}
				for (int i = std::max(20, num_solution / 2); i < std::min(num_solution / 2 + 20, num_solution); ++i) {
					printf("eigen[%d] = %f, %f\n", i, ei[i], S_e_value[i]);
				}
				for (int i = std::max(40, num_solution - 20); i < num_solution; ++i) {
					printf("eigen[%d] = %f, %f\n", i, ei[i], S_e_value[i]);
				}
			}
			fflush(stdout);
		}

#endif

		watch.Record(10);
		
        if (is_escape) {
            if (is_root) {
                printf("LOBPCG escape |residual|^2 = %g < %g\n", max_diff_norm_r, LIMIT_Z_NORM_R_FOR_ESCAPE);
                fflush(stdout);
            }
            break;
        }

#endif
	}

	
	////////////////////////////End of loop of Lanczos to create Symmetric Tridiagonal T// 

#ifdef S_MATRIX_RELOAD
	if ((keep_S_matrix != nullptr)) {
		int num_gemm3_call = 0;
		num_gemm3_call = NextMatrix_z(num_solution, Sa, Sb, S_e_vector, S_e_value, (iter_max>1), false);
		watch.Record(19, num_gemm3_call);
		//normalize XtX ... in H- and S-matrix//
		ScaleHSMatrix_z(num_solution, Sa, Sb, sum_norms, sum_norms_p, sum_norms_r, false);
		SaveSMatrix_d(num_solution, keep_S_matrix, Sb);
	}
#endif

	delete[] l_norms;
	delete[] Sa;
	
	delete[] S_e_value;
	delete[] S_e_vector;

#ifdef USE_BLAS
	delete[] temp_mat33;
#endif


	watch.Record( 10);
}

inline
size_t WorkSize_z_multi_mpi(int local_size, int num_solution) {	
	return 10 * local_size * num_solution;
}

}

#endif
