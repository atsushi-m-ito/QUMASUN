#pragma once

#include <stdio.h>
#include <stdlib.h>
#include "wrap_lapack.h"
#include "transpose.h"

#ifdef LAPACK_NO_HEADER
#ifndef lapack_complex_double
#define lapack_complex_double  std::complex<double>
#endif
extern "C" {
#ifdef _NEC
//#define dgemm
#endif
void dgemm_(const char* transa, const char* transb, const int* m, const int* n, const int* k,
		const double* alpha, const double* a, const int* lda, const double* b, const int* ldb,
		const double* beta, double* c, const int* ldc);
}
#endif




/*
* DGEMM wrapper
* C = alpha*(A^t*B) + beta*C
* matrix A is KxM grids (A[iy + K * ix] where ix in [0,M), iy in [0,K)
* matrix B is KxN grids (A[iy + K * ix] where ix in [0,N), iy in [0,K)
*/
inline
void blas_DGEMM_t(int M, int N, int K, double* A, double* B, double* C, int strideC, double alpha, double beta) {
	cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, M, N, K, alpha, A, K, B, K, beta, C, strideC);
	//cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, M, N, K, alpha, A, M, B, K, beta, C, strideC);
	//cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, M, N, K, alpha, A, K, B, K, beta, C, strideC);
	
}

inline
void blas_ZGEMM_t(int M, int N, int K, void* A, void* B, void* C, int strideC, OneComplex alpha, OneComplex beta) {
    cblas_zgemm(CblasColMajor, CblasTrans, CblasNoTrans, M, N, K, &alpha, A, K, B, K, &beta, C, strideC);
    
}

#ifdef USE_DGEMMT
inline
void blas_DGEMMT_t_up(int N, int K, double* A, double* B, double* C, int strideC, double alpha, double beta) {
	cblas_dgemmt(CblasColMajor, CblasUpper, CblasTrans, CblasNoTrans, N, K, alpha, A, K, B, K, beta, C, strideC);
	//cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, M, N, K, alpha, A, M, B, K, beta, C, strideC);
	//cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, M, N, K, alpha, A, K, B, K, beta, C, strideC);

}
#endif

inline
void blas_DGEMM_t2(int M, int N, int K, double* A, double* B, double* C, int strideC, double alpha, double beta) {
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, M, N, K, alpha, A, K, B, K, beta, C, strideC);
	
}


/*
* DGEMM wrapper
* C = alpha*(A^t*B) + beta*C
* matrix A is KxM grids (A[iy + K * ix] where ix in [0,M), iy in [0,K)
* matrix B is KxN grids (A[iy + K * ix] where ix in [0,N), iy in [0,K)
*/
inline
void blas_DGEMM_n(int M, int N, int K, double* A, double* B, double* C, int strideC, double alpha, double beta) {
	//cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, alpha, A, K, B, N, beta, C, strideC);
	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, M, N, K, alpha, A, M, B, K, beta, C, strideC);
}


inline
void nec_DMM_add(int N, int M, int K, double* A, double* B, double* C, int strideC, double alpha) {
#pragma _NEC ivdep
	for (int m = 0; m < M; ++m) {
#pragma _NEC ivdep
		for (int k = 0; k < K; ++k) {
#pragma _NEC ivdep
			for (int n = 0; n < N; ++n) {
				C[n + strideC * m] += alpha * A[k + K * n] * B[k + K * m];
			}
		}
	}
}

/*
* DGEMV wrapper
* C = alpha*(A^t*V) + beta*C
* matrix A is KxM grids (A[iy + K * ix] where ix in [0,M), iy in [0,K)
* matrix V is K grids (A[iy] where iy in [0,K)
*/
inline
void blas_DGEMV_t(int M, int K, double* A, double* V, double* C, double alpha, double beta) {

	cblas_dgemv(CblasRowMajor, CblasNoTrans, M, K, alpha, A, K, V, 1, beta, C, 1);

}

inline
void blas_ZGEMV_t(int M, int K, void* A, void* V, void* C, OneComplex alpha, OneComplex beta) {

    cblas_zgemv(CblasRowMajor, CblasNoTrans, M, K, &alpha, A, K, V, 1, &beta, C, 1);

}
