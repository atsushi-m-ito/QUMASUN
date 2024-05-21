#pragma once

#include <stdio.h>
#include <stdlib.h>
#include "wrap_lapack.h"
#ifdef LAPACK_NO_HEADER
#ifndef lapack_complex_double
#define lapack_complex_double  std::complex<double>
#endif
extern "C"{
#ifdef _NEC
#define ZHEGVX  zhegvx_
#endif
void ZHEGVX(const int* itype, const char* jobz, const char* range,
	const char* uplo, const int* n, lapack_complex_double* a,
	const int* lda, lapack_complex_double* b, const int* ldb,
	const double* vl, const double* vu, const int* il,
	const int* iu, const double* abstol, int* m, double* w,
	lapack_complex_double* z, const int* ldz, lapack_complex_double* work,
	const int* lwork, double* rwork, int* iwork,
	int* ifail, int* info);
}
#endif

template<class MYCOMPLEX>
inline
int lapack_ZHEGVX(MYCOMPLEX*A, MYCOMPLEX*B, double *eigen_values, MYCOMPLEX*eigen_vectors, int N) {

	/*
	ZHEGVX

	input:  N;
	input:  A[n][n];  matrix A
	input:  B[n][n];  matrix B: S matrix
	output: eigen_vectors[n][n];  eigevectors
	output: eigen_values[n];    eigenvalues
	*/





	INTEGER ITYPE = 1;	// 1 indicate "A v = e B v" type problems//
	char  JOBZ = ((eigen_vectors) ? 'V' : 'N');

	char  RANGE = 'A';
	char  UPLO = 'U';

	INTEGER n = N;
	INTEGER LDA = n;
	INTEGER LDB = n;
	double VL, VU; /* dummy */
	INTEGER IL, IU;
	double ABSTOL = LAPACK_ABSTOL;
	INTEGER M;

	INTEGER LDZ = n;
	INTEGER LIWORK;
	INTEGER *IWORK;
	INTEGER *IFAIL = (INTEGER*)malloc(sizeof(INTEGER)*N);
	INTEGER INFO;

	//int i, j;

	INTEGER LWORK = 4 * N; // > 2*N //
	MYCOMPLEX* WORK = (MYCOMPLEX*)malloc(sizeof(MYCOMPLEX)*LWORK);


	INTEGER LRWORK = 7 * N;
	double* RWORK = (double*)malloc(sizeof(double)*LRWORK);
	//memset(RWORK, 0, sizeof(double)*LRWORK);		/* AITUNE */

	LIWORK = 5 * N;
	IWORK = (INTEGER*)malloc(sizeof(INTEGER)*LIWORK);
	//memset(IWORK, 0, sizeof(INTEGER)*LIWORK);		/* AITUNE */

	IL = 1;
	IU = N;



	ZHEGVX(&ITYPE, &JOBZ, &RANGE, &UPLO, &n, (lapack_complex_double*)A, &LDA, (lapack_complex_double*)B, &LDB,
		&VL, &VU, &IL, &IU, &ABSTOL, &M, 
		eigen_values, (lapack_complex_double*)eigen_vectors, &LDZ, (lapack_complex_double*)WORK, &LWORK, RWORK, IWORK, IFAIL, &INFO);



	if (INFO>0) {
		/* printf("\n%s: error in dsyevd_, info=%d\n\n",name,INFO); */
	} else if (INFO<0) {
		printf(" info=%d\n", INFO);

	} else { /* (INFO==0) */
			/* store eigenvectors */
			/*
			for (i=0;i<EVmax;i++) {
			for (j=0;j<n;j++) {
			a[i*n + j]= A[i*n+j];
			}
			}
			*/
	}
	free(IFAIL);
	free(IWORK); free(RWORK); free(WORK);

	return INFO;
}

