#pragma once

#include <stdio.h>
#include <stdlib.h>
#include "wrap_lapack.h"

#ifdef LAPACK_NO_HEADER
extern "C" {
#ifdef _NEC
#define DSTEVD dstevd_
#endif

void DSTEVD(const char* jobz, const int* n, double* d, double* e,
	double* z, const int* ldz, double* work,
	const int* lwork, int* iwork, const int* liwork,
	int* info );
}
#endif

inline
int lapack_DSTEVD(double *diagonal, double *off_diagonal, double *eigen_values, double *eigen_vectors, int N) {


	char JOBZ = ((eigen_vectors) ? 'V' : 'N');

	INTEGER n = N;
	INTEGER LDZ = n;
	//double *Z = (double*)malloc(sizeof(double)*n*n);

	
	INTEGER LWORK = 1 + 4 * n + n * n;
	double *WORK = (double*)malloc(sizeof(double)*LWORK);

	INTEGER LIWORK = 3 + 5 * n;
	INTEGER* IWORK = (INTEGER*)malloc(sizeof(INTEGER)*LIWORK);

	//double VL, VU;
	//INTEGER IL, IU;
	//double ABSTOL = LAPACK_ABSTOL;
	//INTEGER M;

	INTEGER INFO;


	//A=(double*)malloc(sizeof(double)*n*n);

	




	double* D = ((eigen_values != nullptr) ? eigen_values : new double[N]);
	
	for (int i = 0; i<n; i++) {
		D[i] = diagonal[i];
	}
	
	double* Z = ((eigen_vectors!= nullptr) ? eigen_vectors : new double[1]);

	//D is overriten by the eigen values//
	//Z is overwriten by the eigen vector//
	DSTEVD(&JOBZ, &n, D, off_diagonal, Z, &LDZ, WORK, &LWORK, IWORK, &LIWORK, &INFO);
	//dstevd(&JOBZ, &n, D, off_diagonal, Z, &LDZ, WORK, &LWORK, IWORK, &LIWORK, &INFO);



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
	
	free(IWORK); free(WORK);

	if (eigen_vectors == nullptr) { free(Z); }
	if (eigen_values == nullptr) { free(D); }

	return INFO;
}
