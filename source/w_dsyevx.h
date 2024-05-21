#pragma once
#include <stdio.h>
#include <stdlib.h>
#include "wrap_lapack.h"

inline
int lapack_DSYEVX(double *A, double *eigen_values, double*eigen_vectors, int N) {

	/*
	DSYEVX

	input:  N;
	input:  A[n][n];  matrix A
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

	INTEGER LWORK = 8 * N; // > 8*N //
	double* WORK = (double*)malloc(sizeof(double)*LWORK);


	//INTEGER LRWORK = 7 * N;
	//double* RWORK = (double*)malloc(sizeof(double)*LRWORK);
	//memset(RWORK, 0, sizeof(double)*LRWORK);		/* AITUNE */

	LIWORK = 5 * N;
	IWORK = (INTEGER*)malloc(sizeof(INTEGER)*LIWORK);
	//memset(IWORK, 0, sizeof(INTEGER)*LIWORK);		/* AITUNE */

	IL = 1;
	IU = N;

	
	DSYEVX(&JOBZ, &RANGE, &UPLO, &n, A, &LDA, &VL, &VU, &IL, &IU, &ABSTOL, &M, 
		eigen_values, eigen_vectors, &LDZ, WORK, &LWORK, IWORK, IFAIL, &INFO);



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
	free(IWORK);; free(WORK);

	return INFO;
}
