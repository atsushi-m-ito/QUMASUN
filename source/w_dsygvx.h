#pragma once
#ifdef USE_MPI
#include <mpi.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include "wrap_lapack.h"
#ifdef LAPACK_NO_HEADER
extern "C" {
#ifdef _NEC
#define DSYGVX dsygvx_
#endif
void DSYGVX(const int* itype, const char* jobz, const char* range,
	const char* uplo, const int* n, double* a,
	const int* lda, double* b, const int* ldb,
	const double* vl, const double* vu, const int* il,
	const int* iu, const double* abstol, int* m, double* w,
	double* z, const int* ldz, double* work,
	const int* lwork, int* iwork, int* ifail,
	int* info);
}
#endif

inline
int lapack_DSYGVX(double *A, double *B, double *eigen_values, double*eigen_vectors, int N) {

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



	
	DSYGVX(&ITYPE, &JOBZ, &RANGE, &UPLO, &n, A, &LDA, B, &LDB,
		&VL, &VU, &IL, &IU, &ABSTOL, &M, 
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


inline
int lapack_DSYGVX_range(double* A, double* B, double* eigen_values, double* eigen_vectors, int N, int ilower, int iupper) {

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

	char  RANGE = 'I';
	char  UPLO = 'U';

	INTEGER n = N;
	INTEGER LDA = n;
	INTEGER LDB = n;
	double VL, VU; /* dummy */
	double ABSTOL = LAPACK_ABSTOL;
	INTEGER M;

	INTEGER LDZ = n;
	INTEGER LIWORK;
	INTEGER* IWORK;
	INTEGER* IFAIL = (INTEGER*)malloc(sizeof(INTEGER) * N);
	INTEGER INFO;

	//int i, j;

	INTEGER LWORK = 8 * N; // > 8*N //
	double* WORK = (double*)malloc(sizeof(double) * LWORK);


	//INTEGER LRWORK = 7 * N;
	//double* RWORK = (double*)malloc(sizeof(double)*LRWORK);
	//memset(RWORK, 0, sizeof(double)*LRWORK);		/* AITUNE */

	LIWORK = 5 * N;
	IWORK = (INTEGER*)malloc(sizeof(INTEGER) * LIWORK);
	//memset(IWORK, 0, sizeof(INTEGER)*LIWORK);		/* AITUNE */

	INTEGER IL = ilower;
	INTEGER IU = iupper;

	
	DSYGVX(&ITYPE, &JOBZ, &RANGE, &UPLO, &n, A, &LDA, B, &LDB,
		&VL, &VU, &IL, &IU, &ABSTOL, &M,
		eigen_values, eigen_vectors, &LDZ, WORK, &LWORK, IWORK, IFAIL, &INFO);



	if (INFO > 0) {
		/* printf("\n%s: error in dsyevd_, info=%d\n\n",name,INFO); */
	} else if (INFO < 0) {
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



#ifdef USE_MPI

inline
int lapack_DSYGVX_mpi_split(double* A, double* B, double* eigen_values, double* eigen_vectors, int N, const MPI_Comm& mpi_comm) {
	constexpr int SOLUTION_PER_PROC = 4;
	int num_procs;
	int proc_id;
	MPI_Comm_size(mpi_comm, &num_procs);
	MPI_Comm_rank(mpi_comm, &proc_id);
	

	int used_procs = num_procs;
	if (N < num_procs * SOLUTION_PER_PROC) {
		used_procs = std::max(N / SOLUTION_PER_PROC, 1);		
	}

	if (proc_id == 0) {
		printf("Use %d / %d procs for LAPACK Eigen solver\n", used_procs, num_procs);
	}

	if (used_procs == 1) {
		if (proc_id == 0) {
			
			return lapack_DSYGVX(A, B, eigen_values, eigen_vectors, N);
		} else {
			return 0;
		}
	} 


		

	int ilower = 1 + (proc_id * N) / used_procs;  //1 need for fortran index
	int iupper = ((proc_id + 1) * N) / used_procs;

	int num_local_ev = iupper - ilower + 1;
	int res=0;
	if(proc_id==0){

		//double* evec = new double[num_local_ev*N];

		res = lapack_DSYGVX_range(A, B, eigen_values, eigen_vectors, N, ilower, iupper);


		int* local_size_list = new int[num_procs];
		int* displs = new int[num_procs];
		const int ibegin1 = iupper;
		local_size_list[0] = 0;
		displs[0] = 0;
		for (int pid = 1; pid < used_procs; ++pid) {
			int ibegin = (pid * N) / used_procs;
			int iend = ((pid + 1) * N) / used_procs;

			local_size_list[pid] = (iend - ibegin) * N;
			displs[pid] = (ibegin-ibegin1) * N;
		}
		for (int pid = used_procs; pid < num_procs; ++pid) {
			local_size_list[pid] = 0;
			displs[pid] = (N - ibegin1) * N;
		}

		MPI_Gatherv(nullptr, 0, MPI_DOUBLE, &eigen_vectors[ibegin1 *N], local_size_list, displs, MPI_DOUBLE, 0, mpi_comm);
		delete[] local_size_list;
		delete[] displs;
	} else if(proc_id < used_procs){

		res = lapack_DSYGVX_range(A, B, eigen_values, eigen_vectors, N, ilower, iupper);

		MPI_Gatherv(eigen_vectors, N * num_local_ev, MPI_DOUBLE, nullptr, nullptr, nullptr, MPI_DOUBLE, 0, mpi_comm);
	} else {
		//LAPACKを解かないprocess//
		MPI_Gatherv(nullptr, 0, MPI_DOUBLE, nullptr, nullptr, nullptr, MPI_DOUBLE, 0, mpi_comm);
			
	}
	return res;
	
}
#endif
