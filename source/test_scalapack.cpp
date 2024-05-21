/*******************************

Test for ScaLAPACK in C++

********************************/
#ifdef USE_MPI
#include <mpi.h>
#endif
#include <cstdlib>
#include <cstdio>
#include <cstdint>
#include <cmath>
#include <random>
#include "wrap_chrono.h"
#include "print_matrix.h"
#define DEBUG_PRINT
#include "w_dsygvx.h"
#include "w_dsyevx.h"

#ifdef USE_MPI
#include "w_pdsyevd.h"
#include "w_pdsyevx.h"
#include "w_pdsygvx.h"


//test to cut sub firld and to interporate from field
int main(int argc, char* argv[]) {
	MPI_Init(&argc, &argv);
	MPI_Comm mpi_comm = MPI_COMM_WORLD;
	int proc_id;
	int num_procs;
	MPI_Comm_rank(mpi_comm, &proc_id);
	MPI_Comm_size(mpi_comm, &num_procs);
	
	const int N = 8;
#if 1
	double A[] = { 2533.947598 ,    -2917.065942,    -1869.931743,    -1340.469133,    5044.671746 ,    -1741.604346,    -994.427265 ,    -624.017194 ,
		-2917.065942,    9087.184014 ,    -682.979643 ,    -512.728426 ,    -5683.111391,    5485.292180 ,    -272.934708 ,    -323.883183 ,
		-1869.931743,    -682.979643 ,    9763.949926 ,    -329.645406 ,    -3684.991000,    -337.947310 ,    5098.370988 ,    -191.253179 ,
		-1340.469133,    -512.728426 ,    -329.645406 ,    9963.609335 ,    -2582.361038,    -372.216164 ,    -197.337159 ,    4916.421493 ,
		5044.671746 ,    -5683.111391,    -3684.991000,    -2582.361038,    11784.950849,    -5446.336224,    -3447.425605,    -2234.852550,
		-1741.604346,    5485.292180 ,    -337.947310 ,    -372.216164 ,    -5446.336224,    9046.223708 ,    -539.512903 ,    -725.937482 ,
		-994.427265 ,    -272.934708 ,    5098.370988 ,    -197.337159 ,    -3447.425605,    -539.512903 ,    9936.743885 ,    -415.993280 ,
		-624.017194 ,    -323.883183 ,    -191.253179 ,    4916.421493 ,    -2234.852550,    -725.937482 ,    -415.993280 ,    10349.368937
	};
	double B[] = { 106.720703 ,     0.000000   ,     0.000000   ,     -0.000000  ,     -0.000000  ,     -56.753828 ,     -39.141999,      -29.097548,
				0.000000   ,     106.720703 ,     -0.000000  ,     -0.000000  ,     -61.710919 ,     -0.000000  ,     -14.296344,      -11.129790,
				0.000000   ,     -0.000000  ,     106.720703 ,     0.000000   ,     -39.558655 ,     -13.287910 ,     -0.000000 ,      -7.155609 ,
				-0.000000  ,     -0.000000  ,     0.000000   ,     106.720703 ,     -28.357803 ,     -9.975538  ,     -6.900241 ,      -0.000000 ,
				-0.000000  ,     -61.710919 ,     -39.558655 ,     -28.357803 ,     106.720703 ,     -8.336323  ,     -1.376146 ,      1.414583  ,
				-56.753828 ,     -0.000000  ,     -13.287910 ,     -9.975538  ,     -8.336323  ,     106.720703 ,     18.373812 ,      12.136711 ,
				-39.141999 ,     -14.296344 ,     -0.000000  ,     -6.900241  ,     -1.376146  ,     18.373812  ,     106.720703,      9.700419  ,
				-29.097548 ,     -11.129790 ,     -7.155609  ,     -0.000000  ,     1.414583   ,     12.136711  ,     9.700419  ,      106.720703
	};

#else
	
	double* A = new double[N * N];
	double* B = new double[N * N];
#endif

	double* C = new double[N * N];
	double* triangle = new double[N * N];
	double eigen_value[N];
	//const double dV = 1.0 / (double)num_grids;
#if 0
	for (int k = 0; k < N ; ++k) {
		for (int i = 0; i < N; ++i) {
			A[i + N * k] = (double)(1 + i) / (double)(1 + (i - k) * (i - k));
			B[i + N * k] = (i == k) ? 1.0 : 0.0;
		}
	}
#endif

#if 0
	lapack_PDSYEVD(A, eigen_value, C, N, mpi_comm);

	if (proc_id == 0) {
		OutputMatrix(C, N, N, "test_eigenvector.txt");
		OutputMatrix(eigen_value, N, 1, "test_eigenvalue.txt");
	}

	for (int i = 0; i < N; ++i)eigen_value[i] = 0.0;
	for (int i = 0; i < N*N; ++i)C[i] = 0.0;

	lapack_PDSYEVX(A, eigen_value, C, N, mpi_comm);

	if (proc_id == 0) {
		OutputMatrix(C, N, N, "test_eigenvector_x.txt");
		OutputMatrix(eigen_value, N, 1, "test_eigenvalue_x.txt");
	}


	for (int i = 0; i < N; ++i)eigen_value[i] = 0.0;
	for (int i = 0; i < N * N; ++i)C[i] = 0.0;

	lapack_PDSYGVX(A, B, eigen_value, C, N, mpi_comm);

	if (proc_id == 0) {
		OutputMatrix(C, N, N, "test_eigenvector_g.txt");
		OutputMatrix(eigen_value, N, 1, "test_eigenvalue_g.txt");
	}
#endif
	{
		double* A2 = new double[N * N];
		double* B2 = new double[N * N];
		double* C2 = new double[N * N];
		double* eigen_value2 = new double[N];
		printf("proc(%d):1\n", proc_id); fflush(stdout);
		for (int i = 0; i < N * N; ++i) {
			A2[i] = A[i];			
		}

		for (int k = 0; k < N; ++k) {
			for (int i = 0; i < N; ++i) {
				B2[i + N * k] = (i == k) ? 1.0 : 0.0;
			}
		}

		int sub_procs = (num_procs / 2 > 0) ? num_procs/2 : 1;
		printf("proc(%d):sub_procs = %d\n", proc_id, sub_procs); fflush(stdout);
		int key = proc_id % sub_procs;
		int color = proc_id / sub_procs;
		printf("proc(%d): key = %d, color=%d\n", proc_id, key, color); fflush(stdout);
		MPI_Comm m_ddm_comm;
		MPI_Comm_split(mpi_comm, color, key, &m_ddm_comm);


		int group_id = color;

		BlacsGridInfo blacs_grid = BeginBLACS(m_ddm_comm, group_id);

		if (group_id == 0) {
			lapack_PDSYGVX_split(A, B, eigen_value, C, N, 1, N, m_ddm_comm, blacs_grid);

			if (proc_id == 0) {
				OutputMatrix(C, N, N, "test_eigenvector_g.txt");
				OutputMatrix(eigen_value, N, 1, "test_eigenvalue_g.txt");
			}
		} else {
			lapack_PDSYGVX_split(A2, B2, eigen_value2, C2, N, 1, N, m_ddm_comm, blacs_grid);

			if (key == 0) {
				OutputMatrix(C2, N, N, "test_eigenvector_g2.txt");
				OutputMatrix(eigen_value2, N, 1, "test_eigenvalue_g2.txt");
			}
		}

		EndBLACS(blacs_grid);

		delete[] eigen_value2;
		delete[] A2;
		delete[] B2;
		delete[] C2;
	}



#if 0
	delete[] A;
	delete[] B;
#endif
	delete[] C;

	MPI_Finalize();
	return 0;
}

#else

inline
int lapack_DSYGVX_range(double* A, double* B, double* eigen_values, double* eigen_vectors, double* triangle_factor, int N, int ilower, int iupper) {

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
	INTEGER IL, IU;
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

	IL = ilower;
	IU = iupper;



	//Because lapack input matrix and output eigen vector is same pointor, A should be copied to eiven_vectors buffer//
	double* b = ((triangle_factor) ? triangle_factor : new double[N * N]);

	for (int i = 0; i < n * n; i++) {
		b[i] = B[i];
	}

	DSYGVX(&ITYPE, &JOBZ, &RANGE, &UPLO, &n, A, &LDA, b, &LDB,
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

	if (triangle_factor == NULL) {
		delete[] b;
	}

	return INFO;
}


//test to cut sub firld and to interporate from field
int main(int argc, char* argv[]) {
	const int N = 1000;
	const int STEPS = 4;
	const int SPLIT_COUNT = 4;
	double* A = new double[N * N];
	double* B = new double[N * N];
	double* Aorg = new double[N * N];
	double* Borg = new double[N * N];
	double* C = new double[N * N];
	double* triangle = new double[N * N];
	double eigen_value[N];
	//const double dV = 1.0 / (double)num_grids;

	printf("matrix size : %d^2\n", N);

	for (int k = 0; k < N; ++k) {
		for (int i = 0; i < N; ++i) {
			Aorg[i + N * k] = (double)(1 + i) / (double)(1 + (i - k) * (i - k));
			Borg[i + N * k] = (i == k) ? 1.0 : 0./ (double)(1 + (i - k) * (i - k));
		}
	}

	auto Copy = [](auto* dest, auto* src, size_t N) {
		for (size_t i = 0; i < N; ++i) {
			dest[i] = src[i];
		}
	};

	{

		auto start_tm = std::chrono::high_resolution_clock::now();
		for (int istep = 0; istep < STEPS; ++istep) {
			Copy(A, Aorg, N*N);
			Copy(B, Borg, N * N);
			lapack_DSYEVX(A, eigen_value, C, N);
		}

		auto end_tm = std::chrono::high_resolution_clock::now();
		printf("DSYEVX (V-A): %f [s/step]\n", DoubleSec(end_tm - start_tm) / (double)STEPS);

		OutputMatrix(C, N, N, "test_eigenvector.txt");
		OutputMatrix(eigen_value, N, 1, "test_eigenvalue.txt");
	}

	{
		auto start_tm = std::chrono::high_resolution_clock::now();
		for (int istep = 0; istep < STEPS; ++istep) {
			Copy(A, Aorg, N * N);
			Copy(B, Borg, N * N);
			lapack_DSYGVX(A, B, eigen_value, C, triangle, N);
		}

		auto end_tm = std::chrono::high_resolution_clock::now();
		printf("DSYGVX (V-A): %f [s/step]\n", DoubleSec(end_tm - start_tm) / (double)STEPS);

		OutputMatrix(C, N, N, "test_eigenvector_g.txt");
		OutputMatrix(eigen_value, N, 1, "test_eigenvalue_g.txt");
	}

	
	for(int isplit = 0; isplit < SPLIT_COUNT;isplit++){

		int ilower = 1 + (isplit * N) / SPLIT_COUNT;  //1 need for fortran index
		int iupper = ((isplit+1) * N) / SPLIT_COUNT;
		auto start_tm = std::chrono::high_resolution_clock::now();
		for (int istep = 0; istep < STEPS; ++istep) {
			Copy(A, Aorg, N * N);
			Copy(B, Borg, N * N);
			lapack_DSYGVX_range(A, B, eigen_value, C, triangle, N, ilower, iupper);
		}

		auto end_tm = std::chrono::high_resolution_clock::now();
		printf("DSYGVX (V-I[%d,%d]split): %f [s/step]\n", ilower, iupper, DoubleSec(end_tm - start_tm) / (double)STEPS);

	}
	

	delete[] A;
	delete[] B;
	delete[] Aorg;
	delete[] Borg;
	delete[] C;
}



#endif

