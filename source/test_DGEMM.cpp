/*******************************

QUantum MAterial Simulation UNraveler (QUMASUN)

QUMASUN is numerical simulation code for Density Functional Theory(DFT) and Time-dependent DFT based on the real space grid.

********************************/

#include <cstdlib>
#include <cstdio>
#include <cstdint>
#include <cmath>
#include "wrap_chrono.h"
#include "print_matrix.h"
#include "w_dgemm.h"
#ifdef _NEC
#include <asl.h>
#endif



void BenchDGEMM(
	const int num_grids,// 32 * 32 * 32;
	const int num_rows_A,//100;
    const int num_cols_B,//100;
	const int STEPS,
	const bool is_output){
	double* A = new double[num_grids * num_rows_A];
	double* B = new double[num_grids * num_cols_B];
	double* C = new double[num_rows_A * num_cols_B];
	const double dV = 1.0 / (double)(num_grids);
	const double coef = dV / (double)STEPS;

	printf("Matrix A (in) : %d * %d\n", num_rows_A, num_grids);
	printf("Matrix B (in) : %d * %d\n", num_cols_B, num_grids);
	printf("Matrix C (out): %d * %d\n", num_rows_A, num_cols_B);
	printf("STEPS         : %d\n", STEPS);

	for (int n = 0; n < num_rows_A; ++n) {
		for (int ix = 0; ix < num_grids; ++ix) {
			A[ix + num_grids * n] = (double)(1+n);
		}
	}
    for (int n = 0; n < num_cols_B; ++n) {
        for (int ix = 0; ix < num_grids; ++ix) {
            B[ix + num_grids * n] = 1.0 / (double)(1 + n);
        }
    }

	{
		auto start_tm = std::chrono::high_resolution_clock::now();
		for (int istep = 0; istep < STEPS; ++istep) {
			for (int n = 0; n < num_rows_A* num_cols_B; ++n) {
				C[n] = 0.0;
			}
			blas_DGEMM_t(num_rows_A, num_cols_B, num_grids, A, B, C, num_rows_A, coef, 1.0);
		}
		auto end_tm = std::chrono::high_resolution_clock::now();
		if (is_output)OutputMatrix(C, num_rows_A, num_cols_B, "test_matrix_b.txt");
		printf("M * M dgemm: %f [s/step], %f [GFLOPS]\n", DoubleSec(end_tm - start_tm) / (double)STEPS,
			(double)STEPS * (double)num_rows_A  * (double)num_cols_B * (double)num_grids * 2.0e-9 / DoubleSec(end_tm - start_tm));
		fflush(stdout);
	}

#ifdef USE_DGEMMT
	{
		auto start_tm = std::chrono::high_resolution_clock::now();
		for (int istep = 0; istep < STEPS; ++istep) {
            for (int n = 0; n < num_rows_A * num_cols_B; ++n) {
				C[n] = 0.0;
			}
            blas_DGEMM_t(num_rows_A, num_cols_B, num_grids, A, B, C, num_rows_A, coef, 1.0);
		}
		auto end_tm = std::chrono::high_resolution_clock::now();
		if (is_output)OutputMatrix(C, num_rows_A, num_cols_B, "test_matrix_b.txt");
		printf("M * M dgemmt: %f [s/step], %f [GFLOPS]\n", DoubleSec(end_tm - start_tm) / (double)STEPS,
			(double)STEPS * ((double)num_rows_A * (double)(num_cols_B -1)/2.0) * (double)num_grids * 2.0e-9 / DoubleSec(end_tm - start_tm));
		fflush(stdout);
	}
#endif

	{
		auto start_tm = std::chrono::high_resolution_clock::now();
		for (int istep = 0; istep < STEPS; ++istep) {
			for (int n = 0; n < num_rows_A * num_cols_B; ++n) {
				C[n] = 0.0;
			}
			blas_DGEMM_t2(num_rows_A, num_cols_B, num_grids, A, B, C, num_cols_B, coef, 1.0);
		}
		auto end_tm = std::chrono::high_resolution_clock::now();
		if (is_output)OutputMatrix(C, num_rows_A, num_cols_B, "test_matrix_b2.txt");
		printf("M * M dgemm: %f [s/step], %f [GFLOPS]\n", DoubleSec(end_tm - start_tm) / (double)STEPS,
			(double)STEPS * (double)num_rows_A * (double)num_cols_B * (double)num_grids * 2.0e-9 / DoubleSec(end_tm - start_tm));
		fflush(stdout);
	}


	{
		double* A_t = new double[num_cols_B * num_grids];

		auto start_tm = std::chrono::high_resolution_clock::now();
		for (int istep = 0; istep < STEPS; ++istep) {
			for (int n = 0; n < num_rows_A * num_cols_B; ++n) {
				C[n] = 0.0;
			}
			/*
			int ierr = ASL_dam1tp(B, num_grids, num_grids, num_states, A_t, num_states);
			if (ierr > 0) {
				printf("ierr=%d\n", ierr);
				break;
			}
			*/
			transpose(num_grids, num_cols_B, B, num_grids, A_t, num_cols_B);
			if (is_output)OutputMatrix(A, num_grids, num_rows_A, "test_matrix_src_A.txt");
			if (is_output)OutputMatrix(B, num_grids, num_cols_B, "test_matrix_src_B.txt");
			if (is_output)OutputMatrix(A_t, num_cols_B, num_grids, "test_matrix_src_A_t.txt");
			//blas_DGEMM_n(num_states, num_states, num_grids, A_t, B, C, num_states, coef, 1.0);
			blas_DGEMM_n(num_rows_A, num_cols_B, num_grids, A, A_t, C, num_rows_A, coef, 1.0);
			//blas_DGEMM_t(num_states, num_states, num_grids, A, B, C, num_states, coef, 1.0);
		}
		auto end_tm = std::chrono::high_resolution_clock::now();
		if (is_output)OutputMatrix(C, num_rows_A, num_cols_B, "test_matrix_bn.txt");
		printf("M * M dgemm with dam1tp: %f [s/step], %f [GFLOPS]\n", DoubleSec(end_tm - start_tm) / (double)STEPS,
			(double)STEPS * (double)num_rows_A * (double)num_cols_B * (double)num_grids * 2.0e-9 / DoubleSec(end_tm - start_tm));
		fflush(stdout);
		delete[] A_t;
	}


	{
		auto start_tm = std::chrono::high_resolution_clock::now();

		for (int istep = 0; istep < STEPS; ++istep) {
			for (int n = 0; n < num_rows_A * num_cols_B; ++n) {
				C[n] = 0.0;
			}
			for (int m = 0; m < num_cols_B; ++m) {
				for (int n = 0; n < num_rows_A; ++n) {
					for (int ix = 0; ix < num_grids; ++ix) {
						C[n + num_rows_A * m] += coef * A[ix + num_grids * n] * B[ix + num_grids * m];
					}
				}
			}
		}
		auto end_tm = std::chrono::high_resolution_clock::now();
		if(is_output)OutputMatrix(C, num_rows_A, num_cols_B, "test_matrix.txt");
		printf("M * M handmade: %f [s/step], %f [GFLOPS]\n", DoubleSec(end_tm - start_tm) / (double)STEPS,
			(double)STEPS * (double)num_rows_A* (double)num_cols_B* (double)num_grids * 2.0e-9 / DoubleSec(end_tm - start_tm) );
        fflush(stdout);
	}

	{
		auto start_tm = std::chrono::high_resolution_clock::now();
		for (int istep = 0; istep < STEPS; ++istep) {
			for (int n = 0; n < num_rows_A * num_cols_B; ++n) {
				C[n] = 0.0;
			}
			nec_DMM_add(num_rows_A, num_cols_B, num_grids, A, B, C, num_rows_A, coef);
		}
		auto end_tm = std::chrono::high_resolution_clock::now();
		if(is_output)OutputMatrix(C, num_rows_A, num_cols_B, "test_matrix_n.txt");
		printf("M * M nec_DMM_add: %f [s/step], %f [GFLOPS]\n", DoubleSec(end_tm - start_tm) / (double)STEPS,
			(double)STEPS * (double)num_rows_A* (double)num_cols_B* (double)num_grids * 2.0e-9 / DoubleSec(end_tm - start_tm));
        fflush(stdout);
	}


	
	{
		auto start_tm = std::chrono::high_resolution_clock::now();
		for (int istep = 0; istep < STEPS; ++istep) {
			for (int n = 0; n < num_rows_A * num_cols_B; ++n) {
				C[n] = 0.0;
			}
			for (int m = 0; m < num_cols_B; ++m) {
				blas_DGEMV_t(num_rows_A, num_grids, A, &B[num_grids*m], &C[num_rows_A * m], coef, 1.0);
			}
		}
		auto end_tm = std::chrono::high_resolution_clock::now();
		if(is_output)OutputMatrix(C, num_rows_A, num_cols_B, "test_matrix_v.txt");
		printf("(M * V dgemv) x col: %f [s/step], %f [GFLOPS]\n", DoubleSec(end_tm - start_tm) / (double)STEPS,
			(double)STEPS * (double)num_rows_A * (double)num_cols_B * (double)num_grids * 1.0e-9 / DoubleSec(end_tm - start_tm));
        fflush(stdout);
	}
	
	delete[] A;
	delete[] B;
	delete[] C;
	
}

//test to cut sub firld and to interporate from field
int main(int argc, char* argv[]) {
  //  BenchDGEMM(6800, 1280, 1280, 10, false);
    BenchDGEMM(2184, 8, 8, 100000, false);
    BenchDGEMM(11678, 8, 8, 100000, false);
    BenchDGEMM(2184, 8, 40, 100000, false);
    BenchDGEMM(11678, 8, 40, 100000, false);
    BenchDGEMM(100, 100, 1000, 10000, false);
    BenchDGEMM(1000, 100, 1000, 10000, false);
    BenchDGEMM(4000, 100, 100, 1000, false);
    BenchDGEMM(32 * 32 * 32, 8, 8, 100, false);
    BenchDGEMM(32 * 32 * 32, 16, 16, 100, false);
    BenchDGEMM(32 * 32 * 32, 32, 32, 100, false);
    BenchDGEMM(32 * 32 * 32, 48, 48, 100, false);
    BenchDGEMM(32 * 32 * 32, 64, 64, 100, false);
    BenchDGEMM(32 * 32 * 32, 128, 128, 10, false);
    BenchDGEMM(32 * 32 * 32, 300, 300, 10, false);
    BenchDGEMM(32 * 32 * 32, 1000, 1000, 1, false);
	BenchDGEMM(1000, 1000, 1000, 1, false);
	BenchDGEMM(4000, 4000, 1000, 1, false);
	return 0;
}



