#ifdef USE_MPI
#pragma once
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include "wrap_scalapack.h"

inline
int lapack_PDSYEVD(double* A, double* eigen_values, double* eigen_vectors, int N, const MPI_Comm& mpi_comm) {

	

	//INTEGER ITYPE = 1;	// 1 indicate "A v = e B v" type problems//
	char  JOBZ = 'V';//((eigen_vectors) ? 'V' : 'N');

	//char  RANGE = 'A';
	char  UPLO = 'U';

	int n = N;
	int nblk = 2;  //width of block cyclic//
	

	//printf("ScaLapack: start\n");

	int proc_id, num_procs;
	MPI_Comm_size(mpi_comm, &num_procs);
	MPI_Comm_rank(mpi_comm, &proc_id);



	int np_rows = (int)(sqrt((float)num_procs));
	do {
		if ((num_procs % np_rows) == 0) break;
		np_rows--;
	} while (np_rows >= 2);

	int np_cols = num_procs / np_rows;

	int icontxt = Csys2blacs_handle(mpi_comm);


	//int icontxt;
	//const char order = 'R'; //R means row-major (default)//
	//blacs_gridinit_(&icontxt, "R", &np_rows, &np_cols);
	//note: blacs_gridinit_ is not run, then sl_init_ is used.
	sl_init_(&icontxt, &np_rows, &np_cols);

	printf("Csys2blacs_handle [%d/%d contxt=%d]: cols*rows = %d * %d\n", proc_id, num_procs, icontxt, np_cols, np_rows);

	int myrow, mycol;
	blacs_gridinfo_(&icontxt, &np_rows, &np_cols, &myrow, &mycol);

	printf("blacs_gridinfo_ [%d/%d contxt=%d]: myrow=%d, mycol=%d\n", proc_id, num_procs, icontxt, myrow, mycol);

	const int ZERO = 0;
	int local_row_size = numroc_(&n, &nblk, &myrow, &ZERO, &np_rows);
	int local_col_size = numroc_(&n, &nblk, &mycol, &ZERO, &np_cols);

	printf("numroc_ [%d/%d contxt=%d]: local_row_size=%d, local_col_size=%d\n", proc_id, num_procs, icontxt, local_row_size, local_col_size);



	double* local_a = new double[local_row_size * local_col_size];
	double* local_z = new double[local_row_size * local_col_size];


	for (int j = 0; j < local_col_size; j++) {
		for (int i = 0; i < local_row_size; i++) {
			int ig = np_rows * nblk * ((i) / nblk) + (i) % nblk + ((np_rows + myrow) % np_rows) * nblk;
			int jg = np_cols * nblk * ((j) / nblk) + (j) % nblk + ((np_cols + mycol) % np_cols) * nblk;
			local_a[j * local_row_size + i] = A[ig + N * jg];
		}
	}




	int descA[9] = { 0 };//size is 9//
	int max_a_row = local_row_size;//MAX(1, localA_row);//	
	int info;
	descinit_(descA, &n, &n, &nblk, &nblk, &ZERO, &ZERO, &icontxt, &max_a_row, &info);
	int descZ[9] = { 0 };
	descinit_(descZ, &n, &n, &nblk, &nblk, &ZERO, &ZERO, &icontxt, &max_a_row, &info);

	printf("descinit_ [%d/%d contxt=%d]: descA=%d, %d, %d, %d, %d, %d, %d, %d, %d\n", proc_id, num_procs, icontxt, descA[0], descA[1], descA[2], descA[3], descA[4], descA[5], descA[6], descA[7], descA[8]);


	int ia = 1;  //always 1 is enough//
	int ja = 1;
	int iz = 1;
	int jz = 1;

	int lwork;
	{
		int NP = numroc_(&n, &nblk, &ZERO, &ZERO, &np_rows);
		int NQ = numroc_(&n, &nblk, &ZERO, &ZERO, &np_cols);
		int TRILWMIN = 3 * N + std::max(nblk * (NP + 1), 3 * nblk);
		lwork = std::max(1 + 6 * n + 2 * NP * NQ, TRILWMIN) + 2 * n;			
	}

	int liwork;
	{
		int NNP = std::max(std::max(n, np_rows * np_cols + 1), 4);
		liwork = 6 * NNP;
	}

	double work_size;
	int iwork_size;
/*
	info = 0;
	// get working size
	pdsyevd_(&JOBZ, &UPLO, &n, local_a, &ia, &ja, descA,
		eigen_values, local_z, &iz, &jz, descZ, &work_size, &lwork, &iwork_size, &liwork, &info);

	printf("pdsyevd_ [%d/%d contxt=%d]: work_size=%d, iwork_size=%d, info=%d\n", proc_id, num_procs, icontxt, (int)work_size, iwork_size, info);

	lwork = (int)work_size;
	liwork = iwork_size;
*/
	double* work = new double[lwork];
	int* iwork = new int[liwork];
	pdsyevd_(&JOBZ, &UPLO, &n, local_a, &ia, &ja, descA,
		eigen_values, local_z, &iz, &jz, descZ, work, &lwork, iwork, &liwork, &info);

	printf("pdsyevd_ [%d/%d contxt=%d]: info=%d\n", proc_id, num_procs, icontxt, info);


	blacs_gridexit_(&icontxt);

	printf("Cblacs_gridexit [%d/%d contxt=%d]:\n", proc_id, num_procs, icontxt);

	delete[] work;
	delete[] iwork;

	/*
	if (proc_id == 0) {
		printf("Eigen Values\n");
		int k;
		for (k = 0; k < N; k++) {
			printf("%f\n", eigen_values[k]);
		}
		fflush(stdout);
	}
	*/

	int* xy_inod = nullptr;
	int* local_size_list = nullptr;
	int* displs = nullptr;
	if (proc_id == 0) {
		xy_inod = new int[4 * num_procs];
		local_size_list = new int[num_procs];
		displs = new int[num_procs];
	}


	int my_inod[4];
	my_inod[0] = myrow;
	my_inod[1] = mycol;
	my_inod[2] = local_row_size;
	my_inod[3] = local_col_size;

	MPI_Gather(my_inod, 4, MPI_INT, xy_inod, 4, MPI_INT, 0, mpi_comm);

	int total_grid_size = 0;
	double* z_buffers = nullptr;
	if (proc_id == 0) {

		for (int pid = 0; pid < num_procs; ++pid) {
			//int irow = xy_inod[pid*4 + 0];
			//int icol = xy_inod[pid * 4 + 1];
			int LOCr = xy_inod[pid * 4 + 2];
			int LOCc = xy_inod[pid * 4 + 3];
			local_size_list[pid] = LOCr * LOCc;
			displs[pid] = total_grid_size;
			total_grid_size += LOCr * LOCc;
		}

		z_buffers = new double[total_grid_size];
		MPI_Gatherv(local_z, local_row_size* local_col_size, MPI_DOUBLE, z_buffers, local_size_list, displs, MPI_DOUBLE, 0, mpi_comm);		


		for (int pid = 0; pid < num_procs; ++pid) {
			int irow = xy_inod[pid * 4 + 0];
			int icol = xy_inod[pid * 4 + 1];
			int LOCr = xy_inod[pid * 4 + 2];
			int LOCc = xy_inod[pid * 4 + 3];			
			int offset = displs[pid];

			for (int j = 0; j < LOCc; j++) {
				for (int i = 0; i < LOCr; i++) {
					int ig = np_rows * nblk * ((i) / nblk) + (i) % nblk + ((np_rows + irow) % np_rows) * nblk;
					int jg = np_cols * nblk * ((j) / nblk) + (j) % nblk + ((np_cols + icol) % np_cols) * nblk;
					eigen_vectors[ig + N * jg] = z_buffers[offset + j * LOCr + i];
				}
			}

		}
		delete[] xy_inod;
		delete[] local_size_list;
		delete[] displs;
		delete[] z_buffers;
	} else {
		MPI_Gatherv(local_z, local_row_size* local_col_size, MPI_DOUBLE, nullptr, local_size_list, displs, MPI_DOUBLE, 0, mpi_comm);
	}


	delete[] local_a;
	delete[] local_z;

	return info;

}
#endif
