#ifdef USE_MPI
#pragma once
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include "wrap_scalapack.h"
#include "w_dsygvx.h"



inline
int lapack_PDSYGVX_split(double* A, double* B, double* eigen_values, double* eigen_vectors, int N, int ilower, int iupper, const MPI_Comm& mpi_comm, const BlacsGridInfo& blacs_grid) {

	int ibtype = 1;	// 1 indicate "A v = e B v" type problems//
	char  JOBZ = 'V';//((eigen_vectors) ? 'V' : 'N');

	int il = ilower;
	int iu = iupper;
	char  RANGE = ((il==1)&&(iu==N)) ? 'A' : 'I'; //all eigen values and vectors//
	char  UPLO = 'U';

	int n = N;
	int ZERO = 0;


	//printf("ScaLapack: start\n");

	int proc_id, num_procs;
	MPI_Comm_size(mpi_comm, &num_procs);
	MPI_Comm_rank(mpi_comm, &proc_id);


	int np_rows = blacs_grid.np_rows;
	int np_cols = blacs_grid.np_cols;
	int icontxt = blacs_grid.icontxt;


	int nblk;  //width of block cyclic//
	//check for num_proc is much larger than matrix size//
	{

		auto bits_msb = [](unsigned int v)
			{
				v = v | (v >> 1);
				v = v | (v >> 2);
				v = v | (v >> 4);
				v = v | (v >> 8);
				v = v | (v >> 16);
				return v ^ (v >> 1);
			};

		const int max_col_row = std::max(np_cols, np_rows);
		nblk = bits_msb(N / max_col_row);

#ifdef DEBUG_PRINT
		if (proc_id == 0) {
			printf("Block size in ScaLAPACK [%d/%d]: nblk=%d, np_rows=%d, np_cols=%d\n", proc_id, num_procs, nblk, np_rows, np_cols);
		}
#endif
		if (nblk <= 0) {  //when matrix size if much smaller than num_procs, call LAPACK on root procs//

			if (proc_id == 0) {
				return lapack_DSYGVX_range(A, B, eigen_values, eigen_vectors, N, ilower, iupper);
			} else {
				return 0;
			}
		}
	}



	int myrow, mycol;
	blacs_gridinfo_(&icontxt, &np_rows, &np_cols, &myrow, &mycol);

#ifdef DEBUG_PRINT
		printf("blacs_gridinfo_ [%d/%d contxt=%d]: myrow=%d, mycol=%d\n", proc_id, num_procs, icontxt, myrow, mycol);
		fflush(stdout);
#endif


	int local_row_size = numroc_(&n, &nblk, &myrow, &ZERO, &np_rows);
	int local_col_size = numroc_(&n, &nblk, &mycol, &ZERO, &np_cols);

#ifdef DEBUG_PRINT
	printf("numroc_ [%d/%d contxt=%d]: local_row_size=%d, local_col_size=%d\n", proc_id, num_procs, icontxt, local_row_size, local_col_size);
	fflush(stdout);
#endif
	//#undef DEBUG_PRINT    //NG

	double* local_a = new double[local_row_size * local_col_size * 2];
	double* local_b = new double[local_row_size * local_col_size * 2];
	double* local_z = new double[local_row_size * local_col_size * 2];


	for (int j = 0; j < local_col_size; j++) {
		for (int i = 0; i < local_row_size; i++) {
			int ig = np_rows * nblk * ((i) / nblk) + (i) % nblk + ((np_rows + myrow) % np_rows) * nblk;
			int jg = np_cols * nblk * ((j) / nblk) + (j) % nblk + ((np_cols + mycol) % np_cols) * nblk;
			local_a[j * local_row_size + i] = A[ig + N * jg];
			local_b[j * local_row_size + i] = B[ig + N * jg];
		}
	}




	int descA[9] = { 0 };//size is 9//
	int max_a_row = local_row_size;//MAX(1, localA_row);//	
	int info;
	descinit_(descA, &n, &n, &nblk, &nblk, &ZERO, &ZERO, &icontxt, &max_a_row, &info);
	int descB[9] = { 0 };
	descinit_(descB, &n, &n, &nblk, &nblk, &ZERO, &ZERO, &icontxt, &max_a_row, &info);
	int descZ[9] = { 0 };
	descinit_(descZ, &n, &n, &nblk, &nblk, &ZERO, &ZERO, &icontxt, &max_a_row, &info);

#ifdef DEBUG_PRINT
	printf("descinit_ [%d/%d contxt=%d]: descA=%d, %d, %d, %d, %d, %d, %d, %d, %d\n", proc_id, num_procs, icontxt, descA[0], descA[1], descA[2], descA[3], descA[4], descA[5], descA[6], descA[7], descA[8]);
#endif

	int ia = 1;  //always 1 is enough//
	int ja = 1;
	int ib = 1;
	int jb = 1;
	int iz = 1;
	int jz = 1;
	double vl = -100.0;
	double vu = 1000.0;
	const double ABSTOL = 0.0e-16;
	const double orfac = 1.0;

	int lwork;
	{
		int nn = std::max(std::max(n, nblk), 2);
		int NP0 = numroc_(&n, &nblk, &ZERO, &ZERO, &np_rows);
		lwork = 5 * n + std::max(5 * nn, NP0 * NP0 + 2 * nblk * nblk) + (int)ceil((double)n / (double)(np_rows * np_cols)) * nn;
		//lwork += 1000;//lwork *= 10;
	}
	int liwork;
	{
		int NNP = std::max(std::max(n, np_rows * np_cols + 1), 4);
		liwork = 6 * NNP;
	}
#ifdef DEBUG_PRINT    //OK
	printf("worksize [%d/%d contxt=%d]: lwork=%d, liwork=%d\n", proc_id, num_procs, icontxt, lwork, liwork);
#endif


	int num_found = 0;
	int num_found_z = 0;
	int* ifail = new int[N];
	int* iclustr = new int[2 * np_rows * np_cols];
	double* gap = new double[np_rows * np_cols];

	info = 0;

	double* work = new double[lwork];
	int* iwork = new int[liwork];

	pdsygvx_(&ibtype, &JOBZ, &RANGE, &UPLO, &n, local_a, &ia, &ja, descA,
		local_b, &ib, &jb, descB,
		&vl, &vu, &il, &iu, &ABSTOL, &num_found, &num_found_z,
		eigen_values, &orfac, local_z, &iz, &jz, descZ, work, &lwork, iwork, &liwork, ifail, iclustr, gap, &info);

#ifdef DEBUG_PRINT
	printf("pdsygvx_ [%d/%d contxt=%d]: info=%d, %d, %d\n", proc_id, num_procs, icontxt, info, num_found, num_found_z);
#endif



#ifdef DEBUG_PRINT
	printf("Cblacs_gridexit [%d/%d contxt=%d]:\n", proc_id, num_procs, icontxt); fflush(stdout);
#endif
	delete[] work;
	delete[] iwork;
	delete[] ifail;
	delete[] iclustr;
	delete[] gap;


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
#ifdef DEBUG_PRINT
	printf("MPI_Gather before[%d/%d contxt=%d]: %d x %d\n", proc_id, num_procs, icontxt, local_row_size, local_col_size); fflush(stdout);
	MPI_Barrier(mpi_comm);
#endif
	MPI_Gather(&my_inod[0], 4, MPI_INT, &xy_inod[0], 4, MPI_INT, 0, mpi_comm);
#ifdef DEBUG_PRINT
	printf("MPI_Gather [%d/%d contxt=%d]:\n", proc_id, num_procs, icontxt); fflush(stdout);
#endif
	int total_grid_size = 0;
	double* z_buffers = nullptr;
	if (proc_id == 0) {

		for (int pid = 0; pid < num_procs; ++pid) {
			//int irow = xy_inod[pid*4 + 0];
			//int icol = xy_inod[pid * 4 + 1];
			int LOCr = xy_inod[pid * 4 + 2];
			int LOCc = xy_inod[pid * 4 + 3];
#ifdef DEBUG_PRINT
			printf("[%d/%d contxt=%d]:LOCr = %d, LOCc = %d\n", proc_id, num_procs, icontxt, LOCr, LOCc); fflush(stdout);
#endif
			local_size_list[pid] = LOCr * LOCc;
			displs[pid] = total_grid_size;
			total_grid_size += LOCr * LOCc;
		}
#ifdef DEBUG_PRINT
		printf("total_grid_size [%d/%d contxt=%d]: %d\n", proc_id, num_procs, icontxt, total_grid_size); fflush(stdout);
#endif

		z_buffers = new double[total_grid_size];
		for (int i = 0; i < total_grid_size; ++i) {
			z_buffers[i] = 0.0;
		}

		MPI_Gatherv(local_z, local_row_size * local_col_size, MPI_DOUBLE, z_buffers, local_size_list, displs, MPI_DOUBLE, 0, mpi_comm);
#ifdef DEBUG_PRINT
		printf("MPI_Gatherv [%d/%d contxt=%d]:\n", proc_id, num_procs, icontxt); fflush(stdout);
#endif

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
#ifdef DEBUG_PRINT
		printf("MPI_Gatherv before [%d/%d contxt=%d]: %d x %d\n", proc_id, num_procs, icontxt, local_row_size, local_col_size); fflush(stdout);
#endif
		MPI_Gatherv(local_z, local_row_size * local_col_size, MPI_DOUBLE, nullptr, local_size_list, displs, MPI_DOUBLE, 0, mpi_comm);
#ifdef DEBUG_PRINT
		printf("MPI_Gatherv [%d/%d contxt=%d]:\n", proc_id, num_procs, icontxt); fflush(stdout);
#endif
	}


	delete[] local_a;
	delete[] local_b;
	delete[] local_z;

	return info;

}
//#undef DEBUG_PRINT
#endif
