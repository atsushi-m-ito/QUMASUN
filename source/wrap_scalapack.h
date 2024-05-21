#pragma once
#include <mpi.h>

#if defined(__INTEL_LLVM_COMPILER ) || defined(__INTEL_COMPILER )
#include <complex>
#include <mkl.h>
#include <mkl_scalapack.h>
#elif defined(_NEC)
#include <complex>
#include <cblas.h>
#define MKL_Complex16 std::complex<double>

#else
#define HAVE_LAPACK_CONFIG_H
#define LAPACK_COMPLEX_CPP
#include <complex>
//#define lapack_complex_float std::complex<float>
//#define lapack_complex_double std::complex<double>
#define LAPACK_GLOBAL_PATTERN_UC
#define __EMSCRIPTEN__
#include <cblas.h>
//#include <f77blas.h>
#include <lapack.h>
#define MKL_Complex16 std::complex<double>
#endif

using INTEGER = int;

#ifndef ___dcomplex_definition___
using dcomplex = std::complex<double>;
//typedef struct { double r,i; } dcomplex;
#define ___dcomplex_definition___ 
#endif

#ifndef LAPACK_ABSTOL
#define LAPACK_ABSTOL (1.0E-13)
#endif



extern "C" {
    int Csys2blacs_handle(MPI_Comm comm);

    void sl_init_(int* icontext, int* nprow, int* npcolumn);
    void blacs_get_(int* icontext, const int* what, int* val);

    //void blacs_gridinit_(int* ConTxt, const char* layout, const int* nprow, const int* npcol);

    void blacs_gridinfo_(int* icontext, int* nprow, int* npcolumn, int* myrow,
        int* mycolumn);

    void blacs_gridmap_(int* icontext, int* usermap, int* ldumap, int* nprow, int* npcolumn);
    

    int numroc_(const int* n, const int* nb, const int* iproc,
        const int* isrcproc, const int* nprocs);

    
    void descinit_(int* desc, const int* m, const int* n,
        const int* mb, const int* nb, const int* irsrc,
        const int* icsrc, const int* ictxt, const int* lld,
        int* info);


    void pdsyevx_(const char* jobz, const char* range, const char* uplo, const int* n,
        const double* a, const int* ia, const int* ja, const int* desca,
        const double* vl, const double* vu, const int* il, const int* iu,
        const double* abstol, int* m, int* nz, double* w, const double* orfac,
        double* z, const int* iz, const int* jz, const int* descz,
        double* work, const int* lwork, int* iwork, const int* liwork,
        int* ifail, int* iclustr, double* gap, int* info);
    
    void pdsyevd_(const char* jobz, const char* uplo, const int* n,
        const double* a, const int* ia, const int* ja, const int* desca,
        double* w, double* z, const int* iz, const int* jz, const int* descz, 
        double* work, const int* lwork, int* iwork, const int* liwork, int* info);

    void pdsygvx_(const int* ibtype, const char* jobz, const char* range, const char* uplo, const int* n, double* a, const int* ia, const int* ja, const int* desca, double* b, const int* ib, const int* jb, const int* descb, const double* vl, const double* vu, const int* il, const int* iu, const double* abstol, int* m, int* nz, double* w, const double* orfac, double* z, const int* iz, const int* jz, const int* descz, double* work, const int* lwork, int* iwork, const int* liwork, int* ifail, int* iclustr, double* gap, int* info);
    

    void pzhegvx_(const int* ibtype, const char* jobz, const char* range, const char* uplo, const int* n, lapack_complex_double* a, const int* ia, const int* ja, const int* desca, lapack_complex_double* b, const int* ib, const int* jb, const int* descb, const double* vl, const double* vu, const int* il, const int* iu, const double* abstol, int* m, int* nz, double* w, const double* orfac, lapack_complex_double* z, const int* iz, const int* jz, const int* descz, lapack_complex_double* work, const int* lwork, double* rwork, const int* lrwork, int* iwork, const int* liwork, int* ifail, int* iclustr, double* gap, int* info);


    void blacs_exit_(int* cont);


    void blacs_gridexit_(int* icontext);

}



struct BlacsGridInfo {
	int icontxt;
	int np_rows;
	int np_cols;
};

inline
BlacsGridInfo BeginBLACS(MPI_Comm& mpi_comm, int group_id) {

	int proc_id, num_procs;
	MPI_Comm_size(mpi_comm, &num_procs);
	MPI_Comm_rank(mpi_comm, &proc_id);

	int ZERO = 0;

	int np_rows = (int)(sqrt((float)num_procs));
	do {
		if ((num_procs % np_rows) == 0) break;
		np_rows--;
	} while (np_rows >= 2);

#if 1
	int np_cols = num_procs / np_rows;
#else
	int np_cols = np_rows;
	np_rows = num_procs / np_cols;
#endif


	int orgcontxt;
	blacs_get_(&ZERO, &ZERO, &orgcontxt); //get default blacs context
#ifdef DEBUG_PRINT
	printf("blacs_get_ [%d/%d contxt=%d]:group_id=%d\n", proc_id, num_procs, orgcontxt, group_id); fflush(stdout);
#endif

	//make multi-grid (corresponding MPI_Comm_split)//
	//reference: https://www.ibm.com/docs/en/pessl/5.5?topic=blacs-gridmap-routine
	int icontxt;
	int TEN = 10;
	int* usermap = new int[num_procs];
	for (int i = 0; i < num_procs; ++i) {
		usermap[i] = i + num_procs * group_id;
	}
	//blacs_get_(&orgcontxt, &TEN, &icontxt); //redefinition mode//
	icontxt = orgcontxt;

#ifdef DEBUG_PRINT
	printf("blacs_get_(2) [%d/%d contxt=%d]:icontxt=%d\n", proc_id, num_procs, orgcontxt, icontxt); fflush(stdout);
#endif

	blacs_gridmap_(&icontxt, &usermap[0], &np_rows, &np_rows, &np_cols);
#ifdef DEBUG_PRINT
	printf("blacs_gridmap_ [%d/%d contxt=%d]:icontxt=%d\n", proc_id, num_procs, orgcontxt, icontxt); fflush(stdout);
#endif
	delete[] usermap;

	if (proc_id == 0) {
		printf("Use ScaLAPACK [%d/%d,group:%d]: np_rows=%d, np_cols=%d\n", proc_id, num_procs, group_id, np_rows, np_cols);
	}

	return BlacsGridInfo{ icontxt, np_rows, np_cols };
}

inline
void EndBLACS(const BlacsGridInfo& blacs_grid) {
	int icontxt = blacs_grid.icontxt;
	blacs_gridexit_(&icontxt);
}

