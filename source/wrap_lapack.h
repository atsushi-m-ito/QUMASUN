#pragma once

#if defined(__INTEL_LLVM_COMPILER ) || defined(__INTEL_COMPILER )
#include <complex>
#include <mkl.h>
#include <mkl_lapack.h>
#elif defined(_NEC)
#include <complex>
#include <cblas.h>

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

/*
typedef int INTEGER;

void DSYEVD( const char* jobz, const char* uplo, const INTEGER* n, double* a,
             const INTEGER* lda, double* w, double* work, const INTEGER* lwork,
             INTEGER* iwork, const INTEGER* liwork, INTEGER* info );

void dsygvd_( const INTEGER* itype, const char* jobz, const char* uplo,
              const INTEGER* n, double* a, const INTEGER* lda, double* b,
              const INTEGER* ldb, double* w, double* work,
              const INTEGER* lwork, INTEGER* iwork, const INTEGER* liwork,
              INTEGER* info );
*/

