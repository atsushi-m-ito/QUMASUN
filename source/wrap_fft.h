#pragma once
//#define _USE_MATH_DEFINES
#include <complex>
#include <cmath>
#ifdef _NEC
#include <aslfftw3.h>
#else
#include <fftw3.h>
#endif
#include "wave_function.h"

template<class COMPLEX_T>
inline
void IFFT_3D(COMPLEX_T* dest_r, COMPLEX_T* src_k,
	const int Nx ,const int Ny , const int Nz )
{

	/*
	Forward FFT
	out[k] = \sum_x f[x] exp(-i 2\pi k x / N)
	*/

	const int dimension = 3;

	fftw_complex* in = (fftw_complex*)(src_k);
	fftw_complex* out = (fftw_complex*)(dest_r);

	fftw_plan p = fftw_plan_dft_3d(Nz, Ny, Nx, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_execute(p);
	fftw_destroy_plan(p);

	//inverse FFT in FFTW3 does not include 1/N coefficient//
	const size_t N3D = Nx * Ny * Nz;
	const double invN = 1.0 / (double)(N3D);
	for (size_t i = 0; i < N3D; ++i) {
		out[i][0] *= invN;
		out[i][1] *= invN;
	}
}

template<class COMPLEX_T>
inline
void IFFT_3D_c2r(double* dest_r, COMPLEX_T* src_k,
	const int Nx, const int Ny, const int Nz) {

	/*
	Forward FFT
	out[k] = \sum_x f[x] exp(-i 2\pi k x / N)
	*/

	const int dimension = 3;

	fftw_complex* in = (fftw_complex*)(src_k);
	double* out = (dest_r);

	/* from FFTW3 manual
	An r2c transform produces the same output as a FFTW_FORWARD complex DFT of the same input,
	and a c2r transform is correspondingly equivalent to FFTW_BACKWARD. For more information
	*/
	fftw_plan p = fftw_plan_dft_c2r_3d(Nz, Ny, Nx, in, out, FFTW_ESTIMATE);
	fftw_execute(p);
	fftw_destroy_plan(p);

	//inverse FFT in FFTW3 does not include 1/N coefficient//
	const size_t N3D = Nx * Ny * Nz;
	const double invN = 1.0 / (double)(N3D);
	for (size_t i = 0; i < N3D; ++i) {
		out[i] *= invN;
	}
}



template<class COMPLEX_T>
inline
void FFT_3D(COMPLEX_T* dest_k, COMPLEX_T* src_r,
	const int Nx, const int Ny, const int Nz) {

	/*
	Forward FFT
	out[k] = \sum_x f[x] exp(-i 2\pi k x / N)
	*/

	const int dimension = 3;

	fftw_complex* in = (fftw_complex*)(src_r);
	fftw_complex* out = (fftw_complex*)(dest_k);
	
	fftw_plan p = fftw_plan_dft_3d(Nz, Ny, Nx, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

	fftw_execute(p);

	fftw_destroy_plan(p);

}

template<class COMPLEX_T>
inline
void FFT_3D_r2c(const COMPLEX_T* dest_k, double* src_r,
	const int Nx, const int Ny, const int Nz) {

	/*
	Forward FFT
	out[k] = \sum_x f[x] exp(-i 2\pi k x / N)
	*/

	const int dimension = 3;

	double* in = src_r;
	fftw_complex* out = (fftw_complex*)dest_k;

	/* from FFTW3 manual
	An r2c transform produces the same output as a FFTW_FORWARD complex DFT of the same input, 
	and a c2r transform is correspondingly equivalent to FFTW_BACKWARD. For more information

	For out-of-place transforms, this is the end of the story: 
	the real data is stored as a row-major array of size n0 * n1 * n2 * ... * n_{d-1}
	and the complex data is stored as a row-major array of size n0 * n1 * n2 * ... * (n_{d-1}/2 + 1).
	For in-place transforms, however, extra padding of the real-data array is necessary ...
	*/
	fftw_plan p = fftw_plan_dft_r2c_3d(Nz, Ny, Nx, in, out, FFTW_ESTIMATE);

	fftw_execute(p);

	fftw_destroy_plan(p);

}
