#pragma once
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <cstdint>
#include "wrap_fft.h"
#include "fftw_executor.h"

/*
* necessary work size is (size_x * size_y * size_z * 4).
*/
inline
void SetPotentialByPoissonFFT(double* Vhart, const double* rho, int size_x, int size_y, int size_z,
	double dx, double dy, double dz, double* work) {

	dcomplex* rhok=(dcomplex*)work;
	const int64_t size_3d = (int64_t)size_x * (int64_t)size_y * (int64_t)size_z;

	//#undef USE_R2C_C2R
#ifndef USE_R2C_C2R

	dcomplex* rho_comp = (dcomplex*)(work + size_3d * 2);
	for (size_t i = 0; i < size_3d; ++i) {
		rho_comp[i] = { rho[i],0.0 };
	}

	FFT_3D(rhok, rho_comp, size_x, size_y, size_z);

	const int kz_end = size_z;
	const int kx_end = size_x;

#else
	/*
	problem: this path cannot obtain correct Vhart
	*/
	RspaceFunc<double>& rho_d = m_work_rd;


	{
		vecmath::Copy<double>(rho_d, rho, size_3d);
		//rho_d[m_size_x / 2 + m_size_x * (m_size_y / 2 + m_size_y * (m_size_z / 2))] -= 1.0 / (m_dx * m_dy * m_dz);
	}

	double sumrho = 0.0;
	for (size_t i = 0; i < size_3d; ++i) {
		sumrho += rho_d[i];
	}
	sumrho *= dx * dy * dz;
	printf("rho in poisson = %f\n", sumrho);

	/*fftw_plan_dft_r2c_3d*/
	FFT_3D(rhok, rho_d, m_size_x, m_size_y, m_size_z);

	const int kz_end = m_size_z;
	const int kx_end = m_size_x / 2 + 1;
#endif

	//NOTE: 
	//  V = c rho, where c = 4pi / (kx^2 + ky^2 + kz^2) //
	//    = 1.0/ ((pi*kx^2/box_x^2) + (pi*ky^2/box_y^2) + (pi*kz^2/box_z^2)) 

	const double coef1_x = (2.0 * M_PI / (dx * (double)size_x));
	const double coef1_y = (2.0 * M_PI / (dy * (double)size_y));
	const double coef1_z = (2.0 * M_PI / (dz * (double)size_z));

	auto SQ = [](double x) { return x * x; };

	//#define NO_FOLD

	for (int kz = 0; kz < kz_end; ++kz) {
#ifdef NO_FOLD
		const double kz2 = coef_z * SQ((double)(kz));
#else

		const double kz2 = SQ(coef1_z * (double)(kz < size_z / 2 ? kz : kz - size_z));
#endif
		for (int ky = 0; ky < size_y; ++ky) {
#ifdef NO_FOLD
			const double ky2 = coef_z * SQ((double)(ky));
#else
			const double ky2 = SQ(coef1_y * (double)(ky < size_y / 2 ? ky : ky - size_y));
#endif


			for (int kx = 0; kx < kx_end; ++kx) {
#ifdef NO_FOLD
				const double kx2 = coef_z * SQ((double)(kx));
#else
				const double kx2 = SQ(coef1_x * (double)(kx < size_x / 2 ? kx : kx - size_x));
#endif
				const size_t i = kx + kx_end * (ky + (size_y * kz));

				if (kx == 0 && ky == 0 && kz == 0) {
					rhok[i] = 0.0;
				} else {
					rhok[i] *= 4.0 * M_PI / (kx2 + ky2 + kz2);
				}
			}

		}
	}

#ifndef USE_R2C_C2R 

    dcomplex* Vh =(dcomplex*)(work + size_3d * 2);
	/*
	for (size_t i = 0; i < m_size_3d; ++i) {
		Vh[i] = { 0.0,0.0 };
	}
	*/
	IFFT_3D(Vh, rhok, size_x, size_y, size_z);

	for (size_t i = 0; i < size_3d; ++i) {
		Vhart[i] = Vh[i].real();
	}
#else

	IFFT_3D(Vhart, rhok, size_x, size_y, size_z);

#endif


}


/*
* necessary work size is (size_x * size_y * size_z * 6).
*/
inline
void SetPotentialByPoissonFFT_keepRhok(double* Vhart, const double* rho, int size_x, int size_y, int size_z,
    double dx, double dy, double dz, double* work) {

    const int64_t size_3d = (int64_t)size_x * (int64_t)size_y * (int64_t)size_z;
    dcomplex* rhok = (dcomplex*)work;
    dcomplex* rhok_kk = (dcomplex*)(work + size_3d * 2);

    //#undef USE_R2C_C2R
#ifndef USE_R2C_C2R

    dcomplex* rho_comp = (dcomplex*)(work + size_3d * 4);
    for (size_t i = 0; i < size_3d; ++i) {
        rho_comp[i] = { rho[i],0.0 };
    }

    FFT_3D(rhok, rho_comp, size_x, size_y, size_z);

    const int kz_end = size_z;
    const int kx_end = size_x;

#else
    /*
    problem: this path cannot obtain correct Vhart
    */
    RspaceFunc<double>& rho_d = m_work_rd;


    {
        vecmath::Copy<double>(rho_d, rho, size_3d);
        //rho_d[m_size_x / 2 + m_size_x * (m_size_y / 2 + m_size_y * (m_size_z / 2))] -= 1.0 / (m_dx * m_dy * m_dz);
    }

    double sumrho = 0.0;
    for (size_t i = 0; i < size_3d; ++i) {
        sumrho += rho_d[i];
    }
    sumrho *= dx * dy * dz;
    printf("rho in poisson = %f\n", sumrho);

    /*fftw_plan_dft_r2c_3d*/
    FFT_3D(rhok, rho_d, m_size_x, m_size_y, m_size_z);

    const int kz_end = m_size_z;
    const int kx_end = m_size_x / 2 + 1;
#endif

    //NOTE: 
    //  V = c rho, where c = 4pi / (kx^2 + ky^2 + kz^2) //
    //    = 1.0/ ((pi*kx^2/box_x^2) + (pi*ky^2/box_y^2) + (pi*kz^2/box_z^2)) 

    const double coef1_x = (2.0 * M_PI / (dx * (double)size_x));
    const double coef1_y = (2.0 * M_PI / (dy * (double)size_y));
    const double coef1_z = (2.0 * M_PI / (dz * (double)size_z));

    auto SQ = [](double x) { return x * x; };

    //#define NO_FOLD

    for (int kz = 0; kz < kz_end; ++kz) {
#ifdef NO_FOLD
        const double kz2 = coef_z * SQ((double)(kz));
#else

        const double kz2 = SQ(coef1_z * (double)(kz < size_z / 2 ? kz : kz - size_z));
#endif
        for (int ky = 0; ky < size_y; ++ky) {
#ifdef NO_FOLD
            const double ky2 = coef_z * SQ((double)(ky));
#else
            const double ky2 = SQ(coef1_y * (double)(ky < size_y / 2 ? ky : ky - size_y));
#endif


            for (int kx = 0; kx < kx_end; ++kx) {
#ifdef NO_FOLD
                const double kx2 = coef_z * SQ((double)(kx));
#else
                const double kx2 = SQ(coef1_x * (double)(kx < size_x / 2 ? kx : kx - size_x));
#endif
                const size_t i = kx + kx_end * (ky + (size_y * kz));

                if (kx == 0 && ky == 0 && kz == 0) {
                    rhok_kk[i] = 0.0;
                } else {
                    rhok_kk[i] = rhok[i] * 4.0 * M_PI / (kx2 + ky2 + kz2);
                }
            }

        }
    }

#ifndef USE_R2C_C2R 

    dcomplex* Vh = (dcomplex*)(work + size_3d * 4);
    /*
    for (size_t i = 0; i < m_size_3d; ++i) {
        Vh[i] = { 0.0,0.0 };
    }
    */
    IFFT_3D(Vh, rhok_kk, size_x, size_y, size_z);

    for (size_t i = 0; i < size_3d; ++i) {
        Vhart[i] = Vh[i].real();
    }
#else

    IFFT_3D(Vhart, rhok_kk, size_x, size_y, size_z);

#endif


}


inline 
void SolvePoissonInKspace(fftw_complex* rhok, int size_x, int size_y, int size_z, double dx, double dy, double dz) {
    const int kz_end = size_z;
    const int kx_end = size_x;

    //NOTE: 
    //  V = c rho, where c = 4pi / (kx^2 + ky^2 + kz^2) //
    //    = 1.0/ ((pi*kx^2/box_x^2) + (pi*ky^2/box_y^2) + (pi*kz^2/box_z^2)) 

    const double coef1_x = (2.0 * M_PI / (dx * (double)size_x));
    const double coef1_y = (2.0 * M_PI / (dy * (double)size_y));
    const double coef1_z = (2.0 * M_PI / (dz * (double)size_z));

    auto SQ = [](double x) { return x * x; };

    //#define NO_FOLD

    for (int kz = 0; kz < kz_end; ++kz) {
#ifdef NO_FOLD
        const double kz2 = coef_z * SQ((double)(kz));
#else

        const double kz2 = SQ(coef1_z * (double)(kz < size_z / 2 ? kz : kz - size_z));
#endif
        for (int ky = 0; ky < size_y; ++ky) {
#ifdef NO_FOLD
            const double ky2 = coef_z * SQ((double)(ky));
#else
            const double ky2 = SQ(coef1_y * (double)(ky < size_y / 2 ? ky : ky - size_y));
#endif


            for (int kx = 0; kx < kx_end; ++kx) {
#ifdef NO_FOLD
                const double kx2 = coef_z * SQ((double)(kx));
#else
                const double kx2 = SQ(coef1_x * (double)(kx < size_x / 2 ? kx : kx - size_x));
#endif
                const size_t i = kx + kx_end * (ky + (size_y * kz));

                if (kx == 0 && ky == 0 && kz == 0) {
                    rhok[i][0] = 0.0;
                    rhok[i][1] = 0.0;
                } else {
                    rhok[i][0] *= 4.0 * M_PI / (kx2 + ky2 + kz2);
                    rhok[i][1] *= 4.0 * M_PI / (kx2 + ky2 + kz2);
                }
            }

        }
    }
}


/*
* necessary work size is (size_x * size_y * size_z * 4).
*/
inline
void SetPotentialByPoissonFFT_v2(FFTW_Executor* fftw_e, double* Vhart, const double* rho, int size_x, int size_y, int size_z,
    double dx, double dy, double dz) {


    const int64_t size_3d = (int64_t)size_x * (int64_t)size_y * (int64_t)size_z;

    auto* rho_comp = fftw_e->GetBuffer();
    for (size_t i = 0; i < size_3d; ++i) {
        rho_comp[i][0] = rho[i];
        rho_comp[i][1] = 0.0;
    }

    auto rhok = fftw_e->ForwardExecute(nullptr, nullptr);

    SolvePoissonInKspace(rhok, size_x, size_y, size_z, dx, dy, dz);

    //IFFT_3D(Vh, rhok, size_x, size_y, size_z);
    auto* Vh = fftw_e->BackwardExecute(nullptr, nullptr);
    const double invN = 1.0 / (double)(size_3d);

    for (size_t i = 0; i < size_3d; ++i) {
        Vhart[i] = Vh[i][0] * invN;
    }


}
