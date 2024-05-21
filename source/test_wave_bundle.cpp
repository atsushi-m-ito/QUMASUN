/**********************************************************************
* 
* Wave functionがNs本あるときに、それらを束ねたデータ構造として
* コンベンショナルなphi(x,y,z,n)のオーダーと
* バンドルを先に持ってくるphi(n,x,y,z)のオーダーで
* どちらが速いかの検証
* 
* Intel AVX512系列では、4次精度2階差分に関してphi(x,y,z,n)の方が速い
* 
* 
* 
* 
***********************************************************************/







#ifdef USE_MPI
#include <mpi.h>
#include "mpi_helper.h"
#endif
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cmath>
#include <chrono>
#include "GridRange.h"
#include "GridFor.h"
#include "GridDifference2nd.h"
#include "GridDifference4th.h"
#include "wrap_chrono.h"

// operatort A and wave function p_n
// for all n,
// res_n = <p_n|A|p_n>
//
// where, wave functions is bundled in ordere of x,y,z,n//
void PVP_xyzn(GridRange& grid, int64_t num_waves, double* result, const double* waves_xyzn, const double* V) {
	const int64_t Nxyz = grid.Size3D();
	for (int64_t n = 0; n < num_waves; ++n) {
		const double* wave_n = &waves_xyzn[n * Nxyz];
		double sum = 0.0;
#pragma ivdep
		for (int64_t i = 0; i < Nxyz; ++i) {
			sum += wave_n[i] * V[i] * wave_n[i];
		}
		result[n] = sum;
	}
}


// operatort A and wave function p_n
// for all n,
// res_n = <p_n|A|p_n>
//
// where, wave functions is bundled in ordere of n,x,y,z//
void PVP_nxyz(GridRange& grid, int64_t num_waves, double* result, const double* waves_nxyz, const double* V) {
	const int64_t Nxyz = grid.Size3D();
	for (int64_t n = 0; n < num_waves; ++n) {
		result[n] = .0;
	}
#pragma ivdep
	for (int64_t i = 0; i < Nxyz; ++i) {
		for (int64_t n = 0; n < num_waves; ++n) {
			const double p = waves_nxyz[n + i * num_waves];
			result[n] += p * V[i] * p;
		}
	}
}


// operatort A and wave function p_n, q_n
// for all n,
// res_n = <q_n|A|p_n>
//
// where, wave functions is bundled in ordere of x,y,z,n//
void QVP_xyzn(GridRange& grid, int64_t num_waves, double* result, const double* q_xyzn, const double* p_xyzn, const double* V) {
	const int64_t Nxyz = grid.Size3D();
	for (int64_t n = 0; n < num_waves; ++n) {
		const double* p_n = &p_xyzn[n * Nxyz];
		const double* q_n = &q_xyzn[n * Nxyz];
		double sum = 0.0;
#pragma ivdep
		for (int64_t i = 0; i < Nxyz; ++i) {
			sum += q_n[i] * V[i] * p_n[i];
		}
		result[n] = sum;
	}
}


// operatort A and wave function p_n
// for all n,
// res_n = <p_n|A|p_n>
//
// where, wave functions is bundled in ordere of n,x,y,z//
void QVP_nxyz(GridRange& grid, int64_t num_waves, double* result, const double* q_nxyz, const double* p_nxyz, const double* V) {
	const int64_t Nxyz = grid.Size3D();
	for (int64_t n = 0; n < num_waves; ++n) {
		result[n] = .0;
	}
#pragma ivdep
	for (int64_t i = 0; i < Nxyz; ++i) {
		for (int64_t n = 0; n < num_waves; ++n) {
			const double& p = p_nxyz[n + i * num_waves];
			const double& q = q_nxyz[n + i * num_waves];
			result[n] += q * V[i] * p;
		}
	}
}



// operatort A and wave function p_n
// for all n,
// |r_n> = A|p_n>
// where r_n is also bundled vectors.
// where, wave functions is bundled in ordere of x,y,z,n//
void VP_xyzn(GridRange& grid, int64_t num_waves, double* bundled_result, const double* waves_xyzn, const double* V) {
	const int64_t Nxyz = grid.Size3D();
	for (int64_t n = 0; n < num_waves; ++n) {
		const double* wave_n = &waves_xyzn[n * Nxyz];
		double* res_n = &bundled_result[n * Nxyz];
		double sum = 0.0;
#pragma ivdep
		for (int64_t i = 0; i < Nxyz; ++i) {
			res_n[i] += V[i] * wave_n[i];
		}
	}
}


// operatort A and wave function p_n
// for all n,
// res_n = <p_n|A|p_n>
//
// where, wave functions is bundled in ordere of n,x,y,z//
void VP_nxyz(GridRange& grid, int64_t num_waves, double* bundled_result, const double* waves_nxyz, const double* V) {
	const int64_t Nxyz = grid.Size3D();
#pragma ivdep
	for (int64_t i = 0; i < Nxyz; ++i) {
		for (int64_t n = 0; n < num_waves; ++n) {
			const double& p = waves_nxyz[n + i * num_waves];
			double& r = bundled_result[n + i * num_waves];
			r += V[i] * p;
		}
	}
}



void TestWaveBundle1() {
	constexpr int STEP = 10;
	constexpr int64_t Ns = 1024; //number of wave function (bundle)//
	constexpr int64_t Nx = 32;
	constexpr int64_t Ny = 32;
	constexpr int64_t Nz = 32;
	constexpr int64_t size_xyz = Nx * Ny * Nz;
	const double dx = 0.1;
	const double dy = 0.1;
	const double dz = 0.1;
	const double center = (Nx * dx) / 2.0;

	double* waves = new double[size_xyz * Ns];
	double* Hp = new double[size_xyz * Ns];
	double* V = new double[size_xyz];   //potential//
	double* result_xyzn = new double[Ns];
	double* result_nxyz = new double[Ns];
	GridRange grid{ 0,0,0,Nx,Ny,Nz };

#if 1
	printf("Data order xyzn======================\n");


	ForXYZ(grid, [&](size_t i, size_t ix, size_t iy, size_t iz)
		{
			const double x = ix * dx;
			const double y = iy * dy;
			const double z = iz * dz;
			V[i] = (x - center) * (x - center) + (y - center) * (y - center) + +(z - center) * (z - center);
		});

	//initialize waves ordered x,y,z,n//
	for (int64_t n = 0; n < Ns; ++n) {
		double* wave_n = &waves[n * size_xyz];
		double* Hp_n = &Hp[n * size_xyz];
		double norm = 0.0;
		ForXYZ(grid, [&](size_t i, size_t ix, size_t iy, size_t iz)
			{
				const double x = ix * dx;
				const double y = iy * dy;
				const double z = iz * dz;
				wave_n[i] = exp(-((x - center) * (x - center) + (y - center) * (y - center) + +(z - center) * (z - center))/2.0);
				norm += wave_n[i] * wave_n[i];
				Hp_n[i] = 0.0;
			});
		norm = 1.0 / sqrt(norm);
		For(grid, [&](size_t i)
			{
				wave_n[i] *= norm;
			});

	}
	

	{
		const auto begin_tm = std::chrono::high_resolution_clock::now();
		for (int istep = 0; istep < STEP; ++istep) {
			PVP_xyzn(grid, Ns, result_xyzn, waves, V);
		}
		const auto end_tm = std::chrono::high_resolution_clock::now();
		printf("PVP_xyzn: time = %f [s]\n", DoubleSec(end_tm - begin_tm));
	}


	{
		const auto begin_tm = std::chrono::high_resolution_clock::now();
		for (int istep = 0; istep < STEP; ++istep) {
			VP_xyzn(grid, Ns, Hp, waves, V);
		}
		const auto end_tm = std::chrono::high_resolution_clock::now();
		printf("VP_xyzn: time = %f [s]\n", DoubleSec(end_tm - begin_tm));
	}
	
	{
		double* b_wide = new double[(Nx + 2) * (Ny + 2) * (Nz + 2)];
		const auto begin_tm = std::chrono::high_resolution_clock::now();
		for (int istep = 0; istep < STEP; ++istep) {
			for (int64_t n = 0; n < Ns; ++n) {
				Laplasian2nd_overspace2(grid, Hp + size_xyz * n, waves + size_xyz * n, b_wide, -1.0 / (2.0 * (double)(Ns * 2)), dx, dy, dz);
			}
		}
		const auto end_tm = std::chrono::high_resolution_clock::now();
		printf("Laplasian2nd_overspace: time = %f [s]\n", DoubleSec(end_tm - begin_tm));
		delete[] b_wide;
	}
	{
		double* b_wide = new double[(Nx + 6) * (Ny + 6) * (Nz + 6)];
		const auto begin_tm = std::chrono::high_resolution_clock::now();
		for (int istep = 0; istep < STEP; ++istep) {
			for (int64_t n = 0; n < Ns; ++n) {
				Laplasian4th_overspace2(grid, Hp + size_xyz * n, waves + size_xyz * n, b_wide, -1.0 / (2.0 * (double)(Ns * 2)), dx, dy, dz);
			}
		}
		const auto end_tm = std::chrono::high_resolution_clock::now();
		printf("Laplasian4th_overspace: time = %f [s]\n", DoubleSec(end_tm - begin_tm));
		delete[] b_wide;
	}

	{
		const auto begin_tm = std::chrono::high_resolution_clock::now();
		for (int istep = 0; istep < STEP; ++istep) {
			QVP_xyzn(grid, Ns, result_xyzn, waves, Hp, V);
		}
		const auto end_tm = std::chrono::high_resolution_clock::now();
		printf("QVP_xyzn: time = %f [s]\n", DoubleSec(end_tm - begin_tm));
	}
#endif

#if 1
	printf("Data order nxyz======================\n");
	///////////////////////////////////////////////////////////////////////////////////////////////////////
	//initialize waves ordered n,x,y,z//
	{
		double norm = 0.0;
		ForXYZ(grid, [&](size_t i, size_t ix, size_t iy, size_t iz)
			{
				const double x = ix * dx;
				const double y = iy * dy;
				const double z = iz * dz;
				double val = (x - center) * (x - center) + (y - center) * (y - center) + +(z - center) * (z - center);
				for (int64_t n = 0; n < Ns; ++n) {
					waves[n + Ns * i] = exp(-val / 2.0);
					norm += waves[n + Ns * i] * waves[n + Ns * i];
					Hp[n + Ns * i] = 0.0;
				}
			});
		norm = 1.0 / (sqrt(norm / (double)Ns));
		For(grid, [&](size_t i)
			{
				for (int64_t n = 0; n < Ns; ++n) {
					waves[n + Ns * i] *= norm;
				}
			});
	}


	{
		const auto begin_tm = std::chrono::high_resolution_clock::now();
		for (int istep = 0; istep < STEP; ++istep) {
			PVP_nxyz(grid, Ns, result_nxyz, waves, V);
		}
		const auto end_tm = std::chrono::high_resolution_clock::now();
		printf("PVP_nxyz: time = %f [s]\n", DoubleSec(end_tm - begin_tm));
	}
	{
		const auto begin_tm = std::chrono::high_resolution_clock::now();
		for (int istep = 0; istep < STEP; ++istep) {
			VP_nxyz(grid, Ns, Hp, waves, V);
		}
		const auto end_tm = std::chrono::high_resolution_clock::now();
		printf("VP_nxyz: time = %f [s]\n", DoubleSec(end_tm - begin_tm));
	}

	{
		const auto begin_tm = std::chrono::high_resolution_clock::now();
		for (int istep = 0; istep < STEP; ++istep) {
			Laplasian2nd_bundle(grid, Ns, Hp, waves, -1.0 / (2.0 * (double)(Ns * 2)), dx, dy, dz);
		}
		const auto end_tm = std::chrono::high_resolution_clock::now();
		printf("Laplasian2nd_bundle: time = %f [s]\n", DoubleSec(end_tm - begin_tm));
	}
	{
		const auto begin_tm = std::chrono::high_resolution_clock::now();
		for (int istep = 0; istep < STEP; ++istep) {
			Laplasian4th_bundle(grid, Ns, Hp, waves, -1.0 / (2.0 * (double)(Ns * 2)), dx, dy, dz);
		}
		const auto end_tm = std::chrono::high_resolution_clock::now();
		printf("Laplasian4th_bundle: time = %f [s]\n", DoubleSec(end_tm - begin_tm));
	}

	{
		const auto begin_tm = std::chrono::high_resolution_clock::now();
		for (int istep = 0; istep < STEP; ++istep) {
			QVP_nxyz(grid, Ns, result_nxyz, waves, Hp, V);
		}
		const auto end_tm = std::chrono::high_resolution_clock::now();
		printf("QVP_nxyz: time = %f [s]\n", DoubleSec(end_tm - begin_tm));
	}
#endif
	

	{//comparison//
		int64_t error_count = 0;
		for (int64_t n = 0; n < Ns; ++n) {
			if (fabs(result_xyzn[n] - result_nxyz[n]) > 1.0e-5) {
				printf("ERROR: %zd, %f != %f\n", n, result_xyzn[n], result_nxyz[n]);
				error_count++;
			} else {
				//printf("OK: %zd, %f != %f\n", n, result_xyzn[n], result_nxyz[n]);
			}
		}
		if (error_count == 0) {
			printf("All clear:\n");
		} else {
			printf("ERROR COUNT: %zd\n", error_count);
		}
	}
	
	delete[] V;
	delete[] result_xyzn;
	delete[] result_nxyz;
	delete[] waves;
	delete[] Hp;
}


int main(int argc, char* argv[]) {
	TestWaveBundle1();
}
