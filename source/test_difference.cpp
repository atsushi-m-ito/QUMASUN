

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cmath>
#include <chrono>
#include "wrap_chrono.h"
#include "GridDifference2nd.h"
#include "GridDifference4th.h"

//3次元デカルト座標グリッドの計算のためのループを司る昨日//
//MPI並列時にはリダクションも提供する

inline
void Diff2_1_0(double* a, const double* b, const double* m, int grid_size) {
	constexpr int margin_size = 1;
#pragma ivdep
	for (int ix = margin_size; ix < grid_size - margin_size; ++ix) {
		const int ix_m1 = ix - 1;
		const int ix_p1 = ix + 1;

		a[ix] = (b[ix_m1] - 2.0 * b[ix] + b[ix_p1]) / 2.0;
	}
	a[0] = (m[margin_size*2 - 1] - 2.0 * b[0] + b[1]) / 2.0;
	a[grid_size - 1] = (b[grid_size - 2] - 2.0 * b[grid_size - 1] + m[0]) / 2.0;
}

inline
void Diff2_1_m(double* a, const double* b, const double* m, int grid_size) {
	constexpr int margin_size = 1;
#pragma ivdep
	for (int ix = 0; ix < grid_size; ++ix) {

		const double b_m1 = ((ix - 1 < 0) ? m[ix - 1 + margin_size] : b[ix - 1]);
		const double b_p1 = ((ix + 1 >= grid_size) ? m[ix + 1 - margin_size] : b[ix + 1]);

		a[ix] = (b_m1 - 2.0 * b[ix] + b_p1) / 2.0;
	}
}

inline
void Diff2_1_2(double* a, const double* b, const double* m, int grid_size) {
	constexpr int margin_size = 1;
#pragma ivdep
	for (int ix = 0; ix < grid_size; ++ix) {

		const double b_m1 = ((ix - 1 < 0) ? 0.0 : b[ix - 1]);
		const double b_p1 = ((ix + 1 >= grid_size) ? 0.0 : b[ix + 1]);

		a[ix] = (b_m1 - 2.0 * b[ix] + b_p1) / 2.0;
	}
	a[0] = (m[margin_size*2 - 1] - 2.0 * b[0] + b[1]) / 2.0;
	a[grid_size - 1] = (b[grid_size - 2] - 2.0 * b[grid_size - 1] + m[0]) / 2.0;
}

inline
void Diff2_1_p(double* a, const double* p, const double* m, int size_x, int size_y, int size_z) {
	constexpr int margin_size = 1;
	const double m_dx = 0.025;
	const double m_dy = 0.025;
	const double m_dz = 0.025;

#pragma ivdep
	for (int iz = 0; iz < size_z; ++iz) {
		const int idz_m = ((iz == 0) ? size_z - 1 : -1) * size_x * size_y;
		const int idz_p = ((iz == size_z - 1) ? 1 - size_z : 1) * size_x * size_y;
		
		for (int iy = 0; iy < size_y; ++iy) {
			const int idy_m = ((iy == 0) ? size_y - 1 : -1) * size_x;
			const int idy_p = ((iy == size_y - 1) ? 1 - size_y : 1) * size_x;
			for (int ix = margin_size; ix < size_x - margin_size; ++ix) {
				const int idx_m = ((ix == 0) ? size_x - 1 : -1);
				const int idx_p = ((ix == size_x - 1) ? 1 - size_x : 1);
				const size_t i = ix + size_x * (iy + (size_y * iz));

				double d2psidx2 = ((p[i + idx_p] + p[i + idx_m]) - 2.0* p[i]) / (m_dx * m_dx);
				d2psidx2 += ((p[i + idy_p] + p[i + idy_m]) - 2 * p[i]) / (m_dy * m_dy);
				d2psidx2 += ((p[i + idz_p] + p[i + idz_m]) - 2 * p[i]) / (m_dz * m_dz);
				a[i] += -d2psidx2 ;
			}
		}
	}
}

inline
void Diff2_2_p(double* a, const double* p, const double* m, int size_x, int size_y, int size_z) {
	constexpr int margin_size = 1;
	const double m_dx = 0.025;
	const double m_dy = 0.025;
	const double m_dz = 0.025;

#pragma ivdep
	for (int iz = 0; iz < size_z; ++iz) {
		const int idz_m2 = ((iz <= 1) ? size_z - 2 : -2) * size_x * size_y;
		const int idz_m = ((iz == 0) ? size_z - 1 : -1) * size_x * size_y;
		const int idz_p = ((iz == size_z - 1) ? 1 - size_z : 1) * size_x * size_y;
		const int idz_p2 = ((iz >= size_z - 2) ? 2 - size_z : 2) * size_x * size_y;

		for (int iy = 0; iy < size_y; ++iy) {
			const int idy_m2 = ((iy <= 1) ? size_y - 2 : -2) * size_x;
			const int idy_m = ((iy == 0) ? size_y - 1 : -1) * size_x;
			const int idy_p = ((iy == size_y - 1) ? 1 - size_y : 1) * size_x;
			const int idy_p2 = ((iy >= size_y - 2) ? 2 - size_y : 2) * size_x;

			for (int ix = margin_size; ix < size_x - margin_size; ++ix) {
				const int idx_m2 = ((ix <= 1) ? size_x - 2 : -2);
				const int idx_m = ((ix == 0) ? size_x - 1 : -1);
				const int idx_p = ((ix == size_x - 1) ? 1 - size_x : 1);
				const int idx_p2 = ((ix >= size_x - 2) ? 2 - size_x : 2);
				const size_t i = ix + size_x * (iy + (size_y * iz));

				double d2psidx2 = (-54.0 * (p[i + idx_p2] + p[i + idx_m2]) + 783.0 * (p[i + idx_p] + p[i + idx_m]) - 1460.0 * p[i]) / (m_dx * m_dx);
				d2psidx2 += (-54.0 * (p[i + idy_p2] + p[i + idy_m2]) + 783.0 * (p[i + idy_p] + p[i + idy_m]) - 1460.0 * p[i]) / (m_dy * m_dy);
				d2psidx2 += (-54.0 * (p[i + idz_p2] + p[i + idz_m2]) + 783.0 * (p[i + idz_p] + p[i + idz_m]) - 1460.0 * p[i]) / (m_dz * m_dz);
				a[i] += -d2psidx2 / (2.0 * 24.0 * 24.0);
			}
		}
	}
}

inline
void Diff2_3_p(double* a, const double* p, const double* m, int size_x, int size_y, int size_z) {
	constexpr int margin_size = 1;
	const double m_dx = 0.025;
	const double m_dy = 0.025;
	const double m_dz = 0.025;

#pragma ivdep
	for (int iz = 0; iz < size_z; ++iz) {
		const int idz_m3 = ((iz <= 2) ? size_z - 3 : -3) * size_x * size_y;
		const int idz_m2 = ((iz <= 1) ? size_z - 2 : -2) * size_x * size_y;
		const int idz_m = ((iz == 0) ? size_z - 1 : -1) * size_x * size_y;
		const int idz_p = ((iz == size_z - 1) ? 1 - size_z : 1) * size_x * size_y;
		const int idz_p2 = ((iz >= size_z - 2) ? 2 - size_z : 2) * size_x * size_y;
		const int idz_p3 = ((iz >= size_z - 3) ? 3 - size_z : 3) * size_x * size_y;

		for (int iy = 0; iy < size_y; ++iy) {
			const int idy_m3 = ((iy <= 2) ? size_y - 3 : -3) * size_x;
			const int idy_m2 = ((iy <= 1) ? size_y - 2 : -2) * size_x;
			const int idy_m = ((iy == 0) ? size_y - 1 : -1) * size_x;
			const int idy_p = ((iy == size_y - 1) ? 1 - size_y : 1) * size_x;
			const int idy_p2 = ((iy >= size_y - 2) ? 2 - size_y : 2) * size_x;
			const int idy_p3 = ((iy >= size_y - 3) ? 3 - size_y : 3) * size_x;
			for (int ix = 0; ix < size_x ; ++ix) {
				const int idx_m3 = ((ix <= 2) ? size_x - 3 : -3);
				const int idx_m2 = ((ix <= 1) ? size_x - 2 : -2);
				const int idx_m = ((ix == 0) ? size_x - 1 : -1);
				const int idx_p = ((ix == size_x - 1) ? 1 - size_x : 1);
				const int idx_p2 = ((ix >= size_x - 2) ? 2 - size_x : 2);
				const int idx_p3 = ((ix >= size_x - 3) ? 3 - size_x : 3);
				const size_t i = ix + size_x * (iy + (size_y * iz));

				double d2psidx2 = ((p[i + idx_p3] + p[i + idx_m3]) - 54.0 * (p[i + idx_p2] + p[i + idx_m2]) + 783.0 * (p[i + idx_p] + p[i + idx_m]) - 1460.0 * p[i]) / (m_dx * m_dx);
				d2psidx2 += ((p[i + idy_p3] + p[i + idy_m3]) - 54.0 * (p[i + idy_p2] + p[i + idy_m2]) + 783.0 * (p[i + idy_p] + p[i + idy_m]) - 1460.0 * p[i]) / (m_dy * m_dy);
				d2psidx2 += ((p[i + idz_p3] + p[i + idz_m3]) - 54.0 * (p[i + idz_p2] + p[i + idz_m2]) + 783.0 * (p[i + idz_p] + p[i + idz_m]) - 1460.0 * p[i]) / (m_dz * m_dz);
				a[i] += -d2psidx2 / (2.0 * 24.0 * 24.0);
			}
		}
	}
}

inline
void Diff2_3_0(double* a, const double* p, const double* m, int size_x, int size_y, int size_z) {
	constexpr int margin_size = 3;
	const double m_dx = 0.025;
	const double m_dy = 0.025;
	const double m_dz = 0.025;

#pragma ivdep
	for (int iz = 0; iz < size_z; ++iz) {
		const int idz_m3 = ((iz <= 2) ? size_z - 3 : -3) * size_x * size_y;
		const int idz_m2 = ((iz <= 1) ? size_z - 2 : -2) * size_x * size_y;
		const int idz_m = ((iz == 0) ? size_z - 1 : -1) * size_x * size_y;
		const int idz_p = ((iz == size_z - 1) ? 1 - size_z : 1) * size_x * size_y;
		const int idz_p2 = ((iz >= size_z - 2) ? 2 - size_z : 2) * size_x * size_y;
		const int idz_p3 = ((iz >= size_z - 3) ? 3 - size_z : 3) * size_x * size_y;

		for (int iy = 0; iy < size_y; ++iy) {
			const int idy_m3 = ((iy <= 2) ? size_y - 3 : -3) * size_x;
			const int idy_m2 = ((iy <= 1) ? size_y - 2 : -2) * size_x;
			const int idy_m = ((iy == 0) ? size_y - 1 : -1) * size_x;
			const int idy_p = ((iy == size_y - 1) ? 1 - size_y : 1) * size_x;
			const int idy_p2 = ((iy >= size_y - 2) ? 2 - size_y : 2) * size_x;
			const int idy_p3 = ((iy >= size_y - 3) ? 3 - size_y : 3) * size_x;
			
			const int idx_m3 = -3;
			const int idx_m2 = -2;
			const int idx_m = -1;
			const int idx_p = 1;
			const int idx_p2 = 2;
			const int idx_p3 = 3;
			for (int ix = margin_size; ix < size_x - margin_size; ++ix) {
				const size_t i = ix + size_x * (iy + (size_y * iz));

				double d2psidx2 = ((p[i + idx_p3] + p[i + idx_m3]) - 54.0 * (p[i + idx_p2] + p[i + idx_m2]) + 783.0 * (p[i + idx_p] + p[i + idx_m]) - 1460.0 * p[i]) / (m_dx * m_dx);
				d2psidx2 += ((p[i + idy_p3] + p[i + idy_m3]) - 54.0 * (p[i + idy_p2] + p[i + idy_m2]) + 783.0 * (p[i + idy_p] + p[i + idy_m]) - 1460.0 * p[i]) / (m_dy * m_dy);
				d2psidx2 += ((p[i + idz_p3] + p[i + idz_m3]) - 54.0 * (p[i + idz_p2] + p[i + idz_m2]) + 783.0 * (p[i + idz_p] + p[i + idz_m]) - 1460.0 * p[i]) / (m_dz * m_dz);
				a[i] += -d2psidx2 / (2.0 * 24.0 * 24.0);
			}

			{
				const int ix = 0;
				const size_t i = ix + size_x * (iy + (size_y * iz));
				const size_t im = margin_size*2 * (iy + (size_y * iz));
				double d2psidx2 = ((p[i + idx_p3] + m[margin_size*2-3+ im]) - 54.0 * (p[i + idx_p2] + m[margin_size * 2 - 2 +im]) + 783.0 * (p[i + idx_p] + m[margin_size * 2 - 1 + im]) - 1460.0 * p[i]) / (m_dx * m_dx);
				d2psidx2 += ((p[i + idy_p3] + p[i + idy_m3]) - 54.0 * (p[i + idy_p2] + p[i + idy_m2]) + 783.0 * (p[i + idy_p] + p[i + idy_m]) - 1460.0 * p[i]) / (m_dy * m_dy);
				d2psidx2 += ((p[i + idz_p3] + p[i + idz_m3]) - 54.0 * (p[i + idz_p2] + p[i + idz_m2]) + 783.0 * (p[i + idz_p] + p[i + idz_m]) - 1460.0 * p[i]) / (m_dz * m_dz);
				a[i] += -d2psidx2 / (2.0 * 24.0 * 24.0);
			}

			{
				const int ix = 1;
				const size_t i = ix + size_x * (iy + (size_y * iz));
				const size_t im = margin_size * 2 * (iy + (size_y * iz));
				double d2psidx2 = ((p[i + idx_p3] + m[margin_size * 2 - 2 + im]) - 54.0 * (p[i + idx_p2] + m[margin_size * 2 - 1 + im]) + 783.0 * (p[i + idx_p] + p[i - 1]) - 1460.0 * p[i]) / (m_dx * m_dx);
				d2psidx2 += ((p[i + idy_p3] + p[i + idy_m3]) - 54.0 * (p[i + idy_p2] + p[i + idy_m2]) + 783.0 * (p[i + idy_p] + p[i + idy_m]) - 1460.0 * p[i]) / (m_dy * m_dy);
				d2psidx2 += ((p[i + idz_p3] + p[i + idz_m3]) - 54.0 * (p[i + idz_p2] + p[i + idz_m2]) + 783.0 * (p[i + idz_p] + p[i + idz_m]) - 1460.0 * p[i]) / (m_dz * m_dz);
				a[i] += -d2psidx2 / (2.0 * 24.0 * 24.0);
			}

			{
				const int ix = 2;
				const size_t i = ix + size_x * (iy + (size_y * iz));
				const size_t im = margin_size * 2 * (iy + (size_y * iz));
				double d2psidx2 = ((p[i + idx_p3] + m[margin_size * 2 - 1 + im]) - 54.0 * (p[i + idx_p2] + p[i-2]) + 783.0 * (p[i + idx_p] + p[i - 1]) - 1460.0 * p[i]) / (m_dx * m_dx);
				d2psidx2 += ((p[i + idy_p3] + p[i + idy_m3]) - 54.0 * (p[i + idy_p2] + p[i + idy_m2]) + 783.0 * (p[i + idy_p] + p[i + idy_m]) - 1460.0 * p[i]) / (m_dy * m_dy);
				d2psidx2 += ((p[i + idz_p3] + p[i + idz_m3]) - 54.0 * (p[i + idz_p2] + p[i + idz_m2]) + 783.0 * (p[i + idz_p] + p[i + idz_m]) - 1460.0 * p[i]) / (m_dz * m_dz);
				a[i] += -d2psidx2 / (2.0 * 24.0 * 24.0);
			}

			{
				const int ix = size_x-3;
				const size_t i = ix + size_x * (iy + (size_y * iz));
				const size_t im = margin_size * 2 * (iy + (size_y * iz));
				double d2psidx2 = ((m[0 + im] + p[i + idx_m3]) - 54.0 * (p[i + idx_p2] + p[i + idx_m2]) + 783.0 * (p[i + idx_p] + p[i + idx_m]) - 1460.0 * p[i]) / (m_dx * m_dx);
				d2psidx2 += ((p[i + idy_p3] + p[i + idy_m3]) - 54.0 * (p[i + idy_p2] + p[i + idy_m2]) + 783.0 * (p[i + idy_p] + p[i + idy_m]) - 1460.0 * p[i]) / (m_dy * m_dy);
				d2psidx2 += ((p[i + idz_p3] + p[i + idz_m3]) - 54.0 * (p[i + idz_p2] + p[i + idz_m2]) + 783.0 * (p[i + idz_p] + p[i + idz_m]) - 1460.0 * p[i]) / (m_dz * m_dz);
				a[i] += -d2psidx2 / (2.0 * 24.0 * 24.0);
			}
			{
				const int ix = size_x - 2;
				const size_t i = ix + size_x * (iy + (size_y * iz));
				const size_t im = margin_size * 2 * (iy + (size_y * iz));
				double d2psidx2 = ((m[1 + im] + p[i + idx_m3]) - 54.0 * (m[0 + im] + p[i + idx_m2]) + 783.0 * (p[i + idx_p] + p[i + idx_m]) - 1460.0 * p[i]) / (m_dx * m_dx);
				d2psidx2 += ((p[i + idy_p3] + p[i + idy_m3]) - 54.0 * (p[i + idy_p2] + p[i + idy_m2]) + 783.0 * (p[i + idy_p] + p[i + idy_m]) - 1460.0 * p[i]) / (m_dy * m_dy);
				d2psidx2 += ((p[i + idz_p3] + p[i + idz_m3]) - 54.0 * (p[i + idz_p2] + p[i + idz_m2]) + 783.0 * (p[i + idz_p] + p[i + idz_m]) - 1460.0 * p[i]) / (m_dz * m_dz);
				a[i] += -d2psidx2 / (2.0 * 24.0 * 24.0);
			}
			{
				const int ix = size_x - 2;
				const size_t i = ix + size_x * (iy + (size_y * iz));
				const size_t im = margin_size * 2 * (iy + (size_y * iz));
				double d2psidx2 = ((m[2 + im] + p[i + idx_m3]) - 54.0 * (m[1 + im] + p[i + idx_m2]) + 783.0 * (m[0 + im] + p[i + idx_m]) - 1460.0 * p[i]) / (m_dx * m_dx);
				d2psidx2 += ((p[i + idy_p3] + p[i + idy_m3]) - 54.0 * (p[i + idy_p2] + p[i + idy_m2]) + 783.0 * (p[i + idy_p] + p[i + idy_m]) - 1460.0 * p[i]) / (m_dy * m_dy);
				d2psidx2 += ((p[i + idz_p3] + p[i + idz_m3]) - 54.0 * (p[i + idz_p2] + p[i + idz_m2]) + 783.0 * (p[i + idz_p] + p[i + idz_m]) - 1460.0 * p[i]) / (m_dz * m_dz);
				a[i] += -d2psidx2 / (2.0 * 24.0 * 24.0);
			}
		}
	}
}

inline
void TestDifference() {
	constexpr int STEP = 1000;
	constexpr int Nx = 32;
	constexpr int Ny = 32;
	constexpr int Nz = 32;
	constexpr int margin_size = 3;
	double* a = new double[Nx*Ny*Nz];
	double* b = new double[Nx * Ny * Nz];
	double* c = new double[margin_size*2*Ny * Nz];
	for (int iy = 0; iy < Ny * Nz; ++iy) {
		for (int ix = 0; ix < Nx; ++ix) {
			const int i = ix + Nx * iy;
			a[i] = 0.0;
			b[i] = exp(-(double)(ix - Nx / 2) * (double)(ix - Nx / 2) / 100.0);
		}
		c[0 + margin_size * 2 * iy] = 1.0;
		c[1 + margin_size * 2 * iy] = 2.0;
		c[2 + margin_size * 2 * iy] = 2.0;
		c[3 + margin_size * 2 * iy] = 2.0;
		c[4 + margin_size * 2 * iy] = 2.0;
		c[5 + margin_size * 2 * iy] = 2.0;

	}


	{
		const auto begin_tm = std::chrono::high_resolution_clock::now();
		for (int istep = 0; istep < STEP; ++istep) {
			for (int iy = 0; iy < Ny * Nz; ++iy) {
				Diff2_1_m(a + Nx * iy, b + Nx * iy, c + margin_size * 2 * iy, Nx);
			}
			for (int iy = 0; iy < Ny * Nz; ++iy) {
				Diff2_1_m(b + Nx * iy, a + Nx * iy, c + margin_size * 2 * iy, Nx);
			}
		}
		const auto end_tm = std::chrono::high_resolution_clock::now();
		printf("Diff2_1_m: time = %f [s]\n", DoubleSec(end_tm - begin_tm));
	}


	{
		const auto begin_tm = std::chrono::high_resolution_clock::now();
		for (int istep = 0; istep < STEP; ++istep) {
			for (int iy = 0; iy < Ny * Nz; ++iy) {
				Diff2_1_2(a + Nx * iy, b + Nx * iy, c + margin_size * 2 * iy, Nx);
			}
			for (int iy = 0; iy < Ny * Nz; ++iy) {
				Diff2_1_2(b + Nx * iy, a + Nx * iy, c + margin_size * 2 * iy, Nx);
			}
		}
		const auto end_tm = std::chrono::high_resolution_clock::now();
		printf("Diff2_1_2: time = %f [s]\n", DoubleSec(end_tm - begin_tm));
	}

	{
		const auto begin_tm = std::chrono::high_resolution_clock::now();
		for (int istep = 0; istep < STEP; ++istep) {
			for (int iy = 0; iy < Ny * Nz; ++iy) {
				Diff2_1_0(a + Nx * iy, b + Nx * iy, c + margin_size * 2 * iy, Nx);
			}
			for (int iy = 0; iy < Ny * Nz; ++iy) {
				Diff2_1_0(b + Nx * iy, a + Nx * iy, c + margin_size * 2 * iy, Nx);
			}
		}
		const auto end_tm = std::chrono::high_resolution_clock::now();
		printf("Diff2_1_0: time = %f [s]\n", DoubleSec(end_tm - begin_tm));
	}

	{
		const auto begin_tm = std::chrono::high_resolution_clock::now();
		for (int istep = 0; istep < STEP; ++istep) {
			Diff2_1_p(a, b, c, Nx, Ny, Nz);
			Diff2_1_p(b, a, c, Nx, Ny, Nz);
		}
		const auto end_tm = std::chrono::high_resolution_clock::now();
		printf("Diff2_1_p: time = %f [s]\n", DoubleSec(end_tm - begin_tm));
	}

	{
		const auto begin_tm = std::chrono::high_resolution_clock::now();
		for (int istep = 0; istep < STEP; ++istep) {
			Diff2_2_p(a, b, c, Nx, Ny, Nz);
			Diff2_2_p(b, a, c, Nx, Ny, Nz);
		}
		const auto end_tm = std::chrono::high_resolution_clock::now();
		printf("Diff2_2_p: time = %f [s]\n", DoubleSec(end_tm - begin_tm));
	}

	{
		const auto begin_tm = std::chrono::high_resolution_clock::now();
		for (int istep = 0; istep < STEP; ++istep) {
			Diff2_3_p(a, b, c, Nx, Ny, Nz);
			Diff2_3_p(b, a, c, Nx, Ny, Nz);
		}
		const auto end_tm = std::chrono::high_resolution_clock::now();
		printf("Diff2_3_p: time = %f [s]\n", DoubleSec(end_tm - begin_tm));
	}

	{
		const auto begin_tm = std::chrono::high_resolution_clock::now();
		for (int istep = 0; istep < STEP; ++istep) {
			Diff2_3_0(a, b, c, Nx, Ny, Nz);
			Diff2_3_0(b, a, c, Nx, Ny, Nz);
		}
		const auto end_tm = std::chrono::high_resolution_clock::now();
		printf("Diff2_3_0: time = %f [s]\n", DoubleSec(end_tm - begin_tm));
	}


	{
		const auto begin_tm = std::chrono::high_resolution_clock::now();
		for (int istep = 0; istep < STEP; ++istep) {
			Laplasian2nd_1(GridRange{ 0,0,0, Nx, Ny, Nz }, a, b, -1.0 / 2.0, 0.025, 0.025, 0.025);
			Laplasian2nd_1(GridRange{ 0,0,0, Nx, Ny, Nz }, b, a, -1.0 / 2.0, 0.025, 0.025, 0.025);
		}
		const auto end_tm = std::chrono::high_resolution_clock::now();
		printf("Laplasian2nd: time = %f [s]\n", DoubleSec(end_tm - begin_tm));
	}
	{
		const auto begin_tm = std::chrono::high_resolution_clock::now();
		for (int istep = 0; istep < STEP; ++istep) {
			Laplasian2nd_2(GridRange{ 0,0,0, Nx, Ny, Nz }, a, b, -1.0 / 2.0, 0.025, 0.025, 0.025);
			Laplasian2nd_2(GridRange{ 0,0,0, Nx, Ny, Nz }, b, a, -1.0 / 2.0, 0.025, 0.025, 0.025);
		}
		const auto end_tm = std::chrono::high_resolution_clock::now();
		printf("Laplasian2nd_2: time = %f [s]\n", DoubleSec(end_tm - begin_tm));
	}
	{
		const auto begin_tm = std::chrono::high_resolution_clock::now();
		for (int istep = 0; istep < STEP; ++istep) {
			Laplasian2nd_3(GridRange{ 0,0,0, Nx, Ny, Nz }, a, b, -1.0 / 2.0, 0.025, 0.025, 0.025);
			Laplasian2nd_3(GridRange{ 0,0,0, Nx, Ny, Nz }, b, a, -1.0 / 2.0, 0.025, 0.025, 0.025);
		}
		const auto end_tm = std::chrono::high_resolution_clock::now();
		printf("Laplasian2nd_3: time = %f [s]\n", DoubleSec(end_tm - begin_tm));
	}
	{
		const auto begin_tm = std::chrono::high_resolution_clock::now();
		for (int istep = 0; istep < STEP; ++istep) {
			Laplasian2nd_4(GridRange{ 0,0,0, Nx, Ny, Nz }, a, b, -1.0 / 2.0, 0.025, 0.025, 0.025);
			Laplasian2nd_4(GridRange{ 0,0,0, Nx, Ny, Nz }, b, a, -1.0 / 2.0, 0.025, 0.025, 0.025);
		}
		const auto end_tm = std::chrono::high_resolution_clock::now();
		printf("Laplasian2nd_4: time = %f [s]\n", DoubleSec(end_tm - begin_tm));
	}
	{
		const auto begin_tm = std::chrono::high_resolution_clock::now();
		for (int istep = 0; istep < STEP; ++istep) {
			Laplasian2nd_5(GridRange{ 0,0,0, Nx, Ny, Nz }, a, b, -1.0 / 2.0, 0.025, 0.025, 0.025);
			Laplasian2nd_5(GridRange{ 0,0,0, Nx, Ny, Nz }, b, a, -1.0 / 2.0, 0.025, 0.025, 0.025);
		}
		const auto end_tm = std::chrono::high_resolution_clock::now();
		printf("Laplasian2nd_5: time = %f [s]\n", DoubleSec(end_tm - begin_tm));
	}
	{
		const auto begin_tm = std::chrono::high_resolution_clock::now();
		for (int istep = 0; istep < STEP; ++istep) {
			Laplasian2nd_6(GridRange{ 0,0,0, Nx, Ny, Nz }, a, b, -1.0 / 2.0, 0.025, 0.025, 0.025);
			Laplasian2nd_6(GridRange{ 0,0,0, Nx, Ny, Nz }, b, a, -1.0 / 2.0, 0.025, 0.025, 0.025);
		}
		const auto end_tm = std::chrono::high_resolution_clock::now();
		printf("Laplasian2nd_6: time = %f [s]\n", DoubleSec(end_tm - begin_tm));
	}
	{
		const auto begin_tm = std::chrono::high_resolution_clock::now();
		for (int istep = 0; istep < STEP; ++istep) {
			Laplasian2nd_7(GridRange{ 0,0,0, Nx, Ny, Nz }, a, b, -1.0 / 2.0, 0.025, 0.025, 0.025);
			Laplasian2nd_7(GridRange{ 0,0,0, Nx, Ny, Nz }, b, a, -1.0 / 2.0, 0.025, 0.025, 0.025);
		}
		const auto end_tm = std::chrono::high_resolution_clock::now();
		printf("Laplasian2nd_7: time = %f [s]\n", DoubleSec(end_tm - begin_tm));
	}
	{
		const auto begin_tm = std::chrono::high_resolution_clock::now();
		for (int istep = 0; istep < STEP; ++istep) {
			Laplasian2nd_8(GridRange{ 0,0,0, Nx, Ny, Nz }, a, b, -1.0 / 2.0, 0.025, 0.025, 0.025);
			Laplasian2nd_8(GridRange{ 0,0,0, Nx, Ny, Nz }, b, a, -1.0 / 2.0, 0.025, 0.025, 0.025);
		}
		const auto end_tm = std::chrono::high_resolution_clock::now();
		printf("Laplasian2nd_8: time = %f [s]\n", DoubleSec(end_tm - begin_tm));
	}
	{
		const auto begin_tm = std::chrono::high_resolution_clock::now();
		for (int istep = 0; istep < STEP; ++istep) {
			Laplasian2nd_9(GridRange{ 0,0,0, Nx, Ny, Nz }, a, b, -1.0 / 2.0, 0.025, 0.025, 0.025);
			Laplasian2nd_9(GridRange{ 0,0,0, Nx, Ny, Nz }, b, a, -1.0 / 2.0, 0.025, 0.025, 0.025);
		}
		const auto end_tm = std::chrono::high_resolution_clock::now();
		printf("Laplasian2nd_9: time = %f [s]\n", DoubleSec(end_tm - begin_tm));
	}
	{
		double* b_wide = new double[(Nx + 2) * (Ny + 2) * (Nz + 2)];
		const auto begin_tm = std::chrono::high_resolution_clock::now();
		for (int istep = 0; istep < STEP; ++istep) {
			Laplasian2nd_overspace(GridRange{ 0,0,0, Nx, Ny, Nz }, a, b, b_wide, -1.0 / 2.0, 0.025, 0.025, 0.025);
			Laplasian2nd_overspace(GridRange{ 0,0,0, Nx, Ny, Nz }, b, a, b_wide ,-1.0 / 2.0, 0.025, 0.025, 0.025);
		}
		const auto end_tm = std::chrono::high_resolution_clock::now();
		printf("Laplasian2nd_overspace: time = %f [s]\n", DoubleSec(end_tm - begin_tm));
		delete[] b_wide;
	}
	{
		const auto begin_tm = std::chrono::high_resolution_clock::now();
		for (int istep = 0; istep < STEP; ++istep) {
			Laplasian4th(GridRange{ 0,0,0, Nx, Ny, Nz }, a, b, -1.0 / 2.0, 0.025, 0.025, 0.025);
			Laplasian4th(GridRange{ 0,0,0, Nx, Ny, Nz }, b, a, -1.0 / 2.0, 0.025, 0.025, 0.025);
		}
		const auto end_tm = std::chrono::high_resolution_clock::now();
		printf("Laplasian4th: time = %f [s]\n", DoubleSec(end_tm - begin_tm));
	}
	{
		const auto begin_tm = std::chrono::high_resolution_clock::now();
		for (int istep = 0; istep < STEP; ++istep) {
			Laplasian4th_d(GridRange{ 0,0,0, Nx, Ny, Nz }, a, b, -1.0 / 2.0, 0.025, 0.025, 0.025);
			Laplasian4th_d(GridRange{ 0,0,0, Nx, Ny, Nz }, b, a, -1.0 / 2.0, 0.025, 0.025, 0.025);
		}
		const auto end_tm = std::chrono::high_resolution_clock::now();
		printf("Laplasian4th_d: time = %f [s]\n", DoubleSec(end_tm - begin_tm));
	}
	{
		double* b_wide = new double[(Nx + 6) * (Ny + 6)* (Nz + 6)];
		const auto begin_tm = std::chrono::high_resolution_clock::now();
		for (int istep = 0; istep < STEP; ++istep) {
			Laplasian4th_overspace(GridRange{ 0,0,0, Nx, Ny, Nz }, a, b, b_wide, -1.0 / 2.0, 0.025, 0.025, 0.025);
			Laplasian4th_overspace(GridRange{ 0,0,0, Nx, Ny, Nz }, b, a, b_wide, -1.0 / 2.0, 0.025, 0.025, 0.025);
		}
		const auto end_tm = std::chrono::high_resolution_clock::now();
		printf("Laplasian4th_overspace: time = %f [s]\n", DoubleSec(end_tm - begin_tm));
		delete[] b_wide;
	}
	{
		double* b_wide = new double[(Nx + 6) * (Ny + 6) * (Nz + 6)];
		const auto begin_tm = std::chrono::high_resolution_clock::now();
		for (int istep = 0; istep < STEP; ++istep) {
			Laplasian4th_overspace2(GridRange{ 0,0,0, Nx, Ny, Nz }, a, b, b_wide, -1.0 / 2.0, 0.025, 0.025, 0.025);
			Laplasian4th_overspace2(GridRange{ 0,0,0, Nx, Ny, Nz }, b, a, b_wide, -1.0 / 2.0, 0.025, 0.025, 0.025);
		}
		const auto end_tm = std::chrono::high_resolution_clock::now();
		printf("Laplasian4th_overspace2: time = %f [s]\n", DoubleSec(end_tm - begin_tm));
		delete[] b_wide;
	}
	/*
	{
		const auto begin_tm = std::chrono::high_resolution_clock::now();
		for (int istep = 0; istep < STEP; ++istep) {
			Laplasian4th_8(GridRange{ 0,0,0, Nx, Ny, Nz }, a, b, -1.0 / 2.0, 0.025, 0.025, 0.025);
			Laplasian4th_8(GridRange{ 0,0,0, Nx, Ny, Nz }, b, a, -1.0 / 2.0, 0.025, 0.025, 0.025);
		}
		const auto end_tm = std::chrono::high_resolution_clock::now();
		printf("Laplasian4th_8: time = %f [s]\n", DoubleSec(end_tm - begin_tm));		
	}
	*/
}


int main(int argc, char* argv[]) {
	TestDifference();
	return 0;
}
