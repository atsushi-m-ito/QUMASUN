#pragma once
#ifdef USE_MPI
#include <mpi.h>
#include "mpi_helper.h"
#endif
#include "GridRange.h"
#include <vector>

#include "GridDifference2nd.h"


//3次元デカルト座標グリッドの計算のためのループを司る機能//
//MPI並列時にはリダクションも提供する



template <class T>
void Laplasian4th_base(const GridRange& grid, T* out, const T* p, const double coef, const double dx, const double dy, const double dz) {
	const int size_x = grid.SizeX();
	const int size_y = grid.SizeY();
	const int size_z = grid.SizeZ();

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
			for (int ix = 0; ix < size_x; ++ix) {
				const int idx_m3 = ((ix <= 2) ? size_x - 3 : -3);
				const int idx_m2 = ((ix <= 1) ? size_x - 2 : -2);
				const int idx_m = ((ix == 0) ? size_x - 1 : -1);
				const int idx_p = ((ix == size_x - 1) ? 1 - size_x : 1);
				const int idx_p2 = ((ix >= size_x - 2) ? 2 - size_x : 2);
				const int idx_p3 = ((ix >= size_x - 3) ? 3 - size_x : 3);
				const int64_t i = ix + size_x * (iy + (size_y * iz));

				auto d2psidx2 = ((p[i + idx_p3] + p[i + idx_m3]) - 54.0 * (p[i + idx_p2] + p[i + idx_m2]) + 783.0 * (p[i + idx_p] + p[i + idx_m]) - 1460.0 * p[i]) / (dx * dx);
				d2psidx2 += ((p[i + idy_p3] + p[i + idy_m3]) - 54.0 * (p[i + idy_p2] + p[i + idy_m2]) + 783.0 * (p[i + idy_p] + p[i + idy_m]) - 1460.0 * p[i]) / (dy * dy);
				d2psidx2 += ((p[i + idz_p3] + p[i + idz_m3]) - 54.0 * (p[i + idz_p2] + p[i + idz_m2]) + 783.0 * (p[i + idz_p] + p[i + idz_m]) - 1460.0 * p[i]) / (dz * dz);
				out[i] += coef * d2psidx2 / (24.0 * 24.0);

			}
		}
	}
}

template <class T>
void Laplasian4th_d(const GridRange& grid, T* out, const T* p, const double coef, const double dx, const double dy, const double dz) {
	const int size_x = grid.SizeX();
	const int size_y = grid.SizeY();
	const int size_z = grid.SizeZ();
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
			
			{
				const int ix = 0;
				const int idx_m3 = ((ix <= 2) ? size_x - 3 : -3);
				const int idx_m2 = ((ix <= 1) ? size_x - 2 : -2);
				const int idx_m = ((ix == 0) ? size_x - 1 : -1);
				const int idx_p = ((ix == size_x - 1) ? 1 - size_x : 1);
				const int idx_p2 = ((ix >= size_x - 2) ? 2 - size_x : 2);
				const int idx_p3 = ((ix >= size_x - 3) ? 3 - size_x : 3);
				const int64_t i = ix + size_x * (iy + (size_y * iz));

				auto d2psidx2 = ((p[i + idx_p3] + p[i + idx_m3]) - 54.0 * (p[i + idx_p2] + p[i + idx_m2]) + 783.0 * (p[i + idx_p] + p[i + idx_m]) - 1460.0 * p[i]) / (dx * dx);
				d2psidx2 += ((p[i + idy_p3] + p[i + idy_m3]) - 54.0 * (p[i + idy_p2] + p[i + idy_m2]) + 783.0 * (p[i + idy_p] + p[i + idy_m]) - 1460.0 * p[i]) / (dy * dy);
				d2psidx2 += ((p[i + idz_p3] + p[i + idz_m3]) - 54.0 * (p[i + idz_p2] + p[i + idz_m2]) + 783.0 * (p[i + idz_p] + p[i + idz_m]) - 1460.0 * p[i]) / (dz * dz);
				out[i] += coef * d2psidx2 / (24.0 * 24.0);

			}
			{
				const int ix = 1;
				const int idx_m3 = ((ix <= 2) ? size_x - 3 : -3);
				const int idx_m2 = ((ix <= 1) ? size_x - 2 : -2);
				const int idx_m = ((ix == 0) ? size_x - 1 : -1);
				const int idx_p = ((ix == size_x - 1) ? 1 - size_x : 1);
				const int idx_p2 = ((ix >= size_x - 2) ? 2 - size_x : 2);
				const int idx_p3 = ((ix >= size_x - 3) ? 3 - size_x : 3);
				const int64_t i = ix + size_x * (iy + (size_y * iz));

				auto d2psidx2 = ((p[i + idx_p3] + p[i + idx_m3]) - 54.0 * (p[i + idx_p2] + p[i + idx_m2]) + 783.0 * (p[i + idx_p] + p[i + idx_m]) - 1460.0 * p[i]) / (dx * dx);
				d2psidx2 += ((p[i + idy_p3] + p[i + idy_m3]) - 54.0 * (p[i + idy_p2] + p[i + idy_m2]) + 783.0 * (p[i + idy_p] + p[i + idy_m]) - 1460.0 * p[i]) / (dy * dy);
				d2psidx2 += ((p[i + idz_p3] + p[i + idz_m3]) - 54.0 * (p[i + idz_p2] + p[i + idz_m2]) + 783.0 * (p[i + idz_p] + p[i + idz_m]) - 1460.0 * p[i]) / (dz * dz);
				out[i] += coef * d2psidx2 / (24.0 * 24.0);

			}
			{
				const int ix = 2;
				const int idx_m3 = ((ix <= 2) ? size_x - 3 : -3);
				const int idx_m2 = ((ix <= 1) ? size_x - 2 : -2);
				const int idx_m = ((ix == 0) ? size_x - 1 : -1);
				const int idx_p = ((ix == size_x - 1) ? 1 - size_x : 1);
				const int idx_p2 = ((ix >= size_x - 2) ? 2 - size_x : 2);
				const int idx_p3 = ((ix >= size_x - 3) ? 3 - size_x : 3);
				const int64_t i = ix + size_x * (iy + (size_y * iz));

				auto d2psidx2 = ((p[i + idx_p3] + p[i + idx_m3]) - 54.0 * (p[i + idx_p2] + p[i + idx_m2]) + 783.0 * (p[i + idx_p] + p[i + idx_m]) - 1460.0 * p[i]) / (dx * dx);
				d2psidx2 += ((p[i + idy_p3] + p[i + idy_m3]) - 54.0 * (p[i + idy_p2] + p[i + idy_m2]) + 783.0 * (p[i + idy_p] + p[i + idy_m]) - 1460.0 * p[i]) / (dy * dy);
				d2psidx2 += ((p[i + idz_p3] + p[i + idz_m3]) - 54.0 * (p[i + idz_p2] + p[i + idz_m2]) + 783.0 * (p[i + idz_p] + p[i + idz_m]) - 1460.0 * p[i]) / (dz * dz);
				out[i] += coef * d2psidx2 / (24.0 * 24.0);

			}
			for (int ix = 3; ix < size_x - 3; ++ix) {
				const int idx_m3 = -3;
				const int idx_m2 = -2;
				const int idx_m = -1;
				const int idx_p = 1;
				const int idx_p2 = 2;
				const int idx_p3 = 3;
				const int64_t i = ix + size_x * (iy + (size_y * iz));

				auto d2psidx2 = ((p[i + idx_p3] + p[i + idx_m3]) - 54.0 * (p[i + idx_p2] + p[i + idx_m2]) + 783.0 * (p[i + idx_p] + p[i + idx_m]) - 1460.0 * p[i]) / (dx * dx);
				d2psidx2 += ((p[i + idy_p3] + p[i + idy_m3]) - 54.0 * (p[i + idy_p2] + p[i + idy_m2]) + 783.0 * (p[i + idy_p] + p[i + idy_m]) - 1460.0 * p[i]) / (dy * dy);
				d2psidx2 += ((p[i + idz_p3] + p[i + idz_m3]) - 54.0 * (p[i + idz_p2] + p[i + idz_m2]) + 783.0 * (p[i + idz_p] + p[i + idz_m]) - 1460.0 * p[i]) / (dz * dz);
				out[i] += coef * d2psidx2 / (24.0 * 24.0);

			}
			{
				const int ix = size_z - 3;
				const int idx_m3 = ((ix <= 2) ? size_x - 3 : -3);
				const int idx_m2 = ((ix <= 1) ? size_x - 2 : -2);
				const int idx_m = ((ix == 0) ? size_x - 1 : -1);
				const int idx_p = ((ix == size_x - 1) ? 1 - size_x : 1);
				const int idx_p2 = ((ix >= size_x - 2) ? 2 - size_x : 2);
				const int idx_p3 = ((ix >= size_x - 3) ? 3 - size_x : 3);
				const int64_t i = ix + size_x * (iy + (size_y * iz));

				auto d2psidx2 = ((p[i + idx_p3] + p[i + idx_m3]) - 54.0 * (p[i + idx_p2] + p[i + idx_m2]) + 783.0 * (p[i + idx_p] + p[i + idx_m]) - 1460.0 * p[i]) / (dx * dx);
				d2psidx2 += ((p[i + idy_p3] + p[i + idy_m3]) - 54.0 * (p[i + idy_p2] + p[i + idy_m2]) + 783.0 * (p[i + idy_p] + p[i + idy_m]) - 1460.0 * p[i]) / (dy * dy);
				d2psidx2 += ((p[i + idz_p3] + p[i + idz_m3]) - 54.0 * (p[i + idz_p2] + p[i + idz_m2]) + 783.0 * (p[i + idz_p] + p[i + idz_m]) - 1460.0 * p[i]) / (dz * dz);
				out[i] += coef * d2psidx2 / (24.0 * 24.0);

			}
			{
				const int ix = size_z - 2;
				const int idx_m3 = ((ix <= 2) ? size_x - 3 : -3);
				const int idx_m2 = ((ix <= 1) ? size_x - 2 : -2);
				const int idx_m = ((ix == 0) ? size_x - 1 : -1);
				const int idx_p = ((ix == size_x - 1) ? 1 - size_x : 1);
				const int idx_p2 = ((ix >= size_x - 2) ? 2 - size_x : 2);
				const int idx_p3 = ((ix >= size_x - 3) ? 3 - size_x : 3);
				const int64_t i = ix + size_x * (iy + (size_y * iz));

				auto d2psidx2 = ((p[i + idx_p3] + p[i + idx_m3]) - 54.0 * (p[i + idx_p2] + p[i + idx_m2]) + 783.0 * (p[i + idx_p] + p[i + idx_m]) - 1460.0 * p[i]) / (dx * dx);
				d2psidx2 += ((p[i + idy_p3] + p[i + idy_m3]) - 54.0 * (p[i + idy_p2] + p[i + idy_m2]) + 783.0 * (p[i + idy_p] + p[i + idy_m]) - 1460.0 * p[i]) / (dy * dy);
				d2psidx2 += ((p[i + idz_p3] + p[i + idz_m3]) - 54.0 * (p[i + idz_p2] + p[i + idz_m2]) + 783.0 * (p[i + idz_p] + p[i + idz_m]) - 1460.0 * p[i]) / (dz * dz);
				out[i] += coef * d2psidx2 / (24.0 * 24.0);

			}
			{
				const int ix = size_z-1;
				const int idx_m3 = ((ix <= 2) ? size_x - 3 : -3);
				const int idx_m2 = ((ix <= 1) ? size_x - 2 : -2);
				const int idx_m = ((ix == 0) ? size_x - 1 : -1);
				const int idx_p = ((ix == size_x - 1) ? 1 - size_x : 1);
				const int idx_p2 = ((ix >= size_x - 2) ? 2 - size_x : 2);
				const int idx_p3 = ((ix >= size_x - 3) ? 3 - size_x : 3);
				const int64_t i = ix + size_x * (iy + (size_y * iz));

				auto d2psidx2 = ((p[i + idx_p3] + p[i + idx_m3]) - 54.0 * (p[i + idx_p2] + p[i + idx_m2]) + 783.0 * (p[i + idx_p] + p[i + idx_m]) - 1460.0 * p[i]) / (dx * dx);
				d2psidx2 += ((p[i + idy_p3] + p[i + idy_m3]) - 54.0 * (p[i + idy_p2] + p[i + idy_m2]) + 783.0 * (p[i + idy_p] + p[i + idy_m]) - 1460.0 * p[i]) / (dy * dy);
				d2psidx2 += ((p[i + idz_p3] + p[i + idz_m3]) - 54.0 * (p[i + idz_p2] + p[i + idz_m2]) + 783.0 * (p[i + idz_p] + p[i + idz_m]) - 1460.0 * p[i]) / (dz * dz);
				out[i] += coef * d2psidx2 / (24.0 * 24.0);

			}
		}
	}
}

#if 0
template <class T>
void Laplasian4th_overspace(const GridRange& grid, T* out, const T* src, T* p, const double coef, const double dx, const double dy, const double dz) {
	const double c0 = -1460.0 / (dx * dx) - 1460.0 / (dy * dy) - 1460.0 / (dz * dz); 
	const int Nx = grid.SizeX();
	const int Ny = grid.SizeY();
	const int Nz = grid.SizeZ();

	{
#pragma ivdep
		for (int iz = 0; iz < Nz; ++iz) {
			for (int iy = 0; iy < Ny; ++iy) {
				for (int ix = 0; ix < Nx; ++ix) {
					p[(ix + 3) + (Nx + 6) * (iy + 3 + (Ny + 6) * (iz + 3))] = src[ix + Nx * (iy + Ny * iz)];
				}
			}
		}
#pragma ivdep
		for (int iz = 0; iz < Nz; ++iz) {
			for (int iy = 0; iy < Ny; ++iy) {
				p[0 + (Nx + 6) * (iy + 3 + (Ny + 6) * (iz + 3))] = src[(Nx - 3) + Nx * (iy + Ny * iz)];
				p[1 + (Nx + 6) * (iy + 3 + (Ny + 6) * (iz + 3))] = src[(Nx - 2) + Nx * (iy + Ny * iz)];
				p[2 + (Nx + 6) * (iy + 3 + (Ny + 6) * (iz + 3))] = src[(Nx - 1) + Nx * (iy + Ny * iz)];
				p[Nx + 5 + (Nx + 6) * (iy + 3 + (Ny + 6) * (iz + 3))] = src[2 + Nx * (iy + Ny * iz)];
				p[Nx + 4 + (Nx + 6) * (iy + 3 + (Ny + 6) * (iz + 3))] = src[1 + Nx * (iy + Ny * iz)];
				p[Nx + 3 + (Nx + 6) * (iy + 3 + (Ny + 6) * (iz + 3))] = src[0 + Nx * (iy + Ny * iz)];
			}
		}
#pragma ivdep
		for (int iz = 0; iz < Nz; ++iz) {
			for (int ix = 0; ix < Nx; ++ix) {
				p[(ix + 3) + (Nx + 6) * (0 + (Ny + 6) * (iz + 3))] = src[ix + Nx * (Ny - 3 + Ny * iz)];
				p[(ix + 3) + (Nx + 6) * (1 + (Ny + 6) * (iz + 3))] = src[ix + Nx * (Ny - 2 + Ny * iz)];
				p[(ix + 3) + (Nx + 6) * (2 + (Ny + 6) * (iz + 3))] = src[ix + Nx * (Ny - 1 + Ny * iz)];
				p[(ix + 3) + (Nx + 6) * (Ny + 5 + (Ny + 6) * (iz + 3))] = src[ix + Nx * (2 + Ny * iz)];
				p[(ix + 3) + (Nx + 6) * (Ny + 4 + (Ny + 6) * (iz + 3))] = src[ix + Nx * (1 + Ny * iz)];
				p[(ix + 3) + (Nx + 6) * (Ny + 3 + (Ny + 6) * (iz + 3))] = src[ix + Nx * (0 + Ny * iz)];
			}
		}
#pragma ivdep
		for (int iy = 0; iy < Ny; ++iy) {
			for (int ix = 0; ix < Nx; ++ix) {
				p[(ix + 3) + (Nx + 6) * (iy + 3 + (Ny + 6) * (0))] = src[ix + Nx * (iy + Ny * (Nz - 3))];
				p[(ix + 3) + (Nx + 6) * (iy + 3 + (Ny + 6) * (1))] = src[ix + Nx * (iy + Ny * (Nz - 2))];
				p[(ix + 3) + (Nx + 6) * (iy + 3 + (Ny + 6) * (2))] = src[ix + Nx * (iy + Ny * (Nz - 1))];
			}
		}
		for (int iy = 0; iy < Ny; ++iy) {
			for (int ix = 0; ix < Nx; ++ix) {
				p[(ix + 3) + (Nx + 6) * (iy + 3 + (Ny + 6) * (Nz + 3))] = src[ix + Nx * (iy + Ny * 0)];
				p[(ix + 3) + (Nx + 6) * (iy + 3 + (Ny + 6) * (Nz + 4))] = src[ix + Nx * (iy + Ny * 1)];
				p[(ix + 3) + (Nx + 6) * (iy + 3 + (Ny + 6) * (Nz + 5))] = src[ix + Nx * (iy + Ny * 2)];
			}
		}
	}


	const int size_x = grid.SizeX() + 6;
	const int size_y = grid.SizeY() + 6;
	const int size_z = grid.SizeZ() + 6;

	for (int iz = 3; iz < size_z-3; ++iz) {
		const int idz_m3 = -3 * size_x * size_y;
		const int idz_m2 = -2 * size_x * size_y;
		const int idz_m = -1 * size_x * size_y;
		const int idz_p = 1 * size_x * size_y;
		const int idz_p2 = 2 * size_x * size_y;
		const int idz_p3 = 3 * size_x * size_y;

		for (int iy = 3; iy < size_y-3; ++iy) {
			const int idy_m3 = -3 * size_x;
			const int idy_m2 = -2 * size_x;
			const int idy_m = -1 * size_x;
			const int idy_p = 1 * size_x;
			const int idy_p2 = 2 * size_x;
			const int idy_p3 = 3 * size_x;
#pragma ivdep
			for (int ix = 3; ix < size_x-3; ++ix) {
				const int idx_m3 = -3;
				const int idx_m2 = -2;
				const int idx_m = -1;
				const int idx_p = 1;
				const int idx_p2 = 2;
				const int idx_p3 = 3;
				const int64_t i = ix + size_x * (iy + (size_y * iz));
				const int64_t io = (ix-3) + Nx * (iy - 3 + (Ny * (iz-3)));
				/*
				auto d2psidx2 = ((p[i + idx_p3] + p[i + idx_m3]) - 54.0 * (p[i + idx_p2] + p[i + idx_m2]) + 783.0 * (p[i + idx_p] + p[i + idx_m]) - 1460.0 * p[i]) / (dx * dx);
				d2psidx2 += ((p[i + idy_p3] + p[i + idy_m3]) - 54.0 * (p[i + idy_p2] + p[i + idy_m2]) + 783.0 * (p[i + idy_p] + p[i + idy_m]) - 1460.0 * p[i]) / (dy * dy);
				d2psidx2 += ((p[i + idz_p3] + p[i + idz_m3]) - 54.0 * (p[i + idz_p2] + p[i + idz_m2]) + 783.0 * (p[i + idz_p] + p[i + idz_m]) - 1460.0 * p[i]) / (dz * dz);
				*/
				auto d2psidx2 = c0 * p[i];
				d2psidx2 += ((p[i + idx_p3] + p[i + idx_m3]) - 54.0 * (p[i + idx_p2] + p[i + idx_m2]) + 783.0 * (p[i + idx_p] + p[i + idx_m]) ) / (dx * dx);
				d2psidx2 += ((p[i + idy_p3] + p[i + idy_m3]) - 54.0 * (p[i + idy_p2] + p[i + idy_m2]) + 783.0 * (p[i + idy_p] + p[i + idy_m]) ) / (dy * dy);
				d2psidx2 += ((p[i + idz_p3] + p[i + idz_m3]) - 54.0 * (p[i + idz_p2] + p[i + idz_m2]) + 783.0 * (p[i + idz_p] + p[i + idz_m]) ) / (dz * dz);
				out[io] += coef * d2psidx2 / (24.0 * 24.0);

			}
		}
	}
}

#else
//faster than above code//

template <class T>
void Laplasian4th_overspace(const GridRange& grid, T* out, const T* src, T* p, const double coef, const double dx, const double dy, const double dz) {
	const int Nx = grid.SizeX();
	const int Ny = grid.SizeY();
	const int Nz = grid.SizeZ();

	auto pow2 = [](T a) {return a * a; };
	const double c0 = -1460.0 * coef * (1.0 / pow2(24.0 * dx) + 1.0 / pow2(24.0 * dy) + 1.0 / pow2(24.0 * dz));
	const double c3x = coef / pow2(24.0 * dx);
	const double c2x = -54.0 * coef / pow2(24.0 * dx);
	const double c1x = 783.0 * coef / pow2(24.0 * dx);
	const double c3y = coef / pow2(24.0 * dy);
	const double c2y = -54.0 * coef / pow2(24.0 * dy);
	const double c1y = 783.0 * coef / pow2(24.0 * dy);
	const double c3z = coef / pow2(24.0 * dz);
	const double c2z = -54.0 * coef / pow2(24.0 * dz);
	const double c1z = 783.0 * coef / pow2(24.0 * dz);


	{
#pragma ivdep
		for (int iz = 0; iz < Nz; ++iz) {
			for (int iy = 0; iy < Ny; ++iy) {
				for (int ix = 0; ix < Nx; ++ix) {
					p[(ix + 3) + (Nx + 6) * (iy + 3 + (Ny + 6) * (iz + 3))] = src[ix + Nx * (iy + Ny * iz)];
				}
			}
		}
#pragma ivdep
		for (int iz = 0; iz < Nz; ++iz) {
			for (int iy = 0; iy < Ny; ++iy) {
				p[0 + (Nx + 6) * (iy + 3 + (Ny + 6) * (iz + 3))] = src[(Nx - 3) + Nx * (iy + Ny * iz)];
				p[1 + (Nx + 6) * (iy + 3 + (Ny + 6) * (iz + 3))] = src[(Nx - 2) + Nx * (iy + Ny * iz)];
				p[2 + (Nx + 6) * (iy + 3 + (Ny + 6) * (iz + 3))] = src[(Nx - 1) + Nx * (iy + Ny * iz)];
				p[Nx + 5 + (Nx + 6) * (iy + 3 + (Ny + 6) * (iz + 3))] = src[2 + Nx * (iy + Ny * iz)];
				p[Nx + 4 + (Nx + 6) * (iy + 3 + (Ny + 6) * (iz + 3))] = src[1 + Nx * (iy + Ny * iz)];
				p[Nx + 3 + (Nx + 6) * (iy + 3 + (Ny + 6) * (iz + 3))] = src[0 + Nx * (iy + Ny * iz)];
			}
		}
#pragma ivdep
		for (int iz = 0; iz < Nz; ++iz) {
			for (int ix = 0; ix < Nx; ++ix) {
				p[(ix + 3) + (Nx + 6) * (0 + (Ny + 6) * (iz + 3))] = src[ix + Nx * (Ny - 3 + Ny * iz)];
				p[(ix + 3) + (Nx + 6) * (1 + (Ny + 6) * (iz + 3))] = src[ix + Nx * (Ny - 2 + Ny * iz)];
				p[(ix + 3) + (Nx + 6) * (2 + (Ny + 6) * (iz + 3))] = src[ix + Nx * (Ny - 1 + Ny * iz)];
				p[(ix + 3) + (Nx + 6) * (Ny + 5 + (Ny + 6) * (iz + 3))] = src[ix + Nx * (2 + Ny * iz)];
				p[(ix + 3) + (Nx + 6) * (Ny + 4 + (Ny + 6) * (iz + 3))] = src[ix + Nx * (1 + Ny * iz)];
				p[(ix + 3) + (Nx + 6) * (Ny + 3 + (Ny + 6) * (iz + 3))] = src[ix + Nx * (0 + Ny * iz)];
			}
		}
#pragma ivdep
		for (int iy = 0; iy < Ny; ++iy) {
			for (int ix = 0; ix < Nx; ++ix) {
				p[(ix + 3) + (Nx + 6) * (iy + 3 + (Ny + 6) * (0))] = src[ix + Nx * (iy + Ny * (Nz - 3))];
				p[(ix + 3) + (Nx + 6) * (iy + 3 + (Ny + 6) * (1))] = src[ix + Nx * (iy + Ny * (Nz - 2))];
				p[(ix + 3) + (Nx + 6) * (iy + 3 + (Ny + 6) * (2))] = src[ix + Nx * (iy + Ny * (Nz - 1))];
			}
		}
		for (int iy = 0; iy < Ny; ++iy) {
			for (int ix = 0; ix < Nx; ++ix) {
				p[(ix + 3) + (Nx + 6) * (iy + 3 + (Ny + 6) * (Nz + 3))] = src[ix + Nx * (iy + Ny * 0)];
				p[(ix + 3) + (Nx + 6) * (iy + 3 + (Ny + 6) * (Nz + 4))] = src[ix + Nx * (iy + Ny * 1)];
				p[(ix + 3) + (Nx + 6) * (iy + 3 + (Ny + 6) * (Nz + 5))] = src[ix + Nx * (iy + Ny * 2)];
			}
		}
	}


	const int size_x = grid.SizeX() + 6;
	const int size_y = grid.SizeY() + 6;
	const int size_z = grid.SizeZ() + 6;

	for (int iz = 3; iz < size_z - 3; ++iz) {
		const int idz_m3 = -3 * size_x * size_y;
		const int idz_m2 = -2 * size_x * size_y;
		const int idz_m = -1 * size_x * size_y;
		const int idz_p = 1 * size_x * size_y;
		const int idz_p2 = 2 * size_x * size_y;
		const int idz_p3 = 3 * size_x * size_y;

		for (int iy = 3; iy < size_y - 3; ++iy) {
			const int idy_m3 = -3 * size_x;
			const int idy_m2 = -2 * size_x;
			const int idy_m = -1 * size_x;
			const int idy_p = 1 * size_x;
			const int idy_p2 = 2 * size_x;
			const int idy_p3 = 3 * size_x;
#pragma ivdep
			for (int ix = 3; ix < size_x - 3; ++ix) {
				const int idx_m3 = -3;
				const int idx_m2 = -2;
				const int idx_m = -1;
				const int idx_p = 1;
				const int idx_p2 = 2;
				const int idx_p3 = 3;
				const int64_t i = ix + size_x * (iy + (size_y * iz));
				const int64_t io = (ix - 3) + Nx * (iy - 3 + (Ny * (iz - 3)));

#if 1
				auto d2psidx2 = c3x * ((p[i + idx_p3] + p[i + idx_m3]) - 54.0 * (p[i + idx_p2] + p[i + idx_m2]) + 783.0 * (p[i + idx_p] + p[i + idx_m]) - 1460.0 * p[i])
				+c3y * ((p[i + idy_p3] + p[i + idy_m3]) - 54.0 * (p[i + idy_p2] + p[i + idy_m2]) + 783.0 * (p[i + idy_p] + p[i + idy_m]) - 1460.0 * p[i])
					+ c3z * ((p[i + idz_p3] + p[i + idz_m3]) - 54.0 * (p[i + idz_p2] + p[i + idz_m2]) + 783.0 * (p[i + idz_p] + p[i + idz_m]) - 1460.0 * p[i]);
				out[io] += d2psidx2;
#else
				auto d2psidx2 = c0 * p[i];
				d2psidx2 += c3x*(p[i + idx_p3] + p[i + idx_m3]) +c2x* (p[i + idx_p2] + p[i + idx_m2]) + c1x * (p[i + idx_p] + p[i + idx_m]);
				d2psidx2 += c3y*(p[i + idy_p3] + p[i + idy_m3]) +c2y* (p[i + idy_p2] + p[i + idy_m2]) + c1y * (p[i + idy_p] + p[i + idy_m]);
				d2psidx2 += c3z*(p[i + idz_p3] + p[i + idz_m3]) +c2z* (p[i + idz_p2] + p[i + idz_m2]) + c1z * (p[i + idz_p] + p[i + idz_m]);
				out[io] += d2psidx2;
#endif
			}
		}
	}
}
#endif

#if 1

#define PasteHalo_f3_defined


template <class T, class FUNC_XL, class FUNC_XR, class FUNC_YL, class FUNC_YR, class FUNC_ZL, class FUNC_ZR>
inline
void PasteHalo_f3(const GridRange& grid, T* p, const T* src, FUNC_XL&& func_xl, FUNC_XR&& func_xr, FUNC_YL&& func_yl, FUNC_YR&& func_yr, FUNC_ZL&& func_zl, FUNC_ZR&& func_zr) {
	const int Nx = grid.SizeX();
	const int Ny = grid.SizeY();
	const int Nz = grid.SizeZ();


	{
#pragma ivdep
		for (int iz = 0; iz < Nz; ++iz) {
			for (int iy = 0; iy < Ny; ++iy) {
				for (int ix = 0; ix < Nx; ++ix) {
					p[(ix + 3) + (Nx + 6) * (iy + 3 + (Ny + 6) * (iz + 3))] = src[ix + Nx * (iy + Ny * iz)];
				}
			}
		}
#pragma ivdep
		for (int iz = 0; iz < Nz; ++iz) {
			for (int iy = 0; iy < Ny; ++iy) {
				p[0 + (Nx + 6) * (iy + 3 + (Ny + 6) * (iz + 3))] = func_xl(-3, iy, iz);
				p[1 + (Nx + 6) * (iy + 3 + (Ny + 6) * (iz + 3))] = func_xl(-2, iy, iz);
				p[2 + (Nx + 6) * (iy + 3 + (Ny + 6) * (iz + 3))] = func_xl(-1, iy, iz);
				p[Nx + 5 + (Nx + 6) * (iy + 3 + (Ny + 6) * (iz + 3))] = func_xr(Nx + 2, iy, iz);
				p[Nx + 4 + (Nx + 6) * (iy + 3 + (Ny + 6) * (iz + 3))] = func_xr(Nx + 1, iy, iz);
				p[Nx + 3 + (Nx + 6) * (iy + 3 + (Ny + 6) * (iz + 3))] = func_xr(Nx + 0, iy, iz);
			}
		}
#pragma ivdep
		for (int iz = 0; iz < Nz; ++iz) {
			for (int ix = 0; ix < Nx; ++ix) {
				p[(ix + 3) + (Nx + 6) * (0 + (Ny + 6) * (iz + 3))] = func_yl(ix, -3, iz);
				p[(ix + 3) + (Nx + 6) * (1 + (Ny + 6) * (iz + 3))] = func_yl(ix, -2, iz);
				p[(ix + 3) + (Nx + 6) * (2 + (Ny + 6) * (iz + 3))] = func_yl(ix, -1, iz);
				p[(ix + 3) + (Nx + 6) * (Ny + 5 + (Ny + 6) * (iz + 3))] = func_yr(ix, Ny + 2, iz);
				p[(ix + 3) + (Nx + 6) * (Ny + 4 + (Ny + 6) * (iz + 3))] = func_yr(ix, Ny + 1, iz);
				p[(ix + 3) + (Nx + 6) * (Ny + 3 + (Ny + 6) * (iz + 3))] = func_yr(ix, Ny + 0, iz);
			}
		}
#pragma ivdep
		for (int iy = 0; iy < Ny; ++iy) {
			for (int ix = 0; ix < Nx; ++ix) {
				p[(ix + 3) + (Nx + 6) * (iy + 3 + (Ny + 6) * (0))] = func_zl(ix, iy, -3);
				p[(ix + 3) + (Nx + 6) * (iy + 3 + (Ny + 6) * (1))] = func_zl(ix, iy, -2);
				p[(ix + 3) + (Nx + 6) * (iy + 3 + (Ny + 6) * (2))] = func_zl(ix, iy, -1);
			}
		}
		for (int iy = 0; iy < Ny; ++iy) {
			for (int ix = 0; ix < Nx; ++ix) {
				p[(ix + 3) + (Nx + 6) * (iy + 3 + (Ny + 6) * (Nz + 3))] = func_zr(ix, iy, Nz + 0);
				p[(ix + 3) + (Nx + 6) * (iy + 3 + (Ny + 6) * (Nz + 4))] = func_zr(ix, iy, Nz + 1);
				p[(ix + 3) + (Nx + 6) * (iy + 3 + (Ny + 6) * (Nz + 5))] = func_zr(ix, iy, Nz + 2);
			}
		}
	}
}

template <class T>
void Laplasian4th_halo_f(const GridRange& grid, T* out, const T* p, double coef, double dx, double dy, double dz) {
	
	const int Nx = grid.SizeX();
	const int Ny = grid.SizeY();
	const int Nz = grid.SizeZ();

	auto pow2 = [](double a) {return a * a; };
	const double c0 = -1460.0 * coef * (1.0 / pow2(24.0 * dx) + 1.0 / pow2(24.0 * dy) + 1.0 / pow2(24.0 * dz));
	const double c3x = coef / pow2(24.0 * dx);
	const double c2x = -54.0 * coef / pow2(24.0 * dx);
	const double c1x = 783.0 * coef / pow2(24.0 * dx);
	const double c3y = coef / pow2(24.0 * dy);
	const double c2y = -54.0 * coef / pow2(24.0 * dy);
	const double c1y = 783.0 * coef / pow2(24.0 * dy);
	const double c3z = coef / pow2(24.0 * dz);
	const double c2z = -54.0 * coef / pow2(24.0 * dz);
	const double c1z = 783.0 * coef / pow2(24.0 * dz);

	const int size_x = grid.SizeX() + 6;
	const int size_y = grid.SizeY() + 6;
	const int size_z = grid.SizeZ() + 6;

	for (int iz = 3; iz < size_z - 3; ++iz) {
		const int idz_m3 = -3 * size_x * size_y;
		const int idz_m2 = -2 * size_x * size_y;
		const int idz_m = -1 * size_x * size_y;
		const int idz_p = 1 * size_x * size_y;
		const int idz_p2 = 2 * size_x * size_y;
		const int idz_p3 = 3 * size_x * size_y;

		for (int iy = 3; iy < size_y - 3; ++iy) {
			const int idy_m3 = -3 * size_x;
			const int idy_m2 = -2 * size_x;
			const int idy_m = -1 * size_x;
			const int idy_p = 1 * size_x;
			const int idy_p2 = 2 * size_x;
			const int idy_p3 = 3 * size_x;
#pragma ivdep
			for (int ix = 3; ix < size_x - 3; ++ix) {
				const int idx_m3 = -3;
				const int idx_m2 = -2;
				const int idx_m = -1;
				const int idx_p = 1;
				const int idx_p2 = 2;
				const int idx_p3 = 3;
				const int64_t i = ix + size_x * (iy + (size_y * iz));
				const int64_t io = (ix - 3) + Nx * (iy - 3 + (Ny * (iz - 3)));
#if 1
				const auto d2psidx2 = c3x * ((p[i + idx_p3] + p[i + idx_m3]) - 54.0 * (p[i + idx_p2] + p[i + idx_m2]) + 783.0 * (p[i + idx_p] + p[i + idx_m]) - 1460.0 * p[i])
					+ c3y * ((p[i + idy_p3] + p[i + idy_m3]) - 54.0 * (p[i + idy_p2] + p[i + idy_m2]) + 783.0 * (p[i + idy_p] + p[i + idy_m]) - 1460.0 * p[i])
					+ c3z * ((p[i + idz_p3] + p[i + idz_m3]) - 54.0 * (p[i + idz_p2] + p[i + idz_m2]) + 783.0 * (p[i + idz_p] + p[i + idz_m]) - 1460.0 * p[i]);
				out[io] += d2psidx2;
#else
				auto d2psidx2 = c0 * p[i];
				d2psidx2 += c3x * (p[i + idx_p3] + p[i + idx_m3]) + c2x * (p[i + idx_p2] + p[i + idx_m2]) + c1x * (p[i + idx_p] + p[i + idx_m]);
				d2psidx2 += c3y * (p[i + idy_p3] + p[i + idy_m3]) + c2y * (p[i + idy_p2] + p[i + idy_m2]) + c1y * (p[i + idy_p] + p[i + idy_m]);
				d2psidx2 += c3z * (p[i + idz_p3] + p[i + idz_m3]) + c2z * (p[i + idz_p2] + p[i + idz_m2]) + c1z * (p[i + idz_p] + p[i + idz_m]);
				out[io] += d2psidx2;
#endif
			}
		}
	}
}


template <class T>
void Laplasian4th_halo(const GridRange& grid, T* out, const T* src, T* buf_with_halo, const double coef, const double dx, const double dy, const double dz) {
	const int size_x = grid.SizeX();
	const int size_y = grid.SizeY();
	const int size_z = grid.SizeZ();

	PasteHalo_f3(grid, buf_with_halo,src,
		[&](int ix, int iy, int iz) { return src[(ix + size_x * (iy + 1 + size_y * iz))]; },
		[&](int ix, int iy, int iz) { return src[(ix + size_x * (iy - 1 + size_y * iz))]; },
		[&](int ix, int iy, int iz) { return src[(ix + size_x * (iy + size_y * (iz + 1)))]; },
		[&](int ix, int iy, int iz) { return src[(ix + size_x * (iy + size_y * (iz - 1)))]; },
		[&](int ix, int iy, int iz) { return src[(ix + size_x * (iy + size_y * (iz + size_z)))]; },
		[&](int ix, int iy, int iz) { return src[(ix + size_x * (iy + size_y * (iz - size_z)))]; }	);


	Laplasian4th_halo_f(grid, out, buf_with_halo, coef, dx, dy, dz);


}

#else

template <class T, class FUNC_XL, class FUNC_XR, class FUNC_YL, class FUNC_YR, class FUNC_ZL, class FUNC_ZR>
void Laplasian4th_overspace_f(const GridRange& grid, T* out, const T* src, T* p, FUNC_XL&& func_xl, FUNC_XR&& func_xr, FUNC_YL&& func_yl, FUNC_YR&& func_yr, FUNC_ZL&& func_zl, FUNC_ZR&& func_zr, double coef, double dx, double dy, double dz) {
	const int Nx = grid.SizeX();
	const int Ny = grid.SizeY();
	const int Nz = grid.SizeZ();

	auto pow2 = [](T a) {return a * a; };
	const double c0 = -1460.0 * coef * (1.0 / pow2(24.0 * dx) + 1.0 / pow2(24.0 * dy) + 1.0 / pow2(24.0 * dz));
	const double c3x = coef / pow2(24.0 * dx);
	const double c2x = -54.0 * coef / pow2(24.0 * dx);
	const double c1x = 783.0 * coef / pow2(24.0 * dx);
	const double c3y = coef / pow2(24.0 * dy);
	const double c2y = -54.0 * coef / pow2(24.0 * dy);
	const double c1y = 783.0 * coef / pow2(24.0 * dy);
	const double c3z = coef / pow2(24.0 * dz);
	const double c2z = -54.0 * coef / pow2(24.0 * dz);
	const double c1z = 783.0 * coef / pow2(24.0 * dz);


	{
#pragma ivdep
		for (int iz = 0; iz < Nz; ++iz) {
			for (int iy = 0; iy < Ny; ++iy) {
				for (int ix = 0; ix < Nx; ++ix) {
					p[(ix + 3) + (Nx + 6) * (iy + 3 + (Ny + 6) * (iz + 3))] = src[ix + Nx * (iy + Ny * iz)];
				}
			}
		}
#pragma ivdep
		for (int iz = 0; iz < Nz; ++iz) {
			for (int iy = 0; iy < Ny; ++iy) {
				p[0 + (Nx + 6) * (iy + 3 + (Ny + 6) * (iz + 3))] = func_xl(-3, iy, iz);
				p[1 + (Nx + 6) * (iy + 3 + (Ny + 6) * (iz + 3))] = func_xl(-2, iy, iz);
				p[2 + (Nx + 6) * (iy + 3 + (Ny + 6) * (iz + 3))] = func_xl(-1, iy, iz);
				p[Nx + 5 + (Nx + 6) * (iy + 3 + (Ny + 6) * (iz + 3))] = func_xr(Nx + 2, iy, iz);
				p[Nx + 4 + (Nx + 6) * (iy + 3 + (Ny + 6) * (iz + 3))] = func_xr(Nx + 1, iy, iz);
				p[Nx + 3 + (Nx + 6) * (iy + 3 + (Ny + 6) * (iz + 3))] = func_xr(Nx + 0, iy, iz);
			}
		}
#pragma ivdep
		for (int iz = 0; iz < Nz; ++iz) {
			for (int ix = 0; ix < Nx; ++ix) {
				p[(ix + 3) + (Nx + 6) * (0 + (Ny + 6) * (iz + 3))] = func_yl(ix, -3, iz);
				p[(ix + 3) + (Nx + 6) * (1 + (Ny + 6) * (iz + 3))] = func_yl(ix, -2, iz);
				p[(ix + 3) + (Nx + 6) * (2 + (Ny + 6) * (iz + 3))] = func_yl(ix, -1, iz);
				p[(ix + 3) + (Nx + 6) * (Ny + 5 + (Ny + 6) * (iz + 3))] = func_yr(ix, Ny + 2, iz);
				p[(ix + 3) + (Nx + 6) * (Ny + 4 + (Ny + 6) * (iz + 3))] = func_yr(ix, Ny + 1, iz);
				p[(ix + 3) + (Nx + 6) * (Ny + 3 + (Ny + 6) * (iz + 3))] = func_yr(ix, Ny + 0, iz);
			}
		}
#pragma ivdep
		for (int iy = 0; iy < Ny; ++iy) {
			for (int ix = 0; ix < Nx; ++ix) {
				p[(ix + 3) + (Nx + 6) * (iy + 3 + (Ny + 6) * (0))] = func_zl(ix, iy, -3);
				p[(ix + 3) + (Nx + 6) * (iy + 3 + (Ny + 6) * (1))] = func_zl(ix, iy, -2);
				p[(ix + 3) + (Nx + 6) * (iy + 3 + (Ny + 6) * (2))] = func_zl(ix, iy, -1);
			}
		}
		for (int iy = 0; iy < Ny; ++iy) {
			for (int ix = 0; ix < Nx; ++ix) {
				p[(ix + 3) + (Nx + 6) * (iy + 3 + (Ny + 6) * (Nz + 3))] = func_zr(ix, iy, Nz + 0);
				p[(ix + 3) + (Nx + 6) * (iy + 3 + (Ny + 6) * (Nz + 4))] = func_zr(ix, iy, Nz + 1);
				p[(ix + 3) + (Nx + 6) * (iy + 3 + (Ny + 6) * (Nz + 5))] = func_zr(ix, iy, Nz + 2);
			}
		}
	}


	const int size_x = grid.SizeX() + 6;
	const int size_y = grid.SizeY() + 6;
	const int size_z = grid.SizeZ() + 6;

	for (int iz = 3; iz < size_z - 3; ++iz) {
		const int idz_m3 = -3 * size_x * size_y;
		const int idz_m2 = -2 * size_x * size_y;
		const int idz_m = -1 * size_x * size_y;
		const int idz_p = 1 * size_x * size_y;
		const int idz_p2 = 2 * size_x * size_y;
		const int idz_p3 = 3 * size_x * size_y;

		for (int iy = 3; iy < size_y - 3; ++iy) {
			const int idy_m3 = -3 * size_x;
			const int idy_m2 = -2 * size_x;
			const int idy_m = -1 * size_x;
			const int idy_p = 1 * size_x;
			const int idy_p2 = 2 * size_x;
			const int idy_p3 = 3 * size_x;
#pragma ivdep
			for (int ix = 3; ix < size_x - 3; ++ix) {
				const int idx_m3 = -3;
				const int idx_m2 = -2;
				const int idx_m = -1;
				const int idx_p = 1;
				const int idx_p2 = 2;
				const int idx_p3 = 3;
				const int64_t i = ix + size_x * (iy + (size_y * iz));
				const int64_t io = (ix - 3) + Nx * (iy - 3 + (Ny * (iz - 3)));
#if 1
				const auto d2psidx2 = c3x*((p[i + idx_p3] + p[i + idx_m3]) - 54.0 * (p[i + idx_p2] + p[i + idx_m2]) + 783.0 * (p[i + idx_p] + p[i + idx_m]) - 1460.0 * p[i])
					+ c3y * ((p[i + idy_p3] + p[i + idy_m3]) - 54.0 * (p[i + idy_p2] + p[i + idy_m2]) + 783.0 * (p[i + idy_p] + p[i + idy_m]) - 1460.0 * p[i])
					+ c3z * ((p[i + idz_p3] + p[i + idz_m3]) - 54.0 * (p[i + idz_p2] + p[i + idz_m2]) + 783.0 * (p[i + idz_p] + p[i + idz_m]) - 1460.0 * p[i]);
				out[io] += d2psidx2;
#else
				auto d2psidx2 = c0 * p[i];
				d2psidx2 += c3x * (p[i + idx_p3] + p[i + idx_m3]) + c2x * (p[i + idx_p2] + p[i + idx_m2]) + c1x * (p[i + idx_p] + p[i + idx_m]);
				d2psidx2 += c3y * (p[i + idy_p3] + p[i + idy_m3]) + c2y * (p[i + idy_p2] + p[i + idy_m2]) + c1y * (p[i + idy_p] + p[i + idy_m]);
				d2psidx2 += c3z * (p[i + idz_p3] + p[i + idz_m3]) + c2z * (p[i + idz_p2] + p[i + idz_m2]) + c1z * (p[i + idz_p] + p[i + idz_m]);
				out[io] += d2psidx2;
#endif
			}
		}
	}
}


template <class T>
void Laplasian4th_overspace2(const GridRange& grid, T* out, const T* src, T* p, const double coef, const double dx, const double dy, const double dz) {
	const int size_x = grid.SizeX();
	const int size_y = grid.SizeY();
	const int size_z = grid.SizeZ();

	Laplasian4th_overspace_f(grid, out, src, p,
		[&](int ix, int iy, int iz) { return src[(ix + size_x * (iy + 1 + size_y * iz))]; },
		[&](int ix, int iy, int iz) { return src[(ix + size_x * (iy - 1 + size_y * iz))]; },
		[&](int ix, int iy, int iz) { return src[(ix + size_x * (iy + size_y * (iz + 1)))]; },
		[&](int ix, int iy, int iz) { return src[(ix + size_x * (iy + size_y * (iz - 1)))]; },
		[&](int ix, int iy, int iz) { return src[(ix + size_x * (iy + size_y * (iz + size_z)))]; },
		[&](int ix, int iy, int iz) { return src[(ix + size_x * (iy + size_y * (iz - size_z)))]; },
		coef, dx, dy, dz);

	
}
#endif

template <int MY, int MZ, class T>
void Laplasian4th_7_X(int size_x, int size_y, int iy, int iz, T* out, const T* p, const T* m_y_left, const T* m_y_right, const T* m_z_left, const T* m_z_right, double coef, double dx, double dy, double dz) {
	const double c0 = -1460.0 / (dx * dx) - 1460.0 / (dy * dy) - 1460.0 / (dz * dz);
	const int idy = size_x;
	const int idz = size_x * size_y;

	auto d2dy2 = [&](int64_t i, int ix) {
		if constexpr (MY == 0) {
			return ((p[i - 3*idy] + p[i + 3* idy]) - 54.0 * (p[i - 2* idy] + p[i + 2*idy]) + 783.0 * (p[i - idy] + p[i + idy]) ) / (dy * dy);
		} else if constexpr (MY == 1) {
			return ((p[i - 3 * idy] + m_y_right[ix]) - 54.0 * (p[i - 2 * idy] + p[i + 2 * idy]) + 783.0 * (p[i - idy] + p[i + idy])) / (dy * dy);
		} else if constexpr (MY == 2) {
			return ((p[i - 3 * idy] + m_y_right[ix + size_x]) - 54.0 * (p[i - 2 * idy] + m_y_right[ix]) + 783.0 * (p[i - idy] + p[i + idy])) / (dy * dy);
		} else if constexpr (MY == 3) {
			return ((p[i - 3 * idy] + m_y_right[ix + 2*size_x]) - 54.0 * (p[i - 2 * idy] + m_y_right[ix + size_x]) + 783.0 * (p[i - idy] + m_y_right[ix])) / (dy * dy);
		} else if constexpr (MY == -1) {
			return ((m_y_left[ix + 2 * size_x] + p[i + 3 * idy]) - 54.0 * (p[i - 2 * idy] + p[i + 2 * idy]) + 783.0 * (p[i - idy] + p[i + idy])) / (dy * dy);
		} else if constexpr (MY == -2) {
			return ((m_y_left[ix + 1 * size_x] + p[i + 3 * idy]) - 54.0 * (m_y_left[ix + 2 * size_x] + p[i + 2 * idy]) + 783.0 * (p[i - idy] + p[i + idy])) / (dy * dy);
		} else if constexpr (MY == -3) {
			return ((m_y_left[ix] + p[i + 3 * idy]) - 54.0 * (m_y_left[ix + 1 * size_x] + p[i + 2 * idy]) + 783.0 * (m_y_left[ix + 2 * size_x] + p[i + idy])) / (dy * dy);
		}
	};

	auto d2dz2 = [&](int64_t i, int ix) {
		if constexpr (MZ == 0) {
			return ((p[i - 3*idz] + p[i + 3*idz]) - 54.0 * (p[i - 2*idz] + p[i + 2*idz]) + 783.0 * (p[i - idz] + p[i + idz])) / (dz * dz);
		} else if constexpr (MZ == 1) {
			return ((p[i - 3*idz] + m_z_right[ix]) - 54.0 * (p[i - 2*idz] + p[i + 2*idz]) + 783.0 * (p[i - idz] + p[i + idz])) / (dz * dz);
		} else if constexpr (MZ == 2) {
			return ((p[i - 3*idz] + m_z_right[ix + size_x]) - 54.0 * (p[i - 2*idz] + m_z_right[ix]) + 783.0 * (p[i - idz] + p[i + idz])) / (dz * dz);
		} else if constexpr (MZ == 3) {
			return ((p[i - 3*idz] + m_z_right[ix + 2*size_x]) - 54.0 * (p[i - 2*idz] + m_z_right[ix + size_x]) + 783.0 * (p[i - idz] + m_z_right[ix])) / (dz * dz);
		} else if constexpr (MZ == -1) {
			return ((m_z_left[ix + 2*size_x] + p[i + 3*idz]) - 54.0 * (p[i - 2*idz] + p[i + 2*idz]) + 783.0 * (p[i - idz] + p[i + idz])) / (dz * dz);
		} else if constexpr (MZ == -2) {
			return ((m_z_left[ix + 1 * size_x] + p[i + 3*idz]) - 54.0 * (m_z_left[ix + 2 * size_x] + p[i + 2*idz]) + 783.0 * (p[i - idz] + p[i + idz])) / (dz * dz);
		} else if constexpr (MZ == -3) {
			return ((m_z_left[ix] + p[i + 3*idz]) - 54.0 * (m_z_left[ix + size_x] + p[i + 2*idz]) + 783.0 * (m_z_left[ix + 2 * size_x] + p[i + idz])) / (dz * dz);
		}
	};

	for (int ix = 3; ix < size_x - 3; ++ix) {
		const int64_t i = ix + size_x * (iy + (size_y * iz));

		auto d2psidx2 = c0 * p[i];
		d2psidx2 += ((p[i - 3] + p[i + 3]) - 54.0 * (p[i - 2] + p[i + 2]) + 783.0 * (p[i - 1] + p[i + 1])) / (dx * dx);
		d2psidx2 += d2dy2(i,ix) + d2dz2(i, ix);
		out[i] += coef * d2psidx2 / (24.0 * 24.0);

	}
	{
		int ix = 2;
		const int64_t i = ix + size_x * (iy + (size_y * iz));

		auto d2psidx2 = c0 * p[i];
		d2psidx2 += ((p[i - 3 + size_x] + p[i + 3]) - 54.0 * (p[i - 2] + p[i + 2]) + 783.0 * (p[i - 1] + p[i + 1])) / (dx * dx);
		d2psidx2 += d2dy2(i, ix) + d2dz2(i, ix);
		out[i] += coef * d2psidx2 / (24.0 * 24.0);

	}
	{
		int ix = 1;
		const int64_t i = ix + size_x * (iy + (size_y * iz));

		auto d2psidx2 = c0 * p[i];
		d2psidx2 += ((p[i - 3 + size_x] + p[i + 3]) - 54.0 * (p[i - 2 + size_x] + p[i + 2]) + 783.0 * (p[i - 1] + p[i + 1])) / (dx * dx);
		d2psidx2 += d2dy2(i, ix) + d2dz2(i, ix);
		out[i] += coef * d2psidx2 / (24.0 * 24.0);

	}
	{
		int ix = 0;
		const int64_t i = ix + size_x * (iy + (size_y * iz));

		auto d2psidx2 = c0 * p[i];
		d2psidx2 += ((p[i - 3 + size_x] + p[i + 3]) - 54.0 * (p[i - 2 + size_x] + p[i + 2]) + 783.0 * (p[i - 1 + size_x] + p[i + 1])) / (dx * dx);
		d2psidx2 += d2dy2(i, ix) + d2dz2(i, ix);
		out[i] += coef * d2psidx2 / (24.0 * 24.0);

	}
	{
		int ix = size_x - 3;
		const int64_t i = ix + size_x * (iy + (size_y * iz));
		
		auto d2psidx2 = c0 * p[i];
		d2psidx2 += ((p[i - 3] + p[i + 3 - size_x]) - 54.0 * (p[i - 2] + p[i + 2]) + 783.0 * (p[i - 1] + p[i + 1])) / (dx * dx);
		d2psidx2 += d2dy2(i, ix) + d2dz2(i, ix);
		out[i] += coef * d2psidx2 / (24.0 * 24.0);

	}
	{
		int ix = size_x - 2;
		const int64_t i = ix + size_x * (iy + (size_y * iz));

		auto d2psidx2 = c0 * p[i];
		d2psidx2 += ((p[i - 3] + p[i + 3 - size_x]) - 54.0 * (p[i - 2] + p[i + 2 - size_x]) + 783.0 * (p[i - 1] + p[i + 1])) / (dx * dx);
		d2psidx2 += d2dy2(i, ix) + d2dz2(i, ix);
		out[i] += coef * d2psidx2 / (24.0 * 24.0);

	}
	{
		int ix = size_x - 1;
		const int64_t i = ix + size_x * (iy + (size_y * iz));

		auto d2psidx2 = c0 * p[i];
		d2psidx2 += ((p[i - 3] + p[i + 3 - size_x]) - 54.0 * (p[i - 2] + p[i + 2 - size_x]) + 783.0 * (p[i - 1] + p[i + 1 - size_x])) / (dx * dx);
		d2psidx2 += d2dy2(i, ix) + d2dz2(i, ix);
		out[i] += coef * d2psidx2 / (24.0 * 24.0);

	}

}


template <int MZ, class T>
void Laplasian4th_8_Y(int size_x, int size_y, int size_z, int iz, T* out, const T* p, double coef, double dx, double dy, double dz) {

	for (int iy = 3; iy < size_y - 3; ++iy) {
		Laplasian4th_7_X<0, MZ>(size_x, size_y, iy, iz, out, p,
			&p[size_x * (iy - 3 + size_y)], &p[size_x * (iy + 1 - size_y)],
			&p[size_x * size_y * (iz - 3 + size_z)], &p[size_x * size_y * (iz + 1 - size_z)], coef, dx, dy, dz);
	}
	{
		const int iy = 0;
		Laplasian4th_7_X<-3, MZ>(size_x, size_y, iy, iz, out, p,
			&p[size_x * (iy - 3 + size_y)], &p[size_x * (iy + 1 - size_y)],
			&p[size_x * size_y * (iz - 3 + size_z)], &p[size_x * size_y * (iz + 1 - size_z)], coef, dx, dy, dz);
	}
	{
		const int iy = 1;
		Laplasian4th_7_X<-2, MZ>(size_x, size_y, iy, iz, out, p,
			&p[size_x * (iy - 3 + size_y)], &p[size_x * (iy + 1 - size_y)],
			&p[size_x * size_y * (iz - 3 + size_z)], &p[size_x * size_y * (iz + 1 - size_z)], coef, dx, dy, dz);
	}
	{
		const int iy = 2;
		Laplasian4th_7_X<-1, MZ>(size_x, size_y, iy, iz, out, p,
			&p[size_x * (iy - 3 + size_y)], &p[size_x * (iy + 1 - size_y)],
			&p[size_x * size_y * (iz - 3 + size_z)], &p[size_x * size_y * (iz + 1 - size_z)], coef, dx, dy, dz);
	}

	{
		const int iy = size_y - 3;
		Laplasian4th_7_X<1, MZ>(size_x, size_y, iy, iz, out, p,
			&p[size_x * (iy - 3 + size_y)], &p[size_x * (iy + 1 - size_y)],
			&p[size_x * size_y * (iz - 3 + size_z)], &p[size_x * size_y * (iz + 1 - size_z)], coef, dx, dy, dz);
	}
	{
		const int iy = size_y - 2;
		Laplasian4th_7_X<2, MZ>(size_x, size_y, iy, iz, out, p,
			&p[size_x * (iy - 3 + size_y)], &p[size_x * (iy + 1 - size_y)],
			&p[size_x * size_y * (iz - 3 + size_z)], &p[size_x * size_y * (iz + 1 - size_z)], coef, dx, dy, dz);
	}
	{
		const int iy = size_y - 1;
		Laplasian4th_7_X<3, MZ>(size_x, size_y, iy, iz, out, p,
			&p[size_x * (iy - 3 + size_y)], &p[size_x * (iy + 1 - size_y)],
			&p[size_x * size_y * (iz - 3 + size_z)], &p[size_x * size_y * (iz + 1 - size_z)], coef, dx, dy, dz);
	}

}



template <class T>
void Laplasian4th_8(const GridRange& grid, T* out, const T* p, const double coef, const double dx, const double dy, const double dz) {
	const int size_x = grid.SizeX();
	const int size_y = grid.SizeY();
	const int size_z = grid.SizeZ();
	for (int iz = 3; iz < size_z - 3; ++iz) {
		Laplasian4th_8_Y<0>(size_x, size_y, size_z, iz, out, p, coef, dx, dy, dz);
	}
	{
		const int iz = 0;
		Laplasian4th_8_Y<-3>(size_x, size_y, size_z, iz, out, p, coef, dx, dy, dz);

	}
	{
		const int iz = 1;
		Laplasian4th_8_Y<-2>(size_x, size_y, size_z, iz, out, p, coef, dx, dy, dz);

	}
	{
		const int iz = 2;
		Laplasian4th_8_Y<-1>(size_x, size_y, size_z, iz, out, p, coef, dx, dy, dz);

	}
	{
		const int iz = size_z - 3;
		Laplasian4th_8_Y<1>(size_x, size_y, size_z, iz, out, p, coef, dx, dy, dz);

	}
	{
		const int iz = size_z - 2;
		Laplasian4th_8_Y<2>(size_x, size_y, size_z, iz, out, p, coef, dx, dy, dz);

	}
	{
		const int iz = size_z - 1;
		Laplasian4th_8_Y<3>(size_x, size_y, size_z, iz, out, p, coef, dx, dy, dz);

	}
}


/*
* wrapper of the fastest algorithm of 4th order difference equation
* 
*/
template <class T>
void Laplasian4th(const GridRange& grid, T* out, const T* p, double coef, double dx, double dy, double dz) {
	
	const int margin_width = 3;
	const int over_size = (grid.SizeX() + 2 * margin_width) * (grid.SizeY() + 2 * margin_width) * (grid.SizeZ() + 2 * margin_width);
	static std::vector<T> over_data(over_size);//margin data//

	Laplasian4th_halo(grid, out, p, over_data.data(), coef, dx, dy, dz);
}



template <class T, class FUNC_XL, class FUNC_XR, class FUNC_YL, class FUNC_YR, class FUNC_ZL, class FUNC_ZR>
void Laplasian4th_bundle_X(int64_t Ns, int size_x, int size_y, int size_z, T* out, const T* p, FUNC_XL&& func_xl, FUNC_XR&& func_xr, FUNC_YL&& func_yl, FUNC_YR&& func_yr, FUNC_ZL&& func_zl, FUNC_ZR&& func_zr, double coef, double dx, double dy, double dz) {
	const double c0 = -1460.0 / (dx * dx) - 1460.0 / (dy * dy) - 1460.0 / (dz * dz);
	
	auto GetHeadXL = [&](int ix, int iy, int iz) {
		if (ix < 0) return func_xl(ix, iy, iz);
		return p + Ns * (ix + size_x * (iy + (size_y * iz)));
		};
	auto GetHeadXR = [&](int ix, int iy, int iz) {
		if (ix >= size_x) return func_xr(ix, iy, iz);
		return p + Ns * (ix + size_x * (iy + (size_y * iz)));
		};
	auto GetHeadYL = [&](int ix, int iy, int iz) {
		if (iy < 0) return func_yl(ix, iy, iz);
		return p + Ns * (ix + size_x * (iy + (size_y * iz)));
		};
	auto GetHeadYR = [&](int ix, int iy, int iz) {
		if (iy >= size_y) return func_yr(ix, iy, iz);
		return p + Ns * (ix + size_x * (iy + (size_y * iz)));
		};
	auto GetHeadZL = [&](int ix, int iy, int iz) {
		if (iz < 0) return func_zl(ix, iy, iz);
		return p + Ns * (ix + size_x * (iy + (size_y * iz)));
		};
	auto GetHeadZR = [&](int ix, int iy, int iz) {
		if (iz >= size_z) return func_zr(ix, iy, iz);
		return p + Ns * (ix + size_x * (iy + (size_y * iz)));
		};

	
	for (int iz = 0; iz < size_z; ++iz) {
		for (int iy = 0; iy < size_y; ++iy) {
			for (int ix = 0; ix < size_x; ++ix) {
				const int64_t i = Ns * (ix + size_x * (iy + (size_y * iz)));
				const T* p0 = p + i;
				const T* xm1 = GetHeadXL(ix - 1, iy, iz);
				const T* xm2 = GetHeadXL(ix - 2, iy, iz);
				const T* xm3 = GetHeadXL(ix - 3, iy, iz);
				const T* xp1 = GetHeadXR(ix + 1, iy, iz);
				const T* xp2 = GetHeadXR(ix + 2, iy, iz);
				const T* xp3 = GetHeadXR(ix + 3, iy, iz);
				const T* zm1 = GetHeadZL(ix, iy, iz - 1);
				const T* zm2 = GetHeadZL(ix, iy, iz - 2);
				const T* zm3 = GetHeadZL(ix, iy, iz - 3);
				const T* zp1 = GetHeadZR(ix, iy, iz + 1);
				const T* zp2 = GetHeadZR(ix, iy, iz + 2);
				const T* zp3 = GetHeadZR(ix, iy, iz + 3);
				const T* ym1 = GetHeadYL(ix, iy - 1, iz);
				const T* ym2 = GetHeadYL(ix, iy - 2, iz);
				const T* ym3 = GetHeadYL(ix, iy - 3, iz);
				const T* yp1 = GetHeadYR(ix, iy + 1, iz);
				const T* yp2 = GetHeadYR(ix, iy + 2, iz);
				const T* yp3 = GetHeadYR(ix, iy + 3, iz);
#pragma ivdep
				for (int n = 0; n < Ns; ++n) {
					auto d2psidx2 = c0 * p0[n];
					d2psidx2 += ((xm3[n] + xp3[n]) - 54.0 * (xm2[n] + xp2[n]) + 783.0 * (xm1[n] + xp1[n])) / (dx * dx);
					d2psidx2 += ((ym3[n] + yp3[n]) - 54.0 * (ym2[n] + yp2[n]) + 783.0 * (ym1[n] + yp1[n])) / (dy * dy);
					d2psidx2 += ((zm3[n] + zp3[n]) - 54.0 * (zm2[n] + zp2[n]) + 783.0 * (zm1[n] + zp1[n])) / (dz * dz);
					out[n + i] += coef * d2psidx2 / (24.0 * 24.0);
				}
			}
		}
	}
	
}


template <class T, class FUNC_XL, class FUNC_XR, class FUNC_YL, class FUNC_YR, class FUNC_ZL, class FUNC_ZR>
void Laplasian4th_bundle_X1c(int64_t Ns, int size_x, int size_y, int size_z, T* out, const T* p, FUNC_XL&& func_xl, FUNC_XR&& func_xr, FUNC_YL&& func_yl, FUNC_YR&& func_yr, FUNC_ZL&& func_zl, FUNC_ZR&& func_zr, double coef, double dx, double dy, double dz) {
	auto pow2 = [](T a) {return a * a; };
	const double c0 = -1460.0 * coef * (1.0 / (dx * dx) + 1.0 / (dy * dy) + 1.0 / (dz * dz)) / (24.0 * 24.0);
	const double cx1 = coef / pow2(24.0 * dx);
	const double cx2 = -54.0 * coef / pow2(24.0 * dx);
	const double cx3 = 783.0 * coef / pow2(24.0 * dx);
	const double cy1 = coef / pow2(24.0 * dy);
	const double cy2 = -54.0 * coef / pow2(24.0 * dy);
	const double cy3 = 783.0 * coef / pow2(24.0 * dy);
	const double cz1 = coef / pow2(24.0 * dz);
	const double cz2 = -54.0 * coef / pow2(24.0 * dz);
	const double cz3 = 783.0 * coef / pow2(24.0 * dz);

	auto GetHeadXL = [&](int ix, int iy, int iz) {
		if (ix < 0) return func_xl(ix, iy, iz);
		return p + Ns * (ix + size_x * (iy + (size_y * iz)));
		};
	auto GetHeadXR = [&](int ix, int iy, int iz) {
		if (ix >= size_x) return func_xr(ix, iy, iz);
		return p + Ns * (ix + size_x * (iy + (size_y * iz)));
		};
	auto GetHeadYL = [&](int ix, int iy, int iz) {
		if (iy < 0) return func_yl(ix, iy, iz);
		return p + Ns * (ix + size_x * (iy + (size_y * iz)));
		};
	auto GetHeadYR = [&](int ix, int iy, int iz) {
		if (iy >= size_y) return func_yr(ix, iy, iz);
		return p + Ns * (ix + size_x * (iy + (size_y * iz)));
		};
	auto GetHeadZL = [&](int ix, int iy, int iz) {
		if (iz < 0) return func_zl(ix, iy, iz);
		return p + Ns * (ix + size_x * (iy + (size_y * iz)));
		};
	auto GetHeadZR = [&](int ix, int iy, int iz) {
		if (iz >= size_z) return func_zr(ix, iy, iz);
		return p + Ns * (ix + size_x * (iy + (size_y * iz)));
		};


	for (int iz = 0; iz < size_z; ++iz) {
		for (int iy = 0; iy < size_y; ++iy) {
			for (int ix = 0; ix < size_x; ++ix) {
				const int64_t i = Ns * (ix + size_x * (iy + (size_y * iz)));
				const T* p0 = p + i;
				const T* xm1 = GetHeadXL(ix - 1, iy, iz);
				const T* xm2 = GetHeadXL(ix - 2, iy, iz);
				const T* xm3 = GetHeadXL(ix - 3, iy, iz);
				const T* xp1 = GetHeadXR(ix + 1, iy, iz);
				const T* xp2 = GetHeadXR(ix + 2, iy, iz);
				const T* xp3 = GetHeadXR(ix + 3, iy, iz);
				const T* zm1 = GetHeadZL(ix, iy, iz - 1);
				const T* zm2 = GetHeadZL(ix, iy, iz - 2);
				const T* zm3 = GetHeadZL(ix, iy, iz - 3);
				const T* zp1 = GetHeadZR(ix, iy, iz + 1);
				const T* zp2 = GetHeadZR(ix, iy, iz + 2);
				const T* zp3 = GetHeadZR(ix, iy, iz + 3);
				const T* ym1 = GetHeadYL(ix, iy - 1, iz);
				const T* ym2 = GetHeadYL(ix, iy - 2, iz);
				const T* ym3 = GetHeadYL(ix, iy - 3, iz);
				const T* yp1 = GetHeadYR(ix, iy + 1, iz);
				const T* yp2 = GetHeadYR(ix, iy + 2, iz);
				const T* yp3 = GetHeadYR(ix, iy + 3, iz);
#pragma ivdep
				for (int n = 0; n < Ns; ++n) {
					const auto d2psidx2 = c0 * p0[n]
						+ cx1 * (xm3[n] + xp3[n]) + cx2 * (xm2[n] + xp2[n]) + cx3 * (xm1[n] + xp1[n])
						+ cy1 * (ym3[n] + yp3[n]) + cy2 * (ym2[n] + yp2[n]) + cy3 * (ym1[n] + yp1[n])
						+ cz1 * (zm3[n] + zp3[n]) + cz2 * (zm2[n] + zp2[n]) + cz3 * (zm1[n] + zp1[n]);
					out[n + i] += d2psidx2;
				}
			}
		}
	}


}



template <class T, class FUNC_XL, class FUNC_XR, class FUNC_YL, class FUNC_YR, class FUNC_ZL, class FUNC_ZR>
void Laplasian4th_bundle_X3(int64_t Ns, int size_x, int size_y, int size_z, T* out, const T* p, FUNC_XL&& func_xl, FUNC_XR&& func_xr, FUNC_YL&& func_yl, FUNC_YR&& func_yr, FUNC_ZL&& func_zl, FUNC_ZR&& func_zr, double coef, double dx, double dy, double dz) {
	const double c0 = -1460.0 / (dx * dx) - 1460.0 / (dy * dy) - 1460.0 / (dz * dz);


	auto GetHeadXL = [&](int ix, int iy, int iz) {
		if (ix < 0) return func_xl(ix, iy, iz);
		return p + Ns * (ix + size_x * (iy + (size_y * iz)));
		};
	auto GetHeadXR = [&](int ix, int iy, int iz) {
		if (ix >= size_x) return func_xr(ix, iy, iz);
		return p + Ns * (ix + size_x * (iy + (size_y * iz)));
		};
	auto GetHeadYL = [&](int ix, int iy, int iz) {
		if (iy < 0) return func_yl(ix, iy, iz);
		return p + Ns * (ix + size_x * (iy + (size_y * iz)));
		};
	auto GetHeadYR = [&](int ix, int iy, int iz) {
		if (iy >= size_y) return func_yr(ix, iy, iz);
		return p + Ns * (ix + size_x * (iy + (size_y * iz)));
		};
	auto GetHeadZL = [&](int ix, int iy, int iz) {
		if (iz < 0) return func_zl(ix, iy, iz);
		return p + Ns * (ix + size_x * (iy + (size_y * iz)));
		};
	auto GetHeadZR = [&](int ix, int iy, int iz) {
		if (iz >= size_z) return func_zr(ix, iy, iz);
		return p + Ns * (ix + size_x * (iy + (size_y * iz)));
		};


	for (int iz = 0; iz < size_z; ++iz) {
		for (int iy = 0; iy < size_y; ++iy) {
			for (int ix = 0; ix < size_x; ++ix) {
				const int64_t i = Ns * (ix + size_x * (iy + (size_y * iz)));
				const T* p0 = p + i;
				const T* xm1 = GetHeadXL(ix - 1, iy, iz);
				const T* xm2 = GetHeadXL(ix - 2, iy, iz);
				const T* xm3 = GetHeadXL(ix - 3, iy, iz);
				const T* xp1 = GetHeadXR(ix + 1, iy, iz);
				const T* xp2 = GetHeadXR(ix + 2, iy, iz);
				const T* xp3 = GetHeadXR(ix + 3, iy, iz);
				const T* zm1 = GetHeadZL(ix, iy, iz - 1);
				const T* zm2 = GetHeadZL(ix, iy, iz - 2);
				const T* zm3 = GetHeadZL(ix, iy, iz - 3);
				const T* zp1 = GetHeadZR(ix, iy, iz + 1);
				const T* zp2 = GetHeadZR(ix, iy, iz + 2);
				const T* zp3 = GetHeadZR(ix, iy, iz + 3);
				const T* ym1 = GetHeadYL(ix, iy - 1, iz);
				const T* ym2 = GetHeadYL(ix, iy - 2, iz);
				const T* ym3 = GetHeadYL(ix, iy - 3, iz);
				const T* yp1 = GetHeadYR(ix, iy + 1, iz);
				const T* yp2 = GetHeadYR(ix, iy + 2, iz);
				const T* yp3 = GetHeadYR(ix, iy + 3, iz);
#pragma ivdep
				for (int n = 0; n < Ns; ++n) {
					const auto d2psidx2 = c0 * p0[n]
						+ ((xm3[n] + xp3[n]) - 54.0 * (xm2[n] + xp2[n]) + 783.0 * (xm1[n] + xp1[n])) / (dx * dx);
					//d2psidx2 += ((ym3[n] + yp3[n]) - 54.0 * (ym2[n] + yp2[n]) + 783.0 * (ym1[n] + yp1[n])) / (dy * dy);
					//d2psidx2 += ((zm3[n] + zp3[n]) - 54.0 * (zm2[n] + zp2[n]) + 783.0 * (zm1[n] + zp1[n])) / (dz * dz);
					out[n + i] += coef * d2psidx2 / (24.0 * 24.0);
				}
			}
		}
	}
	for (int iz = 0; iz < size_z; ++iz) {
		for (int iy = 0; iy < size_y; ++iy) {
			for (int ix = 0; ix < size_x; ++ix) {
				const int64_t i = Ns * (ix + size_x * (iy + (size_y * iz)));
				const T* p0 = p + i;
				const T* xm1 = GetHeadXL(ix - 1, iy, iz);
				const T* xm2 = GetHeadXL(ix - 2, iy, iz);
				const T* xm3 = GetHeadXL(ix - 3, iy, iz);
				const T* xp1 = GetHeadXR(ix + 1, iy, iz);
				const T* xp2 = GetHeadXR(ix + 2, iy, iz);
				const T* xp3 = GetHeadXR(ix + 3, iy, iz);
				const T* zm1 = GetHeadZL(ix, iy, iz - 1);
				const T* zm2 = GetHeadZL(ix, iy, iz - 2);
				const T* zm3 = GetHeadZL(ix, iy, iz - 3);
				const T* zp1 = GetHeadZR(ix, iy, iz + 1);
				const T* zp2 = GetHeadZR(ix, iy, iz + 2);
				const T* zp3 = GetHeadZR(ix, iy, iz + 3);
				const T* ym1 = GetHeadYL(ix, iy - 1, iz);
				const T* ym2 = GetHeadYL(ix, iy - 2, iz);
				const T* ym3 = GetHeadYL(ix, iy - 3, iz);
				const T* yp1 = GetHeadYR(ix, iy + 1, iz);
				const T* yp2 = GetHeadYR(ix, iy + 2, iz);
				const T* yp3 = GetHeadYR(ix, iy + 3, iz);
#pragma ivdep
				for (int n = 0; n < Ns; ++n) {
					//const auto d2psidx2 = c0 * p0[n];
//					d2psidx2 += ((xm3[n] + xp3[n]) - 54.0 * (xm2[n] + xp2[n]) + 783.0 * (xm1[n] + xp1[n])) / (dx * dx);
//					d2psidx2 += ((ym3[n] + yp3[n]) - 54.0 * (ym2[n] + yp2[n]) + 783.0 * (ym1[n] + yp1[n])) / (dy * dy);
//					d2psidx2 += ((zm3[n] + zp3[n]) - 54.0 * (zm2[n] + zp2[n]) + 783.0 * (zm1[n] + zp1[n])) / (dz * dz);
					//out[n + i] += coef * d2psidx2 / (24.0 * 24.0);
				}
			}
		}
	}

	for (int iz = 0; iz < size_z; ++iz) {
		for (int iy = 0; iy < size_y; ++iy) {
			for (int ix = 0; ix < size_x; ++ix) {
				const int64_t i = Ns * (ix + size_x * (iy + (size_y * iz)));
				const T* p0 = p + i;
				const T* xm1 = GetHeadXL(ix - 1, iy, iz);
				const T* xm2 = GetHeadXL(ix - 2, iy, iz);
				const T* xm3 = GetHeadXL(ix - 3, iy, iz);
				const T* xp1 = GetHeadXR(ix + 1, iy, iz);
				const T* xp2 = GetHeadXR(ix + 2, iy, iz);
				const T* xp3 = GetHeadXR(ix + 3, iy, iz);
				const T* zm1 = GetHeadZL(ix, iy, iz - 1);
				const T* zm2 = GetHeadZL(ix, iy, iz - 2);
				const T* zm3 = GetHeadZL(ix, iy, iz - 3);
				const T* zp1 = GetHeadZR(ix, iy, iz + 1);
				const T* zp2 = GetHeadZR(ix, iy, iz + 2);
				const T* zp3 = GetHeadZR(ix, iy, iz + 3);
				const T* ym1 = GetHeadYL(ix, iy - 1, iz);
				const T* ym2 = GetHeadYL(ix, iy - 2, iz);
				const T* ym3 = GetHeadYL(ix, iy - 3, iz);
				const T* yp1 = GetHeadYR(ix, iy + 1, iz);
				const T* yp2 = GetHeadYR(ix, iy + 2, iz);
				const T* yp3 = GetHeadYR(ix, iy + 3, iz);
#pragma ivdep
				for (int n = 0; n < Ns; ++n) {
					//auto d2psidx2 = c0 * p0[n];
					//d2psidx2 += ((xm3[n] + xp3[n]) - 54.0 * (xm2[n] + xp2[n]) + 783.0 * (xm1[n] + xp1[n])) / (dx * dx);
					const auto d2psidx2 = ((ym3[n] + yp3[n]) - 54.0 * (ym2[n] + yp2[n]) + 783.0 * (ym1[n] + yp1[n])) / (dy * dy);
					//d2psidx2 += ((zm3[n] + zp3[n]) - 54.0 * (zm2[n] + zp2[n]) + 783.0 * (zm1[n] + zp1[n])) / (dz * dz);
					out[n + i] += coef * d2psidx2 / (24.0 * 24.0);
				}
			}
		}
	}

	for (int iz = 0; iz < size_z; ++iz) {
		for (int iy = 0; iy < size_y; ++iy) {
			for (int ix = 0; ix < size_x; ++ix) {
				const int64_t i = Ns * (ix + size_x * (iy + (size_y * iz)));
				const T* p0 = p + i;
				const T* xm1 = GetHeadXL(ix - 1, iy, iz);
				const T* xm2 = GetHeadXL(ix - 2, iy, iz);
				const T* xm3 = GetHeadXL(ix - 3, iy, iz);
				const T* xp1 = GetHeadXR(ix + 1, iy, iz);
				const T* xp2 = GetHeadXR(ix + 2, iy, iz);
				const T* xp3 = GetHeadXR(ix + 3, iy, iz);
				const T* zm1 = GetHeadZL(ix, iy, iz - 1);
				const T* zm2 = GetHeadZL(ix, iy, iz - 2);
				const T* zm3 = GetHeadZL(ix, iy, iz - 3);
				const T* zp1 = GetHeadZR(ix, iy, iz + 1);
				const T* zp2 = GetHeadZR(ix, iy, iz + 2);
				const T* zp3 = GetHeadZR(ix, iy, iz + 3);
				const T* ym1 = GetHeadYL(ix, iy - 1, iz);
				const T* ym2 = GetHeadYL(ix, iy - 2, iz);
				const T* ym3 = GetHeadYL(ix, iy - 3, iz);
				const T* yp1 = GetHeadYR(ix, iy + 1, iz);
				const T* yp2 = GetHeadYR(ix, iy + 2, iz);
				const T* yp3 = GetHeadYR(ix, iy + 3, iz);
#pragma ivdep
				for (int n = 0; n < Ns; ++n) {
					//auto d2psidx2 = c0 * p0[n];
					//d2psidx2 += ((xm3[n] + xp3[n]) - 54.0 * (xm2[n] + xp2[n]) + 783.0 * (xm1[n] + xp1[n])) / (dx * dx);
					//d2psidx2 += ((ym3[n] + yp3[n]) - 54.0 * (ym2[n] + yp2[n]) + 783.0 * (ym1[n] + yp1[n])) / (dy * dy);
					const auto d2psidx2 = ((zm3[n] + zp3[n]) - 54.0 * (zm2[n] + zp2[n]) + 783.0 * (zm1[n] + zp1[n])) / (dz * dz);
					out[n + i] += coef * d2psidx2 / (24.0 * 24.0);
				}
			}
		}
	}

}


template <class T, class FUNC_XL, class FUNC_XR, class FUNC_YL, class FUNC_YR, class FUNC_ZL, class FUNC_ZR>
void Laplasian4th_bundle_X3c(int64_t Ns, int size_x, int size_y, int size_z, T* out, const T* p, FUNC_XL&& func_xl, FUNC_XR&& func_xr, FUNC_YL&& func_yl, FUNC_YR&& func_yr, FUNC_ZL&& func_zl, FUNC_ZR&& func_zr, double coef, double dx, double dy, double dz) {
	auto pow2 = [](T a) {return a * a; };
	const double c0 = -1460.0 * coef * (1.0/ (dx * dx) + 1.0/ (dy * dy) + 1.0/ (dz * dz)) / (24.0 * 24.0);
	const double cx1 = coef / pow2(24.0 * dx);
	const double cx2 = -54.0 * coef / pow2(24.0 * dx);
	const double cx3 = 783.0 * coef / pow2(24.0 * dx);
	const double cy1 = coef / pow2(24.0 * dy);
	const double cy2 = -54.0 * coef / pow2(24.0 * dy);
	const double cy3 = 783.0 * coef / pow2(24.0 * dy);
	const double cz1 = coef / pow2(24.0 * dz);
	const double cz2 = -54.0 * coef / pow2(24.0 * dz);
	const double cz3 = 783.0 * coef / pow2(24.0 * dz);

	auto GetHeadXL = [&](int ix, int iy, int iz) {
		if (ix < 0) return func_xl(ix, iy, iz);
		return p + Ns * (ix + size_x * (iy + (size_y * iz)));
		};
	auto GetHeadXR = [&](int ix, int iy, int iz) {
		if (ix >= size_x) return func_xr(ix, iy, iz);
		return p + Ns * (ix + size_x * (iy + (size_y * iz)));
		};
	auto GetHeadYL = [&](int ix, int iy, int iz) {
		if (iy < 0) return func_yl(ix, iy, iz);
		return p + Ns * (ix + size_x * (iy + (size_y * iz)));
		};
	auto GetHeadYR = [&](int ix, int iy, int iz) {
		if (iy >= size_y) return func_yr(ix, iy, iz);
		return p + Ns * (ix + size_x * (iy + (size_y * iz)));
		};
	auto GetHeadZL = [&](int ix, int iy, int iz) {
		if (iz < 0) return func_zl(ix, iy, iz);
		return p + Ns * (ix + size_x * (iy + (size_y * iz)));
		};
	auto GetHeadZR = [&](int ix, int iy, int iz) {
		if (iz >= size_z) return func_zr(ix, iy, iz);
		return p + Ns * (ix + size_x * (iy + (size_y * iz)));
		};


	for (int iz = 0; iz < size_z; ++iz) {
		for (int iy = 0; iy < size_y; ++iy) {
			for (int ix = 0; ix < size_x; ++ix) {
				const int64_t i = Ns * (ix + size_x * (iy + (size_y * iz)));
				const T* p0 = p + i;
				const T* xm1 = GetHeadXL(ix - 1, iy, iz);
				const T* xm2 = GetHeadXL(ix - 2, iy, iz);
				const T* xm3 = GetHeadXL(ix - 3, iy, iz);
				const T* xp1 = GetHeadXR(ix + 1, iy, iz);
				const T* xp2 = GetHeadXR(ix + 2, iy, iz);
				const T* xp3 = GetHeadXR(ix + 3, iy, iz);
				const T* zm1 = GetHeadZL(ix, iy, iz - 1);
				const T* zm2 = GetHeadZL(ix, iy, iz - 2);
				const T* zm3 = GetHeadZL(ix, iy, iz - 3);
				const T* zp1 = GetHeadZR(ix, iy, iz + 1);
				const T* zp2 = GetHeadZR(ix, iy, iz + 2);
				const T* zp3 = GetHeadZR(ix, iy, iz + 3);
				const T* ym1 = GetHeadYL(ix, iy - 1, iz);
				const T* ym2 = GetHeadYL(ix, iy - 2, iz);
				const T* ym3 = GetHeadYL(ix, iy - 3, iz);
				const T* yp1 = GetHeadYR(ix, iy + 1, iz);
				const T* yp2 = GetHeadYR(ix, iy + 2, iz);
				const T* yp3 = GetHeadYR(ix, iy + 3, iz);
#pragma ivdep
				for (int n = 0; n < Ns; ++n) {
					const auto d2psidx2 = c0 * p0[n]
						+ cx1 * (xm3[n] + xp3[n]) + cx2 * (xm2[n] + xp2[n]) + cx3 * (xm1[n] + xp1[n]);
					//d2psidx2 += ((ym3[n] + yp3[n]) - 54.0 * (ym2[n] + yp2[n]) + 783.0 * (ym1[n] + yp1[n])) / (dy * dy);
					//d2psidx2 += ((zm3[n] + zp3[n]) - 54.0 * (zm2[n] + zp2[n]) + 783.0 * (zm1[n] + zp1[n])) / (dz * dz);
					out[n + i] += d2psidx2;
				}
			}
		}
	}
	for (int iz = 0; iz < size_z; ++iz) {
		for (int iy = 0; iy < size_y; ++iy) {
			for (int ix = 0; ix < size_x; ++ix) {
				const int64_t i = Ns * (ix + size_x * (iy + (size_y * iz)));
				const T* p0 = p + i;
				const T* xm1 = GetHeadXL(ix - 1, iy, iz);
				const T* xm2 = GetHeadXL(ix - 2, iy, iz);
				const T* xm3 = GetHeadXL(ix - 3, iy, iz);
				const T* xp1 = GetHeadXR(ix + 1, iy, iz);
				const T* xp2 = GetHeadXR(ix + 2, iy, iz);
				const T* xp3 = GetHeadXR(ix + 3, iy, iz);
				const T* zm1 = GetHeadZL(ix, iy, iz - 1);
				const T* zm2 = GetHeadZL(ix, iy, iz - 2);
				const T* zm3 = GetHeadZL(ix, iy, iz - 3);
				const T* zp1 = GetHeadZR(ix, iy, iz + 1);
				const T* zp2 = GetHeadZR(ix, iy, iz + 2);
				const T* zp3 = GetHeadZR(ix, iy, iz + 3);
				const T* ym1 = GetHeadYL(ix, iy - 1, iz);
				const T* ym2 = GetHeadYL(ix, iy - 2, iz);
				const T* ym3 = GetHeadYL(ix, iy - 3, iz);
				const T* yp1 = GetHeadYR(ix, iy + 1, iz);
				const T* yp2 = GetHeadYR(ix, iy + 2, iz);
				const T* yp3 = GetHeadYR(ix, iy + 3, iz);
#pragma ivdep
				for (int n = 0; n < Ns; ++n) {
					//const auto d2psidx2 = c0 * p0[n];
//					d2psidx2 += ((xm3[n] + xp3[n]) - 54.0 * (xm2[n] + xp2[n]) + 783.0 * (xm1[n] + xp1[n])) / (dx * dx);
//					d2psidx2 += ((ym3[n] + yp3[n]) - 54.0 * (ym2[n] + yp2[n]) + 783.0 * (ym1[n] + yp1[n])) / (dy * dy);
//					d2psidx2 += ((zm3[n] + zp3[n]) - 54.0 * (zm2[n] + zp2[n]) + 783.0 * (zm1[n] + zp1[n])) / (dz * dz);
					//out[n + i] += coef * d2psidx2 / (24.0 * 24.0);
				}
			}
		}
	}

	for (int iz = 0; iz < size_z; ++iz) {
		for (int iy = 0; iy < size_y; ++iy) {
			for (int ix = 0; ix < size_x; ++ix) {
				const int64_t i = Ns * (ix + size_x * (iy + (size_y * iz)));
				const T* p0 = p + i;
				const T* xm1 = GetHeadXL(ix - 1, iy, iz);
				const T* xm2 = GetHeadXL(ix - 2, iy, iz);
				const T* xm3 = GetHeadXL(ix - 3, iy, iz);
				const T* xp1 = GetHeadXR(ix + 1, iy, iz);
				const T* xp2 = GetHeadXR(ix + 2, iy, iz);
				const T* xp3 = GetHeadXR(ix + 3, iy, iz);
				const T* zm1 = GetHeadZL(ix, iy, iz - 1);
				const T* zm2 = GetHeadZL(ix, iy, iz - 2);
				const T* zm3 = GetHeadZL(ix, iy, iz - 3);
				const T* zp1 = GetHeadZR(ix, iy, iz + 1);
				const T* zp2 = GetHeadZR(ix, iy, iz + 2);
				const T* zp3 = GetHeadZR(ix, iy, iz + 3);
				const T* ym1 = GetHeadYL(ix, iy - 1, iz);
				const T* ym2 = GetHeadYL(ix, iy - 2, iz);
				const T* ym3 = GetHeadYL(ix, iy - 3, iz);
				const T* yp1 = GetHeadYR(ix, iy + 1, iz);
				const T* yp2 = GetHeadYR(ix, iy + 2, iz);
				const T* yp3 = GetHeadYR(ix, iy + 3, iz);
#pragma ivdep
				for (int n = 0; n < Ns; ++n) {
					//auto d2psidx2 = c0 * p0[n];
					//d2psidx2 += ((xm3[n] + xp3[n]) - 54.0 * (xm2[n] + xp2[n]) + 783.0 * (xm1[n] + xp1[n])) / (dx * dx);
					const auto d2psidx2 = cy1 * (ym3[n] + yp3[n]) + cy2 * (ym2[n] + yp2[n]) + cy3 * (ym1[n] + yp1[n]);
					//d2psidx2 += ((zm3[n] + zp3[n]) - 54.0 * (zm2[n] + zp2[n]) + 783.0 * (zm1[n] + zp1[n])) / (dz * dz);
					out[n + i] += d2psidx2;
				}
			}
		}
	}

	for (int iz = 0; iz < size_z; ++iz) {
		for (int iy = 0; iy < size_y; ++iy) {
			for (int ix = 0; ix < size_x; ++ix) {
				const int64_t i = Ns * (ix + size_x * (iy + (size_y * iz)));
				const T* p0 = p + i;
				const T* xm1 = GetHeadXL(ix - 1, iy, iz);
				const T* xm2 = GetHeadXL(ix - 2, iy, iz);
				const T* xm3 = GetHeadXL(ix - 3, iy, iz);
				const T* xp1 = GetHeadXR(ix + 1, iy, iz);
				const T* xp2 = GetHeadXR(ix + 2, iy, iz);
				const T* xp3 = GetHeadXR(ix + 3, iy, iz);
				const T* zm1 = GetHeadZL(ix, iy, iz - 1);
				const T* zm2 = GetHeadZL(ix, iy, iz - 2);
				const T* zm3 = GetHeadZL(ix, iy, iz - 3);
				const T* zp1 = GetHeadZR(ix, iy, iz + 1);
				const T* zp2 = GetHeadZR(ix, iy, iz + 2);
				const T* zp3 = GetHeadZR(ix, iy, iz + 3);
				const T* ym1 = GetHeadYL(ix, iy - 1, iz);
				const T* ym2 = GetHeadYL(ix, iy - 2, iz);
				const T* ym3 = GetHeadYL(ix, iy - 3, iz);
				const T* yp1 = GetHeadYR(ix, iy + 1, iz);
				const T* yp2 = GetHeadYR(ix, iy + 2, iz);
				const T* yp3 = GetHeadYR(ix, iy + 3, iz);
#pragma ivdep
				for (int n = 0; n < Ns; ++n) {
					//auto d2psidx2 = c0 * p0[n];
					//d2psidx2 += ((xm3[n] + xp3[n]) - 54.0 * (xm2[n] + xp2[n]) + 783.0 * (xm1[n] + xp1[n])) / (dx * dx);
					//d2psidx2 += ((ym3[n] + yp3[n]) - 54.0 * (ym2[n] + yp2[n]) + 783.0 * (ym1[n] + yp1[n])) / (dy * dy);
					const auto d2psidx2 = cz1 * (zm3[n] + zp3[n]) + cz2 * (zm2[n] + zp2[n]) + cz3 * (zm1[n] + zp1[n]);
					out[n + i] += d2psidx2;
				}
			}
		}
	}

}



template <class T, class FUNC_XL, class FUNC_XR, class FUNC_YL, class FUNC_YR, class FUNC_ZL, class FUNC_ZR>
void Laplasian4th_bundle_X4(int64_t Ns, int size_x, int size_y, int size_z, T* out, const T* p, FUNC_XL&& func_xl, FUNC_XR&& func_xr, FUNC_YL&& func_yl, FUNC_YR&& func_yr, FUNC_ZL&& func_zl, FUNC_ZR&& func_zr, double coef, double dx, double dy, double dz) {
	auto pow2 = [](T a) {return a * a; };
	const double c0 = -1460.0 * coef * (1.0 / (dx * dx) + 1.0 / (dy * dy) + 1.0 / (dz * dz)) / (24.0 * 24.0);
	const double cx1 = coef / pow2(24.0 * dx);
	const double cx2 = -54.0 * coef / pow2(24.0 * dx);
	const double cx3 = 783.0 * coef / pow2(24.0 * dx);
	const double cy1 = coef / pow2(24.0 * dy);
	const double cy2 = -54.0 * coef / pow2(24.0 * dy);
	const double cy3 = 783.0 * coef / pow2(24.0 * dy);
	const double cz1 = coef / pow2(24.0 * dz);
	const double cz2 = -54.0 * coef / pow2(24.0 * dz);
	const double cz3 = 783.0 * coef / pow2(24.0 * dz);

	auto GetHeadXL = [&](int ix, int iy, int iz) {
		if (ix < 0) return func_xl(ix, iy, iz);
		return p + Ns * (ix + size_x * (iy + (size_y * iz)));
		};
	auto GetHeadXR = [&](int ix, int iy, int iz) {
		if (ix >= size_x) return func_xr(ix, iy, iz);
		return p + Ns * (ix + size_x * (iy + (size_y * iz)));
		};
	auto GetHeadYL = [&](int ix, int iy, int iz) {
		if (iy < 0) return func_yl(ix, iy, iz);
		return p + Ns * (ix + size_x * (iy + (size_y * iz)));
		};
	auto GetHeadYR = [&](int ix, int iy, int iz) {
		if (iy >= size_y) return func_yr(ix, iy, iz);
		return p + Ns * (ix + size_x * (iy + (size_y * iz)));
		};
	auto GetHeadZL = [&](int ix, int iy, int iz) {
		if (iz < 0) return func_zl(ix, iy, iz);
		return p + Ns * (ix + size_x * (iy + (size_y * iz)));
		};
	auto GetHeadZR = [&](int ix, int iy, int iz) {
		if (iz >= size_z) return func_zr(ix, iy, iz);
		return p + Ns * (ix + size_x * (iy + (size_y * iz)));
		};



	for (int iz = 0; iz < size_z; ++iz) {
		for (int iy = 0; iy < size_y; ++iy) {
			for (int ix = 0; ix < size_x; ++ix) {
				const int64_t i = Ns * (ix + size_x * (iy + (size_y * iz)));
				const T* p0 = p + i;
				const T* xm1 = GetHeadXL(ix - 1, iy, iz);
				const T* xm2 = GetHeadXL(ix - 2, iy, iz);
				const T* xm3 = GetHeadXL(ix - 3, iy, iz);
				const T* xp1 = GetHeadXR(ix + 1, iy, iz);
				const T* xp2 = GetHeadXR(ix + 2, iy, iz);
				const T* xp3 = GetHeadXR(ix + 3, iy, iz);
				const T* zm1 = GetHeadZL(ix, iy, iz - 1);
				const T* zm2 = GetHeadZL(ix, iy, iz - 2);
				const T* zm3 = GetHeadZL(ix, iy, iz - 3);
				const T* zp1 = GetHeadZR(ix, iy, iz + 1);
				const T* zp2 = GetHeadZR(ix, iy, iz + 2);
				const T* zp3 = GetHeadZR(ix, iy, iz + 3);
				const T* ym1 = GetHeadYL(ix, iy - 1, iz);
				const T* ym2 = GetHeadYL(ix, iy - 2, iz);
				const T* ym3 = GetHeadYL(ix, iy - 3, iz);
				const T* yp1 = GetHeadYR(ix, iy + 1, iz);
				const T* yp2 = GetHeadYR(ix, iy + 2, iz);
				const T* yp3 = GetHeadYR(ix, iy + 3, iz);
#pragma ivdep
				for (int n = 0; n < Ns; ++n) {
					const auto d2psidx2 = cz1 * (zm3[n]) + cz2 * (zm2[n]) + cz3 * (zm1[n]);
					out[n + i] += d2psidx2;
				}
			}
		}
	}


	for (int iz = 0; iz < size_z; ++iz) {
		for (int iy = 0; iy < size_y; ++iy) {
			for (int ix = 0; ix < size_x; ++ix) {
				const int64_t i = Ns * (ix + size_x * (iy + (size_y * iz)));
				const T* p0 = p + i;
				const T* xm1 = GetHeadXL(ix - 1, iy, iz);
				const T* xm2 = GetHeadXL(ix - 2, iy, iz);
				const T* xm3 = GetHeadXL(ix - 3, iy, iz);
				const T* xp1 = GetHeadXR(ix + 1, iy, iz);
				const T* xp2 = GetHeadXR(ix + 2, iy, iz);
				const T* xp3 = GetHeadXR(ix + 3, iy, iz);
				const T* zm1 = GetHeadZL(ix, iy, iz - 1);
				const T* zm2 = GetHeadZL(ix, iy, iz - 2);
				const T* zm3 = GetHeadZL(ix, iy, iz - 3);
				const T* zp1 = GetHeadZR(ix, iy, iz + 1);
				const T* zp2 = GetHeadZR(ix, iy, iz + 2);
				const T* zp3 = GetHeadZR(ix, iy, iz + 3);
				const T* ym1 = GetHeadYL(ix, iy - 1, iz);
				const T* ym2 = GetHeadYL(ix, iy - 2, iz);
				const T* ym3 = GetHeadYL(ix, iy - 3, iz);
				const T* yp1 = GetHeadYR(ix, iy + 1, iz);
				const T* yp2 = GetHeadYR(ix, iy + 2, iz);
				const T* yp3 = GetHeadYR(ix, iy + 3, iz);
#pragma ivdep
				for (int n = 0; n < Ns; ++n) {
					const auto d2psidx2 = c0 * p0[n]
						+ cx1 * (xm3[n] + xp3[n]) + cx2 * (xm2[n] + xp2[n]) + cx3 * (xm1[n] + xp1[n]);
					out[n + i] += d2psidx2;
				}
			}
		}
	}


	for (int iz = 0; iz < size_z; ++iz) {
		for (int iy = 0; iy < size_y; ++iy) {
			for (int ix = 0; ix < size_x; ++ix) {
				const int64_t i = Ns * (ix + size_x * (iy + (size_y * iz)));
				const T* p0 = p + i;
				const T* xm1 = GetHeadXL(ix - 1, iy, iz);
				const T* xm2 = GetHeadXL(ix - 2, iy, iz);
				const T* xm3 = GetHeadXL(ix - 3, iy, iz);
				const T* xp1 = GetHeadXR(ix + 1, iy, iz);
				const T* xp2 = GetHeadXR(ix + 2, iy, iz);
				const T* xp3 = GetHeadXR(ix + 3, iy, iz);
				const T* zm1 = GetHeadZL(ix, iy, iz - 1);
				const T* zm2 = GetHeadZL(ix, iy, iz - 2);
				const T* zm3 = GetHeadZL(ix, iy, iz - 3);
				const T* zp1 = GetHeadZR(ix, iy, iz + 1);
				const T* zp2 = GetHeadZR(ix, iy, iz + 2);
				const T* zp3 = GetHeadZR(ix, iy, iz + 3);
				const T* ym1 = GetHeadYL(ix, iy - 1, iz);
				const T* ym2 = GetHeadYL(ix, iy - 2, iz);
				const T* ym3 = GetHeadYL(ix, iy - 3, iz);
				const T* yp1 = GetHeadYR(ix, iy + 1, iz);
				const T* yp2 = GetHeadYR(ix, iy + 2, iz);
				const T* yp3 = GetHeadYR(ix, iy + 3, iz);
#pragma ivdep
				for (int n = 0; n < Ns; ++n) {
					//auto d2psidx2 = c0 * p0[n];
					//d2psidx2 += ((xm3[n] + xp3[n]) - 54.0 * (xm2[n] + xp2[n]) + 783.0 * (xm1[n] + xp1[n])) / (dx * dx);
					const auto d2psidx2 = cy1 * (ym3[n] + yp3[n]) + cy2 * (ym2[n] + yp2[n]) + cy3 * (ym1[n] + yp1[n]);
					//d2psidx2 += ((zm3[n] + zp3[n]) - 54.0 * (zm2[n] + zp2[n]) + 783.0 * (zm1[n] + zp1[n])) / (dz * dz);
					out[n + i] += d2psidx2;
				}
			}
		}
	}



	for (int iz = 0; iz < size_z; ++iz) {
		for (int iy = 0; iy < size_y; ++iy) {
			for (int ix = 0; ix < size_x; ++ix) {
				const int64_t i = Ns * (ix + size_x * (iy + (size_y * iz)));
				const T* p0 = p + i;
				const T* xm1 = GetHeadXL(ix - 1, iy, iz);
				const T* xm2 = GetHeadXL(ix - 2, iy, iz);
				const T* xm3 = GetHeadXL(ix - 3, iy, iz);
				const T* xp1 = GetHeadXR(ix + 1, iy, iz);
				const T* xp2 = GetHeadXR(ix + 2, iy, iz);
				const T* xp3 = GetHeadXR(ix + 3, iy, iz);
				const T* zm1 = GetHeadZL(ix, iy, iz - 1);
				const T* zm2 = GetHeadZL(ix, iy, iz - 2);
				const T* zm3 = GetHeadZL(ix, iy, iz - 3);
				const T* zp1 = GetHeadZR(ix, iy, iz + 1);
				const T* zp2 = GetHeadZR(ix, iy, iz + 2);
				const T* zp3 = GetHeadZR(ix, iy, iz + 3);
				const T* ym1 = GetHeadYL(ix, iy - 1, iz);
				const T* ym2 = GetHeadYL(ix, iy - 2, iz);
				const T* ym3 = GetHeadYL(ix, iy - 3, iz);
				const T* yp1 = GetHeadYR(ix, iy + 1, iz);
				const T* yp2 = GetHeadYR(ix, iy + 2, iz);
				const T* yp3 = GetHeadYR(ix, iy + 3, iz);
#pragma ivdep
				for (int n = 0; n < Ns; ++n) {
					const auto d2psidx2 = cz1 * (zp3[n]) + cz2 * (zp2[n]) + cz3 * (zp1[n]);
					out[n + i] += d2psidx2;
				}
			}
		}
	}

}



//periodic//
template <class T>
void Laplasian4th_bundle(const GridRange& grid, int64_t Ns, T* out, const T* p, const double coef, const double dx, const double dy, const double dz) {
	const int size_x = grid.SizeX();
	const int size_y = grid.SizeY();
	const int size_z = grid.SizeZ();

	
	Laplasian4th_bundle_X1c<T>(Ns, size_x, size_y, size_z, out, p, 
		[&](int ix, int iy, int iz) { return p + Ns * (ix + size_x * (iy + 1 + size_y * iz));},
		[&](int ix, int iy, int iz) { return p + Ns * (ix + size_x * (iy - 1 + size_y * iz));},
		[&](int ix, int iy, int iz) { return p + Ns * (ix + size_x * (iy + size_y * (iz + 1))); },
		[&](int ix, int iy, int iz) { return p + Ns * (ix + size_x * (iy + size_y * (iz - 1))); },
		[&](int ix, int iy, int iz) { return p + Ns * (ix + size_x * (iy + size_y * (iz + size_z))); },
		[&](int ix, int iy, int iz) { return p + Ns * (ix + size_x * (iy + size_y * (iz - size_z))); },
		coef, dx, dy, dz);
}

#ifdef USE_MPI


/*
* 1次元領域分割の場合のLaplasian
* 分割方向はz方向
*/
template <class T>
void Laplasian4th_MPI_1D_overspace(const GridRangeMPI& grid, T* out, const T* src, T* buf_with_halo, double coef, double dx, double dy, double dz) {

	//のり代送受信//
	constexpr int TAG = 0x1;
	const int size_x = grid.SizeX();
	const int size_y = grid.SizeY();
	const int size_xy = size_x * size_y;
	const int size_z = grid.SizeZ();
	const int margin_width = 3;
	static std::vector<T> margin_data(size_xy * margin_width * 2);//margin data//

	TransferMarginZ(grid, margin_data.data(), margin_data.data() + size_xy * margin_width, src, margin_width);

#ifdef TEST_PRINT
	printf("Laplasian4th_MPI_1D_halo\n"); fflush(stdout);
#endif
	PasteHalo_f3(grid, buf_with_halo, src, 
		[&](int ix, int iy, int iz) { return src[(ix + size_x * (iy + 1 + size_y * iz))]; },
		[&](int ix, int iy, int iz) { return src[(ix + size_x * (iy - 1 + size_y * iz))]; },
		[&](int ix, int iy, int iz) { return src[(ix + size_x * (iy + size_y * (iz + 1)))]; },
		[&](int ix, int iy, int iz) { return src[(ix + size_x * (iy + size_y * (iz - 1)))]; },
		[&](int ix, int iy, int iz) { return margin_data[(ix + size_x * (iy + size_y * (iz + margin_width)))]; },
		[&](int ix, int iy, int iz) { return margin_data[size_xy * margin_width + (ix + size_x * (iy + size_y * (iz - size_z)))]; });
		

	Laplasian4th_halo_f(grid, out, buf_with_halo, coef, dx, dy, dz);


}


/*
* 1次元領域分割の場合のLaplasian
* 分割方向はz方向
*/
template <class T>
void Laplasian4th_MPI_3D_overspace(const GridRangeMPI& grid, T* out, const T* src, T* buf_with_halo, double coef, double dx, double dy, double dz) {

	//のり代送受信//
	constexpr int TAG = 0x1;
	const int size_x = grid.SizeX();
	const int size_y = grid.SizeY();
	const int size_xy = size_x * size_y;
	const int size_z = grid.SizeZ();
	const int size_xz = size_x * size_z;
	const int size_yz = size_y * size_z;
	const int margin_width = 3;
	static std::vector<T> margin_data((size_xy + size_xz + size_yz) * margin_width * 2);//margin data//

	T* margin_x_left = margin_data.data();
	T* margin_x_right = margin_x_left + size_yz * margin_width;
	T* margin_y_left = margin_x_right + size_yz * margin_width;
	T* margin_y_right = margin_y_left + size_xz * margin_width;
	T* margin_z_left = margin_y_right + size_xz * margin_width;
	T* margin_z_right = margin_z_left + size_xy * margin_width;

#ifdef TEST_PRINT
	int proc_id;
	MPI_Comm_rank(grid.mpi_comm, &proc_id);
	printf("[%d]test4-1\n", proc_id); fflush(stdout);
#endif
	TransferMarginX(grid, margin_x_left, margin_x_right, src, margin_width);
#ifdef TEST_PRINT
	printf("[%d]test4-2\n", proc_id); fflush(stdout);
#endif
	TransferMarginY(grid, margin_y_left, margin_y_right, src, margin_width);

#ifdef TEST_PRINT
	MPI_Barrier(grid.mpi_comm);
	printf("[%d]test4-3\n", proc_id); fflush(stdout);
#endif
	TransferMarginZ(grid, margin_z_left, margin_z_right, src, margin_width);

#ifdef TEST_PRINT
	MPI_Barrier(grid.mpi_comm);
	printf("[%d]test4-4\n", proc_id); fflush(stdout);
#endif

	PasteHalo_f3(grid, buf_with_halo, src,

		//[&](int ix, int iy, int iz) { return src[(ix + size_x * (iy + 1 + size_y * iz))]; },
		[&](int ix, int iy, int iz) { return margin_x_left[ix + margin_width + margin_width * (iy + size_y * iz)]; },

		//[&](int ix, int iy, int iz) { return src[(ix + size_x * (iy - 1 + size_y * iz))]; },
		[&](int ix, int iy, int iz) { return margin_x_right[ix - size_x + margin_width * (iy + size_y * iz)]; },

		//[&](int ix, int iy, int iz) { return src[(ix + size_x * (iy + size_y * (iz + 1)))]; },
		[&](int ix, int iy, int iz) { return margin_y_left[ix + size_x * (iy + margin_width + margin_width * iz)]; },

		//[&](int ix, int iy, int iz) { return src[(ix + size_x * (iy + size_y * (iz - 1)))]; },
		[&](int ix, int iy, int iz) { return margin_y_right[ix + size_x * (iy - size_y + margin_width * (iz))]; },

		[&](int ix, int iy, int iz) { return margin_z_left[ix + size_x * (iy + size_y * (iz + margin_width))]; },
		[&](int ix, int iy, int iz) { return margin_z_right[ix + size_x * (iy + size_y * (iz - size_z))]; });
		
	Laplasian4th_halo_f(grid, out, buf_with_halo, coef, dx, dy, dz);

#ifdef TEST_PRINT
	MPI_Barrier(grid.mpi_comm);
	printf("[%d]test4-5\n", proc_id); fflush(stdout);
#endif
}


template <class T>
void Laplasian4th_ddm(const GridRangeMPI& grid, T* out, const T* src, double coef, double dx, double dy, double dz) {

	const int margin_width = 3;
	const int over_size = (grid.SizeX() + 2 * margin_width) * (grid.SizeY() + 2 * margin_width) * (grid.SizeZ() + 2 * margin_width);
	static std::vector<T> buf_with_halo(over_size);//margin data//

#ifdef TEST_PRINT
	printf("L4 split_dimension = %d\n", grid.split_dimension);
#endif
	switch (grid.split_dimension) {
	case 1:
		Laplasian4th_MPI_1D_overspace(grid, out, src, buf_with_halo.data(), coef, dx, dy, dz);
		break;
	case 2:

		Laplasian4th_MPI_3D_overspace(grid, out, src, buf_with_halo.data(), coef, dx, dy, dz);
		break;
	case 3:
		
		Laplasian4th_MPI_3D_overspace(grid, out, src, buf_with_halo.data(), coef, dx, dy, dz);
		break;
	default:
		Laplasian4th_halo(grid, out, src, buf_with_halo.data(), coef, dx, dy, dz);
		break;
	}
}

#endif
