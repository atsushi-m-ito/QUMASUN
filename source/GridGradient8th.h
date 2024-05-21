#pragma once
#ifdef USE_MPI
#include <mpi.h>
#include "mpi_helper.h"
#endif
#include "GridRange.h"
#include <vector>
#include <cstdio>

#include "GridDifference8th.h"


//3次元デカルト座標グリッドの計算のためのループを司る機能//
//MPI並列時にはリダクションも提供する




template <class T>
void Gradient8th_x_halo_f(const GridRange& grid, T* out, const T* p, double coef, double dx, double dy, double dz) {
	
	const int Nx = grid.SizeX();
	const int Ny = grid.SizeY();
	const int Nz = grid.SizeZ();
	constexpr int MW = 4;

	const double c4x = -3.0 * coef / (840.0 * dx);
	const double c3x = 32.0*coef / (840.0 * dx);
	const double c2x = -168.0 * coef / (840.0 * dx);
	const double c1x = 672.0 * coef / (840.0 * dx);
	
	const int size_x = grid.SizeX() + MW*2;
	const int size_y = grid.SizeY() + MW * 2;
	const int size_z = grid.SizeZ() + MW * 2;

	for (int iz = MW; iz < size_z - MW; ++iz) {
		const int idz_m4 = -4 * size_x * size_y;
		const int idz_m3 = -3 * size_x * size_y;
		const int idz_m2 = -2 * size_x * size_y;
		const int idz_m = -1 * size_x * size_y;
		const int idz_p = 1 * size_x * size_y;
		const int idz_p2 = 2 * size_x * size_y;
		const int idz_p3 = 3 * size_x * size_y;
		const int idz_p4 = 4 * size_x * size_y;

		for (int iy = MW; iy < size_y - MW; ++iy) {
			const int idy_m4 = -4 * size_x;
			const int idy_m3 = -3 * size_x;
			const int idy_m2 = -2 * size_x;
			const int idy_m = -1 * size_x;
			const int idy_p = 1 * size_x;
			const int idy_p2 = 2 * size_x;
			const int idy_p3 = 3 * size_x;
			const int idy_p4 = 4 * size_x;
#pragma ivdep
			for (int ix = MW; ix < size_x - MW; ++ix) {
				const int idx_m4 = -4;
				const int idx_m3 = -3;
				const int idx_m2 = -2;
				const int idx_m = -1;
				const int idx_p = 1;
				const int idx_p2 = 2;
				const int idx_p3 = 3;
				const int idx_p4 = 4;
				const int64_t i = ix + size_x * (iy + (size_y * iz));
				const int64_t io = (ix - MW) + Nx * (iy - MW + (Ny * (iz - MW)));


				auto dpsi_dx = c4x * (p[i + idx_p4] - p[i + idx_m4]) + c3x * (p[i + idx_p3] - p[i + idx_m3]) + c2x * (p[i + idx_p2] - p[i + idx_m2]) + c1x * (p[i + idx_p] - p[i + idx_m]);
				out[io] = dpsi_dx;
			}
		}
	}
}


template <class T>
void Gradient8th_y_halo_f(const GridRange& grid, T* out, const T* p, double coef, double dx, double dy, double dz) {

	const int Nx = grid.SizeX();
	const int Ny = grid.SizeY();
	const int Nz = grid.SizeZ();
	constexpr int MW = 4;

	const double c4y = -3.0 * coef / (840.0 * dy);
	const double c3y = 32.0 * coef / (840.0 * dy);
	const double c2y = -168.0 * coef / (840.0 * dy);
	const double c1y = 672.0 * coef / (840.0 * dy);

	const int size_x = grid.SizeX() + MW * 2;
	const int size_y = grid.SizeY() + MW * 2;
	const int size_z = grid.SizeZ() + MW * 2;

	for (int iz = MW; iz < size_z - MW; ++iz) {
		const int idz_m4 = -4 * size_x * size_y;
		const int idz_m3 = -3 * size_x * size_y;
		const int idz_m2 = -2 * size_x * size_y;
		const int idz_m = -1 * size_x * size_y;
		const int idz_p = 1 * size_x * size_y;
		const int idz_p2 = 2 * size_x * size_y;
		const int idz_p3 = 3 * size_x * size_y;
		const int idz_p4 = 4 * size_x * size_y;

		for (int iy = MW; iy < size_y - MW; ++iy) {
			const int idy_m4 = -4 * size_x;
			const int idy_m3 = -3 * size_x;
			const int idy_m2 = -2 * size_x;
			const int idy_m = -1 * size_x;
			const int idy_p = 1 * size_x;
			const int idy_p2 = 2 * size_x;
			const int idy_p3 = 3 * size_x;
			const int idy_p4 = 4 * size_x;
#pragma ivdep
			for (int ix = MW; ix < size_x - MW; ++ix) {
				const int idx_m4 = -4;
				const int idx_m3 = -3;
				const int idx_m2 = -2;
				const int idx_m = -1;
				const int idx_p = 1;
				const int idx_p2 = 2;
				const int idx_p3 = 3;
				const int idx_p4 = 4;
				const int64_t i = ix + size_x * (iy + (size_y * iz));
				const int64_t io = (ix - MW) + Nx * (iy - MW + (Ny * (iz - MW)));


				auto dpsi_dy = c4y * (p[i + idy_p4] - p[i + idy_m4]) + c3y * (p[i + idy_p3] - p[i + idy_m3]) + c2y * (p[i + idy_p2] - p[i + idy_m2]) + c1y * (p[i + idy_p] - p[i + idy_m]);
				out[io] = dpsi_dy;
			}
		}
	}
}


template <class T>
void Gradient8th_z_halo_f(const GridRange& grid, T* out, const T* p, double coef, double dx, double dy, double dz) {

	const int Nx = grid.SizeX();
	const int Ny = grid.SizeY();
	const int Nz = grid.SizeZ();
	constexpr int MW = 4;

	const double c4z = -3.0 * coef / (840.0 * dz);
	const double c3z = 32.0 * coef / (840.0 * dz);
	const double c2z = -168.0 * coef / (840.0 * dz);
	const double c1z = 672.0 * coef / (840.0 * dz);

	const int size_x = grid.SizeX() + MW * 2;
	const int size_y = grid.SizeY() + MW * 2;
	const int size_z = grid.SizeZ() + MW * 2;

	for (int iz = MW; iz < size_z - MW; ++iz) {
		const int idz_m4 = -4 * size_x * size_y;
		const int idz_m3 = -3 * size_x * size_y;
		const int idz_m2 = -2 * size_x * size_y;
		const int idz_m = -1 * size_x * size_y;
		const int idz_p = 1 * size_x * size_y;
		const int idz_p2 = 2 * size_x * size_y;
		const int idz_p3 = 3 * size_x * size_y;
		const int idz_p4 = 4 * size_x * size_y;

		for (int iy = MW; iy < size_y - MW; ++iy) {
			const int idy_m4 = -4 * size_x;
			const int idy_m3 = -3 * size_x;
			const int idy_m2 = -2 * size_x;
			const int idy_m = -1 * size_x;
			const int idy_p = 1 * size_x;
			const int idy_p2 = 2 * size_x;
			const int idy_p3 = 3 * size_x;
			const int idy_p4 = 4 * size_x;
#pragma ivdep
			for (int ix = MW; ix < size_x - MW; ++ix) {
				const int idx_m4 = -4;
				const int idx_m3 = -3;
				const int idx_m2 = -2;
				const int idx_m = -1;
				const int idx_p = 1;
				const int idx_p2 = 2;
				const int idx_p3 = 3;
				const int idx_p4 = 4;
				const int64_t i = ix + size_x * (iy + (size_y * iz));
				const int64_t io = (ix - MW) + Nx * (iy - MW + (Ny * (iz - MW)));

				auto dpsi_dz = c4z * (p[i + idz_p4] - p[i + idz_m4]) + c3z * (p[i + idz_p3] - p[i + idz_m3]) + c2z * (p[i + idz_p2] - p[i + idz_m2]) + c1z * (p[i + idz_p] - p[i + idz_m]);
				out[io] = dpsi_dz;
			}
		}
	}
}


template <class T>
void Gradient8th_halo(const GridRange& grid, T* out_x, T* out_y, T* out_z, const T* src, T* buf_with_halo, const double coef, const double dx, const double dy, const double dz) {
	const int size_x = grid.SizeX();
	const int size_y = grid.SizeY();
	const int size_z = grid.SizeZ();

	PasteHalo_f4(grid, buf_with_halo,src,
		[&](int ix, int iy, int iz) { return src[(ix + size_x * (iy + 1 + size_y * iz))]; },
		[&](int ix, int iy, int iz) { return src[(ix + size_x * (iy - 1 + size_y * iz))]; },
		[&](int ix, int iy, int iz) { return src[(ix + size_x * (iy + size_y * (iz + 1)))]; },
		[&](int ix, int iy, int iz) { return src[(ix + size_x * (iy + size_y * (iz - 1)))]; },
		[&](int ix, int iy, int iz) { return src[(ix + size_x * (iy + size_y * (iz + size_z)))]; },
		[&](int ix, int iy, int iz) { return src[(ix + size_x * (iy + size_y * (iz - size_z)))]; }	);

	const size_t local_size = grid.Size3D();
	Gradient8th_x_halo_f(grid, out_x, buf_with_halo, coef, dx, dy, dz);
	Gradient8th_y_halo_f(grid, out_y, buf_with_halo, coef, dx, dy, dz);
	Gradient8th_z_halo_f(grid, out_z, buf_with_halo, coef, dx, dy, dz);


}


/*
* wrapper of the fastest algorithm of 8th order difference equation
* 
*/
template <class T>
void Gradient8th(const GridRange& grid, T* out_x, T* out_y, T* out_z, const T* p, double coef, double dx, double dy, double dz) {
	
	const int margin_width = 4;
	const int over_size = (grid.SizeX() + 2 * margin_width) * (grid.SizeY() + 2 * margin_width) * (grid.SizeZ() + 2 * margin_width);
	std::vector<T> over_data(over_size);//margin data//

	Gradient8th_halo(grid, out_x, out_y, out_z, p, over_data.data(), coef, dx, dy, dz);
}


#ifdef USE_MPI


/*
* 1次元領域分割の場合のGradient
* 分割方向はz方向
*/
template <class T>
void Gradient8th_MPI_1D_overspace(const GridRangeMPI& grid, T* out_x, T* out_y, T* out_z, const T* src, T* buf_with_halo, double coef, double dx, double dy, double dz) {

	MakeOverspace_1D_mpi(grid, buf_with_halo, src);

	Gradient8th_x_halo_f(grid, out_x, buf_with_halo, coef, dx, dy, dz);
	Gradient8th_y_halo_f(grid, out_y, buf_with_halo, coef, dx, dy, dz);
	Gradient8th_z_halo_f(grid, out_z, buf_with_halo, coef, dx, dy, dz);

}


/*
* 3次元領域分割の場合のGradient
* 分割方向はz方向
*/
template <class T>
void Gradient8th_MPI_3D_overspace(const GridRangeMPI& grid, T* out_x, T* out_y, T* out_z, const T* src, T* buf_with_halo, double coef, double dx, double dy, double dz) {

   	MakeOverspace8th_3D_mpi(grid, buf_with_halo, src);

	Gradient8th_x_halo_f(grid, out_x, buf_with_halo, coef, dx, dy, dz);
	Gradient8th_y_halo_f(grid, out_y, buf_with_halo, coef, dx, dy, dz);
	Gradient8th_z_halo_f(grid, out_z, buf_with_halo, coef, dx, dy, dz);
}


template <class T>
void Gradient8th_ddm(const GridRangeMPI& grid, T* out_x, T* out_y, T* out_z, const T* src, double coef, double dx, double dy, double dz) {

	const int margin_width = 4;
	const int over_size = (grid.SizeX() + 2 * margin_width) * (grid.SizeY() + 2 * margin_width) * (grid.SizeZ() + 2 * margin_width);
	std::vector<T> buf_with_halo(over_size);//margin data//

#ifdef TEST_PRINT
    printf("buf_with_halo = %zd\n", buf_with_halo.size()); fflush(stdout);
#endif
	switch (grid.split_dimension) {
	case 1:
		Gradient8th_MPI_1D_overspace(grid, out_x, out_y, out_z, src, buf_with_halo.data(), coef, dx, dy, dz);
		break;
	case 2:
		Gradient8th_MPI_3D_overspace(grid, out_x, out_y, out_z, src, buf_with_halo.data(), coef, dx, dy, dz);
		break;
	case 3:		
		Gradient8th_MPI_3D_overspace(grid, out_x, out_y, out_z, src, buf_with_halo.data(), coef, dx, dy, dz);
		break;
	default:
		Gradient8th_halo(grid, out_x, out_y, out_z, src, buf_with_halo.data(), coef, dx, dy, dz);
		break;
	}
}



#endif
