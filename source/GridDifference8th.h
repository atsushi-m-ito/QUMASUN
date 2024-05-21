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



#ifndef PasteHalo_f4_defined

template <class T, class FUNC_XL, class FUNC_XR, class FUNC_YL, class FUNC_YR, class FUNC_ZL, class FUNC_ZR>
void PasteHalo_f4(const GridRange& grid, T* p, const T* src, FUNC_XL&& func_xl, FUNC_XR&& func_xr, FUNC_YL&& func_yl, FUNC_YR&& func_yr, FUNC_ZL&& func_zl, FUNC_ZR&& func_zr) {
	const int Nx = grid.SizeX();
	const int Ny = grid.SizeY();
	const int Nz = grid.SizeZ();
	constexpr int MW = 4;

	{
#pragma ivdep
		for (int iz = 0; iz < Nz; ++iz) {
			for (int iy = 0; iy < Ny; ++iy) {
				for (int ix = 0; ix < Nx; ++ix) {
					p[(ix + MW) + (Nx + MW*2) * (iy + MW + (Ny + MW*2) * (iz + MW))] = src[ix + Nx * (iy + Ny * iz)];
				}
			}
		}
#pragma ivdep
		for (int iz = 0; iz < Nz; ++iz) {
			for (int iy = 0; iy < Ny; ++iy) {
				p[0 + (Nx + MW*2) * (iy + MW + (Ny + MW*2) * (iz + MW))] = func_xl(-MW, iy, iz);
				p[1 + (Nx + MW * 2) * (iy + MW + (Ny + MW * 2) * (iz + MW))] = func_xl(-MW+1, iy, iz);
				p[2 + (Nx + MW * 2) * (iy + MW + (Ny + MW * 2) * (iz + MW))] = func_xl(-MW+2, iy, iz);
				p[3 + (Nx + MW * 2) * (iy + MW + (Ny + MW * 2) * (iz + MW))] = func_xl(-MW + 3, iy, iz);
				p[Nx + MW + 3 + (Nx + MW * 2) * (iy + MW + (Ny + MW * 2) * (iz + MW))] = func_xr(Nx + 3, iy, iz);
				p[Nx + MW+2 + (Nx + MW * 2) * (iy + MW + (Ny + MW * 2) * (iz + MW))] = func_xr(Nx + 2, iy, iz);
				p[Nx + MW+1 + (Nx + MW * 2) * (iy + MW + (Ny + MW * 2) * (iz + MW))] = func_xr(Nx + 1, iy, iz);
				p[Nx + MW + (Nx + MW * 2) * (iy + MW + (Ny + MW * 2) * (iz + MW))] = func_xr(Nx + 0, iy, iz);
			}
		}
#pragma ivdep
		for (int iz = 0; iz < Nz; ++iz) {
			for (int ix = 0; ix < Nx; ++ix) {
				p[(ix + MW) + (Nx + MW * 2) * (0 + (Ny + MW * 2) * (iz + MW))] = func_yl(ix, -MW, iz);
				p[(ix + MW) + (Nx + MW * 2) * (1 + (Ny + MW * 2) * (iz + MW))] = func_yl(ix, -MW+1, iz);
				p[(ix + MW) + (Nx + MW * 2) * (2 + (Ny + MW * 2) * (iz + MW))] = func_yl(ix, -MW+2, iz);
				p[(ix + MW) + (Nx + MW * 2) * (3 + (Ny + MW * 2) * (iz + MW))] = func_yl(ix, -MW + 3, iz);
				p[(ix + MW) + (Nx + MW * 2) * (Ny + MW + 3 + (Ny + MW * 2) * (iz + MW))] = func_yr(ix, Ny + 3, iz);
				p[(ix + MW) + (Nx + MW * 2) * (Ny + MW+2 + (Ny + MW * 2) * (iz + MW))] = func_yr(ix, Ny + 2, iz);
				p[(ix + MW) + (Nx + MW * 2) * (Ny + MW+1 + (Ny + MW * 2) * (iz + MW))] = func_yr(ix, Ny + 1, iz);
				p[(ix + MW) + (Nx + MW * 2) * (Ny + MW + (Ny + MW * 2) * (iz + MW))] = func_yr(ix, Ny + 0, iz);
			}
		}
#pragma ivdep
		for (int iy = 0; iy < Ny; ++iy) {
			for (int ix = 0; ix < Nx; ++ix) {
				p[(ix + MW) + (Nx + MW * 2) * (iy + MW + (Ny + MW * 2) * (0))] = func_zl(ix, iy, -MW);
				p[(ix + MW) + (Nx + MW * 2) * (iy + MW + (Ny + MW * 2) * (1))] = func_zl(ix, iy, -MW+1);
				p[(ix + MW) + (Nx + MW * 2) * (iy + MW + (Ny + MW * 2) * (2))] = func_zl(ix, iy, -MW+2);
				p[(ix + MW) + (Nx + MW * 2) * (iy + MW + (Ny + MW * 2) * (3))] = func_zl(ix, iy, -MW + 3);
			}
		}
		for (int iy = 0; iy < Ny; ++iy) {
			for (int ix = 0; ix < Nx; ++ix) {
				p[(ix + MW) + (Nx + MW * 2) * (iy + MW + (Ny + MW * 2) * (Nz + MW))] = func_zr(ix, iy, Nz + 0);
				p[(ix + MW) + (Nx + MW * 2) * (iy + MW + (Ny + MW * 2) * (Nz + MW+1))] = func_zr(ix, iy, Nz + 1);
				p[(ix + MW) + (Nx + MW * 2) * (iy + MW + (Ny + MW * 2) * (Nz + MW+2))] = func_zr(ix, iy, Nz + 2);
				p[(ix + MW) + (Nx + MW * 2) * (iy + MW + (Ny + MW * 2) * (Nz + MW + 3))] = func_zr(ix, iy, Nz + 3);
			}
		}
	}
}
#endif

template <class T>
void Laplasian8th_halo_f(const GridRange& grid, T* out, const T* p, double coef, double dx, double dy, double dz) {
	
	const int Nx = grid.SizeX();
	const int Ny = grid.SizeY();
	const int Nz = grid.SizeZ();
	constexpr int MW = 4;

	const double c0 = -14350.0/5040.0 * coef * (1.0 / ( dx * dx) + 1.0 / ( dy * dy) + 1.0 / ( dz * dz));
	const double c4x = -9.0 * coef / (5040.0 * dx * dx); 
	const double c3x = 128.0*coef / (5040.0 * dx * dx);
	const double c2x = -1008.0 * coef / (5040.0 * dx * dx);
	const double c1x = 8064.0 * coef / (5040.0 * dx * dx);
	const double c4y = -9.0 * coef / (5040.0 * dy * dy);
	const double c3y = 128.0 * coef / (5040.0 * dy * dy);
	const double c2y = -1008.0 * coef / (5040.0 * dy * dy);
	const double c1y = 8064.0 * coef / (5040.0 * dy * dy);
	const double c4z = -9.0 * coef / (5040.0 * dz * dz);
	const double c3z = 128.0 * coef / (5040.0 * dz * dz);
	const double c2z = -1008.0 * coef / (5040.0 * dz * dz);
	const double c1z = 8064.0 * coef / (5040.0 * dz * dz);

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


				auto d2psidx2 = c0 * p[i];
				d2psidx2 += c4x * (p[i + idx_p4] + p[i + idx_m4]) + c3x * (p[i + idx_p3] + p[i + idx_m3]) + c2x * (p[i + idx_p2] + p[i + idx_m2]) + c1x * (p[i + idx_p] + p[i + idx_m]);
				d2psidx2 += c4y * (p[i + idy_p4] + p[i + idy_m4]) + c3y * (p[i + idy_p3] + p[i + idy_m3]) + c2y * (p[i + idy_p2] + p[i + idy_m2]) + c1y * (p[i + idy_p] + p[i + idy_m]);
				d2psidx2 += c4z * (p[i + idz_p4] + p[i + idz_m4]) + c3z * (p[i + idz_p3] + p[i + idz_m3]) + c2z * (p[i + idz_p2] + p[i + idz_m2]) + c1z * (p[i + idz_p] + p[i + idz_m]);
				out[io] += d2psidx2;
			}
		}
	}
}


template <class T>
void Laplasian8th_halo(const GridRange& grid, T* out, const T* src, T* buf_with_halo, const double coef, const double dx, const double dy, const double dz) {
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


	Laplasian8th_halo_f(grid, out, buf_with_halo, coef, dx, dy, dz);


}


/*
* wrapper of the fastest algorithm of 8th order difference equation
* 
*/
template <class T>
void Laplasian8th(const GridRange& grid, T* out, const T* p, double coef, double dx, double dy, double dz) {
	
	const int margin_width = 4;
	const int over_size = (grid.SizeX() + 2 * margin_width) * (grid.SizeY() + 2 * margin_width) * (grid.SizeZ() + 2 * margin_width);
	std::vector<T> over_data(over_size);//margin data//

	Laplasian8th_halo(grid, out, p, over_data.data(), coef, dx, dy, dz);
}


#ifdef USE_MPI


/*
* 1次元領域分割の場合のLaplasian
* 分割方向はz方向
*/
template <class T>
void MakeOverspace_1D_mpi(const GridRangeMPI& grid, T* buf_with_halo, const T* src) {

	//のり代送受信//
	constexpr int TAG = 0x1;
	const int size_x = grid.SizeX();
	const int size_y = grid.SizeY();
	const int size_xy = size_x * size_y;
	const int size_z = grid.SizeZ();
	const int margin_width = 4;
	std::vector<T> margin_data(size_xy * margin_width * 2);//margin data//

	TransferMarginZ(grid, margin_data.data(), margin_data.data() + size_xy * margin_width, src, margin_width);

#ifdef TEST_PRINT
	printf("Laplasian6th_MPI_1D_halo\n"); fflush(stdout);
#endif
	PasteHalo_f4(grid, buf_with_halo, src,
		[&](int ix, int iy, int iz) { return src[(ix + size_x * (iy + 1 + size_y * iz))]; },
		[&](int ix, int iy, int iz) { return src[(ix + size_x * (iy - 1 + size_y * iz))]; },
		[&](int ix, int iy, int iz) { return src[(ix + size_x * (iy + size_y * (iz + 1)))]; },
		[&](int ix, int iy, int iz) { return src[(ix + size_x * (iy + size_y * (iz - 1)))]; },
		[&](int ix, int iy, int iz) { return margin_data[(ix + size_x * (iy + size_y * (iz + margin_width)))]; },
		[&](int ix, int iy, int iz) { return margin_data[size_xy * margin_width + (ix + size_x * (iy + size_y * (iz - size_z)))]; });

}

/*
* 1次元領域分割の場合のLaplasian
* 分割方向はz方向
*/
template <class T>
void Laplasian8th_MPI_1D_overspace(const GridRangeMPI& grid, T* out, const T* src, T* buf_with_halo, double coef, double dx, double dy, double dz) {

	MakeOverspace_1D_mpi(grid, buf_with_halo, src);

	Laplasian8th_halo_f(grid, out, buf_with_halo, coef, dx, dy, dz);

}


/*
* 3次元領域分割の場合のLaplasian
* 分割方向はz方向
*/
template <class T>
void MakeOverspace8th_3D_mpi(const GridRangeMPI& grid, T* buf_with_halo, const T* src) {

	//のり代送受信//
	constexpr int TAG = 0x1;
	const int size_x = grid.SizeX();
	const int size_y = grid.SizeY();
	const int size_xy = size_x * size_y;
	const int size_z = grid.SizeZ();
	const int size_xz = size_x * size_z;
	const int size_yz = size_y * size_z;
	const int margin_width = 4;
	std::vector<T> margin_data((size_xy + size_xz + size_yz) * margin_width * 2);//margin data//

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

	PasteHalo_f4(grid, buf_with_halo, src,

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

}

/*
* 3次元領域分割の場合のLaplasian
* 分割方向はz方向
*/
template <class T>
void Laplasian8th_MPI_3D_overspace(const GridRangeMPI & grid, T * out, const T * src, T * buf_with_halo, double coef, double dx, double dy, double dz) {

	MakeOverspace8th_3D_mpi(grid, buf_with_halo, src);

	Laplasian8th_halo_f(grid, out, buf_with_halo, coef, dx, dy, dz);

#ifdef TEST_PRINT
	MPI_Barrier(grid.mpi_comm);
	printf("[%d]test4-5\n", proc_id); fflush(stdout);
#endif
}


template <class T>
void Laplasian8th_ddm(const GridRangeMPI& grid, T* out, const T* src, double coef, double dx, double dy, double dz) {

	const int margin_width = 4;
	const int over_size = (grid.SizeX() + 2 * margin_width) * (grid.SizeY() + 2 * margin_width) * (grid.SizeZ() + 2 * margin_width);
	std::vector<T> buf_with_halo(over_size);//margin data//

#ifdef TEST_PRINT
	printf("L4 split_dimension = %d\n", grid.split_dimension);
#endif
	switch (grid.split_dimension) {
	case 1:
		Laplasian8th_MPI_1D_overspace(grid, out, src, buf_with_halo.data(), coef, dx, dy, dz);
		break;
	case 2:
		Laplasian8th_MPI_3D_overspace(grid, out, src, buf_with_halo.data(), coef, dx, dy, dz);
		break;
	case 3:		
		Laplasian8th_MPI_3D_overspace(grid, out, src, buf_with_halo.data(), coef, dx, dy, dz);
		break;
	default:
		Laplasian8th_halo(grid, out, src, buf_with_halo.data(), coef, dx, dy, dz);
		break;
	}
}



#endif
