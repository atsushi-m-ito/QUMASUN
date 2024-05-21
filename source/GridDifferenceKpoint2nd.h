#pragma once
#ifdef USE_MPI
#include <mpi.h>
#include "mpi_helper.h"
#endif
#include "GridRange.h"
#include <vector>

#include "GridDifference2nd.h"
#include "soacomplex.h"


//3次元デカルト座標グリッドの計算のためのループを司る機能//
//MPI並列時にはリダクションも提供する



template <class T, int KX, int KY, int KZ >
void LaplasianKpoint2nd_halo_f(const GridRange& grid, T* out, T* out2, const T* p, double coef, double sign1st, double dx, double dy, double dz, double gx, double gy, double gz) {

	constexpr int MW = 1;
	const int Nx = grid.SizeX();
	const int Ny = grid.SizeY();
	const int Nz = grid.SizeZ();


	const double c0 = -coef * (2.0 / (dx * dx) + 2.0 / (dy * dy) + 2.0 / (dz * dz)) - coef*(gx*gx + gy*gy+ gz*gz);
	const double c1x = coef / (dx * dx);
	const double c1y = coef / (dy * dy);
	const double c1z = coef / (dz * dz);

	const double coef2 = sign1st * coef * 2.0;
	const double a1x = coef2 * gx / (2.0 * dx);
	const double a1y = coef2 * gy / (2.0 * dy);
	const double a1z = coef2 * gz / (2.0 * dz);


	const int size_x = grid.SizeX() + MW * 2;
	const int size_y = grid.SizeY() + MW * 2;
	const int size_z = grid.SizeZ() + MW * 2;

	for (int iz = MW; iz < size_z - MW; ++iz) {
		//const int idz_m3 = -3 * size_x * size_y;
		//const int idz_m2 = -2 * size_x * size_y;
		const int idz_m = -1 * size_x * size_y;
		const int idz_p = 1 * size_x * size_y;
		//const int idz_p2 = 2 * size_x * size_y;
		//const int idz_p3 = 3 * size_x * size_y;

		for (int iy = MW; iy < size_y - MW; ++iy) {
			//const int idy_m3 = -3 * size_x;
			//const int idy_m2 = -2 * size_x;
			const int idy_m = -1 * size_x;
			const int idy_p = 1 * size_x;
			//const int idy_p2 = 2 * size_x;
			//const int idy_p3 = 3 * size_x;
#pragma ivdep
			for (int ix = MW; ix < size_x - MW; ++ix) {
				//const int idx_m3 = -3;
				//const int idx_m2 = -2;
				const int idx_m = -1;
				const int idx_p = 1;
				//const int idx_p2 = 2;
				//const int idx_p3 = 3;
				const int64_t i = ix + size_x * (iy + (size_y * iz));
				const int64_t io = (ix - MW) + Nx * (iy - MW + (Ny * (iz - MW)));

				//this is bad accuracy, and speed is slower than //
				const auto d2psidx2 = c0 * p[i]
					+ c1x * (p[i + idx_p] + p[i + idx_m])
					+ c1y * (p[i + idy_p] + p[i + idy_m])
					+ c1z * (p[i + idz_p] + p[i + idz_m]);
				out[io] += d2psidx2;


				auto dpsi_dx = 0.0;
				if constexpr (KX != 0) {
					dpsi_dx += a1x * (p[i + idx_p] - p[i + idx_m]);
				}
				if constexpr (KY != 0) {
					dpsi_dx += a1y * (p[i + idy_p] - p[i + idy_m]);
				}
				if constexpr (KZ != 0) {
					dpsi_dx += a1z * (p[i + idz_p] - p[i + idz_m]);

				}
				if constexpr ((KX != 0) || (KY != 0) || (KZ != 0)) {
					out2[io] += dpsi_dx;
				}
			}
		}
	}

}


template<int KX, int KY, int KZ>
void LaplasianKpoint2nd_halo(const GridRange& grid, SoAComplex& out, const SoAComplex& src, double* re_with_halo, double* im_with_halo, const double coef, const double dx, const double dy, const double dz, const double gx, const double gy, const double gz) {
	const int size_x = grid.SizeX();
	const int size_y = grid.SizeY();
	const int size_z = grid.SizeZ();

	PasteHalo_f1(grid, re_with_halo, src.re,
		[&](int ix, int iy, int iz) { return src.re[(ix + size_x * (iy + 1 + size_y * iz))]; },
		[&](int ix, int iy, int iz) { return src.re[(ix + size_x * (iy - 1 + size_y * iz))]; },
		[&](int ix, int iy, int iz) { return src.re[(ix + size_x * (iy + size_y * (iz + 1)))]; },
		[&](int ix, int iy, int iz) { return src.re[(ix + size_x * (iy + size_y * (iz - 1)))]; },
		[&](int ix, int iy, int iz) { return src.re[(ix + size_x * (iy + size_y * (iz + size_z)))]; },
		[&](int ix, int iy, int iz) { return src.re[(ix + size_x * (iy + size_y * (iz - size_z)))]; });

	LaplasianKpoint2nd_halo_f<double, KX, KY, KZ>(grid, out.re, out.im, re_with_halo, coef, +1.0, dx, dy, dz, gx, gy, gz);


	PasteHalo_f1(grid, im_with_halo, src.im,
		[&](int ix, int iy, int iz) { return src.im[(ix + size_x * (iy + 1 + size_y * iz))]; },
		[&](int ix, int iy, int iz) { return src.im[(ix + size_x * (iy - 1 + size_y * iz))]; },
		[&](int ix, int iy, int iz) { return src.im[(ix + size_x * (iy + size_y * (iz + 1)))]; },
		[&](int ix, int iy, int iz) { return src.im[(ix + size_x * (iy + size_y * (iz - 1)))]; },
		[&](int ix, int iy, int iz) { return src.im[(ix + size_x * (iy + size_y * (iz + size_z)))]; },
		[&](int ix, int iy, int iz) { return src.im[(ix + size_x * (iy + size_y * (iz - size_z)))]; });

	LaplasianKpoint2nd_halo_f<double, KX, KY, KZ>(grid, out.im, out.re, im_with_halo, coef, -1.0, dx, dy, dz, gx, gy, gz);


}


/*
* 1次元領域分割の場合のLaplasian
* 分割方向はz方向
*/
template <int KX, int KY, int KZ>
void LaplasianKpoint2nd_MPI_1D_overspace(const GridRangeMPI& grid, SoAComplex& out, const SoAComplex& src, double* re_with_halo, double* im_with_halo, double coef, double dx, double dy, double dz, const double gx, const double gy, const double gz) {

	//のり代送受信//
	constexpr int TAG = 0x1;
	const int size_x = grid.SizeX();
	const int size_y = grid.SizeY();
	const int size_xy = size_x * size_y;
	const int size_z = grid.SizeZ();
	const int margin_width = 1;
	static std::vector<double> margin_data(size_xy * margin_width * 2);//margin data//
	
#ifdef TEST_PRINT
	printf("Laplasian6th_MPI_1D_halo\n"); fflush(stdout);
#endif

	TransferMarginZ(grid, margin_data.data(), margin_data.data() + size_xy * margin_width, src.re, margin_width);
	
	PasteHalo_f1(grid, re_with_halo, src.re,
		[&](int ix, int iy, int iz) { return src.re[(ix + size_x * (iy + 1 + size_y * iz))]; },
		[&](int ix, int iy, int iz) { return src.re[(ix + size_x * (iy - 1 + size_y * iz))]; },
		[&](int ix, int iy, int iz) { return src.re[(ix + size_x * (iy + size_y * (iz + 1)))]; },
		[&](int ix, int iy, int iz) { return src.re[(ix + size_x * (iy + size_y * (iz - 1)))]; },
		[&](int ix, int iy, int iz) { return margin_data[(ix + size_x * (iy + size_y * (iz + margin_width)))]; },
		[&](int ix, int iy, int iz) { return margin_data[size_xy * margin_width + (ix + size_x * (iy + size_y * (iz - size_z)))]; });
	
	LaplasianKpoint2nd_halo_f<double, KX, KY, KZ>(grid, out.re, out.im, re_with_halo, coef, +1.0, dx, dy, dz, gx, gy, gz);



	TransferMarginZ(grid, margin_data.data(), margin_data.data() + size_xy * margin_width, src.im, margin_width);

	PasteHalo_f1(grid, im_with_halo, src.im,
		[&](int ix, int iy, int iz) { return src.im[(ix + size_x * (iy + 1 + size_y * iz))]; },
		[&](int ix, int iy, int iz) { return src.im[(ix + size_x * (iy - 1 + size_y * iz))]; },
		[&](int ix, int iy, int iz) { return src.im[(ix + size_x * (iy + size_y * (iz + 1)))]; },
		[&](int ix, int iy, int iz) { return src.im[(ix + size_x * (iy + size_y * (iz - 1)))]; },
		[&](int ix, int iy, int iz) { return margin_data[(ix + size_x * (iy + size_y * (iz + margin_width)))]; },
		[&](int ix, int iy, int iz) { return margin_data[size_xy * margin_width + (ix + size_x * (iy + size_y * (iz - size_z)))]; });

	LaplasianKpoint2nd_halo_f<double, KX, KY, KZ>(grid, out.im, out.re, im_with_halo, coef, -1.0, dx, dy, dz, gx, gy, gz);

}


/*
* 1次元領域分割の場合のLaplasian
* 分割方向はz方向
*/
template <int KX, int KY, int KZ>
void LaplasianKpoint2nd_MPI_3D_overspace(const GridRangeMPI& grid, SoAComplex& out, const SoAComplex& src, double* re_with_halo, double* im_with_halo, double coef, double dx, double dy, double dz, double gx, double gy, double gz) {

	//のり代送受信//
	constexpr int TAG = 0x1;
	const int size_x = grid.SizeX();
	const int size_y = grid.SizeY();
	const int size_xy = size_x * size_y;
	const int size_z = grid.SizeZ();
	const int size_xz = size_x * size_z;
	const int size_yz = size_y * size_z;
	const int margin_width = 1;
	static std::vector<double> margin_data((size_xy + size_xz + size_yz) * margin_width * 2);//margin data//
	
	auto* margin_x_left = margin_data.data();
	auto* margin_x_right = margin_x_left + size_yz * margin_width;
	auto* margin_y_left = margin_x_right + size_yz * margin_width;
	auto* margin_y_right = margin_y_left + size_xz * margin_width;
	auto* margin_z_left = margin_y_right + size_xz * margin_width;
	auto* margin_z_right = margin_z_left + size_xy * margin_width;

#ifdef TEST_PRINT
	int proc_id;
	MPI_Comm_rank(grid.mpi_comm, &proc_id);
	printf("[%d]test4-1\n", proc_id); fflush(stdout);
#endif
	TransferMarginX(grid, margin_x_left, margin_x_right, src.re, margin_width);
#ifdef TEST_PRINT
	printf("[%d]test4-2\n", proc_id); fflush(stdout);
#endif
	TransferMarginY(grid, margin_y_left, margin_y_right, src.re, margin_width);

#ifdef TEST_PRINT
	MPI_Barrier(grid.mpi_comm);
	printf("[%d]test4-3\n", proc_id); fflush(stdout);
#endif
	TransferMarginZ(grid, margin_z_left, margin_z_right, src.re, margin_width);

#ifdef TEST_PRINT
	MPI_Barrier(grid.mpi_comm);
	printf("[%d]test4-4\n", proc_id); fflush(stdout);
#endif

	PasteHalo_f1(grid, re_with_halo, src.re,
		[&](int ix, int iy, int iz) { return margin_x_left[ix + margin_width + margin_width * (iy + size_y * iz)]; },
		[&](int ix, int iy, int iz) { return margin_x_right[ix - size_x + margin_width * (iy + size_y * iz)]; },
		[&](int ix, int iy, int iz) { return margin_y_left[ix + size_x * (iy + margin_width + margin_width * iz)]; },
		[&](int ix, int iy, int iz) { return margin_y_right[ix + size_x * (iy - size_y + margin_width * (iz))]; },
		[&](int ix, int iy, int iz) { return margin_z_left[ix + size_x * (iy + size_y * (iz + margin_width))]; },
		[&](int ix, int iy, int iz) { return margin_z_right[ix + size_x * (iy + size_y * (iz - size_z))]; });

	LaplasianKpoint2nd_halo_f<double, KX, KY, KZ>(grid, out.re, out.im, re_with_halo, coef, +1.0, dx, dy, dz, gx, gy, gz);

	//Laplasian6th_halo_f<double>(grid, out.re, re_with_halo, coef,  dx, dy, dz);

#ifdef TEST_PRINT
	int proc_id;
	MPI_Comm_rank(grid.mpi_comm, &proc_id);
	printf("[%d]test4-1\n", proc_id); fflush(stdout);
#endif
	TransferMarginX(grid, margin_x_left, margin_x_right, src.im, margin_width);
#ifdef TEST_PRINT
	printf("[%d]test4-2\n", proc_id); fflush(stdout);
#endif
	TransferMarginY(grid, margin_y_left, margin_y_right, src.im, margin_width);

#ifdef TEST_PRINT
	MPI_Barrier(grid.mpi_comm);
	printf("[%d]test4-3\n", proc_id); fflush(stdout);
#endif
	TransferMarginZ(grid, margin_z_left, margin_z_right, src.im, margin_width);

#ifdef TEST_PRINT
	MPI_Barrier(grid.mpi_comm);
	printf("[%d]test4-4\n", proc_id); fflush(stdout);
#endif

	PasteHalo_f1(grid, im_with_halo, src.im,
		[&](int ix, int iy, int iz) { return margin_x_left[ix + margin_width + margin_width * (iy + size_y * iz)]; },
		[&](int ix, int iy, int iz) { return margin_x_right[ix - size_x + margin_width * (iy + size_y * iz)]; },
		[&](int ix, int iy, int iz) { return margin_y_left[ix + size_x * (iy + margin_width + margin_width * iz)]; },
		[&](int ix, int iy, int iz) { return margin_y_right[ix + size_x * (iy - size_y + margin_width * (iz))]; },
		[&](int ix, int iy, int iz) { return margin_z_left[ix + size_x * (iy + size_y * (iz + margin_width))]; },
		[&](int ix, int iy, int iz) { return margin_z_right[ix + size_x * (iy + size_y * (iz - size_z))]; });


	LaplasianKpoint2nd_halo_f<double, KX, KY, KZ>(grid, out.im, out.re, im_with_halo, coef, -1.0, dx, dy, dz, gx, gy, gz);
	//Laplasian6th_halo_f<double>(grid, out.im, im_with_halo, coef, dx, dy, dz);
#ifdef TEST_PRINT
	MPI_Barrier(grid.mpi_comm);
	printf("[%d]test4-5\n", proc_id); fflush(stdout);
#endif
}



template<int KX, int KY, int KZ>
void LaplasianKpoint2nd_t(const GridRangeMPI& grid, SoAComplex& out, const SoAComplex& src, double* re_with_halo, double* im_with_halo, 
	double coef, double dx, double dy, double dz, double gx, double gy, double gz) {
	

	switch (grid.split_dimension) {
	case 1:
		LaplasianKpoint2nd_MPI_1D_overspace<KX,KY,KZ>(grid, out, src, re_with_halo, im_with_halo, coef, dx, dy, dz, gx,gy,gz);
		break;
	case 2:
		LaplasianKpoint2nd_MPI_3D_overspace<KX, KY, KZ>(grid, out, src, re_with_halo, im_with_halo, coef, dx, dy, dz, gx, gy, gz);
		break;
	case 3:
		LaplasianKpoint2nd_MPI_3D_overspace<KX, KY, KZ>(grid, out, src, re_with_halo, im_with_halo, coef, dx, dy, dz, gx, gy, gz);
		break;
	default:
		LaplasianKpoint2nd_halo<KX, KY, KZ>(grid, out, src, re_with_halo, im_with_halo, coef, dx, dy, dz, gx, gy, gz);
		break;
	}
}

void LaplasianKpoint2nd_ddm(const GridRangeMPI& grid, SoAComplex& out, const SoAComplex& src, double coef, double dx, double dy, double dz, double gx, double gy, double gz) {
	const int margin_width = 1;
	const int over_size = (grid.SizeX() + 2 * margin_width) * (grid.SizeY() + 2 * margin_width) * (grid.SizeZ() + 2 * margin_width);
	static std::vector<double> re_with_halo(over_size);//margin data//
	static std::vector<double> im_with_halo(over_size);//margin data//

	constexpr double threshold = 1.0e-8;

	//int proc_id = GetProcessID(grid.mpi_comm);
	//printf("[%d] call LaplasianKpoint6th: gx=%f, %f, %f, %f\n", proc_id, gx, gy, gz, coef); fflush(stdout);

	if (fabs(gx)< threshold) {
		if (fabs(gy) < threshold) {
			if (fabs(gz) < threshold) {
				//printf("[%d] call LaplasianKpoint6th_t<0, 0, 0>\n", proc_id); fflush(stdout);
				LaplasianKpoint2nd_t<0, 0, 0>(grid, out, src, re_with_halo.data(), im_with_halo.data(), coef, dx, dy, dz, gx,gy,gz);
			} else {
				LaplasianKpoint2nd_t<0, 0, 1>(grid, out, src, re_with_halo.data(), im_with_halo.data(), coef, dx, dy, dz, gx, gy, gz);
			}
		} else {
			if (fabs(gz) < threshold) {
				LaplasianKpoint2nd_t<0, 1, 0>(grid, out, src, re_with_halo.data(), im_with_halo.data(), coef, dx, dy, dz, gx, gy, gz);
			} else {
				LaplasianKpoint2nd_t<0, 1, 1>(grid, out, src, re_with_halo.data(), im_with_halo.data(), coef, dx, dy, dz, gx, gy, gz);
			}
		}
	} else {
		if (fabs(gy) < threshold) {
			if (fabs(gz) < threshold) {
				LaplasianKpoint2nd_t<1, 0, 0>(grid, out, src, re_with_halo.data(), im_with_halo.data(), coef, dx, dy, dz, gx, gy, gz);
			} else {
				LaplasianKpoint2nd_t<1, 0, 1>(grid, out, src, re_with_halo.data(), im_with_halo.data(), coef, dx, dy, dz, gx, gy, gz);
			}
		} else {
			if (fabs(gz) < threshold) {
				LaplasianKpoint2nd_t<1, 1, 0>(grid, out, src, re_with_halo.data(), im_with_halo.data(), coef, dx, dy, dz, gx, gy, gz);
			} else {
				LaplasianKpoint2nd_t<1, 1, 1>(grid, out, src, re_with_halo.data(), im_with_halo.data(), coef, dx, dy, dz, gx, gy, gz);
			}
		}
	}

}

