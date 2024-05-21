#pragma once
#ifdef USE_MPI
#include <mpi.h>
#include "mpi_helper.h"
#endif
#include "GridRange.h"
#include "GridSubgrid.h"
#include <vector>



template <class T>
void Laplasian2nd_1(const GridRange& grid, T* out, const T* p, double coef, double dx, double dy, double dz) {
	const int size_x = grid.SizeX();
	const int size_y = grid.SizeY();
	const int size_z = grid.SizeZ();

	for (int iz = 0; iz < size_z; ++iz) {
		const int idz_m = ((iz == 0) ? size_z - 1 : -1) * size_x * size_y;
		const int idz_p = ((iz == size_z - 1) ? 1 - size_z : 1) * size_x * size_y;
		for (int iy = 0; iy < size_y; ++iy) {
			const int idy_m = ((iy == 0) ? size_y - 1 : -1) * size_x;
			const int idy_p = ((iy == size_y - 1) ? 1 - size_y : 1) * size_x;
			for (int ix = 0; ix < size_x; ++ix) {
				const int idx_m = ((ix == 0) ? size_x - 1 : -1);
				const int idx_p = ((ix == size_x - 1) ? 1 - size_x : 1);
				const int64_t i = ix + size_x * (iy + (size_y * iz));

				const auto d2psidx2 = (p[i + idx_p] + p[i + idx_m] - 2.0 * p[i]) / (dx * dx)
					+ (p[i + idy_p] + p[i + idy_m] - 2.0 * p[i]) / (dy * dy)
					+ (p[i + idz_p] + p[i + idz_m] - 2.0 * p[i]) / (dz * dz);
				out[i] += coef*d2psidx2;

			}
		}
	}
}


template <class T>
void Laplasian2nd_2(const GridRange& grid, T* out, const T* p, double coef, double dx, double dy, double dz) {
	const int size_x = grid.SizeX();
	const int size_y = grid.SizeY();
	const int size_z = grid.SizeZ();

	for (int iz = 0; iz < size_z; ++iz) {
		const int idz_m = ((iz == 0) ? size_z - 1 : -1) * size_x * size_y;
		const int idz_p = ((iz == size_z - 1) ? 1 - size_z : 1) * size_x * size_y;
		for (int iy = 0; iy < size_y; ++iy) {
			const int idy_m = ((iy == 0) ? size_y - 1 : -1) * size_x;
			const int idy_p = ((iy == size_y - 1) ? 1 - size_y : 1) * size_x;
			for (int ix = 1; ix < size_x-1; ++ix) {
				const int64_t i = ix + size_x * (iy + (size_y * iz));

				const auto d2psidx2 = (p[i + 1] + p[i - 1] - 2.0 * p[i]) / (dx * dx)
					+ (p[i + idy_p] + p[i + idy_m] - 2.0 * p[i]) / (dy * dy)
					+ (p[i + idz_p] + p[i + idz_m] - 2.0 * p[i]) / (dz * dz);
				out[i] += coef * d2psidx2;

			}
			{
				const int ix = 0;
				const int64_t i = ix + size_x * (iy + (size_y * iz));

				const auto d2psidx2 = (p[i + 1] + p[i - 1 + size_x] - 2.0 * p[i]) / (dx * dx)
					+ (p[i + idy_p] + p[i + idy_m] - 2.0 * p[i]) / (dy * dy)
					+ (p[i + idz_p] + p[i + idz_m] - 2.0 * p[i]) / (dz * dz);
				out[i] += coef * d2psidx2;
			}
			{
				const int ix = size_x - 1;
				const int64_t i = ix + size_x * (iy + (size_y * iz));

				const auto d2psidx2 = (p[i + 1 - size_x] + p[i - 1] - 2.0 * p[i]) / (dx * dx)
					+ (p[i + idy_p] + p[i + idy_m] - 2.0 * p[i]) / (dy * dy)
					+ (p[i + idz_p] + p[i + idz_m] - 2.0 * p[i]) / (dz * dz);
				out[i] += coef * d2psidx2;
			}
		}
	}
}


template <class T>
void Laplasian2nd_3(const GridRange & grid, T * out, const T * p, double coef, double dx, double dy, double dz) {
	const int size_x = grid.SizeX();
	const int size_y = grid.SizeY();
	const int size_z = grid.SizeZ();

	for (int iz = 0; iz < size_z; ++iz) {
		//const int idz_m = ((iz == 0) ? size_z - 1 : -1) * size_x * size_y;
		const int idz_p = ((iz == size_z - 1) ? 1 - size_z : 1) * size_x * size_y;
		for (int iy = 0; iy < size_y; ++iy) {
			const int idy_m = ((iy == 0) ? size_y - 1 : -1) * size_x;
			const int idy_p = ((iy == size_y - 1) ? 1 - size_y : 1) * size_x;
			for (int ix = 0; ix < size_x; ++ix) {
				const int idx_m = ((ix == 0) ? size_x - 1 : -1);
				const int idx_p = ((ix == size_x - 1) ? 1 - size_x : 1);
				const int64_t i = ix + size_x * (iy + (size_y * iz));

				const auto d2psidx2 = (p[i + idx_p] + p[i + idx_m] - 2.0 * p[i]) / (dx * dx)
					+ (p[i + idy_p] + p[i + idy_m] - 2.0 * p[i]) / (dy * dy)
					+ (p[i + idz_p] - 2.0 * p[i]) / (dz * dz);
				out[i] += coef * d2psidx2;

				const auto d2psidx2_up = (p[i]) / (dz * dz);
				out[i + idz_p] += coef * d2psidx2_up;
			}
		}
	}
}



template <class T>
void Laplasian2nd_4(const GridRange& grid, T* out, const T* p, double coef, double dx, double dy, double dz) {
	const int size_x = grid.SizeX();
	const int size_y = grid.SizeY();
	const int size_z = grid.SizeZ();

	for (int iz = 0; iz < size_z; ++iz) {
		const int idz_m = ((iz == 0) ? size_z - 1 : -1) * size_x * size_y;
		const int idz_p = ((iz == size_z - 1) ? 1 - size_z : 1) * size_x * size_y;
		for (int iy = 1; iy < size_y-1; ++iy) {
			const int idy_m = - size_x;
			const int idy_p = size_x;
			for (int ix = 1; ix < size_x - 1; ++ix) {
				const int64_t i = ix + size_x * (iy + (size_y * iz));

				const auto d2psidx2 = (p[i + 1] + p[i - 1] - 2.0 * p[i]) / (dx * dx)
					+ (p[i + idy_p] + p[i + idy_m] - 2.0 * p[i]) / (dy * dy)
					+ (p[i + idz_p] + p[i + idz_m] - 2.0 * p[i]) / (dz * dz);
				out[i] += coef * d2psidx2;

			}
			{
				const int ix = 0;
				const int64_t i = ix + size_x * (iy + (size_y * iz));

				const auto d2psidx2 = (p[i + 1] + p[i - 1 + size_x] - 2.0 * p[i]) / (dx * dx)
					+ (p[i + idy_p] + p[i + idy_m] - 2.0 * p[i]) / (dy * dy)
					+ (p[i + idz_p] + p[i + idz_m] - 2.0 * p[i]) / (dz * dz);
				out[i] += coef * d2psidx2;
			}
			{
				const int ix = size_x - 1;
				const int64_t i = ix + size_x * (iy + (size_y * iz));

				const auto d2psidx2 = (p[i + 1 - size_x] + p[i - 1] - 2.0 * p[i]) / (dx * dx)
					+ (p[i + idy_p] + p[i + idy_m] - 2.0 * p[i]) / (dy * dy)
					+ (p[i + idz_p] + p[i + idz_m] - 2.0 * p[i]) / (dz * dz);
				out[i] += coef * d2psidx2;
			}
		}
		{
			const int iy = 0;
			const int idy_m = (size_y - 1) * size_x;
			const int idy_p = size_x;
			for (int ix = 1; ix < size_x - 1; ++ix) {
				const int64_t i = ix + size_x * (iy + (size_y * iz));

				const auto d2psidx2 = (p[i + 1] + p[i - 1] - 2.0 * p[i]) / (dx * dx)
					+ (p[i + idy_p] + p[i + idy_m] - 2.0 * p[i]) / (dy * dy)
					+ (p[i + idz_p] + p[i + idz_m] - 2.0 * p[i]) / (dz * dz);
				out[i] += coef * d2psidx2;

			}
			{
				const int ix = 0;
				const int64_t i = ix + size_x * (iy + (size_y * iz));

				const auto d2psidx2 = (p[i + 1] + p[i - 1 + size_x] - 2.0 * p[i]) / (dx * dx)
					+ (p[i + idy_p] + p[i + idy_m] - 2.0 * p[i]) / (dy * dy)
					+ (p[i + idz_p] + p[i + idz_m] - 2.0 * p[i]) / (dz * dz);
				out[i] += coef * d2psidx2;
			}
			{
				const int ix = size_x - 1;
				const int64_t i = ix + size_x * (iy + (size_y * iz));

				const auto d2psidx2 = (p[i + 1 - size_x] + p[i - 1] - 2.0 * p[i]) / (dx * dx)
					+ (p[i + idy_p] + p[i + idy_m] - 2.0 * p[i]) / (dy * dy)
					+ (p[i + idz_p] + p[i + idz_m] - 2.0 * p[i]) / (dz * dz);
				out[i] += coef * d2psidx2;
			}
		}
		{
			const int iy = size_y - 1;
			const int idy_m = - size_x;
			const int idy_p = (1 - size_y) * size_x;
			for (int ix = 1; ix < size_x - 1; ++ix) {
				const int64_t i = ix + size_x * (iy + (size_y * iz));

				const auto d2psidx2 = (p[i + 1] + p[i - 1] - 2.0 * p[i]) / (dx * dx)
					+ (p[i + idy_p] + p[i + idy_m] - 2.0 * p[i]) / (dy * dy)
					+ (p[i + idz_p] + p[i + idz_m] - 2.0 * p[i]) / (dz * dz);
				out[i] += coef * d2psidx2;

			}
			{
				const int ix = 0;
				const int64_t i = ix + size_x * (iy + (size_y * iz));

				const auto d2psidx2 = (p[i + 1] + p[i - 1 + size_x] - 2.0 * p[i]) / (dx * dx)
					+ (p[i + idy_p] + p[i + idy_m] - 2.0 * p[i]) / (dy * dy)
					+ (p[i + idz_p] + p[i + idz_m] - 2.0 * p[i]) / (dz * dz);
				out[i] += coef * d2psidx2;
			}
			{
				const int ix = size_x - 1;
				const int64_t i = ix + size_x * (iy + (size_y * iz));

				const auto d2psidx2 = (p[i + 1 - size_x] + p[i - 1] - 2.0 * p[i]) / (dx * dx)
					+ (p[i + idy_p] + p[i + idy_m] - 2.0 * p[i]) / (dy * dy)
					+ (p[i + idz_p] + p[i + idz_m] - 2.0 * p[i]) / (dz * dz);
				out[i] += coef * d2psidx2;
			}
		}
	}
}


template <class T>
void Laplasian2nd_5_X(int size_x, int size_y, int iy, int iz, T* out, const T* p, const T* m_y_left, const T* m_y_right, const T* m_z_left, const T* m_z_right, double coef, double dx, double dy, double dz) {
	for (int ix = 1; ix < size_x - 1; ++ix) {
		const int64_t i = ix + size_x * (iy + (size_y * iz));

		const auto d2psidx2 = (p[i + 1] + p[i - 1] - 2.0 * p[i]) / (dx * dx)
			+ (m_y_left[ix] + m_y_right[ix] - 2.0 * p[i]) / (dy * dy)
			+ (m_z_left[ix] + m_z_right[ix] - 2.0 * p[i]) / (dz * dz);
		out[i] += coef * d2psidx2;

	}
	{
		int ix = 0;
		const int64_t i = ix + size_x * (iy + (size_y * iz));

		const auto d2psidx2 = (p[i + 1] + p[i - 1 + size_x] - 2.0 * p[i]) / (dx * dx)
			+ (m_y_left[ix] + m_y_right[ix] - 2.0 * p[i]) / (dy * dy)
			+ (m_z_left[ix] + m_z_right[ix] - 2.0 * p[i]) / (dz * dz);
		out[i] += coef * d2psidx2;
	}
	{
		int ix = size_x - 1;
		const int64_t i = ix + size_x * (iy + (size_y * iz));

		const auto d2psidx2 = (p[i + 1 - size_x] + p[i - 1] - 2.0 * p[i]) / (dx * dx)
			+ (m_y_left[ix] + m_y_right[ix] - 2.0 * p[i]) / (dy * dy)
			+ (m_z_left[ix] + m_z_right[ix] - 2.0 * p[i]) / (dz * dz);
		out[i] += coef * d2psidx2;
	}

}

template <class T>
void Laplasian2nd_5(const GridRange& grid, T* out, const T* p, double coef, double dx, double dy, double dz) {
	const int size_x = grid.SizeX();
	const int size_y = grid.SizeY();
	const int size_z = grid.SizeZ();

	for (int iz = 0; iz < size_z; ++iz) {
		const auto* mz_l = (iz == 0) ? &p[size_x * size_y * (iz - 1 + size_z)] : &p[size_x * size_y * (iz - 1)];
		const auto* mz_r = (iz == size_z - 1) ? &p[size_x * size_y * (iz + 1 - size_z)] : &p[size_x * size_y * (iz + 1)];
		for (int iy = 1; iy < size_y - 1; ++iy) {
			const int idy_m = -size_x;
			const int idy_p = size_x;
			Laplasian2nd_5_X(size_x, size_y, iy,iz, out,p, &p[size_x * (iy-1)], &p[size_x * (iy + 1)],
				mz_l, mz_r, -1.0/2.0, dx,dy,dz);
		}
		{
			const int iy = 0;
			Laplasian2nd_5_X(size_x, size_y, iy, iz, out, p, &p[size_x * (iy - 1 + size_y)], &p[size_x * (iy + 1)],
				mz_l, mz_r, -1.0 / 2.0, dx, dy, dz);
		}
		{
			const int iy = size_y - 1;
			Laplasian2nd_5_X(size_x, size_y, iy, iz, out, p, &p[size_x * (iy - 1)], &p[size_x * (iy + 1 - size_y)],
				mz_l, mz_r, -1.0 / 2.0, dx, dy, dz);
		}
	}
}

template <class T>
void Laplasian2nd_6_X(int size_x, int size_y, int iy, int iz, T* out, const T* p, const T* m_y_left, const T* m_y_right, const T* m_z_left, const T* m_z_right, double coef, double dx, double dy, double dz) {
	const double c0 = -2.0 / (dx * dx) - 2.0 / (dy * dy) - 2.0 / (dz * dz);
	for (int ix = 1; ix < size_x - 1; ++ix) {
		const int64_t i = ix + size_x * (iy + (size_y * iz));

		auto d2psidx2 = c0 * p[i];
		d2psidx2 += (p[i + 1] + p[i - 1]) / (dx * dx);
		d2psidx2 += (m_y_left[ix] + m_y_right[ix]) / (dy * dy);
		d2psidx2 += (m_z_left[ix] + m_z_right[ix]) / (dz * dz);
		out[i] += coef * d2psidx2;

	}
	{
		int ix = 0;
		const int64_t i = ix + size_x * (iy + (size_y * iz));

		auto d2psidx2 = c0 * p[i];
		d2psidx2 += (p[i + 1] + p[i - 1 + size_x]) / (dx * dx);
		d2psidx2 += (m_y_left[ix] + m_y_right[ix]) / (dy * dy);
		d2psidx2 += (m_z_left[ix] + m_z_right[ix]) / (dz * dz);
		out[i] += coef * d2psidx2;
	}
	{
		int ix = size_x - 1;
		const int64_t i = ix + size_x * (iy + (size_y * iz));

		auto d2psidx2 = c0 * p[i]; 
		d2psidx2 += (p[i + 1 - size_x] + p[i - 1]) / (dx * dx);
		d2psidx2 += (m_y_left[ix] + m_y_right[ix]) / (dy * dy);
		d2psidx2 += (m_z_left[ix] + m_z_right[ix]) / (dz * dz);
		out[i] += coef * d2psidx2;
	}

}


template <class T>
void Laplasian2nd_6(const GridRange& grid, T* out, const T* p, double coef, double dx, double dy, double dz) {
	const int size_x = grid.SizeX();
	const int size_y = grid.SizeY();
	const int size_z = grid.SizeZ();

	for (int iz = 0; iz < size_z; ++iz) {
		const auto* mz_l = (iz == 0) ? &p[size_x * size_y * (iz - 1 + size_z)] : &p[size_x * size_y * (iz - 1)];
		const auto* mz_r = (iz == size_z - 1) ? &p[size_x * size_y * (iz + 1 - size_z)] : &p[size_x * size_y * (iz + 1)];
		for (int iy = 1; iy < size_y - 1; ++iy) {
			const int idy_m = -size_x;
			const int idy_p = size_x;
			Laplasian2nd_6_X(size_x, size_y, iy, iz, out, p, &p[size_x * (iy - 1)], &p[size_x * (iy + 1)],
				mz_l, mz_r, -1.0 / 2.0, dx, dy, dz);
		}
		{
			const int iy = 0;
			Laplasian2nd_6_X(size_x, size_y, iy, iz, out, p, &p[size_x * (iy - 1 + size_y)], &p[size_x * (iy + 1)],
				mz_l, mz_r, -1.0 / 2.0, dx, dy, dz);
		}
		{
			const int iy = size_y - 1;
			Laplasian2nd_6_X(size_x, size_y, iy, iz, out, p, &p[size_x * (iy - 1)], &p[size_x * (iy + 1 - size_y)],
				mz_l, mz_r, -1.0 / 2.0, dx, dy, dz);
		}
	}
}


template <int MY, int MZ, class T>
void Laplasian2nd_7_X(int size_x, int size_y, int iy, int iz, T* out, const T* p, const T* m_y_left, const T* m_y_right, const T* m_z_left, const T* m_z_right, double coef, double dx, double dy, double dz) {
	const double c0 = -2.0 / (dx * dx) - 2.0 / (dy * dy) - 2.0 / (dz * dz);
	const int idy = size_x;
	const int idz = size_x * size_y;

	auto d2dy2 = [&](int64_t i, int ix) {
		if constexpr (MY == 0) {
			return (p[i - idy] + p[i + idy]) / (dy * dy);
		} else if constexpr (MY == 1) {
			return (p[i - idy] + m_y_right[ix]) / (dy * dy);
		} else if constexpr (MY == -1) {
			return (m_y_left[ix] + p[i + idy]) / (dy * dy);
		}
	};

	auto d2dz2 = [&](int64_t i, int ix) {
		if constexpr (MZ == 0) {
			return (p[i - idz] + p[i + idz]) / (dz * dz);
		} else if constexpr (MZ == 1) {
			return (p[i - idz] + m_z_right[ix]) / (dz * dz);
		} else if constexpr (MZ == -1) {
			return  (m_z_left[ix] + p[i + idz]) / (dz * dz);
		}
	};

	for (int ix = 1; ix < size_x - 1; ++ix) {
		const int64_t i = ix + size_x * (iy + (size_y * iz));

		auto d2psidx2 = c0 * p[i];
		d2psidx2 += (p[i-1] + p[i + 1]) / (dx * dx);

		d2psidx2 += d2dy2(i, ix) + d2dz2(i, ix);
		/*
		if constexpr (MY == 0) {
			d2psidx2 += (p[i-idy] + p[i+ idy]) / (dy * dy);
		} else if constexpr (MY == 1) {
			d2psidx2 += (p[i - idy] + m_y_right[ix]) / (dy * dy);
		} else if constexpr (MY == -1) {
			d2psidx2 += (m_y_left[ix] + p[i + idy]) / (dy * dy);
		}

		if constexpr (MZ == 0) {
			d2psidx2 += (p[i - idz] + p[i + idz]) / (dz * dz);
		} else if constexpr (MZ == 1) {
			d2psidx2 += (p[i - idz] + m_z_right[ix]) / (dz * dz);
		} else if constexpr (MZ == -1) {
			d2psidx2 += (m_z_left[ix] + p[i + idz]) / (dz * dz);
		}
		*/
		out[i] += coef * d2psidx2;

	}
	{
		int ix = 0;
		const int64_t i = ix + size_x * (iy + (size_y * iz));

		auto d2psidx2 = c0 * p[i];
		d2psidx2 += (p[i + 1] + p[i - 1 + size_x]) / (dx * dx);

		d2psidx2 += d2dy2(i, ix) + d2dz2(i, ix);
		/*
		if constexpr (MY == 0) {
			d2psidx2 += (p[i - idy] + p[i + idy]) / (dy * dy);
		} else if constexpr (MY == 1) {
			d2psidx2 += (p[i - idy] + m_y_right[ix]) / (dy * dy);
		} else if constexpr (MY == -1) {
			d2psidx2 += (m_y_left[ix] + p[i + idy]) / (dy * dy);
		}

		if constexpr (MZ == 0) {
			d2psidx2 += (p[i - idz] + p[i + idz]) / (dz * dz);
		} else if constexpr (MZ == 1) {
			d2psidx2 += (p[i - idz] + m_z_right[ix]) / (dz * dz);
		} else if constexpr (MZ == -1) {
			d2psidx2 += (m_z_left[ix] + p[i + idz]) / (dz * dz);
		}
		*/
		out[i] += coef * d2psidx2;

	}
	{
		int ix = size_x - 1;
		const int64_t i = ix + size_x * (iy + (size_y * iz));

		auto d2psidx2 = c0 * p[i];
		d2psidx2 += (p[i + 1 - size_x] + p[i - 1]) / (dx * dx);

		d2psidx2 += d2dy2(i, ix) + d2dz2(i, ix);
		/*
		if constexpr (MY == 0) {
			d2psidx2 += (p[i - idy] + p[i + idy]) / (dy * dy);
		} else if constexpr (MY == 1) {
			d2psidx2 += (p[i - idy] + m_y_right[ix]) / (dy * dy);
		} else if constexpr (MY == -1) {
			d2psidx2 += (m_y_left[ix] + p[i + idy]) / (dy * dy);
		}

		if constexpr (MZ == 0) {
			d2psidx2 += (p[i - idz] + p[i + idz]) / (dz * dz);
		} else if constexpr (MZ == 1) {
			d2psidx2 += (p[i - idz] + m_z_right[ix]) / (dz * dz);
		} else if constexpr (MZ == -1) {
			d2psidx2 += (m_z_left[ix] + p[i + idz]) / (dz * dz);
		}
		*/
		out[i] += coef * d2psidx2;

	}

}



template <class T>
void Laplasian2nd_7(const GridRange& grid, T* out, const T* p, double coef, double dx, double dy, double dz) {
	const int size_x = grid.SizeX();
	const int size_y = grid.SizeY();
	const int size_z = grid.SizeZ();


	for (int iz = 1; iz < size_z-1; ++iz) {
		for (int iy = 1; iy < size_y - 1; ++iy) {
			Laplasian2nd_7_X<0, 0>(size_x, size_y, iy, iz, out, p, 
				&p[size_x * (iy - 1 + size_y + size_y * iz)], &p[size_x * (iy + 1 - size_y + size_y * iz)],
				&p[size_x * (iy + size_y * (iz - 1 + size_z))], &p[size_x * (iy + size_y * (iz + 1 - size_z))], -1.0 / 2.0, dx, dy, dz);
		}
		{
			const int iy = 0;
			Laplasian2nd_7_X<-1, 0>(size_x, size_y, iy, iz, out, p, 
				&p[size_x * (iy - 1 + size_y + size_y * iz)], &p[size_x * (iy + 1 - size_y + size_y * iz)],
				&p[size_x * (iy + size_y * (iz - 1 + size_z))], &p[size_x * (iy + size_y * (iz + 1 - size_z))], -1.0 / 2.0, dx, dy, dz);
		}
		{
			const int iy = size_y - 1;
			Laplasian2nd_7_X<1, 0>(size_x, size_y, iy, iz, out, p, 
				&p[size_x * (iy - 1 + size_y + size_y * iz)], &p[size_x * (iy + 1 - size_y + size_y * iz)],
				&p[size_x * (iy + size_y * (iz - 1 + size_z))], &p[size_x * (iy + size_y * (iz + 1 - size_z))], -1.0 / 2.0, dx, dy, dz);
		}

	}
	{
		int iz = 0;
		for (int iy = 1; iy < size_y - 1; ++iy) {
			Laplasian2nd_7_X<0, -1>(size_x, size_y, iy, iz, out, p, 
				&p[size_x * (iy - 1 + size_y + size_y * iz)], &p[size_x * (iy + 1 - size_y + size_y * iz)],
				&p[size_x * (iy + size_y * (iz - 1 + size_z))], &p[size_x * (iy + size_y * (iz + 1 - size_z))], -1.0 / 2.0, dx, dy, dz);
		}
		{
			const int iy = 0;
			Laplasian2nd_7_X<-1, -1>(size_x, size_y, iy, iz, out, p, 
				&p[size_x * (iy - 1 + size_y + size_y * iz)], &p[size_x * (iy + 1 - size_y + size_y * iz)],
				&p[size_x * (iy + size_y * (iz - 1 + size_z))], &p[size_x * (iy + size_y * (iz + 1 - size_z))], -1.0 / 2.0, dx, dy, dz);
		}
		{
			const int iy = size_y - 1;
			Laplasian2nd_7_X<1, -1>(size_x, size_y, iy, iz, out, p, 
				&p[size_x * (iy - 1 + size_y + size_y * iz)], &p[size_x * (iy + 1 - size_y + size_y * iz)],
				&p[size_x * (iy + size_y * (iz - 1 + size_z))], &p[size_x * (iy + size_y * (iz + 1 - size_z))], -1.0 / 2.0, dx, dy, dz);
		}

	}
	{
		int iz = size_z-1;
		for (int iy = 1; iy < size_y - 1; ++iy) {
			Laplasian2nd_7_X<0, 1>(size_x, size_y, iy, iz, out, p, 
				&p[size_x * (iy - 1 + size_y + size_y * iz)], &p[size_x * (iy + 1 - size_y + size_y * iz)],
				&p[size_x * (iy + size_y * (iz - 1 + size_z))], &p[size_x * (iy + size_y * (iz + 1 - size_z))], -1.0 / 2.0, dx, dy, dz);
		}
		{
			const int iy = 0;
			Laplasian2nd_7_X<-1, 1>(size_x, size_y, iy, iz, out, p, 
				&p[size_x * (iy - 1 + size_y + size_y * iz)], &p[size_x * (iy + 1 - size_y + size_y * iz)],
				&p[size_x * (iy + size_y * (iz - 1 + size_z))], &p[size_x * (iy + size_y * (iz + 1 - size_z))], -1.0 / 2.0, dx, dy, dz);
		}
		{
			const int iy = size_y - 1;
			Laplasian2nd_7_X<1, 1>(size_x, size_y, iy, iz, out, p, 
				&p[size_x * (iy - 1 + size_y + size_y * iz)], &p[size_x * (iy + 1 - size_y + size_y * iz)],
				&p[size_x * (iy + size_y * (iz - 1 + size_z))], &p[size_x * (iy + size_y * (iz + 1 - size_z))], -1.0 / 2.0, dx, dy, dz);
		}

	}
}


template <int MZ, class T>
void Laplasian2nd_8_Y(int size_x, int size_y, int size_z, int iz, T* out, const T* p, double coef, double dx, double dy, double dz) {
	for (int iy = 1; iy < size_y - 1; ++iy) {
		Laplasian2nd_7_X<0, MZ>(size_x, size_y, iy, iz, out, p,
			&p[size_x * (iy - 1 + size_y + size_y * iz)], & p[size_x * (iy + 1 - size_y + size_y * iz)],
			& p[size_x * (iy + size_y * (iz - 1 + size_z))], & p[size_x * (iy + size_y * (iz + 1 - size_z))], -1.0 / 2.0, dx, dy, dz);
	}
	{
		const int iy = 0;
		Laplasian2nd_7_X<-1, MZ>(size_x, size_y, iy, iz, out, p,
			&p[size_x * (iy - 1 + size_y + size_y * iz)], &p[size_x * (iy + 1 - size_y + size_y * iz)],
			&p[size_x * (iy + size_y * (iz - 1 + size_z))], &p[size_x * (iy + size_y * (iz + 1 - size_z))], -1.0 / 2.0, dx, dy, dz);
	}
	{
		const int iy = size_y - 1;
		Laplasian2nd_7_X<1, MZ>(size_x, size_y, iy, iz, out, p,
			&p[size_x * (iy - 1 + size_y + size_y * iz)], &p[size_x * (iy + 1 - size_y + size_y * iz)],
			&p[size_x * (iy + size_y * (iz - 1 + size_z))], &p[size_x * (iy + size_y * (iz + 1 - size_z))], -1.0 / 2.0, dx, dy, dz);
	}
}


template <class T>
void Laplasian2nd_8(const GridRange& grid, T* out, const T* p, double coef, double dx, double dy, double dz) {
	const int size_x = grid.SizeX();
	const int size_y = grid.SizeY();
	const int size_z = grid.SizeZ();


	for (int iz = 1; iz < size_z - 1; ++iz) {
		Laplasian2nd_8_Y<0>(size_x, size_y, size_z, iz, out, p, coef, dx, dy, dz);
	}
	{
		int iz = 0;
		Laplasian2nd_8_Y<-1>(size_x, size_y, size_z, iz, out, p, coef, dx, dy, dz);

	}
	{
		int iz = size_z - 1;
		Laplasian2nd_8_Y<1>(size_x, size_y, size_z, iz, out, p, coef, dx, dy, dz);

	}
}



template <int MY, int MZ, class T, class FUNC_YL, class FUNC_YR, class FUNC_ZL, class FUNC_ZR>
void Laplasian2nd_9_X(int size_x, int size_y, int iy, int iz, T* out, const T* p, FUNC_YL&& f_y_left, FUNC_YR&& f_y_right, FUNC_ZL&& f_z_left, FUNC_ZR&& f_z_right, double coef, double dx, double dy, double dz) {
	const double c0 = -2.0 / (dx * dx) - 2.0 / (dy * dy) - 2.0 / (dz * dz);
	const int idy = size_x;
	const int idz = size_x * size_y;

	auto d2dy2 = [&](int64_t i, int ix) {
		if constexpr (MY == 0) {
			return (p[i - idy] + p[i + idy]) / (dy * dy);
		} else if constexpr (MY == -1) {
			return (f_y_left(ix,iy-1,iz) + p[i + idy]) / (dy * dy);
		} else if constexpr (MY == 1) {
			return (p[i - idy] + f_y_right(ix, iy + 1, iz)) / (dy * dy);
		}
		};

	auto d2dz2 = [&](int64_t i, int ix) {
		if constexpr (MZ == 0) {
			return (p[i - idz] + p[i + idz]) / (dz * dz);
		} else if constexpr (MZ == 1) {
			return (p[i - idz] + f_z_right(ix,iy,iz+1)) / (dz * dz);
		} else if constexpr (MZ == -1) {
			return  (f_z_left(ix,iy,iz-1) + p[i + idz]) / (dz * dz);
		}
		};

	for (int ix = 1; ix < size_x - 1; ++ix) {
		const int64_t i = ix + size_x * (iy + (size_y * iz));

		auto d2psidx2 = c0 * p[i];
		d2psidx2 += (p[i - 1] + p[i + 1]) / (dx * dx);

		d2psidx2 += d2dy2(i, ix) + d2dz2(i, ix);
		out[i] += coef * d2psidx2;

	}
	{
		int ix = 0;
		const int64_t i = ix + size_x * (iy + (size_y * iz));

		auto d2psidx2 = c0 * p[i];
		d2psidx2 += (p[i + 1] + p[i - 1 + size_x]) / (dx * dx);

		d2psidx2 += d2dy2(i, ix) + d2dz2(i, ix);
		
		out[i] += coef * d2psidx2;

	}
	{
		int ix = size_x - 1;
		const int64_t i = ix + size_x * (iy + (size_y * iz));

		auto d2psidx2 = c0 * p[i];
		d2psidx2 += (p[i + 1 - size_x] + p[i - 1]) / (dx * dx);

		d2psidx2 += d2dy2(i, ix) + d2dz2(i, ix);
		
		out[i] += coef * d2psidx2;

	}

}


template <int MZ, class T, class FUNC_YL, class FUNC_YR, class FUNC_ZL, class FUNC_ZR>
void Laplasian2nd_9_Y(int size_x, int size_y, int size_z, int iz, T* out, const T* p, FUNC_YL&& f_y_left, FUNC_YR&& f_y_right, FUNC_ZL&& f_z_left, FUNC_ZR&& f_z_right, double coef, double dx, double dy, double dz) {


	for (int iy = 1; iy < size_y - 1; ++iy) {
		Laplasian2nd_9_X<0, MZ>(size_x, size_y, iy, iz, out, p,
			f_y_left, f_y_right, f_z_left, f_z_right, -1.0 / 2.0, dx, dy, dz);
	}
	{
		const int iy = 0;
		Laplasian2nd_9_X<-1, MZ>(size_x, size_y, iy, iz, out, p,
			f_y_left, f_y_right, f_z_left, f_z_right, -1.0 / 2.0, dx, dy, dz);
	}
	{
		const int iy = size_y - 1;
		Laplasian2nd_9_X<1, MZ>(size_x, size_y, iy, iz, out, p,
			f_y_left, f_y_right, f_z_left, f_z_right, -1.0 / 2.0, dx, dy, dz);
	}
}


template <class T>
void Laplasian2nd_9(const GridRange& grid, T* out, const T* p, double coef, double dx, double dy, double dz) {
	const int size_x = grid.SizeX();
	const int size_y = grid.SizeY();
	const int size_z = grid.SizeZ();

	auto f_y_left = [&](int ix, int iy, int iz) {//periodic accessor when iy < 0//
		return p[ix + size_x * (iy + size_y * (iz + 1))];
		};
	auto f_y_right = [&](int ix, int iy, int iz) {//periodic accessor when iy >= Ny//
		return p[ix + size_x * (iy + size_y * (iz - 1))];
		};
	auto f_z_left = [&](int ix, int iy, int iz) {//periodic accessor when iz < 0//
		return p[ix + size_x * (iy + size_y * (iz + size_z))];
		};
	auto f_z_right = [&](int ix, int iy, int iz) {//periodic accessor when iz >= Nz//
		return p[ix + size_x * (iy + size_y * (iz - size_z))];
		};


	for (int iz = 1; iz < size_z - 1; ++iz) {
		Laplasian2nd_9_Y<0>(size_x, size_y, size_z, iz, out, p, f_y_left, f_y_right, f_z_left, f_z_right, coef, dx, dy, dz);
	}
	{
		int iz = 0;
		Laplasian2nd_9_Y<-1>(size_x, size_y, size_z, iz, out, p, f_y_left, f_y_right, f_z_left, f_z_right, coef, dx, dy, dz);

	}
	{
		int iz = size_z - 1;
		Laplasian2nd_9_Y<1>(size_x, size_y, size_z, iz, out, p, f_y_left, f_y_right, f_z_left, f_z_right, coef, dx, dy, dz);

	}
}

#if 0

template <class T>
void Laplasian2nd_overspace(const GridRange& grid, T* out, const T* src, T* p, const double coef, const double dx, const double dy, const double dz) {

	constexpr int MW = 1; //margin_width
	//constexpr int MW2 = MW; //margin_width

	const int Nx = grid.SizeX();
	const int Ny = grid.SizeY();
	const int Nz = grid.SizeZ();

	{//inside data//
#pragma ivdep
		for (int iz = 0; iz < Nz; ++iz) {
			for (int iy = 0; iy < Ny; ++iy) {
				for (int ix = 0; ix < Nx; ++ix) {
					p[(ix + MW) + (Nx + MW*2) * (iy + MW + (Ny + MW*2) * (iz + MW))] = src[ix + Nx * (iy + Ny * iz)];
				}
			}
		}
		//periodic condition: margin is copied from opposite side inside data //
#pragma ivdep
		for (int iz = 0; iz < Nz; ++iz) {
			for (int iy = 0; iy < Ny; ++iy) {
				p[0 + (Nx + MW * 2) * (iy + MW + (Ny + MW * 2) * (iz + MW))] = src[(Nx - MW) + Nx * (iy + Ny * iz)];
				p[Nx + MW + (Nx + MW * 2) * (iy + MW + (Ny + MW * 2) * (iz + MW))] = src[0 + Nx * (iy + Ny * iz)];
			}
		}
#pragma ivdep
		for (int iz = 0; iz < Nz; ++iz) {
			for (int ix = 0; ix < Nx; ++ix) {
				p[(ix + MW) + (Nx + MW * 2) * (0 + (Ny + MW * 2) * (iz + MW))] = src[ix + Nx * (Ny - MW + Ny * iz)];
				p[(ix + MW) + (Nx + MW * 2) * (Ny + MW + (Ny + MW * 2) * (iz + MW))] = src[ix + Nx * (0 + Ny * iz)];
			}
		}
#pragma ivdep
		for (int iy = 0; iy < Ny; ++iy) {
			for (int ix = 0; ix < Nx; ++ix) {
				p[(ix + MW) + (Nx + MW*2) * (iy + MW + (Ny + MW*2) * (0))] = src[ix + Nx * (iy + Ny * (Nz - MW))];
			}
		}
		for (int iy = 0; iy < Ny; ++iy) {
			for (int ix = 0; ix < Nx; ++ix) {
				p[(ix + MW) + (Nx + MW*2) * (iy + MW + (Ny + MW*2) * (Nz + MW))] = src[ix + Nx * (iy + Ny * 0)];				
			}
		}
	}


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

				const auto d2psidx2 = (p[i + idx_p] + p[i + idx_m] - 2.0 * p[i]) / (dx * dx)
					+ (p[i + idy_p] + p[i + idy_m] - 2.0 * p[i]) / (dy * dy)
					+ (p[i + idz_p] + p[i + idz_m] - 2.0 * p[i]) / (dz * dz);
				out[io] += coef * d2psidx2;
			}
		}
	}
}


template <class T, class FUNC_XL, class FUNC_XR, class FUNC_YL, class FUNC_YR, class FUNC_ZL, class FUNC_ZR>
void Laplasian2nd_bundle_X(int64_t Ns, int size_x, int size_y, int size_z, T* out, const T* p, FUNC_XL&& func_xl, FUNC_XR&& func_xr, FUNC_YL&& func_yl, FUNC_YR&& func_yr, FUNC_ZL&& func_zl, FUNC_ZR&& func_zr, double coef, double dx, double dy, double dz) {
	const double c0 = -2.0 / (dx * dx) - 2.0 / (dy * dy) - 2.0 / (dz * dz);

	auto GetHead = [&](int ix, int iy, int iz) {
		if (ix < 0) return func_xl(ix, iy, iz);
		if (ix >= size_x) return func_xr(ix, iy, iz);
		if (iy < 0) return func_yl(ix, iy, iz);
		if (iy >= size_y) return func_yr(ix, iy, iz);
		if (iz < 0) return func_zl(ix, iy, iz);
		if (iz >= size_z) return func_zr(ix, iy, iz);
		return p + Ns * (ix + size_x * (iy + (size_y * iz)));
		};


	for (int iz = 0; iz < size_z; ++iz) {
		for (int iy = 0; iy < size_y; ++iy) {
			for (int ix = 0; ix < size_x; ++ix) {
				const int64_t i = Ns * (ix + size_x * (iy + (size_y * iz)));
				const T* p0 = p + i;
				const T* xm1 = GetHead(ix - 1, iy, iz);
				const T* xp1 = GetHead(ix + 1, iy, iz);
				const T* zm1 = GetHead(ix, iy, iz - 1);
				const T* zp1 = GetHead(ix, iy, iz + 1);
				const T* ym1 = GetHead(ix, iy - 1, iz);
				const T* yp1 = GetHead(ix, iy + 1, iz);

				for (int n = 0; n < Ns; ++n) {
					auto d2psidx2 = c0 * p0[n];
					d2psidx2 += ((xm1[n] + xp1[n])) / (dx * dx);
					d2psidx2 += ((ym1[n] + yp1[n])) / (dy * dy);
					d2psidx2 += ((zm1[n] + zp1[n])) / (dz * dz);
					out[n + i] += coef * d2psidx2;
				}
			}
		}
	}

}

#else
//faster than above code//

template <class T>
void Laplasian2nd_overspace(const GridRange& grid, T* out, const T* src, T* p, const double coef, const double dx, const double dy, const double dz) {

	constexpr int MW = 1; //margin_width
	//constexpr int MW2 = MW; //margin_width

	const int Nx = grid.SizeX();
	const int Ny = grid.SizeY();
	const int Nz = grid.SizeZ();
	auto pow2 = [](T a) {return a * a; };
	const double c0 = -2.0*coef * (1.0 / pow2(dx) + 1.0 / pow2(dy) + 1.0 / pow2(dz));
	const double c1x = coef / pow2(dx);
	const double c1y = coef / pow2(dy);
	const double c1z = coef / pow2(dz);
	

	{//inside data//
#pragma ivdep
		for (int iz = 0; iz < Nz; ++iz) {
			for (int iy = 0; iy < Ny; ++iy) {
				for (int ix = 0; ix < Nx; ++ix) {
					p[(ix + MW) + (Nx + MW * 2) * (iy + MW + (Ny + MW * 2) * (iz + MW))] = src[ix + Nx * (iy + Ny * iz)];
				}
			}
		}
		//periodic condition: margin is copied from opposite side inside data //
#pragma ivdep
		for (int iz = 0; iz < Nz; ++iz) {
			for (int iy = 0; iy < Ny; ++iy) {
				p[0 + (Nx + MW * 2) * (iy + MW + (Ny + MW * 2) * (iz + MW))] = src[(Nx - MW) + Nx * (iy + Ny * iz)];
				p[Nx + MW + (Nx + MW * 2) * (iy + MW + (Ny + MW * 2) * (iz + MW))] = src[0 + Nx * (iy + Ny * iz)];
			}
		}
#pragma ivdep
		for (int iz = 0; iz < Nz; ++iz) {
			for (int ix = 0; ix < Nx; ++ix) {
				p[(ix + MW) + (Nx + MW * 2) * (0 + (Ny + MW * 2) * (iz + MW))] = src[ix + Nx * (Ny - MW + Ny * iz)];
				p[(ix + MW) + (Nx + MW * 2) * (Ny + MW + (Ny + MW * 2) * (iz + MW))] = src[ix + Nx * (0 + Ny * iz)];
			}
		}
#pragma ivdep
		for (int iy = 0; iy < Ny; ++iy) {
			for (int ix = 0; ix < Nx; ++ix) {
				p[(ix + MW) + (Nx + MW * 2) * (iy + MW + (Ny + MW * 2) * (0))] = src[ix + Nx * (iy + Ny * (Nz - MW))];
			}
		}
		for (int iy = 0; iy < Ny; ++iy) {
			for (int ix = 0; ix < Nx; ++ix) {
				p[(ix + MW) + (Nx + MW * 2) * (iy + MW + (Ny + MW * 2) * (Nz + MW))] = src[ix + Nx * (iy + Ny * 0)];
			}
		}
	}


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

				const auto d2psidx2 = c0 * p[i]
					+ c1x * (p[i + idx_p] + p[i + idx_m])
					+ c1y * (p[i + idy_p] + p[i + idy_m])
					+ c1z * (p[i + idz_p] + p[i + idz_m]);
				out[io] += d2psidx2;
			}
		}
	}
}

#if 1

template <class T, class FUNC_XL, class FUNC_XR, class FUNC_YL, class FUNC_YR, class FUNC_ZL, class FUNC_ZR>
void PasteHalo_f1(const GridRange& grid, T* p, const T* src, FUNC_XL&& func_xl, FUNC_XR&& func_xr, FUNC_YL&& func_yl, FUNC_YR&& func_yr, FUNC_ZL&& func_zl, FUNC_ZR&& func_zr) {

	constexpr int MW = 1; //margin_width
	//constexpr int MW2 = MW; //margin_width

	const int Nx = grid.SizeX();
	const int Ny = grid.SizeY();
	const int Nz = grid.SizeZ();
	

	{//inside data//
#pragma ivdep
		for (int iz = 0; iz < Nz; ++iz) {
			for (int iy = 0; iy < Ny; ++iy) {
				for (int ix = 0; ix < Nx; ++ix) {
					p[(ix + MW) + (Nx + MW * 2) * (iy + MW + (Ny + MW * 2) * (iz + MW))] = src[ix + Nx * (iy + Ny * iz)];
				}
			}
		}

		//periodic condition: margin is copied from opposite side inside data //
#pragma ivdep
		for (int iz = 0; iz < Nz; ++iz) {
			for (int iy = 0; iy < Ny; ++iy) {
				p[0 + (Nx + MW * 2) * (iy + MW + (Ny + MW * 2) * (iz + MW))] = func_xl(-1, iy, iz);
				p[Nx + MW + (Nx + MW * 2) * (iy + MW + (Ny + MW * 2) * (iz + MW))] = func_xr(Nx, iy, iz);
			}
		}

#pragma ivdep
		for (int iz = 0; iz < Nz; ++iz) {
			for (int ix = 0; ix < Nx; ++ix) {
				p[(ix + MW) + (Nx + MW * 2) * (0 + (Ny + MW * 2) * (iz + MW))] = func_yl(ix, -1, iz);
				p[(ix + MW) + (Nx + MW * 2) * (Ny + MW + (Ny + MW * 2) * (iz + MW))] = func_yr(ix, Ny, iz);
			}
		}

#pragma ivdep
		for (int iy = 0; iy < Ny; ++iy) {
			for (int ix = 0; ix < Nx; ++ix) {
				p[(ix + MW) + (Nx + MW * 2) * (iy + MW + (Ny + MW * 2) * (0))] = func_zl(ix, iy, -1);
			}
		}

		for (int iy = 0; iy < Ny; ++iy) {
			for (int ix = 0; ix < Nx; ++ix) {
				p[(ix + MW) + (Nx + MW * 2) * (iy + MW + (Ny + MW * 2) * (Nz + MW))] = func_zr(ix, iy, Nz + 0);
			}
		}
	}


}


template <class T>
void Laplasian2nd_halo_f(const GridRange& grid, T* out, const T* p, double coef, double dx, double dy, double dz) {

	constexpr int MW = 1; //margin_width
	//constexpr int MW2 = MW; //margin_width

	const int Nx = grid.SizeX();
	const int Ny = grid.SizeY();
	const int Nz = grid.SizeZ();
	
	const double c0 = -coef * (2.0 / (dx* dx) + 2.0 / (dy* dy) + 2.0 / (dz* dz));
	const double c1x = coef / (dx * dx);
	const double c1y = coef / (dy * dy);
	const double c1z = coef / (dz * dz);



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
#if 1
				const auto d2psidx2 =
					+c1x * (p[i + idx_p] + p[i + idx_m] - 2.0 * p[i])
					+ c1y * (p[i + idy_p] + p[i + idy_m] - 2.0 * p[i])
					+ c1z * (p[i + idz_p] + p[i + idz_m] - 2.0 * p[i]);
#else
				//this is bad accuracy, and speed is slower than //
				const auto d2psidx2 = c0 * p[i]
					+ c1x * (p[i + idx_p] + p[i + idx_m])
					+ c1y * (p[i + idy_p] + p[i + idy_m])
					+ c1z * (p[i + idz_p] + p[i + idz_m]);
#endif
				out[io] += d2psidx2;
			}
		}
	}
}


template <class T>
void Laplasian2nd_halo(const GridRange& grid, T* out, const T* src, T* buf_with_halo, const double coef, const double dx, const double dy, const double dz) {
	const int size_x = grid.SizeX();
	const int size_y = grid.SizeY();
	const int size_z = grid.SizeZ();

	PasteHalo_f1(grid, buf_with_halo, src,
		[&](int ix, int iy, int iz) { return src[(ix + size_x * (iy + 1 + size_y * iz))]; },
		[&](int ix, int iy, int iz) { return src[(ix + size_x * (iy - 1 + size_y * iz))]; },
		[&](int ix, int iy, int iz) { return src[(ix + size_x * (iy + size_y * (iz + 1)))]; },
		[&](int ix, int iy, int iz) { return src[(ix + size_x * (iy + size_y * (iz - 1)))]; },
		[&](int ix, int iy, int iz) { return src[(ix + size_x * (iy + size_y * (iz + size_z)))]; },
		[&](int ix, int iy, int iz) { return src[(ix + size_x * (iy + size_y * (iz - size_z)))]; }	);

	Laplasian2nd_halo_f(grid, out, buf_with_halo, coef, dx, dy, dz);

}


#else

template <class T, class FUNC_XL, class FUNC_XR, class FUNC_YL, class FUNC_YR, class FUNC_ZL, class FUNC_ZR>
void Laplasian2nd_overspace_f(const GridRange& grid, T* out, const T* src, T* p, FUNC_XL&& func_xl, FUNC_XR&& func_xr, FUNC_YL&& func_yl, FUNC_YR&& func_yr, FUNC_ZL&& func_zl, FUNC_ZR&& func_zr, double coef, double dx, double dy, double dz) {

	constexpr int MW = 1; //margin_width
	//constexpr int MW2 = MW; //margin_width

	const int Nx = grid.SizeX();
	const int Ny = grid.SizeY();
	const int Nz = grid.SizeZ();
	auto pow2 = [](T a) {return a * a; };
	const double c0 = -coef * (2.0 / pow2(dx) + 2.0 / pow2(dy) + 2.0 / pow2(dz));
	const double c1x = coef / pow2(dx);
	const double c1y = coef / pow2(dy);
	const double c1z = coef / pow2(dz);


	{//inside data//
#pragma ivdep
		for (int iz = 0; iz < Nz; ++iz) {
			for (int iy = 0; iy < Ny; ++iy) {
				for (int ix = 0; ix < Nx; ++ix) {
					p[(ix + MW) + (Nx + MW * 2) * (iy + MW + (Ny + MW * 2) * (iz + MW))] = src[ix + Nx * (iy + Ny * iz)];
				}
			}
		}

		//periodic condition: margin is copied from opposite side inside data //
#pragma ivdep
		for (int iz = 0; iz < Nz; ++iz) {
			for (int iy = 0; iy < Ny; ++iy) {
				p[0 + (Nx + MW * 2) * (iy + MW + (Ny + MW * 2) * (iz + MW))] = func_xl(-1, iy, iz);
				p[Nx + MW + (Nx + MW * 2) * (iy + MW + (Ny + MW * 2) * (iz + MW))] = func_xr(Nx, iy, iz);
			}
		}

#pragma ivdep
		for (int iz = 0; iz < Nz; ++iz) {
			for (int ix = 0; ix < Nx; ++ix) {
				p[(ix + MW) + (Nx + MW * 2) * (0 + (Ny + MW * 2) * (iz + MW))] = func_yl(ix, -1, iz);
				p[(ix + MW) + (Nx + MW * 2) * (Ny + MW + (Ny + MW * 2) * (iz + MW))] = func_yr(ix, Ny, iz);
			}
		}

#pragma ivdep
		for (int iy = 0; iy < Ny; ++iy) {
			for (int ix = 0; ix < Nx; ++ix) {
				p[(ix + MW) + (Nx + MW * 2) * (iy + MW + (Ny + MW * 2) * (0))] = func_zl(ix, iy, -1);
			}
		}

		for (int iy = 0; iy < Ny; ++iy) {
			for (int ix = 0; ix < Nx; ++ix) {
				p[(ix + MW) + (Nx + MW * 2) * (iy + MW + (Ny + MW * 2) * (Nz + MW))] = func_zr(ix, iy, Nz + 0);
			}
		}
	}


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
#if 1
				const auto d2psidx2 =
					+c1x * (p[i + idx_p] + p[i + idx_m] - 2.0 * p[i])
					+ c1y * (p[i + idy_p] + p[i + idy_m] - 2.0 * p[i])
					+ c1z * (p[i + idz_p] + p[i + idz_m] - 2.0 * p[i]);
#else
				//this is bad accuracy, and speed is slower than //
				const auto d2psidx2 = c0 * p[i]
					+ c1x * (p[i + idx_p] + p[i + idx_m])
					+ c1y * (p[i + idy_p] + p[i + idy_m])
					+ c1z * (p[i + idz_p] + p[i + idz_m]);
#endif
				out[io] += d2psidx2;
			}
		}
	}
}


template <class T>
void Laplasian2nd_overspace2(const GridRange& grid, T* out, const T* src, T* p, const double coef, const double dx, const double dy, const double dz) {
	const int size_x = grid.SizeX();
	const int size_y = grid.SizeY();
	const int size_z = grid.SizeZ();

	Laplasian2nd_overspace_f(grid, out, src, p,
		[&](int ix, int iy, int iz) { return src[(ix + size_x * (iy + 1 + size_y * iz))]; },
		[&](int ix, int iy, int iz) { return src[(ix + size_x * (iy - 1 + size_y * iz))]; },
		[&](int ix, int iy, int iz) { return src[(ix + size_x * (iy + size_y * (iz + 1)))]; },
		[&](int ix, int iy, int iz) { return src[(ix + size_x * (iy + size_y * (iz - 1)))]; },
		[&](int ix, int iy, int iz) { return src[(ix + size_x * (iy + size_y * (iz + size_z)))]; },
		[&](int ix, int iy, int iz) { return src[(ix + size_x * (iy + size_y * (iz - size_z)))]; },
		coef, dx, dy, dz);


}

#endif


template <class T>
void Laplasian2nd(const GridRange& grid, T* out, const T* src, double coef, double dx, double dy, double dz) {

	const int margin_width = 1;
	const int over_size = (grid.SizeX() + 2 * margin_width) * (grid.SizeY() + 2 * margin_width) * (grid.SizeZ() + 2 * margin_width);
	std::vector<T> over_data(over_size);//margin data//
	
	Laplasian2nd_halo(grid, out, src, over_data.data(), coef, dx, dy, dz);
}


template <class T, class FUNC_XL, class FUNC_XR, class FUNC_YL, class FUNC_YR, class FUNC_ZL, class FUNC_ZR>
void Laplasian2nd_bundle_X(int64_t Ns, int size_x, int size_y, int size_z, T* out, const T* p, FUNC_XL&& func_xl, FUNC_XR&& func_xr, FUNC_YL&& func_yl, FUNC_YR&& func_yr, FUNC_ZL&& func_zl, FUNC_ZR&& func_zr, double coef, double dx, double dy, double dz) {

	auto pow2 = [](T a) {return a * a; };
	const double c0 = -2.0 * coef * (1.0 / pow2(dx) + 1.0 / pow2(dy) + 1.0 / pow2(dz));
	const double c1x = coef / pow2(dx);
	const double c1y = coef / pow2(dy);
	const double c1z = coef / pow2(dz);

	auto GetHead = [&](int ix, int iy, int iz) {
		if (ix < 0) return func_xl(ix, iy, iz);
		if (ix >= size_x) return func_xr(ix, iy, iz);
		if (iy < 0) return func_yl(ix, iy, iz);
		if (iy >= size_y) return func_yr(ix, iy, iz);
		if (iz < 0) return func_zl(ix, iy, iz);
		if (iz >= size_z) return func_zr(ix, iy, iz);
		return p + Ns * (ix + size_x * (iy + (size_y * iz)));
		};


	for (int iz = 0; iz < size_z; ++iz) {
		for (int iy = 0; iy < size_y; ++iy) {
			for (int ix = 0; ix < size_x; ++ix) {
				const int64_t i = Ns * (ix + size_x * (iy + (size_y * iz)));
				const T* p0 = p + i;
				const T* xm1 = GetHead(ix - 1, iy, iz);
				const T* xp1 = GetHead(ix + 1, iy, iz);
				const T* zm1 = GetHead(ix, iy, iz - 1);
				const T* zp1 = GetHead(ix, iy, iz + 1);
				const T* ym1 = GetHead(ix, iy - 1, iz);
				const T* yp1 = GetHead(ix, iy + 1, iz);
#pragma ivdep
				for (int n = 0; n < Ns; ++n) {
					auto d2psidx2 = c0 * p0[n];
					d2psidx2 += c1x * (xm1[n] + xp1[n]);
					d2psidx2 += c1y * (ym1[n] + yp1[n]);
					d2psidx2 += c1z * (zm1[n] + zp1[n]);
					out[n + i] += d2psidx2;
				}
			}
		}
	}

}

#endif




//periodic//
template <class T>
void Laplasian2nd_bundle(const GridRange& grid, int64_t Ns, T* out, const T* p, const double coef, const double dx, const double dy, const double dz) {
	const int size_x = grid.SizeX();
	const int size_y = grid.SizeY();
	const int size_z = grid.SizeZ();


	Laplasian2nd_bundle_X<T>(Ns, size_x, size_y, size_z, out, p,
		[&](int ix, int iy, int iz) { return p + Ns * (ix + size_x * (iy + 1 + size_y * iz)); },
		[&](int ix, int iy, int iz) { return p + Ns * (ix + size_x * (iy - 1 + size_y * iz)); },
		[&](int ix, int iy, int iz) { return p + Ns * (ix + size_x * (iy + size_y * (iz + 1))); },
		[&](int ix, int iy, int iz) { return p + Ns * (ix + size_x * (iy + size_y * (iz - 1))); },
		[&](int ix, int iy, int iz) { return p + Ns * (ix + size_x * (iy + size_y * (iz + size_z))); },
		[&](int ix, int iy, int iz) { return p + Ns * (ix + size_x * (iy + size_y * (iz - size_z))); },
		coef, dx, dy, dz);
}



#ifdef USE_MPI

/*
* Laplasianの演算の内、X,Y方向のループの計算.
* X,Y方向は周期境界
* Z方向の隣のバッファもp自体.
* Please compare Laplasian2nd_calc_XY_margin
*/
template <class T>
void Laplasian2nd_calc_XY(const GridRangeMPI& grid, int iz, T* out, const T* p, double coef, double dx, double dy, double dz) {

	const int size_x = grid.SizeX();
	const int size_y = grid.SizeY();
	const int size_z = grid.SizeZ();

	const int idz_m = - size_x * size_y;
	const int idz_p = size_x * size_y;
	for (int iy = 0; iy < size_y; ++iy) {
		const int idy_m = ((iy == 0) ? size_y - 1 : -1) * size_x;
		const int idy_p = ((iy == size_y - 1) ? 1 - size_y : 1) * size_x;
		for (int ix = 1; ix < size_x - 1; ++ix) {
			const int64_t i = ix + size_x * (iy + (size_y * iz));

			const double d2psidx2 = (p[i + 1] + p[i - 1] - 2.0 * p[i]) / (dx * dx)
				+ (p[i + idy_p] + p[i + idy_m] - 2.0 * p[i]) / (dy * dy)
				+ (p[i + idz_p] + p[i + idz_m] - 2.0 * p[i]) / (dz * dz);
			out[i] += coef * d2psidx2;

		}
		{
			const int ix = 0;
			const int64_t i = ix + size_x * (iy + (size_y * iz));

			const double d2psidx2 = (p[i + 1] + p[i - 1 + size_x] - 2.0 * p[i]) / (dx * dx)
				+ (p[i + idy_p] + p[i + idy_m] - 2.0 * p[i]) / (dy * dy)
				+ (p[i + idz_p] + p[i + idz_m] - 2.0 * p[i]) / (dz * dz);
			out[i] += coef * d2psidx2;
		}
		{
			const int ix = size_x - 1;
			const int64_t i = ix + size_x * (iy + (size_y * iz));

			const double d2psidx2 = (p[i + 1 - size_x] + p[i - 1] - 2.0 * p[i]) / (dx * dx)
				+ (p[i + idy_p] + p[i + idy_m] - 2.0 * p[i]) / (dy * dy)
				+ (p[i + idz_p] + p[i + idz_m] - 2.0 * p[i]) / (dz * dz);
			out[i] += coef * d2psidx2;
		}
	}
}

/*
* Laplasianの演算の内、XY方向のループの計算
* X,Y方向は周期境界
* Z方向の隣のバッファを任意に指定できる
* pと同じでも良いし、marginでもよい。
*/
template <class T>
void Laplasian2nd_calc_XY_margin(const GridRangeMPI& grid, int iz, T* out, const T* p, const T* m_z_left, const T* m_z_right, double coef, double dx, double dy, double dz) {

	const int size_x = grid.SizeX();
	const int size_y = grid.SizeY();
	const int size_z = grid.SizeZ();

	for (int iy = 0; iy < size_y; ++iy) {
		const int idy_m = ((iy == 0) ? size_y - 1 : -1) * size_x;
		const int idy_p = ((iy == size_y - 1) ? 1 - size_y : 1) * size_x;
		for (int ix = 1; ix < size_x - 1; ++ix) {
			const int64_t i = ix + size_x * (iy + (size_y * iz));
			const int64_t ixy = ix + size_x * (iy);

			const double d2psidx2 = (p[i + 1] + p[i - 1] - 2.0 * p[i]) / (dx * dx)
				+ (p[i + idy_p] + p[i + idy_m] - 2.0 * p[i]) / (dy * dy)
				+ (m_z_left[ixy] + m_z_right[ixy] - 2.0 * p[i]) / (dz * dz);
			out[i] += coef * d2psidx2;

		}
		{
			const int ix = 0;
			const int64_t i = ix + size_x * (iy + (size_y * iz));
			const int64_t ixy = ix + size_x * (iy);

			const double d2psidx2 = (p[i + 1] + p[i - 1 + size_x] - 2.0 * p[i]) / (dx * dx)
				+ (p[i + idy_p] + p[i + idy_m] - 2.0 * p[i]) / (dy * dy)
				+ (m_z_left[ixy] + m_z_right[ixy] - 2.0 * p[i]) / (dz * dz);
			out[i] += coef * d2psidx2;
		}
		{
			const int ix = size_x - 1;
			const int64_t i = ix + size_x * (iy + (size_y * iz));
			const int64_t ixy = ix + size_x * (iy);

			const double d2psidx2 = (p[i + 1 - size_x] + p[i - 1] - 2.0 * p[i]) / (dx * dx)
				+ (p[i + idy_p] + p[i + idy_m] - 2.0 * p[i]) / (dy * dy)
				+ (m_z_left[ixy] + m_z_right[ixy] - 2.0 * p[i]) / (dz * dz);
			out[i] += coef * d2psidx2;
		}
	}
}


/*
* 1次元領域分割の場合のLaplasian
* 分割方向はz方向
*/
template <class T>
void Laplasian2nd_MPI_1D(const GridRangeMPI& grid, T* out, const T* p, double coef, double dx, double dy, double dz) {
	
	//のり代送受信//
	constexpr int TAG = 0x1;
	const int64_t size_xy = grid.SizeX() * grid.SizeY();
	const int size_z = grid.SizeZ();
	std::vector<double> margin_data(size_xy*2);//margin data//
	MPI_Status status;
	MPI_Sendrecv(p, size_xy, GetDatatype<T>(), grid.proc_z_left, TAG, margin_data.data() + size_xy * 1, size_xy, GetDatatype<T>(), grid.proc_z_right, TAG, grid.mpi_comm, &status);
	MPI_Sendrecv(p + size_xy * (size_z - 1), size_xy, GetDatatype<T>(), grid.proc_z_right, TAG, margin_data.data(), size_xy, GetDatatype<T>(), grid.proc_z_left, TAG, grid.mpi_comm, &status);

	for (int iz = 1; iz < size_z - 1; ++iz) {	
		Laplasian2nd_calc_XY(grid, iz, out, p, coef, dx, dy, dz);
	}

	{
		Laplasian2nd_calc_XY_margin(grid, 0, out, p, margin_data.data(), p + size_xy , coef, dx, dy, dz);
		Laplasian2nd_calc_XY_margin(grid, size_z-1, out, p, p + size_xy * (size_z-2), margin_data.data()+ size_xy * 1, coef, dx, dy, dz);
	}

}

template <class T>
void TransferMarginZ(const GridRangeMPI& grid, T* margin_z_left, T* margin_z_right, const T* src, int margin_width) {
	//のり代送受信//
	constexpr int TAG = 0x1;
	const int size_xy = grid.SizeX() * grid.SizeY();
	const int size_z = grid.SizeZ();
	const int send_size = size_xy * margin_width;
	
	MPI_Status status;
	MPI_Sendrecv(src, send_size, GetDatatype<T>(), grid.proc_z_left, TAG, margin_z_right, send_size, GetDatatype<T>(), grid.proc_z_right, TAG, grid.mpi_comm, &status);
	MPI_Sendrecv(src + size_xy * (size_z - margin_width), send_size, GetDatatype<T>(), grid.proc_z_right, TAG, margin_z_left, send_size, GetDatatype<T>(), grid.proc_z_left, TAG, grid.mpi_comm, &status);

}

//#define WITH_MPI_SENDRECV_REPLACE 

template <class T>
void TransferMarginX(const GridRangeMPI& grid, T* margin_x_left, T* margin_x_right, const T* src, int margin_width) {
	//のり代送受信//
	constexpr int TAG = 0x1;
	const int size_yz = grid.SizeY() * grid.SizeZ();
	const int size_x = grid.SizeX();
	const int send_size = size_yz * margin_width;

#ifdef WITH_MPI_SENDRECV_REPLACE
	T* send_buf_left = margin_x_right;
	T* send_buf_right = margin_x_left;
#else
	std::vector<T> send_buf(send_size * 2);
	T* send_buf_left = send_buf.data();
	T* send_buf_right = send_buf.data() + send_size;
#endif

	GridRange grid_src{ 0,0,0,grid.SizeX(),grid.SizeY() ,grid.SizeZ() };


	CutSubgrid(GridRange{ 0,0,0,margin_width,grid.SizeY() ,grid.SizeZ() }, send_buf_left, grid_src, src);
	CutSubgrid(GridRange{ size_x- margin_width,0,0,size_x,grid.SizeY() ,grid.SizeZ() }, send_buf_right, grid_src, src);

#ifdef WITH_MPI_SENDRECV_REPLACE
	MPI_Status status;
	MPI_Sendrecv_replace(margin_x_right, send_size, GetDatatype<T>(), grid.proc_x_left, TAG, grid.proc_x_right, TAG, grid.mpi_comm, &status);

	MPI_Sendrecv_replace(margin_x_left, send_size, GetDatatype<T>(), grid.proc_x_right, TAG, grid.proc_x_left, TAG, grid.mpi_comm, &status);
#else
	MPI_Status status;
	MPI_Sendrecv(send_buf_left, send_size, GetDatatype<T>(), grid.proc_x_left, TAG, margin_x_right, send_size, GetDatatype<T>(), grid.proc_x_right, TAG, grid.mpi_comm, &status);

	MPI_Sendrecv(send_buf_right, send_size, GetDatatype<T>(), grid.proc_x_right, TAG, margin_x_left, send_size, GetDatatype<T>(), grid.proc_x_left, TAG, grid.mpi_comm, &status);

#endif

}


template <class T>
void TransferMarginY(const GridRangeMPI& grid, T* margin_y_left, T* margin_y_right, const T* src, int margin_width) {
	//のり代送受信//
	constexpr int TAG = 0x1;
	const int size_xz = grid.SizeX() * grid.SizeZ();
	const int size_y = grid.SizeY();
	const int send_size = size_xz * margin_width;

#ifdef WITH_MPI_SENDRECV_REPLACE
	T* send_buf_left = margin_y_right;
	T* send_buf_right = margin_y_left;
#else
	std::vector<T> send_buf(send_size * 2);
	T* send_buf_left = send_buf.data();
	T* send_buf_right = send_buf.data() + send_size;
#endif
	
	GridRange grid_src{ 0,0,0,grid.SizeX(),grid.SizeY() ,grid.SizeZ() };
	CutSubgrid(GridRange{ 0,0,0,grid.SizeX(), margin_width ,grid.SizeZ() }, send_buf_left, grid_src, src);
	CutSubgrid(GridRange{ 0,size_y - margin_width,0,grid.SizeX(),size_y ,grid.SizeZ() }, send_buf_right, grid_src, src);


#ifdef WITH_MPI_SENDRECV_REPLACE
	MPI_Status status;
	MPI_Sendrecv_replace(margin_y_right, send_size, GetDatatype<T>(), grid.proc_y_left, TAG, grid.proc_y_right, TAG, grid.mpi_comm, &status);

	MPI_Sendrecv_replace(margin_y_left, send_size, GetDatatype<T>(), grid.proc_y_right, TAG, grid.proc_y_left, TAG, grid.mpi_comm, &status);
#else
	MPI_Status status;
	MPI_Sendrecv(send_buf_left, send_size, GetDatatype<T>(), grid.proc_y_left, TAG, margin_y_right, send_size, GetDatatype<T>(), grid.proc_y_right, TAG, grid.mpi_comm, &status);

	MPI_Sendrecv(send_buf_right, send_size, GetDatatype<T>(), grid.proc_y_right, TAG, margin_y_left, send_size, GetDatatype<T>(), grid.proc_y_left, TAG, grid.mpi_comm, &status);
#endif


}

/*
* 1次元領域分割の場合のLaplasian
* 分割方向はz方向
*/
template <class T>
void Laplasian2nd_MPI_1D_overspace(const GridRangeMPI& grid, T* out, const T* src, T* buf_with_halo, double coef, double dx, double dy, double dz) {

	//のり代送受信//
	constexpr int TAG = 0x1;
	const int size_x = grid.SizeX();
	const int size_y = grid.SizeY();
	const int size_xy = size_x * size_y;
	const int size_z = grid.SizeZ();
	const int margin_width = 1;
	std::vector<T> margin_data(size_xy * margin_width * 2);//margin data//
	
	TransferMarginZ(grid, margin_data.data(), margin_data.data() + size_xy * margin_width, src, margin_width);


	PasteHalo_f1(grid, buf_with_halo, src,
		[&](int ix, int iy, int iz) { return src[(ix + size_x * (iy + 1 + size_y * iz))]; },
		[&](int ix, int iy, int iz) { return src[(ix + size_x * (iy - 1 + size_y * iz))]; },
		[&](int ix, int iy, int iz) { return src[(ix + size_x * (iy + size_y * (iz + 1)))]; },
		[&](int ix, int iy, int iz) { return src[(ix + size_x * (iy + size_y * (iz - 1)))]; },
		[&](int ix, int iy, int iz) { return margin_data[(ix + size_x * (iy + size_y * (iz + margin_width)))]; },
		[&](int ix, int iy, int iz) { return margin_data[size_xy * margin_width+(ix + size_x * (iy + size_y * (iz - size_z)))]; }	);


	Laplasian2nd_halo_f(grid, out, buf_with_halo, coef, dx, dy, dz);

}


/*
* 1次元領域分割の場合のLaplasian
* 分割方向はz方向
*/
template <class T>
void Laplasian2nd_MPI_3D_overspace(const GridRangeMPI& grid, T* out, const T* src, T* buf_with_halo, double coef, double dx, double dy, double dz) {

	//のり代送受信//
	constexpr int TAG = 0x1;
	const int size_x = grid.SizeX();
	const int size_y = grid.SizeY();
	const int size_xy = size_x * size_y;
	const int size_z = grid.SizeZ();
	const int size_xz = size_x * size_z;
	const int size_yz = size_y * size_z;
	const int margin_width = 1;
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
	printf("[%d]test2-1\n", proc_id); fflush(stdout);
#endif
	TransferMarginX(grid, margin_x_left, margin_x_right, src, margin_width);
#ifdef TEST_PRINT
	printf("[%d]test2-2\n", proc_id); fflush(stdout);
#endif
	TransferMarginY(grid, margin_y_left, margin_y_right, src, margin_width);

#ifdef TEST_PRINT
	MPI_Barrier(grid.mpi_comm);
	printf("[%d]test2-3\n", proc_id); fflush(stdout);
#endif
	TransferMarginZ(grid, margin_z_left, margin_z_right, src, margin_width);

#ifdef TEST_PRINT
	MPI_Barrier(grid.mpi_comm);
	printf("[%d]test2-4\n", proc_id); fflush(stdout);
#endif

	PasteHalo_f1(grid, buf_with_halo, src,

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
	Laplasian2nd_halo_f(grid, out, buf_with_halo, coef, dx, dy, dz);

#ifdef TEST_PRINT
	MPI_Barrier(grid.mpi_comm);
	printf("[%d]test2-5\n", proc_id); fflush(stdout);
#endif

}


template <class T>
void Laplasian2nd_ddm(const GridRangeMPI& grid, T* out, const T* src, double coef, double dx, double dy, double dz) {

	const int margin_width = 1; 
	const int over_size = (grid.SizeX() + 2 * margin_width) * (grid.SizeY() + 2 * margin_width) * (grid.SizeZ() + 2 * margin_width);
	std::vector<T> buf_with_halo(over_size);//margin data//

#ifdef TEST_PRINT
	printf("L2 split_dimension = %d\n", grid.split_dimension);
#endif
	switch (grid.split_dimension) {
	case 1:
		Laplasian2nd_MPI_1D_overspace(grid, out, src, buf_with_halo.data(), coef, dx, dy, dz);
		break;
	case 2:
		Laplasian2nd_MPI_3D_overspace(grid, out, src, buf_with_halo.data(), coef, dx, dy, dz);
		break;
	case 3:
		Laplasian2nd_MPI_3D_overspace(grid, out, src, buf_with_halo.data(), coef, dx, dy, dz);
		break;
	default:		
		Laplasian2nd_halo(grid, out, src, buf_with_halo.data(), coef, dx, dy, dz);
		break;
	}
}

#endif

