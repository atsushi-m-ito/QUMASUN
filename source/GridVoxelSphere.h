#pragma once
#include <vector>
#include "GridRange.h"

/********
* XYX, cartecian gridから球状のvoxel領域を切り出す際の領域データを作成する.
* 切り出すのはx軸に沿ったラインで,始点(xyz)と幅(w=x_end - x)のデータである.
* 
*********/
inline
int MakeVoxelSphereRange(double Rx, double Ry, double Rz, double cutoff_r,
	double dx, double dy, double dz, std::vector<int>& voxel_lines){

	const int iz_begin = (int)ceil((Rz - cutoff_r) / dz - 0.5);
	const int iz_end = (int)floor((Rz + cutoff_r) / dz + 0.5) + 1;

	const double rr_max = cutoff_r;

	auto pow2 = [](double x) { return x * x; };

	int grid_size = 0;

	for (int iz = iz_begin; iz < iz_end; ++iz) {
		const double rz2 = pow2(dz * (double)(iz) - Rz);
		if (rz2 > rr_max) continue;

		const double curoff_r2 = sqrt(rr_max - rz2);
		const int iy_begin = (int)ceil((Ry - curoff_r2) / dy - 0.5);
		const int iy_end = (int)floor((Ry + curoff_r2) / dy + 0.5) + 1;

		for (int iy = iy_begin; iy < iy_end; ++iy) {
			const double ry2 = pow2(dy * (double)(iy) - Ry);
			if (rz2 + ry2 > rr_max) continue;

			const double curoff_r3 = sqrt(rr_max - rz2 - ry2);
			const int ix_begin = (int)ceil((Rx - curoff_r3) / dx - 0.5);
			const int ix_end = (int)floor((Rx + curoff_r3) / dx + 0.5) + 1;

			voxel_lines.push_back(iz);
			voxel_lines.push_back(iy);
			voxel_lines.push_back(ix_begin);
			voxel_lines.push_back(ix_end - ix_begin);

			grid_size += ix_end - ix_begin;
		}
	}

	return grid_size;
}


inline
int MakeVoxelSphereRangeOverrap(const GridRange& grid, double Rx, double Ry, double Rz, double cutoff_r,
	double dx, double dy, double dz, std::vector<int>& voxel_lines, int* pnum_line) {

	/*auto is_overlap = [](int begin1, int end1, int begin2, int end2) {
		return (begin2 < end1) && (begin1 < end2);
		};
		*/

	const int iz_begin = (int)ceil((Rz - cutoff_r) / dz - 0.5);
	const int iz_end = (int)floor((Rz + cutoff_r) / dz + 0.5) + 1;

	const int iz_begin2 = std::max(iz_begin, grid.begin_z);
	const int iz_end2 = std::min(iz_end, grid.end_z);
	//if(!is_overlap(grid.begin_z, grid.end_z, iz_begin, iz_end)) return 0;	

	const double rr_max = cutoff_r* cutoff_r;

	auto pow2 = [](double x) { return x * x; };

	int num_lines = 0;
	int num_grids = 0;

	for (int iz = iz_begin2; iz < iz_end2; ++iz) {
		const double rz2 = pow2(dz * (double)(iz)-Rz);
		if (rz2 > rr_max) continue;

		const double curoff_r2 = sqrt(rr_max - rz2);
		const int iy_begin = (int)ceil((Ry - curoff_r2) / dy - 0.5);
		const int iy_end = (int)floor((Ry + curoff_r2) / dy + 0.5) + 1;

		const int iy_begin2 = std::max(iy_begin, grid.begin_y);
		const int iy_end2 = std::min(iy_end, grid.end_y);

		for (int iy = iy_begin2; iy < iy_end2; ++iy) {
			const double ry2 = pow2(dy * (double)(iy)-Ry);
			if (rz2 + ry2 > rr_max) continue;

			const double curoff_r3 = sqrt(rr_max - rz2 - ry2);
			const int ix_begin = (int)ceil((Rx - curoff_r3) / dx - 0.5);
			const int ix_end = (int)floor((Rx + curoff_r3) / dx + 0.5) + 1;

			const int ix_begin2 = std::max(ix_begin, grid.begin_x);
			const int ix_end2 = std::min(ix_end, grid.end_x);
			const int width = ix_end2 - ix_begin2;
			if (width > 0) {
				voxel_lines.push_back(iz);
				voxel_lines.push_back(iy);
				voxel_lines.push_back(ix_begin2);
				voxel_lines.push_back(width);
				++num_lines;
				num_grids += width;
			}
		}
		
	}
	*pnum_line = num_lines;
	return num_grids;
}

struct ShiftedPos{
	double x;
	double y;
	double z;
	int num_lines;
};

/*
* ForOverlapPeriodicと同様の領域(複数の領域の束)を返す
*/
inline
int GetVsOverlapPeriodic(const GridRange& grid1, 
	double Rx, double Ry, double Rz, double cutoff_r,
	double dx, double dy, double dz, 
	int peri_width_x, int peri_width_y, int peri_width_z,
	std::vector<int>& voxel_lines1, std::vector<ShiftedPos>& shifted_center) {

	const int ix_begin = (int)ceil((Rx - cutoff_r) / dx - 0.5);
	const int ix_end = (int)floor((Rx + cutoff_r) / dx + 0.5) + 1;
	const int iy_begin = (int)ceil((Ry - cutoff_r) / dy - 0.5);
	const int iy_end = (int)floor((Ry + cutoff_r) / dy + 0.5) + 1;
	const int iz_begin = (int)ceil((Rz - cutoff_r) / dz - 0.5);
	const int iz_end = (int)floor((Rz + cutoff_r) / dz + 0.5) + 1;

	auto n_gt_div = [](int x, int y) {
		return (x >= 0) ? (x / y) + 1 : ((x + 1) / y);
		};
	const int ns_x = n_gt_div(grid1.begin_x - ix_end, peri_width_x);
	const int ne_x = -n_gt_div(ix_begin - grid1.end_x, peri_width_x);
	const int ns_y = n_gt_div(grid1.begin_y - iy_end, peri_width_y);
	const int ne_y = -n_gt_div(iy_begin - grid1.end_y, peri_width_y);
	const int ns_z = n_gt_div(grid1.begin_z - iz_end, peri_width_z);
	const int ne_z = -n_gt_div(iz_begin - grid1.end_z, peri_width_z);

	//printf("ns,ne = %d, %d, %d, %d, %d, %d\n", ns_x, ne_x, ns_y, ne_y, ns_z, ne_z);

	
	int total_size = 0;
	for (int nz = ns_z; nz <= ne_z; ++nz) {
		const double Rz_p = Rz + dz*(double)(nz * peri_width_z);
		
		for (int ny = ns_y; ny <= ne_y; ++ny) {
			const double Ry_p = Ry + dy * (double)(ny * peri_width_y);
		
			for (int nx = ns_x; nx <= ne_x; ++nx) {
				const double Rx_p = Rx + dx * (double)(nx * peri_width_x);

				int num_line;
				int grid_size = MakeVoxelSphereRangeOverrap(grid1, Rx_p, Ry_p, Rz_p, cutoff_r, dx, dy, dz, voxel_lines1, &num_line);
			
				if (num_line > 0) {
					shifted_center.push_back(ShiftedPos{ Rx_p, Ry_p, Rz_p, num_line });
				}
				total_size += grid_size;
			}
		}
	}
	return total_size;
}



template <class T>
inline
int CutVsgridByRanges(const std::vector<int>& voxel_lines1, T* vs_dest, const GridRange& grid1, const T* src) {
	const size_t num_lines4 = voxel_lines1.size();
	int total_i = 0;
	for (size_t il = 0; il < num_lines4; il += 4) {
		int iz_begin = voxel_lines1[il + 0];
		int iy_begin = voxel_lines1[il + 1];
		int ix_begin = voxel_lines1[il + 2];
		int ix_width = voxel_lines1[il + 3];
		size_t src_offset = (ix_begin - grid1.begin_x) + grid1.SizeX() * ((iy_begin - grid1.begin_y) + grid1.SizeY() * (iz_begin - grid1.begin_z));
		for (int ix = 0; ix < ix_width; ++ix) {			
			vs_dest[total_i + ix] = src[ix + src_offset];
		}
		total_i += ix_width;
	}
	return total_i;
}

template <class T>
inline
int AddVsgridByRanges(const std::vector<int>& voxel_lines1, const T* vs_src, const GridRange& grid1, T* dest) {
	const size_t num_lines4 = voxel_lines1.size();
	int total_i = 0;
	for (size_t il = 0; il < num_lines4; il += 4) {
		int iz_begin = voxel_lines1[il + 0];
		int iy_begin = voxel_lines1[il + 1];
		int ix_begin = voxel_lines1[il + 2];
		int ix_width = voxel_lines1[il + 3];
		size_t src_offset = (ix_begin - grid1.begin_x) + grid1.SizeX() * ((iy_begin - grid1.begin_y) + grid1.SizeY() * (iz_begin - grid1.begin_z));
		for (int ix = 0; ix < ix_width; ++ix) {
			dest[ix + src_offset]+=vs_src[total_i + ix];
		}
		total_i += ix_width;
	}
	return total_i;
}
