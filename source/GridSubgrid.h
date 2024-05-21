#pragma once
#include "GridRange.h"
#ifdef USE_MPI
#include "mpi_helper.h"
#endif
#include <algorithm>

/*
* global_gridの範囲のデータから、gridで指定される範囲のデータを抜き出す
* 
* 使用条件：
* grid が完全にglobal_gridに収まっていること. すなわち, 
* (global_grid.begin_x <= grid.begin_x) && (grid.end_x <= global_grid.end_x)
*/
template <class T>
void CutSubgrid(const GridRange& grid, T* dest, const GridRange& global_grid, const T* src) {
	const int ix_end = grid.SizeX();
	const int iy_end = grid.SizeY();
	const int iz_end = grid.SizeZ();
	const int ix_offset = grid.begin_x - global_grid.begin_x;
	const int iy_offset = grid.begin_y - global_grid.begin_y;
	const int iz_offset = grid.begin_z - global_grid.begin_z;
	const int gwx = global_grid.SizeX();
	const int gwy = global_grid.SizeY();
	//const int gwz = global_grid.SizeZ();

	//printf("%d, %d, %d,   %d, %d, %d,  %d, %d", ix_end, iy_end, iz_end, ix_offset, iy_offset, iz_offset, gwx, gwy); fflush(stdout);

	for (int iz = 0; iz < iz_end; ++iz) {
		for (int iy = 0; iy < iy_end; ++iy) {
			for (int ix = 0; ix < ix_end; ++ix) {
				dest[ix + ix_end * (iy + iy_end * iz)] = src[(ix + ix_offset) + gwx * (iy + iy_offset + gwy * (iz + iz_offset))];
			}
		}
	}

}


/*
* global_gridの範囲のデータから、gridで指定される範囲のデータを抜き出す
*
* 使用条件：
* grid が完全にglobal_gridに収まっていること. すなわち,
* (global_grid.begin_x <= grid.begin_x) && (grid.end_x <= global_grid.end_x)
*/
template <class T, int STRIDE_D, int STRIDE_S>
void CutSubgrid_s(const GridRange& grid, T* dest, const GridRange& global_grid, const T* src) {
    const int ix_end = grid.SizeX();
    const int iy_end = grid.SizeY();
    const int iz_end = grid.SizeZ();
    const int ix_offset = grid.begin_x - global_grid.begin_x;
    const int iy_offset = grid.begin_y - global_grid.begin_y;
    const int iz_offset = grid.begin_z - global_grid.begin_z;
    const int gwx = global_grid.SizeX();
    const int gwy = global_grid.SizeY();
    //const int gwz = global_grid.SizeZ();

    //printf("%d, %d, %d,   %d, %d, %d,  %d, %d", ix_end, iy_end, iz_end, ix_offset, iy_offset, iz_offset, gwx, gwy); fflush(stdout);

    for (int iz = 0; iz < iz_end; ++iz) {
        for (int iy = 0; iy < iy_end; ++iy) {
            for (int ix = 0; ix < ix_end; ++ix) {
                dest[(ix + ix_end * (iy + iy_end * iz))* STRIDE_D] = src[((ix + ix_offset) + gwx * (iy + iy_offset + gwy * (iz + iz_offset)))* STRIDE_S];
            }
        }
    }

}


/*
* global_gridの範囲のデータに、gridで指定される範囲のデータを書き込む
*
* 使用条件：
* grid が完全にglobal_gridに収まっていること. すなわち,
* (global_grid.begin_x <= grid.begin_x) && (grid.end_x <= global_grid.end_x)
*/
template <class T>
void PasteSubgrid(const GridRange& global_grid, T* dest, const GridRange& grid, const T* sub_src) {
	const int ix_end = grid.SizeX();
	const int iy_end = grid.SizeY();
	const int iz_end = grid.SizeZ();
	const int ix_offset = grid.begin_x - global_grid.begin_x;
	const int iy_offset = grid.begin_y - global_grid.begin_y;
	const int iz_offset = grid.begin_z - global_grid.begin_z;
	const int gwx = global_grid.SizeX();
	const int gwy = global_grid.SizeY();
	//const int gwz = global_grid.SizeZ();
#pragma ivdep
	for (int iz = 0; iz < iz_end; ++iz) {
		for (int iy = 0; iy < iy_end; ++iy) {
			for (int ix = 0; ix < ix_end; ++ix) {
				dest[(ix + ix_offset) + gwx * (iy + iy_offset + gwy * (iz + iz_offset))] = sub_src[ix + ix_end * (iy + iy_end * iz)];
			}
		}
	}
}

template <class T>
void AddSubgrid(const GridRange& global_grid, T* dest, const GridRange& grid, const T* sub_src) {
	const int ix_end = grid.SizeX();
	const int iy_end = grid.SizeY();
	const int iz_end = grid.SizeZ();
	const int ix_offset = grid.begin_x - global_grid.begin_x;
	const int iy_offset = grid.begin_y - global_grid.begin_y;
	const int iz_offset = grid.begin_z - global_grid.begin_z;
	const int gwx = global_grid.SizeX();
	const int gwy = global_grid.SizeY();
	//const int gwz = global_grid.SizeZ();
#pragma ivdep
	for (int iz = 0; iz < iz_end; ++iz) {
		for (int iy = 0; iy < iy_end; ++iy) {
			for (int ix = 0; ix < ix_end; ++ix) {
				dest[(ix + ix_offset) + gwx * (iy + iy_offset + gwy * (iz + iz_offset))] += sub_src[ix + ix_end * (iy + iy_end * iz)];
			}
		}
	}
}


/*
* global_gridの範囲のデータから、gridで指定される範囲のデータを抜き出す
*
* 使用条件：
* grid が完全にglobal_gridに収まっていること. すなわち,
* (global_grid.begin_x <= grid.begin_x) && (grid.end_x <= global_grid.end_x)
*/
template <class T>
void CutSubgrid_bundle(int num_bundle, const GridRange& grid, T* dest, const GridRange& global_grid, const T* src) {
	const int ix_end = grid.SizeX();
	const int iy_end = grid.SizeY();
	const int iz_end = grid.SizeZ();
	const int ix_offset = grid.begin_x - global_grid.begin_x;
	const int iy_offset = grid.begin_y - global_grid.begin_y;
	const int iz_offset = grid.begin_z - global_grid.begin_z;
	const int gwx = global_grid.SizeX();
	const int gwy = global_grid.SizeY();
	//const int gwz = global_grid.SizeZ();

	//printf("%d, %d, %d,   %d, %d, %d,  %d, %d", ix_end, iy_end, iz_end, ix_offset, iy_offset, iz_offset, gwx, gwy); fflush(stdout);
#pragma ivdep
	for (int iz = 0; iz < iz_end; ++iz) {
		for (int iy = 0; iy < iy_end; ++iy) {
			for (int ix = 0; ix < ix_end; ++ix) {
				for (int ib = 0; ib < num_bundle; ++ib) {
					dest[ib + num_bundle * (int64_t)(ix + ix_end * (iy + iy_end * iz))] = src[ib + num_bundle * (int64_t)((ix + ix_offset) + gwx * (iy + iy_offset + gwy * (iz + iz_offset)))];
				}
			}
		}
	}

}



template <class T>
void AddSubgrid_bundle(int num_bundle, const GridRange& global_grid, T* dest, const GridRange& grid, const T* sub_src) {
	const int ix_end = grid.SizeX();
	const int iy_end = grid.SizeY();
	const int iz_end = grid.SizeZ();
	const int ix_offset = grid.begin_x - global_grid.begin_x;
	const int iy_offset = grid.begin_y - global_grid.begin_y;
	const int iz_offset = grid.begin_z - global_grid.begin_z;
	const int gwx = global_grid.SizeX();
	const int gwy = global_grid.SizeY();
	//const int gwz = global_grid.SizeZ();
#pragma ivdep
	for (int iz = 0; iz < iz_end; ++iz) {
		for (int iy = 0; iy < iy_end; ++iy) {
			for (int ix = 0; ix < ix_end; ++ix) {
				for (int ib = 0; ib < num_bundle; ++ib) {
					dest[ib + num_bundle * (int64_t)((ix + ix_offset) + gwx * (iy + iy_offset + gwy * (iz + iz_offset)))] += sub_src[ib + num_bundle * (int64_t)(ix + ix_end * (iy + iy_end * iz))];
				}
			}
		}
	}
}



/*Overlapの折り畳み回数算出用

#include <cstdio>
#include <cstdlib>
#include <vector>
#include <random>

int main()
{
	const int begin_g = 0;
	const int end_g = 10;

	const int width = end_g - begin_g;

	const int begin_sub = -13;
	const int end_sub = -11;


	auto n_gt_div = [](int x, int y) {
		return (x >= 0) ? (x / y) + 1 : ((x + 1) / y);
		};

	const int n = n_gt_div(begin_g - end_sub, width);

	printf("n = %d, range = [%d, %d)\n", n, begin_sub + n * width, end_sub + n * width);
	return n;
}
*/


/*
grid1 と grid2のoverlapする範囲内でループを回す。
ここで、overlap範囲内は、
max(grid1.begin_x, grid2.begin_x) to min(grid1.end_x, grid2.end_x) 
である。
*/
template<class FUNC>
void ForOverlapRange(const GridRange& grid1, const GridRange& grid2, FUNC&& func){
	const int ix_begin = std::max(grid1.begin_x, grid2.begin_x);
	const int iy_begin = std::max(grid1.begin_y, grid2.begin_y);
	const int iz_begin = std::max(grid1.begin_z, grid2.begin_z);
	const int ix_end = std::min(grid1.end_x, grid2.end_x);
	const int iy_end = std::min(grid1.end_y, grid2.end_y);
	const int iz_end = std::min(grid1.end_z, grid2.end_z);

	const int wx1 = grid1.SizeX();
	const int wy1 = grid1.SizeY();
	const int wz1 = grid1.SizeZ();
	const int wx2 = grid2.SizeX();
	const int wy2 = grid2.SizeY();
	const int wz2 = grid2.SizeZ();

#ifdef __DEBUG
	const int proc_id = GetProcessID(MPI_COMM_WORLD);
	printf("[%d] ibegin-end = %d, %d, %d, %d, %d, %d\n", proc_id, ix_begin, iy_begin, iz_begin, ix_end, iy_end, iz_end);
	//printf("[%d] wx1,wx2 = %d, %d, %d, %d, %d, %d\n", proc_id, wx1, wy1, wz1, wx2, wy2, wz2);
	fflush(stdout);
#endif

	int64_t sum = 0;
	for(int iz = iz_begin; iz < iz_end; ++iz){
		for(int iy = iy_begin; iy < iy_end; ++iy){
			for(int ix = ix_begin; ix < ix_end; ++ix){
				const int64_t i1 = (ix - grid1.begin_x) + wx1 * ((iy - grid1.begin_y) + wy1 * (iz - grid1.begin_z));
				const int64_t i2 = (ix - grid2.begin_x) + wx2 * ((iy - grid2.begin_y) + wy2 * (iz - grid2.begin_z));
				func(i1, i2);
				++sum;
	
			}
		}
	}
	
}


/*
grid1 と grid2のoverlapする範囲内でループを回す。
ただし、周期境界条件としてgrid2はperi_width_{x,y,z}で自由にシフトできる
このとき、シフト回数は、負方向にns回、正方向にne回のシフトをしても、grid1とoverlapできる。
言い換えれば、負方向にns+1回以上、正方向にne+1回以上シフトした場合はgrid1とoverlapしなくなる。
各overlap範囲毎に、ForOverlapRange関数を読んで、グリッドの内部のループを行う.

例としては、grid1はfixされたglobal_gridや、ddm分割されたgrid,
grid2はpseudo-potentialのnon-local項などの原子に紐づいた部分波をイメージしている.
*/
template<class FUNC>
void ForOverlapPeriodic(const GridRange& grid1, const GridRange& grid2,
	int peri_width_x, int peri_width_y, int peri_width_z, FUNC&& func){
	
	auto n_gt_div = [](int x, int y) {
		return (x >= 0) ? (x / y) + 1 : ((x + 1) / y);
		};
	const int ns_x = n_gt_div(grid1.begin_x - grid2.end_x, peri_width_x);
    const int ne_x = -n_gt_div(grid2.begin_x - grid1.end_x, peri_width_x);
	const int ns_y = n_gt_div(grid1.begin_y - grid2.end_y, peri_width_y);
    const int ne_y = -n_gt_div(grid2.begin_y - grid1.end_y, peri_width_y);
	const int ns_z = n_gt_div(grid1.begin_z - grid2.end_z, peri_width_z);
    const int ne_z = -n_gt_div(grid2.begin_z - grid1.end_z, peri_width_z);
	
	//printf("ns,ne = %d, %d, %d, %d, %d, %d\n", ns_x, ne_x, ns_y, ne_y, ns_z, ne_z);

	GridRange grid2_shift;
	for(int nz = ns_z; nz <= ne_z;++nz){
		grid2_shift.begin_z = grid2.begin_z + nz * peri_width_z;
		grid2_shift.end_z = grid2.end_z + nz * peri_width_z;
		for(int ny = ns_y; ny <= ne_y;++ny){
			grid2_shift.begin_y = grid2.begin_y + ny * peri_width_y;
			grid2_shift.end_y = grid2.end_y + ny * peri_width_y;
		
			for(int nx = ns_x; nx <= ne_x;++nx){
				grid2_shift.begin_x = grid2.begin_x + nx * peri_width_x;
				grid2_shift.end_x = grid2.end_x + nx * peri_width_x;
		
				ForOverlapRange(grid1, grid2_shift, func);
			}
		}
	}
}

inline
GridRange GetOverlapRange(const GridRange& grid1, const GridRange& grid2) {
	const int ix_begin = std::max(grid1.begin_x, grid2.begin_x);
	const int iy_begin = std::max(grid1.begin_y, grid2.begin_y);
	const int iz_begin = std::max(grid1.begin_z, grid2.begin_z);
	const int ix_end = std::min(grid1.end_x, grid2.end_x);
	const int iy_end = std::min(grid1.end_y, grid2.end_y);
	const int iz_end = std::min(grid1.end_z, grid2.end_z);

	return GridRange{ ix_begin, iy_begin, iz_begin, ix_end, iy_end, iz_end };
}


template<class T>
void CutSubgridPeriodic(const GridRange& grid1, const T* src, const GridRange& grid2, T* sub_dest,
	int peri_width_x, int peri_width_y, int peri_width_z) {

	ForOverlapPeriodic(grid1, grid2,
		peri_width_x, peri_width_y, peri_width_z,
		[&sub_dest, &src](const int64_t i1, const int64_t i2) {
			sub_dest[i2] = src[i1];
		});
}

template<class T>
void PasteSubgridPeriodic(const GridRange& grid1, T* dest, const GridRange& grid2, const T* sub_src,
	int peri_width_x, int peri_width_y, int peri_width_z) {

	ForOverlapPeriodic(grid1, grid2,
		peri_width_x, peri_width_y, peri_width_z,
		[&dest, &sub_src](const int64_t i1, const int64_t i2) {
			dest[i1] = sub_src[i2];
		});
}

template<class T>
void AddSubgridPeriodic(const GridRange& grid1, T* dest, const GridRange& grid2, const T* sub_src,
	int peri_width_x, int peri_width_y, int peri_width_z) {

	ForOverlapPeriodic(grid1, grid2,
		peri_width_x, peri_width_y, peri_width_z,
		[&dest, &sub_src](const int64_t i1, const int64_t i2) {
			dest[i1] += sub_src[i2];
		});
}

struct ShiftedGrid {
	int x;
	int y;
	int z;
};


struct SubgridBlock {

	std::vector<GridRange> range_blocks;
	std::vector<ShiftedGrid> shifted_grids;//ブロックが周期境界を跨いで切り取られたときに、空間中の座標として何周期分シフトしているかをしめしたもの//
	int grid_sizes;//blockの総グリッド数//
	double cutoff;
	/********************************
	* note: block化
	* 擬ポテンシャルの存在する部分領域と、
	* 系全体の領域分割によって本mpi-processが保持する波動関数の領域の
	* オーバーラップ部分だけに相当する領域のデータをblock化されたデータと呼ぶことにする.
	* 周期境界の為に、1つの擬ポテンシャルが複数のブロックで構成されることがあり得るが,
	* 配列内では連続的につないでしまって一つの配列として保持する.
	* 内積を取る演算時も一続きにアクセスして良い.
	*****************************/
	/*
	~PP_Vlocal_block() {
		delete[] rho_vlocal;
	}
	*/
};


/*
* ForOverlapPeriodicと同様の領域(複数の領域の束)を返す
*/
inline
int GetOverlapPeriodic(const GridRange& grid1, const GridRange& grid2,
	int peri_width_x, int peri_width_y, int peri_width_z, 
	std::vector<GridRange>& ranges_in_grid1, std::vector<ShiftedGrid>& shifted_grid2) {

	auto n_gt_div = [](int x, int y) {
		return (x >= 0) ? (x / y) + 1 : ((x + 1) / y);
		};
	const int ns_x = n_gt_div(grid1.begin_x - grid2.end_x, peri_width_x);
	const int ne_x = -n_gt_div(grid2.begin_x - grid1.end_x, peri_width_x);
	const int ns_y = n_gt_div(grid1.begin_y - grid2.end_y, peri_width_y);
	const int ne_y = -n_gt_div(grid2.begin_y - grid1.end_y, peri_width_y);
	const int ns_z = n_gt_div(grid1.begin_z - grid2.end_z, peri_width_z);
	const int ne_z = -n_gt_div(grid2.begin_z - grid1.end_z, peri_width_z);

	//printf("ns,ne = %d, %d, %d, %d, %d, %d\n", ns_x, ne_x, ns_y, ne_y, ns_z, ne_z);

	int total_grids = 0;
	for (int nz = ns_z; nz <= ne_z; ++nz) {
		GridRange grid2_shift;
		grid2_shift.begin_z = grid2.begin_z + nz * peri_width_z;
		grid2_shift.end_z = grid2.end_z + nz * peri_width_z;
		for (int ny = ns_y; ny <= ne_y; ++ny) {
			grid2_shift.begin_y = grid2.begin_y + ny * peri_width_y;
			grid2_shift.end_y = grid2.end_y + ny * peri_width_y;

			for (int nx = ns_x; nx <= ne_x; ++nx) {
				grid2_shift.begin_x = grid2.begin_x + nx * peri_width_x;
				grid2_shift.end_x = grid2.end_x + nx * peri_width_x;

				GridRange range_in_1 = GetOverlapRange(grid1, grid2_shift);
				ranges_in_grid1.emplace_back(range_in_1);
				total_grids += range_in_1.Size3D();
#ifdef DEBUG_PRINT
				printf("range_1: [%d, %d) * [%d, %d) * [%d, %d)\n",
					range_in_1.begin_x, range_in_1.end_x, range_in_1.begin_y, range_in_1.end_y, range_in_1.begin_z, range_in_1.end_z);
				fflush(stdout);
#endif

				shifted_grid2.push_back(ShiftedGrid{ nx * peri_width_x, ny * peri_width_y, nz * peri_width_z });

#ifdef DEBUG_PRINT
				printf("range_2: [%d, %d) * [%d, %d) * [%d, %d)\n",
					range_in_2.begin_x, range_in_2.end_x, range_in_2.begin_y, range_in_2.end_y, range_in_2.begin_z, range_in_2.end_z);
				fflush(stdout);
#endif
			}
		}
	}
	return total_grids;
}



template<class T>
void CutSubgridByRanges(const std::vector<GridRange>& ranges_in_1, T* dest, const GridRange& grid1, const T* src) {
	size_t offset = 0;
	for (const auto& range : ranges_in_1) {
		CutSubgrid(range, dest + offset, grid1, src);
		offset += range.Size3D();
	}	
}

template<class T>
void AddSubgridByRanges(const GridRange& grid1, T* dest, const std::vector<GridRange>& ranges_in_1, const T* src) {
	size_t offset = 0;
	for (const auto& range : ranges_in_1) {
		AddSubgrid(grid1, dest, range, src + offset);
		offset += range.Size3D();
	}
}


template<class T>
void CutSubgridByRanges_SOAtoCOMPLEX(const std::vector<GridRange>& ranges_in_1, T* dest, const GridRange& grid1, const T* src_re, const T* src_im) {
    size_t offset = 0;
    for (const auto& range : ranges_in_1) {
        CutSubgrid_s<double, 2, 1>(range, dest + offset * 2, grid1, src_re);
        CutSubgrid_s<double, 2, 1>(range, dest + offset*2+1, grid1, src_im);
        offset += range.Size3D();
    }
}

template<class T>
void CutSubgridByRanges_bundle(int num_bundle, const std::vector<GridRange>& ranges_in_1, T* dest, const GridRange& grid1, const T* src) {
	size_t offset = 0;
	for (const auto& range : ranges_in_1) {
		CutSubgrid_bundle(num_bundle, range, dest + offset, grid1, src);
		offset += range.Size3D()* num_bundle;
	}
}

template<class T>
void AddSubgridByRanges_bundle(int num_bundle, const GridRange& grid1, T* dest, const std::vector<GridRange>& ranges_in_1, const T* src) {
	size_t offset = 0;
	for (const auto& range : ranges_in_1) {
		AddSubgrid_bundle(num_bundle, grid1, dest, range, src + offset);
		offset += range.Size3D()* num_bundle;
	}
}


inline
GridRange GridForSphere(double Rx, double Ry, double Rz, double cutoff_r, double dx, double dy, double dz) {
	GridRange subgrid;
	subgrid.begin_x = (int)ceil((Rx - cutoff_r) / dx - 0.5);
	subgrid.end_x = (int)floor((Rx + cutoff_r) / dx + 0.5) + 1;
	subgrid.begin_y = (int)ceil((Ry - cutoff_r) / dy - 0.5);
	subgrid.end_y = (int)floor((Ry + cutoff_r) / dy + 0.5) + 1;
	subgrid.begin_z = (int)ceil((Rz - cutoff_r) / dz - 0.5);
	subgrid.end_z = (int)floor((Rz + cutoff_r) / dz + 0.5) + 1;
	return subgrid;
}
