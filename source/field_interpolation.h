#pragma once


#include <cmath>
#include <cstdio>
#include <cstdint>
#include <random>
#include "GridRange.h"


inline
bool IsOverlap(const GridI3D& begin1, const GridI3D& end1, const GridI3D& begin2, const GridI3D& end2)
{
	return (begin2.x < end1.x) && (begin1.x < end2.x) &&
		(begin2.y < end1.y) && (begin1.y < end2.y) &&
		(begin2.z < end1.z) && (begin1.z < end2.z);
}

/*
* srcから[begin, begin + width)で指定された区間を切り取りdestにコピーする//
* 指定区間が範囲外の時は何もしない
*/
template <class T>
void CutoutField(T* dest, const GridRange& cut, 
	const T* src, const GridRange& src_range)
{
	const int ix_begin = std::max(cut.begin_x - src_range.begin_x, 0);
	const int ix_end = std::min(cut.end_x, src_range.end_x) - src_range.begin_x;
	const int iy_begin = std::max(cut.begin_y - src_range.begin_y, 0);
	const int iy_end = std::min(cut.end_y, src_range.end_y) - src_range.begin_y;
	const int iz_begin = std::max(cut.begin_z - src_range.begin_z, 0);
	const int iz_end = std::min(cut.end_z, src_range.end_z) - src_range.begin_z;
	const int i2k_x = src_range.begin_x - cut.begin_x;
	const int i2k_y = src_range.begin_y - cut.begin_y;
	const int i2k_z = src_range.begin_z - cut.begin_z;
	const int dest_wx = cut.end_x - cut.begin_x;
	const int dest_wy = cut.end_y - cut.begin_y;
	const int size_x = src_range.end_x - src_range.begin_x;
	const int size_y = src_range.end_y - src_range.begin_y;

	for (int iz = iz_begin; iz < iz_end; ++iz) {
		const int kz = (iz + i2k_z) * dest_wx * dest_wy;
		for (int iy = iy_begin; iy < iy_end; ++iy) {
			const int ky = (iy  + i2k_y) * dest_wx;
			for (int ix = ix_begin; ix < ix_end; ++ix) {
				const int k = (ix + i2k_x) + ky + kz;
				const int i = ix + size_x * (iy + size_y * iz);
				dest[k] = src[i];
			}
		}
	}
}



/*
* 整数を対象とした計算
* ある点xが周期cycle_lengthで繰り返される時に(繰り返し点をミラー)、
* 範囲( >= range_begin)の中に納まる最小のミラーx_nを返す
* ただし、ある整数nを用いて
* x_n = x - cycle_length * n 
* であり、
* x_n >= range_begin
* かつ
* x_n - cycle_length < range_begin
* となる
*/
inline
int BeginCycle(int x, int cycle_length, int range_begin) {
	auto imodpositive = [](int x, int N) {
		return (x >= 0) ? x % N : (N - ((-x - 1) % N + 1));
		};
	const int x_n = imodpositive(x - range_begin, cycle_length) + range_begin;
	return x_n;
}


/*
* 整数を対象とした計算
* ある範囲[x0,x1)が周期cycle_lengthで繰り返される時に(繰り返し点をミラー)、
* 範囲[r0,r1)とオーバーラップする最小のx0_beginと, 
* x0_begin以上でオーバーラップしなくなる最小のx0_endを返す
* つまり、ある整数nを用いて
* x_n = x1 - 1 + n * cycle_length
* としたとき、
* x_n >= r0
* となる最小のnを用いて
* x0_begin = x0 + n * cycle_length
* となる. 
* ある整数mを用いて
* x_m = x0 + m * cycle_length
* としたとき、
* x_m >= r1
* となる最小のmを用いて
* x0_end = x0 + m * cycle_length
* とする.
* 
* 使い方
* int x0 = 3;       //subspace in field
* int x1 = 13;      //subspace in field
* int cycle_length = 10;  //global box_size
* int r0 = 5;       // local box_region
* int r1 = 10;      // local box_region
*
* Example of usage:
* int x0_end;
* int x0_begin = OverlapRangeInCycle(x0, x1, cycle_length, r0, r1, &x0_end);
* int width = x1 - x0;
* for(int x0i = x0_begin; x0i = x0_end; ++x0i){
*   CutoutField(dest, x0i, x0i + width, src, r0, r1);
* }
*/
inline
int OverlapRangeInCycle(int x0, int x1, int cycle_length, int r0, int r1, int* x0_end) {
	auto imodpositive = [](int x, int N) {
		return (x >= 0) ? x % N : (N - ((-x - 1) % N + 1));
		};
	const int x0_begin = BeginCycle(x1 - 1, cycle_length, r0) - (x1-1-x0);
	*x0_end = BeginCycle(x0, cycle_length, r1);
	return x0_begin;
}

/*
* srcから[begin, begin + width)で指定された区間を切り取りdestにコピーする//
* 指定区間が範囲外の時はperiodic boundaryとして範囲内に折り畳みを行ってsrcを読み込む
*
*/
template <class T>
void CutoutField_periodic(T* dest, const GridRange& cut,
	const T* src, const GridRange& src_region, const GridI3D& global_size)
{
	auto imodpositive = [](int x, int N) {
		return (x >= 0) ? x % N : (N - ((-x - 1) % N + 1));
		};

	int ib_x_end, ib_y_end, ib_z_end;
	int ib_x = OverlapRangeInCycle(cut.begin_x, cut.end_x, global_size.x, src_region.begin_x, src_region.end_x, &ib_x_end);
	int ib_y = OverlapRangeInCycle(cut.begin_y, cut.end_y, global_size.y, src_region.begin_y, src_region.end_y, &ib_y_end);
	int ib_z = OverlapRangeInCycle(cut.begin_z, cut.end_z, global_size.z, src_region.begin_z, src_region.end_z, &ib_z_end);


	for (int iz0 = ib_z; iz0 < ib_z_end; iz0 += global_size.z) {
		for (int iy0 = ib_y; iy0 < ib_y_end; iy0 += global_size.y) {
			for (int ix0 = ib_x; ix0 < ib_x_end; ix0 += global_size.x) {
				GridRange cut_local = { ix0, iy0, iz0, ix0 + 4, iy0 + 4, iz0 + 4 };
				CutoutField(dest, cut_local, src, src_region);								
			}
		}
	}
}



/*
グリッドデータから近似的にエルミート補間する

The following transform should be calculated before input into argument
	scaled_Rx = Rx / dx

	field: field data with the size of (size_x * size_y * size_z)

	org_x: offset of field index;
*/
template <class T>
T Interpolation(double scaled_Rx, double scaled_Ry, double scaled_Rz,
	const T* src,
	const GridRange& src_range, const GridI3D& global_size)
{
	const int ix = (int)floor(scaled_Rx);
	const int iy = (int)floor(scaled_Ry);
	const int iz = (int)floor(scaled_Rz);


	T data[4 * 4 * 4] = { 0.0 };//補完に必要な4x4x4領域を切り取ったデータ//
	//4x4x4領域の切り取り(領域外の場合は0が埋まる)
	
	CutoutField_periodic(data, { ix - 1, iy - 1, iz - 1, ix + 3, iy + 3, iz + 3 },
		src, src_range, global_size);



	const double x = scaled_Rx - (double)ix;
	const double y = scaled_Ry - (double)iy;
	const double z = scaled_Rz - (double)iz;

	const double x2 = x * x;
	const double x3 = x * x * x;
	const double y2 = y * y;
	const double y3 = y * y * y;
	const double z2 = z * z;
	const double z3 = z * z * z;

	//Cubic Hermite Interpolation//
	const double g1x = 3.0 * x2 - 2.0 * x3;
	const double g0x = 1.0 - g1x;
	const double g3x = -x2 + x3;
	const double g2x = x - x2 + g3x;

	const double g1y = 3.0 * y2 - 2.0 * y3;
	const double g0y = 1.0 - g1y;
	const double g3y = -y2 + y3;
	const double g2y = y - y2 + g3y;

	const double g1z = 3.0 * z2 - 2.0 * z3;
	const double g0z = 1.0 - g1z;
	const double g3z = -z2 + z3;
	const double g2z = z - z2 + g3z;


	T Fyz[4 * 4];	//x方向の補間後のデータ//
	T Fz[4 * 4];	//x,y方向の補間後のデータ//

	for (int byz = 0; byz < 16; ++byz) {
		Fyz[byz] = -g2x / 2.0 * data[0 + byz * 4] + (g0x - g3x / 2.0) * data[1 + byz * 4] + (g1x + g2x / 2.0) * data[2 + byz * 4] + g3x / 2.0 * data[3 + byz * 4];
	}

	for (int bz = 0; bz < 4; ++bz) {
		Fz[bz] = -g2y / 2.0 * Fyz[0 + bz * 4] + (g0y - g3y / 2.0) * Fyz[1 + bz * 4] + (g1y + g2y / 2.0) * Fyz[2 + bz * 4] + g3y / 2.0 * Fyz[3 + bz * 4];
	}

	const double f_inter = -g2z / 2.0 * Fz[0] + (g0z - g3z / 2.0) * Fz[1] + (g1z + g2z / 2.0) * Fz[2] + g3z / 2.0 * Fz[3];

	return f_inter;


}

