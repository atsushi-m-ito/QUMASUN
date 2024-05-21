#pragma once
#include "GridFor.h"
#include "GridSubgrid.h"
#include "GridSubgrid_arithmetic.h"
#include "GridVoxelSphere.h"
#include "RealSphericalHarmonics.h"


template<class YLM>
void SetBlockYlmAroundAny(YLM Y, double value_org, double* buf, const std::vector<GridRange>& range_blocks,
	const std::vector<ShiftedGrid>& shifted_grid,
	double Rx, double Ry, double Rz,
	double dx, double dy, double dz)
{

	size_t num_blocks = range_blocks.size();
	constexpr double rr_min = 1.0e-14;
	int offset_i = 0;
	for (size_t ib = 0; ib < num_blocks; ++ib) {
		const int ix_begin = range_blocks[ib].begin_x;
		const int iy_begin = range_blocks[ib].begin_y;
		const int iz_begin = range_blocks[ib].begin_z;
		const int ix_end = range_blocks[ib].end_x;
		const int iy_end = range_blocks[ib].end_y;
		const int iz_end = range_blocks[ib].end_z;

		const int size_x = ix_end - ix_begin;
		const int size_y = iy_end - iy_begin;

		for (int iz = iz_begin; iz < iz_end; ++iz) {
			const double z = dz * (double)(iz - shifted_grid[ib].z) - Rz;
			const double zz = z * z;
			for (int iy = iy_begin; iy < iy_end; ++iy) {
				const double y = dy * (double)(iy - shifted_grid[ib].y) - Ry;
				const double yy_zz = y * y + zz;
				for (int ix = ix_begin; ix < ix_end; ++ix) {
					const double x = dx * (double)(ix - shifted_grid[ib].x) - Rx;
					const double rr = x * x + yy_zz;
					const int i = offset_i + (ix - ix_begin) + size_x * ((iy - iy_begin) + size_y * (iz - iz_begin));

					if (rr_min > rr) {
						buf[i] = value_org;
					} else {
						const double r = sqrt(rr);
						buf[i] = Y(x / r, y / r, z / r);
					}
				}
			}
		}
		offset_i += range_blocks[ib].Size3D();
	}

}

inline
void SetBlockYlmAround(int L, int M, double* buf, const std::vector<GridRange>& range_blocks,
	const std::vector<ShiftedGrid>& shifted_grid,
	double Rx, double Ry, double Rz,
	double dx, double dy, double dz)
{
	const double Y00 = 1.0 / sqrt(4.0 * M_PI);
	switch (L * L + L + M) {
	case 0:
	{
		Ylm<0, 0> Y;
		SetBlockYlmAroundAny(Y, Y00, buf, range_blocks, shifted_grid, Rx, Ry, Rz, dx, dy, dz);
		break;
	}
	case 1:
	{
		Ylm<1, -1> Y;
		SetBlockYlmAroundAny(Y, 0.0, buf, range_blocks, shifted_grid, Rx, Ry, Rz, dx, dy, dz);
		break;
	}
	case 2:
	{
		Ylm<1, 0> Y;
		SetBlockYlmAroundAny(Y, 0.0, buf, range_blocks, shifted_grid, Rx, Ry, Rz, dx, dy, dz);
		break;
	}
	case 3:
	{
		Ylm<1, 1> Y;
		SetBlockYlmAroundAny(Y, 0.0, buf, range_blocks, shifted_grid, Rx, Ry, Rz, dx, dy, dz);
		break;
	}
	case 4:
	{
		Ylm<2, -2> Y;
		SetBlockYlmAroundAny(Y, 0.0, buf, range_blocks, shifted_grid, Rx, Ry, Rz, dx, dy, dz);
		break;
	}
	case 5:
	{
		Ylm<2, -1> Y;
		SetBlockYlmAroundAny(Y, 0.0, buf, range_blocks, shifted_grid, Rx, Ry, Rz, dx, dy, dz);
		break;
	}
	case 6:
	{
		Ylm<2, 0> Y;
		SetBlockYlmAroundAny(Y, 0.0, buf, range_blocks, shifted_grid, Rx, Ry, Rz, dx, dy, dz);
		break;
	}
	case 7:
	{
		Ylm<2, 1> Y;
		SetBlockYlmAroundAny(Y, 0.0, buf, range_blocks, shifted_grid, Rx, Ry, Rz, dx, dy, dz);
		break;
	}
	case 8:
	{
		Ylm<2, 2> Y;
		SetBlockYlmAroundAny(Y, 0.0, buf, range_blocks, shifted_grid, Rx, Ry, Rz, dx, dy, dz);
		break;
	}
	case 9:
	{
		Ylm<3, -3> Y;
		SetBlockYlmAroundAny(Y, 0.0, buf, range_blocks, shifted_grid, Rx, Ry, Rz, dx, dy, dz);
		break;
	}
	case 10:
	{
		Ylm<3, -2> Y;
		SetBlockYlmAroundAny(Y, 0.0, buf, range_blocks, shifted_grid, Rx, Ry, Rz, dx, dy, dz);
		break;
	}
	case 11:
	{
		Ylm<3, -1> Y;
		SetBlockYlmAroundAny(Y, 0.0, buf, range_blocks, shifted_grid, Rx, Ry, Rz, dx, dy, dz);
		break;
	}
	case 12:
	{
		Ylm<3, 0> Y;
		SetBlockYlmAroundAny(Y, 0.0, buf, range_blocks, shifted_grid, Rx, Ry, Rz, dx, dy, dz);
		break;
	}
	case 13:
	{
		Ylm<3, 1> Y;
		SetBlockYlmAroundAny(Y, 0.0, buf, range_blocks, shifted_grid, Rx, Ry, Rz, dx, dy, dz);
		break;
	}
	case 14:
	{
		Ylm<3, 2> Y;
		SetBlockYlmAroundAny(Y, 0.0, buf, range_blocks, shifted_grid, Rx, Ry, Rz, dx, dy, dz);
		break;
	}
	case 15:
	{
		Ylm<3, 3> Y;
		SetBlockYlmAroundAny(Y, 0.0, buf, range_blocks, shifted_grid, Rx, Ry, Rz, dx, dy, dz);
		break;
	}
	default:
	{
		YlmSimple Y(L, M);
		SetBlockYlmAroundAny(Y, 0.0, buf, range_blocks, shifted_grid, Rx, Ry, Rz, dx, dy, dz);
		break;
	}
	}
}

template<class YLM>
void SetVsYlmAroundAny(YLM& Y, double value_org, double* buf, const std::vector<int>& voxel_lines,
	const std::vector<ShiftedPos>& shifted_pos,
	double dx, double dy, double dz)
{


	constexpr double rr_min = 1.0e-14;

	int total_lines = 0;
	int total_i = 0;
	for (const auto& pos : shifted_pos) {
		for (int il = 0; il < pos.num_lines; ++il) {
			const int index = (il + total_lines) * 4;
			const int iz = voxel_lines[index];
			const int iy = voxel_lines[index + 1];
			const int ix_begin = voxel_lines[index + 2];
			const int ix_end = ix_begin + voxel_lines[index + 3];

			const double z = dz * (double)iz - pos.z;
			const double y = dy * (double)iy - pos.y;
			const double yy_zz = y * y + z * z;
			for (int ix = ix_begin; ix < ix_end; ++ix) {
				const double x = dx * (double)ix - pos.x;
				const double rr = x * x + yy_zz;

				if (rr_min > rr) {
					buf[total_i + ix - ix_begin] = value_org;
				} else {
					const double r = sqrt(rr);
					buf[total_i + ix - ix_begin] = Y(x / r, y / r, z / r);
				}
			}
			total_i += voxel_lines[index + 3];
		}
		total_lines += pos.num_lines;
	}

}

inline
void SetVsYlmAround(int L, int M, double* buf, const std::vector<int>& voxel_lines,
	const std::vector<ShiftedPos>& shifted_pos,
	double dx, double dy, double dz)
{
	const double Y00 = 1.0 / sqrt(4.0 * M_PI);
	switch (L * L + L + M) {
	case 0:
	{
		Ylm<0, 0> Y;
		SetVsYlmAroundAny(Y, Y00, buf, voxel_lines, shifted_pos, dx, dy, dz);
		break;
	}
	case 1:
	{
		Ylm<1, -1> Y;
		SetVsYlmAroundAny(Y, 0.0, buf, voxel_lines, shifted_pos, dx, dy, dz);
		break;
	}
	case 2:
	{
		Ylm<1, 0> Y;
		SetVsYlmAroundAny(Y, 0.0, buf, voxel_lines, shifted_pos, dx, dy, dz);
		break;
	}
	case 3:
	{
		Ylm<1, 1> Y;
		SetVsYlmAroundAny(Y, 0.0, buf, voxel_lines, shifted_pos, dx, dy, dz);
		break;
	}
	case 4:
	{
		Ylm<2, -2> Y;
		SetVsYlmAroundAny(Y, 0.0, buf, voxel_lines, shifted_pos, dx, dy, dz);
		break;
	}
	case 5:
	{
		Ylm<2, -1> Y;
		SetVsYlmAroundAny(Y, 0.0, buf, voxel_lines, shifted_pos, dx, dy, dz);
		break;
	}
	case 6:
	{
		Ylm<2, 0> Y;
		SetVsYlmAroundAny(Y, 0.0, buf, voxel_lines, shifted_pos, dx, dy, dz);
		break;
	}
	case 7:
	{
		Ylm<2, 1> Y;
		SetVsYlmAroundAny(Y, 0.0, buf, voxel_lines, shifted_pos, dx, dy, dz);
		break;
	}
	case 8:
	{
		Ylm<2, 2> Y;
		SetVsYlmAroundAny(Y, 0.0, buf, voxel_lines, shifted_pos, dx, dy, dz);
		break;
	}
	case 9:
	{
		Ylm<3, -3> Y;
		SetVsYlmAroundAny(Y, 0.0, buf, voxel_lines, shifted_pos, dx, dy, dz);
		break;
	}
	case 10:
	{
		Ylm<3, -2> Y;
		SetVsYlmAroundAny(Y, 0.0, buf, voxel_lines, shifted_pos, dx, dy, dz);
		break;
	}
	case 11:
	{
		Ylm<3, -1> Y;
		SetVsYlmAroundAny(Y, 0.0, buf, voxel_lines, shifted_pos, dx, dy, dz);
		break;
	}
	case 12:
	{
		Ylm<3, 0> Y;
		SetVsYlmAroundAny(Y, 0.0, buf, voxel_lines, shifted_pos, dx, dy, dz);
		break;
	}
	case 13:
	{
		Ylm<3, 1> Y;
		SetVsYlmAroundAny(Y, 0.0, buf, voxel_lines, shifted_pos, dx, dy, dz);
		break;
	}
	case 14:
	{
		Ylm<3, 2> Y;
		SetVsYlmAroundAny(Y, 0.0, buf, voxel_lines, shifted_pos, dx, dy, dz);
		break;
	}
	case 15:
	{
		Ylm<3, 3> Y;
		SetVsYlmAroundAny(Y, 0.0, buf, voxel_lines, shifted_pos, dx, dy, dz);
		break;
	}
	default:
	{
		YlmSimple Y(L, M);
		SetVsYlmAroundAny(Y, 0.0, buf, voxel_lines, shifted_pos, dx, dy, dz);
		break;
	}
	}
}

