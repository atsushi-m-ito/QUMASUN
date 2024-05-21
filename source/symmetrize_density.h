#pragma once
#include "qumasun_input.h"

inline 
void SymmetrizeDensity(double* rho, int size_x, int size_y, int size_z, uint32_t kpoint_symmetry) {
	using namespace QUMASUN::KPOINT_SYMMETRY;
	/*
	k=222でNEGAPOSIの処理を入れると,local minimumに引っ掛かる//
	if (kpoint_symmetry & NEGAPOSI_X) {
		for (int iz = 0; iz < size_z; ++iz) {
			for (int iy = 0; iy < size_y; ++iy) {
				for (int ix = 1; ix < (size_x+1)/2; ++ix) {
					const int64_t i1 = ix + size_x * (iy + size_y * iz);
					const int64_t i2 = (size_x - ix) + size_x * (iy + size_y * iz);
					double ave = (rho[i1] + rho[i2])/2.0;
					rho[i1] = ave;
					rho[i2] = ave;
				}
			}
		}
	}

	if (kpoint_symmetry & NEGAPOSI_Y) {
		for (int iz = 0; iz < size_z; ++iz) {
			for (int iy = 1; iy < (size_y+1)/2; ++iy) {
				for (int ix = 0; ix < size_x; ++ix) {
					const int64_t i1 = ix + size_x * (iy + size_y * iz);
					const int64_t i2 = ix + size_x * ((size_y - iy) + size_y * iz);
					double ave = (rho[i1] + rho[i2]) / 2.0;
					rho[i1] = ave;
					rho[i2] = ave;
				}
			}
		}
	}

	if (kpoint_symmetry & NEGAPOSI_Z) {
		for (int iz = 1; iz < (size_z + 1) / 2; ++iz) {
			for (int iy = 0; iy < size_y; ++iy) {
				for (int ix = 0; ix < size_x; ++ix) {
					const int64_t i1 = ix + size_x * (iy + size_y * iz);
					const int64_t i2 = ix + size_x * (iy + size_y * (size_z - iz));
					double ave = (rho[i1] + rho[i2]) / 2.0;
					rho[i1] = ave;
					rho[i2] = ave;
				}
			}
		}
	}
	*/
	if (kpoint_symmetry & MIRROR_XY) {
		for (int iz = 0; iz < size_z; ++iz) {
			for (int iy = 0; iy < size_y; ++iy) {
				for (int ix = iy+1; ix < size_x; ++ix) {
					const int64_t i1 = ix + size_x * (iy + size_y * iz);
					const int64_t i2 = iy + size_x * (ix + size_y * iz);
					double ave = (rho[i1] + rho[i2]) / 2.0;
					rho[i1] = ave;
					rho[i2] = ave;
				}
			}
		}
	}

	if (kpoint_symmetry & MIRROR_YZ) {
		for (int iz = 0; iz < size_z; ++iz) {
			for (int iy = iz+1; iy < size_y; ++iy) {
				for (int ix = 0; ix < size_x; ++ix) {
					const int64_t i1 = ix + size_x * (iy + size_y * iz);
					const int64_t i2 = ix + size_x * (iz + size_y * iy);
					double ave = (rho[i1] + rho[i2]) / 2.0;
					rho[i1] = ave;
					rho[i2] = ave;
				}
			}
		}
	}

	if (kpoint_symmetry & MIRROR_ZX) {
		for (int iz = 0; iz < size_z; ++iz) {
			for (int iy = 0; iy < size_y; ++iy) {
				for (int ix = iz + 1; ix < size_x; ++ix) {
					const int64_t i1 = ix + size_x * (iy + size_y * iz);
					const int64_t i2 = iz + size_x * (iy + size_y * ix);
					double ave = (rho[i1] + rho[i2]) / 2.0;
					rho[i1] = ave;
					rho[i2] = ave;
				}
			}
		}
	}

}
