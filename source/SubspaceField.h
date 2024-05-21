#pragma once


#include <cmath>
#include <cstdio>
#include "GridSubgrid.h"


/*
動径方向の座標において、グリッド上の半径rの値と、グリッド間隔等に関する情報を保持するクラス

他のクラスに参照させて、メモリは共有で使いまわす
*/
class SubspaceField{
private:
	double* buffer = nullptr;
public:
	int begin_grid_x = 0;
	int begin_grid_y = 0;
	int begin_grid_z = 0;
	int num_grid_x = 0;
	int num_grid_y = 0;
	int num_grid_z = 0;
	size_t num_total_grid = 0;
	
	
	SubspaceField(int begin_grid_x_, int begin_grid_y_, int begin_grid_z_,
		int num_grid_x_, int num_grid_y_, int num_grid_z_) :
		begin_grid_x(begin_grid_x_),
		begin_grid_y(begin_grid_y_),
		begin_grid_z(begin_grid_z_),
		num_grid_x(num_grid_x_),
		num_grid_y(num_grid_y_),
		num_grid_z(num_grid_z_),
		num_total_grid((size_t)num_grid_x_ * (size_t)num_grid_y_ * (size_t)num_grid_z_)
	{
	};


	double* Reserve() {
		buffer = new double[num_total_grid];
		return buffer;
	}

	void Release() {
		delete[] buffer;
	}

	double* GetGridPointer() {
		return buffer;
	}

	const double* GetGridPointer() const {
		return buffer;
	}
};



/*
* subspaceとglobal spaceのオーバーラップ部分の格子点に対して任意の関数FUNCを実行しながらループする
* 前提として周期境界になっている
* 領域分割MPI並列ではForOverlapPeriodicを利用すること
*/
#if 1
template <class FUNC>
inline void OverlapLoop_periodic(int sub_begin_x, int sub_begin_y, int sub_begin_z,
	int sub_size_x, int sub_size_y, int sub_size_z,
	int whole_x, int whole_y, int whole_z, FUNC&& func) {

	ForOverlapPeriodic(
		GridRange{ 0, 0, 0, whole_x, whole_y, whole_z },
		GridRange{ sub_begin_x, sub_begin_y, sub_begin_z, sub_begin_x + sub_size_x, sub_begin_y + sub_size_y, sub_begin_z + sub_size_z },
		whole_x, whole_y, whole_z, func);
}

#else
template <class FUNC>
inline void OverlapLoop_periodic(int sub_begin_x, int sub_begin_y, int sub_begin_z,
	int sub_size_x, int sub_size_y, int sub_size_z,
	int whole_x, int whole_y, int whole_z, FUNC ope) {


	auto imodpositive = [](int x, int N) {
		return (x >= 0) ? x % N : (N - ((-x - 1) % N + 1));
		};

	int ib_x = imodpositive(sub_begin_x, whole_x);
	int ib_y = imodpositive(sub_begin_y, whole_y);
	int ib_z = imodpositive(sub_begin_z, whole_z);

	for (int iz0 = ib_z, offset_z = 0; offset_z < sub_size_z; iz0 = 0) {
		//subspaceがglobal spaceを周期境界で跨ぐ場合には2度目のループへ入る//

		const int iz1 = std::min<int>(iz0 + sub_size_z - offset_z, whole_z);
		//printf("iz_begin-end = %d, %d\n", iz0, iz1);

		for (int iz = iz0; iz < iz1; ++iz) {
			const size_t kz = iz - iz0 + offset_z;

			for (int iy0 = ib_y, offset_y = 0; offset_y < sub_size_y; iy0 = 0) {
				//subspaceがglobal spaceを周期境界で跨ぐ場合には2度目のループへ入る//

				const int iy1 = std::min<int>(iy0 + sub_size_y - offset_y, whole_y);
				for (int iy = iy0; iy < iy1; ++iy) {
					const size_t ky = iy - iy0 + offset_y;


					for (int ix0 = ib_x, offset_x = 0; offset_x < sub_size_x; ix0 = 0) {
						//subspaceがglobal spaceを周期境界で跨ぐ場合には2度目のループへ入る//

						const int ix1 = std::min<int>(ix0 + sub_size_x - offset_x, whole_x);
						for (int ix = ix0; ix < ix1; ++ix) {
							const size_t i = ix + whole_x * (iy + whole_y * iz);
							const size_t kx = ix - ix0 + offset_x;
							const size_t k = kx + sub_size_x * (ky + sub_size_y * kz);
							ope(i, k);
						}
						offset_x += ix1 - ix0;
					}
				}
				offset_y += iy1 - iy0;
			}
		}
		offset_z += iz1 - iz0;
	}

}
#endif

/*
* calculate <p(x,y,z)|psi(x,y,z)>
* <p| is projector and <p| = <w(r)Y(theta,phi)|
* In addition, radial function w(r) and spherical harminic function Y(theta,phi) are converted to
* the data w(x,y,z) and Y(x,y,z) in real subspace grid.
*/
inline
double InnerProd_2S_1R(const SubspaceField& w, const SubspaceField& Y, const RspaceFunc<double>& psi, int Nx, int Ny, int Nz) {

	double sum = 0.0;
	const double* wbuf = w.GetGridPointer();
	const double* ybuf = Y.GetGridPointer();
	
	OverlapLoop_periodic(w.begin_grid_x, w.begin_grid_y, w.begin_grid_z, w.num_grid_x, w.num_grid_y, w.num_grid_z,
		Nx, Ny, Nz,
		[&sum, &wbuf, &ybuf, &psi](size_t whole_i, size_t sub_k) {
			sum += wbuf[sub_k] * ybuf[sub_k] * psi[whole_i];
		});
		
	/*
	ForOverlapPeriodic(
		GridRange{ 0, 0, 0, Nx, Ny, Nz },
		GridRange{ w.begin_grid_x, w.begin_grid_y, w.begin_grid_z,
		w.begin_grid_x + w.num_grid_x, w.begin_grid_y + w.num_grid_y, w.begin_grid_z + w.num_grid_z },		
		Nx, Ny, Nz,
		[&sum, &wbuf, &ybuf, &psi](size_t whole_i, size_t sub_k) {
			sum += wbuf[sub_k] * ybuf[sub_k] * psi[whole_i];
		});*/
	return sum;
}

/*
* calculate <p(x,y,z)|psi(x,y,z)>
* <p| is projector and <p| = <w(r)Y(theta,phi)|
* In addition, radial function w(r) is converted to
* the data w(x,y,z) in real subspace grid.
* The spherical harminic function Y(theta,phi) is Y_{00} = 1
*/
#if 1
inline
double InnerProd_1S_1R(const SubspaceField& w, const RspaceFunc<double>& psi, int Nx, int Ny, int Nz) {

	double sum = 0.0;
	const double* wbuf = w.GetGridPointer();

	OverlapLoop_periodic(w.begin_grid_x, w.begin_grid_y, w.begin_grid_z, w.num_grid_x, w.num_grid_y, w.num_grid_z,
		Nx, Ny, Nz,
		[&sum, &wbuf, &psi](size_t whole_i, size_t sub_k) {
			sum += wbuf[sub_k] * psi[whole_i];
		});
	return sum;
}

#else

inline
double InnerProd_1S_1R(const SubspaceField & w, const RspaceFunc<double>&psi, int Nx, int Ny, int Nz) {

	int ib_x = imodpositive(w.begin_grid_x, Nx);
	int ib_y = imodpositive(w.begin_grid_y, Ny);
	int ib_z = imodpositive(w.begin_grid_z, Nz);
	
	double sum = 0.0;
	const double* wbuf = w.GetGridPointer();
	

	for (int iz0 = ib_z, offset_z = 0; offset_z < w.num_grid_z; iz0 = 0) {
		const int iz1 = std::min<int>(iz0 + w.num_grid_z - offset_z, Nz);
		for (int iz = iz0; iz < iz1; ++iz) {
			size_t kz = iz - iz0 + offset_z;
		
			for (int iy0 = ib_y, offset_y = 0; offset_y < w.num_grid_y; iy0 = 0) {
				const int iy1 = std::min<int>(iy0 + w.num_grid_y - offset_y, Ny);
				for (int iy = iy0; iy < iy1; ++iy) {
					size_t ky = iy - iy0 + offset_y;


					for (int ix0 = ib_x, offset_x = 0; offset_x < w.num_grid_x; ix0 = 0) {
						const int ix1 = std::min<int>(ix0 + w.num_grid_x - offset_x, Nx);
						for (int ix = ix0; ix < ix1; ++ix) {
							size_t i = ix + Nx * (iy + Ny * iz);
							size_t kx = ix - ix0 + offset_x;
							size_t k = kx + w.num_grid_x * (ky + w.num_grid_y * kz);
							sum += wbuf[k] * psi[i];
						}
						offset_x += ix1 - ix0;
					}
				}
				offset_y += iy1 - iy0;
			}			
		}
		offset_z += iz1 - iz0;
	}
	return sum;
}
#endif

inline
void Add_1S(RspaceFunc<double>& v, const SubspaceField& w, double coef, int Nx, int Ny, int Nz) {

	
	const double* wbuf = w.GetGridPointer();

	OverlapLoop_periodic(w.begin_grid_x, w.begin_grid_y, w.begin_grid_z, w.num_grid_x, w.num_grid_y, w.num_grid_z,
		Nx, Ny, Nz,
		[coef, &v, &wbuf](size_t whole_i, size_t sub_k) {
			v[whole_i] += coef * wbuf[sub_k];
		});
	
}

inline
void Add_2S(RspaceFunc<double>& v, const SubspaceField& w, const SubspaceField& Y, double coef, int Nx, int Ny, int Nz) {


	const double* wbuf = w.GetGridPointer();
	const double* ybuf = Y.GetGridPointer();

	OverlapLoop_periodic(w.begin_grid_x, w.begin_grid_y, w.begin_grid_z, w.num_grid_x, w.num_grid_y, w.num_grid_z,
		Nx, Ny, Nz,
		[coef, &v, &wbuf, &ybuf](size_t whole_i, size_t sub_k) {
			v[whole_i] += coef * wbuf[sub_k] * ybuf[sub_k];
		});

}

