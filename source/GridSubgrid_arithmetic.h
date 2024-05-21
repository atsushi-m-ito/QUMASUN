#pragma once
#include "GridRange.h"
#include "GridSubgrid.h"
#include "SubspaceField.h"
#include "soacomplex.h"



/*
* calculate <p(x,y,z)|psi(x,y,z)>
* <p| is projector and <p| = <w(r)Y(theta,phi)|
* In addition, radial function w(r) is converted to
* the data w(x,y,z) in real subspace grid.
* The spherical harminic function Y(theta,phi) is Y_{00} = 1
*/
inline
double InnerProd_1S_1R(const SubspaceField& w, const double* l_psi, const GridRangeMPI& grid, int Nx, int Ny, int Nz) {

	double sum = 0.0;
	const double* wbuf = w.GetGridPointer();

#if 0
	const int proc_id = GetProcessID(grid.mpi_comm);
	printf("[%d] Nx, Ny, Nz = %d, %d, %d\n", proc_id, Nx, Ny, Nz); fflush(stdout);
	printf("[%d] grid1: %d, %d, %d, %d, %d, %d\n", proc_id, grid.begin_x, grid.begin_y, grid.begin_z, grid.end_x, grid.end_y, grid.end_z); fflush(stdout);
	printf("[%d] grid2: %d, %d, %d, %d, %d, %d\n", proc_id, w.begin_grid_x, w.begin_grid_y, w.begin_grid_z, w.begin_grid_x + w.num_grid_x, w.begin_grid_y + w.num_grid_y, w.begin_grid_z + w.num_grid_z); fflush(stdout);
#endif

	ForOverlapPeriodic(
		*(GridRange*)&grid, 
		GridRange{ w.begin_grid_x, w.begin_grid_y, w.begin_grid_z, w.begin_grid_x + w.num_grid_x, w.begin_grid_y + w.num_grid_y, w.begin_grid_z + w.num_grid_z },
		Nx, Ny, Nz, 
		[&sum, &wbuf, &l_psi](size_t whole_i, size_t sub_k) {
			sum += wbuf[sub_k] * l_psi[whole_i];
		});

	return sum;
}

/*
* calculate <p(x,y,z)|f(x,y,z) psi(x,y,z)>
* Here, f(x,y,z) and psi(x,y,z) is Complex number.
* And, f imply exp(ikx) in Bloch theorem
* <p| is projector and <p| = <w(r)Y(theta,phi)|
* In addition, radial function w(r) is converted to
* the data w(x,y,z) in real subspace grid.
* The spherical harminic function Y(theta,phi) is Y_{00} = 1
*
*/
inline
OneComplex InnerProd_1SZ_1R(const SubspaceField& w, const double* f_re, const double* f_im, const double* l_psi_re, const double* l_psi_im, const GridRangeMPI& grid, int Nx, int Ny, int Nz) {

	OneComplex sum{ 0.0, 0.0 };
	const double* wbuf = w.GetGridPointer();


	ForOverlapPeriodic(
		*(GridRange*)&grid,
		GridRange{ w.begin_grid_x, w.begin_grid_y, w.begin_grid_z, w.begin_grid_x + w.num_grid_x, w.begin_grid_y + w.num_grid_y, w.begin_grid_z + w.num_grid_z },
		Nx, Ny, Nz,
		[&sum, &wbuf, &f_re, &f_im, &l_psi_re, &l_psi_im](size_t whole_i, size_t sub_k) {
			sum.r += wbuf[sub_k] * (f_re[sub_k] * l_psi_re[whole_i] - f_im[sub_k] * l_psi_im[whole_i]);
			sum.i += wbuf[sub_k] * (f_im[sub_k] * l_psi_re[whole_i] + f_re[sub_k] * l_psi_im[whole_i]);
		});

	return sum;
}


/*
* calculate <p(x,y,z)|psi(x,y,z)>
* <p| is projector and <p| = <w(r)Y(theta,phi)|
* In addition, radial function w(r) and spherical harminic function Y(theta,phi) are converted to
* the data w(x,y,z) and Y(x,y,z) in real subspace grid.
*/
inline
double InnerProd_2S_1R(const SubspaceField& w, const SubspaceField& Y, const double* l_psi, const GridRangeMPI& grid, int Nx, int Ny, int Nz) {

	double sum = 0.0;
	const double* wbuf = w.GetGridPointer();
	const double* ybuf = Y.GetGridPointer();

	ForOverlapPeriodic(
		*(GridRange*)&grid,
		GridRange{ w.begin_grid_x, w.begin_grid_y, w.begin_grid_z, w.begin_grid_x + w.num_grid_x, w.begin_grid_y + w.num_grid_y, w.begin_grid_z + w.num_grid_z },
		Nx, Ny, Nz,
		[&sum, &wbuf, &ybuf, &l_psi](size_t whole_i, size_t sub_k) {
			sum += wbuf[sub_k] * ybuf[sub_k] * l_psi[whole_i];
		});

	return sum;
}

inline
OneComplex InnerProd_2SZ_1R(const SubspaceField& w, const SubspaceField& Y, const double* f_re, const double* f_im, const double* l_psi_re, const double* l_psi_im, const GridRangeMPI& grid, int Nx, int Ny, int Nz) {

	OneComplex sum{ 0.0,0.0 };
	const double* wbuf = w.GetGridPointer();
	const double* ybuf = Y.GetGridPointer();

	ForOverlapPeriodic(
		*(GridRange*)&grid,
		GridRange{ w.begin_grid_x, w.begin_grid_y, w.begin_grid_z, w.begin_grid_x + w.num_grid_x, w.begin_grid_y + w.num_grid_y, w.begin_grid_z + w.num_grid_z },
		Nx, Ny, Nz,
		[&sum, &wbuf, &ybuf, &l_psi_re, &l_psi_im, &f_re, &f_im](size_t whole_i, size_t sub_k) {
			sum.r += wbuf[sub_k] * ybuf[sub_k] * (f_re[sub_k] * l_psi_re[whole_i] - f_im[sub_k] * l_psi_im[whole_i]);
			sum.i += wbuf[sub_k] * ybuf[sub_k] * (f_im[sub_k] * l_psi_re[whole_i] + f_re[sub_k] * l_psi_im[whole_i]);
		});

	return sum;
}


inline
void Add_1S(double* v, const GridRangeMPI& grid, const SubspaceField& w, double coef, int Nx, int Ny, int Nz) {


	const double* wbuf = w.GetGridPointer();
	
	ForOverlapPeriodic(
		*(GridRange*)&grid,
		GridRange{ w.begin_grid_x, w.begin_grid_y, w.begin_grid_z, w.begin_grid_x + w.num_grid_x, w.begin_grid_y + w.num_grid_y, w.begin_grid_z + w.num_grid_z },
		Nx, Ny, Nz,
		[coef, &v, &wbuf](size_t whole_i, size_t sub_k) {
			v[whole_i] += coef * wbuf[sub_k];
		});
}



inline
void Add_1SZ(double* v_re, double* v_im, const GridRangeMPI& grid, const SubspaceField& w, const double* f_re, const double* f_im, double coef_re, double coef_im, int Nx, int Ny, int Nz) {


	const double* wbuf = w.GetGridPointer();
	

	ForOverlapPeriodic(
		*(GridRange*)&grid,
		GridRange{ w.begin_grid_x, w.begin_grid_y, w.begin_grid_z, w.begin_grid_x + w.num_grid_x, w.begin_grid_y + w.num_grid_y, w.begin_grid_z + w.num_grid_z },
		Nx, Ny, Nz,
		[coef_re, coef_im, &v_re, &v_im, &wbuf, &f_re, &f_im](size_t whole_i, size_t sub_k) {
			v_re[whole_i] += wbuf[sub_k] * (f_re[sub_k] * coef_re - f_im[sub_k] * coef_im);
			v_im[whole_i] += wbuf[sub_k] * (f_im[sub_k] * coef_re + f_re[sub_k] * coef_im);
		});
}


inline
void Add_2S(double* v, const GridRangeMPI& grid, const SubspaceField& w, const SubspaceField& Y, double coef, int Nx, int Ny, int Nz) {


	const double* wbuf = w.GetGridPointer();
	const double* ybuf = Y.GetGridPointer();


	ForOverlapPeriodic(
		*(GridRange*)&grid,
		GridRange{ w.begin_grid_x, w.begin_grid_y, w.begin_grid_z, w.begin_grid_x + w.num_grid_x, w.begin_grid_y + w.num_grid_y, w.begin_grid_z + w.num_grid_z },
		Nx, Ny, Nz,
		[coef, &v, &wbuf, &ybuf](size_t whole_i, size_t sub_k) {
			v[whole_i] += coef * wbuf[sub_k] * ybuf[sub_k];
		});

}

inline
void Add_2SZ(double* v_re, double* v_im, const GridRangeMPI& grid, const SubspaceField& w, const SubspaceField& Y, const double* f_re, const double* f_im, double coef_re, double coef_im, int Nx, int Ny, int Nz) {


	const double* wbuf = w.GetGridPointer();
	const double* ybuf = Y.GetGridPointer();

	ForOverlapPeriodic(
		*(GridRange*)&grid,
		GridRange{ w.begin_grid_x, w.begin_grid_y, w.begin_grid_z, w.begin_grid_x + w.num_grid_x, w.begin_grid_y + w.num_grid_y, w.begin_grid_z + w.num_grid_z },
		Nx, Ny, Nz,
		[coef_re, coef_im, &v_re, &v_im, &wbuf, &ybuf, &f_re, &f_im](size_t whole_i, size_t sub_k) {
			v_re[whole_i] += wbuf[sub_k] * ybuf[sub_k] * (f_re[sub_k] * coef_re - f_im[sub_k] * coef_im);
			v_im[whole_i] += wbuf[sub_k] * ybuf[sub_k] * (f_im[sub_k] * coef_re + f_re[sub_k] * coef_im);
		});

}


inline
void Add_1ScZ(double* v_re, double* v_im, const GridRangeMPI& grid, const SubspaceField& w, const double* f_re, const double* f_im, double coef_re, double coef_im, int Nx, int Ny, int Nz) {


	const double* wbuf = w.GetGridPointer();

	ForOverlapPeriodic(
		*(GridRange*)&grid,
		GridRange{ w.begin_grid_x, w.begin_grid_y, w.begin_grid_z, w.begin_grid_x + w.num_grid_x, w.begin_grid_y + w.num_grid_y, w.begin_grid_z + w.num_grid_z },
		Nx, Ny, Nz,
		[coef_re, coef_im, &v_re, &v_im, &wbuf, &f_re, &f_im](size_t whole_i, size_t sub_k) {
			v_re[whole_i] += wbuf[sub_k] * (f_re[sub_k] * coef_re + f_im[sub_k] * coef_im);
			v_im[whole_i] += wbuf[sub_k] * (-f_im[sub_k] * coef_re + f_re[sub_k] * coef_im);
		});

}
