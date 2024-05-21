/*******************************

QUantum MAterial Simulation UNraveler (QUMASUN)

QUMASUN is numerical simulation code for Density Functional Theory(DFT) and Time-dependent DFT based on the real space grid.

********************************/
//#define TEST_PRINT
#include "mpi.h"

#include <cstdlib>
#include <cstdio>
#include "wave_function.h"
#include "GridRange.h"
#include "GridRange.h"
#include "GridFor.h"
#include "GetArg.h"
#include "GridScatterGather.h"
#include "GridDifference2nd.h"
#include "GridDifference4th.h"


//test to cut sub firld and to interporate from field
int main(int argc, char* argv[]) {
	MPI_Init(&argc, &argv);
	MPI_Comm mpi_comm = MPI_COMM_WORLD;
	int proc_id;
	int num_procs;
	MPI_Comm_rank(mpi_comm, &proc_id);
	MPI_Comm_size(mpi_comm, &num_procs);

	const int Nx = 256;
	const int Ny = 256;
	const int Nz = 256;
	double dx = 20.0 / Nx;
	double dy = 20.0 / Ny;
	double dz = 20.0 / Nz;

	size_t global_size_3d = (size_t)Nx * (size_t)Ny * (size_t)Nz;
	GridRange global_grid{ 0,0,0,Nx,Ny,Nz };

	int split_num[3] = { 1,1,1 };
	int res = GetArgumentNumList<int>(argc, argv, "-ddm", 3, split_num);

	if (IsRoot(mpi_comm)) {
		printf("split grid by %d x %d x %d\n", split_num[0], split_num[1], split_num[2]);
		fflush(stdout);
	}
	MPI_Barrier(mpi_comm);
	GridRangeMPI grid = MakeRange(Nx, Ny, Nz, mpi_comm, split_num);
	printf("[%d] grid_size_3d = %zd / %zd\n", proc_id, grid.Size3D(), global_grid.Size3D()); fflush(stdout);

	double* gphi = IsRoot(mpi_comm) ? new double[global_grid.Size3D()] : nullptr;
	double* gd2phi_dx2 = IsRoot(mpi_comm) ? new double[global_grid.Size3D()] : nullptr;
	double* gd2phi_dx2_4th = IsRoot(mpi_comm) ? new double[global_grid.Size3D()] : nullptr;
	double* b_gath = IsRoot(mpi_comm) ? new double[global_grid.Size3D()] : nullptr;
	double* a = new double[grid.Size3D()];
	double* d2a_dx2_ref = new double[grid.Size3D()];
	double* d2a_dx2_ref4th = new double[grid.Size3D()];
	double* b = new double[grid.Size3D()];
	double* c = new double[grid.Size3D()];
	double* c_gath = IsRoot(mpi_comm) ? new double[global_grid.Size3D()] : nullptr;

	
	size_t min_a = INT_MAX;
	size_t max_a = 0;

	{//data initialize//
		double norm_g = 0.0;
		if (IsRoot(mpi_comm)) {
			ForXYZ(global_grid, [&](int64_t i, int64_t gix, int64_t giy, int64_t giz) {
				const double val = (double)(gix - Nx / 2) * (double)(gix - Nx / 2) * dx * dx / 2.0;
				+(double)(giy - Ny / 2) * (double)(giy - Ny / 2) * dy * dy / 2.0
					+ (double)(giz - Nz / 2) * (double)(giz - Nz / 2) * dz * dz / 2.0;
				gphi[i] = exp(-val);
				//gphi[i] = (double)(giz);
				norm_g += gphi[i]* gphi[i];
				gd2phi_dx2[i]=0.0;
				gd2phi_dx2_4th[i] = 0.0;
				});

			norm_g *= dx * dy * dz;
			const double coef = 1.0 / sqrt(norm_g);
			norm_g = 0.0;
			For(global_grid, [&](int64_t i) {
				gphi[i] *= coef;
				norm_g += gphi[i] * gphi[i];
				});

			norm_g *= dx * dy * dz;
			printf("[%d]initialize: norm=%f\n", proc_id, norm_g); fflush(stdout);

			Laplasian2nd(global_grid, gd2phi_dx2, gphi, -1.0 / 2.0, dx, dy, dz);
			Laplasian4th(global_grid, gd2phi_dx2_4th, gphi, -1.0 / 2.0, dx, dy, dz);
			
			
		}

		For(grid, [&](int64_t i) {
			b[i] = 0.0;
			c[i] = 0.0;
			});


		ScatterGrid(grid, a, global_grid, gphi, 0);
		ScatterGrid(grid, d2a_dx2_ref, global_grid, gd2phi_dx2, 0);
		ScatterGrid(grid, d2a_dx2_ref4th, global_grid, gd2phi_dx2_4th, 0);
		
		{
			double sum_a = 0.0;
			For(grid, [&](int64_t i) {
				sum_a += a[i];
			});
			printf("[%d]local_sum_a = %f\n", proc_id, sum_a); fflush(stdout);
		}

		double norm = ForReduce2(grid, 0.0, [&](int64_t i) {
			return a[i]* a[i];
			});
		norm *= dx * dy * dz;
		
		/*
		For(grid, [&](int64_t i) {
			return a[i] /= sum2;
			});
			*/
		if (IsRoot(mpi_comm)) {
			if (fabs(norm_g - norm) < 1.0e-7) {
				printf("[%d]sum = %f == %f\n", proc_id, norm_g, norm); fflush(stdout);
			} else {
				printf("ERROR: [%d]sum = %f != %f\n", proc_id, norm_g, norm); fflush(stdout);
			}
		}
	}

	printf("[%d]local_size: %d, %d, %d\n", proc_id, grid.SizeX(), grid.SizeY(), grid.SizeZ()); fflush(stdout);
	Laplasian2nd(grid, b, a, -1.0 / 2.0, dx, dy, dz);
	Laplasian4th(grid, c, a, -1.0 / 2.0, dx, dy, dz);

	GatherGrid(global_grid, b_gath, grid, b, 0);
	GatherGrid(global_grid, c_gath, grid, c, 0);
	

	{
		
		double sum_bc_hand[2]{ 0.0, 0.0 };
		{
			double local_b_hand[2]{ 0.0, 0.0 };
			const int64_t i_end = grid.Size3D();
			for (int64_t i = 0; i < i_end; ++i) {
				local_b_hand[0] += a[i]*b[i];
				local_b_hand[1] += a[i] * c[i];
			}
			local_b_hand[0] *= dx * dy * dz;
			local_b_hand[1] *= dx * dy * dz;
			printf("[%d]local_b = %.15g\n", proc_id, local_b_hand[0]); fflush(stdout);
			printf("[%d]local_c = %.15g\n", proc_id, local_b_hand[1]); fflush(stdout);
			MPI_Reduce(local_b_hand, sum_bc_hand, 2, MPI_DOUBLE_PRECISION, MPI_SUM, 0, grid.mpi_comm);
		}
		double sum_b = ForReduce2(grid, 0.0, [&](int64_t i) {
			return a[i] * b[i];
			});
		sum_b *= dx * dy * dz;

		double sum_c = ForReduce2(grid, 0.0, [&](int64_t i) {
			return a[i] * c[i];
			});
		sum_c *= dx * dy * dz;

		//printf("[%d]test6: %f\n", proc_id, sum_b); fflush(stdout);

		if (IsRoot(mpi_comm)) {
			double sum_a_ref = 0.0;
			double sum_a2_ref = 0.0;
			double sum_b2_ref = 0.0;
			double sum_b_ref = 0.0;
			double sum_b_gath = 0.0;
			double sum_c_ref = 0.0;
			double sum_c_gath = 0.0;
			For(global_grid, [&](int64_t i) {
				sum_b_ref += gphi[i] * gd2phi_dx2[i];
				sum_b2_ref += gd2phi_dx2[i] * gd2phi_dx2[i];
				sum_b_gath += gphi[i] * b_gath[i];

				sum_a_ref += gphi[i];
				sum_a2_ref += gphi[i]* gphi[i];

				sum_c_ref += gphi[i] * gd2phi_dx2_4th[i];
				sum_c_gath += gphi[i]  * c_gath[i];
			});
			sum_b_ref *= dx * dy * dz;
			sum_b2_ref *= dx * dy * dz;
			sum_b_gath *= dx * dy * dz;

			sum_a_ref *= dx * dy * dz;
			sum_a2_ref *= dx * dy * dz;

			sum_c_ref *= dx * dy * dz;
			sum_c_gath *= dx * dy * dz;

			//printf("[%d]test7: %f\n", proc_id, sum_b_ref); fflush(stdout);

			if (fabs(sum_b - sum_b_ref) < 1.0e-10) {
				printf("[%d]out = %g ~ %g, %g, %g\n", proc_id, sum_b_ref, sum_b, sum_b_gath, sum_bc_hand[0]); fflush(stdout);
			} else {
				printf("[%d]sum = %g != %g, %g, %g\n", proc_id, sum_b_ref, sum_b, sum_b_gath, sum_bc_hand[0]); fflush(stdout);
			}
			printf("[%d]b^2out = %g\n", proc_id, sum_b2_ref); fflush(stdout);
			printf("[%d]a, a^2 = %g, %g\n", proc_id, sum_a_ref, sum_a2_ref); fflush(stdout);
			
			if (fabs(sum_c - sum_c_ref) < 1.0e-10) {
				printf("[%d]out = %g ~ %g, %g, %g\n", proc_id, sum_c_ref, sum_c, sum_c_gath, sum_bc_hand[1]); fflush(stdout);
			} else {
				printf("ERROR: [%d]sum = %g != %g, %g, %g\n", proc_id, sum_c_ref, sum_c, sum_c_gath, sum_bc_hand[1]); fflush(stdout);
			}
			
		}
	}



	delete[]gphi;
	delete[]gd2phi_dx2;
	delete[]a;
	delete[]d2a_dx2_ref;
	delete[]b;
	delete[]b_gath;
	delete[] gd2phi_dx2_4th;
	delete[] c;
	delete[] c_gath;

	MPI_Finalize();
}


