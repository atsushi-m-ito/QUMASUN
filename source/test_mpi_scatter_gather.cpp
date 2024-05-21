/*******************************

QUantum MAterial Simulation UNraveler (QUMASUN)

QUMASUN is numerical simulation code for Density Functional Theory(DFT) and Time-dependent DFT based on the real space grid.

********************************/
#include <mpi.h>

#include <cstdlib>
#include <cstdio>
#include <ccomplex>
#include "mpi_helper.h"
#include "wave_function.h"
#include "GridRange.h"
#include "GridRange.h"
#include "GridFor.h"
#include "GetArg.h"
#include "GridScatterGather.h"



//test to cut sub firld and to interporate from field
int main(int argc, char* argv[]) {
	MPI_Init(&argc, &argv);
	MPI_Comm mpi_comm = MPI_COMM_WORLD;
	const int proc_id = GetProcessID(mpi_comm);
	const int num_procs = GetNumProcess(mpi_comm);

	int size_x = 32;
	int size_y = 32;
	int size_z = 32;
	size_t global_size_3d = (size_t)size_x * (size_t)size_y * (size_t)size_z;
	GridRange global_grid{ 0,0,0,size_x,size_y,size_z };

	int split_num[3] = { 1,1,1 };
	int res = GetArgumentNumList<int>(argc, argv, "-ddm", 3, split_num);

	if (IsRoot(mpi_comm)) {
		printf("split grid by %d x %d x %d\n", split_num[0], split_num[1], split_num[2]);
		fflush(stdout);
	}
	MPI_Barrier(mpi_comm);
	GridRangeMPI grid = MakeRange(size_x, size_y, size_z, mpi_comm, split_num);
	printf("[%d] grid_size_3d = %zd\n", proc_id, grid.Size3D());
	fflush(stdout);

	RspaceFunc<double> a(grid.Size3D());
	RspaceFunc<double> b(grid.Size3D());

	For(grid, [&](size_t i) {
		a[i] = 1.0;
		});

	For(grid, [&](size_t i) {
		b[i] = a[i] * 2.0 + 1.0;
		});

	double sum = 0.0;
	ForReduce(grid, &sum, [&](size_t i) {
		return b[i];
		});

	
	sum /= (double)(global_size_3d);
	if (IsRoot(mpi_comm)) {
		printf("[%d]sum = %f should be 3.0\n", proc_id, sum);
	}

	size_t min_a=INT_MAX;
	size_t max_a = 0;

	ForXYZ(grid, [&](size_t i, size_t gix, size_t giy, size_t giz) {
		a[i] = ((double)gix) * 2.0;
		min_a = std::min(min_a, gix);
		max_a = std::max(max_a, gix);
		});

	double sum2 = ForReduce2(grid, 0.0, [&](size_t i) {
		return a[i];
		});
	sum2 /= (double)(global_size_3d);

	if (IsRoot(mpi_comm)) {
		printf("[%d]sum = %f should be %d, min=%zu, max=%zu\n", proc_id, sum2, size_x - 1, min_a, max_a);
	}

	{//scatter and gather of data grid//
		MPI_Barrier(mpi_comm);

		double* global_psi = nullptr;
		double* global_res = nullptr;
		double res_sum0;
		if (IsRoot(mpi_comm)) {
			global_psi = new double[global_grid.Size3D()];
			global_res = new double[global_grid.Size3D()];
			const int center = size_x / 2;
			ForXYZ(global_grid, [&](int64_t i, int64_t ix, int64_t iy, int64_t iz) {
				const int64_t xx = (ix - center) * (ix - center) + (iy - center) * (iy - center) + (iz - center) * (iz - center);
				global_psi[i] = exp(-(double)xx / (2.0 * ((double)center) * ((double)center)));
				});

			res_sum0 = ForReduce2(global_grid, 0.0, [&](int64_t i) {return global_psi[i]; });
			For(global_grid, [&](int64_t i) {				
				global_psi[i] /= res_sum0;
				//global_psi[i] = 1.0;
				});
			res_sum0 = ForReduce2(global_grid, 0.0, [&](int64_t i) {return global_psi[i]; });
			printf("res_sum0 = %f\n", res_sum0);

		}

		printf("[%d] local_pnt = 0x%zx\n", proc_id, a.Pointer());
		fflush(stdout);


		ScatterGrid(grid, a.Pointer(), global_grid, global_psi, 0);

		{
			double sum = 0.0;
			const int size_local = grid.Size3D();
			for (int i = 0; i < size_local; ++i) {
				sum += a[i];
			}
			printf("[%d] res_sum1_local = %f\n", proc_id, sum); fflush(stdout);
		}
		const double res_sum1 = ForReduce2(grid, 0.0, [&](int64_t i) {return a[i]; });

		GatherGrid(global_grid, global_res, grid, a.Pointer(), 0);

		if (IsRoot(mpi_comm)) {
			printf("res_sum1 = %f\n", res_sum1);
		
			const double res_sum2 = ForReduce2(global_grid, 0.0, [&](int64_t i) {return global_res[i]; });
			printf("res_sum2 = %f\n", res_sum2);

			if (fabs(res_sum1 - res_sum0) > 1.0e-6) {
				printf("ERROR: res_sum1 != res_sum0\n");
			} else {
				printf("OK: res_sum1 - res_sum0 = %f\n", res_sum1 - res_sum0);
			}
			if (fabs(res_sum2 - res_sum0) > 1.0e-6) {
				printf("ERROR: res_sum2 != res_sum0\n");
			} else {
				printf("OK: res_sum2 - res_sum0 = %f\n", res_sum2 - res_sum0);
			}

			delete[] global_psi;
			delete[] global_res;
			
		}
	}
	fflush(stdout);


	{//Broadcast any variables
		
		MPI_Barrier(mpi_comm);

		struct A {
			int a;
			int b;
			int c;
			double d1;
			double d2;
			std::complex<double> z;
		};

		A val;
		if (proc_id == 0) {
			val.a = 32;
			val.b = 128;
			val.c = 256;
			val.d1 = -3456.0;
			val.d2 = -0.001;
			val.z = { 3.0, -4.0 };
		}

		
		BroadcastAny(mpi_comm, 0, val.a, val.b, val.c, val.d1, val.d2, val.z);

		printf("[%d]BcastAny: %d, %d, %d, %f, %f, %f, %f\n", proc_id, val.a, val.b, val.c, val.d1, val.d2, val.z.real(), val.z.imag());

	}


	MPI_Finalize();
}


