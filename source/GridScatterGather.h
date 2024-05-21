#pragma once
#ifdef USE_MPI

#include <mpi.h>
#include "mpi_helper.h"
#include "GridRange.h"
#include "GridSubgrid.h"
#include <vector>

template <class T>
void ScatterGrid(const GridRangeMPI& grid, T* each_dest, const GridRange& global_grid, const T* src, int root_id, int HR_ratio_x=1, int HR_ratio_y=1, int HR_ratio_z=1, T* work = nullptr) {
	int proc_id;
	MPI_Comm_rank(grid.mpi_comm, &proc_id);
	int num_procs;
	MPI_Comm_size(grid.mpi_comm, &num_procs);

	const int local_size = grid.Size3D();
    const int HR_ratio = HR_ratio_x * HR_ratio_y * HR_ratio_z;

	//data preparation//
	if (grid.split_dimension == 1) {
		if (local_size * num_procs == global_grid.Size3D()) {
			//全領域が同じグリッドサイズ//			
			//printf("[%d]MPI_Scatter, size=%d, pnt=0x%x\n", proc_id, local_size, each_dest); fflush(stdout);
			//MPI_Barrier(grid.mpi_comm);
			MPI_Scatter(src, local_size* HR_ratio, GetDatatype<T>(), each_dest, local_size* HR_ratio, GetDatatype<T>(), root_id, grid.mpi_comm);

			//printf("[%d]MPI_Scatter: fin.\n", proc_id); fflush(stdout);
			//MPI_Barrier(grid.mpi_comm);

		} else {
			//全領域が異なるグリッドサイズ//
			int* send_sizes = nullptr;
			int* send_heads = nullptr;
			if (proc_id == root_id) {
				send_sizes = new int[num_procs * 2];
				send_heads = send_sizes + num_procs;
				const int Nz = global_grid.SizeZ();
				const int Nxy = global_grid.SizeX() * global_grid.SizeY();
				for (int n = 0; n < num_procs; ++n) {
					send_heads[n] = (Nxy * ((Nz * n) / num_procs))* HR_ratio;
					send_sizes[n] = (Nxy * ((Nz * (n + 1)) / num_procs))* HR_ratio - send_heads[n];
				//	printf("send_head, size=%d, %d\n", send_heads[n], send_sizes[n]); fflush(stdout);
				}
			}
			//printf("[%d](3)MPI_Scatter, size=%d, pnt=0x%x\n", proc_id, local_size, each_dest); fflush(stdout);
			MPI_Scatterv(src, send_sizes, send_heads, GetDatatype<T>(), each_dest, local_size* HR_ratio, GetDatatype<T>(), root_id, grid.mpi_comm);

			delete[] send_sizes;
		}
	} else {



		int* send_sizes = nullptr;
		int* send_heads = nullptr;

		T* sub_buf = nullptr;
		bool is_allocated = false;
		if (root_id == proc_id) {
			if (work) {
				sub_buf = work;
			} else {
				sub_buf = new T[global_grid.Size3D()* HR_ratio];
				//memset(sub_buf, 0, sizeof(T) * (global_grid.Size3D()-1));
				//printf("zero: %f, %f\n", sub_buf[0], sub_buf[global_grid.Size3D() -1]); fflush(stdout);
				is_allocated = true;
			}
			

			send_sizes = new int[num_procs * 2];
			send_heads = send_sizes + num_procs;

			const int Nx = global_grid.SizeX();
			const int Ny = global_grid.SizeY();
			const int Nz = global_grid.SizeZ();

            auto hr_global_grid = ScaleRange(global_grid, HR_ratio_x, HR_ratio_y, HR_ratio_z);

			int head = 0;
			int pid = 0;
			GridRange grid_p;
			for (int pz = 0; pz < grid.num_split_z; ++pz) {
				grid_p.begin_z = ((Nz * pz) / grid.num_split_z) * HR_ratio_z;
				grid_p.end_z = ((Nz * (pz + 1)) / grid.num_split_z)* HR_ratio_z;

				for (int py = 0; py < grid.num_split_y; ++py) {
					grid_p.begin_y = ((Ny * py) / grid.num_split_y) * HR_ratio_y;
					grid_p.end_y = ((Ny * (py + 1)) / grid.num_split_y)* HR_ratio_y;

					for (int px = 0; px < grid.num_split_x; ++px) {						
						grid_p.begin_x = ((Nx * px) / grid.num_split_x) * HR_ratio_x;
						grid_p.end_x = ((Nx * (px+1)) / grid.num_split_x)* HR_ratio_x;
						
#if 0
						printf("head = sub_buf + %zd\n", head);
#endif
						CutSubgrid(grid_p, sub_buf + head, hr_global_grid, src);
						send_heads[pid] = head;
						send_sizes[pid] = grid_p.Size3D();
						head += send_sizes[pid];
						++pid;
					}
				}
			}

		}

		//printf("test33: %d\n", num_procs); fflush(stdout);

		MPI_Scatterv(sub_buf, send_sizes, send_heads, GetDatatype<T>(), each_dest, local_size*HR_ratio, GetDatatype<T>(), root_id, grid.mpi_comm);


		if (is_allocated == true) {
			delete[] sub_buf;
		}
		delete[]send_sizes;
	}

}


template <class T>
void GatherGrid(const GridRange& global_grid, T* dest, const GridRangeMPI& grid, const T* each_src, int root_id, int HR_ratio_x=1, int HR_ratio_y=1, int HR_ratio_z=1, T* work = nullptr) {
	

	int proc_id;
	MPI_Comm_rank(grid.mpi_comm, &proc_id);
	int num_procs;
	MPI_Comm_size(grid.mpi_comm, &num_procs);

	const int local_size = grid.Size3D();
    const int HR_ratio = HR_ratio_x * HR_ratio_y * HR_ratio_z;

	//data preparation//
	if (grid.split_dimension == 1) {
		if (local_size * num_procs == global_grid.Size3D()) {			
			MPI_Gather(each_src, local_size*HR_ratio, GetDatatype<T>(), dest, local_size * HR_ratio, GetDatatype<T>(), root_id, grid.mpi_comm);
		} else {
			//全領域が異なるグリッドサイズ//
			int* send_sizes = nullptr;
			int* send_heads = nullptr;
			if (proc_id == root_id) {
				send_sizes = new int[num_procs * 2];
				send_heads = send_sizes + num_procs;
				const int Nz = global_grid.SizeZ();
				const int Nxy = global_grid.SizeX() * global_grid.SizeY();
				for (int n = 0; n < num_procs; ++n) {
					send_heads[n] = (Nxy * ((Nz * n) / num_procs))* HR_ratio;
					send_sizes[n] = (Nxy * ((Nz * (n + 1)) / num_procs))* HR_ratio - send_heads[n];
				}
			}
			MPI_Gatherv(each_src, local_size* HR_ratio, GetDatatype<T>(), dest, send_sizes, send_heads, GetDatatype<T>(), root_id, grid.mpi_comm);

			delete[] send_sizes;
		}
	} else {
		int* send_sizes = nullptr;
		int* send_heads = nullptr;

		T* gather_buf = nullptr;
		bool is_allocated = false;
		if (root_id == proc_id) {
			if (work) {
				gather_buf = work;
			} else {
				gather_buf = new T[global_grid.Size3D()* HR_ratio];
				is_allocated = true;
			}

			send_sizes = new int[num_procs * 2];
			send_heads = send_sizes + num_procs;

			const int Nx = global_grid.SizeX();
			const int Ny = global_grid.SizeY();
			const int Nz = global_grid.SizeZ();

			int head = 0;
			int pid = 0;
			GridRange grid_p;
			for (int pz = 0; pz < grid.num_split_z; ++pz) {
				grid_p.begin_z = ((Nz * pz) / grid.num_split_z)*HR_ratio_z;
				grid_p.end_z = ((Nz * (pz + 1)) / grid.num_split_z)* HR_ratio_z;

				for (int py = 0; py < grid.num_split_y; ++py) {
					grid_p.begin_y = ((Ny * py) / grid.num_split_y) * HR_ratio_y;
					grid_p.end_y = ((Ny * (py + 1)) / grid.num_split_y)* HR_ratio_y;

					for (int px = 0; px < grid.num_split_x; ++px) {
						grid_p.begin_x = ((Nx * px) / grid.num_split_x) * HR_ratio_x;
						grid_p.end_x = ((Nx * (px + 1)) / grid.num_split_x)* HR_ratio_x;

#if 0
						printf("head = sub_buf + %zd\n", head);
#endif
						//CutSubgrid(grid_p, sub_buf + head, global_grid, src);
						send_heads[pid] = head;
						send_sizes[pid] = grid_p.Size3D();
						head += send_sizes[pid];
						++pid;
					}
				}
			}
		}

		MPI_Gatherv(each_src, local_size*HR_ratio, GetDatatype<T>(), gather_buf, send_sizes, send_heads, GetDatatype<T>(), root_id, grid.mpi_comm);


		//printf("each_src = 0x%zx, 0x%zx\n", each_src, gather_buf); fflush(stdout);

		if (IsRoot(grid.mpi_comm)) {



			const int Nx = global_grid.SizeX();
			const int Ny = global_grid.SizeY();
			const int Nz = global_grid.SizeZ();

            auto hr_global_grid = ScaleRange(global_grid, HR_ratio_x, HR_ratio_y, HR_ratio_z);

			int head = 0;
			GridRange grid_p;
			for (int pz = 0; pz < grid.num_split_z; ++pz) {
				grid_p.begin_z = ((Nz * pz) / grid.num_split_z)*HR_ratio_z;
				grid_p.end_z = ((Nz * (pz + 1)) / grid.num_split_z)* HR_ratio_z;

				for (int py = 0; py < grid.num_split_y; ++py) {
					grid_p.begin_y = ((Ny * py) / grid.num_split_y) * HR_ratio_y;
					grid_p.end_y = ((Ny * (py + 1)) / grid.num_split_y)* HR_ratio_y;

					for (int px = 0; px < grid.num_split_x; ++px) {
						grid_p.begin_x = ((Nx * px) / grid.num_split_x) * HR_ratio_x;
						grid_p.end_x = ((Nx * (px + 1)) / grid.num_split_x)* HR_ratio_x;

#if 0 
						printf("head = sub_buf + %zd\n", head);
#endif
						//CutSubgrid(grid_p, sub_buf + head, global_grid, src);
						PasteSubgrid(hr_global_grid, dest, grid_p,	gather_buf + head);
						head += grid_p.Size3D();						
					}
				}
			}



			if (is_allocated == true) {
				delete[] gather_buf;
			}
		}
		
		delete[] send_sizes;
	}

}

#endif
