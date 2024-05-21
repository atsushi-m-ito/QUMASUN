#pragma once
#ifdef USE_MPI
#include <mpi.h>
#endif


struct GridI3D {
	int x;
	int y;
	int z;
};


struct GridRange {
public:
	int begin_x;
	int begin_y;
	int begin_z;
	int end_x;
	int end_y;
	int end_z;

public:
	size_t Size3D() const {
		return (size_t)(end_x - begin_x)
			* (size_t)(end_y - begin_y)
			* (size_t)(end_z - begin_z);
	}


	int SizeX() const {
		return (end_x - begin_x);
	}

	int SizeY() const {
		return (end_y - begin_y);
	}

	int SizeZ() const {
		return (end_z - begin_z);
	}
        
};

inline 
GridRange ScaleRange(const GridRange& src, int ratio_x, int ratio_y, int ratio_z) {
    GridRange res{ src };
    res.begin_x*= ratio_x;
    res.begin_y*= ratio_y;
    res.begin_z*= ratio_z;
    res.end_x*= ratio_x;
    res.end_y*= ratio_y;
    res.end_z*= ratio_z;
    return res;
}

inline
GridRange MakeRange(int global_x, int global_y, int global_z) {
	return GridRange{ 0,0,0, global_x, global_y, global_z };
}

#ifdef USE_MPI

struct GridRangeMPI : public GridRange {
public:
	MPI_Comm mpi_comm;
	int split_dimension;
	int num_split_x;
	int num_split_y;
	int num_split_z;
	int proc_x_left;
	int proc_x_right;
	int proc_y_left;
	int proc_y_right;
	int proc_z_left;
	int proc_z_right;


};

inline
GridRangeMPI ScaleRange(const GridRangeMPI& src, int ratio_x, int ratio_y, int ratio_z) {
    GridRangeMPI res{ src };
    res.begin_x *= ratio_x;
    res.begin_y *= ratio_y;
    res.begin_z *= ratio_z;
    res.end_x *= ratio_x;
    res.end_y *= ratio_y;
    res.end_z *= ratio_z;
    return res;
}

inline
GridRangeMPI MakeRange(int global_x, int global_y, int global_z, MPI_Comm mpi_comm, const int split_num[3]) {
	//constexpr int split_dimension = 1;
	
	int proc_id;
	int num_procs;
	MPI_Comm_rank(mpi_comm, &proc_id);
	MPI_Comm_size(mpi_comm, &num_procs);
	const int num_split_all = split_num[0] * split_num[1] * split_num[2];
	if ((num_split_all != num_procs)&&(num_split_all*2 != num_procs)) {
		printf("ERROR: split num is not difference from the num of total procs.\n");
		return { 0,0,0,0,0,0 };
	}
	if (split_num[0] * split_num[1] * split_num[2] == 0) {
		printf("ERROR: split num includes zero.\n");
		return { 0,0,0,0,0,0 };
	}

	const int proc_x = proc_id % split_num[0];
	const int proc_y = (proc_id / split_num[0]) % split_num[1];
	const int proc_z = proc_id / (split_num[0] * split_num[1]);

	//preliminary 1D split//
	
	GridRangeMPI grid;
	grid.begin_x = (global_x * proc_x) / split_num[0];
	grid.end_x = (global_x * (proc_x+1)) / split_num[0];
	grid.begin_y = (global_y * proc_y) / split_num[1];
	grid.end_y = (global_y * (proc_y+1)) / split_num[1];
	grid.begin_z = (global_z * proc_z) / split_num[2];
	grid.end_z = (global_z * (proc_z + 1)) / split_num[2];
	grid.mpi_comm = mpi_comm;

	auto MakeProcID = [&](int px, int py, int pz) {
		return ((px + split_num[0]) % split_num[0])
			+ split_num[0] * (((py + split_num[1]) % split_num[1])
				+ split_num[1] * ((pz + split_num[2]) % split_num[2]));

		};

	grid.proc_x_left = MakeProcID(proc_x - 1,proc_y,proc_z);
	grid.proc_x_right = MakeProcID(proc_x + 1, proc_y, proc_z);
	grid.proc_y_left = MakeProcID(proc_x, proc_y - 1, proc_z);
	grid.proc_y_right = MakeProcID(proc_x, proc_y + 1, proc_z);
	grid.proc_z_left = MakeProcID(proc_x, proc_y, proc_z - 1);
	grid.proc_z_right = MakeProcID(proc_x, proc_y, proc_z + 1);
	//grid.split_dimension = ((split_num[0] > 1) ? 1 : 0) + ((split_num[1] > 1) ? 1 : 0) + ((split_num[2] > 1) ? 1 : 0);
	//grid.split_dimension = 3;
	grid.split_dimension = (split_num[0] > 1) ? 3 : (split_num[1] > 1) ? 2 : (split_num[2] > 1) ? 1 : 0;
	grid.num_split_x = split_num[0];
	grid.num_split_y = split_num[1];
	grid.num_split_z = split_num[2];

	//printf("split_dimension = %d\n", grid.split_dimension);
	return grid;
	
}

#endif

