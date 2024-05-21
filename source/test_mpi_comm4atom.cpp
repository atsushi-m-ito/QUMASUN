/*******************************

QUantum MAterial Simulation UNraveler (QUMASUN)

QUMASUN is numerical simulation code for Density Functional Theory(DFT) and Time-dependent DFT based on the real space grid.

********************************/
#include <mpi.h>

#include <cstdlib>
#include <cstdio>
#include <vector>
#include <random>
#include "mpi_helper.h"
#include "nucleus.h"
#include "GridRange.h"
#include "GetArg.h"

#define TEST_WITH_CLASS
#ifdef TEST_WITH_CLASS
#include "comm4atom.h"
#endif
//#define USE_GROUP

//test to cut sub firld and to interporate from field
int main(int argc, char* argv[]) {
	MPI_Init(&argc, &argv);
	MPI_Comm mpi_comm = MPI_COMM_WORLD;
	const int proc_id = GetProcessID(mpi_comm);
	const int num_procs = GetNumProcess(mpi_comm);

	const int num_test_nuclei = 32;
	const double cutoff_r = 2.0;
	const double cutoff_r2 = cutoff_r* cutoff_r;
	const double dx = 0.2;
	const double dy = 0.2;
	const double dz = 0.2;
	const int size_x = 128;
	const int size_y = 128;
	const int size_z = 128;
	const double box_x = (double)size_x * dx;
	const double box_y = (double)size_y * dy;
	const double box_z = (double)size_z * dz;
	const double box_half_x = box_x / 2.0;
	const double box_half_y = box_y / 2.0;
	const double box_half_z = box_z / 2.0;

	size_t global_size_3d = (size_t)size_x * (size_t)size_y * (size_t)size_z;
	GridRange global_grid{ 0,0,0,size_x,size_y,size_z };

	int split_num[3] = { 1,1,1 };
	int res = GetArgumentNumList<int>(argc, argv, "-ddm", 3, split_num);

	if (IsRoot(mpi_comm)) {
		printf("split grid by %d x %d x %d\n", split_num[0], split_num[1], split_num[2]);
		fflush(stdout);
	}
	MPI_Barrier(mpi_comm);
	GridRangeMPI l_grid = MakeRange(size_x, size_y, size_z, mpi_comm, split_num);
	printf("[%d] grid_size_3d = %zd\n", proc_id, l_grid.Size3D());
	fflush(stdout);

	//for communicator///////////////
	uint32_t bit_width = sizeof(uint8_t);
	int req_bit_size = num_procs;
	int req_uint_size = (req_bit_size + bit_width-1) / bit_width;
	uint32_t my_proc_bit = 0x1 >> (proc_id % bit_width);
	int my_proc_offset = proc_id / bit_width;

	//set positions
	
	
	std::vector<Nucleus> nuclei;
	nuclei.push_back(Nucleus{ 1, 1.0, 1.0, 1.0 });
	
	std::mt19937 mt(123456789);
	std::uniform_real_distribution<double> dist(.0, 1.0);
	for (int i = 0; i < num_test_nuclei; ++i) {
		nuclei.push_back(Nucleus{ 1, dist(mt) * box_x, dist(mt)* box_y, dist(mt)* box_z});
		
	}

#ifdef TEST_WITH_CLASS
	std::vector<CommForAtoms::AtomInfo> atoms;
	for (const auto& a : nuclei) {
		atoms.push_back({ a.Rx,a.Ry,a.Rz, cutoff_r });
	}
	int num_nuclei = nuclei.size();

	CommForAtoms comm4atoms;
	comm4atoms.CreateComms(&atoms[0], num_nuclei, l_grid, box_x, box_y, box_z, dx, dy, dz);


	printf("TEST MPI_Bcast=============================\n"); fflush(stdout);
	MPI_Barrier(mpi_comm);

	//test transfer//

	for (int i = 0; i < num_nuclei; ++i) {
		MPI_Comm my_comm = comm4atoms.GetComm(i);
		if (my_comm != MPI_COMM_NULL) {
			int proc_id_4_nucl;
			MPI_Comm_rank(my_comm, &proc_id_4_nucl);
			int value = proc_id_4_nucl;
			MPI_Bcast(&value, 1, MPI_INT, 0, my_comm);

			printf("nucleus[%d] -- proc[%d] as [%d]\n", i, proc_id, proc_id_4_nucl); fflush(stdout);

		}
		MPI_Barrier(mpi_comm);
	}

#else
	const double domain_begin_x = (double)l_grid.begin_x * dx;
	const double domain_end_x = (double)l_grid.end_x * dx;
	const double domain_begin_y = (double)l_grid.begin_y * dy;
	const double domain_end_y = (double)l_grid.end_y * dy;
	const double domain_begin_z = (double)l_grid.begin_z * dz;
	const double domain_end_z = (double)l_grid.end_z * dz;
	const double domain_center_x = (domain_end_x + domain_begin_x) / 2.0;
	const double domain_center_y = (domain_end_y + domain_begin_y) / 2.0;
	const double domain_center_z = (domain_end_z + domain_begin_z) / 2.0;
	const double domain_half_w_x = (domain_end_x - domain_begin_x) / 2.0;
	const double domain_half_w_y = (domain_end_y - domain_begin_y) / 2.0;
	const double domain_half_w_z = (domain_end_z - domain_begin_z) / 2.0;

	
	//check overlap atoms for each region/////////////


	const int64_t num_nuclei = nuclei.size();
	std::vector<uint8_t> interact_flags(num_nuclei * req_uint_size);
	std::vector<MPI_Comm> comm_4_nuclei(num_nuclei);

	for (int64_t i = 0; i < num_nuclei; ++i) {
		double x = fabs(nuclei[i].Rx - domain_center_x);
		double y = fabs(nuclei[i].Ry - domain_center_y);
		double z = fabs(nuclei[i].Rz - domain_center_z);
		if (x > box_half_x) x = box_x - x;
		if (y > box_half_y) y = box_y - y;
		if (z > box_half_z) z = box_z - z;

		bool is_interact = false;

		auto length2 = [](double x, double y) {
			return x * x + y * y;
			};
		auto length3 = [](double x, double y, double z) {
			return x * x + y * y + z * z;
			};

		if (x < domain_half_w_x + cutoff_r) {
			if (y < domain_half_w_y + cutoff_r) {
				if (z < domain_half_w_z + cutoff_r) {
					//以降は2方向でボックス範囲内なら確実に相互作用する
					if (x < domain_half_w_x) {
						if (y < domain_half_w_y) {
							//粒子は相互作用圏内//
							is_interact = true;
						} else if (z < domain_half_w_z) {
							//粒子は相互作用圏内//
							is_interact = true;
						} else {
							//y,zで範囲外なのでボックスの辺との距離を測る//
							double rr = length2(y - domain_half_w_y, z - domain_half_w_z);
							if (rr < cutoff_r2) {
								//粒子は相互作用圏内//
								is_interact = true;
							}
						}
					} else if (y < domain_half_w_y) {
						if (z < domain_half_w_z) {
							//粒子は相互作用圏内//
							is_interact = true;
						} else {
							//x,zで範囲外なのでボックスの辺との距離を測る//
							double rr = length2(x - domain_half_w_x, z - domain_half_w_z);
							if (rr < cutoff_r2) {
								//粒子は相互作用圏内//
								is_interact = true;
							}
						}
					} else if (z < domain_half_w_z) {
						//x,yで範囲外なのでボックスの辺との距離を測る//
						double rr = length2(x - domain_half_w_x, y - domain_half_w_y);
						if (rr < cutoff_r2) {
							//粒子は相互作用圏内//
							is_interact = true;
						}
					} else {
						//x,y,zで範囲外なのでボックスの辺との距離を測る//
						double rr = length3(x - domain_half_w_x, y - domain_half_w_y, z - domain_half_w_z);
						if (rr < cutoff_r2) {
							//粒子は相互作用圏内//
							is_interact = true;
						}
					}
				}
			}
		}

		if (is_interact) {
#if 1
			printf("[%d] nucleus[%d](%f,%f,%f) interacts with domain [%f,%f)x[%f,%f)x[%f,%f)\n",
				proc_id, i, nuclei[i].Rx, nuclei[i].Ry, nuclei[i].Rz,
				domain_begin_x, domain_end_x, domain_begin_y, domain_end_y, domain_begin_z, domain_end_z);
#endif
		}
#ifdef USE_GROUP
		if (is_interact) {
			interact_flags[i * req_uint_size + my_proc_offset] = my_proc_bit;
		}
#else
		if (is_interact) {
			MPI_Comm_split(mpi_comm, i, proc_id, &comm_4_nuclei[i]);
		} else {
			MPI_Comm_split(mpi_comm, MPI_UNDEFINED, proc_id, &comm_4_nuclei[i]);
		}
#endif
		
	}

#ifdef USE_GROUP
	//communicator creation//
	std::vector<uint8_t> interact_flags_red(num_nuclei * req_uint_size);
	MPI_Allreduce(&interact_flags[0], &interact_flags_red[0], num_nuclei * req_uint_size, MPI_BYTE, MPI_BOR, mpi_comm);

	for (int i = 0; i < num_nuclei; ++i) {
		std::vector<int> proc_list;
		int num = 0;
		for (int k = 0; k < req_uint_size; ++k) {
			for (int b = 0; b < bit_width; ++b) {
				if (interact_flags[i * req_uint_size + my_proc_offset] & (0x1 >> b)) {
					proc_list.push_back(b + bit_width * k);
					num++;
				}
			}
		}
		MPI_Group group_owner;
		MPI_Comm_group(mpi_comm, &group_owner);
		MPI_Group group;
		MPI_Group_incl(group_owner, num, &proc_list[0], &group);

		MPI_Comm_create(mpi_comm, group, &comm_4_nuclei[i]);
		MPI_Group_free(&group);

		

	}

#endif

	printf("TEST MPI_Bcast=============================\n"); fflush(stdout);
	MPI_Barrier(mpi_comm);

	//test transfer//

	for (int i = 0; i < num_nuclei; ++i) {
		if (comm_4_nuclei[i] != MPI_COMM_NULL) {
			int proc_id_4_nucl;
			MPI_Comm_rank(comm_4_nuclei[i], &proc_id_4_nucl );
			int value = proc_id_4_nucl;
			MPI_Bcast(&value, 1, MPI_INT, 0, comm_4_nuclei[i]);

			printf("nucleus[%d] -- proc[%d] as [%d]\n", i, proc_id, proc_id_4_nucl); fflush(stdout);

		}
		MPI_Barrier(mpi_comm);
	}

#endif

	MPI_Finalize();
}


