#pragma once
#include <mpi.h>
#include <vector>
#include "nucleus.h"
#include "GridRange.h"

class CommForAtoms {
public:
	struct AtomInfo {
		double x;
		double y;
		double z;
		double cutoff;
	};
private:
	std::vector<MPI_Comm> comm_4_nuclei;

public:
	~CommForAtoms() {
		DeleteComms();
	}

	int CreateComms(const AtomInfo* nuclei, int num_atoms, const GridRangeMPI& l_grid,
		double box_x, double box_y, double box_z,
		double dx, double dy, double dz	) {
		const MPI_Comm& mpi_comm = l_grid.mpi_comm;
		int proc_id;
		MPI_Comm_rank(mpi_comm, &proc_id);

		const double box_half_x = box_x / 2.0;
		const double box_half_y = box_y / 2.0;
		const double box_half_z = box_z / 2.0;
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

		auto length2 = [](double x, double y) {
			return x * x + y * y;
			};
		auto length3 = [](double x, double y, double z) {
			return x * x + y * y + z * z;
			};


		comm_4_nuclei.resize(num_atoms);
		for (int64_t i = 0; i < num_atoms; ++i) {
			double x = fabs(nuclei[i].x - domain_center_x);
			double y = fabs(nuclei[i].y - domain_center_y);
			double z = fabs(nuclei[i].z - domain_center_z);
			if (x > box_half_x) x = box_x - x;
			if (y > box_half_y) y = box_y - y;
			if (z > box_half_z) z = box_z - z;

			bool is_interact = false;

			const double cutoff_r = nuclei[i].cutoff;
			const double cutoff_r2 = cutoff_r* cutoff_r;

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
				MPI_Comm_split(mpi_comm, i, proc_id, &comm_4_nuclei[i]);
			} else {
				MPI_Comm_split(mpi_comm, MPI_UNDEFINED, proc_id, &comm_4_nuclei[i]);
			}


		}
		return 0;
	}

	MPI_Comm GetComm(int i) {
		if (i < comm_4_nuclei.size()) {
			return comm_4_nuclei[i];
		} else {
			return MPI_COMM_NULL;
		}
	}

	void DeleteComms() {
		for (auto&& a : comm_4_nuclei) {
			if (a != MPI_COMM_NULL) {
				MPI_Comm_free(&a);
			}
		}
	}

	int CountValidComms() {
		int num = 0;
		for (auto&& a : comm_4_nuclei) {
			if (a != MPI_COMM_NULL) {
				++num;
			}
		}
		return num;
	}
};
