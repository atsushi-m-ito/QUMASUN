#pragma once
#include "qumasun_mpi.h"
#include <cstdio>
#include "cube_writer.h"

inline
void QUMASUN_MPI::OutputDensity(QUMASUN::OUTPUT_TARGET target, const char* filepath) {
	using namespace QUMASUN;
	if (!is_construction_successful)return;

	if (IsRoot(m_mpi_comm)) {
		CubeWriter::Frame frame;
		frame.grid_x = m_size_x;
		frame.grid_y = m_size_y;
		frame.grid_z = m_size_z;
		frame.boxaxis[0] = m_box_x;
		frame.boxaxis[1] = 0.0;
		frame.boxaxis[2] = 0.0;
		frame.boxaxis[3] = 0.0;
		frame.boxaxis[4] = m_box_y;
		frame.boxaxis[5] = 0.0;
		frame.boxaxis[6] = 0.0;
		frame.boxaxis[7] = 0.0;
		frame.boxaxis[8] = m_box_z;
		frame.boxorg[0] = 0.0;
		frame.boxorg[1] = 0.0;
		frame.boxorg[2] = 0.0;

		CubeWriter writer;


		switch (target) {
		case OUTPUT_TARGET::Density:
			writer.SaveCube(filepath, frame, m_num_nuclei, m_nuclei, m_rho);
			break;
		case OUTPUT_TARGET::Vhart:
			writer.SaveCube(filepath, frame, m_num_nuclei, m_nuclei, m_Vhart);
			break;
		}
	}
}



inline
void QUMASUN_MPI::OutputEigenValue(const char* filepath) {
	using namespace QUMASUN;
	if (!is_construction_successful)return;

	FILE* fp = nullptr;
	const bool is_root_global = IsRoot(m_mpi_comm);
	const bool is_root_each_ddm = IsRoot(m_ddm_comm);
	const int num_procs_global = GetNumProcess(m_mpi_comm);
	const int proc_id = GetProcessID(m_mpi_comm);

	const int TAG = 40000;
	const int TAGO = 50000;

	if (is_root_global) {
		fp = fopen(filepath, "w");
	

		fprintf(fp, "Eigen values in spin: up\n");
		for (int i = 0; i < num_solution; ++i) {
			fprintf(fp, "%d\t%f\t%f\n", i, m_eigen_values[i], m_occupancy[i]);
		}

		if (is_spin_on) {
			fprintf(fp, "Eigen values in spin: down\n");
			/*
			* すでにoccupationの計算時にrootプロセスに集約されている//
			if (m_num_spin == 1) {
				//spin_parallel//
				const int down_proc_id = num_procs_global / 2;
				MPI_Status status;
				MPI_Recv(m_eigen_values_down, num_solution, MPI_DOUBLE, down_proc_id, TAG, m_mpi_comm, &status);
				MPI_Recv(m_occupancy_down, num_solution, MPI_DOUBLE, down_proc_id, TAGO, m_mpi_comm, &status);
			}*/
			for (int i = 0; i < num_solution; ++i) {
				fprintf(fp, "%d\t%f\t%f\n", i, m_eigen_values_down[i], m_occupancy_down[i]);
			}
		}
		fclose(fp);
	} /*else {
		//すでにoccupationの計算時にrootプロセスに送っている//
		const int down_proc_id = num_procs_global / 2;
		if (is_spin_on && (proc_id == down_proc_id)) {
			MPI_Send(m_eigen_values, num_solution, MPI_DOUBLE, 0, TAG, m_mpi_comm);
			MPI_Send(m_occupancy, num_solution, MPI_DOUBLE, 0, TAGO, m_mpi_comm);

		}
		
	}*/

	
}

