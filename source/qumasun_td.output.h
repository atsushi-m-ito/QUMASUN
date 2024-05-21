#pragma once
#include "qumasun_td.h"
#include <cstdio>
#include <string>
#include "cube_writer.h"

inline
void QUMASUN_TD::OutputDensity(QUMASUN::OUTPUT_TARGET target, const char* filepath) {
	using namespace QUMASUN;

	if (!is_construction_successful)return;

	if (IsRoot(m_mpi_comm)) {
        watch.Restart();

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
        watch.Record(44);
	}
}

inline
void QUMASUN_TD::OutputEigenValue(const char* filepath) {
	using namespace QUMASUN;
	if (!is_construction_successful)return;

	FILE* fp = nullptr;
	const bool is_root_global = IsRoot(m_mpi_comm);
	const bool is_root_each_ddm = IsRoot(m_ddm_comm);
	
	if (is_root_each_ddm) {
		double* eigen_buf = nullptr;
		double* occupancy_buf = nullptr;
		if (is_root_global) {
			fp = fopen(filepath, "w");
			eigen_buf = new double[num_solution*2];
			occupancy_buf = eigen_buf + num_solution;
		}
		const int num_same_place_procs = GetNumProcess(m_same_ddm_place_comm);
		

		const int TAG = 40000;
		const int TAGO = 50000;
	
        /*
        std::vector<int> num_list(num_same_place_procs);
        MPI_Gather(&num_solution, 1, MPI_INT, &num_list[0], num_same_place_procs, 0, m_same_ddm_place_comm);

        std::vector<int> heads(num_same_place_procs + 1);
        heads[0] = 0;
        for (int p = 0; p < m_num_procs_sd; ++p) {
            heads[p+1] = heads[p] + num_list[p];
        }
        int total = heads[m_num_procs_sd];
        */

        std::vector<Kpoint3D> kpoint_list;
        int all_kinds_kpoint = ListupKpoints(m_kpoint_sampling[0], m_kpoint_sampling[1], m_kpoint_sampling[2], kpoint_symmetry, kpoint_list);



		for (int sk = 0; sk < m_all_kinds_spin_kpoint; ++sk) {
            if (is_root_global) {
                if (all_kinds_kpoint > sk) {
                    fprintf(fp, "Eigen values in kpoint and spin: %d/%d, %d/%d, %d/%d, %s\n", kpoint_list[sk].kx, m_kpoint_sampling[0], kpoint_list[sk].ky, m_kpoint_sampling[1], kpoint_list[sk].kz, m_kpoint_sampling[2], "up");
                } else {
                    fprintf(fp, "Eigen values in kpoint and spin: %d/%d, %d/%d, %d/%d, %s\n", kpoint_list[sk - all_kinds_kpoint].kx, m_kpoint_sampling[0], kpoint_list[sk - all_kinds_kpoint].ky, m_kpoint_sampling[1], kpoint_list[sk - all_kinds_kpoint].kz, m_kpoint_sampling[2], "down");
                }
            }

            std::vector<int> num_list(num_same_place_procs);
            bool has_kpoint = false;
            if ((m_having_spin_kpoint_begin <= sk) && (sk < m_having_spin_kpoint_begin + m_num_having_spin_kpoint)) {
                has_kpoint = true;
            }

            int num = has_kpoint ? num_solution : 0;

            MPI_Gather(&num, 1, MPI_INT, &num_list[0], num_same_place_procs, MPI_INT,0, m_same_ddm_place_comm);

            std::vector<int> heads(num_same_place_procs + 1);
            heads[0] = 0;
            for (int p = 0; p < m_num_procs_sd; ++p) {
                heads[p + 1] = heads[p] + num_list[p];
            }
            int total = heads[m_num_procs_sd];

            std::vector<double> occ_all(m_total_state);
            std::vector<double> eigen_all(m_total_state);
            double* ep = (has_kpoint) ? ml_wave_set[sk - m_having_spin_kpoint_begin].eigen_values : nullptr;
            double* eo = (has_kpoint) ? ml_wave_set[sk - m_having_spin_kpoint_begin].occupancy : nullptr;
                
            MPI_Gatherv(ep, num_solution, MPI_DOUBLE, &eigen_all[0], &num_list[0], &heads[0], MPI_DOUBLE, 0, m_same_ddm_place_comm);
            MPI_Gatherv(eo, num_solution, MPI_DOUBLE, &occ_all[0], &num_list[0], &heads[0], MPI_DOUBLE, 0, m_same_ddm_place_comm);
        
            if (is_root_global) {
                for (int i = 0; i < m_total_state; ++i) {
                    fprintf(fp, "%d\t%f\t%f\n", i, eigen_all[i], occ_all[i]);
                }
            }
		}
	
		if (is_root_global) {
			fclose(fp);
			delete[] eigen_buf;
		}
	}
}



inline
void QUMASUN_TD::OutputEigenVector(const char* filepath) {
    using namespace QUMASUN;
    if (!is_construction_successful)return;

    watch.Restart();
    
    const bool is_root_global = IsRoot(m_mpi_comm);
    const bool is_root_each_ddm = IsRoot(m_ddm_comm);


    double* state_vector_re = nullptr;
    double* state_vector_im = nullptr;

    if (is_root_each_ddm) {
        state_vector_re = new double[m_size_3d * 2];
        state_vector_im = state_vector_re + m_size_3d;
    }
    const int num_same_place_procs = GetNumProcess(m_same_ddm_place_comm);

    const int TAG = 40000;
    const int TAGO = 50000;

    const int all_kpoints = is_spin_on ? m_all_kinds_spin_kpoint / 2 : m_all_kinds_spin_kpoint;
    auto FileName = [&filepath, &all_kpoints](int sk, int n, const char* tail) {
        std::string filename(filepath);
        filename += (sk < all_kpoints) ? "_up" : "_down";
        filename += "_k" + std::to_string(sk % all_kpoints) + "_n" + std::to_string(n) + tail + ".cube";
        return filename;
        };

    if (is_root_global) {
        FILE* fp = fopen("state_file_list.txt", "w");
        for (int sk = 0; sk < m_all_kinds_spin_kpoint; ++sk) {
            for (int n = 0; n < m_total_state; ++n) {
                fprintf(fp, "%s\n", FileName(sk, n, "_re").c_str());
                fprintf(fp, "%s\n", FileName(sk, n, "_im").c_str());
            }
        }
        fclose(fp);
    }

    watch.Record(42);

    for (int sk = 0; sk < m_all_kinds_spin_kpoint; ++sk) {


        if ((m_having_spin_kpoint_begin <= sk) && (sk < m_having_spin_kpoint_begin + m_num_having_spin_kpoint)) {
            const auto& ws = ml_wave_set[sk - m_having_spin_kpoint_begin];



            for (int n = 0; n < num_solution; ++n) {
                mGatherField(state_vector_re, ws.l_psi_set[n].re);
                mGatherField(state_vector_im, ws.l_psi_set[n].im);
                watch.Record(43);

                if (is_root_each_ddm) {

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
                    sprintf(frame.comment, "eigenvalue=%.10f, kpoint=%d/%d,%d/%d,%d/%d", ws.eigen_values[n], ws.kpoint_x, m_kpoint_sampling[0], ws.kpoint_y, m_kpoint_sampling[1], ws.kpoint_z, m_kpoint_sampling[2]);

                    CubeWriter writer;
                    auto filepath_re = FileName(sk, n + m_begin_state, "_re");
                    auto filepath_im = FileName(sk, n + m_begin_state, "_im");
                    printf("save %s\n", filepath_re.c_str()); fflush(stdout);
                    writer.SaveCube(filepath_re.c_str(), frame, m_num_nuclei, m_nuclei, state_vector_re);
                    printf("save %s\n", filepath_im.c_str()); fflush(stdout);
                    writer.SaveCube(filepath_im.c_str(), frame, m_num_nuclei, m_nuclei, state_vector_im);


                    watch.Record(44);
                }
            }

        }

    }


    if (is_root_each_ddm) {
        delete[] state_vector_re;
    }


}

