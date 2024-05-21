#pragma once
#include "qumasun_kpoint.h"
#include <cstdio>
#include <string>
#include "cube_writer.h"

inline
void QUMASUN_KPOINT::OutputDensity(QUMASUN::OUTPUT_TARGET target, const char* filepath) {
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
        case OUTPUT_TARGET::DensityHR:
            frame.grid_x *= m_HR_ratio_x;
            frame.grid_y *= m_HR_ratio_y;
            frame.grid_z *= m_HR_ratio_z;
            writer.SaveCube(filepath, frame, m_num_nuclei, m_nuclei, m_hr_rho);
            break;
        }
	}
}

inline
void QUMASUN_KPOINT::OutputEigenValue(const char* filepath) {
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
		const int id_of_num_same_place_procs = GetNumProcess(m_same_ddm_place_comm);

		const int TAG = 40000;
		const int TAGO = 50000;
	
		for (int sk = 0; sk < m_all_kinds_spin_kpoint; ++sk) {
			if (is_root_global) {
				if ((m_having_spin_kpoint_begin <= sk) && (sk < m_having_spin_kpoint_begin + m_num_having_spin_kpoint)) {
					const auto& ws = ml_wave_set[sk - m_having_spin_kpoint_begin];
					fprintf(fp, "Eigen values in kpoint and spin: %d/%d, %d/%d, %d/%d, %s\n", ws.kpoint_x, m_kpoint_sampling[0], ws.kpoint_y, m_kpoint_sampling[1], ws.kpoint_z, m_kpoint_sampling[2], (ws.spin == SPIN::UP) ? "up" : "down");
					for (int i = 0; i < num_solution; ++i) {
						fprintf(fp, "%d\t%f\t%f\n", i, ws.eigen_values[i], ws.occupancy[i]);
					}
				} else {
					MPI_Status status;
					int kpoint_spin[4];
					MPI_Recv(kpoint_spin, 4, MPI_INT, MPI_ANY_SOURCE, TAG + sk, m_same_ddm_place_comm, &status);
					MPI_Recv(eigen_buf, num_solution, MPI_DOUBLE, MPI_ANY_SOURCE, TAG + sk, m_same_ddm_place_comm, &status);
					MPI_Recv(occupancy_buf, num_solution, MPI_DOUBLE, MPI_ANY_SOURCE, TAGO + sk, m_same_ddm_place_comm, &status);

					int kx = kpoint_spin[0];
					int ky = kpoint_spin[1];
					int kz = kpoint_spin[2];
					int spin = kpoint_spin[3];
                    fprintf(fp, "Eigen values in kpoint and spin: %d/%d, %d/%d, %d/%d, %s\n", kx, m_kpoint_sampling[0], ky, m_kpoint_sampling[1], kz, m_kpoint_sampling[2], (spin == 0) ? "up" : "down");
					for (int i = 0; i < num_solution; ++i) {
						fprintf(fp, "%d\t%f\t%f\n", i, eigen_buf[i], occupancy_buf[i]);
					}
				}
			} else {//no-root//
				if ((m_having_spin_kpoint_begin <= sk) && (sk < m_having_spin_kpoint_begin + m_num_having_spin_kpoint)) {
					const auto& ws = ml_wave_set[sk - m_having_spin_kpoint_begin];
					int kpoint_spin[4]{ ws.kpoint_x, ws.kpoint_y, ws.kpoint_z, ws.spin == SPIN::UP ? 0 : 1 };
					MPI_Send(kpoint_spin, 4, MPI_INT, 0, TAG + sk, m_same_ddm_place_comm);
					MPI_Send(ws.eigen_values, num_solution, MPI_DOUBLE, 0, TAG + sk, m_same_ddm_place_comm);
					MPI_Send(ws.occupancy, num_solution, MPI_DOUBLE, 0, TAGO + sk, m_same_ddm_place_comm);
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
void QUMASUN_KPOINT::OutputEigenVector(const char* filepath) {
	using namespace QUMASUN;
	if (!is_construction_successful)return;

    
	const bool is_root_global = IsRoot(m_mpi_comm);
	const bool is_root_each_ddm = IsRoot(m_ddm_comm);
    //const bool proc_id = GetProcessID(m_mpi_comm);
    //DEBUG_PRINTF("[%d]OutputEigenVector1\n", proc_id);



	double* state_vector_re = nullptr;
	double* state_vector_im = nullptr;

	if (is_root_each_ddm) {
        state_vector_re = new double[m_size_3d * 2];
        state_vector_im = state_vector_re + m_size_3d;
	}
	const int num_same_place_procs = GetNumProcess(m_same_ddm_place_comm);
	const int id_of_num_same_place_procs = GetNumProcess(m_same_ddm_place_comm);

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
            for (int n = 0; n < num_solution; ++n) {
                fprintf(fp, "%s\n", FileName(sk, n, "_re").c_str());
                fprintf(fp, "%s\n", FileName(sk, n, "_im").c_str());
            }
        }
        fclose(fp);
    }

	for (int sk = 0; sk < m_all_kinds_spin_kpoint; ++sk) {

        
		if ((m_having_spin_kpoint_begin <= sk) && (sk < m_having_spin_kpoint_begin + m_num_having_spin_kpoint)) {
			const auto& ws = ml_wave_set[sk - m_having_spin_kpoint_begin];

			

			for (int n = 0; n < num_solution; ++n) {
				mGatherField(state_vector_re, ws.l_psi_set[n].re);
				mGatherField(state_vector_im, ws.l_psi_set[n].im);

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
                    writer.SaveCube(FileName(sk, n, "_re").c_str(), frame, m_num_nuclei, m_nuclei, state_vector_re);
                    writer.SaveCube(FileName(sk, n, "_im").c_str(), frame, m_num_nuclei, m_nuclei, state_vector_im);

                    /*
					double gkx = m_dkx * (double)ws.kpoint_x;
					fprintf(fp, "# %d th eigen vector\n", n);
					for (int i = 0; i < m_size_x; ++i) {
						double cos1 = cos(gkx * m_dx * (double)i);
						double sin1 = sin(gkx * m_dx * (double)i);
						double re_1 = cos1 * eigen_vector_re[i] - sin1 * state_vector_im[i];
						double im_1 = sin1 * eigen_vector_re[i] + cos1 * state_vector_im[i];
						double ang = atan2(eigen_vector_re[i], eigen_vector_im[i]) / (2.0 * M_PI);
						double ang_1 = atan2(re_1, im_1) / (2.0 * M_PI);
						fprintf(fp, "%d\t%f\t%f\t%f\t%f\t%f\t%f\n", i, eigen_vector_re[i], eigen_vector_im[i], ang, re_1,im_1, ang_1);
					}
                    */
				}
			}

		}

	}


    if (is_root_each_ddm) {
        delete[] state_vector_re;
    }

}

