#pragma once
#ifdef USE_MPI
#include "mpi_helper.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <random>
#include "qumasun_kpoint.h"
#include "vecmath.h"
#include "vps_loader.h"
#include "pao_loader.h"
#include "cube_reader2.h"
#include "GramSchmidt_mpi.h"
#include "GridFor.h"

//#define READ_INITIAL_CUBE    "H_psi_ob1_30.cube"
//#define READ_INITIAL_PAO     1


inline
void QUMASUN_KPOINT::mInitializeStateRandom() {


	{
		//random state
		const int proc_id = GetProcessID(m_ddm_comm);
		const bool is_root_ddm = (proc_id == 0);
		const double dVol = m_dx * m_dy * m_dz;

		const double phase = 2.0 * M_PI * 0.0;// 0.0625;  //random phase rotation//
		const double cos_p = cos(phase);
		const double sin_p = sin(phase);
		const int local_size = ml_grid.Size3D();

		if (is_root_ddm) {

			double* psi_r = new double[m_size_3d];
			double* psi_i = new double[m_size_3d];

			const uint32_t seed = 123456789;
			std::mt19937 mt(seed);
			std::uniform_real_distribution<> distribution(0.0, 1.0);
			
			for (int sk = 0; sk < m_num_having_spin_kpoint ; ++sk) {
				double val = 1.0;

				for (int n = 0; n < num_solution; ++n) {

					for (size_t i = 0; i < m_size_3d; ++i) {
						psi_r[i] = distribution(mt);
					}
					vecmath::Normalize(psi_r, m_size_3d, dVol);
					for (size_t i = 0; i < m_size_3d; ++i) {
						size_t ix = i % m_size_x;
						double angle = 2.0 * M_PI * (double)ml_wave_set[sk].kpoint_x * m_dkx * ((double)ix * m_dx);

						psi_i[i] = sin( - angle)* psi_r[i];
						psi_r[i] *= cos( - angle);

					}
					
					//SoAC::SetZero(ml_wave_set[sk].l_psi_set[n], ml_grid.Size3D());
					mScatterField(ml_wave_set[sk].l_psi_set[n].re, psi_r);
					mScatterField(ml_wave_set[sk].l_psi_set[n].im, psi_i);
					/*
					for (size_t i = 0; i < local_size; ++i) {
						size_t ix = i % m_size_x;
						double angle = 2.0 * M_PI * (double)ml_wave_set[sk].kpoint_x * m_dkx * ((double)ix*m_dx);
						ml_wave_set[sk].l_psi_set[n].im[i] = (-angle) * ml_wave_set[sk].l_psi_set[n].re[i];
						ml_wave_set[sk].l_psi_set[n].re[i] *= cos_p;
					}
					*/

					val += distribution(mt);
					ml_wave_set[sk].eigen_values[n] = -1.0 / val;
					
				}

				MPI_Bcast(ml_wave_set[sk].eigen_values, num_solution, MPI_DOUBLE, 0, ml_grid.mpi_comm);
				
#ifdef INIT_GRAM_SCHMIDT
				GramSchmidt_ddm(ml_grid, num_solution, ml_wave_set[sk].l_psi_set, dVol);
#endif
			}

			delete[] psi_r;
			delete[] psi_i;

		} else {
			//slave process//

			for (int sk = 0; sk < m_num_having_spin_kpoint; ++sk) {
				for (int n = 0; n < num_solution; ++n) {
					//SoAC::SetZero(ml_wave_set[sk].l_psi_set[n], ml_grid.Size3D());
					mScatterField(ml_wave_set[sk].l_psi_set[n].re, nullptr);
					mScatterField(ml_wave_set[sk].l_psi_set[n].im, nullptr);
					/*
					for (size_t i = 0; i < local_size; ++i) {
						ml_wave_set[sk].l_psi_set[n].im[i] = sin_p * ml_wave_set[sk].l_psi_set[n].re[i];
						ml_wave_set[sk].l_psi_set[n].re[i] *= cos_p;
					}*/

				}

				MPI_Bcast(ml_wave_set[sk].eigen_values, num_solution, MPI_DOUBLE, 0, ml_grid.mpi_comm);
				
#ifdef INIT_GRAM_SCHMIDT
				GramSchmidt_ddm(ml_grid, num_solution, ml_wave_set[sk].l_psi_set, dVol);
#endif
			}
		}

		
	}

}

inline
void QUMASUN_KPOINT::mInitializeStateFile() {
	//using namespace QUMASUN;
	//const size_t num_file = m_initial_state.size();

	//FILE* fp = nullptr;
	const bool is_root_global = IsRoot(m_mpi_comm);
	const bool is_root_each_ddm = IsRoot(m_ddm_comm);



	
	
	const int TAG = 30000;
	const int TAGO = 310000;

	for (int sk = 0; sk < m_all_kinds_spin_kpoint; ++sk) {
		if ((m_having_spin_kpoint_begin <= sk) && (sk < m_having_spin_kpoint_begin + m_num_having_spin_kpoint)) {
			const auto& ws = ml_wave_set[sk - m_having_spin_kpoint_begin];

			for (int n = 0; n < num_solution; ++n) {
                double* state_vector_re = nullptr;
                double* state_vector_im = nullptr;
				if (is_root_each_ddm) {
					CubeReader2::Frame frame;
                    state_vector_re = CubeReader2::LoadCube(m_initial_state[sk * num_solution*2 + n*2].c_str(), &frame);
                    {
                        const char word[] = "eigenvalue=";
                        const int length = strlen(word);
                        //if (std::strncmp(frame.comment, word, length) == 0) {
                        const char* p_ev = std::strstr(frame.comment, word);
                        if (p_ev) {
                            ws.eigen_values[n] = strtod(p_ev+ length, nullptr);
                        } else {
                            ws.eigen_values[n] = 0.0;
                        }
                    }

                    long kpoint_x, kpoint_y, kpoint_z, kx_size, ky_size, kz_size;
                    if (m_initial_expand_from_kpoint) {
                        const char kword[] = "kpoint=";
                        const char* pk = std::strstr(frame.comment, kword);
                        if (pk) {
                            pk += strlen(kword);
                            char* pnext = nullptr;
                            kpoint_x = strtol(pk, &pnext, 10);
                            pk = pnext + 1;
                            kx_size = strtol(pk, &pnext, 10);
                            pk = pnext + 1;
                            kpoint_y = strtol(pk, &pnext, 10);
                            pk = pnext + 1;
                            ky_size = strtol(pk, &pnext, 10);
                            pk = pnext + 1;
                            kpoint_z = strtol(pk, &pnext, 10);
                            pk = pnext + 1;
                            kz_size = strtol(pk, &pnext, 10);
                            pk = pnext + 1;
                        }
                    }

                    state_vector_im = CubeReader2::LoadCube(m_initial_state[sk * num_solution*2 + n*2+1].c_str(), &frame);

                    //kpoint sampleした状態を、拡張してガンマ点に読み込み
                    if (m_initial_expand_from_kpoint) {
                        if ((frame.grid_x * kx_size != m_size_x) ||
                            (frame.grid_y * ky_size != m_size_y) ||
                            (frame.grid_z * kz_size != m_size_z)) {
                            printf("ERROR: at the expansion of initial state from kpoint, grid size is inconsistent\n"); fflush(stdout);
                            delete[] state_vector_re;
                            delete[] state_vector_im;
                            return;
                        }

                        double* expanded_re = new double[m_size_x * m_size_y * m_size_z];
                        double* expanded_im = new double[m_size_x * m_size_y * m_size_z];
                        const double gx_dx = (double)kpoint_x * m_dkx * m_dx;
                        const double gy_dy = (double)kpoint_y * m_dky * m_dy;
                        const double gz_dz = (double)kpoint_z * m_dkz * m_dz;

                        const int src_size_x = frame.grid_x;
                        const int src_size_y = frame.grid_y;
                        const int src_size_z = frame.grid_z;
                        
                        const double scale = 1.0 / std::sqrt((double)((m_size_x / src_size_x) * (m_size_y / src_size_y) * (m_size_z / src_size_z)));


                        for (int iz = 0; iz < m_size_z; ++iz) {
                            const double kzz = gz_dz * (double)iz;
                            for (int iy = 0; iy < m_size_y; ++iy) {
                                const double kyy = gy_dy * (double)iy;
                                for (int ix = 0; ix < m_size_x; ++ix) {
                                    const double kxx = gx_dx * (double)ix;
                                    const int i = ix + m_size_x * (iy + m_size_y * iz);
                                    const int si = (ix% src_size_x) + src_size_x * ((iy % src_size_y) + src_size_y * (iz % src_size_z));

                                    double cosikx = cos(kxx + kyy + kzz);
                                    double sinikx = sin(kxx + kyy + kzz);

                                    expanded_re[i] = scale*(cosikx * state_vector_re[si] - sinikx * state_vector_im[si]);
                                    expanded_im[i] = scale*(sinikx * state_vector_re[si] + cosikx * state_vector_im[si]);

                                }
                            }
                        }
                        delete[] state_vector_re;
                        delete[] state_vector_im;
                        state_vector_re = expanded_re;
                        state_vector_im = expanded_im;
                    }
				}

				mScatterField(ws.l_psi_set[n].re, state_vector_re);
				mScatterField(ws.l_psi_set[n].im, state_vector_im);
				

				if (is_root_each_ddm) {
					delete[] state_vector_re;
					delete[] state_vector_im;
				}
			}
		}

	}

}

inline
void QUMASUN_KPOINT::mInitializeState() {
	if (m_initial_state.empty()) {
		mInitializeStateRandom();
	} else {
		mInitializeStateFile();
	}
}

inline
void QUMASUN_KPOINT::mInitializeDensity() {
	if (m_initial_density == "none") {

		if (IsRoot(m_mpi_comm)) {
			vecmath::SetZero<double>(m_rho, m_size_3d);
		}
		mSetDensity(false);
	} else if (m_initial_density == "atom") {
		//rootで計算////////////////////////////////
		const bool is_root_global = IsRoot(m_mpi_comm);
		if (is_root_global) {
			RspaceFunc<double>& rho = m_rho;

			vecmath::SetZero<double>(rho, m_size_3d);

			std::map<int, Orbital_PAO> pao_list;
			for (const auto& p : m_atomic_wave_set) {
				pao_list.emplace(p.first, LoadPAO(p.second.c_str()));
			}


			for (int n = 0; n < m_num_nuclei; ++n) {
				const auto& pao = pao_list[m_nuclei[n].Z];
				const double rr_min = pow2(exp(pao.xi_min));

				const double R0_x = m_nuclei[n].Rx;
				const double R0_y = m_nuclei[n].Ry;
				const double R0_z = m_nuclei[n].Rz;

				//const int pao_id = (READ_INITIAL_PAO - 1 > 0) ? (READ_INITIAL_PAO - 1) : 0;
				const int pao_id = n;
				ForXYZ(GridRange{ 0,0,0,m_size_x,m_size_y,m_size_z },
					[&](int64_t i, int64_t ix, int64_t iy, int64_t iz) {
						const double rz2 = pow2(fold(m_dz * (double)iz - R0_z, m_box_z));
						const double ry2 = pow2(fold(m_dy * (double)iy - R0_y, m_box_y));
						const double rx2 = pow2(fold(m_dx * (double)ix - R0_x, m_box_x));

						const double rr = rx2 + ry2 + rz2;
						if (rr_min > rr) {
							double val = pao.valence_charge[0];
							rho[i] += val;
						} else {
							double val = RadialGrid2::GetValueBySquare(rr, pao.valence_charge, pao.num_radial_grids, pao.xi_min, pao.xi_delta);
							rho[i] += val;
						}

					});
			}


			double sum = 0.0;
			for (size_t i = 0; i < m_size_3d; ++i) {
				sum += rho[i];
			}
			sum *= m_dx * m_dy * m_dz;
			printf("init_rho=%f\n", sum);
            const double scale = (double)num_electrons / sum;
            for (size_t i = 0; i < m_size_3d; ++i) {
                rho[i] *= scale;
            }
            printf("init_rho(corrected)=%f\n", sum);

			for (auto&& p : pao_list) {
				p.second.Release();
			}


			if (is_spin_on) {
				const double ratio_diff = 1.0 / sum;
				for (size_t i = 0; i < m_size_3d; ++i) {
					m_rho_diff[i] = rho[i] * ratio_diff;
				}
			}
		}

		mHierarchyScatterField(ml_rho, m_rho.Pointer());
		if (is_spin_on) {
			mHierarchyScatterField(ml_rho_diff, m_rho_diff.Pointer());
		}
	} else {
        //file read//////////////////////////////////////////

		const bool is_root_spin = IsRoot(m_mpi_comm);
        if (is_root_spin) {

            CubeReader2::Frame frame;
            {
                double* density = CubeReader2::LoadCube(m_initial_density.c_str(), &frame);
                if (m_initial_expand_from_kpoint) {

                    const int src_size_x = frame.grid_x;
                    const int src_size_y = frame.grid_y;
                    const int src_size_z = frame.grid_z;

                    const double scale = 1.0 / (double)((m_size_x / src_size_x) * (m_size_y / src_size_y) * (m_size_z / src_size_z));

                    for (int iz = 0; iz < m_size_z; ++iz) {
                        for (int iy = 0; iy < m_size_y; ++iy) {
                            for (int ix = 0; ix < m_size_x; ++ix) {
                                const int i = ix + m_size_x * (iy + m_size_y * iz);
                                const int si = (ix % src_size_x) + src_size_x * ((iy % src_size_y) + src_size_y * (iz % src_size_z));
                                m_rho[i] = density[si];
                            }
                        }
                    }
                } else {
                    for (size_t i = 0; i < m_size_3d; ++i) {
                        m_rho[i] = density[i];
                    }
                }

                delete[] density;
            }

			if (is_spin_on) {
                if (m_initial_density_difference != "none" ) {
                    double* density_diff = CubeReader2::LoadCube(m_initial_density_difference.c_str(), &frame);
                    if (m_initial_expand_from_kpoint) {

                        const int src_size_x = frame.grid_x;
                        const int src_size_y = frame.grid_y;
                        const int src_size_z = frame.grid_z;

                        const double scale = 1.0 / (double)((m_size_x / src_size_x) * (m_size_y / src_size_y) * (m_size_z / src_size_z));

                        for (int iz = 0; iz < m_size_z; ++iz) {
                            for (int iy = 0; iy < m_size_y; ++iy) {
                                for (int ix = 0; ix < m_size_x; ++ix) {
                                    const int i = ix + m_size_x * (iy + m_size_y * iz);
                                    const int si = (ix % src_size_x) + src_size_x * ((iy % src_size_y) + src_size_y * (iz % src_size_z));
                                    m_rho_diff[i] = density_diff[si];
                                }
                            }
                        }
                    } else {
                        for (size_t i = 0; i < m_size_3d; ++i) {
                            m_rho_diff[i] = density_diff[i];
                        }
                    }


                    delete[] density_diff;
                } else {
                    for (size_t i = 0; i < m_size_3d; ++i) {
                        m_rho_diff[i] = 0.0;
                    }
                }
			}
		}

		mHierarchyScatterField(ml_rho, m_rho.Pointer());
		if (is_spin_on) {
			mHierarchyScatterField(ml_rho_diff, m_rho_diff.Pointer());
		}
	}


	if (IsRoot(m_mpi_comm)) {
		vecmath::Copy<double>(m_rho_prev, m_rho, m_size_3d);
	}
}


#endif
