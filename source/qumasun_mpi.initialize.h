#pragma once
#ifdef USE_MPI
#include "mpi_helper.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <random>
#include "qumasun_mpi.h"
#include "vecmath.h"
#include "vps_loader.h"
#include "pao_loader.h"
#include "cube_reader2.h"
#include "GramSchmidt_mpi.h"
#include "GridFor.h"

//#define READ_INITIAL_CUBE    "H_psi_ob1_30.cube"
//#define READ_INITIAL_PAO     1

//#define INIT_GRAM_SCHMIDT


void QUMASUN_MPI::mInitializeState() {


#ifdef READ_INITIAL_CUBE
	CubeReader2::Frame frame;
	double* density = CubeReader2::LoadCube(m_density_filepath, &frame);

	RspaceFunc<double>& psi_r = m_psi_set[0];
	double sum = 0.0;
	for (size_t i = 0; i < m_size_3d; ++i) {
		psi_r[i] = sqrt(density[i]);
		sum += density[i];
	}
	sum *= m_dx * m_dy * m_dz;
	for (size_t i = 0; i < m_size_3d; ++i) {
		density[i] /= sum;
		psi_r[i] = sqrt(density[i]);
		//sum += density[i];
	}
	printf("init_rho=%f\n", sum);
	/*
	{

		//double* density = CubeReader2::LoadCube(m_density_filepath, &frame);
		double* psi0 = CubeReader2::LoadCube("PP/H_psi_ob1_30.cube", &frame);
		double* vhart = CubeReader2::LoadCube("PP/H_vhart.cube", &frame);

		RspaceFunc<double>& psi_r = m_psi_set[0];
		double sum = 0.0;
		double sum2 = 0.0;
		double Ehart = 0.0;
		double Ehart2 = 0.0;
		double Ehart3 = 0.0;
		double sum_vhart = 0.0;
		for (size_t i = 0; i < m_size_3d; ++i) {
			psi_r[i] = psi0[i];
			density[i]= psi0[i]* psi0[i];
			sum += density[i];
			sum2 += psi_r[i] * psi_r[i];
			Ehart += (psi_r[i]* psi_r[i])* vhart[i];
			m_Vhart[i] = vhart[i];
			sum_vhart += vhart[i];
		}

		RspaceFunc<double>& den2 = m_work_rd;
		mSetLaplasianFFT(den2, m_Vhart, 1.0 / (4.0 * M_PI));
		double sum3 = 0.0;
		double sum_vhart2 = 0.0;
		for (size_t i = 0; i < m_size_3d; ++i) {
			sum3 += den2[i];
		}
		sum3 = (1.0/ (m_dx * m_dy * m_dz) - sum3) / (double)m_size_3d;
		{
			double dif_rho_max = -DBL_MAX;
			double dif_rho_min = -DBL_MIN;
			for (size_t i = 0; i < m_size_3d; ++i) {
				den2[i] = den2[i] + sum3;
				Ehart2 += den2[i] * vhart[i];

				if (dif_rho_max < den2[i] - density[i]) dif_rho_max = den2[i] - density[i];
				if (dif_rho_min > den2[i] - density[i]) dif_rho_min = den2[i] - density[i];
			}
			printf("dif_rho %f, %f\n", dif_rho_min, dif_rho_max);
		}
		RspaceFunc<double> den(density);
		mSetPotentialVhart(m_Vhart, den);
		for (size_t i = 0; i < m_size_3d; ++i) {
			sum_vhart2 += m_Vhart[i];
			Ehart3 += (psi_r[i] * psi_r[i]) * m_Vhart[i];
		}

		RspaceFunc<double> den3(m_size_3d);
		mSetLaplasianFFT(den3, m_Vhart, 1.0 / (4.0 * M_PI));
		sum3 = 0.0;
		for (size_t i = 0; i < m_size_3d; ++i) {
			sum3 += den3[i];
		}
		sum3 = (1.0 / (m_dx * m_dy * m_dz) - sum3) / (double)m_size_3d;
		double Ehart4 = 0.0;
		{
			double dif_rho_max = -DBL_MAX;
			double dif_rho_min = -DBL_MIN;
			for (size_t i = 0; i < m_size_3d; ++i) {
				den3[i] = den3[i] + sum3;
				Ehart4 += den3[i] * vhart[i];

				if (dif_rho_max < den3[i] - density[i]) dif_rho_max = den3[i] - density[i];
				if (dif_rho_min > den3[i] - density[i]) dif_rho_min = den3[i] - density[i];
			}
			printf("dif_rho2 %f, %f\n", dif_rho_min, dif_rho_max);
		}

		sum *= m_dx * m_dy * m_dz;
		sum2 *= m_dx * m_dy * m_dz;
		Ehart *= m_dx * m_dy * m_dz / 2.0;
		Ehart2 *= m_dx * m_dy * m_dz / 2.0;
		Ehart3 *= m_dx * m_dy * m_dz / 2.0;
		Ehart4 *= m_dx * m_dy * m_dz / 2.0;
		sum_vhart *= m_dx * m_dy * m_dz;
		sum_vhart2 *= m_dx * m_dy * m_dz;
		printf("init_rho2=%f, %f\n", sum, sum2);
		printf("Ehart=%f, %f, %f, %f\n", Ehart, Ehart2, Ehart3, Ehart4);
		printf("sum_vhart=%.15f, %.15f\n", sum_vhart, sum_vhart2);
		delete[] psi0;
		delete[] vhart;
	}
	*/

	delete[] density;


	//dummy values are set in eigen values//
	for (int n = 0; n < num_solution; ++n) {
		m_eigen_values[n] = -1.0 / (double)((n + 1) * (n + 1));
	}

#elif defined( READ_INITIAL_PAO)
	if (mode_initial_state == 10) {
		auto pao = LoadPAO("PP/Be7.0.pao");

		const double rr_min = pow2(exp(pao.xi_min));
		const double R0_x = m_nuclei[0].Rx;
		const double R0_y = m_nuclei[0].Ry;
		const double R0_z = m_nuclei[0].Rz;
		for (int n = 0; n < num_solution; ++n) {
			//const int pao_id = (READ_INITIAL_PAO - 1 > 0) ? (READ_INITIAL_PAO - 1) : 0;
			const int pao_id = n;
			RspaceFunc<double>& psi_r = m_psi_set[n];
			for (int iz = 0; iz < m_size_z; ++iz) {
				const double rz2 = pow2(fold(m_dz * (double)iz - R0_z, m_box_z));
				for (int iy = 0; iy < m_size_y; ++iy) {
					const double ry2 = pow2(fold(m_dy * (double)iy - R0_y, m_box_y));
					for (int ix = 0; ix < m_size_x; ++ix) {
						const double rx2 = pow2(fold(m_dx * (double)ix - R0_x, m_box_x));
						const size_t i = (size_t)ix + (size_t)m_size_x * ((size_t)iy + ((size_t)m_size_y * (size_t)iz));

						const double rr = rx2 + ry2 + rz2;
						if (rr_min > rr) {
							double val = pao.orbitals[pao_id][0] / sqrt(4.0 * M_PI);
							psi_r[i] = val;
						} else {
							double val = RadialGrid2::GetValueBySquare(rr, pao.orbitals[pao_id], pao.num_radial_grids, pao.xi_min, pao.xi_delta) / sqrt(4.0 * M_PI);
							psi_r[i] = val;
						}
					}
				}
			}


			//dummy values are set in eigen values//
			m_eigen_values[n] = -1.0 / (double)((n + 1) * (n + 1));

			double sum = 0.0;
			for (size_t i = 0; i < m_size_3d; ++i) {
				sum += psi_r[i] * psi_r[i];
			}
			sum *= m_dx * m_dy * m_dz;
			printf("init_rho=%f\n", sum);
		}

		const double dVol = m_dx * m_dy * m_dz;
		GramSchmidt(m_size_3d, num_solution, m_psi_set, dVol);


	} else


#elif 1
	{
		//random state
		const int proc_id = GetProcessID(m_ddm_comm);
		const bool is_root_ddm = (proc_id == 0);
		const double dVol = m_dx * m_dy * m_dz;

		if (is_root_ddm) {
			double* psi_r = new double[m_size_3d];

			const uint32_t seed = 123456789;
			std::mt19937 mt(seed);
			std::uniform_real_distribution<> distribution(0.0, 1.0);
			{
				double val = 1.0;


				for (int n = 0; n < num_solution; ++n) {
					for (size_t i = 0; i < m_size_3d; ++i) {
						psi_r[i] = distribution(mt);
					}
					vecmath::Normalize(psi_r, m_size_3d, dVol);
					mScatterField(ml_psi[n], psi_r);

					val += distribution(mt);
					m_eigen_values[n] = -1.0 / val;					
				}
				MPI_Bcast(m_eigen_values, num_solution, MPI_DOUBLE, 0, ml_grid.mpi_comm);
			}

			if (m_num_spin == 2) {
				double val = 1.0;

				for (int n = 0; n < num_solution; ++n) {


					for (size_t i = 0; i < m_size_3d; ++i) {
						psi_r[i] = distribution(mt);
					}
					vecmath::Normalize(psi_r, m_size_3d, dVol);
					mScatterField(ml_psi[n + num_solution], psi_r);//for down spin//

					val += distribution(mt);
					m_eigen_values_down[n] = -1.0 / val;					
				}
				MPI_Bcast(m_eigen_values_down, num_solution, MPI_DOUBLE, 0, ml_grid.mpi_comm);
			}
			delete[] psi_r;

		} else {
			//slave process//
			for (int n = 0; n < num_solution; ++n) {
				mScatterField(ml_psi[n], nullptr);				
			}
			MPI_Bcast(m_eigen_values, num_solution, MPI_DOUBLE, 0, ml_grid.mpi_comm);

			if (m_num_spin == 2) {
				for (int n = 0; n < num_solution; ++n) {
					mScatterField(ml_psi[n + num_solution], nullptr);//for down spin//					
				}
				MPI_Bcast(m_eigen_values_down, num_solution, MPI_DOUBLE, 0, ml_grid.mpi_comm);
			}
		}
#ifdef INIT_GRAM_SCHMIDT 
		GramSchmidt_ddm(ml_grid, num_solution, ml_psi, dVol);
		if (m_num_spin == 2) {
			GramSchmidt_ddm(ml_grid, num_solution, ml_psi + num_solution, dVol);
		}
#endif
	}

#elif 1
	
	if (is_root) {
		double* psi_r = new double[m_size_3d];
		for (size_t i = 0; i < m_size_3d; ++i) {
			psi_r[i] = 0.0;
		}
		double sum = 0.0;
		for (int ni = 0; ni < m_num_nuclei; ++ni) {
			//const double Ze = (double)m_nuclei[ni].Z;
			const double R0_x = m_nuclei[ni].Rx;
			const double R0_y = m_nuclei[ni].Ry;
			const double R0_z = m_nuclei[ni].Rz;

			for (int iz = 0; iz < m_size_z; ++iz) {
				const double rz2 = SQ(Length(m_dz * (double)iz - R0_z, m_box_z));
				for (int iy = 0; iy < m_size_y; ++iy) {
					const double ry2 = SQ(Length(m_dy * (double)iy - R0_y, m_box_y));
					for (int ix = 0; ix < m_size_x; ++ix) {
						const double rx2 = SQ(Length(m_dx * (double)ix - R0_x, m_box_x));
						const size_t i = ix + m_size_x * (iy + (m_size_y * iz));

						const double val = exp(-sqrt(rx2 + ry2 + rz2));
						psi_r[i] += val;
						sum += val * val;
					}
				}
			}
		}
		{
			sum = 1.0 / sqrt(sum * m_dx * m_dy * m_dz);
			//sum = 1.0 / sqrt(sum );
			double sum2 = 0.0;
			for (size_t i = 0; i < m_size_3d; ++i) {
				psi_r[i] *= sum;
				sum2 += psi_r[i] * psi_r[i];
			}
			sum2 *= m_dx * m_dy * m_dz;
			printf("sum2=%f\n", sum2); fflush(stdout);
		}
		m_eigen_values[0] = 0.0;

		mScatterField(ml_psi[0], psi_r);

		delete[] psi_r;
	} else {
		mScatterField(ml_psi[0], nullptr);
	}

#else
	{
		RspaceFunc<double>& psi_r = m_psi_set[0];
		for (size_t i = 0; i < m_size_3d; ++i) {
			psi_r[i] = 1.0;
		}
		psi_r[0] =  2.0;
	}
#endif
}


void QUMASUN_MPI::mInitializeDensity() {
	const int proc_id = GetProcessID(m_mpi_comm);

	if (m_initial_density == "none") {

		if (IsRoot(m_mpi_comm)) {
			vecmath::SetZero<double>(m_rho_prev, m_size_3d);
		}
		mSetDensity(false);
	} else if (m_initial_density == "atom") {
		//rootで計算////////////////////////////////	
		const bool is_root_spin = IsRoot(m_mpi_comm);
		if (is_root_spin) {
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


		if (is_spin_on) {
			if (m_num_spin == 1) {
				const int num_procs = GetNumProcess(m_mpi_comm);
				const int TAG = 2300;
				if (proc_id == 0) {
					MPI_Send(m_rho.Pointer(), m_size_3d, MPI_DOUBLE, num_procs / 2, TAG, m_mpi_comm);					
				} else if (proc_id == num_procs / 2) {
					MPI_Status status;
					MPI_Recv(m_rho.Pointer(), m_size_3d, MPI_DOUBLE, 0, TAG, m_mpi_comm, &status);
				}
			}
		}
		mScatterField(ml_rho, m_rho.Pointer());

	} else{

		const bool is_root_spin = IsRoot(m_mpi_comm);
		if (is_root_spin) {

			CubeReader2::Frame frame;
			double* density = CubeReader2::LoadCube(m_initial_density.c_str(), &frame);
			for (size_t i = 0; i < m_size_3d; ++i) {
				m_rho[i] = density[i];
			}

			if (is_spin_on) {
				for (size_t i = 0; i < m_size_3d; ++i) {
					m_rho_diff[i] = 0.0;
				}
			}

			delete[] density;
		}


		if (is_spin_on) {
			if (m_num_spin == 1) {
				const int num_procs = GetNumProcess(m_mpi_comm);
				const int TAG = 2301;
				if (proc_id == 0) {
					MPI_Send(m_rho.Pointer(), m_size_3d, MPI_DOUBLE, num_procs / 2, TAG, m_mpi_comm);
				} else if (proc_id == num_procs / 2) {
					MPI_Status status;
					MPI_Recv(m_rho.Pointer(), m_size_3d, MPI_DOUBLE, 0, TAG, m_mpi_comm, &status);
				}
			}
		}
		mScatterField(ml_rho, m_rho.Pointer());
	}

	if (IsRoot(m_mpi_comm)) {
		vecmath::Copy<double>(m_rho_prev, m_rho, m_size_3d);
	}
}


#endif
