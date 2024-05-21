#pragma once
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include "qumasun_td.h"
#include "poisson_fft.h"
#include "Vxc.h"
#include "nucleus.h"
#include "vecmath.h"

//#define NORMALIZE_CORE_CHARGE_TD

/*
* Perform the following calculations
* - Vext or Vlocal potential due to nuclei
* - Core-Core repulsive energy
* where they only need to be calculated once before the SCF calculation
*/
inline
void QUMASUN_TD::mPrepareCore() {
	const bool is_root = IsRoot(m_mpi_comm);
	

	if (IsRoot(m_same_ddm_place_comm)) {
		if (m_hamiltonian_type == HAMILTONIAN::KohnSham_PP) {//DFT//

            //Vext calculation on HR (high resolution) grid///////////////////////////////////////////////////////

            m_pp_integrator.SetChargeVlocal_HR(ml_hr_nucl_rho, ml_grid, m_nuclei, m_num_nuclei);
            watch.Record(31);

            if (IsRoot(ml_grid.mpi_comm)) {
                double* hr_nucl_rho = m_work;
                mGatherField_HR(hr_nucl_rho, ml_hr_nucl_rho);

                double sumrho = 0.0;
                const size_t hr_size_3d = m_size_3d * m_HR_ratio_x * m_HR_ratio_y * m_HR_ratio_z;
                for (size_t i = 0; i < hr_size_3d; ++i) {
                    sumrho += hr_nucl_rho[i];
                }
                sumrho *= m_dx * m_dy * m_dz / (double)(m_HR_ratio_x * m_HR_ratio_y * m_HR_ratio_z);
                printf("ChargeForVlocal(HR)=%f\n", sumrho);
                watch.Record(33);
#if 1
                SetPotentialByPoissonFFT_v2(m_HR_fftw, m_hr_Vext, hr_nucl_rho, m_size_x * m_HR_ratio_x, m_size_y * m_HR_ratio_y, m_size_z * m_HR_ratio_z, m_dx / (double)m_HR_ratio_x, m_dy / (double)m_HR_ratio_y, m_dz / (double)m_HR_ratio_z);
#else
                SetPotentialByPoissonFFT(m_hr_Vext, hr_nucl_rho, m_size_x * m_HR_ratio_x, m_size_y* m_HR_ratio_y, m_size_z* m_HR_ratio_z, m_dx/ (double)m_HR_ratio_x, m_dy / (double)m_HR_ratio_y, m_dz / (double)m_HR_ratio_z, m_work + hr_size_3d);
#endif
                watch.Record(34);
            } else {
                mGatherField_HR(nullptr, ml_hr_nucl_rho);
            }

            mScatterField_HR(ml_hr_Vext, m_hr_Vext);
            //printf("[%d]mGatherField\n", proc_id); fflush(stdout);
            watch.Record(35);


            //Vext calculation on low resolution grid///////////////////////////////////////////////////////

			m_pp_integrator.SetChargeVlocal_v2(ml_nucl_rho, ml_grid, m_nuclei, m_num_nuclei);
            watch.Record(21);
			//printf("[%d]SetChargeVlocal_v2\n", proc_id); fflush(stdout);

#ifdef VEXT_FROM_HR

            if (IsRoot(ml_grid.mpi_comm)) {
                DownConvert_Realspace(m_Vext.Pointer(), m_size_x, m_size_y, m_size_z, m_hr_Vext.Pointer(), m_HR_ratio_x, m_HR_ratio_y, m_HR_ratio_z);
                watch.Record(36);
            }
            mScatterField(ml_Vext, m_Vext);
            watch.Record(25);
#else

			if (IsRoot(ml_grid.mpi_comm)) {
				RspaceFunc<double> nucl_rho(m_work + m_size_3d * 5);
				mGatherField(nucl_rho, ml_nucl_rho);

				double sumrho = 0.0;
				for (size_t i = 0; i < m_size_3d; ++i) {
					sumrho += nucl_rho[i];
				}
				sumrho *= m_dx * m_dy * m_dz;
				printf("ChargeForVlocal=%f\n", sumrho);
                watch.Record(23);
				SetPotentialByPoissonFFT(m_Vext, nucl_rho, m_size_x, m_size_y, m_size_z, m_dx, m_dy, m_dz, m_work + m_size_3d);
                watch.Record(24);


                printf("Vlocal[center]=%f\n", m_Vext[m_size_x / 2 + m_size_x * (m_size_y / 2 + m_size_y * (m_size_z / 2))]);

			} else {
				mGatherField(nullptr, ml_nucl_rho);
			}

			mScatterField(ml_Vext, m_Vext);
			//printf("[%d]mGatherField\n", proc_id); fflush(stdout);
            watch.Record(25);
#endif
		}
	}

#if 1

    if (IsRoot(m_same_ddm_place_comm)) {
        if (m_hamiltonian_type == HAMILTONIAN::KohnSham_PP) {//DFT//

            m_pp_integrator.SetPccCharge_v2(ml_pcc_rho, ml_grid, m_nuclei, m_num_nuclei);
            watch.Record(26);
            mGatherField(m_pcc_rho, ml_pcc_rho);
            watch.Record(27);
        } else {
            if (IsRoot(ml_grid.mpi_comm)) {
                mSetPotentialVext(m_Vext, m_nuclei, m_num_nuclei);
            }
        }
    }

#else
	if (is_root) {

		//Vext or Vlocal/////////////////////////////////////////
		vecmath::SetZero(m_pcc_rho.Pointer(), m_size_3d);

		if (m_hamiltonian_type == HAMILTONIAN::KohnSham_PP) {//DFT//

			//Vxcの計算にのみPCC chargeを加味する//
			m_pp_integrator.SetPccCharge(m_pcc_rho, m_nuclei, m_num_nuclei);
			mScatterField(ml_pcc_rho, m_pcc_rho);
		} else {
			mSetPotentialVext(m_Vext, m_nuclei, m_num_nuclei);
		}
	} else {
		if(IsRoot(m_same_ddm_place_comm)){
			if (m_hamiltonian_type == HAMILTONIAN::KohnSham_PP) {//DFT//
				
				mScatterField(ml_pcc_rho, nullptr);
			}
		}
	}
#endif
	

	if (is_root) {
		//Valence Electron/////////////////////////////////////////
		{

			for (int i = 0; i < m_num_nuclei; ++i) {
				m_nuclei_valence_elecron[i] = m_pp_integrator.NumValenceElectron(m_nuclei[i].Z);
			}

            watch.Record(28);
		}

		//E_core_core and PP内側の補正(VlocalとZZ/rの差)//////////////////////////////////////////////
		{
			constexpr double cutoff_length_PP = 2.0;// PPのVlocalとZZ/rがズレ始める距離
			double EV = 0.0;
			double E_self = 0.0;

			auto Folding = [](double x, double box_width) {
				return x - floor(x / box_width) * box_width;
				};

			auto Distance = [](double x, double box_width) {
				return (fabs(x) < box_width * 0.5) ? x : -copysign(box_width - fabs(x), x);
				};


			for (int i = 0; i < m_num_nuclei; ++i) {
				const double Qi = m_nuclei_valence_elecron[i];
				double xi = Folding(m_nuclei[i].Rx, m_box_x);
				double yi = Folding(m_nuclei[i].Ry, m_box_y);
				double zi = Folding(m_nuclei[i].Rz, m_box_z);

				//自己相互作用：原子核自身の発するVlocalの原点でのエネルギー//
				const double Vlocal_0 = m_pp_integrator.GetVlocalRadial(m_nuclei[i].Z, 0.0);
				E_self += -(Qi * (-Vlocal_0));


				//近接の原子核との相互作用のCoulombとVlocalの差//
				//擬ポテンシャルのカットオフ長よりも内側におけるVlocalと理想のZZ/rとの差を補正//
				for (int k = 0; k < m_num_nuclei; ++k) {
					if (i == k) continue;
					const double Qk = m_nuclei_valence_elecron[k];
					double xk = Folding(m_nuclei[k].Rx, m_box_x);
					double yk = Folding(m_nuclei[k].Ry, m_box_y);
					double zk = Folding(m_nuclei[k].Rz, m_box_z);
					double dx = Distance(xk - xi, m_box_x);
					double dy = Distance(yk - yi, m_box_y);
					double dz = Distance(zk - zi, m_box_z);
					double r = sqrt(dx * dx + dy * dy + dz * dz);
					if (r < cutoff_length_PP) {
						const double Vlocal_i = m_pp_integrator.GetVlocalRadial(m_nuclei[i].Z, r);
						const double diff_Coulomb_Vlocal = (Qi * Qk / r) - (Qk * (-Vlocal_i));
						EV += diff_Coulomb_Vlocal;
					}
				}
			}

			m_Ecore_self = E_self / 2.0;
			m_diff_Coulomb_Vlocal = EV / 2.0;

            watch.Record(29);
		}


        //(v2)E_core_core and PP内側の補正(VlocalとZZ/rの差)//////////////////////////////////////////////
        {
            //カットオフの記憶//
            std::vector<double> cutoff_vlocal(m_num_nuclei);
            int prev_Z = 0;
            double prev_cut = 0.0;
            for (int i = 0; i < m_num_nuclei; ++i) {
                if (m_nuclei[i].Z == prev_Z) {
                    cutoff_vlocal[i] = prev_cut;
                } else {
                    prev_cut = m_pp_integrator.GetCutoffVlocal(m_nuclei[i].Z);
                    cutoff_vlocal[i] = prev_cut;
                    prev_Z = m_nuclei[i].Z;
                }
            }

            //force clear//
            m_force_nn_correction.resize(m_num_nuclei);
            m_force_nn_TF.resize(m_num_nuclei);
            for (int i = 0; i < m_num_nuclei; ++i) {
                m_force_nn_correction[i] = { 0.0,0.0,0.0 };
                m_force_nn_TF[i] = { 0.0,0.0,0.0 };
            }

            double EV1 = 0.0;
            double EV2 = 0.0;

            auto Folding = [](double x, double box_width) {
                return x - floor(x / box_width) * box_width;
                };

            auto Distance = [](double x, double box_width) {
                return (fabs(x) < box_width * 0.5) ? x : -copysign(box_width - fabs(x), x);
                };

            auto E_TF = [](double R, double Z_a, double Z_b, double Q_a, double Q_b, double a0, double b0) {
                return (Q_a + (Z_a - Q_a) * exp(-R / a0)) * (Q_b + (Z_b - Q_b) * exp(-R / b0)) / R;
                };

            auto Force_TF = [](double R, double Z_a, double Z_b, double Q_a, double Q_b, double a0, double b0) {
                double F = (Q_a + (Z_a - Q_a) * exp(-R / a0)) * (Q_b + (Z_b - Q_b) * exp(-R / b0)) / (R);
                F += ((Z_a - Q_a) * exp(-R / a0)) * (Q_b + (Z_b - Q_b) * exp(-R / b0)) / (a0);
                F += (Q_a + (Z_a - Q_a) * exp(-R / a0)) * ((Z_b - Q_b) * exp(-R / b0)) / (b0);
                return F / R;
                };


            for (int i = 0; i < m_num_nuclei; ++i) {
                const double Qi = m_nuclei_valence_elecron[i];
                const double Z_a = m_nuclei[i].Z;
                const double a0 = cutoff_vlocal[i] / (2.0 * Z_a / Qi);
                double xi = Folding(m_nuclei[i].Rx, m_box_x);
                double yi = Folding(m_nuclei[i].Ry, m_box_y);
                double zi = Folding(m_nuclei[i].Rz, m_box_z);



                //近接の原子核との相互作用のCoulombとVlocalの差//
                //擬ポテンシャルのカットオフ長よりも内側におけるVlocalと理想のZZ/rとの差を補正//
                for (int k = 0; k < m_num_nuclei; ++k) {
                    if (i == k) continue;
                    const double cutoff_length_PP = cutoff_vlocal[i] + cutoff_vlocal[k];// PPのVlocalとZZ/rがズレ始める距離

                    const double Qk = m_nuclei_valence_elecron[k];
                    const double Z_b = m_nuclei[k].Z;
                    const double b0 = cutoff_vlocal[k] / (2.0 * Z_b / Qk);
                    double xk = Folding(m_nuclei[k].Rx, m_box_x);
                    double yk = Folding(m_nuclei[k].Ry, m_box_y);
                    double zk = Folding(m_nuclei[k].Rz, m_box_z);
                    double dx = Distance(xk - xi, m_box_x);
                    double dy = Distance(yk - yi, m_box_y);
                    double dz = Distance(zk - zi, m_box_z);
                    double r = sqrt(dx * dx + dy * dy + dz * dz);
                    if (r < cutoff_length_PP) {
                        double force_nn = 0.0;
                        const double Enn_correct = m_pp_integrator.CoreCoreEnergyCorrection(r, m_nuclei[i].Z, m_nuclei[k].Z, &force_nn);
                        EV1 += Enn_correct;
                        EV2 += E_TF(r, Z_a, Z_b, Qi, Qk, a0, b0);
                        double force_tf = Force_TF(r, Z_a, Z_b, Qi, Qk, a0, b0);
                        {
                            const double fx = force_nn * dx / (r * 2.0);
                            const double fy = force_nn * dy / (r * 2.0);
                            const double fz = force_nn * dz / (r * 2.0);
                            m_force_nn_correction[i].x -= fx;
                            m_force_nn_correction[i].y -= fy;
                            m_force_nn_correction[i].z -= fz;
                            m_force_nn_correction[k].x += fx;
                            m_force_nn_correction[k].y += fy;
                            m_force_nn_correction[k].z += fz;
                        }
                        {
                            const double fx = force_tf * dx / (r * 2.0);
                            const double fy = force_tf * dy / (r * 2.0);
                            const double fz = force_tf * dz / (r * 2.0);
                            m_force_nn_TF[i].x -= fx;
                            m_force_nn_TF[i].y -= fy;
                            m_force_nn_TF[i].z -= fz;
                            m_force_nn_TF[k].x += fx;
                            m_force_nn_TF[k].y += fy;
                            m_force_nn_TF[k].z += fz;
                        }
                    }
                }
            }


            m_Enn_close_correction = (EV2 - EV1) / 2.0;

            watch.Record(29);
        }

        //自己相互作用補正//
        m_Ecore_self_v2 = m_pp_integrator.GetSelfCoreEnergy()/2.0;
        m_Ecore_self_HR = m_pp_integrator.GetSelfCoreEnergy_HR() / 2.0;
	}
}

