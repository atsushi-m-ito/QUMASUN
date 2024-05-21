#pragma once
//#define _USE_MATH_DEFINES
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include "qumasun_td.h"
#include "wrap_fft.h"
#include "Vxc.h"
#include "nucleus.h"
#include "vecmath.h"
#include "convertHR.h"

//#define USE_R2C_C2R


/*
* V_hart, Vx, Vc V_totを計算する
* 
* V_extはmPrepareCoreで最初に一度だけ計算
* is_calc_Vhart_xcがfalseの時は電子密度が変化していないのでVhartとVxcの計算をスキップ
*/
inline
void QUMASUN_TD::mSetPotential(bool is_calc_Vhart_xc) {
	const bool is_root_spin = IsRoot(m_mpi_comm);
	if (is_root_spin) {

		if (m_hamiltonian_type == HAMILTONIAN::Schrodinger) {
            /*
			if (is_calc_Vhart_xc) {
            	SetPotentialByPoissonFFT(m_Vhart, m_rho, m_size_x, m_size_y, m_size_z, m_dx, m_dy, m_dz, m_work);
			}
			for (size_t i = 0; i < m_size_3d; ++i) {
				m_Vtot[i] = m_Vext[i] + m_Vhart[i];
			}
            */
		} else if (m_hamiltonian_type == HAMILTONIAN::KohnSham_PP) {//DFT//
			if(is_calc_Vhart_xc){

#if 1
                watch.Record(4);
                {

                    const int64_t size_3d = m_size_3d;

                    auto* rho_comp = m_fftw->GetBuffer();
                    for (size_t i = 0; i < size_3d; ++i) {
                        rho_comp[i][0] = m_rho[i];
                        rho_comp[i][1] = 0.0;
                    }
                }
                auto rhok = m_fftw->ForwardExecute(nullptr, nullptr);
                watch.Record(50);
                
                //calculate HR rho
                {
                    //OneComplex* rhok = (OneComplex*)rhok;
                    OneComplex* hr_rhok = (OneComplex*)(m_HR_fftw->GetBuffer());
                    UpConvert_Kspace((OneComplex*)rhok, m_size_x, m_size_y, m_size_z, hr_rhok, m_HR_ratio_x, m_HR_ratio_y, m_HR_ratio_z);
                    watch.Record(51);
                    const size_t hr_size_3d = m_size_3d * m_HR_ratio_x * m_HR_ratio_y * m_HR_ratio_z;
                    //OneComplex* hr_rho_c = hr_rhok + hr_size_3d;
                    //IFFT_3D(hr_rho_c, hr_rhok, m_size_x * m_HR_ratio_x, m_size_y * m_HR_ratio_y, m_size_z * m_HR_ratio_z);
                    auto* res = m_HR_fftw->BackwardExecute(nullptr, nullptr);
                    OneComplex* hr_rho_c = (OneComplex*)res;
                    watch.Record(52);
                    double sum_rho = 0.0;
                    //const double scaling_factor = (double)(m_HR_ratio_x * m_HR_ratio_y * m_HR_ratio_z)/(double);
                    const double scaling_factor = 1.0/(double)m_size_3d;
                    for (size_t i = 0; i < hr_size_3d; ++i) {
                        m_hr_rho[i] = hr_rho_c[i].r * scaling_factor;
                        sum_rho += hr_rho_c[i].r * scaling_factor;
                    }
                    printf("Rho(HR) = %f\n", sum_rho * m_dx * m_dy * m_dz / (double)(m_HR_ratio_x * m_HR_ratio_y * m_HR_ratio_z));
                    watch.Record(53);

                    double diff = 0.0;
                    for (int iz = 0; iz < m_size_z; ++iz) {
                        for (int iy = 0; iy < m_size_y; ++iy) {
                            for (int ix = 0; ix < m_size_x; ++ix) {
                                int64_t i = ix + m_size_x * (iy + m_size_y * iz);
                                int64_t oi = m_HR_ratio_x * (ix + m_size_x * (m_HR_ratio_y * (iy + m_size_y * m_HR_ratio_z * iz)));
                                //int64_t oi = (ix + hr_ratio_x * size_x * (iy + hr_ratio_y * size_y * iz));
                                diff += (m_hr_rho[oi] - m_rho[i]) * (m_hr_rho[oi] - m_rho[i]);
                            }
                        }
                    }
                    printf("diff_Rho(HR) = %f\n", diff * m_dx * m_dy * m_dz);
                    watch.Record(54);
                }

                SolvePoissonInKspace(rhok, m_size_x, m_size_y, m_size_z, m_dx, m_dy, m_dz);
                watch.Record(50);

                //calculate HR Vhart
                {
                    OneComplex* vhartk = (OneComplex*)(rhok);
                    OneComplex* hr_vhartk = (OneComplex*)(m_HR_fftw->GetBuffer());
                    UpConvert_Kspace(vhartk, m_size_x, m_size_y, m_size_z, hr_vhartk, m_HR_ratio_x, m_HR_ratio_y, m_HR_ratio_z);
                    watch.Record(51);
                    const size_t hr_size_3d = m_size_3d * m_HR_ratio_x * m_HR_ratio_y * m_HR_ratio_z;
                    //OneComplex* hr_vhart_c = hr_vhartk + hr_size_3d;
                    //IFFT_3D(hr_vhart_c, hr_vhartk, m_size_x * m_HR_ratio_x, m_size_y * m_HR_ratio_y, m_size_z * m_HR_ratio_z);
                    auto* res = m_HR_fftw->BackwardExecute(nullptr, nullptr);
                    OneComplex* hr_vhart_c = (OneComplex*)res;
                    watch.Record(52);
                    //double sum_rho = 0.0;
                    //const double scaling_factor = (double)(m_HR_ratio_x * m_HR_ratio_y * m_HR_ratio_z);
                    const double scaling_factor = 1.0 / (double)m_size_3d;
                    for (size_t i = 0; i < hr_size_3d; ++i) {
                        m_hr_Vhart[i] = hr_vhart_c[i].r * scaling_factor;
                    }
                    watch.Record(53);

                }

                auto Vh = m_fftw->BackwardExecute(nullptr, nullptr);
                const double invN = 1.0 / (double)(m_size_3d);
                for (size_t i = 0; i < m_size_3d; ++i) {
                    m_Vhart[i] = Vh[i][0] * invN;
                }
                watch.Record(50);

                {
                    double diff = 0.0;
                    for (int iz = 0; iz < m_size_z; ++iz) {
                        for (int iy = 0; iy < m_size_y; ++iy) {
                            for (int ix = 0; ix < m_size_x; ++ix) {
                                int64_t i = ix + m_size_x * (iy + m_size_y * iz);
                                int64_t oi = m_HR_ratio_x * (ix + m_size_x * (m_HR_ratio_y * (iy + m_size_y * m_HR_ratio_z * iz)));
                                //int64_t oi = (ix + hr_ratio_x * size_x * (iy + hr_ratio_y * size_y * iz));
                                diff += (m_hr_Vhart[oi] - m_Vhart[i]) * (m_hr_Vhart[oi] - m_Vhart[i]);
                            }
                        }
                    }
                    printf("diff_Vhart(HR) = %f\n", diff * m_dx * m_dy * m_dz);
                    watch.Record(54);
                }

#else
                watch.Record(4);
				SetPotentialByPoissonFFT_keepRhok(m_Vhart, m_rho, m_size_x, m_size_y, m_size_z, m_dx, m_dy, m_dz, m_work);
                watch.Record(50);
				//calculate HR rho
				{
					OneComplex* rhok = (OneComplex*)m_work;
					OneComplex* hr_rhok = (OneComplex*)(m_work + m_size_3d * 4);
					UpConvert_Kspace(rhok, m_size_x, m_size_y, m_size_z, hr_rhok, m_HR_ratio_x, m_HR_ratio_y, m_HR_ratio_z);
                    watch.Record(51);
                    const size_t hr_size_3d = m_size_3d * m_HR_ratio_x * m_HR_ratio_y * m_HR_ratio_z;
					OneComplex* hr_rho_c = hr_rhok + hr_size_3d;
					IFFT_3D(hr_rho_c, hr_rhok, m_size_x * m_HR_ratio_x, m_size_y * m_HR_ratio_y, m_size_z * m_HR_ratio_z);
                    watch.Record(52);
                    double sum_rho = 0.0;
					const double scaling_factor = (double)(m_HR_ratio_x * m_HR_ratio_y * m_HR_ratio_z);
					for (size_t i = 0; i < hr_size_3d; ++i) {
						m_hr_rho[i] = hr_rho_c[i].r * scaling_factor;
						sum_rho += hr_rho_c[i].r * scaling_factor;
					}
					printf("Rho(HR) = %f\n", sum_rho * m_dx * m_dy * m_dz / (double)(m_HR_ratio_x * m_HR_ratio_y * m_HR_ratio_z));
                    watch.Record(53);

					double diff = 0.0;
					for (int iz = 0; iz < m_size_z; ++iz) {
						for (int iy = 0; iy < m_size_y; ++iy) {
							for (int ix = 0; ix < m_size_x; ++ix) {
								int64_t i = ix + m_size_x * (iy + m_size_y * iz);
								int64_t oi = m_HR_ratio_x * (ix + m_size_x * (m_HR_ratio_y * (iy + m_size_y * m_HR_ratio_z * iz)));
								//int64_t oi = (ix + hr_ratio_x * size_x * (iy + hr_ratio_y * size_y * iz));
								diff += (m_hr_rho[oi] - m_rho[i]) * (m_hr_rho[oi] - m_rho[i]);
							}
						}
					}
					printf("diff_Rho(HR) = %f\n", diff * m_dx * m_dy * m_dz);
                    watch.Record(54);
				}


				//calculate HR Vhart
				{
					OneComplex* vhartk = (OneComplex*)(m_work + m_size_3d * 2);
					OneComplex* hr_vhartk = (OneComplex*)(m_work + m_size_3d * 4);
					UpConvert_Kspace(vhartk, m_size_x, m_size_y, m_size_z, hr_vhartk, m_HR_ratio_x, m_HR_ratio_y, m_HR_ratio_z);
                    watch.Record(51);
                    const size_t hr_size_3d = m_size_3d * m_HR_ratio_x * m_HR_ratio_y * m_HR_ratio_z;
					OneComplex* hr_vhart_c = hr_vhartk + hr_size_3d;
					IFFT_3D(hr_vhart_c, hr_vhartk, m_size_x * m_HR_ratio_x, m_size_y * m_HR_ratio_y, m_size_z * m_HR_ratio_z);
                    watch.Record(52);
                    //double sum_rho = 0.0;
					const double scaling_factor = (double)(m_HR_ratio_x * m_HR_ratio_y * m_HR_ratio_z);
					for (size_t i = 0; i < hr_size_3d; ++i) {
						m_hr_Vhart[i] = hr_vhart_c[i].r * scaling_factor;                    
					}
                    watch.Record(53);

					double diff = 0.0;
					for (int iz = 0; iz < m_size_z; ++iz) {
						for (int iy = 0; iy < m_size_y; ++iy) {
							for (int ix = 0; ix < m_size_x; ++ix) {
								int64_t i = ix + m_size_x * (iy + m_size_y * iz);
								int64_t oi = m_HR_ratio_x * (ix + m_size_x * (m_HR_ratio_y * (iy + m_size_y * m_HR_ratio_z * iz)));
								//int64_t oi = (ix + hr_ratio_x * size_x * (iy + hr_ratio_y * size_y * iz));
								diff += (m_hr_Vhart[oi] - m_Vhart[i]) * (m_hr_Vhart[oi] - m_Vhart[i]);
							}
						}
					}
					printf("diff_Vhart(HR) = %f\n", diff * m_dx * m_dy * m_dz);
                    watch.Record(54);
				}
#endif

			}

		} else {//KohnSham_AE//
            /*
            if (is_calc_Vhart_xc) {
				SetPotentialByPoissonFFT(m_Vhart, m_rho, m_size_x, m_size_y, m_size_z, m_dx, m_dy, m_dz, m_work);				 
                //m_pcc_rho is set to zero//
                mSetPotentialVxc(m_Vx, m_Vc, m_rho, m_pcc_rho);
            }

			for (size_t i = 0; i < m_size_3d; ++i) {
				m_Vtot[i] = m_Vext[i] + m_Vhart[i] + m_Vx[i] + m_Vc[i];
			}

			mCheckPotentialVhart(m_Vhart, m_rho);
            */
		}
	}

    if (m_hamiltonian_type == HAMILTONIAN::KohnSham_PP) {//DFT//
        if (IsRoot(m_same_ddm_place_comm)) {
            if (is_calc_Vhart_xc) {
                mScatterField(ml_Vhart, m_Vhart);
                mScatterField_HR(ml_hr_rho, m_hr_rho);
                mScatterField_HR(ml_hr_Vhart, m_hr_Vhart);

                watch.Record(57);

                if (!is_spin_on) {
                    //Vxcの計算にのみPCC chargeを加味する//                    
                    mSetPotentialVxc_local(ml_Vx, ml_Vc, ml_rho, ml_pcc_rho);                    
                } else {
                    //Vxcの計算にのみPCC chargeを加味する//		
                    mSetPotentialVxcSpin_local(ml_Vx, ml_Vc, ml_Vx_down, ml_Vc_down, ml_rho, ml_pcc_rho, ml_rho_diff);                    
                }
                watch.Record(55);
            }


            const size_t local_size = ml_grid.Size3D();
            for (size_t i = 0; i < local_size; ++i) {
                ml_Vtot[i] = ml_Vext[i] + ml_Vhart[i] + ml_Vx[i] + ml_Vc[i];
            }
            MPI_Bcast(ml_Vtot, (int)local_size, MPI_DOUBLE, 0, m_same_ddm_place_comm);


            if (is_spin_on) {
                for (size_t i = 0; i < local_size; ++i) {
                    ml_Vtot_down[i] = ml_Vext[i] + ml_Vhart[i] + ml_Vx_down[i] + ml_Vc_down[i];
                }
                MPI_Bcast(ml_Vtot_down, (int)local_size, MPI_DOUBLE, 0, m_same_ddm_place_comm);

            }
            watch.Record(56);
        } else {

            const size_t local_size = ml_grid.Size3D();
            MPI_Bcast(ml_Vtot, (int)local_size, MPI_DOUBLE, 0, m_same_ddm_place_comm);
            if (is_spin_on) {
                MPI_Bcast(ml_Vtot_down, (int)local_size, MPI_DOUBLE, 0, m_same_ddm_place_comm);

            }
            watch.Record(56);

        }
    }

}

inline
void QUMASUN_TD::mSetPotentialVext(RspaceFunc<double>& V, const Nucleus* nuclei, int num_nuclei) {
	
	for (size_t i = 0; i < m_size_3d; ++i) {
		V[i] = 0.0;
	}

	for (int ni = 0; ni < num_nuclei; ++ni) {
		const double Qe = m_nuclei_valence_elecron[ni];
		const double R0_x = nuclei[ni].Rx;
		const double R0_y = nuclei[ni].Ry;
		const double R0_z = nuclei[ni].Rz;
		for (int iz = 0; iz < m_size_z; ++iz) {
			const double rz2 = SQ(Length(m_dz * (double)iz - R0_z, m_box_z));
			for (int iy = 0; iy < m_size_y; ++iy) {
				const double ry2 = SQ(Length(m_dy * (double)iy - R0_y, m_box_y));
				for (int ix = 0; ix < m_size_x; ++ix) {
					const double rx2 = SQ(Length(m_dx * (double)ix - R0_x, m_box_x));
					const size_t i = ix + m_size_x * (iy + (m_size_y * iz));

					V[i] += std::max<double>(-Qe / sqrt(rx2 + ry2 + rz2), -1.0e4);

				}
			}
		}
	}

}

#if 0
inline
void QUMASUN_TD::mSetPotentialVhart(RspaceFunc<double>& Vhart, RspaceFunc<double>& rho) {
	SetPotentialByPoissonFFT(Vhart.Pointer(), rho.Pointer(),
		m_size_x, m_size_y, m_size_z, m_dx, m_dy, m_dz, m_work);
}

inline
void QUMASUN_TD::mSetPotentialVlocalFromCorrespondingRho(RspaceFunc<double>& Vext, RspaceFunc<double>& nucl_rho) {
	SetPotentialByPoissonFFT(Vext.Pointer(), nucl_rho.Pointer(),
		m_size_x, m_size_y, m_size_z, m_dx, m_dy, m_dz, m_work);
}
#endif


inline
void QUMASUN_TD::mSetPotentialVxc(RspaceFunc<double>& Vx, RspaceFunc<double>& Vc, RspaceFunc<double>& rho, RspaceFunc<double>& pcc_rho) {

	for (size_t i = 0; i < m_size_3d; ++i) {
		auto ret = Calc_XC_LDA(rho[i] + pcc_rho[i]);
		Vx[i] = ret.V_x;
		Vc[i] = ret.V_c;
	}
}

inline
void QUMASUN_TD::mSetPotentialVxcSpin(double* Vx_up, double* Vc_up, double* Vx_down, double* Vc_down, const double* rho, const double* pcc_rho, const double* rho_diff) {
	for (size_t i = 0; i < m_size_3d; ++i) {
		auto ret = Calc_XC_LSDA(rho[i] + pcc_rho[i], rho_diff[i] / (rho[i] + pcc_rho[i]));
		Vx_up[i] = ret.V_x_up;
		Vc_up[i] = ret.V_c_up;
		Vx_down[i] = ret.V_x_down;
		Vc_down[i] = ret.V_c_down;
	}

}

inline
void QUMASUN_TD::mSetPotentialVxc_local(double* Vx, double* Vc, const double* rho, const double* pcc_rho) {
    const size_t size_3d = ml_grid.Size3D();
    for (size_t i = 0; i < size_3d; ++i) {
        auto ret = Calc_XC_LDA(rho[i] + pcc_rho[i]);
        Vx[i] = ret.V_x;
        Vc[i] = ret.V_c;
    }
}

inline
void QUMASUN_TD::mSetPotentialVxcSpin_local(double* Vx_up, double* Vc_up, double* Vx_down, double* Vc_down, const double* rho, const double* pcc_rho, const double* rho_diff) {
    const size_t size_3d = ml_grid.Size3D();
    for (size_t i = 0; i < size_3d; ++i) {
        auto ret = Calc_XC_LSDA(rho[i] + pcc_rho[i], rho_diff[i] / (rho[i] + pcc_rho[i]));
        Vx_up[i] = ret.V_x_up;
        Vc_up[i] = ret.V_c_up;
        Vx_down[i] = ret.V_x_down;
        Vc_down[i] = ret.V_c_down;
    }

}


inline
void QUMASUN_TD::mCheckPotentialVhart(RspaceFunc<double>& Vhart, RspaceFunc<double>& rho) {
	//RspaceFunc<double>& d2Vdx2 = m_work_rd;////*m_Vext;//
	RspaceFunc<double> rho2(m_size_3d);	
	mSetLaplasianFFT(rho2, Vhart, 1.0/(4.0*M_PI), m_work);


	printf("out: iz=  , Vhart, -(1/4pi)d2Vh/dx^2, rho, (rho2:byIFFT)\n");

	double diff = 0.0;
	double var = 0.0;
	double rho_tot = 0.0;
	double V_tot = 0.0;
	double d2V_tot = 0.0;
	double Ehart = 0.0;
	double Ehart2 = 0.0;

	double delta_rho = 0.0;
	double delta2_rho = 0.0;
	
	for (int iz = 0; iz < m_size_z; ++iz) {
		const int idz_m = ((iz == 0) ? m_size_z - 1 : - 1) * m_size_x * m_size_y;
		const int idz_p = ((iz == m_size_z - 1) ? 1 - m_size_z : 1) * m_size_x * m_size_y;
		for (int iy = 0; iy < m_size_y; ++iy) {
			const int idy_m = ((iy == 0) ? m_size_y - 1 : -1) * m_size_x;
			const int idy_p = ((iy == m_size_y - 1) ? 1 - m_size_y : 1) * m_size_x;
			for (int ix = 0; ix < m_size_x; ++ix) {
				const int idx_m = ((ix == 0) ? m_size_x - 1 : -1);
				const int idx_p = ((ix == m_size_x - 1) ? 1 - m_size_x : 1);
				const size_t i = ix + m_size_x * (iy + (m_size_y * iz));

				const double d2Vdx2 = (Vhart[i + idx_p] + Vhart[i + idx_m] - 2.0 * Vhart[i]) / (m_dx*m_dx)
					+ (Vhart[i + idy_p] + Vhart[i + idy_m] - 2.0 * Vhart[i]) / (m_dy*m_dy)
					+ (Vhart[i + idz_p] + Vhart[i + idz_m] - 2.0 * Vhart[i]) / (m_dz*m_dz);
				const double rh = rho[i];
				
				double d = -d2Vdx2 / (4.0 * M_PI);//
				diff += d - rh;
				var += (d - rh)*(d - rh);
				rho_tot += rh;
				double vh = Vhart[i];
				V_tot += Vhart[i];
				Ehart += Vhart[i] * rh;
				Ehart2 += Vhart[i] * rho2[i];
				d2V_tot += d;

				if (d > 0.0) {
					if (iz >1) {
						rho_tot += rh*0.;
					}
				}

				if (ix == m_size_x / 2) {
					if (iy == m_size_y / 2) {
						//if (iz == m_size_z / 2) 
						{
							double d2V_dx2 = -(Vhart[i + idx_p] + Vhart[i + idx_m] - 2.0 * Vhart[i]) / (m_dx * m_dx) / (4.0 * M_PI);
							double d2V_dy2 = -(Vhart[i + idy_p] + Vhart[i + idy_m] - 2.0 * Vhart[i]) / (m_dy * m_dy) / (4.0 * M_PI);
							double d2V_dz2 = -(Vhart[i + idz_p] + Vhart[i + idz_m] - 2.0 * Vhart[i]) / (m_dz * m_dz) / (4.0 * M_PI);
							

							printf("out: iz=%d, %f, %f, %f, (rho2:%f)\n", iz, vh, d, m_rho[i], rho2[i]);

						}
					}
				}

				double rh2 = rho2[i];
				double difrho = rh2;// -rh;
				delta_rho += difrho;
				delta2_rho += difrho * difrho;
			}
		}
	}
	diff *= m_dx * m_dy * m_dz;
	var = sqrt(var) * m_dx * m_dy * m_dz;
	rho_tot *= m_dx * m_dy * m_dz;
	V_tot *= m_dx * m_dy * m_dz;
	d2V_tot *= m_dx * m_dy * m_dz;
	Ehart *= m_dx * m_dy * m_dz * 0.5;
	Ehart2 *= m_dx * m_dy * m_dz * 0.5;
	delta_rho *= m_dx * m_dy * m_dz;
	delta2_rho = sqrt(delta2_rho)*m_dx * m_dy * m_dz;
	

	const double Ehart3 = mGetEnergyHartree();

	printf("\\int {-d2Vdx2 / (4.0*M_PI) - rho}dx = %f, %f\n", diff, var);
	printf("\\int {-d2Vdx2 / (4.0*M_PI) }dx = %f\n", d2V_tot);
	printf("\\int rho dx = %f\n", rho_tot);
	printf("\\int Vhart dx = %f\n", V_tot);
	printf("\\int (1/2)Vhart * rho dx = %f, %f, %f\n", Ehart, Ehart2, Ehart3);
	printf("\\int rho2-rho dx = %f, %f\n", delta_rho, delta2_rho);


}

inline
void QUMASUN_TD::mSetLaplasianFFT(RspaceFunc<double>& out, RspaceFunc<double>& in, double coef, double* work) {
	dcomplex* fk = (dcomplex*)(work);

#undef USE_R2C_C2R
#ifndef USE_R2C_C2R

	dcomplex* rho_comp = (dcomplex*)(work + m_size_3d * 2);
	for (size_t i = 0; i < m_size_3d; ++i) {
		rho_comp[i] = { in[i],0.0 };
	}

	FFT_3D(fk, rho_comp, m_size_x, m_size_y, m_size_z);

	const int kz_end = m_size_z;
	const int kx_end = m_size_x;

#else
	/*
	problem: this path cannot obtain correct Vhart
	*/
	RspaceFunc<double>& rho_d = work + m_size_3d * 2;
	vecmath::Copy<double>(rho_d, in, m_size_3d);
	/*fftw_plan_dft_r2c_3d*/
	FFT_3D(fk, rho_d, m_size_x, m_size_y, m_size_z);

	const int kz_end = m_size_z;
	const int kx_end = m_size_x / 2 + 1;
#endif

	//V = c rho, where c = 4pi / (kx^2 + ky^2 + kz^2) //
	//  = 1.0/ ((pi*kx^2/box_x^2) + (pi*ky^2/box_y^2) + (pi*kz^2/box_z^2)) 
	//const double coef_x = M_PI / (m_box_x * m_box_x);
	//const double coef_y = M_PI / (m_box_y * m_box_y);
	//const double coef_z = M_PI / (m_box_z * m_box_z);
	const double coef_x = SQ(2.0 * M_PI / m_box_x);
	const double coef_y = SQ(2.0 * M_PI / m_box_y);
	const double coef_z = SQ(2.0 * M_PI / m_box_z);

	for (int kz = 0; kz < kz_end; ++kz) {
		const double kz2 = coef_z * SQ((double)(kz < m_size_z / 2 ? kz : kz - m_size_z));
		for (int ky = 0; ky < m_size_y; ++ky) {
			const double ky2 = coef_y * SQ((double)(ky < m_size_y / 2 ? ky : ky - m_size_y));

			for (int kx = 0; kx < kx_end; ++kx) {
				const double kx2 = coef_x * SQ((double)(kx < m_size_x / 2 ? kx : kx - m_size_x));
				const size_t i = kx + kx_end * (ky + (m_size_y * kz));

				fk[i] *= coef *(kx2 + ky2 + kz2);
			}

		}
	}

#ifndef USE_R2C_C2R 

	dcomplex* Fr = (dcomplex*)(work);
	/*
	for (size_t i = 0; i < m_size_3d; ++i) {
		Fr[i] = { 0.0,0.0 };
	}
	*/
	IFFT_3D(Fr, fk, m_size_x, m_size_y, m_size_z);
	
	for (size_t i = 0; i < m_size_3d; ++i) {
		out[i] = Fr[i].real();
		/*
		if (i / (128 * 128) == 64) {
			if ((i / 128) % 128 == 64) {
				printf("out: %d, %f, (ref:%f)\n", i % 128, out[i], m_rho[i]);
			}
		}*/
	}
	
#else


	IFFT_3D(Fr, fk, m_size_x, m_size_y, m_size_z);


#endif
}
