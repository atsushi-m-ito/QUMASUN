#pragma once
//#define _USE_MATH_DEFINES
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include "qumasun_single.h"
#include "wrap_fft.h"
#include "Vxc.h"
#include "nucleus.h"
#include "vecmath.h"



//#define USE_R2C_C2R


/*
* V_hart, Vx, Vc V_totを計算する
* 
* V_extはmPrepareCoreで最初に一度だけ計算
*/
inline
void QUMASUN_SINGLE::mSetPotential() {


	if (m_hamiltonian_type == HAMILTONIAN::Schrodinger) {
		mSetPotentialVhart(m_Vhart, m_rho);
		for (size_t i = 0; i < m_size_3d; ++i) {
			m_Vtot[i] = m_Vext[i] + m_Vhart[i];
		}

	} else if(m_hamiltonian_type == HAMILTONIAN::KohnSham_PP) {//DFT//

		mSetPotentialVhart(m_Vhart, m_rho);
		
		if (!is_spin_on) {
			//Vxcの計算にのみPCC chargeを加味する//		
			mSetPotentialVxc(m_Vx, m_Vc, m_rho, m_pcc_rho);

			for (size_t i = 0; i < m_size_3d; ++i) {
				m_Vtot[i] = m_Vext[i] + m_Vhart[i] + m_Vx[i] + m_Vc[i];
			}
		} else {
			//Vxcの計算にのみPCC chargeを加味する//		
			mSetPotentialVxcSpin(m_Vx, m_Vc, m_Vx_down, m_Vc_down, m_rho, m_pcc_rho, m_rho_diff);
			
			for (size_t i = 0; i < m_size_3d; ++i) {
				m_Vtot[i] = m_Vext[i] + m_Vhart[i] + m_Vx[i] + m_Vc[i];
			}
			for (size_t i = 0; i < m_size_3d; ++i) {
				m_Vtot_down[i] = m_Vext[i] + m_Vhart[i] +m_Vx_down[i] + m_Vc_down[i];
			}
		}
		//mCheckPotentialVhart(m_Vhart, m_rho);
	}else{//KohnSham_AE//

		mSetPotentialVhart(m_Vhart, m_rho);
		//m_pcc_rho is set to zero//
		mSetPotentialVxc(m_Vx, m_Vc, m_rho, m_pcc_rho);

		for (size_t i = 0; i < m_size_3d; ++i) {
			m_Vtot[i] = m_Vext[i] + m_Vhart[i] + m_Vx[i] + m_Vc[i];
		}

		mCheckPotentialVhart(m_Vhart, m_rho);
	}
}

inline
void QUMASUN_SINGLE::mSetPotentialVext(RspaceFunc<double>& V, const Nucleus* nuclei, int num_nuclei) {
	
	for (size_t i = 0; i < m_size_3d; ++i) {
		V[i] = 0.0;
	}

	for (int ni = 0; ni < num_nuclei; ++ni) {
		const double Qe = (double)m_nuclei_valence_elecron[ni];
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


inline
void QUMASUN_SINGLE::mSetPotentialVhart(RspaceFunc<double>& Vhart, RspaceFunc<double>& rho) {



	KspaceFunc<dcomplex> rhok((dcomplex*)m_work);

//#undef USE_R2C_C2R
#ifndef USE_R2C_C2R

	RspaceFunc<dcomplex> rho_comp((dcomplex*)(m_work + m_size_3d * 2));
	for (size_t i = 0; i < m_size_3d; ++i) {
		rho_comp[i] = { rho[i],0.0 };
	}

	FFT_3D(rhok, rho_comp, m_size_x, m_size_y, m_size_z);

	const int kz_end = m_size_z;
	const int kx_end = m_size_x;

#else
	/*
	problem: this path cannot obtain correct Vhart
	*/
	RspaceFunc<double>& rho_d = m_work_rd;

	{
		vecmath::Copy<double>(rho_d, rho, m_size_3d);
		//rho_d[m_size_x / 2 + m_size_x * (m_size_y / 2 + m_size_y * (m_size_z / 2))] -= 1.0 / (m_dx * m_dy * m_dz);
	}
	
	double sumrho = 0.0;
	for (size_t i = 0; i < m_size_3d; ++i) {
		sumrho += rho_d[i];
	}
	sumrho *= m_dx * m_dy * m_dz;
	printf("rho in poisson = %f\n", sumrho);

	/*fftw_plan_dft_r2c_3d*/
	FFT_3D(rhok, rho_d, m_size_x, m_size_y, m_size_z);

	const int kz_end = m_size_z;
	const int kx_end = m_size_x/2+1;
#endif

	//NOTE: 
	//  V = c rho, where c = 4pi / (kx^2 + ky^2 + kz^2) //
	//    = 1.0/ ((pi*kx^2/box_x^2) + (pi*ky^2/box_y^2) + (pi*kz^2/box_z^2)) 
	
	const double coef1_x = (2.0 * M_PI / m_box_x);
	const double coef1_y = (2.0 * M_PI / m_box_y);
	const double coef1_z = (2.0 * M_PI / m_box_z);

//#define NO_FOLD

	for (int kz = 0; kz < kz_end; ++kz) {
#ifdef NO_FOLD
		const double kz2 = coef_z * SQ((double)(kz));
#else
		const double kz2 = SQ(coef1_z * (double)(kz < m_size_z / 2 ? kz : kz - m_size_z));
#endif
		for (int ky = 0; ky < m_size_y; ++ky) {
#ifdef NO_FOLD
			const double ky2 = coef_z * SQ((double)(ky));
#else
			const double ky2 = SQ(coef1_y * (double)(ky < m_size_y / 2 ? ky : ky - m_size_y));
#endif


			for (int kx = 0; kx < kx_end; ++kx) {
#ifdef NO_FOLD
				const double kx2 = coef_z * SQ((double)(kx));
#else
				const double kx2 = SQ(coef1_x * (double)(kx < m_size_x / 2 ? kx : kx - m_size_x));
#endif
				const size_t i = kx + kx_end * (ky + (m_size_y * kz));

				if (kx == 0 && ky == 0 && kz == 0) {
					rhok[i] = 0.0;
				} else {
					rhok[i] *= 4.0 * M_PI / (kx2 + ky2 + kz2);
				}
			}

		}
	}
	
#ifndef USE_R2C_C2R 

	RspaceFunc<dcomplex> Vh((dcomplex*)(m_work + m_size_3d * 2));

	IFFT_3D(Vh, rhok, m_size_x, m_size_y, m_size_z);
	
	for (size_t i = 0; i < m_size_3d; ++i) {
		Vhart[i] = Vh[i].real();
	}
#else

	
	IFFT_3D(Vhart, rhok, m_size_x, m_size_y, m_size_z);


#endif



}

inline
void QUMASUN_SINGLE::mSetPotentialVlocalFromCorrespondingRho(RspaceFunc<double>& Vext, RspaceFunc<double>& nucl_rho, const Nucleus* nuclei, int num_nuclei) {

	if (m_hamiltonian_type != HAMILTONIAN::KohnSham_PP) {
		return;
	}
	
	KspaceFunc<dcomplex> rhok((dcomplex*)m_work);
	RspaceFunc<double> work_rd2(m_work + m_size_3d * 2);
	
	{
		for (size_t i = 0; i < m_size_3d; ++i) {
			nucl_rho[i] = 0.0;
		}
		
		m_pp_integrator.SetChargeVlocal(nucl_rho, work_rd2, nuclei, num_nuclei);
		double sumrho = 0.0;
		for (size_t i = 0; i < m_size_3d; ++i) {
			sumrho += nucl_rho[i];
		}
		sumrho *= m_dx * m_dy * m_dz;

		printf("ChargeForVlocal=%f\n", sumrho);

	}

#ifdef USE_R2C_C2R
	/*fftw_plan_dft_r2c_3d*/
	FFT_3D(rhok, rho_d, m_size_x, m_size_y, m_size_z);

	const int kz_end = m_size_z;
	const int kx_end = m_size_x / 2 + 1;
#else

	RspaceFunc<dcomplex> rho_r((dcomplex*)(m_work + m_size_3d * 3));
	for (size_t i = 0; i < m_size_3d; ++i) {
		rho_r[i] = {nucl_rho[i],0.0};
	}
	/*fftw_plan_dft_3d for complex*/
	FFT_3D(rhok, rho_r, m_size_x, m_size_y, m_size_z);

	const int kz_end = m_size_z;
	const int kx_end = m_size_x;
#endif

	//V = c rho, where c = 4pi / (kx^2 + ky^2 + kz^2) //
	//  = 1.0/ ((pi*kx^2/box_x^2) + (pi*ky^2/box_y^2) + (pi*kz^2/box_z^2)) 

	const double coef1_x = (2.0 * M_PI / m_box_x);
	const double coef1_y = (2.0 * M_PI / m_box_y);
	const double coef1_z = (2.0 * M_PI / m_box_z);

	//#define NO_FOLD

	for (int kz = 0; kz < kz_end; ++kz) {
#ifdef NO_FOLD
		const double kz2 = coef_z * SQ((double)(kz));
#else
		//const double kz2 = coef_z * SQ((double)(kz<m_size_z / 2 ? kz : kz - m_size_z));
		const double kz2 = SQ(coef1_z * (double)(kz < m_size_z / 2 ? kz : kz - m_size_z));
#endif
		for (int ky = 0; ky < m_size_y; ++ky) {
#ifdef NO_FOLD
			const double ky2 = coef_z * SQ((double)(ky));
#else
			//const double ky2 = coef_y * SQ((double)(ky<m_size_y / 2 ? ky : ky - m_size_y));
			const double ky2 = SQ(coef1_y * (double)(ky < m_size_y / 2 ? ky : ky - m_size_y));
#endif


			for (int kx = 0; kx < kx_end; ++kx) {
#ifdef NO_FOLD
				const double kx2 = coef_z * SQ((double)(kx));
#else
				const double kx2 = SQ(coef1_x * (double)(kx < m_size_x / 2 ? kx : kx - m_size_x));
#endif
				const size_t i = kx + kx_end * (ky + (m_size_y * kz));
				//const size_t i = kz + kz_end * (ky + (m_size_y * kx));

				//rhok[i] /= (kx2 + ky2 + kz2);
				if (kx == 0 && ky == 0 && kz == 0) {
					rhok[i] = 0.0;
				} else {
					rhok[i] *= 4.0 * M_PI / (kx2 + ky2 + kz2);
				}
			}
		}
	}


#ifdef USE_R2C_C2R
	IFFT_3D(Vext, rhok, m_size_x, m_size_y, m_size_z);
#else
	IFFT_3D(rho_r, rhok, m_size_x, m_size_y, m_size_z);
	for (size_t i = 0; i < m_size_3d; ++i) {
		Vext[i] = rho_r[i].real();
	}
#endif

}


inline
void QUMASUN_SINGLE::mSetPotentialVxc(RspaceFunc<double>& Vx, RspaceFunc<double>& Vc, RspaceFunc<double>& rho, RspaceFunc<double>& pcc_rho) {

	for (size_t i = 0; i < m_size_3d; ++i) {
		auto ret = Calc_XC_LDA(rho[i] + pcc_rho[i]);	
		Vx[i] = ret.V_x;
		Vc[i] = ret.V_c;
	}

//#define DEBUG_PRINT
#ifdef DEBUG_PRINT
	auto ret = Calc_XC_LDA(rho[0] + pcc_rho[0]);
	printf("PCC1: %f, %f, %f (%f, %f), %f (%f, %f)\n",
		rho[0], pcc_rho[0], ret.E_den_x + ret.E_den_c, ret.E_den_x, ret.E_den_c,
		ret.V_x + ret.V_c, ret.V_x, ret.V_c);
#endif
//#undef DEBUG_PRINT
}

inline
void QUMASUN_SINGLE::mSetPotentialVxcSpin(double* Vx_up, double* Vc_up, double* Vx_down, double* Vc_down, const double* rho, const double* pcc_rho, const double* rho_diff) {
	for (size_t i = 0; i < m_size_3d; ++i) {
		auto ret = Calc_XC_LSDA(rho[i] + pcc_rho[i], rho_diff[i] / (rho[i] + pcc_rho[i]));
		Vx_up[i] = ret.V_x_up;
		Vc_up[i] = ret.V_c_up;
		Vx_down[i] = ret.V_x_down;
		Vc_down[i] = ret.V_c_down;
	}

}

inline
void QUMASUN_SINGLE::mCheckPotentialVhart(RspaceFunc<double>& Vhart, RspaceFunc<double>& rho) {
	//RspaceFunc<double>& d2Vdx2 = m_work_rd;////*m_Vext;//
	RspaceFunc<double> rho2(m_size_3d);	
	mSetLaplasianFFT(rho2, Vhart, 1.0/(4.0*M_PI));


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
void QUMASUN_SINGLE::mSetLaplasianFFT(RspaceFunc<double>& out, RspaceFunc<double>& in, double coef) {
	KspaceFunc<dcomplex> fk((dcomplex*)m_work);

#undef USE_R2C_C2R
#ifndef USE_R2C_C2R

	RspaceFunc<dcomplex> rho_comp((dcomplex*)(m_work + m_size_3d*2));
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
	RspaceFunc<double>& rho_d = m_work_rd;
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

	RspaceFunc<dcomplex>Fr((dcomplex*)(m_work + m_size_3d * 2));
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
