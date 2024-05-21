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


/*
* Perform the following calculations
* - Vext or Vlocal potential due to nuclei
* - Core-Core repulsive energy
* where they only need to be calculated once before the SCF calculation
*/
inline
void QUMASUN_SINGLE::mPrepareCore() {

	//Vext or Vlocal/////////////////////////////////////////
	vecmath::SetZero(m_pcc_rho.Pointer(), m_size_3d);

	if (m_hamiltonian_type == HAMILTONIAN::KohnSham_PP) {//DFT//
		mSetPotentialVlocalFromCorrespondingRho(m_Vext, m_nucl_rho, m_nuclei, m_num_nuclei);
		m_pp_integrator.SetPccCharge(m_pcc_rho, m_nuclei, m_num_nuclei);

	} else {
		mSetPotentialVext(m_Vext, m_nuclei, m_num_nuclei);
	}

#if 0
	{//for debug
		
		FILE* fp = fopen("Vext.txt","w");
		fprintf(fp, "#r\tVext\n");
		for (int ix = 0; ix < m_size_x; ++ix) {
			size_t i = ix + m_size_x * (ix + m_size_y * ix);
			const double x = (double)ix * m_dx;
			const double dx = m_nuclei[0].Rx - x;
			const double dy = m_nuclei[0].Ry - x;
			const double dz = m_nuclei[0].Rz - x;
			const double r = sqrt(dx * dx + dy * dy + dz * dz);
						
			fprintf(fp, "%.12f\t%.12f\n", r, m_Vext[i]);
		}
		fclose(fp);
	}
#endif


	//Valence Electron/////////////////////////////////////////
	{

		for (int i = 0; i < m_num_nuclei; ++i) {
			m_nuclei_valence_elecron[i] = m_pp_integrator.NumValenceElectron(m_nuclei[i].Z);
		}
	}

	//E_core_core//////////////////////////////////////////////
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

			//自己相互作用の補正：原子核自身の発生するVlcoalの原点での相互作用
			const double Vlocal_0 = m_pp_integrator.GetVlocalRadial(m_nuclei[i].Z, 0.0);
			E_self += -(Qi * (-Vlocal_0));


			//近接の相互作用の補正：CoulombとVlocalの差
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
	}

}
