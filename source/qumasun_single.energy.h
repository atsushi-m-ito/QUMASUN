#pragma once
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include "qumasun_single.h"
#include "Vxc.h"
#include "vecmath.h"
#include "field_interpolation.h"
#include "GridDifference2nd.h"
#include "GridDifference4th.h"

inline
double QUMASUN_SINGLE::mGetTotalEnergy() {
	
	const double Ekin = mGetEnergyKinetic();
	const double Eext = mGetEnergyExtByVextRho(); 
	const double Ehart = mGetEnergyHartree();
	const double Eext_nucl_density = mGetEnergyExtByVhartNuclRho();
	const double Eext_nucl_point = mGetEnergyExtByVhartAtPoint();
	const double E_corecore_nucl_density = mGetEnergyCoreCoreNuclDensity();
	double Ex = 0.0;
	double Ec = 0.0;
	double VxRho = 0.0;
	double VcRho = 0.0;
	double VxRho_down = 0.0;
	double VcRho_down = 0.0;

	double Exc;
	if (!is_spin_on) {
		Exc=mGetEnergyXC(&Ex, &Ec, &VxRho, &VcRho);
	} else {
		Exc=mGetEnergyXCSpin(&Ex, &Ec, &VxRho, &VcRho, &VxRho_down, &VcRho_down);
	}
	double E_pp = 0.0;
	if (m_hamiltonian_type == HAMILTONIAN::KohnSham_PP) {
		E_pp = mGetEnergyPseudoNonlocal();
	}

	const double E_corecore_direct = mGetEnergyCoreCore_direct();
	const double E_corecore_Vext = mGetEnergyCoreCoreByVextAtPoint();

	const double Etot = Ekin + Eext + Ehart + Exc + E_pp + E_corecore_direct;
	const double Etot2 = Ekin + Eext + Ehart + Exc + E_pp + E_corecore_Vext;
	const double Etot3 = Ekin + Eext + Ehart + Exc + E_pp + E_corecore_Vext + m_diff_Coulomb_Vlocal;
	const double Etot4 = Ekin + Eext + Ehart + Exc + E_pp + E_corecore_Vext + m_Ecore_self + m_diff_Coulomb_Vlocal;
	const double Etot5 = Ekin + Eext + Ehart + Exc + E_pp + E_corecore_nucl_density + m_Ecore_self + m_diff_Coulomb_Vlocal;
	printf("Etot = %f  (test: %f, %f, %f, %f)\n", Etot5, Etot, Etot2, Etot3, Etot4);
	printf("Ekin = %f\n", Ekin);
	if (m_hamiltonian_type == HAMILTONIAN::KohnSham_PP) {
		printf("E_PP_local = %f\n", Eext);
		printf("  Eext_nucl_density = %f (non-use)\n", Eext_nucl_density);
		printf("  Eext_nucl_point = %f (non-use)\n", Eext_nucl_point);
		printf("E_PP_nonlocal = %f\n", E_pp);	
	} else {
		printf("Eext = %f\n", Eext);
	}
	printf("Ehart = %f\n", Ehart);
	printf("Exc = %f\n", Exc);
	printf("  Ex = %f\n", Ex);
	printf("  Ec = %f\n", Ec);
	if (!is_spin_on) {
		printf("  VxRho = %f (non-use)\n", VxRho);
		printf("  VcRho = %f (non-use)\n", VcRho);
	} else {
		printf("  VxRho = %f, %f (non-use)\n", VxRho, VxRho_down);
		printf("  VcRho = %f, %f (non-use)\n", VcRho, VcRho_down);
	}
	printf("Enn(Vlocal*rho_nucl-self) = %f (non-use)\n", E_corecore_nucl_density + m_Ecore_self);
	printf("  Enn(Vlocal*rho_nucl) = %f\n", E_corecore_nucl_density);
	printf("  Enn(self) = %f\n", -m_Ecore_self);
	printf("  Enn(Vlocal) = %f (non-use)\n", E_corecore_Vext);
	printf("  Enn(simpleCoulomb) = %f (non-use)\n", E_corecore_direct);
	printf("correct_Coulomb_Vlocal = %f\n", m_diff_Coulomb_Vlocal);
	return Etot;
}


inline
double QUMASUN_SINGLE::mGetEnergyKinetic() {
	double EK = 0.0;
	for (int s = 0; s < num_solution; ++s) {
		if (m_occupancy[s] <= 0.0) continue;
		EK += mGetOneEnergyKinetic(m_psi_set[s]) * m_occupancy[s];
	}
	
	if (is_spin_on) {
		for (int s = 0; s < num_solution; ++s) {
			if (m_occupancy_down[s] <= 0.0) continue;
			EK += mGetOneEnergyKinetic(m_psi_set[s+num_solution]) * m_occupancy_down[s];
		}
	}

	return EK;
}

inline
double QUMASUN_SINGLE::mGetOneEnergyKinetic(RspaceFunc<double >& p) {
	RspaceFunc<double > Hp(m_work);
	vecmath::SetZero<double>(Hp, m_size_3d);
	mKineticMatrixAdd(Hp, p);
	return vecmath::InnerProd<double>(p, Hp, m_size_3d) * m_dx*m_dy*m_dz;
}






/*
* calculate
* -\int n_nucl V_hart(\rho) dx
*/
inline
double QUMASUN_SINGLE::mGetEnergy_V_rho(const double* V, const double* rho) {
	double sum2 = 0.0;
	for (int i = 0; i < m_size_3d; ++i) {
		sum2 += V[i] * rho[i];
	}	
	return sum2 * m_dx * m_dy * m_dz;
}

inline
double QUMASUN_SINGLE::mGetEnergyExtByVextRho() {
	const double EV = mGetEnergy_V_rho(m_Vext.Pointer(), m_rho.Pointer());
	return EV;
}

/*
* calculate
* (1/2) \int \rho V_hart(\rho) dx
*/
inline
double QUMASUN_SINGLE::mGetEnergyHartree() {
	const double EV = mGetEnergy_V_rho(m_Vhart.Pointer(), m_rho.Pointer());
	return EV/2.0;
}

/*
* calculate
* -\int n_nucl V_hart(\rho) dx
*/
inline
double QUMASUN_SINGLE::mGetEnergyExtByVhartNuclRho() {
	const double EV = mGetEnergy_V_rho(m_Vhart.Pointer(), m_nucl_rho.Pointer());
	return EV;
}

inline
double QUMASUN_SINGLE::mGetEnergyCoreCoreNuclDensity() {
	const double EV = mGetEnergy_V_rho(m_Vext.Pointer(), m_nucl_rho.Pointer());
	return EV / 2.0;
}


inline
double QUMASUN_SINGLE::mGetEnergyXC(double* pEx, double* pEc, double* pVxRho, double* pVcRho) {

	RspaceFunc<double>& rho = m_rho;
	RspaceFunc<double>& pcc_rho = m_pcc_rho;
	
	double E_x = 0.0;
	double E_c = 0.0;
	double VxRho = 0.0;
	double VcRho = 0.0;
	
	for (size_t i = 0; i < m_size_3d; ++i) {
		//pcc_charge should be added//
		const double rho_i = rho[i] + pcc_rho[i];
		auto ret = Calc_XC_LDA(rho_i);
		
		E_x += ret.E_den_x * rho_i;
		E_c += ret.E_den_c * rho_i;
		VxRho += ret.V_x * rho_i;
		VcRho += ret.V_c * rho_i;
	}
	E_x *= m_dx * m_dy * m_dz;
	E_c *= m_dx * m_dy * m_dz;
	VxRho *= m_dx * m_dy * m_dz;
	VcRho *= m_dx * m_dy * m_dz;
	if (pVxRho)*pVxRho = VxRho;
	if (pVcRho)*pVcRho = VcRho;


	*pEx = E_x;
	*pEc = E_c;
	return E_x + E_c;
}

inline
double QUMASUN_SINGLE::mGetEnergyXCSpin(double* pEx, double* pEc, double* pVxRho_up, double* pVcRho_up, double* pVxRho_down, double* pVcRho_down) {

	RspaceFunc<double>& rho = m_rho;
	RspaceFunc<double>& pcc_rho = m_pcc_rho;

	double E_x = 0.0;
	double E_c = 0.0;
	double VxRho_up = 0.0;
	double VcRho_up = 0.0;
	double VxRho_down = 0.0;
	double VcRho_down = 0.0;

	for (size_t i = 0; i < m_size_3d; ++i) {
		//pcc_charge should be added//
		const double rho_i = rho[i] + pcc_rho[i];
		auto ret = Calc_XC_LSDA(rho_i, m_rho_diff[i]/rho_i);

		E_x += ret.E_den_x * rho_i;
		E_c += ret.E_den_c * rho_i;
		VxRho_up += ret.V_x_up * rho_i;
		VcRho_up += ret.V_c_up * rho_i;
		VxRho_down += ret.V_x_down * rho_i;
		VcRho_down += ret.V_c_down * rho_i;
	}
	E_x *= m_dx * m_dy * m_dz;
	E_c *= m_dx * m_dy * m_dz;
	VxRho_up *= m_dx * m_dy * m_dz;
	VcRho_up *= m_dx * m_dy * m_dz;
	VxRho_down *= m_dx * m_dy * m_dz;
	VcRho_down *= m_dx * m_dy * m_dz;
	if (pVxRho_up)*pVxRho_up = VxRho_up;
	if (pVcRho_up)*pVcRho_up = VcRho_up;
	if (pVxRho_down)*pVxRho_down = VxRho_down;
	if (pVcRho_down)*pVcRho_down = VcRho_down;


	*pEx = E_x;
	*pEc = E_c;
	return E_x + E_c;
}

inline
double QUMASUN_SINGLE::mGetEnergyPseudoNonlocal() {
	double EV = 0.0;
	{
		for (int s = 0; s < num_solution; ++s) {

			if (m_occupancy[s] <= 0.0) continue;

			double pp_ene =m_pp_integrator.EnergyPPnonlocal(m_psi_set[s], m_nuclei, m_num_nuclei);

			EV += pp_ene * m_occupancy[s];
		}		
	}

	if(is_spin_on) {
		for (int s = 0; s < num_solution; ++s) {

			if (m_occupancy_down[s] <= 0.0) continue;

			double pp_ene = m_pp_integrator.EnergyPPnonlocal(m_psi_set[s+num_solution], m_nuclei, m_num_nuclei);

			EV += pp_ene * m_occupancy_down[s];
		}
	}
	return EV;
}

inline
double QUMASUN_SINGLE::mGetEnergyCoreCore_direct() {
	double EV = 0.0;

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
		for (int k = i + 1; k < m_num_nuclei; ++k) {
			
			const double Qk = m_nuclei_valence_elecron[k];
			double xk = Folding(m_nuclei[k].Rx, m_box_x);
			double yk = Folding(m_nuclei[k].Ry, m_box_y);
			double zk = Folding(m_nuclei[k].Rz, m_box_z);
			double dx = Distance(xk - xi, m_box_x);
			double dy = Distance(yk - yi, m_box_y);
			double dz = Distance(zk - zi, m_box_z);
			double r = sqrt(dx * dx + dy * dy + dz * dz);

			EV += Qi * Qk / r;
		}
	}

	return EV;
}


inline
double QUMASUN_SINGLE::mGetEnergyPotentialAtNucl(const double* V) {
	double EV = 0.0;
	
	auto Folding = [](double x, double box_width) {
		return x - floor(x / box_width) * box_width;
		};

	for (int i = 0; i < m_num_nuclei; ++i) {

		const double Qi = m_nuclei_valence_elecron[i];
		const double x = Folding(m_nuclei[i].Rx, m_box_x);
		const double y = Folding(m_nuclei[i].Ry, m_box_y);
		const double z = Folding(m_nuclei[i].Rz, m_box_z);
		

		const double val_V_ext = Interpolation<double>(x / m_dx, y / m_dy, z / m_dz, V, { 0, 0, 0, m_size_x, m_size_y, m_size_z },
			{ m_size_x, m_size_y, m_size_z });

		EV += val_V_ext * (-Qi);

	}

	return EV;
}

inline
double QUMASUN_SINGLE::mGetEnergyCoreCoreByVextAtPoint() {
	double EV = mGetEnergyPotentialAtNucl(m_Vext.Pointer());	
	return EV/2.0;
}

inline
double QUMASUN_SINGLE::mGetEnergyExtByVhartAtPoint() {
	double EV = mGetEnergyPotentialAtNucl(m_Vhart.Pointer());
	return EV;
}
