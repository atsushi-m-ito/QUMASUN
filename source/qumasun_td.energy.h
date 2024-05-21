#pragma once
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include "qumasun_td.h"
#include "Vxc.h"
#include "vecmath.h"
#include "field_interpolation.h"
#include "GridFor.h"
#include "debug_print.h"


//partially mpi supported//
inline
double QUMASUN_TD::mGetTotalEnergy(bool with_nonlocal_force) {
	const bool is_root_global = IsRoot(m_mpi_comm);

	const int proc_id = GetProcessID(m_mpi_comm);
	const int num_procs = GetNumProcess(m_mpi_comm);

    DEBUG_PRINTF("[%d] mGetTotalEnergy 1\n", proc_id); 
	const double Ekin = mGetEnergyKinetic();
    watch.Record(60);
	DEBUG_PRINTF("[%d] mGetTotalEnergy 2\n", proc_id); 

	double E_pp = 0.0;
	if (m_hamiltonian_type == HAMILTONIAN::KohnSham_PP) {
        if (with_nonlocal_force) {
            m_force_nonlocal.resize(m_num_nuclei);
            E_pp = mGetForcePseudoNonlocal((double*)&m_force_nonlocal[0]);
            //mGetForcePseudoNonlocal2((double*)&forces_n_nonlocal[0]);
            watch.Record(70);
        }else{
            E_pp = mGetEnergyPseudoNonlocal();
            watch.Record(61);
        }
	}
    DEBUG_PRINTF("[%d] mGetTotalEnergy 3\n", proc_id); 

	if (IsRoot(m_same_ddm_place_comm)) {

        const double Eext = mGetEnergyExtByVextRho();
        const double Ehart = mGetEnergyHartree();
        const double Eext_nucl_density = mGetEnergyExtByVhartNuclRho();
        const double E_corecore_nucl_density = mGetEnergyCoreCoreNuclDensity();
        watch.Record(62);
        const double Eext_HR = mGetEnergyExtByVextRho_HR();
        const double Eext_nucl_density_HR = mGetEnergyExtByVhartNuclRho_HR();
        const double E_corecore_nucl_density_HR = mGetEnergyCoreCoreNuclDensity_HR();
        watch.Record(63);

		double Ex = 0.0;
		double Ec = 0.0;
		double VxRho = 0.0;
		double VcRho = 0.0;
		double VxRho_down = 1.0;
		double VcRho_down = 0.0;

		double Exc;
		if (!is_spin_on) {
			Exc = mGetEnergyXC(&Ex, &Ec, &VxRho, &VcRho);
		} else {
			Exc = mGetEnergyXCSpin(&Ex, &Ec, &VxRho, &VcRho, &VxRho_down, &VcRho_down);
		}
        watch.Record(64);

        DEBUG_PRINTF("[%d] mGetTotalEnergy 8\n", proc_id); 

		if (is_root_global) {

            DEBUG_PRINTF("[%d] mGetTotalEnergy 9\n", proc_id); 
			const double Eext_nucl_point = mGetEnergyExtByVhartAtPoint();
            watch.Record(65);
            DEBUG_PRINTF("[%d] mGetTotalEnergy 10\n", proc_id);
			const double E_corecore_direct = mGetEnergyCoreCore_direct();
            watch.Record(66);
            DEBUG_PRINTF("[%d] mGetTotalEnergy 11\n", proc_id); 
			const double E_corecore_Vext = mGetEnergyCoreCoreByVextAtPoint();
            watch.Record(65);
            DEBUG_PRINTF("[%d] mGetTotalEnergy 12\n", proc_id); 

			const double Etot = Ekin + Eext + Ehart + Exc + E_pp + E_corecore_direct;
            const double Etot2 = Ekin + Eext + Ehart + Exc + E_pp + E_corecore_Vext + m_diff_Coulomb_Vlocal;
            const double Etot3 = Ekin + Eext + Ehart + Exc + E_pp + E_corecore_Vext + m_Ecore_self + m_diff_Coulomb_Vlocal;
            const double Etot4 = Ekin + Eext + Ehart + Exc + E_pp + E_corecore_nucl_density + m_Ecore_self + m_diff_Coulomb_Vlocal;
            const double Etot5 = Ekin + Eext + Ehart + Exc + E_pp + E_corecore_nucl_density - m_Ecore_self_v2 + m_diff_Coulomb_Vlocal;
            const double Etot6 = Ekin + Eext + Ehart + Exc + E_pp + E_corecore_nucl_density - m_Ecore_self_v2 + m_Enn_close_correction;
            const double Etot7 = Ekin + Eext + Ehart + Exc + E_pp + E_corecore_nucl_density_HR - m_Ecore_self_HR + m_Enn_close_correction;
            const double Etot8 = Ekin + Eext_HR + Ehart + Exc + E_pp + E_corecore_nucl_density_HR - m_Ecore_self_HR + m_Enn_close_correction;
            printf("Etot = %f (test: %f, %f, %f, %f, %f)\n", Etot8, Etot7, Etot6, Etot2, Etot3, Etot4);
            printf("Ekin = %f\n", Ekin);
			if (m_hamiltonian_type == HAMILTONIAN::KohnSham_PP) {
                printf("E_PP_local_HR = %f\n", Eext_HR); 
                printf("E_PP_local    = %f (non-use)\n", Eext);
                printf("  Eext_nucl_density_HR = %f (non-use)\n", Eext_nucl_density_HR);
                printf("  Eext_nucl_density    = %f (non-use)\n", Eext_nucl_density);
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
            printf("Enn(Vlocal*rho_nucl-self_v2)HR = %f\n", E_corecore_nucl_density_HR - m_Ecore_self_HR);
            printf("Enn(Vlocal*rho_nucl-self_v2) = %f (non-use)\n", E_corecore_nucl_density - m_Ecore_self_v2);
            printf("Enn(Vlocal*rho_nucl-self) = %f (non-use)\n", E_corecore_nucl_density + m_Ecore_self);
			printf("  Enn(Vlocal*rho_nucl)HR = %f\n", E_corecore_nucl_density_HR);
            printf("  Enn(Vlocal*rho_nucl) = %f (non-use)\n", E_corecore_nucl_density);
            printf("  Enn(self_v2)HR = %f\n", m_Ecore_self_HR);
            printf("  Enn(self_v2) = %f (non-use)\n", m_Ecore_self_v2);
            printf("  Enn(self) = %f (non-use)\n", -m_Ecore_self);
			printf("  Enn(Vlocal) = %f (non-use)\n", E_corecore_Vext);
			printf("  Enn(simpleCoulomb) = %f (non-use)\n", E_corecore_direct);
			printf("correct_Coulomb_Vlocal = %f (non-use)\n", m_diff_Coulomb_Vlocal);
            printf("correct_Enn_v2 = %f\n", m_Enn_close_correction);
            fflush(stdout);

            DEBUG_PRINTF("[%d] mGetTotalEnergy 13\n", proc_id); 

			return Etot8;
		}
	}

    DEBUG_PRINTF("[%d] mGetTotalEnergy 14\n", proc_id);
	return 0.0;
	
}



//mpi supported//
inline
double QUMASUN_TD::mGetEnergyKinetic() {

	const int proc_id = GetProcessID(m_mpi_comm);
	const int local_size = ml_grid.Size3D();

	double EK = 0.0;

	for (int sk = 0; sk < m_num_having_spin_kpoint; ++sk) {

		for (int s = 0; s < num_solution; ++s) {
			const double occupancy = ml_wave_set[sk].occupancy[s];
			if (occupancy > 0.0) {
				auto& p = ml_wave_set[sk].l_psi_set[s];

				
				SoAComplex Kp{ m_work, m_work + local_size };
				SoAC::SetZero(Kp, local_size);

				DEBUG_PRINTF("[%d]mGetEnergyKinetic 1\n",proc_id); 
#if 1
				mKineticMatrixAdd_ddm_kpint(Kp, p, ml_wave_set[sk].kpoint_x, ml_wave_set[sk].kpoint_y, ml_wave_set[sk].kpoint_z);
#else
				mKineticMatrixAdd_ddm(Kp.re, p.re);
				mKineticMatrixAdd_ddm(Kp.im, p.im);
#endif
				EK += SoAC::InnerProdReal(p, Kp, local_size) * occupancy;
			}
		}
	}

	EK *= m_dx * m_dy * m_dz;
	//printf("EK=%f\n"); fflush(stdout);
	
	double tot_EK = 0.0;
	MPI_Reduce(&EK, &tot_EK, 1, MPI_DOUBLE, MPI_SUM, 0, m_mpi_comm);
	return tot_EK;
	
}

inline
double QUMASUN_TD::mGetEnergy_V_rho(const double* V, const double* rho) {

	const int64_t size_3d = ml_grid.Size3D();
	double sum2 = 0.0;
	for (int i = 0; i < size_3d; ++i) {
		sum2 += V[i] * rho[i];
	}
	sum2 *= m_dx * m_dy * m_dz;

	double sum2g = 0.0;
	MPI_Reduce(&sum2, &sum2g, 1, MPI_DOUBLE, MPI_SUM, m_root_id, ml_grid.mpi_comm);
	return sum2g;
}

inline
double QUMASUN_TD::mGetEnergy_V_rho_HR(const double* V, const double* rho) {

    const int64_t size_3d = ml_grid.Size3D() * m_HR_ratio_x * m_HR_ratio_y * m_HR_ratio_z;
    double sum2 = 0.0;
    for (int i = 0; i < size_3d; ++i) {
        sum2 += V[i] * rho[i];
    }
    sum2 *= m_dx * m_dy * m_dz / (double)(m_HR_ratio_x * m_HR_ratio_y * m_HR_ratio_z);

    double sum2g = 0.0;
    MPI_Reduce(&sum2, &sum2g, 1, MPI_DOUBLE, MPI_SUM, m_root_id, ml_grid.mpi_comm);
    return sum2g;
}


//mpi supported//
inline
double QUMASUN_TD::mGetEnergyExtByVextRho() {
	const double EV = mGetEnergy_V_rho(ml_Vext, ml_rho);
	return EV;
}

inline
double QUMASUN_TD::mGetEnergyExtByVextRho_HR() {
    const double EV = mGetEnergy_V_rho_HR(ml_hr_Vext, ml_hr_rho);
    return EV;
}

//mpi supported//
/*
* calculate
* (1/2) \int \rho V_hart(\rho) dx
*/
inline
double QUMASUN_TD::mGetEnergyHartree() {
	const double EV = mGetEnergy_V_rho(ml_Vhart, ml_rho);
	return EV / 2.0;
}

/*
* calculate
* -\int n_nucl V_hart(\rho) dx
*/
inline
double QUMASUN_TD::mGetEnergyExtByVhartNuclRho() {
	const double EV = mGetEnergy_V_rho(ml_Vhart, ml_nucl_rho);
	return EV;
}

inline
double QUMASUN_TD::mGetEnergyExtByVhartNuclRho_HR() {
    const double EV = mGetEnergy_V_rho_HR(ml_hr_Vhart, ml_hr_nucl_rho);
    return EV;
}

inline
double QUMASUN_TD::mGetEnergyCoreCoreNuclDensity() {
	const double EV = mGetEnergy_V_rho(ml_Vext, ml_nucl_rho);
	return EV / 2.0;
}


inline
double QUMASUN_TD::mGetEnergyCoreCoreNuclDensity_HR() {
    const double EV = mGetEnergy_V_rho_HR(ml_hr_Vext, ml_hr_nucl_rho);
    return EV / 2.0;
}


//mpi supported//
inline
double QUMASUN_TD::mGetEnergyXC(double* pEx, double* pEc, double* pVxRho, double* pVcRho) {
	const int64_t size_3d = ml_grid.Size3D();
	const auto& rho = ml_rho;
	const auto& pcc_rho = ml_pcc_rho;
	double sum_buf[4];
	auto& E_x = sum_buf[0];
	auto& E_c = sum_buf[1];
	auto& VxRho  = sum_buf[2];
	auto& VcRho = sum_buf[3];

	E_x = 0.0;
	E_c = 0.0;
	VxRho = 0.0;
	VcRho = 0.0;
	for (size_t i = 0; i < size_3d; ++i) {
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

	
	double sum_buf_g[4];
	MPI_Reduce(sum_buf, sum_buf_g, 4, MPI_DOUBLE, MPI_SUM, m_root_id, ml_grid.mpi_comm);
	for (int i = 0; i < 4; ++i) {
		sum_buf[i] = sum_buf_g[i];
	}
	
	if (pVxRho)*pVxRho = VxRho;
	if (pVcRho)*pVcRho = VcRho;
	*pEx = E_x;
	*pEc = E_c;
	return E_x + E_c;
}

//mpi supported//
inline
double QUMASUN_TD::mGetEnergyXCSpin(double* pEx, double* pEc, double* pVxRho_up, double* pVcRho_up, double* pVxRho_down, double* pVcRho_down) {
	const int64_t size_3d = ml_grid.Size3D();
	const auto& rho = ml_rho;
	const auto& pcc_rho = ml_pcc_rho;
	const auto& rho_diff = ml_rho_diff;
	double sum_buf[6]{ 0.0 };
	auto& E_x = sum_buf[0];
	auto& E_c = sum_buf[1];
	auto& VxRho = sum_buf[2];
	auto& VcRho = sum_buf[3];
	auto& VxRho_down = sum_buf[4];
	auto& VcRho_down = sum_buf[5];

	E_x = 0.0;
	E_c = 0.0;
	VxRho = 0.0;
	VcRho = 0.0;
	for (size_t i = 0; i < size_3d; ++i) {
		//pcc_charge should be added//
		const double rho_i = rho[i] + pcc_rho[i];
		auto ret = Calc_XC_LSDA(rho_i, rho_diff[i]/ rho_i);

		E_x += ret.E_den_x * rho_i;
		E_c += ret.E_den_c * rho_i;
		VxRho += ret.V_x_up * rho_i;
		VcRho += ret.V_c_up * rho_i;
		VxRho_down += ret.V_x_down * rho_i;
		VcRho_down += ret.V_c_down * rho_i;
	}
	E_x *= m_dx * m_dy * m_dz;
	E_c *= m_dx * m_dy * m_dz;
	VxRho *= m_dx * m_dy * m_dz;
	VcRho *= m_dx * m_dy * m_dz;
	VxRho_down *= m_dx * m_dy * m_dz;
	VcRho_down *= m_dx * m_dy * m_dz;


	double sum_buf_g[6];
	MPI_Reduce(sum_buf, sum_buf_g, 6, MPI_DOUBLE, MPI_SUM, m_root_id, ml_grid.mpi_comm);
	for (int i = 0; i < 6; ++i) {
		sum_buf[i] = sum_buf_g[i];
	}

	if (pVxRho_up)*pVxRho_up = VxRho;
	if (pVcRho_up)*pVcRho_up = VcRho;
	if (pVxRho_down)*pVxRho_down = VxRho_down;
	if (pVcRho_down)*pVcRho_down = VcRho_down;
	*pEx = E_x;
	*pEc = E_c;
	return E_x + E_c;
}

//mpi supported//
inline
double QUMASUN_TD::mGetEnergyPseudoNonlocal() {
	double EV = 0.0;
	for (int sk = 0; sk < m_num_having_spin_kpoint; ++sk) {

		for (int s = 0; s < num_solution; ++s) {
			const double occupancy = ml_wave_set[sk].occupancy[s];

			if (occupancy <= 0.0) continue;

			double pp_ene = m_pp_integrator.EnergyPPnonlocal(ml_wave_set[sk].l_psi_set[s], ml_grid, m_nuclei, m_num_nuclei, sk);

			EV += pp_ene * occupancy;
		}
		
	}

	if (IsRoot(m_ddm_comm)) {
		double tot_EV = 0.0;
		MPI_Reduce(&EV, &tot_EV, 1, MPI_DOUBLE, MPI_SUM, m_root_id, m_same_ddm_place_comm);
		return tot_EV;
	} else {
		return 0.0;
	}
	return EV;
}

inline
double QUMASUN_TD::mGetEnergyCoreCore_direct() {
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
double QUMASUN_TD::mGetEnergyPotentialAtNucl(const double* V) {
	double EV = 0.0;

	auto Folding = [](double x, double box_width) {
		return x - floor(x / box_width) * box_width;
		};

	for (int i = 0; i < m_num_nuclei; ++i) {
		const double Qi = m_nuclei_valence_elecron[i];
		double x = Folding(m_nuclei[i].Rx, m_box_x);
		double y = Folding(m_nuclei[i].Ry, m_box_y);
		double z = Folding(m_nuclei[i].Rz, m_box_z);
		

		const double val_V_ext = Interpolation<double>(x / m_dx, y / m_dy, z / m_dz, V, { 0, 0, 0, m_size_x, m_size_y, m_size_z },
			{ m_size_x, m_size_y, m_size_z });

		EV += val_V_ext * (-Qi);

	}


	return EV;
}

inline
double QUMASUN_TD::mGetEnergyCoreCoreByVextAtPoint() {
	double EV = mGetEnergyPotentialAtNucl(m_Vext.Pointer());
	return EV / 2.0;
}

inline
double QUMASUN_TD::mGetEnergyExtByVhartAtPoint() {
	double EV = mGetEnergyPotentialAtNucl(m_Vhart.Pointer());
	return EV;
}

#undef DEBUG_PRINT_ON
