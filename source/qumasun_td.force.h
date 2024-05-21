#pragma once
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include "qumasun_td.h"
#include "Vxc.h"
#include "vecmath.h"
#include "field_interpolation.h"
#include "GridGradient8th.h"
#include "GridFor.h"
#include "vec3.h"

//partially mpi supported//
/*
* If nonlocal force was already estimated in mGetTotalEnergy,
* we can skip the estimation of nonlocal force here by stand the flag "is_nonlocal_already_".
* 
*/
inline
void QUMASUN_TD::mGetForce(bool is_nonlocal_already) {
    const bool is_root_spin = IsRoot(m_mpi_comm);

    const int proc_id = GetProcessID(m_mpi_comm);
    const int num_procs = GetNumProcess(m_mpi_comm);

    m_nucl_forces.resize(m_num_nuclei);

    //nuclei-electron term from nonlocal potential//
    if(!is_nonlocal_already){
        if (m_hamiltonian_type == HAMILTONIAN::KohnSham_PP) {
            watch.Record(8);
            std::vector<vec3d> m_force_nonlocal(m_num_nuclei);
            mGetForcePseudoNonlocal((double*)&m_force_nonlocal[0]);
            //mGetForcePseudoNonlocal2((double*)&forces_n_nonlocal[0]);
                watch.Record(70);
        }
    }


	if (IsRoot(m_same_ddm_place_comm)) {
		std::vector<vec3d> forces_nn(m_num_nuclei); //nuclei-nuclei term from Hartree potential//
		std::vector<vec3d> forces_n_hart(m_num_nuclei); //nuclei-electron term from Hartree potential//
        std::vector<vec3d> forces_hr_nn(m_num_nuclei); //nuclei-nuclei term from Hartree potential//
        std::vector<vec3d> forces_hr_n_hart(m_num_nuclei); //nuclei-electron term from Hartree potential//
		std::vector<vec3d> forces_tot(m_num_nuclei); //nuclei-nuclei term from Hartree potential//
		//std::vector<double> force
		//nuclei-electron term from Hartree potential//
		mGetForceCoreCoreByVextNuclRho((double*)&forces_nn[0]);
        mGetForceHartreeNuclRho((double*)&forces_n_hart[0]);
        mGetForceCoreCoreByVextNuclRho_HR((double*)&forces_hr_nn[0]);
        mGetForceHartreeNuclRho_HR((double*)&forces_hr_n_hart[0]);

        //accept Enn correction when close distance//
        if (is_root_spin) {
            if (m_hamiltonian_type == HAMILTONIAN::KohnSham_PP) {
                for (int i = 0; i < m_num_nuclei; ++i) {
                    forces_nn[i] += m_force_nn_TF[i] - m_force_nn_correction[i];
                    forces_hr_nn[i] += m_force_nn_TF[i] - m_force_nn_correction[i];
                }
            }
        }

		for (int i = 0; i < m_num_nuclei; ++i) {
			//forces_tot[i] = forces_nn[i] + forces_n_hart[i];
		    forces_tot[i] = forces_hr_nn[i] + forces_hr_n_hart[i];
		}
		if (m_hamiltonian_type == HAMILTONIAN::KohnSham_PP) {
			for (int i = 0; i < m_num_nuclei; ++i) {
				forces_tot[i] += m_force_nonlocal[i];
			}
		}

        
        for (int i = 0; i < m_num_nuclei; ++i) {
            m_nucl_forces[i] = forces_tot[i];
        }

		if (is_root_spin) {

			vec3d total_force{ 0.0,0.0,0.0 };
			for (int ni = 0; ni < m_num_nuclei; ++ni) {
				total_force += forces_tot[ni];
			}
			printf("\nForce =========================\n");
			printf("Ftot: %8.5f  %8.5f  %8.5f\n\n", total_force.x, total_force.y, total_force.z);
			for (int ni = 0; ni < m_num_nuclei; ++ni) {
				printf("   %d: %8.5f  %8.5f  %8.5f\n", ni, forces_tot[ni].x, forces_tot[ni].y, forces_tot[ni].z);
			}

            printf("\nForce n-n (HR) with correction for close distance =====================\n");
            for (int ni = 0; ni < m_num_nuclei; ++ni) {
                printf("   %d: %8.5f  %8.5f  %8.5f\n", ni, forces_hr_nn[ni].x, forces_hr_nn[ni].y, forces_hr_nn[ni].z);
            }

            printf("\nForce n-n with correction for close distance =====================\n");
            for (int ni = 0; ni < m_num_nuclei; ++ni) {
                printf("   %d: %8.5f  %8.5f  %8.5f\n", ni, forces_nn[ni].x, forces_nn[ni].y, forces_nn[ni].z);
            }

#if 1
            printf("\nForce n-n correct =====================\n");
            for (int ni = 0; ni < m_num_nuclei; ++ni) {
                printf("   %d: %8.5f  %8.5f  %8.5f\n", ni, m_force_nn_correction[ni].x, m_force_nn_correction[ni].y, m_force_nn_correction[ni].z);
            }

            printf("\nForce n-n TF =====================\n");
            for (int ni = 0; ni < m_num_nuclei; ++ni) {
                printf("   %d: %8.5f  %8.5f  %8.5f\n", ni, m_force_nn_TF[ni].x, m_force_nn_TF[ni].y, m_force_nn_TF[ni].z);
            }
#endif

            printf("\nForce n-hart (HR) =====================\n");
            for (int ni = 0; ni < m_num_nuclei; ++ni) {
                printf("   %d: %8.5f  %8.5f  %8.5f\n", ni, forces_hr_n_hart[ni].x, forces_hr_n_hart[ni].y, forces_hr_n_hart[ni].z);
            }

			printf("\nForce n-hart =====================\n");
			for (int ni = 0; ni < m_num_nuclei; ++ni) {
				printf("   %d: %8.5f  %8.5f  %8.5f\n", ni, forces_n_hart[ni].x, forces_n_hart[ni].y, forces_n_hart[ni].z);
			}

			if (m_hamiltonian_type == HAMILTONIAN::KohnSham_PP) {
				printf("\nForce n-nonlocal =====================\n");
				for (int ni = 0; ni < m_num_nuclei; ++ni) {
					printf("   %d: %8.5f  %8.5f  %8.5f\n", ni, m_force_nonlocal[ni].x, m_force_nonlocal[ni].y, m_force_nonlocal[ni].z);
				}
			}

			fflush(stdout);

			
		}
	}

	//must do Bcast forces////
    MPI_Bcast(&m_nucl_forces[0], m_num_nuclei * 3, MPI_DOUBLE, 0, m_mpi_comm);
}


/*
* calculate froce from
* -\int n_nucl V_any(\rho) dx
*/
inline
void QUMASUN_TD::mGetForceOnNuclRho(double* force, const double* l_Vpot) {
	size_t local_size = ml_grid.Size3D();
	double* l_dV_dx = m_work;
	double* l_dV_dy = m_work + local_size;
	double* l_dV_dz = m_work + local_size*2;

	Gradient8th_ddm(ml_grid, l_dV_dx, l_dV_dy, l_dV_dz, l_Vpot, -1.0, m_dx, m_dy, m_dz);
    watch.Record(71);

	const int stride = 3;
	m_pp_integrator.InnerForChargeVlocal(force, stride, l_dV_dx, ml_grid, m_nuclei, m_num_nuclei);
	m_pp_integrator.InnerForChargeVlocal(force+1, stride, l_dV_dy, ml_grid, m_nuclei, m_num_nuclei);
	m_pp_integrator.InnerForChargeVlocal(force+2, stride, l_dV_dz, ml_grid, m_nuclei, m_num_nuclei);
    watch.Record(72);
}

/*
* calculate froce from
* -\int n_nucl V_any(\rho) dx
*/
inline
void QUMASUN_TD::mGetForceOnNuclRho_HR(double* force, const double* l_hr_Vpot) {
    size_t local_size = ml_grid.Size3D() * m_HR_ratio_x * m_HR_ratio_y * m_HR_ratio_z;
    double* l_dV_dx = m_work;
    double* l_dV_dy = m_work + local_size;
    double* l_dV_dz = m_work + local_size * 2;
    auto hr_grid = ScaleRange(ml_grid, m_HR_ratio_x, m_HR_ratio_y, m_HR_ratio_z);
    DEBUG_PRINTF("[%d]before Gradient8th_ddm: %zd\n", GetProcessID(ml_grid.mpi_comm),local_size);
    Gradient8th_ddm(hr_grid, l_dV_dx, l_dV_dy, l_dV_dz, l_hr_Vpot, -1.0, m_dx/(double)m_HR_ratio_x, m_dy / (double)m_HR_ratio_y, m_dz / (double)m_HR_ratio_z);
    DEBUG_PRINTF("after Gradient8th_ddm\n");
    watch.Record(73);
    const int stride = 3;
    m_pp_integrator.InnerForChargeVlocal(force, stride, l_dV_dx, ml_grid, m_nuclei, m_num_nuclei, true);
    DEBUG_PRINTF("InnerForChargeVlocal x\n");
    m_pp_integrator.InnerForChargeVlocal(force + 1, stride, l_dV_dy, ml_grid, m_nuclei, m_num_nuclei, true);
    DEBUG_PRINTF("InnerForChargeVlocal y\n"); 
    m_pp_integrator.InnerForChargeVlocal(force + 2, stride, l_dV_dz, ml_grid, m_nuclei, m_num_nuclei, true);
    DEBUG_PRINTF("InnerForChargeVlocal z\n");
    watch.Record(74);
}

/*
* calculate froce from
* -\int n_nucl V_nucl(\rho) dx
*/
inline
void QUMASUN_TD::mGetForceCoreCoreByVextNuclRho(double* force) {
    mGetForceOnNuclRho(force, ml_Vext);
#if 0
	size_t local_size = ml_grid.Size3D();
	double* l_dV_dx = m_work;
	double* l_dV_dy = m_work + local_size;
	double* l_dV_dz = m_work + local_size*2;

	Gradient8th_ddm(ml_grid, l_dV_dx, l_dV_dy, l_dV_dz, ml_Vext, -1.0, m_dx, m_dy, m_dz);
	
	const int stride = 3;
	m_pp_integrator.InnerForChargeVlocal(force, stride, l_dV_dx, ml_grid, m_nuclei, m_num_nuclei);
	m_pp_integrator.InnerForChargeVlocal(force+1, stride, l_dV_dy, ml_grid, m_nuclei, m_num_nuclei);
	m_pp_integrator.InnerForChargeVlocal(force+2, stride, l_dV_dz, ml_grid, m_nuclei, m_num_nuclei);
#endif
}

/*
* calculate froce from
* -\int n_nucl V_nucl(\rho) dx
*/
inline
void QUMASUN_TD::mGetForceCoreCoreByVextNuclRho_HR(double* force) {
    mGetForceOnNuclRho_HR(force, ml_hr_Vext);
}

/*
* calculate froce from
* -\int n_nucl V_nucl(\rho) dx
*/
inline
void QUMASUN_TD::mGetForceHartreeNuclRho(double* force) {
	mGetForceOnNuclRho(force, ml_Vhart);
#if 0	
	size_t local_size = ml_grid.Size3D();
	double* l_dV_dx = m_work;
	double* l_dV_dy = m_work + local_size;
	double* l_dV_dz = m_work + local_size * 2;

	Gradient8th_ddm(ml_grid, l_dV_dx, l_dV_dy, l_dV_dz, ml_Vhart, -1.0, m_dx, m_dy, m_dz);

	const int stride = 3;
	m_pp_integrator.InnerForChargeVlocal(force, stride, l_dV_dx, ml_grid, m_nuclei, m_num_nuclei);
	m_pp_integrator.InnerForChargeVlocal(force + 1, stride, l_dV_dy, ml_grid, m_nuclei, m_num_nuclei);
	m_pp_integrator.InnerForChargeVlocal(force + 2, stride, l_dV_dz, ml_grid, m_nuclei, m_num_nuclei);
#endif
}

inline
void QUMASUN_TD::mGetForceHartreeNuclRho_HR(double* force) {
    mGetForceOnNuclRho_HR(force, ml_hr_Vhart);
}

//mpi supported//
inline
double QUMASUN_TD::mGetForcePseudoNonlocal(double* force) {
	
	size_t local_size = ml_grid.Size3D();

	for (int i = 0; i < m_num_nuclei*3; ++i) {
		force[i] = 0;
	}
	SoAComplex l_dpsi_dx{ m_work, m_work + local_size };
	SoAComplex l_dpsi_dy{ m_work + local_size * 2, m_work + local_size * 3 };
	SoAComplex l_dpsi_dz{ m_work + local_size * 4, m_work + local_size * 5 };

	double EV = 0.0;
	for (int sk = 0; sk < m_num_having_spin_kpoint; ++sk) {
		const double gx = (double)(ml_wave_set[sk].kpoint_x) * m_dkx;
		const double gy = (double)(ml_wave_set[sk].kpoint_y) * m_dky;
		const double gz = (double)(ml_wave_set[sk].kpoint_z) * m_dkz;

		for (int s = 0; s < num_solution; ++s) {
			const double occupancy = ml_wave_set[sk].occupancy[s];

			if (occupancy <= 0.0) continue;


			const double coef = -2.0 * occupancy;
			Gradient8th_ddm(ml_grid, l_dpsi_dx.re, l_dpsi_dy.re, l_dpsi_dz.re, ml_wave_set[sk].l_psi_set[s].re, coef, m_dx, m_dy, m_dz);
			Gradient8th_ddm(ml_grid, l_dpsi_dx.im, l_dpsi_dy.im, l_dpsi_dz.im, ml_wave_set[sk].l_psi_set[s].im, coef, m_dx, m_dy, m_dz);
            watch.Record(75);
			double pp_ene = m_pp_integrator.ForceNonlocal(force, ml_wave_set[sk].l_psi_set[s], l_dpsi_dx, l_dpsi_dy, l_dpsi_dz, ml_grid, m_nuclei, m_num_nuclei, sk, gx* coef, gy* coef, gz* coef);
            watch.Record(76);
			EV += pp_ene * occupancy;
		}

	}



	if (IsRoot(m_ddm_comm)) {
		double tot_EV = 0.0;
		MPI_Reduce(&EV, &tot_EV, 1, MPI_DOUBLE, MPI_SUM, m_root_id, m_same_ddm_place_comm);


		std::vector<double> force_sum(m_num_nuclei * 3);
		MPI_Reduce(force, &force_sum[0], m_num_nuclei * 3, MPI_DOUBLE, MPI_SUM, m_root_id, m_same_ddm_place_comm);

		if (IsRoot(m_same_ddm_place_comm)) {
			for (int i = 0; i < m_num_nuclei * 3; ++i) {
				force[i] = force_sum[i];
			}
		}

		return tot_EV;
	} else {
		return 0.0;
	}
	
}

#if 0

//slower than original ForceNonlocal that uses just one psi and its diferentials //
//mpi supported//
inline
double QUMASUN_TD::mGetForcePseudoNonlocal2(double* force) {

    size_t local_size = ml_grid.Size3D();

    for (int i = 0; i < m_num_nuclei * 3; ++i) {
        force[i] = 0;
    }

    std::vector<double> valid_occupancy(num_solution);
    std::vector<SoAComplex> valid_psi(num_solution);

    double EV = 0.0;
    for (int sk = 0; sk < m_num_having_spin_kpoint; ++sk) {
        const double gx = (double)(ml_wave_set[sk].kpoint_x) * m_dkx;
        const double gy = (double)(ml_wave_set[sk].kpoint_y) * m_dky;
        const double gz = (double)(ml_wave_set[sk].kpoint_z) * m_dkz;

        auto* p_buffer = m_work;

        int num_valid = 0;
        for (int s = 0; s < num_solution; ++s) {
            const double occupancy = ml_wave_set[sk].occupancy[s];

            if (occupancy <= 0.0) continue;
            valid_occupancy[num_valid] = occupancy;
            valid_psi[num_valid].re = ml_wave_set[sk].l_psi_set[s].re;
            valid_psi[num_valid].im = ml_wave_set[sk].l_psi_set[s].im;

            SoAComplex l_dpsi_dx{ p_buffer, p_buffer + local_size };
            SoAComplex l_dpsi_dy{ p_buffer + local_size * 2, p_buffer + local_size * 3 };
            SoAComplex l_dpsi_dz{ p_buffer + local_size * 4, p_buffer + local_size * 5 };
            p_buffer += local_size * 6;
            ++num_valid;

            const double coef = -2.0 * occupancy;
            Gradient8th_ddm(ml_grid, l_dpsi_dx.re, l_dpsi_dy.re, l_dpsi_dz.re, ml_wave_set[sk].l_psi_set[s].re, coef, m_dx, m_dy, m_dz);
            Gradient8th_ddm(ml_grid, l_dpsi_dx.im, l_dpsi_dy.im, l_dpsi_dz.im, ml_wave_set[sk].l_psi_set[s].im, coef, m_dx, m_dy, m_dz);

        }

        //ここでp_buffer以降はまだworkとして使える
        watch.Record(75);

        {
            //if (occupancy <= 0.0) continue;
            
            const double coef = -2.0;
            double pp_ene = m_pp_integrator.ForceNonlocal2(force, num_valid, &valid_psi[0], m_work, p_buffer, &valid_occupancy[0], ml_grid, m_nuclei, m_num_nuclei, sk, gx * coef, gy * coef, gz * coef);
            watch.Record(76);
            EV += pp_ene ;
        }
    }



    if (IsRoot(m_ddm_comm)) {
        double tot_EV = 0.0;
        MPI_Reduce(&EV, &tot_EV, 1, MPI_DOUBLE, MPI_SUM, m_root_id, m_same_ddm_place_comm);


        std::vector<double> force_sum(m_num_nuclei * 3);
        MPI_Reduce(force, &force_sum[0], m_num_nuclei * 3, MPI_DOUBLE, MPI_SUM, m_root_id, m_same_ddm_place_comm);

        if (IsRoot(m_same_ddm_place_comm)) {
            for (int i = 0; i < m_num_nuclei * 3; ++i) {
                force[i] = force_sum[i];
            }
        }

        return tot_EV;
    } else {
        return 0.0;
    }

}
#endif

void QUMASUN_TD::GetForce(vec3d* forces) {
    for (int i = 0; i < m_num_nuclei; ++i) {
        forces[i] = m_nucl_forces[i];
    }
}
