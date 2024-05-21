#pragma once
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include "qumasun_mpi.h"
#include "Vxc.h"
#include "vecmath.h"
#include "field_interpolation.h"
#include "GridGradient8th.h"
#include "GridFor.h"
#include "vec3.h"

//partially mpi supported//
inline
void QUMASUN_MPI::mGetForce() {
	const bool is_root_spin = IsRoot(m_mpi_comm);

	const int proc_id = GetProcessID(m_mpi_comm);
	const int num_procs = GetNumProcess(m_mpi_comm);

	std::vector<vec3d> forces_n_nonlocal(m_num_nuclei); //nuclei-electron term from nonlocal potential//
	if (m_hamiltonian_type == HAMILTONIAN::KohnSham_PP) {
		mGetForcePseudoNonlocal((double*)&forces_n_nonlocal[0]);
	}

	if (has_up_spin) {
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
				forces_tot[i] += forces_n_nonlocal[i];
			}
		}

        m_nucl_forces.resize(m_num_nuclei);
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

#if 0
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
					printf("   %d: %8.5f  %8.5f  %8.5f\n", ni, forces_n_nonlocal[ni].x, forces_n_nonlocal[ni].y, forces_n_nonlocal[ni].z);
				}
			}

			fflush(stdout);

			
		}
	}

	//must do Bcast forces////
}


/*
* calculate froce from
* -\int n_nucl V_any(\rho) dx
*/
inline
void QUMASUN_MPI::mGetForceOnNuclRho(double* force, const double* l_Vpot) {
	size_t local_size = ml_grid.Size3D();
	double* l_dV_dx = m_work;
	double* l_dV_dy = m_work + local_size;
	double* l_dV_dz = m_work + local_size*2;

	Gradient8th_ddm(ml_grid, l_dV_dx, l_dV_dy, l_dV_dz, l_Vpot, -1.0, m_dx, m_dy, m_dz);
	
	const int stride = 3;
	m_pp_integrator.InnerForChargeVlocal(force, stride, l_dV_dx, ml_grid, m_nuclei, m_num_nuclei);
	m_pp_integrator.InnerForChargeVlocal(force+1, stride, l_dV_dy, ml_grid, m_nuclei, m_num_nuclei);
	m_pp_integrator.InnerForChargeVlocal(force+2, stride, l_dV_dz, ml_grid, m_nuclei, m_num_nuclei);
	
}

/*
* calculate froce from
* -\int n_nucl V_any(\rho) dx
*/
inline
void QUMASUN_MPI::mGetForceOnNuclRho_HR(double* force, const double* l_hr_Vpot) {
    size_t local_size = ml_grid.Size3D() * m_HR_ratio_x * m_HR_ratio_y * m_HR_ratio_z;
    double* l_dV_dx = m_work;
    double* l_dV_dy = m_work + local_size;
    double* l_dV_dz = m_work + local_size * 2;
    auto hr_grid = ScaleRange(ml_grid, m_HR_ratio_x, m_HR_ratio_y, m_HR_ratio_z);
    
    Gradient8th_ddm(hr_grid, l_dV_dx, l_dV_dy, l_dV_dz, l_hr_Vpot, -1.0, m_dx/(double)m_HR_ratio_x, m_dy / (double)m_HR_ratio_y, m_dz / (double)m_HR_ratio_z);
    
    const int stride = 3;
    m_pp_integrator.InnerForChargeVlocal(force, stride, l_dV_dx, ml_grid, m_nuclei, m_num_nuclei, true);
    m_pp_integrator.InnerForChargeVlocal(force + 1, stride, l_dV_dy, ml_grid, m_nuclei, m_num_nuclei, true);
    m_pp_integrator.InnerForChargeVlocal(force + 2, stride, l_dV_dz, ml_grid, m_nuclei, m_num_nuclei, true);

}

/*
* calculate froce from
* -\int n_nucl V_nucl(\rho) dx
*/
inline
void QUMASUN_MPI::mGetForceCoreCoreByVextNuclRho(double* force) {
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
void QUMASUN_MPI::mGetForceCoreCoreByVextNuclRho_HR(double* force) {
    mGetForceOnNuclRho_HR(force, ml_hr_Vext);
}

/*
* calculate froce from
* -\int n_nucl V_nucl(\rho) dx
*/
inline
void QUMASUN_MPI::mGetForceHartreeNuclRho(double* force) {
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
void QUMASUN_MPI::mGetForceHartreeNuclRho_HR(double* force) {
    mGetForceOnNuclRho_HR(force, ml_hr_Vhart);
}

//mpi supported//
inline
double QUMASUN_MPI::mGetForcePseudoNonlocal(double* force) {
	double EV = 0.0;

	size_t local_size = ml_grid.Size3D();

	for (int i = 0; i < m_num_nuclei*3; ++i) {
		force[i] = 0;
	}

	double* l_dpsi_dx = m_work;
	double* l_dpsi_dy = m_work + local_size;
	double* l_dpsi_dz = m_work + local_size * 2;
	{

		for (int s = 0; s < num_solution; ++s) {

			if (m_occupancy[s] <= 0.0) continue;

			const double coef = -2.0 * m_occupancy[s];
			Gradient8th_ddm(ml_grid, l_dpsi_dx, l_dpsi_dy, l_dpsi_dz, ml_psi[s], coef, m_dx, m_dy, m_dz);

			double pp_ene = m_pp_integrator.ForceNonlocal(force, ml_psi[s], l_dpsi_dx, ml_grid, m_nuclei, m_num_nuclei);

			EV += pp_ene * m_occupancy[s];
		}

	}
	if (m_num_spin == 2) {

		for (int s = 0; s < num_solution; ++s) {

			if (m_occupancy_down[s] <= 0.0) continue;

			const double coef = -2.0 * m_occupancy_down[s];
			Gradient8th_ddm(ml_grid, l_dpsi_dx, l_dpsi_dy, l_dpsi_dz, ml_psi[s+ num_solution], coef, m_dx, m_dy, m_dz);

			double pp_ene = m_pp_integrator.ForceNonlocal(force, ml_psi[s+ num_solution], l_dpsi_dx, ml_grid, m_nuclei, m_num_nuclei);

			
			EV += pp_ene * m_occupancy_down[s];
		}

	} else {
		if (is_spin_on) {
			if (IsRoot(m_ddm_comm)) {
				int proc_id = GetProcessID(m_mpi_comm);
				int num_procs = GetNumProcess(m_mpi_comm);
				const int TAG = 0x123;

				if (proc_id == 0) {
					double E_from_down;
					MPI_Status status;
					MPI_Recv(&E_from_down, 1, MPI_DOUBLE, num_procs / 2, TAG, m_mpi_comm, &status);
					EV += E_from_down;


					std::vector<double> force_down(m_num_nuclei * 3);
					MPI_Recv(&force_down[0], m_num_nuclei * 3, MPI_DOUBLE, num_procs / 2, TAG, m_mpi_comm, &status);

					for (int i = 0; i < m_num_nuclei * 3; ++i) {
						force[i] += force_down[i];
					}
				} else {
					MPI_Send(&EV, 1, MPI_DOUBLE, 0, TAG, m_mpi_comm);
					MPI_Send(&force[0], m_num_nuclei * 3, MPI_DOUBLE, 0, TAG, m_mpi_comm);
				}
			}
		}
	}
	return EV;
}


void QUMASUN_MPI::GetForce(vec3d* forces) {
    for (int i = 0; i < m_num_nuclei; ++i) {
        forces[i] = m_nucl_forces[i];
    }
}

