#pragma once
#ifdef USE_MPI
#include "mpi_helper.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <random>
#include "qumasun_mpi.h"
#include "vecmath.h"
#include "wrap_fft.h"
#include "vps_loader.h"
#include "pao_loader.h"
#include "cube_reader2.h"
#include "GramSchmidt_mpi.h"
#include "lobpcg_d_multi_mpi.h"

//#define READ_INITIAL_CUBE    "H_psi_ob1_30.cube"
//#define READ_INITIAL_PAO     1



inline
QUMASUN_MPI::QUMASUN_MPI(const Input& input, const MPI_Comm& mpi_comm_, const int ddm_num[3]) :
	//condition of SCF solver///////////////
	m_solver(SOLVER::LOBPCG),		//Lanczos or LOBPCG shoule be chosen//
	num_solution(input.num_solutions),
	LIMIT_SCF_STEP(input.scf_step),
	LOBPCG_STEP_PER_SCF(input.eigen_step_per_scf),
	LOBPCG_STEP_FIRST(input.eigen_step_initial),
	m_hamiltonian_type(HAMILTONIAN::KohnSham_PP),    //Schrodinger, KohnSham_AE, KohnSham_PAW, KohnSham_PP//
	//condition of system size///////////////
	m_size_x(input.grid_size[0]), m_size_y(input.grid_size[1]), m_size_z(input.grid_size[2]),
	m_size_3d(m_size_x*m_size_y*m_size_z),
	m_box_x(input.box_axis[0]), m_box_y(input.box_axis[4]), m_box_z(input.box_axis[8]),
	m_dx(m_box_x/(double)m_size_x), m_dy(m_box_y / (double)m_size_y), m_dz(m_box_z / (double)m_size_z),
	m_volume(m_box_x*m_box_y*m_box_z),
    m_HR_ratio_x(input.HR_ratio), m_HR_ratio_y(input.HR_ratio), m_HR_ratio_z(input.HR_ratio),
	//condition of electrons and nuclei///////////////
	num_electrons(input.num_electrons),
	m_num_nuclei(input.num_nuclei),
	m_pseudo_pot_set(input.pseudo_pot_set),
	m_atomic_wave_set(input.atomic_wave_set),
	m_initial_density(input.initial_density), 
	m_kbT(input.temperature_K* KbHartree),
	m_mixing_ratio(input.scf_mixing_ratio),
	//optional and test////////////////////////////////////////
	is_spin_on(input.spin_polarization),

	//mpi////////////////////////////
	m_mpi_comm(mpi_comm_)
{

	const int num_procs = GetNumProcess(m_mpi_comm);
	const int proc_id = GetProcessID(m_mpi_comm);

	const int num_procs_ddm = ddm_num[0] * ddm_num[1] * ddm_num[2];



	if (is_spin_on) {
		if (num_procs_ddm * 2 == num_procs) {
			m_num_spin = 1;

			int key = proc_id % (num_procs / 2);
			int color = proc_id/ (num_procs / 2);
			MPI_Comm_split(m_mpi_comm, color, key, &m_ddm_comm);
			m_mpi_split_color = color;

			if (color == 0) {
				has_up_spin = true;
			} else {
				has_up_spin = false;
			}
			//printf("test0-1: %d / %d\n", ddm_num[0] * ddm_num[1] * ddm_num[2], num_procs); fflush(stdout);
		} else if (num_procs == num_procs_ddm) {
		
			m_num_spin = 2;
			m_ddm_comm = m_mpi_comm;
			has_up_spin = true;
			m_mpi_split_color = 0;
			//printf("test0-2: %d / %d\n", ddm_num[0] * ddm_num[1] * ddm_num[2], num_procs); fflush(stdout);
		}else{
			if (proc_id == m_root_id) {
				printf("ERROR: The number of total MPI processes must be a same as the ddm parallel number.\n"); fflush(stdout);
			}
			return;
		}
	} else {

		if (num_procs != num_procs_ddm) {
			if (proc_id == m_root_id) {
				printf("ERROR: The number of total MPI processes must be a same as the ddm parallel number.\n"); fflush(stdout);
			}
			return;
		}


		m_num_spin = 1;
		m_ddm_comm = m_mpi_comm;
		has_up_spin = true;
		m_mpi_split_color = 0;
	}

	m_global_grid = MakeRange(m_size_x, m_size_y, m_size_z);
	ml_grid = MakeRange(m_size_x, m_size_y, m_size_z, m_ddm_comm, ddm_num);

	mInitializeBuffer(input);


	m_work_size = LOBPCG::WorkSize_d_multi_mpi(ml_grid.Size3D(), num_solution);
    m_work_size = std::max<size_t>(m_work_size, ml_grid.Size3D() * (3 * m_HR_ratio_x * m_HR_ratio_y * m_HR_ratio_z));
	if (IsRoot(m_ddm_comm)) {
		m_work_size = std::max<size_t>(m_work_size, m_size_3d * (6 + 5 * m_HR_ratio_x * m_HR_ratio_y * m_HR_ratio_z));
	}
	m_work = new double[m_work_size];

	is_construction_successful = true;
};

inline
QUMASUN_MPI::~QUMASUN_MPI() {
	delete[] ml_psi[0];
	delete[] ml_psi;
	delete[] ml_lobpcg_keep_p;
	delete[] m_lobpcg_keep_S_matrix;

	delete[] ml_rho;
	delete[] ml_Vtot;
	delete[] ml_Vext;
	delete[] ml_Vhart;
	delete[] ml_Vtot_down;
	delete[] ml_Vx;
	delete[] ml_Vc;
    delete[] ml_Vx_down;
    delete[] ml_Vc_down;
	delete[] ml_pcc_rho;
	delete[] ml_nucl_rho;
	delete[] ml_rho_diff;

    delete[] ml_hr_nucl_rho;
    delete[] ml_hr_Vext;
    delete[] ml_hr_rho;
    delete[] ml_hr_Vhart;

	delete[] m_eigen_values;
	delete[] m_eigen_values_down;
	delete[] m_occupancy;
	delete[] m_occupancy_down;

	delete[] m_nuclei;
	delete[] m_nuclei_valence_elecron;
	delete[] m_work;

	
};

/*
* MPI並列中にcallしてよい
* process independent
*/
inline
bool QUMASUN_MPI::mCheckConditions() {
	if (IsRoot(m_mpi_comm)) {
		printf("sizeof(fftw_complex) = %d\n", (int)sizeof(fftw_complex));
		printf("sizeof(std::complex<double>) = %d\n", (int)sizeof(std::complex<double>));
		fflush(stdout);
	}
	if (sizeof(fftw_complex) != sizeof(std::complex<double>)) {
		if (IsRoot(m_mpi_comm)) {
			printf("Error: different size, sizeof(fftw_complex) != sizeof(std::complex<double>)\n");
		}
		return false;
	}

	return true;
}

inline
void QUMASUN_MPI::mInitializeBuffer(const Input& input) {

	//memory allocation for domain decomposed sim. on MPI//
	const size_t local_size = ml_grid.Size3D();
	ml_psi = new Field[num_solution * m_num_spin];
	ml_psi[0] = new double[num_solution * m_num_spin * local_size];
	for (size_t s = 1; s < num_solution * m_num_spin; ++s) {
		ml_psi[s] = ml_psi[0] + s * local_size;
	}

	ml_lobpcg_keep_p = new double[num_solution * m_num_spin * local_size];
	vecmath::SetZero(ml_lobpcg_keep_p, num_solution * m_num_spin * local_size);
	m_lobpcg_keep_S_matrix = new double[num_solution * num_solution * 3 * m_num_spin];

	ml_rho = new double[local_size];
	vecmath::SetZero(ml_rho, local_size);
	ml_Vtot = new double[local_size];
	ml_Vext = new double[local_size];
	ml_Vhart = new double[local_size];
	ml_Vx = new double[local_size];
	ml_Vc = new double[local_size];
	ml_pcc_rho = new double[local_size];
	ml_nucl_rho = new double[local_size];
    
    const size_t hr_ratio = m_HR_ratio_x * m_HR_ratio_y * m_HR_ratio_z;
    ml_hr_nucl_rho = new double[local_size * hr_ratio];
    ml_hr_Vext = new double[local_size * hr_ratio];
    ml_hr_rho = new double[local_size * hr_ratio];
    ml_hr_Vhart = new double[local_size * hr_ratio];

	//ml_work_rd = new double[local_size];
	//ml_work_rd2 = new double[local_size];

	m_occupancy = new double[num_solution];
	m_eigen_values = new double[num_solution];


	m_nuclei = new Nucleus[m_num_nuclei];
	for (int i = 0; i < m_num_nuclei; ++i) {
		m_nuclei[i] = input.nuclei[i];
	}

	m_nuclei_valence_elecron = new double[m_num_nuclei]; 


	if (m_hamiltonian_type == HAMILTONIAN::KohnSham_PP) {
		
		m_pp_integrator.InitializeGrid(GridInfo{ (size_t)m_size_3d, m_size_x, m_size_y, m_size_z }, m_dx, m_dy, m_dz);
        m_pp_integrator.SetHighResolution(m_HR_ratio_x, m_HR_ratio_y, m_HR_ratio_z);

		for (const auto& p : m_pseudo_pot_set) {
			m_pp_integrator.Load(p.first, p.second.c_str(), m_mpi_comm);
			//m_pp_integrator.Load(p.first, p.second.c_str());
		}

        watch.Restart();
		m_pp_integrator.UpdatePosition(m_nuclei, m_num_nuclei, ml_grid);
        watch.Record(20);

	}


	const bool is_root_spin = IsRoot(m_mpi_comm);
	const bool is_root_ddm = IsRoot(m_ddm_comm);

	if (is_spin_on) {

		ml_rho_diff = new double[local_size];

		if ((m_num_spin == 2) || (is_root_spin)) {
			//spin is not parallelized//
			m_occupancy_down = new double[num_solution];
			m_eigen_values_down = new double[num_solution];
		}

		ml_Vtot_down = new double[local_size];
        ml_Vx_down = new double[local_size];
        ml_Vc_down = new double[local_size];
	}

	if (!is_root_ddm) return;//////////////////////////////////////////////

	//initialize////////////////////////////////////////////

	//m_Vtot.Renew(m_size_3d);
	m_Vext.Renew(m_size_3d);
	m_rho.Renew(m_size_3d);
	m_rho_prev.Renew(m_size_3d);
	m_Vhart.Renew(m_size_3d);
	//m_Vx.Renew(m_size_3d);
	//m_Vc.Renew(m_size_3d);

    m_hr_Vext.Renew(m_size_3d* hr_ratio);
    m_hr_rho.Renew(m_size_3d * hr_ratio);
    m_hr_Vhart.Renew(m_size_3d * hr_ratio);

	//for spin polarization///////////////////
	if (is_spin_on) {
		m_rho_diff.Renew(m_size_3d);
		

	}



	//load pseudo_potential//
	if (m_hamiltonian_type == HAMILTONIAN::KohnSham_PP) {
		
		m_pcc_rho.Renew(m_size_3d);


	}
}

inline
void QUMASUN_MPI::mFinalize() {

}


#endif
