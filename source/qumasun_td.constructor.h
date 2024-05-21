#pragma once
#ifdef USE_MPI
#include "mpi_helper.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <random>
#include "qumasun_td.h"
#include "vecmath.h"
#include "wrap_fft.h"
#include "vps_loader.h"
#include "pao_loader.h"
#include "cube_reader2.h"
#include "GramSchmidt_mpi.h"
#include "actual_kpoint.h"
#include "lobpcg_z_multi_mpi.h"

//#define READ_INITIAL_CUBE    "H_psi_ob1_30.cube"
//#define READ_INITIAL_PAO     1



inline
QUMASUN_TD::QUMASUN_TD(const Input& input, const MPI_Comm& mpi_comm_, const int ddm_num[4]) :
	//condition of SCF solver///////////////
	m_solver(SOLVER::LOBPCG),		//Lanczos or LOBPCG shoule be chosen//
	m_total_state(input.num_solutions),
	LIMIT_TDDFT_STEP(input.dynamics_step),
	//LOBPCG_STEP_PER_SCF(input.eigen_step_per_scf),
	//LOBPCG_STEP_FIRST(input.eigen_step_initial),
	m_hamiltonian_type(HAMILTONIAN::KohnSham_PP),    //Schrodinger, KohnSham_AE, KohnSham_PAW, KohnSham_PP//
	//condition of system size///////////////
	m_size_x(input.grid_size[0]), m_size_y(input.grid_size[1]), m_size_z(input.grid_size[2]),
	m_size_3d(m_size_x* m_size_y* m_size_z),
	m_box_x(input.box_axis[0]), m_box_y(input.box_axis[4]), m_box_z(input.box_axis[8]),
	m_dx(m_box_x / (double)m_size_x), m_dy(m_box_y / (double)m_size_y), m_dz(m_box_z / (double)m_size_z),
	m_volume(m_box_x* m_box_y* m_box_z),
    m_HR_ratio_x(input.HR_ratio), m_HR_ratio_y(input.HR_ratio), m_HR_ratio_z(input.HR_ratio),
	//kpoint///////////////////////////////////////
	m_dkx(2.0 * M_PI / (input.box_axis[0] * (double)(input.kpoint_sample[0]))),
	m_dky(2.0 * M_PI / (input.box_axis[4] * (double)(input.kpoint_sample[1]))),
	m_dkz(2.0 * M_PI / (input.box_axis[8] * (double)(input.kpoint_sample[2]))),
	//condition of electrons and nuclei///////////////
	num_electrons(input.num_electrons),
	m_num_nuclei(input.num_nuclei),
	m_pseudo_pot_set(input.pseudo_pot_set),
	m_atomic_wave_set(input.atomic_wave_set),
	m_kbT(input.temperature_K* KbHartree),
	m_mixing_ratio(input.scf_mixing_ratio),
	//optional and test////////////////////////////////////////
	is_spin_on(input.spin_polarization > 0),
    m_initial_spin_differnce(input.spin_polarization > 0 ? input.initial_spin_difference : 0),
	kpoint_symmetry(input.kpoint_symmetry),
    //initialconditions/////////////////////////
	m_initial_density(input.initial_density),
    m_initial_density_difference(input.initial_density_difference),
	m_initial_state(input.initial_state),
    m_initial_expand_from_kpoint(input.initial_expand_from_kpoint != 0),
    m_velocity_for_wave(input.velocity_for_wave),
    //mpi////////////////////////////
	m_mpi_comm(mpi_comm_)
{
//	printf("%f, %f, %f\n", m_dkx, m_dx, 1.0/m_dx);

	const int num_procs = GetNumProcess(m_mpi_comm);
	const int proc_id = GetProcessID(m_mpi_comm);
	m_num_procs_ddm = ddm_num[0] * ddm_num[1] * ddm_num[2];
    m_num_procs_sd = ddm_num[3];
	const int num_procs_spin_kpoint = num_procs / (m_num_procs_ddm * m_num_procs_sd);

	if (num_procs_spin_kpoint * m_num_procs_ddm * m_num_procs_sd != num_procs) {
		if (proc_id == m_root_id) {
			printf("The number of total MPI processes must be a multiplier of the ddm parallel number.\n"); fflush(stdout);
		}
		return;
	}

	m_kpoint_sampling[0] = input.kpoint_sample[0];
	m_kpoint_sampling[1] = input.kpoint_sample[1];
	m_kpoint_sampling[2] = input.kpoint_sample[2];
	

	std::vector<Kpoint3D> kpoint_list;
	int all_kinds_kpoint = ListupKpoints(m_kpoint_sampling[0], m_kpoint_sampling[1], m_kpoint_sampling[2], kpoint_symmetry, kpoint_list);


	//	int all_kinds_spin_kpoint = input.kpoint_sample[0] * input.kpoint_sample[1] * input.kpoint_sample[2];
	m_all_kinds_spin_kpoint = all_kinds_kpoint;
	if (is_spin_on) {
		m_all_kinds_spin_kpoint *= 2;
	}

	if (proc_id == m_root_id) {
		printf("kpoint-sample: %d, %d, %d\n", m_kpoint_sampling[0], m_kpoint_sampling[1], m_kpoint_sampling[2]); fflush(stdout);
		printf("Number of active k-point = %d\n", all_kinds_kpoint);
		for (const auto& a : kpoint_list) {
			printf("k-point(%d,%d,%d): weight=%d\n", a.kx, a.ky, a.kz, a.weight);
		}
		fflush(stdout);
		printf("all-spin-kpoint-sample: %d, %d\n", m_all_kinds_spin_kpoint, num_procs_spin_kpoint); fflush(stdout);
	}
	
	if (m_all_kinds_spin_kpoint < num_procs_spin_kpoint) {
		if (proc_id == m_root_id) {
			printf("ERROR: The number of total MPI processes is greater than the required processes.\n");
			printf("Please check the process balance of DDM and k-point parallelization.\n"); 
			fflush(stdout);
		}
		MPI_Barrier(m_mpi_comm);
		return;
	}
	
	if ((m_all_kinds_spin_kpoint / num_procs_spin_kpoint) * num_procs_spin_kpoint != m_all_kinds_spin_kpoint) {
		if (proc_id == m_root_id) {
			printf("Recommendation: it is better performance\n"
				"  when the number of total MPI processes is a multiplier of the kpoint and spin parallel number.\n"); fflush(stdout);
		}		
	}
	
	{
		
		int key_ddm_place = proc_id % m_num_procs_ddm;
        const int color_ddm = proc_id / m_num_procs_ddm;
        MPI_Comm_split(m_mpi_comm, color_ddm, key_ddm_place, &m_ddm_comm);
        //m_mpi_split_color = color_ddm;
        
        int key_state_spin_kpoint = proc_id / m_num_procs_ddm;
		const int color_place = key_ddm_place;
		MPI_Comm_split(m_mpi_comm, color_place, key_state_spin_kpoint, &m_same_ddm_place_comm);

		m_global_grid = MakeRange(m_size_x, m_size_y, m_size_z);
		ml_grid = MakeRange(m_size_x, m_size_y, m_size_z, m_ddm_comm, ddm_num);
		const size_t local_size = ml_grid.Size3D();

		
        int key_sd_rank = (proc_id / m_num_procs_ddm) % m_num_procs_sd;
        int key_spin_kpoint = proc_id / (m_num_procs_ddm * m_num_procs_sd);

		const int begin_spin_kpoint = (m_all_kinds_spin_kpoint * key_spin_kpoint) / num_procs_spin_kpoint;
		const int end_spin_kpoint = (m_all_kinds_spin_kpoint * (key_spin_kpoint+1)) / num_procs_spin_kpoint;
		//printf("[%d] spin-kpoint-range = %d, %d\n", proc_id, begin_spin_kpoint, end_spin_kpoint); fflush(stdout);

        m_begin_state = (m_total_state * key_sd_rank) / m_num_procs_sd;
        const int end_state = (m_total_state * (key_sd_rank + 1)) / m_num_procs_sd;
        num_solution = end_state - m_begin_state;

		m_having_spin_kpoint_begin = begin_spin_kpoint;
		m_num_having_spin_kpoint = end_spin_kpoint - begin_spin_kpoint;
		ml_wave_set = new WaveSet[m_num_having_spin_kpoint];
		ml_keep_lobpcg = new KeepLOBPCGSet[m_num_having_spin_kpoint];

		for (int sk = begin_spin_kpoint; sk < end_spin_kpoint; ++sk) {
			int spin = (sk >= all_kinds_kpoint) ? 1 : 0;
			const auto& a = kpoint_list[sk % all_kinds_kpoint];
			

			ml_wave_set[sk - begin_spin_kpoint].Allocate(local_size, num_solution, a.kx, a.ky, a.kz, ((spin == 0) ? SPIN::UP : SPIN::DOWN), a.weight);
			//printf("[%d] kpoint = %d,%d,%d, spin=%s\n", proc_id, a.kx, a.ky, a.kz, (spin == 0) ? "up" : "down"); fflush(stdout);

			ml_keep_lobpcg[sk - begin_spin_kpoint].Allocate(local_size, num_solution);
			vecmath::SetZero(ml_keep_lobpcg[sk - begin_spin_kpoint].keep_p, local_size * num_solution * 2);
		}

	}
	
	
	
	
	mInitializeBuffer(input);

	m_work_size = mWorkSizeTimeEvoH4(ml_grid.Size3D(), num_solution);
    m_work_size = std::max(m_work_size, ml_grid.Size3D() * 6);
    //m_work_size = std::max(m_work_size, ml_grid.Size3D() * (6+8)* num_solution);//for ForceNonlocal2
    m_work_size = std::max(m_work_size, ml_grid.Size3D() * (3 * m_HR_ratio_x * m_HR_ratio_y * m_HR_ratio_z));
	if (IsRoot(m_ddm_comm)) {
		m_work_size = std::max<size_t>(m_work_size, m_size_3d * (6 + 5 * m_HR_ratio_x * m_HR_ratio_y * m_HR_ratio_z));
	}
	m_work = new double[m_work_size];
    if (IsRoot(m_ddm_comm)) {
        printf("[%d] work_size: %zd, %p\n", proc_id, m_work_size * 8, m_work);
    }

    if (proc_id == m_root_id) {
        m_fftw = new FFTW_Executor;
        m_fftw->Initialize(m_size_x, m_size_y, m_size_z, FFTW_ESTIMATE);
        m_HR_fftw = new FFTW_Executor;
        m_HR_fftw->Initialize(m_size_x * m_HR_ratio_x, m_size_y * m_HR_ratio_y, m_size_z * m_HR_ratio_z, FFTW_ESTIMATE);
    }

	is_construction_successful = true;
};

inline
QUMASUN_TD::~QUMASUN_TD() {

    const int proc_id = GetProcessID(m_mpi_comm);

	//delete[] ml_psi[0];
	//delete[] ml_psi;
	delete[] ml_wave_set;
	delete[] ml_keep_lobpcg;
	

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


    DEBUG_PRINTF("[%d]destructor1: 0x%zx\n", proc_id, m_nuclei);

	delete[] m_nuclei;
    DEBUG_PRINTF("[%d]destructor2: 0x%zx\n", proc_id, m_nuclei_valence_elecron);
    delete[] m_nuclei_valence_elecron;
    DEBUG_PRINTF("[%d]destructor3: 0x%zx\n", proc_id, m_work);
    delete[] m_work;
    DEBUG_PRINTF("[%d]destructor4\n", proc_id);

    delete m_fftw;
    delete m_HR_fftw;
	
};

/*
* MPI並列中にcallしてよい
* process independent
*/
inline
bool QUMASUN_TD::mCheckConditions() {
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
void QUMASUN_TD::mInitializeBuffer(const Input& input) {

	//memory allocation for domain decomposed sim. on MPI//
	const size_t local_size = ml_grid.Size3D();


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
		}

        watch.Restart();

		m_pp_integrator.UpdatePosition(m_nuclei, m_num_nuclei, ml_grid);

		for (int sk = 0; sk < m_num_having_spin_kpoint; ++sk) {

			const int kpoint_x = ml_wave_set[sk].kpoint_x;
			const int kpoint_y = ml_wave_set[sk].kpoint_y;
			const int kpoint_z = ml_wave_set[sk].kpoint_z;

			const double gx_dx = (double)kpoint_x * m_dkx * m_dx;
			const double gy_dy = (double)kpoint_y * m_dky * m_dy;
			const double gz_dz = (double)kpoint_z * m_dkz * m_dz;

			m_pp_integrator.UpdateNonlocalBlock(m_nuclei, m_num_nuclei, ml_grid, sk, gx_dx, gy_dy, gz_dz);

		}
        watch.Record(20);
	}


	const bool is_root_spin = IsRoot(m_mpi_comm);
	const bool is_root_ddm = IsRoot(m_ddm_comm);

	if (is_spin_on) {

		ml_rho_diff = new double[local_size];
		vecmath::SetZero(ml_rho_diff, local_size);
		ml_Vtot_down = new double[local_size];
        ml_Vx_down = new double[local_size];
        ml_Vc_down = new double[local_size];

	}

	if (!is_root_ddm) return;//////////////////////////////////////////////

	//initialize////////////////////////////////////////////

	m_Vext.Renew(m_size_3d);
	m_rho.Renew(m_size_3d);
	m_rho_prev.Renew(m_size_3d);
	m_Vhart.Renew(m_size_3d);

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
void QUMASUN_TD::mFinalize() {

}

#endif
