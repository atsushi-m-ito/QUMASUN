#pragma once
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <random>
#include "qumasun_single.h"
#include "vecmath.h"
#include "wrap_fft.h"
#include "vps_loader.h"
#include "pao_loader.h"
#include "cube_reader2.h"
#include "GramSchmidt.h"
#include "lobpcg_d_multi.h"

//#define READ_INITIAL_CUBE    "H_psi_ob1_30.cube"
//#define READ_INITIAL_PAO     1



inline
QUMASUN_SINGLE::QUMASUN_SINGLE(const Input& input) :
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
	//condition of electrons and nuclei///////////////
	num_electrons(input.num_electrons),
	m_num_nuclei(input.num_nuclei),
	m_pseudo_pot_set(input.pseudo_pot_set),
	m_atomic_wave_set(input.atomic_wave_set),
	m_initial_density(input.initial_density),
	//optional and test////////////////////////////////////////
	is_spin_on(input.spin_polarization==1),
	m_num_spin(is_spin_on ? 2 : 1)
{

	m_nuclei = new Nucleus[m_num_nuclei];
	for (int i = 0; i < m_num_nuclei; ++i) {
		m_nuclei[i] = input.nuclei[i];
	}

	m_nuclei_valence_elecron = new double[m_num_nuclei];

	//initialize////////////////////////////////////////////

	m_Vtot.Renew(m_size_3d);
	m_Vext.Renew(m_size_3d);
	m_rho.Renew(m_size_3d);
	m_rho_prev.Renew(m_size_3d);
	m_Vhart.Renew(m_size_3d);
	m_Vx.Renew(m_size_3d);
	m_Vc.Renew(m_size_3d);
	m_psi_set = new RspaceFunc<double>[num_solution * m_num_spin];
	m_psi_buffer.Renew(m_size_3d * num_solution* m_num_spin);
	for (size_t s = 0; s < num_solution* m_num_spin; ++s) {
		m_psi_set[s].Refer(m_psi_buffer.Pointer() + m_size_3d * s);
	}
	m_lobpcg_keep_p = new double[m_size_3d * num_solution * m_num_spin];
	m_occupancy = new double[num_solution];
	m_eigen_values = new double[num_solution];

	//for spin polarization///////////////////
	if (is_spin_on) {
		m_occupancy_down = new double[num_solution];
		m_eigen_values_down = new double[num_solution];
		m_Vtot_down.Renew(m_size_3d);
		m_Vx_down.Renew(m_size_3d);
		m_Vc_down.Renew(m_size_3d);
		m_rho_diff.Renew(m_size_3d);
	}
	
	
	
	//load pseudo_potential//
	if (m_hamiltonian_type == HAMILTONIAN::KohnSham_PP) {
		m_pp_integrator.InitializeGrid(GridInfo{(size_t)m_size_3d, m_size_x, m_size_y, m_size_z}, m_dx, m_dy, m_dz);
		
		for (const auto& p : m_pseudo_pot_set) {
			m_pp_integrator.Load(p.first, p.second.c_str());
		}
		m_pp_integrator.UpdatePosition(m_nuclei, m_num_nuclei);
		m_pcc_rho.Renew(m_size_3d);
		m_nucl_rho.Renew(m_size_3d);

		/* test output//
		{
			FILE* fp = fopen("read.PP.txt", "w");
			fprintf(fp, "#xi, r, V_local, p1, p2, ,,,\n");
			const auto* pp = m_pp_integrator.mFindPseudoPot(m_nuclei[0].Z);
			const int OUT_N = 500;
			const double xi_min = -7.8;
			const double xi_max = 2.177230128181548; //from PAO//
			const double xi_delta = (xi_max - xi_min) / (double)(OUT_N - 1);
			for (int i = 0; i < OUT_N; ++i) {
				const double xi = xi_delta * (double)i + xi_min;
				const double r = exp(xi);

				const double vloc = RadialGrid2::GetValue(r, pp->V_local, pp->num_radial_grids, pp->xi_min, pp->xi_delta);
				const double nonloc_s1 = RadialGrid2::GetValue(r, pp->projector[0], pp->num_radial_grids, pp->xi_min, pp->xi_delta);
				const double nonloc_s2 = RadialGrid2::GetValue(r, pp->projector[2], pp->num_radial_grids, pp->xi_min, pp->xi_delta);
				fprintf(fp, "%.12f\t%.12f\t%.12f\t%.12f\t%.12f\n", xi,r,vloc, nonloc_s1, nonloc_s2);
			}
			fclose(fp);
		}
		*/
	}

	size_t work_size = WorkSizeLOBPCG_d_multi(m_size_3d, num_solution);
	work_size = std::max(work_size, m_size_3d * (sizeof(double) * 2 + sizeof(dcomplex) * 2));
	m_work = new double[work_size];

	
};

inline
QUMASUN_SINGLE::~QUMASUN_SINGLE() {
	delete[] m_psi_set;	
	
	delete[] m_eigen_values;
	delete[] m_eigen_values_down;
	delete[] m_occupancy;
	delete[] m_occupancy_down;
	
	delete[] m_nuclei;
	delete[] m_nuclei_valence_elecron;
	delete[] m_work;
	delete[] m_lobpcg_keep_p;
	/*
	if (m_hamiltonian_type == HAMILTONIAN::KohnSham_PP) {
		pp_mbk.Release();
	}
	*/
};

inline
bool QUMASUN_SINGLE::mCheckConditions() {
	printf("sizeof(fftw_complex) = %d\n", (int)sizeof(fftw_complex));
	printf("sizeof(std::complex<double>) = %d\n", (int)sizeof(std::complex<double>));
	if (sizeof(fftw_complex) != sizeof(std::complex<double>)) {
		printf("Error: different size, sizeof(fftw_complex) != sizeof(std::complex<double>)\n");
		return false;
	}

	return true;
}


inline
void QUMASUN_SINGLE::mPrintTime() {

	printf("\nCalculation time====================\n");
	const double total_tm = watch.Total({ 0,1,2,3,4,5,6,7 });
	printf("Total time      : %f [s]\n", total_tm);
	watch.Print("mInitialState   :", 0);
	watch.Print("mPrepareCore    :", 1);
	watch.Print("InitialDensity  :", 7);
	watch.Print("mSetOccupancy   :", 2);
	watch.Print("mSetDensity     :", 3);
	watch.Print("mSetPotential   :", 4);
	printf("EigenSolver     :--\n");
	//watch.Print("EigenSolver     :", 5);
	watch.Print("--K-operation   :", 11);
	watch.Print("--V-operation   :", 12);
	watch.Print("--PPNonlocal    :", 13);
	watch.Print("--Matrix        :", 17);
	//watch.Print("--MPI           :", 14);
	watch.Print("--LAPACK        :", 15);
	watch.Print("--LinearComb.   :", 16);
	watch.Print("--others        :", 10);
	watch.Print("mGetTotalEnergy :", 6);
}

inline
void QUMASUN_SINGLE::mFinalize() {

#if 1
	if (m_num_nuclei == 1) {//for debug

		auto fold = [](double x, double box_w) { return (x * 2.0 > box_w) ? x - box_w : (x * 2.0 < -box_w) ? x + box_w : x; };

		FILE* fp = fopen("Vext.txt", "w");
		fprintf(fp, "#r\tVext\n");
		for (int ix = 0; ix < m_size_x; ++ix) {
			size_t i = ix + m_size_x * (ix + m_size_y * ix);
			const double x = (double)ix * m_dx;
			{
				const double dx = fold(m_nuclei[0].Rx - x, m_box_x);
				const double dy = fold(m_nuclei[0].Ry - x, m_box_y);
				const double dz = fold(m_nuclei[0].Rz - x, m_box_z);
				const double r = sqrt(dx * dx + dy * dy + dz * dz);
				fprintf(fp, "%.12f\t%.12f\n", r, m_Vext[i]);
			}
			{
				const double dx = fold(m_nuclei[0].Rx - x, m_box_x);
				const double dy = fold(m_nuclei[0].Ry - 0.0, m_box_y);
				const double dz = fold(m_nuclei[0].Rz - 0.0, m_box_z);
				const double r = sqrt(dx * dx + dy * dy + dz * dz);
				fprintf(fp, "%.12f\t%.12f\n", r, m_Vext[ix]);
			}
						

		}
		fclose(fp);
	}
#endif
#if 0
	if(m_num_nuclei==1){//for debug
		
		FILE* fp = fopen("psi_fin.txt","w");
		fprintf(fp, "#r\tpsi[0]\tVtot\n");
		for (int ix = 0; ix < m_size_x; ++ix) {
			size_t i = ix + m_size_x * (ix + m_size_y * ix);
			const double x = (double)ix * m_dx;
			const double dx = m_nuclei[0].Rx - x;
			const double dy = m_nuclei[0].Ry - x;
			const double dz = m_nuclei[0].Rz - x;
			const double r = sqrt(dx * dx + dy * dy + dz * dz);
						
			fprintf(fp, "%.12f\t%.12f\t%.12f\n", r, m_psi_set[0][i], m_Vtot[i]);
		}
		fclose(fp);
	}
#endif
}
