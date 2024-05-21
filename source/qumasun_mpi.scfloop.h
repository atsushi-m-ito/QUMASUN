#pragma once
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <thread>
#include <chrono>
#include "qumasun_mpi.h"
#include "mpi_helper.h"
#include "qumasun_note.h"

//#define DEBUG_PRINT_ON
#include "debug_print.h"


inline
void QUMASUN_MPI::Execute(int dynamics_step) {

	constexpr double THRESHOLD_DIFF_RHO = 1.0e-7;

	if (!is_construction_successful) return;

	const int proc_id = GetProcessID(m_mpi_comm);

	//initialize////////////////////
	PrintCondition();
	if (! mCheckConditions()) {
		return;
	}
	
	const bool is_root_spin = IsRoot(m_mpi_comm);
	watch.Restart();

    if (dynamics_step == 0) {
        DEBUG_PRINTF("[%d] test0\n", proc_id);
        mInitializeState();
        watch.Record(0);
        DEBUG_PRINTF("[%d] test1\n", proc_id);

        mSetOccupancy();
        watch.Record(2);
        DEBUG_PRINTF("[%d] test2\n", proc_id);
        DEBUG_BARRIER(m_mpi_comm);
        DEBUG_PRINTF("[%d] test2-2\n", proc_id);

        //calculate electron density in real space///////////////
        mInitializeDensity();
        watch.Record(7);
        DEBUG_BARRIER(m_mpi_comm);
        DEBUG_PRINTF("[%d] test4\n", proc_id);
        DEBUG_BARRIER(m_mpi_comm);
    }

    mPrepareCore();
    watch.Record(1);
    DEBUG_PRINTF("[%d] test3\n", proc_id);

	//calculate potential in real space//////////////////////
	mSetPotential();
	watch.Record(4);
	DEBUG_BARRIER(m_mpi_comm);
	DEBUG_PRINTF("[%d] test5\n", proc_id); 
	DEBUG_BARRIER(m_mpi_comm);

#if 0
	//Arnoldi (Lanczos)から初期値を作る方法, 上手くいかない//
	mInitStateArnoldi();
	mSetDensity(true);
	mSetPotential();
	watch.Record(0);
	DEBUG_BARRIER(m_mpi_comm);
	DEBUG_PRINTF("[%d] test0B\n", proc_id);
	DEBUG_BARRIER(m_mpi_comm);
#endif

	double E_tot = 0.0;

	

	//SCF Loop//////////////////////////////////////////////////////
	for (int scf_step = 0; scf_step < LIMIT_SCF_STEP; ++scf_step) {

		if (is_root_spin) {
			printf("Begin SCF Step %d ==============================\n", scf_step + 1); fflush(stdout);
		}
		
		//Solve the Kohn-Sham equation ////////////////
		switch (m_solver) {
		case SOLVER::Lanczos:
			if (is_root_spin) {
				printf("ERROR: Lanczos solver is not supported in mpi run.\n");
			}				
			//mSolveLanczos();
				
			break;
		case SOLVER::LOBPCG:			
			{
				bool is_first_step=(scf_step == 0);
				mSolveLOBPCG(is_first_step ? LOBPCG_STEP_FIRST : LOBPCG_STEP_PER_SCF, scf_step);
				break;
			}
		}
		

		watch.Record(5);
		DEBUG_PRINTF("[%d] test5\n", proc_id); 

		//calculate occupancy of states//////////////////////
		mSetOccupancy();
		watch.Record(2);
		DEBUG_PRINTF("[%d] test6\n", proc_id); 
		double diff_rho = mSetDensity(true);
		diff_rho /= (double)num_electrons;
		watch.Record(3);
		DEBUG_PRINTF("[%d] test7\n", proc_id); 

		//calculate potential in real space//////////////////////
		mSetPotential();
		watch.Record(4);
		DEBUG_PRINTF("[%d] test8\n", proc_id); 
		
		const double prev_E_tot = E_tot;
		//calculate total energy//////////////////////////////
		E_tot = mGetTotalEnergy();
		watch.Record(6);
		DEBUG_PRINTF("[%d] test9\n", proc_id); 
		
		
		//judgement convergence
		bool is_convergence = false;
		if (is_root_spin) {
			if (scf_step > 0) {
				printf("delta E_tot = %.15f\n", prev_E_tot - E_tot);
				printf("\\int |rho-rho_prev| dr / Ne = %.15f\n", diff_rho); fflush(stdout);
			}

			//Check convergence//////////////////////////////////
#if 1
			if (diff_rho < THRESHOLD_DIFF_RHO) {

				is_convergence = true;

				printf("SCF is convergence because of delta rho is smaller than threshold.\n");
			}

#else
			if (fabs(prev_E_tot - E_tot) < 1.0e-7) {

				is_convergence = true;

				printf( "SCF is convergence because of delta Etot is smaller than threshold.\n");
			}
#endif

			printf("End SCF Step %d ==============================\n\n", scf_step + 1); fflush(stdout);
		}

		MPI_Bcast(&is_convergence, 1, MPI_C_BOOL, m_root_id, m_mpi_comm);
		if(is_convergence){
			break;
		}
	}
	//////////////////////////////////////End of SCF Loop//

	//print result/////////////////////////////////////
	if (is_root_spin) {
		printf("%s\n", QUMASUN::note_def_energy); fflush(stdout);
	}

    watch.Restart();
	mGetForce();
    watch.Record(8);

	mPrintTime();


	mFinalize();
	//end of DFT////////////////////////////////////////
}

