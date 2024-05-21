#pragma once
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include "qumasun_td.h"
#include "mpi_helper.h"
#include "qumasun_note.h"

inline
double QUMASUN_TD::Evolve(int dynamics_step, int incremental_steps, double time_step_dt) {
	constexpr double THRESHOLD_DIFF_RHO = 1.0e-7;
	
	if (!is_construction_successful) return 0.0;

	const int proc_id = GetProcessID(m_mpi_comm);

	//initialize////////////////////
	PrintCondition();
	if (! mCheckConditions()) {
		return 0.0;
	}

	const bool is_root_global = IsRoot(m_mpi_comm);
    watch.Restart();

    if (!is_state_initialized) {
        //printf("[%d] test0\n", proc_id); fflush(stdout);
        mInitializeState();
        watch.Record(0);
        //printf("[%d] test1\n", proc_id); fflush(stdout);

        mSetOccupancy();
        watch.Record(2);
        //printf("[%d] test2\n", proc_id); fflush(stdout);
        //MPI_Barrier(m_mpi_comm);
        //printf("[%d] test2-2\n", proc_id); fflush(stdout);
    
        //calculate electron density in real space///////////////
        mInitializeDensity();
        watch.Record(7);
        //MPI_Barrier(m_mpi_comm);
        //printf("[%d] test3-2\n", proc_id); fflush(stdout);


        mPrepareCore();
        watch.Record(1);
        //printf("[%d] test3\n", proc_id); fflush(stdout);

        //calculate potential in real space//////////////////////
        mSetPotential(true);
        watch.Record(4);
    } else {

        mPrepareCore();
        watch.Record(1);
        //printf("[%d] test3\n", proc_id); fflush(stdout);

        //calculate potential in real space//////////////////////
        mSetPotential(false); //no-update of V_hart and V_xc
        watch.Record(4);
    }

    //MPI_Barrier(m_mpi_comm);
    //printf("[%d] test3-3\n", proc_id); fflush(stdout);

	double E_tot = 0.0;

    if (!is_state_initialized) {
        if (is_root_global) {
            printf("Initial Energy ==============================\n"); fflush(stdout);
        }
        E_tot = mGetTotalEnergy(false);
        m_current_Etot = E_tot;
        watch.Record(6);

        //mGetForce();
        //watch.Record(8);

        if (is_root_global) {
            printf("==============================\n\n"); fflush(stdout);
        }
        
        is_state_initialized = true;
    }


    if (is_root_global) {
        printf("Begin TDDFT Step %d ==============================\n", dynamics_step + 1); fflush(stdout);
    }
	//SCF Loop//////////////////////////////////////////////////////
	for (int tddft_step = 0; tddft_step < incremental_steps; ++tddft_step) {

        watch.Restart();
        mTimeEvolutionH4(time_step_dt);
		watch.Record(5);
		//printf("[%d] test5\n", proc_id); fflush(stdout);
		//calculate occupancy of states//////////////////////
		//mSetOccupancy();
		//watch.Record(2);
		//printf("[%d] test6\n", proc_id); fflush(stdout);

		double diff_rho = mSetDensity(false);//no-mixing
		diff_rho /= (double)num_electrons;
		watch.Record(3);
		//printf("[%d] test7\n", proc_id); fflush(stdout);

		//calculate potential in real space//////////////////////
		mSetPotential(true);
		watch.Record(4);
		//printf("[%d] test8\n", proc_id); fflush(stdout);
		
	}
	//////////////////////////////////////End of SCF Loop//


    const double prev_E_tot = m_current_Etot;
    //calculate total energy//////////////////////////////
    E_tot = mGetTotalEnergy(true);
    m_current_Etot = E_tot;
    watch.Record(6);
    //printf("[%d] test9\n", proc_id); fflush(stdout);


    if (is_root_global) {        
        printf("delta E_tot = %.15f\n", prev_E_tot - E_tot);
        printf("End TDDFT Step %d ==============================\n\n", dynamics_step + incremental_steps); fflush(stdout);
    }

	//print result/////////////////////////////////////
	//if (is_root_global) {
	//	printf("%s\n", QUMASUN::note_def_energy); fflush(stdout);
	//}

    //MPI_Barrier(m_mpi_comm);
    DEBUG_PRINTF("[%d]mGetForce:before\n", proc_id);
	mGetForce(true);
    watch.Record(8);
    DEBUG_PRINTF("[%d]mGetForce:after\n", proc_id);
    //MPI_Barrier(m_mpi_comm);
    DEBUG_PRINTF("[%d]mPrintTime\n", proc_id);
	mPrintTime();

    DEBUG_PRINTF("[%d]end mPrintTime\n", proc_id);
	
	mFinalize();
	//end of DFT////////////////////////////////////////

    DEBUG_PRINTF("[%d]end mFinalize\n", proc_id);
    return m_current_Etot;
}

#undef DEBUG_PRINT_ON
