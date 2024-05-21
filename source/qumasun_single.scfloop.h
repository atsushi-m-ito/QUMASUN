#pragma once
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include "qumasun_single.h"
#include "qumasun_note.h"

inline
void QUMASUN_SINGLE::Execute() {
	watch.Restart();
	//initialize////////////////////
	if (! mCheckConditions()) {
		return;
	}

	mInitializeState();
	watch.Record(0);

	mSetOccupancy();
	watch.Record(2);

	mPrepareCore();
	watch.Record(1);

	//calculate electron density in real space///////////////
	mInitializeDensity();
	watch.Record(7);


	double E_tot = 0.0;

	//SCF Loop//////////////////////////////////////////////////////
	for (int scf_step = 0; scf_step < LIMIT_SCF_STEP; ++scf_step) {

		printf("Begin SCF Step %d ==============================\n", scf_step + 1);


		//calculate potential in real space//////////////////////
		mSetPotential();
		watch.Record(4);

		//Solve the Kohn-Sham equation ////////////////
		switch (m_solver) {
		case SOLVER::Lanczos:
			mSolveLanczos();
			break;
		case SOLVER::LOBPCG:
		{
			bool is_first_step = (scf_step == 0);
			mSolveLOBPCG((is_first_step) ? LOBPCG_STEP_FIRST : LOBPCG_STEP_PER_SCF, is_first_step);
			break;
		}
		}
		watch.Record(5);

		//calculate occupancy of states//////////////////////
		mSetOccupancy();
		watch.Record(2);

		//calculate electron density in real space///////////////
		double diff_rho = mSetDensity(true);
		diff_rho /= (double)num_electrons;
		watch.Record(3);
		
		//print log////////////////////////
		if (!is_spin_on) {
			for (int i = 0; i < num_solution; ++i) {
				printf("eigen[%d] = %f\n", i, m_eigen_values[i]); fflush(stdout);
			}
		} else {
			for (int i = 0; i < num_solution; ++i) {
				printf("eigen[%d] = %f, %f\n", i, m_eigen_values[i], m_eigen_values_down[i]); fflush(stdout);
			}
		}


		//calculate total energy//////////////////////////////
		const double prev_E_tot = E_tot;
		E_tot = mGetTotalEnergy();
		watch.Record(6);

		if (scf_step > 0) {
			printf("delta E_tot = %.15f\n", prev_E_tot - E_tot);
			printf("\\int |rho-rho_prev| dr = %.15f\n", diff_rho);
		}

		printf( "End SCF Step %d ==============================\n\n", scf_step + 1);

#if 1
		if (diff_rho < 1.0e-7) {
			//Check convergence//////////////////////////////////
			if (true) {
				printf("SCF is convergence because of delta rho is smaller than threshold.\n");
				break;
			}
		}
#else
		if (fabs(prev_E_tot - E_tot) < 1.0e-7) {
			//Check convergence//////////////////////////////////
			if (true) {
				printf("SCF is convergence because of delta Etot is smaller than threshold.\n");
				break;
			}
		}
#endif
	}
	//////////////////////////////////////End of SCF Loop//

	//print result/////////////////////////////////////
	printf("%s\n", QUMASUN::note_def_energy); fflush(stdout);
	mPrintTime();

	mFinalize();
	//end of DFT////////////////////////////////////////
}
