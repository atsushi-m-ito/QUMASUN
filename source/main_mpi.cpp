/*******************************

QUantum MAterial Simulation UNraveler (QUMASUN)

QUMASUN is numerical simulation code for Density Functional Theory(DFT) and Time-dependent DFT based on the real space grid.

********************************/
#include <mpi.h>
#include "mpi_helper.h"

#include <cstdlib>
#include <cstdio>
#define TIME_PP_MPI
#include "qumasun_mpi.h"
#include "qumasun_kpoint.h"
#include "qumasun_td.h"
#include "atomic_number.h"
#include "field_interpolation.h"
#include "GridRange.h"
#include "GridFor.h"
#include "GetArg.h"
#include "GridScatterGather.h"
#include "qumasun_make_input.h"
#include "dynamics.h"
#include "dynamics_td.h"



//simple DFT

int main(int argc, char* argv[]) {

	MPI_Init(&argc, &argv);
	MPI_Comm mpi_comm = MPI_COMM_WORLD;
	const int proc_id = GetProcessID(mpi_comm);
	const int num_procs = GetNumProcess(mpi_comm);
	if (argc < 2) {
		if (proc_id == 0) {
			printf("ERROR: no input file.\n");
		}
		MPI_Finalize();
		return -1;
	}


	//for MPI/////////////////////////////////////
	int ddm_num[4] = { 1,1,1, 1 }; //4th element is for state decomposition//
	GetArgumentNumList<int>(argc, argv, "-ddm", 3, ddm_num);
    ddm_num[3] = GetArgumentNum<int>(argc, argv, "-sd", 1);
	

	
	//begin calculation///////////////////////////////
	QUMASUN::Input input = QUMASUN::MakeInput(argv[1]);

    if (input.dynamics_mode >= QUMASUN::DynamicsMode::TDDFT) {

        if (IsRoot(mpi_comm)) {
            printf("Begin calculation TDDFT on QUMASUN\n");
        }

        QUMASUN_TD qumasun(input, mpi_comm, ddm_num);
        TimeEvoQUMASUN(qumasun, input, mpi_comm, input.dynamics_mode, ddm_num);
        

    }else if (GetArgumentText(argc, argv, "-kpoint", nullptr) || (input.kpoint_sample[0] * input.kpoint_sample[1] * input.kpoint_sample[2] > 1)) {

		if (IsRoot(mpi_comm)) {
			printf("Begin calculation on QUMASUN with kpoint\n");
		}

        QUMASUN_KPOINT qumasun(input, mpi_comm, ddm_num);
        ExecuteQUMASUN<QUMASUN_KPOINT>(qumasun, input, mpi_comm, input.dynamics_mode, ddm_num);
        

	} else {
		if (IsRoot(mpi_comm)) {
			printf("Begin calculation on QUMASUN on gamma point\n");fflush(stdout);
		}

		QUMASUN_MPI qumasun(input, mpi_comm, ddm_num);
        ExecuteQUMASUN<QUMASUN_MPI>(qumasun, input, mpi_comm, input.dynamics_mode, ddm_num);
    }

    MPI_Barrier(mpi_comm);
	MPI_Finalize();

}

