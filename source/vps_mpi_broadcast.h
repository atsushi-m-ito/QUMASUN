#pragma once
#ifdef USE_MPI
#include <mpi.h>
#include "mpi_helper.h"
#include "vps_loader.h"

/****************
* rootでファイル読み込みした擬ポテンシャルを、slaveプロセスにコピーする
* 
*********************/
inline
void BroadcastPseudoPot(PseudoPot_MBK* pp, MPI_Comm& mpi_comm, int root_id) {
	const int proc_id = GetProcessID(mpi_comm);

	
	BroadcastAny(mpi_comm, 0, 
		pp->valence_electron,
		pp->num_projectors,
		pp->num_radial_grids,
		pp->cutoff_vlocal,
		pp->cutoff_r[0],
		pp->cutoff_r[1],
		pp->cutoff_r[2],
		pp->cutoff_r[3],
		pp->xi_min,
		pp->xi_delta,
        pp->has_pcc_charge);

#ifdef _DEBUG
	printf("[%d]: pp-bcast: %d, %d, %d, %f, %f, %f, %f, %f, %f\n", 
		proc_id,
		pp->valence_electron,
		pp->num_projectors,
		pp->num_radial_grids,
		pp->cutoff_r[0],
		pp->cutoff_r[1],
		pp->cutoff_r[2],
		pp->cutoff_r[3],
		pp->xi_min,
		pp->xi_delta);
#endif
	

	const int num_grids = pp->num_radial_grids;

	if (proc_id != root_id) {		
		//allocate memory in slave processes//
		pp->projector_quantum_l = new int[pp->num_projectors];
		pp->projector_energy_up = new double[pp->num_projectors * 2];

		
		double* ptr = new double[num_grids * (pp->num_projectors * 2 + 3)];
        pp->buffer_all = ptr;
		pp->radius = ptr; ptr += num_grids;
		pp->V_local = ptr; ptr += num_grids;
		for (int i = 0; i < pp->num_projectors * 2; ++i) {
			pp->projector.push_back( ptr );
            ptr += num_grids;
		}
        pp->pcc_charge = ptr;
	}

	MPI_Bcast((pp->projector_quantum_l), pp->num_projectors, MPI_INT, root_id, mpi_comm);
	MPI_Bcast((pp->projector_energy_up), pp->num_projectors*2, MPI_DOUBLE, root_id, mpi_comm);
	MPI_Bcast((pp->buffer_all), num_grids * (pp->num_projectors * 2 + 3), MPI_DOUBLE, root_id, mpi_comm);

	/*
	pp->xi_min = pp->radius[0];
	pp->xi_delta = (pp->radius[num_grids-1] - pp->xi_min) / (double)(num_grids - 1);
	*/
}

#endif
