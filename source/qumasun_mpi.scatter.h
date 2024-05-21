#pragma once
#include "qumasun_mpi.h"
#include "GridRange.h"
#include "GridScatterGather.h"

inline
void QUMASUN_MPI::mScatterField(double* local_dest, const double* global_src) {
	ScatterGrid(ml_grid, local_dest, m_global_grid, global_src, m_root_id);
}

inline
void QUMASUN_MPI::mGatherField(double* global_dest, const double* local_src) {	
	GatherGrid(m_global_grid, global_dest, ml_grid, local_src, m_root_id);
}


inline
void QUMASUN_MPI::mScatterField_HR(double* local_dest, const double* global_src) {
    ScatterGrid(ml_grid, local_dest, m_global_grid, global_src, m_root_id, m_HR_ratio_x, m_HR_ratio_y, m_HR_ratio_z);
}

inline
void QUMASUN_MPI::mGatherField_HR(double* global_dest, const double* local_src) {
    GatherGrid(m_global_grid, global_dest, ml_grid, local_src, m_root_id, m_HR_ratio_x, m_HR_ratio_y, m_HR_ratio_z);
}



/*
* kpoint=0,0,0, spin=UPのマスタープロセス(場所の異なるものが複数ある)でDDM通信でscatter
* 次に、kpoint, spinは異なるがDDMの担当領域は同じプロセスへとコピー.
*/

inline
void QUMASUN_MPI::mHierarchyScatterField(double* local_dest, const double* global_src) {
    if (is_spin_on && (m_num_spin == 1)) {
        const int proc_id = GetProcessID(m_mpi_comm);
        const int num_procs = GetNumProcess(m_mpi_comm);
        const int local_size = ml_grid.Size3D();
        const int TAG = 2201;

        if (proc_id < num_procs / 2) {
            ScatterGrid(ml_grid, local_dest, m_global_grid, global_src, m_root_id);
            MPI_Send(local_dest, local_size, MPI_DOUBLE, proc_id + num_procs / 2, TAG, m_mpi_comm);
        } else {
            MPI_Status status;
            MPI_Recv(local_dest, local_size, MPI_DOUBLE, proc_id - num_procs / 2, TAG, m_mpi_comm, &status);
        }
    } else {
        ScatterGrid(ml_grid, local_dest, m_global_grid, global_src, m_root_id);
    }
}
