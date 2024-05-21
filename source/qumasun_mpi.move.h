#pragma once
#include "qumasun_mpi.h"
#include "nucleus.h"

inline
void QUMASUN_MPI::MoveNuclei(int num_nuclei, const Nucleus* next_nucleis) {
    watch.Restart();

    for (int ni = 0; ni < num_nuclei; ++ni) {
        m_nuclei[ni] = next_nucleis[ni];
    }

    if (m_hamiltonian_type == HAMILTONIAN::KohnSham_PP) {

        m_pp_integrator.UpdatePosition(m_nuclei, m_num_nuclei, ml_grid);

    }

    watch.Record(20);

}

