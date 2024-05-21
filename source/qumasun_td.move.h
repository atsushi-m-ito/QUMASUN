#pragma once
#include "qumasun_td.h"
#include "nucleus.h"

inline
void QUMASUN_TD::MoveNuclei(int num_nuclei, const Nucleus* next_nucleis) {
    
    watch.Restart();

    for (int ni = 0; ni < num_nuclei; ++ni) {
        m_nuclei[ni] = next_nucleis[ni];
    }

    if (m_hamiltonian_type == HAMILTONIAN::KohnSham_PP) {

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

    }
    watch.Record(20);

}

