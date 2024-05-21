
#pragma once
//#define _USE_MATH_DEFINES
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <complex>
#include "qumasun_td.h"
#include "lobpcg_z_multi_mpi.h"


namespace QUMASUN {
    template <class OperationA, class WATCH>
    inline
    void TimeEvoH4(const GridRangeMPI& l_grid, const double dVol, const double dt, 
            const int num_solution, SoAComplex* x,
            OperationA OpeA, WATCH& watch)
    {

        const int proc_id = GetProcessID(l_grid.mpi_comm);
        const bool is_root = (proc_id == 0);
        const size_t local_size = l_grid.Size3D();

        double* buffer = new double[local_size * 4];

        SoAComplex Hx{ buffer, buffer + local_size };
        SoAComplex H2x{ buffer + local_size * 2, buffer + local_size * 3 };
        
        for (int n = 0; n < num_solution; ++n) {
            OpeA(Hx, x[n]);

            for (int i = 0; i < local_size; ++i) {
                x[n].re[i] += dt * Hx.im[i];
            }
            for (int i = 0; i < local_size; ++i) {
                x[n].im[i] += -dt * Hx.re[i];
            }

            OpeA(H2x, Hx);
            const double coef2 = dt * dt / 2.0;
            for (int i = 0; i < local_size; ++i) {
                x[n].re[i] += -coef2 * H2x.re[i];
            }
            for (int i = 0; i < local_size; ++i) {
                x[n].im[i] += -coef2 * H2x.im[i];
            }

            auto& H3x = Hx;
            OpeA(H3x, H2x);
            const double coef3 = coef2 * dt / 3.0;
            for (int i = 0; i < local_size; ++i) {
                x[n].re[i] += -coef3 * H3x.im[i];
            }
            for (int i = 0; i < local_size; ++i) {
                x[n].im[i] += coef3 * H3x.re[i];
            }

            auto& H4x = H2x;
            OpeA(H4x, H3x);
            const double coef4 = coef3 * dt / 4.0;
            for (int i = 0; i < local_size; ++i) {
                x[n].re[i] += coef4 * H4x.re[i];
            }
            for (int i = 0; i < local_size; ++i) {
                x[n].im[i] += coef4 * H4x.im[i];
            }

        }
        
        delete[] buffer;
    }
}

inline
size_t QUMASUN_TD::mWorkSizeTimeEvoH4(int local_size, int num_solution) {
    return 0; // means for complex //
}

//mpi supported//
//四次精度の時間発展演算子を作用する//
inline
void QUMASUN_TD::mTimeEvolutionH4(double time_step_dt) {
	const int proc_id = GetProcessID(m_mpi_comm);
	const bool is_ddm_root = IsRoot(m_ddm_comm);
	const size_t local_size = ml_grid.Size3D();


	for(int sk = 0; sk < m_num_having_spin_kpoint;++sk){
		const int kpoint_x = ml_wave_set[sk].kpoint_x;
		const int kpoint_y = ml_wave_set[sk].kpoint_y;
		const int kpoint_z = ml_wave_set[sk].kpoint_z;
		auto& V_tot = (ml_wave_set[sk].spin == SPIN::UP) ? ml_Vtot : ml_Vtot_down;

        
        
        QUMASUN::TimeEvoH4(ml_grid, m_dx * m_dy * m_dz, time_step_dt, num_solution,
			ml_wave_set[sk].l_psi_set, 
            [&V_tot, &kpoint_x, &kpoint_y, &kpoint_z, &sk, this](SoAComplex& Ax, const SoAComplex& x) {
                this->mHamiltonianMatrix_ddm(Ax, V_tot, x, kpoint_x, kpoint_y, kpoint_z, sk);
            }, watch);
	}
	
}


