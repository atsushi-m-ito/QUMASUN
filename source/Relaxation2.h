#pragma once
#include "mpi.h"
#include <float.h>
#include <vector>
#include "vec3.h"
#include "nucleus.h"

namespace glips {

    /*
    * This class is based on I2_Relaxation2 class in GLIPS,
    * and the class is modified to call QUMASUN
    */
    class Relaxation2Q
    {
    private:
        double m_lower_force_limit = 1.0e-2;
        double m_upper_move_limit = 0.1;

        MPI_Comm m_mpi_comm;
        Relaxation2Q() = delete;
    public:

        Relaxation2Q(MPI_Comm mpi_comm) :
            m_mpi_comm(mpi_comm)
        {

        };

        void Folding(double& rx, double& ry, double& rz, const double* box_axis) {
            double ix = rx / box_axis[0];
            ix -= std::floor(ix);
            rx = ix * box_axis[0];

            double iy = ry / box_axis[4];
            iy -= std::floor(iy);
            ry = iy * box_axis[4];

            double iz = rz / box_axis[8];
            iz -= std::floor(iz);
            rz = iz * box_axis[8];
        }

        /*
        指定したステップだけ構造緩和を行う
        */
        template<class QUMASUN_T>
        void Evolve(QUMASUN_T& qumasun, const int num_steps, double init_dt, int num_nuclei, Nucleus* nuclei, double* box_axis, const int fix_count, double* total_energy)
        {
            const bool is_root = IsRoot(m_mpi_comm);


            //dt *= m_dt_scale;            
            //const double dt_half = init_dt / 2.0;
            double prev_U = 0.0;//DBL_MAX
            double max_f = 0.0;
            double dU = 0.0;
            double dt = init_dt;


            double prev_abs_f = 1.0;
            std::vector<vec3d> prev_f(num_nuclei);
            std::vector<vec3d> f(num_nuclei);

            {
                for (int i = 0; i < num_nuclei; ++i) {
                    prev_f[i].Clear();
                }

                /*
                note: DDM使用時はForceの計算によって粒子のソートや入れ替えが起こるため、
                m_prev_fもソートする必要がある。困った.
                解決:運動量pを前回の力を格納するバッファにすることで解決
                */
            }



            for (int istep = 0; istep < num_steps; istep++) {

                if (is_root) {
                    printf("====================================\n");
                    printf("Relaxation Step = 0\n");
                    fflush(stdout);
                }

                qumasun.Execute(istep);

                if (is_root) {
                    qumasun.GetForce(&f[0]);				//力の計算.
#if 0
                    if (fix_count > 0) {
                        atom_container->ClearFByState();
                    }
#endif

                    double U = 0.0;
                    //potential->GetU(&U);
                    dU = U - prev_U;

                    prev_U = U;
                    //const double* m = atom_container->M(); mass is 1.0, tentatively.

                    double inner_f_prev_f = 0.0;
                    double abs_f = 0.0;
                    max_f = 0.0;
                    for (int i = 0; i < num_nuclei; i++) {
                        //if (m[i] != 0.0) 
                        {
                            double f2 = f[i] * f[i];
                            if (max_f < f2) {
                                max_f = f2;
                            }
                            abs_f += f2;
                            inner_f_prev_f += f[i] * prev_f[i];
                        }
                    }


                    max_f = sqrt(max_f);
                    abs_f = sqrt(abs_f);

                    uint8_t check_convergence = 0;

                    printf("\nRelaxation info===================================\n");
                    printf("istep = %d, max_f = %g, dU = ---\n", istep, max_f);
                    if (max_f < m_lower_force_limit) {
                        printf("force-convergence: max_f = %g < lower_force_limit = %g\n", max_f, m_lower_force_limit);
                        if (istep == 0) {
                            dU = 0.0;
                        }
                        check_convergence = 1;
                    }

                    
                    MPI_Bcast(&check_convergence, 1, MPI_UINT8_T, 0, m_mpi_comm);
                    if (check_convergence) {
                        break;
                    }


                    /*dtの増減/////////////////////
                      前回のforceと向きが同じならdtを増やし、逆なら減らす。直行なら変化なし。
                      これを内積から決める
                      */
                    const double x = inner_f_prev_f / (abs_f * prev_abs_f);
                    const double f1 = 0.5;
                    const double fm1 = -0.5;
                    const double a2 = (f1 + fm1) / 2.0;
                    const double a1 = (f1 - fm1) / 2.0;
                    const double ratio = 1.0 + a1 * x + a2 * x * x;
                    dt *= ratio;
                    prev_abs_f = abs_f;

                    dt = ((max_f * dt > m_upper_move_limit) ? m_upper_move_limit / max_f : dt);

#if 0
                    if (use_state == UseState::FixAtom) {
                        for (int i = fix_count; i < num_nuclei; i++) {
                            if (atom_container->State(i) ^ IAtomContainer::STATE_FLAG_FIX) {
                                r[i] += dt * f[i];
                                prev_f[i] = f[i];
                            }
                        }
                    } else
#endif
                    {
                        for (int i = 0; i < num_nuclei; i++) {
                            //if (m[i] != 0.0) 
                            {
                                nuclei[i].Rx += dt * f[i].x;
                                nuclei[i].Ry += dt * f[i].y;
                                nuclei[i].Rz += dt * f[i].z;


                                //should be folding 
                                Folding(nuclei[i].Rx, nuclei[i].Ry, nuclei[i].Rz, box_axis);

                                prev_f[i] = f[i];
                            }
                        }
                    }
                } else { //slave process//
                    uint8_t check_convergence = 0;
                    MPI_Bcast( &check_convergence, 1, MPI_UINT8_T, 0, m_mpi_comm);
                    if (check_convergence) {
                        break;
                    }
                }

                //nest step//
                MPI_Bcast(nuclei, sizeof(Nucleus) * num_nuclei, MPI_BYTE, 0, m_mpi_comm);
                qumasun.MoveNuclei(num_nuclei, &nuclei[0]);
                

            }//end of loop "istep"

            if (total_energy) {

                double U=0.0;
                //potential->GetU(&U);
                total_energy[0] = U;
                total_energy[1] = (num_steps <= 1) ? 0.0 : dU;//dt//
                total_energy[2] = max_f;
            }


        }


        struct Parameters {
            double lower_force_limit;
        };

        /*
        Integratorの再初期化
        */
        void Reset(const void* si_vars) {
            const Parameters& params = *(const Parameters*)si_vars;
            m_lower_force_limit = params.lower_force_limit;
        }
        void OutputLog(const char* filepath) {};


    
    };

}//namespace//
