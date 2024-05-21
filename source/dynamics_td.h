#pragma once
#include <mpi.h>
#include "mpi_helper.h"
#include "qumasun_input.h"
#include "Relaxation2.h"
#include "folding.h"
#include "md3_writer2.h"


void FprintRVM(const char* filepath, int num_nuclei, Nucleus* nuclei, const vec3d* p, const double* m, const double* box_axis )
{
    FILE* fp = fopen(filepath, "w");

    fprintf(fp, "#Simulation Box===========================\n\n");
    
    fprintf(fp, "Material.SuperCell     1 1 1\n");
    fprintf(fp, "Material.UnitLatticeVector.Unit    Bohr\n");
    fprintf(fp, "Material.UnitLatticeVector.Begin\n");
    fprintf(fp, "  %.15f  %.15f  %.15f\n", box_axis[0], box_axis[1], box_axis[2]);
    fprintf(fp, "  %.15f  %.15f  %.15f\n", box_axis[3], box_axis[4], box_axis[5]);
    fprintf(fp, "  %.15f  %.15f  %.15f\n", box_axis[6], box_axis[7], box_axis[8]);
    fprintf(fp, "Material.UnitLatticeVector.End\n");
    fprintf(fp, "\n\n");

    fprintf(fp, "#Last Positons of Nuclei===========================\n\n");   

    fprintf(fp, "Material.UnitCell.Unit    Bohr\n");
    fprintf(fp, "Material.UnitCell.Begin\n");
    for (int i = 0; i < num_nuclei; ++i) {
        fprintf(fp, "  %s  %.15f  %.15f  %.15f\n", msz::GetAtomicSymbol(nuclei[i].Z), nuclei[i].Rx, nuclei[i].Ry, nuclei[i].Rz);
    }
    fprintf(fp, "Material.UnitCell.End\n");
    fprintf(fp, "\n\n");

    fprintf(fp, "#Last Mass and Velocity of Nuclei===========================\n\n");

    fprintf(fp, "Dynamics.Mass.Unit        u  # u, Da, me, or a.u.\n");
    fprintf(fp, "Dynamics.Velocity.Unit    a.u.  # m / s, or a.u.\n");
    fprintf(fp, "\n");
    fprintf(fp, "#note: velocity of 1[eV] and 1[u] atom corresponds to v = 13891.38738921143[m / s] = 0.0063497933278141[a.u.]\n");
    fprintf(fp, "#note : if mass is m[u] and energy is E[eV], v' = v * sqrt(E/m)\n");
    fprintf(fp, "\n");
    fprintf(fp, "Dynamics.Velocity.Begin\n");
    for (int i = 0; i < num_nuclei; ++i) {
        fprintf(fp, "  %.15f  %.15f  %.15f  %.15f\n", m[i] / mass_u_as_au, (p[i].x/m[i]), (p[i].y / m[i]), (p[i].z / m[i]));
    }
    fprintf(fp, "Dynamics.Velocity.End\n");

    fprintf(fp, "\n");

    fclose(fp);

}



template<class QUMASUN_T>
int TimeEvoQUMASUN(QUMASUN_T& qumasun, QUMASUN::Input& input, MPI_Comm& mpi_comm, QUMASUN::DynamicsMode mode, const int* ddm_num) {
	using namespace QUMASUN;
	
    const bool is_root = IsRoot(mpi_comm);
    int proc_id = GetProcessID(mpi_comm);


    auto RhoFileName = [](int istep, int max_step) {
        int64_t length = std::to_string(max_step).length();
        std::string filename("rho_t");
        length -= std::to_string(istep).length();
        if (length > 0) {
            filename += std::string(length, '0');
        }
        filename += std::to_string(istep) + ".cube";
        return filename;
        };

	switch (mode) {
		case DynamicsMode::TDDFT:
		{

            const double dt = input.dynamics_timestep_dt;
            const int total_steps = input.dynamics_step;
            const int incremental_step = input.dynamics_output_step;

            for (int istep = 0; istep < total_steps; istep += incremental_step) {
                qumasun.Evolve(istep, incremental_step, dt);
                qumasun.OutputDensity(QUMASUN::OUTPUT_TARGET::Density, RhoFileName(istep+ incremental_step, total_steps).c_str());
            }
			

            break;
		}		
        case DynamicsMode::EhrenfestMDv1:
        {

            const double* box_axis =input.box_axis;

            const double dt = input.dynamics_timestep_dt;
            const int start_steps = input.dynamics_start_step;
            const int total_steps = input.dynamics_step + start_steps;
            const int incremental_step = input.dynamics_output_step;

            
            msz::MD3_Writer2<double> md3;
            md3.Open("nucl_position.md3");
            auto WriteMD3_one = [&md3, &box_axis](int num_nuclei, Nucleus* nuclei) {
                std::vector<int> atoms_Z(num_nuclei);
                std::vector<vec3d> r(num_nuclei);
                for (int i = 0; i < num_nuclei; ++i) {
                    atoms_Z[i] = nuclei[i].Z;
                    r[i].x = nuclei[i].Rx;
                    r[i].y = nuclei[i].Ry;
                    r[i].z = nuclei[i].Rz;
                }

                md3.BeginFrame(num_nuclei);
                md3.WriteZ(&atoms_Z[0]);
                md3.WriteR(&r[0]);
                md3.WriteBox(box_axis);
                md3.EndFrame();
                };
            WriteMD3_one(input.num_nuclei, &input.nuclei[0]);

            //show velocity//
            if (is_root) {
                printf("Mass and Velocity [a.u.] ===========================\n");
                for (int i = 0; i < input.num_nuclei; ++i) {
                    printf("%d: %f, %f, %f, %f\n", i, input.mass[i], input.velocity[i*3], input.velocity[i * 3+1], input.velocity[i * 3+2]);
                }
                printf("====================================================\n");
                fflush(stdout);
            }
            std::vector<vec3d> P(input.num_nuclei );
            std::vector<vec3d> F(input.num_nuclei );
            for (int i = 0; i < input.num_nuclei; ++i) {
                P[i].x = input.velocity[i * 3] * input.mass[i];
                P[i].y = input.velocity[i * 3+1] * input.mass[i];
                P[i].z = input.velocity[i * 3+2] * input.mass[i];
            }

            auto EvolveR = [&box_axis](int num_nuclei, Nucleus* nuclei, const vec3d* p, const double* m, double dt) {
                for (int i = 0; i < num_nuclei; ++i) {
                    nuclei[i].Rx += p[i].x / m[i] * dt;
                    nuclei[i].Ry += p[i].y / m[i] * dt;
                    nuclei[i].Rz += p[i].z / m[i] * dt;

                    Folding(nuclei[i].Rx, nuclei[i].Ry, nuclei[i].Rz, box_axis);
                }
                };
            auto EvolveP = [](int num_nuclei, vec3d* p, const vec3d* f, double dt) {
                for (int i = 0; i < num_nuclei; ++i) {
                    p[i].x += f[i].x * dt;
                    p[i].y += f[i].y * dt;
                    p[i].z += f[i].z * dt;
                }
                };
            
            auto KineticEnergy = [](int num_nuclei, const vec3d* p, const double* m) {
                double K = 0.0;
                for (int i = 0; i < num_nuclei; ++i) {
                    K += ((p[i].x * p[i].x) + (p[i].y * p[i].y) + (p[i].z * p[i].z)) / (2.0 * m[i]);
                }
                return K;
                };

            for (int istep = start_steps; istep < total_steps; ++istep) {
                EvolveR(input.num_nuclei, &input.nuclei[0], &P[0], &input.mass[0], dt/2.0);
                MPI_Bcast(&input.nuclei[0], sizeof(Nucleus) * input.num_nuclei, MPI_BYTE, 0, mpi_comm);
                qumasun.MoveNuclei(input.num_nuclei, &input.nuclei[0]);
                
                const double E_ele = qumasun.Evolve(istep, 1, dt);
                //MPI_Barrier(mpi_comm);
                //printf("[%d]emd:4\n", proc_id); fflush(stdout);
                //MPI_Barrier(mpi_comm);

                qumasun.GetForce(&F[0]);				//力の計算.
                //MPI_Barrier(mpi_comm);
                //printf("[%d]emd:5\n", proc_id); fflush(stdout);
                //MPI_Barrier(mpi_comm);

                EvolveP(input.num_nuclei, &P[0], &F[0], dt);
                const double E_kin_nucl = KineticEnergy(input.num_nuclei, &P[0], &input.mass[0]);
                if (is_root) {
                    printf("MD_E_kin_nucl: %d, %f\n", istep + 1, E_kin_nucl);
                    printf("MD_E_electron: %d, %f\n", istep + 1, E_ele);
                    printf("MD_E_total   : %d, %f\n", istep + 1, E_kin_nucl+ E_ele);
                    fflush(stdout);
                }

                EvolveR(input.num_nuclei, &input.nuclei[0], &P[0], &input.mass[0], dt / 2.0);


                if (incremental_step > 0) {
                    if ((istep + 1) % incremental_step == 0) {
                        qumasun.OutputDensity(QUMASUN::OUTPUT_TARGET::Density, RhoFileName(istep + 1, total_steps).c_str());

                        WriteMD3_one(input.num_nuclei, &input.nuclei[0]);
                    }
                }
                //MPI_Barrier(mpi_comm);
                //printf("[%d]emd:8\n", proc_id); fflush(stdout);
                //MPI_Barrier(mpi_comm);
            }
            MPI_Bcast(&input.nuclei[0], sizeof(Nucleus) * input.num_nuclei, MPI_BYTE, 0, mpi_comm);
            qumasun.MoveNuclei(input.num_nuclei, &input.nuclei[0]);


            FprintRVM("final_positions.txt", input.num_nuclei, &input.nuclei[0], &P[0], &input.mass[0], box_axis);
            qumasun.OutputDensity(QUMASUN::OUTPUT_TARGET::Density, "final_rho.cube");




            break;
        }
        default:
        {
            return -1;
        }
        
	}

    
    if (input.output_state_vector != "none") {
        qumasun.OutputEigenVector(input.output_state_vector.c_str());
    }
    if (is_root) {
        qumasun.PrintTime();
    }

    return 0;

}
