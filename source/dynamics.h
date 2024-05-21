#pragma once
#include <mpi.h>
#include "mpi_helper.h"
#include "qumasun_input.h"
#include "Relaxation2.h"
#include "folding.h"
#include "atomic_number.h"

void PrintPositions(int num_nuclei, Nucleus* nuclei)
{    
    printf("Positons of Nuclei===========================\n");
    for (int i = 0; i < num_nuclei; ++i) {
        printf("  %s\t%.15f\t%.15f\t%.15f\n", msz::GetAtomicSymbol(nuclei[i].Z), nuclei[i].Rx, nuclei[i].Ry, nuclei[i].Rz);
    }
    printf("\n");
}

template<class QUMASUN_T>
int ExecuteQUMASUN(QUMASUN_T& qumasun, QUMASUN::Input& input, MPI_Comm& mpi_comm, QUMASUN::DynamicsMode mode, const int* ddm_num) {
	using namespace QUMASUN;
	
    const bool is_root = IsRoot(mpi_comm);

	switch (mode) {
		case DynamicsMode::None:
		{
			qumasun.Execute();
			qumasun.OutputDensity(QUMASUN::OUTPUT_TARGET::Density, "test.cube");
            qumasun.OutputDensity(QUMASUN::OUTPUT_TARGET::DensityHR, "test_HR.cube");
			qumasun.OutputEigenValue("test_eigenvalue.txt");
            break;
		}
		case DynamicsMode::Test:
		{
            if (is_root) {
                printf("====================================\n"); 
                printf("Dynamics Step = 0\n");
                fflush(stdout);
            }

			qumasun.Execute();

			const int num_nuclei = input.num_nuclei;
			std::vector<Nucleus> nuclei = input.nuclei; //copy

            const int STEPS = input.dynamics_step;
            double delta_x = input.box_axis[0] / (double)(input.grid_size[0] * (STEPS));

			for (int istep = 1; istep <= STEPS; ++istep) {
                if (is_root) {
                    printf("====================================\n");
                    printf("Dynamics Step = %d\n", istep);
                    fflush(stdout);
                }

				for (int i = 0; i < num_nuclei; ++i) {
					nuclei[i].Rx += delta_x;
                    Folding(nuclei[i].Rx, nuclei[i].Ry, nuclei[i].Rz, input.box_axis);
				}
				qumasun.MoveNuclei(num_nuclei, &nuclei[0]);
				qumasun.Execute(istep);
			}
			qumasun.OutputDensity(QUMASUN::OUTPUT_TARGET::Density, "test.cube");
            qumasun.OutputEigenValue("test_eigenvalue.txt");

            break;
		}

        case DynamicsMode::Relaxation:
        {

            const int num_nuclei = input.num_nuclei;
            std::vector<Nucleus> nuclei = input.nuclei; //copy
            glips::Relaxation2Q relax(mpi_comm);
            glips::Relaxation2Q::Parameters params;
            params.lower_force_limit = input.dynamics_force_threshold;
            relax.Reset(&params);
            const int steps = input.dynamics_step;
            double init_dt = 0.1;
            relax.Evolve(qumasun, steps, init_dt, num_nuclei, &nuclei[0], input.box_axis, 0, nullptr);
            qumasun.OutputDensity(QUMASUN::OUTPUT_TARGET::Density, "test.cube");
            qumasun.OutputEigenValue("test_eigenvalue.txt");

            if (is_root) {
                PrintPositions(num_nuclei, &nuclei[0]);
            }
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

    return 0;

}
