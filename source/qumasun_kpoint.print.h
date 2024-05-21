#pragma once
#include <mpi.h>
#include "qumasun_kpoint.h"

inline
void QUMASUN_KPOINT::mPrintTime() {
    const bool is_root = IsRoot(m_mpi_comm);

    //for performance////////////////////////////////////////////////////////////
    auto GemmGFLOPS = [](int N, int M, int K, int64_t count, double time) {
        return (double)count * (double)N * (double)M * (double)K * 2.0e-9 / time;
        };

    double gemm_flops[3] = {
        GemmGFLOPS(2 * num_solution, 2 * num_solution, ml_grid.Size3D(), watch.GetCount(17), watch.GetTime(17)),
        GemmGFLOPS(ml_grid.Size3D(), num_solution, num_solution * 2, watch.GetCount(16), watch.GetTime(16)),
        GemmGFLOPS(num_solution, num_solution, num_solution, 4 * watch.GetCount(19), watch.GetTime(19))
    };
    double total_gemm_flops[3]{ 0.0 };
    MPI_Reduce(gemm_flops, total_gemm_flops, 3, MPI_DOUBLE, MPI_SUM, 0, m_mpi_comm);

    if (!is_root) return;

    printf("\nPeak performance in GEMM=====================\n");
    printf("  DGEMM1(x2HS)   : %f GFLOPS  (%f GFLOPS/root_proc)\n", total_gemm_flops[0], gemm_flops[0]);
    printf("  DGEMM2(x2x)    : %f GFLOPS  (%f GFLOPS/root_proc)\n", total_gemm_flops[1], gemm_flops[1]);
    if (watch.GetCount(19) > 0) {
        printf("  ZGEMM3(HS2HS)  : %f GFLOPS  (%f GFLOPS/root_proc)\n", total_gemm_flops[2], gemm_flops[2]);
    }
    else {
        printf("  ZGEMM3(HS2HS)  : no-called\n");
    }

    printf("\nCalculation time====================\n");
    const double total_tm = watch.Total();
    //const double calc_tm = watch.Total({ 2,3,4,6,10,11,12,13,14,15,16,17,18,19 });
    const double prepare_core_tm = watch.Total({ 1,21,22,23,24,25,26,27,28,29 });
    const double prepare_core_hr_tm = watch.Total({ 31,33,34,35 });
    const double init_tm = watch.Total({ 0,7 }) + prepare_core_tm + prepare_core_hr_tm;
    const double calc_tm = total_tm - init_tm;

    printf("Total time        : %f [s]\n", total_tm);
    watch.Print("mInitialState   :", 0);
    printf("mPrepareCore      : %f [s]\n", prepare_core_tm);
    //watch.Print("mPrepareCore    :", 1);
    watch.Print("--ChargeVlocal  :", 21);
    watch.Print("--GatherRho     :", 23);
    watch.Print("--PoissonVlocal :", 24);
    watch.Print("--ScatterVlocal :", 25);
    watch.Print("--SetPccCharge  :", 26);
    watch.Print("--GatherPCC     :", 27);
    watch.Print("--GetValence    :", 28);
    watch.Print("--EnnCorrection :", 29);

    printf("mPrepareCore(HR)  : %f [s]\n", prepare_core_hr_tm);
    watch.Print("--ChargeVlocal  :", 31);
    watch.Print("--GatherRho     :", 33);
    watch.Print("--PoissonVlocal :", 34);
    watch.Print("--ScatterVlocal :", 35);
    watch.Print("--DownConvVlocal:", 36);

    watch.Print("InitialDensity  :", 7);
    watch.Print("MoveNuclei      :", 20);

    printf("\nCalculation time: %f [s]\n", calc_tm);
    watch.Print("mSetOccupancy   :", 2);
    watch.Print("mSetDensity     :", 3);
    watch.Print("mSetPotential   :", 4);
    watch.Print("mGetTotalEnergy :", 6);
    watch.Print("mGetForce       :", 8);
    printf("EigenSolver       :--\n");
    //watch.Print("EigenSolver     :", 5);
    watch.Print("--K-operation   :", 11);
    watch.Print("--V-operation   :", 12);
    watch.Print("--PPNonlocal    :", 13);
    watch.Print("--GEMM1(x2HS)   :", 17);
    watch.Print("--GEMM2(x2x)    :", 16);
    watch.Print("--GEMM3(HS2HS)  :", 19);
    watch.Print("--LAPACK        :", 15);
    watch.Print("--MPI           :", 14);
    watch.Print("--others        :", 10);
    watch.Print("--wait imbalance:", 18);

}



inline
void QUMASUN_KPOINT::PrintCondition() {
    const bool is_root = IsRoot(m_mpi_comm);
    if (!is_root)return;

    printf("\nSimulation Condition====================\n");
    printf("Box width: %f, %f, %f [Bohr]\n", m_box_x, m_box_y, m_box_z);
    printf("         : %f, %f, %f [Ang]\n", m_box_x * Bohr_as_Ang, m_box_y * Bohr_as_Ang, m_box_z * Bohr_as_Ang);
    printf("Grid size: %d, %d, %d\n", m_global_grid.SizeX(), m_global_grid.SizeY(), m_global_grid.SizeZ());
    printf("Grid width(dx,dy,dz): %f, %f, %f [Bohr]\n", m_dx, m_dy, m_dz);
    printf("                    : %f, %f, %f [Ang]\n", m_dx * Bohr_as_Ang, m_dy * Bohr_as_Ang, m_dz * Bohr_as_Ang);

    double Ecut_x = (2.0 * M_PI / m_dx) * (2.0 * M_PI / m_dx) / 2.0;
    double Ecut_y = (2.0 * M_PI / m_dy) * (2.0 * M_PI / m_dy) / 2.0;
    double Ecut_z = (2.0 * M_PI / m_dz) * (2.0 * M_PI / m_dz) / 2.0;
    Ecut_x /= 4.0;  //half of k grid is negative//
    Ecut_y /= 4.0;  //half of k grid is negative//
    Ecut_z /= 4.0;  //half of k grid is negative//
    printf("Ecut = 0.5*(2pi/2dx)^2\n\n"
        "  (Cutoff energy corresponding to plane wave basis)\n");
    printf("Ecut(x): %f [Hartree] = %f [Ry] = %f [eV]\n", Ecut_x, Ecut_x * 2.0, Ecut_x * Hartree_as_eV);
    printf("Ecut(y): %f [Hartree] = %f [Ry] = %f [eV]\n", Ecut_y, Ecut_y * 2.0, Ecut_y * Hartree_as_eV);
    printf("Ecut(z): %f [Hartree] = %f [Ry] = %f [eV]\n", Ecut_z, Ecut_z * 2.0, Ecut_z * Hartree_as_eV);
    printf("Spin polarization: %s\n", (is_spin_on ? "on" : "off"));
    printf("k-point sample: %d, %d, %d\n", m_kpoint_sampling[0], m_kpoint_sampling[1], m_kpoint_sampling[2]);    
    printf("High-resolution (HR) ratio for core charge: %d, %d, %d\n", m_HR_ratio_x, m_HR_ratio_y, m_HR_ratio_z);
    printf("====================Simulation Condition\n\n");
    fflush(stdout);
}
