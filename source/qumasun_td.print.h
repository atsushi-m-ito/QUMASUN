#pragma once
#include <mpi.h>
#include "qumasun_td.h"

inline
void QUMASUN_TD::mPrintTime() {
    const bool is_root = IsRoot(m_mpi_comm);
    if (!is_root) return;


    printf("\nCalculation time====================\n");
    const double total_tm = watch.Total();
    //const double calc_tm = watch.Total({ 2,3,4,6,10,11,12,13,14,15,16,17,18,19 });
    const double prepare_core_tm = watch.Total({ 1,21,22,23,24,25,26,27,28,29 });
    const double prepare_core_hr_tm = watch.Total({ 31,33,34,35 });
    const double init_state_tm = watch.Total({ 0,40,41 }); 
    const double init_tm = watch.Total({ 0,40,41,7 });
    const double output_tm = watch.Total({ 42,43,44 });
    const double potential_tm = watch.Total({ 4,50,51,52,53,54,55,56,57 });
    const double force_tm = watch.Total({ 8,70,71,72,73,74,75,76 });
    const double energy_tm = watch.Total({ 6,60,61,62,63,64,65,66,6 });
    const double calc_tm = total_tm - init_tm - output_tm;

    printf("Total time      : %f [s]\n", total_tm);
    printf("mInitialState   : %f [s]\n", init_state_tm);
    watch.Print("--LoadState     :", 40);
    watch.Print("--ScatterState  :", 41);
    watch.Print("InitialDensity  :", 7);
    printf("\nCalculation time: %f [s]\n", calc_tm);
    printf("mPrepareCore    : %f [s]\n", prepare_core_tm);
    //watch.Print("mPrepareCore    :", 1);
    watch.Print("--ChargeVlocal  :", 21);
    watch.Print("--Scaling       :", 22);
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

    watch.Print("MoveNuclei      :", 20);
    watch.Print("mSetOccupancy   :", 2);
    watch.Print("mSetDensity     :", 3);
    printf("mSetPotential   : %f [s]\n", potential_tm);
    watch.Print("--PoissonVhart  :", 50);
    watch.Print("--UpConvert(HR) :", 51);
    watch.Print("--IFFT(HR)      :", 52);
    watch.Print("--Normalize(HR) :", 53);
    watch.Print("--Check(HR)     :", 54);
    watch.Print("--XC            :", 55);
    watch.Print("--SumTotal      :", 56);
    watch.Print("--Scatter       :", 57);
    
    printf("mGetTotalEnergy : %f [s]\n", energy_tm);
    watch.Print("--Kinetic       :", 60);
    watch.Print("--Nonlocal      :", 61);
    watch.Print("--InnerVRho     :", 62);
    watch.Print("--InnerVRho(HR) :", 63);
    watch.Print("--XC            :", 64);
    watch.Print("--VhartAtPoint  :", 65);
    watch.Print("--CoreCoreDirect:", 66);

    printf("mGetForce       : %f [s]\n", force_tm);
    watch.Print("--Nonlocal(Diff):", 75); 
    watch.Print("--Nonlocal(Inner:", 76);
    watch.Print("--Nonlocal(other:", 70);
    watch.Print("--DiffV         :", 71);
    watch.Print("--InnerNuclRho  :", 72);
    watch.Print("--DiffV(HR)     :", 73);
    watch.Print("--InnerNucl(HR) :", 74);

    printf("TimeEvolusionState     :--\n");
    //watch.Print("EigenSolver     :", 5);
    watch.Print("--K-operation   :", 11);
    watch.Print("--V-operation   :", 12);
    watch.Print("--PPNonlocal    :", 13);
    watch.Print("--others        :", 10);

    printf("Output                 :--\n");
    watch.Print("--Gather        :", 43);
    watch.Print("--WriteCube     :", 44);
    watch.Print("--others        :", 42);
}



inline
void QUMASUN_TD::PrintCondition() {
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


inline
void QUMASUN_TD::PrintTime() {
    mPrintTime();
}
