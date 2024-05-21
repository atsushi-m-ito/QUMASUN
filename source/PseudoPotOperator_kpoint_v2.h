#ifdef USE_MPI
#pragma once
#include <mpi.h>
#include "mpi_helper.h"
#include "PseudoPotOperator_vlocal.h"
#include "soacomplex.h"
#include "SetBlockYlm.h"
#include "w_dgemm.h"
#include "StopWatch.h"

//fast calculation without imaginary part of exp(ikx) for Bloch theorem when gamma point.
#define FAST_GAMMA_POINT 


#define EXPAND_PROJECTOR_YLM_Z
//#define USE_COMPLEX_PROJECTOR  //complex type slower than SoAComplex



class PseudoPotIntegrator_kpoint : public PseudoPotIntegrator_vlocal {
private:

    //Yml * exp(ikx), where k is k-sampling point//
    struct PP_Proj_Ylm_each_atom_kpoint {
        int max_l;   //max of angular momentum l.
        int num_projectors;
        std::vector < std::vector<SoAComplex>> projector_expikx; //block化されたprojector, nonlocal項の数だけ存在//
        std::vector<double*> Ylm; //block化されたYlm, (max_l+1)^2の数だけ存在//
#ifdef USE_COMPLEX_PROJECTOR
        std::vector < std::vector<OneComplex*>> c_projector_expikx; //block化されたprojector, nonlocal項の数だけ存在//
#endif
        int* quantum_l = nullptr;   //projectorごとのangular momentum l, nonlocal項の数だけ存在//
        double* projector_energy = nullptr;   //projectorごとのangular momentum l, nonlocal項の数だけ存在//


#ifdef EXPAND_PROJECTOR_YLM_Z
#ifdef USE_COMPLEX_PROJECTOR
        std::vector < OneComplex*> projector_RY_expikx; //block化されたprojector(Ylmも作用済み), nonlocal項の数だけ存在, Lの種類ごとにまとまる//
#else
        std::vector < double*> projector_RY_expikx; //block化されたprojector(Ylmも作用済み), nonlocal項の数だけ存在, Lの種類ごとにまとまる//
#endif
        int head_projector_l[5]{ 0 };  //4 is 1 + maximum quntam number l;
#endif


        ~PP_Proj_Ylm_each_atom_kpoint() {
            for (auto&& projector : projector_expikx) {
                for (auto&& p : projector) {
                    delete[] p.re;
                    //delete[] p.im;
                }
            }
            for (auto&& y : Ylm) {
                delete[] y;
            }
#ifdef USE_COMPLEX_PROJECTOR
            for (auto&& projector : c_projector_expikx) {
                for (auto&& p : projector) {
                    delete[] p;
                }
            }

#endif
#ifdef EXPAND_PROJECTOR_YLM_Z
            for (auto&& p : projector_RY_expikx) {
                delete[] p;
            }
#endif
        }
    };


    std::vector< PP_Proj_Ylm_each_atom_kpoint >  m_nonlocal_projectors;
    int m_num_all_nonlocal = 0;
    std::vector<int> m_num_list_nonlocal;

#ifdef FAST_GAMMA_POINT
    std::vector<bool> m_is_gamma_point_list;
#endif

#ifdef TIME_PP_MPI
    StopWatch<10, true> watch_pp;
#else
    StopWatch<10, false> watch_pp;
#endif


public:
    ~PseudoPotIntegrator_kpoint() {
#ifdef TIME_PP_MPI
        if (is_root) {
            mPrintTime();
        }
#endif
    }


    void mSetNonlocalProjectorInfo(PP_Proj_Ylm_each_atom_kpoint& nonlocal_projectors,
        const Nucleus nucleus) {

        const PseudoPot_MBK* pp = mFindPseudoPot(nucleus.Z);
        if (pp == nullptr) {
            printf("ERROR: VPS is not loaded: Z = %d\n", nucleus.Z);
            return;
        }

        const int num_projectors = pp->num_projectors;
        nonlocal_projectors.num_projectors = num_projectors;

        const int MAX_L = pp->MaxProjectorL();
        nonlocal_projectors.max_l = MAX_L;
        nonlocal_projectors.quantum_l = pp->projector_quantum_l;
        nonlocal_projectors.projector_energy = pp->projector_energy_up;
    }

#ifdef FAST_GAMMA_POINT
    void mCreateProjectorRadialOnBlock_gamma(NonlocalBlocks& nonlocal_blocks,
        PP_Proj_Ylm_each_atom_kpoint& nonlocal_projectors,
        const Nucleus nucleus, const GridRange& l_grid,
        int id_spin_kpoint) {


        if (nonlocal_projectors.num_projectors == 0) {
            return;
        }

        const int MAX_L = nonlocal_projectors.max_l;


        const int num_projectors = nonlocal_projectors.num_projectors;
        auto& projector = nonlocal_projectors.projector_expikx.emplace_back();
        projector.resize(num_projectors);

        const PseudoPot_MBK* pp = mFindPseudoPot(nucleus.Z);

        for (int k = 0; k < num_projectors; ++k) {
            const int l = nonlocal_projectors.quantum_l[k];
            const size_t total_block_size = nonlocal_blocks[l].grid_sizes;

            auto& proj_block = projector[k];
            proj_block.re = new double[total_block_size * 2];
            proj_block.im = proj_block.re + total_block_size;

            //auto& expikx_l = expikx[l];

            const double rr_min = (pp->radius[0]) * (pp->radius[0]);
            const double rr_max = (pp->cutoff_r[l]) * (pp->cutoff_r[l]);
            const auto& range_blocks = nonlocal_blocks[l].range_blocks;
            const auto& shifted_grid = nonlocal_blocks[l].shifted_grids;
            size_t num_blocks = range_blocks.size();
            int offset_i = 0;
            for (size_t ib = 0; ib < num_blocks; ++ib) {
                const int ix_begin = range_blocks[ib].begin_x;
                const int iy_begin = range_blocks[ib].begin_y;
                const int iz_begin = range_blocks[ib].begin_z;
                const int ix_end = range_blocks[ib].end_x;
                const int iy_end = range_blocks[ib].end_y;
                const int iz_end = range_blocks[ib].end_z;

                const int size_x = ix_end - ix_begin;
                const int size_y = iy_end - iy_begin;

                for (int iz = iz_begin; iz < iz_end; ++iz) {
                    const double z = m_dz * (double)(iz - shifted_grid[ib].z) - nucleus.Rz;
                    const double zz = z * z;
                    for (int iy = iy_begin; iy < iy_end; ++iy) {
                        const double y = m_dy * (double)(iy - shifted_grid[ib].y) - nucleus.Ry;
                        const double yy_zz = y * y + zz;
                        for (int ix = ix_begin; ix < ix_end; ++ix) {
                            const double x = m_dx * (double)(ix - shifted_grid[ib].x) - nucleus.Rx;
                            const double rr = x * x + yy_zz;
                            const int i = offset_i + (ix - ix_begin) + size_x * ((iy - iy_begin) + size_y * (iz - iz_begin));


                            if (rr_min > rr) {
                                const double proj_spin_orbit_up = pp->projector[k * 2][0];
                                const double proj_spin_orbit_dn = pp->projector[k * 2 + 1][0];
                                const double p = (proj_spin_orbit_up + proj_spin_orbit_dn) / 2.0;
                                proj_block.re[i] = p;
                                //proj_block.im[i] = 0.0;

                            } else if (rr_max < rr) {
                                proj_block.re[i] = 0.0;
                                //proj_block.im[i] = 0.0;
                            } else {
                                //NOTE: if it is not relative DFT with spin-orbit interaction, projectors (j+1/2) and (j-1/2) should be averaged.
                                const double proj_spin_orbit_up = RadialGrid2::GetValueBySquare(rr, pp->projector[k * 2], pp->num_radial_grids, pp->xi_min, pp->xi_delta);
                                const double proj_spin_orbit_dn = RadialGrid2::GetValueBySquare(rr, pp->projector[k * 2 + 1], pp->num_radial_grids, pp->xi_min, pp->xi_delta);
                                const double p = (proj_spin_orbit_up + proj_spin_orbit_dn) / 2.0;
                                proj_block.re[i] = p;
                                //proj_block.im[i] = 0.0;

                            }

                        }
                    }
                }
                offset_i += range_blocks[ib].Size3D();

            }

            for (int i = 0; i < offset_i; ++i) {
                proj_block.im[i] = 0.0;
            }
        }

#ifdef USE_COMPLEX_PROJECTOR
        {
            //const int num_projectors = nonlocal_projectors.num_projectors;
            auto& projector = nonlocal_projectors.c_projector_expikx.emplace_back();
            projector.resize(num_projectors);

            const PseudoPot_MBK* pp = mFindPseudoPot(nucleus.Z);

            for (int k = 0; k < num_projectors; ++k) {
                const int l = nonlocal_projectors.quantum_l[k];
                const size_t total_block_size = nonlocal_blocks[l].grid_sizes;

                auto& proj_block = projector[k];
                proj_block = new OneComplex[total_block_size];

                //auto& expikx_l = expikx[l];

                const double rr_min = (pp->radius[0]) * (pp->radius[0]);
                const double rr_max = (pp->cutoff_r[l]) * (pp->cutoff_r[l]);
                const auto& range_blocks = nonlocal_blocks[l].range_blocks;
                const auto& shifted_grid = nonlocal_blocks[l].shifted_grids;
                size_t num_blocks = range_blocks.size();
                int offset_i = 0;
                for (size_t ib = 0; ib < num_blocks; ++ib) {
                    const int ix_begin = range_blocks[ib].begin_x;
                    const int iy_begin = range_blocks[ib].begin_y;
                    const int iz_begin = range_blocks[ib].begin_z;
                    const int ix_end = range_blocks[ib].end_x;
                    const int iy_end = range_blocks[ib].end_y;
                    const int iz_end = range_blocks[ib].end_z;

                    const int size_x = ix_end - ix_begin;
                    const int size_y = iy_end - iy_begin;

                    for (int iz = iz_begin; iz < iz_end; ++iz) {
                        const double z = m_dz * (double)(iz - shifted_grid[ib].z) - nucleus.Rz;
                        const double zz = z * z;
                        for (int iy = iy_begin; iy < iy_end; ++iy) {
                            const double y = m_dy * (double)(iy - shifted_grid[ib].y) - nucleus.Ry;
                            const double yy_zz = y * y + zz;
                            for (int ix = ix_begin; ix < ix_end; ++ix) {
                                const double x = m_dx * (double)(ix - shifted_grid[ib].x) - nucleus.Rx;
                                const double rr = x * x + yy_zz;
                                const int i = offset_i + (ix - ix_begin) + size_x * ((iy - iy_begin) + size_y * (iz - iz_begin));


                                if (rr_min > rr) {
                                    const double proj_spin_orbit_up = pp->projector[k * 2][0];
                                    const double proj_spin_orbit_dn = pp->projector[k * 2 + 1][0];
                                    const double p = (proj_spin_orbit_up + proj_spin_orbit_dn) / 2.0;
                                    proj_block[i].r = p;
                                    proj_block[i].i = 0.0;

                                } else if (rr_max < rr) {
                                    proj_block[i].r = 0.0;
                                    proj_block[i].i = 0.0;
                                } else {
                                    //NOTE: if it is not relative DFT with spin-orbit interaction, projectors (j+1/2) and (j-1/2) should be averaged.
                                    const double proj_spin_orbit_up = RadialGrid2::GetValueBySquare(rr, pp->projector[k * 2], pp->num_radial_grids, pp->xi_min, pp->xi_delta);
                                    const double proj_spin_orbit_dn = RadialGrid2::GetValueBySquare(rr, pp->projector[k * 2 + 1], pp->num_radial_grids, pp->xi_min, pp->xi_delta);
                                    const double p = (proj_spin_orbit_up + proj_spin_orbit_dn) / 2.0;
                                    proj_block[i].r = p;
                                    proj_block[i].i = 0.0;

                                }

                            }
                        }
                    }
                    offset_i += range_blocks[ib].Size3D();

                }
            }
        }
#endif


    }
#endif

    void mCreateProjectorRadialOnBlock(NonlocalBlocks& nonlocal_blocks,
        PP_Proj_Ylm_each_atom_kpoint& nonlocal_projectors,
        const Nucleus nucleus, const GridRange& l_grid,
        int id_spin_kpoint, double gx_dx, double gy_dy, double gz_dz) {


        if (nonlocal_projectors.num_projectors == 0) {
            return;
        }

        const int MAX_L = nonlocal_projectors.max_l;
        //auto& range_block_in_1 = nonlocal_block.range_blocks;
        //auto& shifted_grids = nonlocal_block.shifted_grids;
        //auto& grid_sizes = nonlocal_block.grid_sizes;
        //auto& blocks = nonlocal_block.blocks;


        //prepare exp(ikx)//////////////////////
        std::vector<SoAComplex> expikx(MAX_L + 1);
        for (int l = 0; l <= MAX_L; ++l) {
            const auto& range_blocks = nonlocal_blocks[l].range_blocks;
            const auto& shifted_grid = nonlocal_blocks[l].shifted_grids;
            size_t num_blocks = range_blocks.size();
            const size_t total_block_size = nonlocal_blocks[l].grid_sizes;

            auto& expikx_l = expikx[l];
            expikx_l.re = new double[total_block_size * 2];
            expikx_l.im = expikx[l].re + total_block_size;

            int offset_i = 0;
            for (size_t ib = 0; ib < num_blocks; ++ib) {
                const int ix_begin = range_blocks[ib].begin_x;
                const int iy_begin = range_blocks[ib].begin_y;
                const int iz_begin = range_blocks[ib].begin_z;
                const int ix_end = range_blocks[ib].end_x;
                const int iy_end = range_blocks[ib].end_y;
                const int iz_end = range_blocks[ib].end_z;

                const int size_x = ix_end - ix_begin;
                const int size_y = iy_end - iy_begin;

                for (int iz = iz_begin; iz < iz_end; ++iz) {
                    const double kzz = gz_dz * (double)(iz - shifted_grid[ib].z);
                    for (int iy = iy_begin; iy < iy_end; ++iy) {
                        const double kyy = gy_dy * (double)(iy - shifted_grid[ib].y);
                        for (int ix = ix_begin; ix < ix_end; ++ix) {
                            const double kxx = gx_dx * (double)(ix - shifted_grid[ib].x);
                            const int i = offset_i + (ix - ix_begin) + size_x * ((iy - iy_begin) + size_y * (iz - iz_begin));

                            double cosikx = cos(kxx + kyy + kzz);
                            double sinikx = sin(kxx + kyy + kzz);

                            expikx_l.re[i] = cosikx;
                            expikx_l.im[i] = sinikx;

                        }
                    }
                }
                offset_i += range_blocks[ib].Size3D();

            }
        }


        const int num_projectors = nonlocal_projectors.num_projectors;
        auto& projector = nonlocal_projectors.projector_expikx.emplace_back();
        projector.resize(num_projectors);

        const PseudoPot_MBK* pp = mFindPseudoPot(nucleus.Z);

        for (int k = 0; k < num_projectors; ++k) {
            const int l = nonlocal_projectors.quantum_l[k];
            const size_t total_block_size = nonlocal_blocks[l].grid_sizes;

            auto& proj_block = projector[k];
            proj_block.re = new double[total_block_size * 2];
            proj_block.im = proj_block.re + total_block_size;

            auto& expikx_l = expikx[l];

            const double rr_min = (pp->radius[0]) * (pp->radius[0]);
            const double rr_max = (pp->cutoff_r[l]) * (pp->cutoff_r[l]);
            const auto& range_blocks = nonlocal_blocks[l].range_blocks;
            const auto& shifted_grid = nonlocal_blocks[l].shifted_grids;
            size_t num_blocks = range_blocks.size();
            int offset_i = 0;
            for (size_t ib = 0; ib < num_blocks; ++ib) {
                const int ix_begin = range_blocks[ib].begin_x;
                const int iy_begin = range_blocks[ib].begin_y;
                const int iz_begin = range_blocks[ib].begin_z;
                const int ix_end = range_blocks[ib].end_x;
                const int iy_end = range_blocks[ib].end_y;
                const int iz_end = range_blocks[ib].end_z;

                const int size_x = ix_end - ix_begin;
                const int size_y = iy_end - iy_begin;

                for (int iz = iz_begin; iz < iz_end; ++iz) {
                    const double z = m_dz * (double)(iz - shifted_grid[ib].z) - nucleus.Rz;
                    const double zz = z * z;
                    for (int iy = iy_begin; iy < iy_end; ++iy) {
                        const double y = m_dy * (double)(iy - shifted_grid[ib].y) - nucleus.Ry;
                        const double yy_zz = y * y + zz;
                        for (int ix = ix_begin; ix < ix_end; ++ix) {
                            const double x = m_dx * (double)(ix - shifted_grid[ib].x) - nucleus.Rx;
                            const double rr = x * x + yy_zz;
                            const int i = offset_i + (ix - ix_begin) + size_x * ((iy - iy_begin) + size_y * (iz - iz_begin));


                            if (rr_min > rr) {
                                const double proj_spin_orbit_up = pp->projector[k * 2][0];
                                const double proj_spin_orbit_dn = pp->projector[k * 2 + 1][0];
                                const double p = (proj_spin_orbit_up + proj_spin_orbit_dn) / 2.0;
                                proj_block.re[i] = p * expikx_l.re[i];
                                proj_block.im[i] = p * expikx_l.im[i];

                            } else if (rr_max < rr) {
                                proj_block.re[i] = 0.0;
                                proj_block.im[i] = 0.0;
                            } else {
                                //NOTE: if it is not relative DFT with spin-orbit interaction, projectors (j+1/2) and (j-1/2) should be averaged.
                                const double proj_spin_orbit_up = RadialGrid2::GetValueBySquare(rr, pp->projector[k * 2], pp->num_radial_grids, pp->xi_min, pp->xi_delta);
                                const double proj_spin_orbit_dn = RadialGrid2::GetValueBySquare(rr, pp->projector[k * 2 + 1], pp->num_radial_grids, pp->xi_min, pp->xi_delta);
                                const double p = (proj_spin_orbit_up + proj_spin_orbit_dn) / 2.0;
                                proj_block.re[i] = p * expikx_l.re[i];
                                proj_block.im[i] = p * expikx_l.im[i];

                            }

                        }
                    }
                }
                offset_i += range_blocks[ib].Size3D();

            }
        }

#ifdef USE_COMPLEX_PROJECTOR
        {
            //const int num_projectors = nonlocal_projectors.num_projectors;
            auto& projector = nonlocal_projectors.c_projector_expikx.emplace_back();
            projector.resize(num_projectors);

            const PseudoPot_MBK* pp = mFindPseudoPot(nucleus.Z);

            for (int k = 0; k < num_projectors; ++k) {
                const int l = nonlocal_projectors.quantum_l[k];
                const size_t total_block_size = nonlocal_blocks[l].grid_sizes;

                auto& proj_block = projector[k];
                proj_block = new OneComplex[total_block_size];

                auto& expikx_l = expikx[l];

                const double rr_min = (pp->radius[0]) * (pp->radius[0]);
                const double rr_max = (pp->cutoff_r[l]) * (pp->cutoff_r[l]);
                const auto& range_blocks = nonlocal_blocks[l].range_blocks;
                const auto& shifted_grid = nonlocal_blocks[l].shifted_grids;
                size_t num_blocks = range_blocks.size();
                int offset_i = 0;
                for (size_t ib = 0; ib < num_blocks; ++ib) {
                    const int ix_begin = range_blocks[ib].begin_x;
                    const int iy_begin = range_blocks[ib].begin_y;
                    const int iz_begin = range_blocks[ib].begin_z;
                    const int ix_end = range_blocks[ib].end_x;
                    const int iy_end = range_blocks[ib].end_y;
                    const int iz_end = range_blocks[ib].end_z;

                    const int size_x = ix_end - ix_begin;
                    const int size_y = iy_end - iy_begin;

                    for (int iz = iz_begin; iz < iz_end; ++iz) {
                        const double z = m_dz * (double)(iz - shifted_grid[ib].z) - nucleus.Rz;
                        const double zz = z * z;
                        for (int iy = iy_begin; iy < iy_end; ++iy) {
                            const double y = m_dy * (double)(iy - shifted_grid[ib].y) - nucleus.Ry;
                            const double yy_zz = y * y + zz;
                            for (int ix = ix_begin; ix < ix_end; ++ix) {
                                const double x = m_dx * (double)(ix - shifted_grid[ib].x) - nucleus.Rx;
                                const double rr = x * x + yy_zz;
                                const int i = offset_i + (ix - ix_begin) + size_x * ((iy - iy_begin) + size_y * (iz - iz_begin));


                                if (rr_min > rr) {
                                    const double proj_spin_orbit_up = pp->projector[k * 2][0];
                                    const double proj_spin_orbit_dn = pp->projector[k * 2 + 1][0];
                                    const double p = (proj_spin_orbit_up + proj_spin_orbit_dn) / 2.0;
                                    proj_block[i].r = p * expikx_l.re[i];
                                    proj_block[i].i = p * expikx_l.im[i];

                                } else if (rr_max < rr) {
                                    proj_block[i].r = 0.0;
                                    proj_block[i].i = 0.0;
                                } else {
                                    //NOTE: if it is not relative DFT with spin-orbit interaction, projectors (j+1/2) and (j-1/2) should be averaged.
                                    const double proj_spin_orbit_up = RadialGrid2::GetValueBySquare(rr, pp->projector[k * 2], pp->num_radial_grids, pp->xi_min, pp->xi_delta);
                                    const double proj_spin_orbit_dn = RadialGrid2::GetValueBySquare(rr, pp->projector[k * 2 + 1], pp->num_radial_grids, pp->xi_min, pp->xi_delta);
                                    const double p = (proj_spin_orbit_up + proj_spin_orbit_dn) / 2.0;
                                    proj_block[i].r = p * expikx_l.re[i];
                                    proj_block[i].i = p * expikx_l.im[i];

                                }

                            }
                        }
                    }
                    offset_i += range_blocks[ib].Size3D();

                }
            }
        }
#endif

        //delete memory//
        for (auto& a : expikx) {
            delete[] a.re;
        }

    }

    void CountNonlocalProjector(const Nucleus* nuclei, int num_nuclei){
        
        int num_all_nonlocal = 0;
        m_num_list_nonlocal.clear();
        m_num_list_nonlocal.push_back(0);
        for (int ni = 0; ni < num_nuclei; ++ni) {
            const PseudoPot_MBK* pp = mFindPseudoPot(nuclei[ni].Z);
            if (pp == nullptr) {
                printf("ERROR: VPS is not loaded: Z = %d\n", nuclei[ni].Z);
                return;
            }

            num_all_nonlocal += pp->TotalProjectorLM();
            m_num_list_nonlocal.push_back(num_all_nonlocal);
        }
        
        m_num_all_nonlocal = num_all_nonlocal;
    }

#ifdef EXPAND_PROJECTOR_YLM_Z
    void mExpandProjectorYlm(NonlocalBlocks& nonlocal_blocks,
        PP_Proj_Ylm_each_atom_kpoint& nonlocal_projectors, int id_spin_kpoint) {

        const double Y_00 = 1.0 / sqrt(4.0 * M_PI);

        const int num_projectors = nonlocal_projectors.num_projectors;
        auto& head_l = nonlocal_projectors.head_projector_l;
        {
            int l = 0;
            //int num_l = 0;
            head_l[0] = 0;
            for (int k = 0; k < num_projectors; ++k) {
                while (l < nonlocal_projectors.quantum_l[k]) {
                    head_l[l + 1] = k;
                    ++l;
                }
            }
            while (l < 4) {
                head_l[l + 1] = num_projectors;
                ++l;
            }
        }

        {
            int prev = 0;
            for (int l = 0; l < 4; ++l) {
                int p = head_l[l + 1];
                head_l[l + 1] = (p - prev) * (2 * l + 1) + head_l[l];
                prev = p;
            }
        }

        int64_t total_buf_size = 0;
        for (int l = 0; l <= nonlocal_projectors.max_l; ++l) {
            total_buf_size += nonlocal_blocks[l].grid_sizes * (head_l[l + 1] - head_l[l]);// *(2 * l + 1);
        }

        auto& projector_RY_expikx = nonlocal_projectors.projector_RY_expikx;
        auto*& buffer = projector_RY_expikx.emplace_back(nullptr);
#ifdef USE_COMPLEX_PROJECTOR
        buffer = new OneComplex[total_buf_size];
#else
        buffer = new double [total_buf_size*2];
#endif

        int64_t offset = 0;

#ifdef FAST_GAMMA_POINT
        if(m_is_gamma_point_list[id_spin_kpoint])
        {
            for (int k = 0; k < nonlocal_projectors.num_projectors; ++k) {
                const int l = nonlocal_projectors.quantum_l[k];
                //total_buf_size += nonlocal_blocks[l].grid_sizes * (head_l[l + 1] - head_l[l]);


                const auto* proj_re = nonlocal_projectors.projector_expikx[id_spin_kpoint][k].re;
                //const auto* proj_im = nonlocal_projectors.projector_expikx[id_spin_kpoint][k].im;
                const int i_end = nonlocal_blocks[l].grid_sizes;
                if (l == 0) {
#ifdef USE_COMPLEX_PROJECTOR
                    for (int i = 0; i < i_end; ++i) {
                        buffer[offset + i].r = Y_00 * proj_re[i];
                        buffer[offset + i].i = 0.0;
                    }
                    offset += i_end;

#else
                    for (int i = 0; i < i_end; ++i) {
                        buffer[offset + i] = Y_00 * proj_re[i];
                    }
                    offset += i_end;

#endif
                } else {
                    for (int m = 0; m <= 2 * l; ++m) {
                        const auto* Ylm = nonlocal_projectors.Ylm[l * l + m];
#ifdef USE_COMPLEX_PROJECTOR
                        for (int i = 0; i < i_end; ++i) {
                            buffer[offset + i].r = Ylm[i] * proj_re[i];
                            buffer[offset + i].i = 0.0;
                        }
                        offset += i_end;

#else
                        for (int i = 0; i < i_end; ++i) {
                            buffer[offset + i] = Ylm[i] * proj_re[i];
                        }
                        offset += i_end;

#endif
                    }
                }
            }
        }
        else
#endif
        {
            for (int k = 0; k < nonlocal_projectors.num_projectors; ++k) {
                const int l = nonlocal_projectors.quantum_l[k];
                //total_buf_size += nonlocal_blocks[l].grid_sizes * (head_l[l + 1] - head_l[l]);


                const auto* proj_re = nonlocal_projectors.projector_expikx[id_spin_kpoint][k].re;
                const auto* proj_im = nonlocal_projectors.projector_expikx[id_spin_kpoint][k].im;
                const int i_end = nonlocal_blocks[l].grid_sizes;
                if (l == 0) {
#ifdef USE_COMPLEX_PROJECTOR
                    for (int i = 0; i < i_end; ++i) {
                        buffer[offset + i].r = Y_00 * proj_re[i];
                        buffer[offset + i].i = Y_00 * proj_im[i];
                    }
                    offset += i_end;

#else
                    for (int i = 0; i < i_end; ++i) {
                        buffer[offset + i] = Y_00 * proj_re[i];
                    }
                    offset += i_end;
                    for (int i = 0; i < i_end; ++i) {
                        buffer[offset + i] = Y_00 * proj_im[i];
                    }
                    offset += i_end;
#endif
                } else {
                    for (int m = 0; m <= 2 * l; ++m) {
                        const auto* Ylm = nonlocal_projectors.Ylm[l * l + m];
#ifdef USE_COMPLEX_PROJECTOR
                        for (int i = 0; i < i_end; ++i) {
                            buffer[offset + i].r = Ylm[i] * proj_re[i];
                            buffer[offset + i].i = Ylm[i] * proj_im[i];
                        }
                        offset += i_end;

#else
                        for (int i = 0; i < i_end; ++i) {
                            buffer[offset + i] = Ylm[i] * proj_re[i];
                        }
                        offset += i_end;
                        for (int i = 0; i < i_end; ++i) {
                            buffer[offset + i] = Ylm[i] * proj_im[i];
                        }
                        offset += i_end;
#endif
                    }
                }
            }
        }


    }
#endif

    /*
    * 擬ポテンシャルを原子核位置に合わせて必要な情報を用意.
    */
public:
    void UpdatePosition(const Nucleus* nuclei, int num_nuclei, const GridRangeMPI& l_grid) {

        //base-class//
        UpdatePositionLocal(nuclei, num_nuclei, l_grid);

        CountNonlocalProjector(nuclei, num_nuclei);

        m_nonlocal_projectors.clear();
        m_nonlocal_projectors.resize(num_nuclei);

        for (int ni = 0; ni < num_nuclei; ++ni) {
            //set m_vlocal_block and m_nonlocal_block//
            mSetNonlocalProjectorInfo(m_nonlocal_projectors[ni], nuclei[ni]);
            tCreateProjectorYlmOnBlock(m_nonlocal_blocks[ni], m_nonlocal_projectors[ni], nuclei[ni], m_dx, m_dy, m_dz, l_grid);
        }
    }

    /*
    * 擬ポテンシャルのnonlocal項を原子核位置に合わせて切り取ったsubグリッド領域に用意する.
    */
    void UpdateNonlocalBlock(const Nucleus* nuclei, int num_nuclei, const GridRangeMPI& l_grid,
        int id_spin_kpoint, double gx_dx, double gy_dy, double gz_dz) {

#ifdef FAST_GAMMA_POINT
        if (m_is_gamma_point_list.size() <= id_spin_kpoint) {
            m_is_gamma_point_list.resize(id_spin_kpoint + 1);
        }
#endif

#ifdef FAST_GAMMA_POINT
        if ((gx_dx == 0.0) && (gy_dy == 0.0) && (gz_dz == 0.0)) {

            m_is_gamma_point_list[id_spin_kpoint] = true;

            for (int ni = 0; ni < num_nuclei; ++ni) {
                //fast calculation without imaginary part of exp(ikx) for Bloch theorem//
                mCreateProjectorRadialOnBlock_gamma(m_nonlocal_blocks[ni], m_nonlocal_projectors[ni], nuclei[ni], l_grid,
                    id_spin_kpoint);

#ifdef EXPAND_PROJECTOR_YLM_Z
                mExpandProjectorYlm(m_nonlocal_blocks[ni], m_nonlocal_projectors[ni], id_spin_kpoint);
#endif
            }
        } 
        else 
        {
            m_is_gamma_point_list[id_spin_kpoint] = false;

            for (int ni = 0; ni < num_nuclei; ++ni) {
                mCreateProjectorRadialOnBlock(m_nonlocal_blocks[ni], m_nonlocal_projectors[ni], nuclei[ni], l_grid,
                    id_spin_kpoint, gx_dx, gy_dy, gz_dz);

#ifdef EXPAND_PROJECTOR_YLM_Z
                mExpandProjectorYlm(m_nonlocal_blocks[ni], m_nonlocal_projectors[ni], id_spin_kpoint);
#endif
            }
        }
#else

        for (int ni = 0; ni < num_nuclei; ++ni) {
            mCreateProjectorRadialOnBlock(m_nonlocal_blocks[ni], m_nonlocal_projectors[ni], nuclei[ni], l_grid,
                id_spin_kpoint, gx_dx, gy_dy, gz_dz);

#ifdef EXPAND_PROJECTOR_YLM_Z
            mExpandProjectorYlm(m_nonlocal_blocks[ni], m_nonlocal_projectors[ni], id_spin_kpoint);
#endif
        }

#endif
    }


    void mInnerNonlocalPsi(const SoAComplex& l_psi, const GridRangeMPI& l_grid, const Nucleus* nuclei, int num_nuclei, OneComplex* l_inner, int id_spin_kpoint, int num_all_nonlocal, double dV) {


        const double Y_00 = 1.0 / sqrt(4.0 * M_PI);


        int index = 0;
        for (int ni = 0; ni < num_nuclei; ++ni) {
            const int max_l = m_nonlocal_projectors[ni].max_l;
            std::vector<SoAComplex> cut_psi(max_l + 1);

            const auto& blocks = m_nonlocal_blocks[ni];

            for (int l = 0; l <= max_l; ++l) {

                cut_psi[l].re = new double[blocks[l].grid_sizes * 2];
                cut_psi[l].im = cut_psi[l].re + blocks[l].grid_sizes;

                CutSubgridByRanges(blocks[l].range_blocks, cut_psi[l].re, l_grid, l_psi.re);
                CutSubgridByRanges(blocks[l].range_blocks, cut_psi[l].im, l_grid, l_psi.im);
            }


            for (int k = 0; k < m_nonlocal_projectors[ni].num_projectors; ++k) {
                const int l = m_nonlocal_projectors[ni].quantum_l[k];
                const auto* proj_re = m_nonlocal_projectors[ni].projector_expikx[id_spin_kpoint][k].re;
                const auto* proj_im = m_nonlocal_projectors[ni].projector_expikx[id_spin_kpoint][k].im;
                const auto* psi_re = cut_psi[l].re;
                const auto* psi_im = cut_psi[l].im;

                if (l == 0) {
                    /*
                    double inner = ForReduce2(cut_grid[l], 0.0,
                    [&proj, &l_psi](int64_t i) { return proj[i] * l_psi[i]; });
                    */
                    size_t nl_size = blocks[l].grid_sizes;
                    double inner_re = 0.0;
                    double inner_im = 0.0;
                    for (int i = 0; i < nl_size; ++i) {
                        inner_re += proj_re[i] * psi_re[i] - proj_im[i] * psi_im[i];
                        inner_im += proj_im[i] * psi_re[i] + proj_re[i] * psi_im[i];
                    }
                    l_inner[index].r = inner_re * dV * Y_00;
                    l_inner[index].i = inner_im * dV * Y_00;
                    ++index;

                } else {
                    for (int m1 = 0; m1 < 2 * l + 1; ++m1) {
                        const auto* Ylm = m_nonlocal_projectors[ni].Ylm[l * l + m1];
                        //double inner = ForReduce2(cut_grid[l], 0.0,
                        //	[&proj, &Ylm, &l_psi](int64_t i) { return proj[i] * Ylm[i] * l_psi[i]; });
                        size_t nl_size = blocks[l].grid_sizes;
                        double inner_re = 0.0;
                        double inner_im = 0.0;
                        for (int i = 0; i < nl_size; ++i) {
                            inner_re += Ylm[i] * (proj_re[i] * psi_re[i] - proj_im[i] * psi_im[i]);
                            inner_im += Ylm[i] * (proj_im[i] * psi_re[i] + proj_re[i] * psi_im[i]);
                        }
                        l_inner[index].r = inner_re * dV;
                        l_inner[index].i = inner_im * dV;
                        ++index;
                    }
                }
            }

            for (int l = 0; l <= max_l; ++l) {
                delete[] cut_psi[l].re;
            }


        }


        if (num_all_nonlocal != index) {
            printf("ERROR: number of nonlocal is invalid"); fflush(stdout);
        }
        //printf("[%d] ProjectionPP 2: 0x%zx\n", GetProcessID(l_grid.mpi_comm), l_inner); fflush(stdout);
        //MPI_Barrier(l_grid.mpi_comm);
    }


    void mInnerNonlocalPsi_v2(const SoAComplex& l_psi, const GridRangeMPI& l_grid, const Nucleus* nuclei, int num_nuclei, OneComplex* l_inner, int id_spin_kpoint, std::vector<int>& info_range, double dV) {

        //const auto& nonlocal_block = m_nonlocal_block[id_spin_kpoint];

        const double Y_00 = 1.0 / sqrt(4.0 * M_PI);

        const int num_all_nonlocal = info_range[num_nuclei];
        OneComplex* sum_inner = l_inner + num_all_nonlocal;

        std::vector<SoAComplex> cut_psi(MAX_L_SYSTEM);
        double* cut_psi_0 = new double[m_max_grid_size * MAX_L_SYSTEM * 2];
        for (int l = 0; l < MAX_L_SYSTEM; ++l) {
            cut_psi[l].re = cut_psi_0 + m_max_grid_size * l * 2;
            cut_psi[l].im = cut_psi_0 + m_max_grid_size * (l * 2 + 1);
        }

        watch_pp.Record(2);


        int num_valid = 0;
        for (int ni = 0; ni < num_nuclei; ++ni) {
            auto mycomm = m_comm4atoms.GetComm(ni);
            if (mycomm == MPI_COMM_NULL) {
                continue;
            }

            int index = info_range[ni];

            const int max_l = m_nonlocal_projectors[ni].max_l;

            //const auto& cut_range_block_l = m_nonlocal_block[ni].range_blocks;
            //const auto& nonlocal_size = m_nonlocal_block[ni].grid_sizes;
            const auto& blocks = m_nonlocal_blocks[ni];

            for (int l = 0; l <= max_l; ++l) {

                CutSubgridByRanges(blocks[l].range_blocks, cut_psi[l].re, l_grid, l_psi.re);
                CutSubgridByRanges(blocks[l].range_blocks, cut_psi[l].im, l_grid, l_psi.im);
            }

            watch_pp.Record(3);

            for (int k = 0; k < m_nonlocal_projectors[ni].num_projectors; ++k) {
                const int l = m_nonlocal_projectors[ni].quantum_l[k];
                const auto* proj_re = m_nonlocal_projectors[ni].projector_expikx[id_spin_kpoint][k].re;
                const auto* proj_im = m_nonlocal_projectors[ni].projector_expikx[id_spin_kpoint][k].im;
                const auto* psi_re = cut_psi[l].re;
                const auto* psi_im = cut_psi[l].im;

                if (l == 0) {

                    const size_t nl_size = blocks[l].grid_sizes;
                    double inner_re = 0.0;
                    double inner_im = 0.0;
                    #pragma ivdep
                    for (int i = 0; i < nl_size; ++i) {
                        inner_re += proj_re[i] * psi_re[i] - proj_im[i] * psi_im[i];
                        inner_im += proj_im[i] * psi_re[i] + proj_re[i] * psi_im[i];
                    }
                    l_inner[index].r = inner_re * dV * Y_00;
                    l_inner[index].i = inner_im * dV * Y_00;
                    ++index;

                } else {
                    const size_t nl_size = blocks[l].grid_sizes;

                    for (int m1 = 0; m1 < 2 * l + 1; ++m1) {
                        const auto* Ylm = m_nonlocal_projectors[ni].Ylm[l * l + m1];
                    
                        double inner_re = 0.0;
                        double inner_im = 0.0;
                        #pragma ivdep
                        for (int i = 0; i < nl_size; ++i) {
                            inner_re += Ylm[i] * (proj_re[i] * psi_re[i] - proj_im[i] * psi_im[i]);
                            inner_im += Ylm[i] * (proj_im[i] * psi_re[i] + proj_re[i] * psi_im[i]);
                        }
                        l_inner[index].r = inner_re * dV;
                        l_inner[index].i = inner_im * dV;
                        ++index;
                    }
                }
            }

            watch_pp.Record(4);

            MPI_Iallreduce(l_inner + info_range[ni], sum_inner + info_range[ni],
                info_range[ni + 1] - info_range[ni], MPI_DOUBLE_COMPLEX, MPI_SUM, mycomm, &m_request4atoms[num_valid]);
            ++num_valid;


        }

        delete[]cut_psi_0;

        MPI_Waitall(num_valid, &m_request4atoms[0], &m_status4atoms[0]);
        watch_pp.Record(5);


    }

#ifdef FAST_GAMMA_POINT

    int mInnerNonlocalPsi_v2_gamma(const SoAComplex& l_psi, const GridRangeMPI& l_grid, const Nucleus* nuclei, int num_nuclei, OneComplex* l_inner, int id_spin_kpoint, std::vector<int>& info_range, double dV) {

        //const auto& nonlocal_block = m_nonlocal_block[id_spin_kpoint];

        const double Y_00 = 1.0 / sqrt(4.0 * M_PI);

        const int num_all_nonlocal = info_range[num_nuclei];
        OneComplex* sum_inner = l_inner + num_all_nonlocal;

        std::vector<SoAComplex> cut_psi(MAX_L_SYSTEM);
        double* cut_psi_0 = new double[m_max_grid_size * MAX_L_SYSTEM * 2];
        for (int l = 0; l < MAX_L_SYSTEM; ++l) {
            cut_psi[l].re = cut_psi_0 + m_max_grid_size * l * 2;
            cut_psi[l].im = cut_psi_0 + m_max_grid_size * (l * 2 + 1);
        }

        watch_pp.Record(2);


        int num_valid = 0;
        for (int ni = 0; ni < num_nuclei; ++ni) {
            auto mycomm = m_comm4atoms.GetComm(ni);
            if (mycomm == MPI_COMM_NULL) {
                continue;
            }

            int index = info_range[ni];

            const int max_l = m_nonlocal_projectors[ni].max_l;

            //const auto& cut_range_block_l = m_nonlocal_block[ni].range_blocks;
            //const auto& nonlocal_size = m_nonlocal_block[ni].grid_sizes;
            const auto& blocks = m_nonlocal_blocks[ni];

            for (int l = 0; l <= max_l; ++l) {
                CutSubgridByRanges(blocks[l].range_blocks, cut_psi[l].re, l_grid, l_psi.re);
            }
            for (int l = 0; l <= max_l; ++l) {
                CutSubgridByRanges(blocks[l].range_blocks, cut_psi[l].im, l_grid, l_psi.im);
            }

            watch_pp.Record(3);

            for (int k = 0; k < m_nonlocal_projectors[ni].num_projectors; ++k) {
                const int l = m_nonlocal_projectors[ni].quantum_l[k];
                const auto* proj_re = m_nonlocal_projectors[ni].projector_expikx[id_spin_kpoint][k].re;
                //const auto* proj_im = m_nonlocal_projectors[ni].projector_expikx[id_spin_kpoint][k].im;
                const auto* psi_re = cut_psi[l].re;
                const auto* psi_im = cut_psi[l].im;

                if (l == 0) {

                    const size_t nl_size = blocks[l].grid_sizes;
                    double inner_re = 0.0;
                    double inner_im = 0.0;
                #pragma ivdep
                    for (int i = 0; i < nl_size; ++i) {
                        inner_re += proj_re[i] * psi_re[i];
                        inner_im += proj_re[i] * psi_im[i];
                    }
                    l_inner[index].r = inner_re * dV * Y_00;
                    l_inner[index].i = inner_im * dV * Y_00;
                    ++index;

                } else {
                    const size_t nl_size = blocks[l].grid_sizes;

                    for (int m1 = 0; m1 < 2 * l + 1; ++m1) {
                        const auto* Ylm = m_nonlocal_projectors[ni].Ylm[l * l + m1];

                        double inner_re = 0.0;
                        double inner_im = 0.0;
                        #pragma ivdep
                        for (int i = 0; i < nl_size; ++i) {
                            const double YP = Ylm[i] * proj_re[i];
                            inner_re += YP * psi_re[i];
                            inner_im += YP * psi_im[i];
                        }
                        l_inner[index].r = inner_re * dV;
                        l_inner[index].i = inner_im * dV;
                        ++index;
                    }
                }
            }

            watch_pp.Record(4);

            MPI_Iallreduce(l_inner + info_range[ni], sum_inner + info_range[ni],
                info_range[ni + 1] - info_range[ni], MPI_DOUBLE_COMPLEX, MPI_SUM, mycomm, &m_request4atoms[num_valid]);
            ++num_valid;


        }

        delete[]cut_psi_0;

        return num_valid;

    }

#endif

    void mInnerNonlocalPsi_v2_1(const SoAComplex& l_psi, const GridRangeMPI& l_grid, const Nucleus* nuclei, int num_nuclei, OneComplex* l_inner, int id_spin_kpoint, std::vector<int>& info_range, double dV) {

        //const auto& nonlocal_block = m_nonlocal_block[id_spin_kpoint];

        const double Y_00 = 1.0 / sqrt(4.0 * M_PI);

        const int num_all_nonlocal = info_range[num_nuclei];
        OneComplex* sum_inner = l_inner + num_all_nonlocal;

        std::vector<SoAComplex> cut_psi(MAX_L_SYSTEM);
        double* cut_psi_0 = new double[m_max_grid_size * MAX_L_SYSTEM * 2];
        for (int l = 0; l < MAX_L_SYSTEM; ++l) {
            cut_psi[l].re = cut_psi_0 + m_max_grid_size * l * 2;
            cut_psi[l].im = cut_psi_0 + m_max_grid_size * (l * 2 + 1);
        }

        OneComplex* R_psi = new OneComplex[m_max_grid_size ];

        watch_pp.Record(2);


        int num_valid = 0;
        for (int ni = 0; ni < num_nuclei; ++ni) {
            auto mycomm = m_comm4atoms.GetComm(ni);
            if (mycomm == MPI_COMM_NULL) {
                continue;
            }

            int index = info_range[ni];

            const int max_l = m_nonlocal_projectors[ni].max_l;

            //const auto& cut_range_block_l = m_nonlocal_block[ni].range_blocks;
            //const auto& nonlocal_size = m_nonlocal_block[ni].grid_sizes;
            const auto& blocks = m_nonlocal_blocks[ni];

            for (int l = 0; l <= max_l; ++l) {

                CutSubgridByRanges(blocks[l].range_blocks, cut_psi[l].re, l_grid, l_psi.re);
                CutSubgridByRanges(blocks[l].range_blocks, cut_psi[l].im, l_grid, l_psi.im);
            }

            watch_pp.Record(3);

            for (int k = 0; k < m_nonlocal_projectors[ni].num_projectors; ++k) {
                const int l = m_nonlocal_projectors[ni].quantum_l[k];
                const auto* proj_re = m_nonlocal_projectors[ni].projector_expikx[id_spin_kpoint][k].re;
                const auto* proj_im = m_nonlocal_projectors[ni].projector_expikx[id_spin_kpoint][k].im;
                const auto* psi_re = cut_psi[l].re;
                const auto* psi_im = cut_psi[l].im;

                if (l == 0) {
                    /*
                    double inner = ForReduce2(cut_grid[l], 0.0,
                    [&proj, &l_psi](int64_t i) { return proj[i] * l_psi[i]; });
                    */
                    size_t nl_size = blocks[l].grid_sizes;
                    double inner_re = 0.0;
                    double inner_im = 0.0;
                    #pragma ivdep
                    for (int i = 0; i < nl_size; ++i) {
                        inner_re += proj_re[i] * psi_re[i] - proj_im[i] * psi_im[i];
                        inner_im += proj_im[i] * psi_re[i] + proj_re[i] * psi_im[i];
                    }
                    l_inner[index].r = inner_re * dV * Y_00;
                    l_inner[index].i = inner_im * dV * Y_00;
                    ++index;

                } else {
                    const size_t nl_size = blocks[l].grid_sizes;
                    #pragma ivdep
                    for (int i = 0; i < nl_size; ++i) {
                        R_psi[i].r = (proj_re[i] * psi_re[i] - proj_im[i] * psi_im[i]);
                        R_psi[i].i = (proj_im[i] * psi_re[i] + proj_re[i] * psi_im[i]);
                    }

                    #pragma ivdep
                    for (int m1 = 0; m1 < 2 * l + 1; ++m1) {
                        const auto* Ylm = m_nonlocal_projectors[ni].Ylm[l * l + m1];
                        //double inner = ForReduce2(cut_grid[l], 0.0,
                        //	[&proj, &Ylm, &l_psi](int64_t i) { return proj[i] * Ylm[i] * l_psi[i]; });
                        double inner_re = 0.0;
                        double inner_im = 0.0;
                        for (int i = 0; i < nl_size; ++i) {
                            inner_re += Ylm[i] * R_psi[i].r;
                            inner_im += Ylm[i] * R_psi[i].i;
                        }
                        l_inner[index].r = inner_re * dV;
                        l_inner[index].i = inner_im * dV;
                        ++index;
                    }
                }
            }

            watch_pp.Record(4);

            MPI_Iallreduce(l_inner + info_range[ni], sum_inner + info_range[ni],
                info_range[ni + 1] - info_range[ni], MPI_DOUBLE_COMPLEX, MPI_SUM, mycomm, &m_request4atoms[num_valid]);
            ++num_valid;


        }

        delete[]cut_psi_0;
        delete[] R_psi;

        MPI_Waitall(num_valid, &m_request4atoms[0], &m_status4atoms[0]);
        watch_pp.Record(5);


    }

#ifdef USE_COMPLEX_PROJECTOR
    //convert complex
    void mInnerNonlocalPsi_v2_2(const SoAComplex& l_psi, const GridRangeMPI& l_grid, const Nucleus* nuclei, int num_nuclei, OneComplex* l_inner, int id_spin_kpoint, std::vector<int>& info_range, double dV) {

        //const auto& nonlocal_block = m_nonlocal_block[id_spin_kpoint];

        const double Y_00 = 1.0 / sqrt(4.0 * M_PI);

        const int num_all_nonlocal = info_range[num_nuclei];
        OneComplex* sum_inner = l_inner + num_all_nonlocal;

        std::vector<OneComplex*> cut_psi(MAX_L_SYSTEM);
        double* cut_psi_0 = new double[m_max_grid_size * MAX_L_SYSTEM * 2];
        for (int l = 0; l < MAX_L_SYSTEM; ++l) {
            cut_psi[l] = (OneComplex*)(cut_psi_0 + m_max_grid_size * l * 2);            
        }

        watch_pp.Record(2);


        int num_valid = 0;
        for (int ni = 0; ni < num_nuclei; ++ni) {
            auto mycomm = m_comm4atoms.GetComm(ni);
            if (mycomm == MPI_COMM_NULL) {
                continue;
            }

            int index = info_range[ni];

            const int max_l = m_nonlocal_projectors[ni].max_l;

            //const auto& cut_range_block_l = m_nonlocal_block[ni].range_blocks;
            //const auto& nonlocal_size = m_nonlocal_block[ni].grid_sizes;
            const auto& blocks = m_nonlocal_blocks[ni];

            for (int l = 0; l <= max_l; ++l) {

                CutSubgridByRanges_SOAtoCOMPLEX(blocks[l].range_blocks, (double*)cut_psi[l], l_grid, l_psi.re, l_psi.im);
                //CutSubgridByRanges(blocks[l].range_blocks, cut_psi[l].im, l_grid, l_psi.im);
            }

            watch_pp.Record(3);

            for (int k = 0; k < m_nonlocal_projectors[ni].num_projectors; ++k) {
                const int l = m_nonlocal_projectors[ni].quantum_l[k];
                const auto* proj = m_nonlocal_projectors[ni].c_projector_expikx[id_spin_kpoint][k];
                //const auto* proj_im = m_nonlocal_projectors[ni].projector_expikx[id_spin_kpoint][k].im;
                const auto* psi = cut_psi[l];
                //const auto* psi_im = cut_psi[l];

                if (l == 0) {

                    const size_t nl_size = blocks[l].grid_sizes;
                    double inner_re = 0.0;
                    double inner_im = 0.0;
#pragma ivdep
                    for (int i = 0; i < nl_size; ++i) {
                        inner_re += real(proj[i]) * real(psi[i]) - imag(proj[i]) * imag(psi[i]);
                        inner_im += imag(proj[i]) * real(psi[i]) + real(proj[i]) * imag(psi[i]);
                    }
                    l_inner[index].r = inner_re * dV * Y_00;
                    l_inner[index].i = inner_im * dV * Y_00;
                    ++index;

                } else {
                    const size_t nl_size = blocks[l].grid_sizes;

                    for (int m1 = 0; m1 < 2 * l + 1; ++m1) {
                        const auto* Ylm = m_nonlocal_projectors[ni].Ylm[l * l + m1];

                        double inner_re = 0.0;
                        double inner_im = 0.0;
#pragma ivdep
                        for (int i = 0; i < nl_size; ++i) {
                            inner_re += Ylm[i] * (real(proj[i]) * real(psi[i]) - imag(proj[i]) * imag(psi[i]));
                            inner_im += Ylm[i] * (imag(proj[i]) * real(psi[i]) + real(proj[i]) * imag(psi[i]));
                        }
                        l_inner[index].r = inner_re * dV;
                        l_inner[index].i = inner_im * dV;
                        ++index;
                    }
                }
            }

            watch_pp.Record(4);

            MPI_Iallreduce(l_inner + info_range[ni], sum_inner + info_range[ni],
                info_range[ni + 1] - info_range[ni], MPI_DOUBLE_COMPLEX, MPI_SUM, mycomm, &m_request4atoms[num_valid]);
            ++num_valid;


        }

        delete[]cut_psi_0;

        MPI_Waitall(num_valid, &m_request4atoms[0], &m_status4atoms[0]);
        watch_pp.Record(5);


    }
#endif

#ifdef EXPAND_PROJECTOR_YLM_Z

#ifdef USE_COMPLEX_PROJECTOR
    void mInnerNonlocalPsi_v3(const SoAComplex& l_psi, const GridRangeMPI& l_grid, const Nucleus* nuclei, int num_nuclei, OneComplex* l_inner, int id_spin_kpoint, std::vector<int>& info_range, double dV) {



        const int num_all_nonlocal = info_range[num_nuclei];
        OneComplex* sum_inner = l_inner + num_all_nonlocal;

        auto* cut_psi_0 = new OneComplex[m_max_grid_size];

        auto* inner = new OneComplex[num_all_nonlocal];

        int num_valid = 0;
        for (int ni = 0; ni < num_nuclei; ++ni) {
            auto mycomm = m_comm4atoms.GetComm(ni);
            if (mycomm == MPI_COMM_NULL) {
                continue;
            }

            int index = info_range[ni];

            const int max_l = m_nonlocal_projectors[ni].max_l;
            std::vector<OneComplex*> cut_psi(max_l + 1);


            auto* proj_p = m_nonlocal_projectors[ni].projector_RY_expikx[id_spin_kpoint];
            int64_t offset = 0;
            const auto& blocks = m_nonlocal_blocks[ni];

            for (int l = 0; l <= max_l; ++l) {

                const auto grid_size = blocks[l].grid_sizes;

                //size_t nl_size = blocks[l].grid_sizes;
                const int width = (m_nonlocal_projectors[ni].head_projector_l[l + 1] - m_nonlocal_projectors[ni].head_projector_l[l]);

                if (grid_size == 0) {

                    for (int j = 0; j < width; ++j) {
                        l_inner[index + j].r = 0.0;
                        l_inner[index + j].i = 0.0;
                    }

                    index += width;
                    //offset += grid_size * width * 2; //because grid_size==0
                    continue;
                }

                cut_psi[l] = cut_psi_0;
                

                CutSubgridByRanges_SOAtoCOMPLEX(blocks[l].range_blocks, (double*)cut_psi[l], l_grid, l_psi.re, l_psi.im);
                

                watch_pp.Record(3);

                auto* proj = proj_p + offset;
                blas_ZGEMV_t(width, grid_size, proj, cut_psi[l], inner, { dV,0.0 }, { 0.0,0.0 });
                //blas_DGEMV_t(width * 2, grid_size, proj_re, cut_psi[l].im, inner + width * 2, dV, 0.0);

                for (int j = 0; j < width; ++j) {
                    l_inner[index + j] = inner[j];
                }

                index += width;
                offset += grid_size * width;

                //delete[] inner;
                //delete[] cut_psi[l].re;
                watch_pp.Record(4);
            }


            MPI_Iallreduce(l_inner + info_range[ni], sum_inner + info_range[ni],
                info_range[ni + 1] - info_range[ni], MPI_DOUBLE_COMPLEX, MPI_SUM, mycomm, &m_request4atoms[num_valid]);
            ++num_valid;


        }

        delete[]cut_psi_0;
        delete[] inner;

        MPI_Waitall(num_valid, &m_request4atoms[0], &m_status4atoms[0]);
        watch_pp.Record(5);


    }

#else
    void mInnerNonlocalPsi_v3(const SoAComplex& l_psi, const GridRangeMPI& l_grid, const Nucleus* nuclei, int num_nuclei, OneComplex* l_inner, int id_spin_kpoint, std::vector<int>& info_range, double dV) {



        const int num_all_nonlocal = info_range[num_nuclei];
        OneComplex* sum_inner = l_inner + num_all_nonlocal;

        double* cut_psi_0 = new double[m_max_grid_size * 2];

        double* inner = new double[num_all_nonlocal * 4];

        int num_valid = 0;
        for (int ni = 0; ni < num_nuclei; ++ni) {
            auto mycomm = m_comm4atoms.GetComm(ni);
            if (mycomm == MPI_COMM_NULL) {
                continue;
            }

            int index = info_range[ni];

            const int max_l = m_nonlocal_projectors[ni].max_l;
            std::vector<SoAComplex> cut_psi(max_l + 1);


            auto* proj_p = m_nonlocal_projectors[ni].projector_RY_expikx[id_spin_kpoint];
            int64_t offset = 0;
            const auto& blocks = m_nonlocal_blocks[ni];

            for (int l = 0; l <= max_l; ++l) {

                const auto grid_size = blocks[l].grid_sizes;

                //size_t nl_size = blocks[l].grid_sizes;
                const int width = (m_nonlocal_projectors[ni].head_projector_l[l + 1] - m_nonlocal_projectors[ni].head_projector_l[l]);

                if (grid_size == 0) {

                    for (int j = 0; j < width; ++j) {
                        l_inner[index + j].r = 0.0;
                        l_inner[index + j].i = 0.0;
                    }

                    index += width;
                    //offset += grid_size * width * 2; //because grid_size==0
                    continue;
                }

                cut_psi[l].re = cut_psi_0;
                cut_psi[l].im = cut_psi[l].re + grid_size;

                CutSubgridByRanges(blocks[l].range_blocks, cut_psi[l].re, l_grid, l_psi.re);
                CutSubgridByRanges(blocks[l].range_blocks, cut_psi[l].im, l_grid, l_psi.im);

                watch_pp.Record(3);

                auto* proj_re = proj_p + offset;
#if 1
                blas_DGEMV_t(width * 2, grid_size, proj_re, cut_psi[l].re, inner, dV, 0.0);
                blas_DGEMV_t(width * 2, grid_size, proj_re, cut_psi[l].im, inner + width * 2, dV, 0.0);

                for (int j = 0; j < width; ++j) {
                    l_inner[index + j].r = inner[2 * j + 0] - inner[2 * j + width * 2 + 1];
                    l_inner[index + j].i = inner[2 * j + 1] + inner[2 * j + width * 2 + 0];
                }
#else
                auto* psi_re = cut_psi[l].re;
                blas_DGEMM_t(2, width * 2, grid_size, psi_re, proj_re, inner, 2, dV, 0.0);

                for (int j = 0; j < width; ++j) {                    
                    l_inner[index + j].r = inner[4 * j + 0] - inner[4 * j + 2 + 1];
                    l_inner[index + j].i = inner[4 * j + 1] + inner[4 * j + 2 + 0];
                }
#endif
                index += width;
                offset += grid_size * width * 2;

                //delete[] inner;
                //delete[] cut_psi[l].re;
                watch_pp.Record(4);
            }


            MPI_Iallreduce(l_inner + info_range[ni], sum_inner + info_range[ni],
                info_range[ni + 1] - info_range[ni], MPI_DOUBLE_COMPLEX, MPI_SUM, mycomm, &m_request4atoms[num_valid]);
            ++num_valid;


        }

        delete[]cut_psi_0;
        delete[] inner;

        MPI_Waitall(num_valid, &m_request4atoms[0], &m_status4atoms[0]);
        watch_pp.Record(5);


    }

#ifdef FAST_GAMMA_POINT
#if 1
    
    void mInnerNonlocalPsi_v3_gamma(const SoAComplex& l_psi, const GridRangeMPI& l_grid, const Nucleus* nuclei, int num_nuclei, OneComplex* l_inner, int id_spin_kpoint, std::vector<int>& info_range, double dV) {



        const int num_all_nonlocal = info_range[num_nuclei];
        OneComplex* sum_inner = l_inner + num_all_nonlocal;


        std::vector<SoAComplex> cut_psi(MAX_L_SYSTEM);
        double* cut_psi_0 = new double[m_max_grid_size * MAX_L_SYSTEM * 2];
        for (int l = 0; l < MAX_L_SYSTEM; ++l) {
            cut_psi[l].re = cut_psi_0 + m_max_grid_size * l * 2;
            cut_psi[l].im = cut_psi_0 + m_max_grid_size * (l * 2 + 1);
        }

        
        int num_valid = 0;
        for (int ni = 0; ni < num_nuclei; ++ni) {
            auto mycomm = m_comm4atoms.GetComm(ni);
            if (mycomm == MPI_COMM_NULL) {
                continue;
            }

            int index = info_range[ni];

            const int max_l = m_nonlocal_projectors[ni].max_l;


            auto* proj_p = m_nonlocal_projectors[ni].projector_RY_expikx[id_spin_kpoint];
            int64_t offset = 0;
            const auto& blocks = m_nonlocal_blocks[ni];

            for (int l = 0; l <= max_l; ++l) {
                CutSubgridByRanges(blocks[l].range_blocks, cut_psi[l].re, l_grid, l_psi.re);
            }
            for (int l = 0; l <= max_l; ++l) {
                CutSubgridByRanges(blocks[l].range_blocks, cut_psi[l].im, l_grid, l_psi.im);
            }
            watch_pp.Record(3);


            for (int l = 0; l <= max_l; ++l) {
                const auto grid_size = blocks[l].grid_sizes;

                //size_t nl_size = blocks[l].grid_sizes;
                const int width = (m_nonlocal_projectors[ni].head_projector_l[l + 1] - m_nonlocal_projectors[ni].head_projector_l[l]);

                if (grid_size == 0) {

                    for (int j = 0; j < width; ++j) {
                        l_inner[index + j].r = 0.0;
                        l_inner[index + j].i = 0.0;
                    }

                    index += width;
                    //offset += grid_size * width * 2; //because grid_size==0
                    continue;
                }

#if 1
                const auto* psi_re = cut_psi[l].re;
                const auto* psi_im = cut_psi[l].im;
                for (int w = 0; w < width; ++w) {
                    auto* proj_re = proj_p + offset;
                    double inner_re = 0.0;
                    double inner_im = 0.0;
#pragma ivdep
                    for (int i = 0; i < grid_size; ++i) {
                        inner_re += proj_re[i] * psi_re[i];
                        inner_im += proj_re[i] * psi_im[i];
                    }
                    l_inner[index].r = inner_re * dV;
                    l_inner[index].i = inner_im * dV;
                    ++index;
                    offset += grid_size;
                }
#elif 1
                auto* proj_re = proj_p + offset;
                blas_DGEMV_t(width, grid_size, proj_re, cut_psi[l].re, inner, dV, 0.0);
                blas_DGEMV_t(width, grid_size, proj_re, cut_psi[l].im, inner + width, dV, 0.0);

#pragma ivdep
                for (int j = 0; j < width; ++j) {
                    l_inner[index + j].r = inner[j + 0];
                    l_inner[index + j].i = inner[j + width + 0];
                }
                index += width;
                offset += grid_size * width;
#else
                auto* proj_re = proj_p + offset;
                auto* psi_re = cut_psi[l].re;
                blas_DGEMM_t(2, width, grid_size, psi_re, proj_re, (double*)&l_inner[index], 2, dV, 0.0);


                index += width;
                offset += grid_size * width;
#endif

                //delete[] inner;
                //delete[] cut_psi[l].re;
                
            }
            watch_pp.Record(4);


            MPI_Iallreduce(l_inner + info_range[ni], sum_inner + info_range[ni],
                info_range[ni + 1] - info_range[ni], MPI_DOUBLE_COMPLEX, MPI_SUM, mycomm, &m_request4atoms[num_valid]);
            ++num_valid;


        }

        delete[]cut_psi_0;
        //delete[] inner;

        MPI_Waitall(num_valid, &m_request4atoms[0], &m_status4atoms[0]);
        watch_pp.Record(5);


    }
#else
    void mInnerNonlocalPsi_v3_gamma(const SoAComplex& l_psi, const GridRangeMPI& l_grid, const Nucleus* nuclei, int num_nuclei, OneComplex* l_inner, int id_spin_kpoint, std::vector<int>& info_range, double dV) {



        const int num_all_nonlocal = info_range[num_nuclei];
        OneComplex* sum_inner = l_inner + num_all_nonlocal;

        double* cut_psi_0 = new double[m_max_grid_size * 2];

        //double* inner = new double[num_all_nonlocal * 2];

        int num_valid = 0;
        for (int ni = 0; ni < num_nuclei; ++ni) {
            auto mycomm = m_comm4atoms.GetComm(ni);
            if (mycomm == MPI_COMM_NULL) {
                continue;
            }

            int index = info_range[ni];

            const int max_l = m_nonlocal_projectors[ni].max_l;
            std::vector<SoAComplex> cut_psi(max_l + 1);


            auto* proj_p = m_nonlocal_projectors[ni].projector_RY_expikx[id_spin_kpoint];
            int64_t offset = 0;
            const auto& blocks = m_nonlocal_blocks[ni];

            for (int l = 0; l <= max_l; ++l) {

                const auto grid_size = blocks[l].grid_sizes;

                //size_t nl_size = blocks[l].grid_sizes;
                const int width = (m_nonlocal_projectors[ni].head_projector_l[l + 1] - m_nonlocal_projectors[ni].head_projector_l[l]);

                if (grid_size == 0) {

                    for (int j = 0; j < width; ++j) {
                        l_inner[index + j].r = 0.0;
                        l_inner[index + j].i = 0.0;
                    }

                    index += width;
                    //offset += grid_size * width * 2; //because grid_size==0
                    continue;
                }

                cut_psi[l].re = cut_psi_0;
                cut_psi[l].im = cut_psi[l].re + grid_size;

                CutSubgridByRanges(blocks[l].range_blocks, cut_psi[l].re, l_grid, l_psi.re);
                CutSubgridByRanges(blocks[l].range_blocks, cut_psi[l].im, l_grid, l_psi.im);

                watch_pp.Record(3);

                
#if 1
                const auto* psi_re = cut_psi[l].re;
                const auto* psi_im = cut_psi[l].im;
                for (int w = 0; w < width; ++w) {
                    auto* proj_re = proj_p + offset;
                    double inner_re = 0.0;
                    double inner_im = 0.0;
                    #pragma ivdep
                    for (int i = 0; i < grid_size; ++i) {
                        inner_re += proj_re[i] * psi_re[i];
                        inner_im += proj_re[i] * psi_im[i];
                    }
                    l_inner[index].r = inner_re * dV;
                    l_inner[index].i = inner_im * dV;
                    ++index;
                    offset += grid_size;
                }
#elif 1
                auto* proj_re = proj_p + offset;
                blas_DGEMV_t(width, grid_size, proj_re, cut_psi[l].re, inner, dV, 0.0);
                blas_DGEMV_t(width, grid_size, proj_re, cut_psi[l].im, inner + width, dV, 0.0);

                #pragma ivdep
                for (int j = 0; j < width; ++j) {
                    l_inner[index + j].r = inner[j + 0];
                    l_inner[index + j].i = inner[j + width + 0];
                }
                index += width;
                offset += grid_size * width;
#else
                auto* proj_re = proj_p + offset;
                auto* psi_re = cut_psi[l].re;
                blas_DGEMM_t(2, width, grid_size, psi_re, proj_re, (double*)&l_inner[index], 2, dV, 0.0);
                

                index += width;
                offset += grid_size * width;
#endif

                //delete[] inner;
                //delete[] cut_psi[l].re;
                watch_pp.Record(4);
            }


            MPI_Iallreduce(l_inner + info_range[ni], sum_inner + info_range[ni],
                info_range[ni + 1] - info_range[ni], MPI_DOUBLE_COMPLEX, MPI_SUM, mycomm, &m_request4atoms[num_valid]);
            ++num_valid;


        }

        delete[]cut_psi_0;
        //delete[] inner;

        MPI_Waitall(num_valid, &m_request4atoms[0], &m_status4atoms[0]);
        watch_pp.Record(5);


    }
#endif
#endif
#endif

#if 0
    void mInnerNonlocalPsi_b1(const SoAComplex& l_psi, const GridRangeMPI& l_grid, const Nucleus* nuclei, int num_nuclei, OneComplex* l_inner, int id_spin_kpoint, int num_all_nonlocal, double dV) {


        const double Y_00 = 1.0 / sqrt(4.0 * M_PI);


        int index = 0;
        for (int ni = 0; ni < num_nuclei; ++ni) {
            const int max_l = m_nonlocal_projectors[ni].max_l;
            std::vector<SoAComplex> cut_psi(max_l + 1);

            const auto& blocks = m_nonlocal_blocks[ni];

            for (int l = 0; l <= max_l; ++l) {

                cut_psi[l].re = new double[blocks[l].grid_sizes * 2];
                cut_psi[l].im = cut_psi[l].re + blocks[l].grid_sizes;

                CutSubgridByRanges(blocks[l].range_blocks, cut_psi[l].re, l_grid, l_psi.re);
                CutSubgridByRanges(blocks[l].range_blocks, cut_psi[l].im, l_grid, l_psi.im);
            }


            auto* proj_p = m_nonlocal_projectors[ni].projector_RY_expikx[id_spin_kpoint];
            int64_t offset = 0;
            for (int k = 0; k < m_nonlocal_projectors[ni].num_projectors; ++k) {
                const int l = m_nonlocal_projectors[ni].quantum_l[k];
                size_t nl_size = blocks[l].grid_sizes;
                const auto* psi_re = cut_psi[l].re;
                const auto* psi_im = cut_psi[l].im;

                if (l == 0) {
                    const auto* proj_re = proj_p + offset;
                    const auto* proj_im = proj_re + nl_size;
                    /*
                    double inner = ForReduce2(cut_grid[l], 0.0,
                    [&proj, &l_psi](int64_t i) { return proj[i] * l_psi[i]; });
                    */
                    size_t nl_size = blocks[l].grid_sizes;
                    double inner_re = 0.0;
                    double inner_im = 0.0;
                    for (int i = 0; i < nl_size; ++i) {
                        inner_re += proj_re[i] * psi_re[i] - proj_im[i] * psi_im[i];
                        inner_im += proj_im[i] * psi_re[i] + proj_re[i] * psi_im[i];
                    }
                    l_inner[index].r = inner_re * dV;
                    l_inner[index].i = inner_im * dV;
                    ++index;
                    offset += nl_size * 2;

                } else {
                    for (int m1 = 0; m1 < 2 * l + 1; ++m1) {
                        const auto* proj_re = proj_p + offset;
                        const auto* proj_im = proj_re + nl_size;
                        //const auto* Ylm = m_nonlocal_projectors[ni].Ylm[l * l + m1];
                        //double inner = ForReduce2(cut_grid[l], 0.0,
                        //	[&proj, &Ylm, &l_psi](int64_t i) { return proj[i] * Ylm[i] * l_psi[i]; });
                        //size_t nl_size = blocks[l].grid_sizes;
                        double inner_re = 0.0;
                        double inner_im = 0.0;
                        for (int i = 0; i < nl_size; ++i) {
                            inner_re += (proj_re[i] * psi_re[i] - proj_im[i] * psi_im[i]);
                            inner_im += (proj_im[i] * psi_re[i] + proj_re[i] * psi_im[i]);
                        }
                        l_inner[index].r = inner_re * dV;
                        l_inner[index].i = inner_im * dV;
                        ++index;
                        offset += nl_size * 2;

                    }
                }
            }

            for (int l = 0; l <= max_l; ++l) {
                delete[] cut_psi[l].re;
            }


        }


        if (num_all_nonlocal != index) {
            printf("ERROR: number of nonlocal is invalid"); fflush(stdout);
        }
        //printf("[%d] ProjectionPP 2: 0x%zx\n", GetProcessID(l_grid.mpi_comm), l_inner); fflush(stdout);
        //MPI_Barrier(l_grid.mpi_comm);
    }

    void mInnerNonlocalPsi_b2(const SoAComplex& l_psi, const SoAComplex& l_dpsi_dx, const SoAComplex& l_dpsi_dy, const SoAComplex& l_dpsi_dz, const GridRangeMPI& l_grid, const Nucleus* nuclei, int num_nuclei, OneComplex* l_inner, int id_spin_kpoint, int num_all_nonlocal, double dV) {


        const double Y_00 = 1.0 / sqrt(4.0 * M_PI);


        int index = 0;
        for (int ni = 0; ni < num_nuclei; ++ni) {
            const int max_l = m_nonlocal_projectors[ni].max_l;
            std::vector<SoAComplex> cut_psi(max_l + 1);

            const auto& blocks = m_nonlocal_blocks[ni];

            for (int l = 0; l <= max_l; ++l) {

                const auto grid_size = blocks[l].grid_sizes;
                cut_psi[l].re = new double[grid_size * 2 * 4];
                cut_psi[l].im = cut_psi[l].re + grid_size;

                CutSubgridByRanges(blocks[l].range_blocks, cut_psi[l].re, l_grid, l_psi.re);
                CutSubgridByRanges(blocks[l].range_blocks, cut_psi[l].im, l_grid, l_psi.im);
                CutSubgridByRanges(blocks[l].range_blocks, cut_psi[l].re+ grid_size * 2, l_grid, l_dpsi_dx.re);
                CutSubgridByRanges(blocks[l].range_blocks, cut_psi[l].im+ grid_size * 2, l_grid, l_dpsi_dx.im);
                CutSubgridByRanges(blocks[l].range_blocks, cut_psi[l].re+ grid_size * 4, l_grid, l_dpsi_dy.re);
                CutSubgridByRanges(blocks[l].range_blocks, cut_psi[l].im+ grid_size * 4, l_grid, l_dpsi_dy.im);
                CutSubgridByRanges(blocks[l].range_blocks, cut_psi[l].re+ grid_size * 6, l_grid, l_dpsi_dz.re);
                CutSubgridByRanges(blocks[l].range_blocks, cut_psi[l].im+ grid_size * 6, l_grid, l_dpsi_dz.im);
            }


            auto* proj_p = m_nonlocal_projectors[ni].projector_RY_expikx[id_spin_kpoint];
            int64_t offset = 0;
            for (int k = 0; k < m_nonlocal_projectors[ni].num_projectors; ++k) {
                const int l = m_nonlocal_projectors[ni].quantum_l[k];
                size_t nl_size = blocks[l].grid_sizes;
                const auto* psi_re = cut_psi[l].re;
                const auto* psi_im = cut_psi[l].im;

                const auto* dpsi_dx_re = cut_psi[l].re + nl_size * 2;
                const auto* dpsi_dx_im = cut_psi[l].im + nl_size * 2;

                const auto* dpsi_dy_re = cut_psi[l].re + nl_size * 4;
                const auto* dpsi_dy_im = cut_psi[l].im + nl_size * 4;

                const auto* dpsi_dz_re = cut_psi[l].re + nl_size * 6;
                const auto* dpsi_dz_im = cut_psi[l].im + nl_size * 6;

                if (l == 0) {
                    const auto* proj_re = proj_p + offset;
                    const auto* proj_im = proj_re + nl_size;
                    
                    double inner[8]{ 0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0 };
                    
                    for (int i = 0; i < nl_size; ++i) {
                        inner[0] += proj_re[i] * psi_re[i] - proj_im[i] * psi_im[i];
                        inner[1] += proj_im[i] * psi_re[i] + proj_re[i] * psi_im[i];
                        inner[2] += proj_re[i] * dpsi_dx_re[i] - proj_im[i] * dpsi_dx_im[i];
                        inner[3] += proj_im[i] * dpsi_dx_re[i] + proj_re[i] * dpsi_dx_im[i];
                        inner[4] += proj_re[i] * dpsi_dy_re[i] - proj_im[i] * dpsi_dy_im[i];
                        inner[5] += proj_im[i] * dpsi_dy_re[i] + proj_re[i] * dpsi_dy_im[i];
                        inner[6] += proj_re[i] * dpsi_dz_re[i] - proj_im[i] * dpsi_dz_im[i];
                        inner[7] += proj_im[i] * dpsi_dz_re[i] + proj_re[i] * dpsi_dz_im[i];
                    }
                    l_inner[index].r = inner[0] * dV;
                    l_inner[index].i = inner[1] * dV;
                    l_inner[index+1].r = inner[2] * dV;
                    l_inner[index+1].i = inner[3] * dV;
                    l_inner[index + 2].r = inner[4] * dV;
                    l_inner[index + 2].i = inner[5] * dV;
                    l_inner[index + 3].r = inner[6] * dV;
                    l_inner[index + 3].i = inner[7] * dV;
                    index += 4;;
                    offset += nl_size * 2;

                } else {
                    for (int m1 = 0; m1 < 2 * l + 1; ++m1) {
                        const auto* proj_re = proj_p + offset;
                        const auto* proj_im = proj_re + nl_size;
                        
                        double inner[8]{ 0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0 };
                        for (int i = 0; i < nl_size; ++i) {
                            inner[0] += proj_re[i] * psi_re[i] - proj_im[i] * psi_im[i];
                            inner[1] += proj_im[i] * psi_re[i] + proj_re[i] * psi_im[i];
                            inner[2] += proj_re[i] * dpsi_dx_re[i] - proj_im[i] * dpsi_dx_im[i];
                            inner[3] += proj_im[i] * dpsi_dx_re[i] + proj_re[i] * dpsi_dx_im[i];
                            inner[4] += proj_re[i] * dpsi_dy_re[i] - proj_im[i] * dpsi_dy_im[i];
                            inner[5] += proj_im[i] * dpsi_dy_re[i] + proj_re[i] * dpsi_dy_im[i];
                            inner[6] += proj_re[i] * dpsi_dz_re[i] - proj_im[i] * dpsi_dz_im[i];
                            inner[7] += proj_im[i] * dpsi_dz_re[i] + proj_re[i] * dpsi_dz_im[i];
                        }
                        l_inner[index].r = inner[0] * dV;
                        l_inner[index].i = inner[1] * dV;
                        l_inner[index + 1].r = inner[2] * dV;
                        l_inner[index + 1].i = inner[3] * dV;
                        l_inner[index + 2].r = inner[4] * dV;
                        l_inner[index + 2].i = inner[5] * dV;
                        l_inner[index + 3].r = inner[6] * dV;
                        l_inner[index + 3].i = inner[7] * dV;
                        index += 4;;
                        offset += nl_size * 2;

                    }
                }
            }

            for (int l = 0; l <= max_l; ++l) {
                delete[] cut_psi[l].re;
            }


        }


        if (num_all_nonlocal*4 != index) {
            printf("ERROR: number of nonlocal is invalid"); fflush(stdout);
        }
        //printf("[%d] ProjectionPP 2: 0x%zx\n", GetProcessID(l_grid.mpi_comm), l_inner); fflush(stdout);
        //MPI_Barrier(l_grid.mpi_comm);
    }
#endif


#ifdef USE_COMPLEX_PROJECTOR

    void mInnerNonlocalPsi_b3(const SoAComplex& l_psi, const SoAComplex& l_dpsi_dx, const SoAComplex& l_dpsi_dy, const SoAComplex& l_dpsi_dz, const GridRangeMPI& l_grid, const Nucleus* nuclei, int num_nuclei, OneComplex* l_inner, int id_spin_kpoint, int num_all_nonlocal, double dV) {


        OneComplex* cut_psi_0 = new OneComplex[m_max_grid_size * 4];

        int index = 0;
        for (int ni = 0; ni < num_nuclei; ++ni) {
            const int max_l = m_nonlocal_projectors[ni].max_l;
            std::vector<OneComplex*> cut_psi(max_l + 1);


            auto* proj_p = m_nonlocal_projectors[ni].projector_RY_expikx[id_spin_kpoint];
            int64_t offset = 0;
            const auto& blocks = m_nonlocal_blocks[ni];

            for (int l = 0; l <= max_l; ++l) {

                const auto grid_size = blocks[l].grid_sizes;

                //size_t nl_size = blocks[l].grid_sizes;
                const int width = (m_nonlocal_projectors[ni].head_projector_l[l + 1] - m_nonlocal_projectors[ni].head_projector_l[l]);

                if (grid_size == 0) {

                    for (int j = 0; j < width * 4; ++j) {
                        l_inner[index + j].r = 0.0;
                        l_inner[index + j].i = 0.0;
                    }

                    index += 4 * width;
                    //offset += grid_size * width * 2; //because grid_size==0
                    continue;
                }

                cut_psi[l] = cut_psi_0;
                //cut_psi[l].im = cut_psi[l].re + grid_size;

                CutSubgridByRanges_SOAtoCOMPLEX(blocks[l].range_blocks, (double*)cut_psi[l], l_grid, l_psi.re, l_psi.im);
                CutSubgridByRanges_SOAtoCOMPLEX(blocks[l].range_blocks, (double*)(cut_psi[l] + grid_size), l_grid, l_dpsi_dx.re, l_dpsi_dx.im);
                CutSubgridByRanges_SOAtoCOMPLEX(blocks[l].range_blocks, (double*)(cut_psi[l] + grid_size * 2), l_grid, l_dpsi_dy.re, l_dpsi_dy.im);
                CutSubgridByRanges_SOAtoCOMPLEX(blocks[l].range_blocks, (double*)(cut_psi[l] + grid_size * 3), l_grid, l_dpsi_dz.re, l_dpsi_dz.im);
                

                auto* psi_re = cut_psi[l];


                auto* proj_re = proj_p + offset;


                
#if 1
                blas_ZGEMM_t(4, width, grid_size, psi_re, proj_re, &l_inner[index], 4, { dV,0.0 }, { 0.0, 0.0 });

#else
                OneComplex* inner = new OneComplex[4 * width];
                blas_DGEMM_t(width * 2, 8, grid_size, proj_re, psi_re, inner, width * 2, dV, 0.0);
#pragma ivdep
                for (int j = 0; j < width; ++j) {
#pragma ivdep
                    for (int m = 0; m < 4; ++m) {
                        l_inner[index + j * 4 + m].r = inner[2 * j + m * width * 4 + 0] - inner[2 * j + m * width * 4 + width * 2 + 1];
                        l_inner[index + j * 4 + m].i = inner[2 * j + m * width * 4 + 1] + inner[2 * j + m * width * 4 + width * 2 + 0];

                    }
                }
                delete[] inner;
#endif
                index += 4 * width;
                offset += grid_size * width;

            }


        }

        delete[] cut_psi_0;

        if (num_all_nonlocal * 4 != index) {
            printf("ERROR: number of nonlocal is invalid"); fflush(stdout);
        }

    }

#else
    void mInnerNonlocalPsi_b3(const SoAComplex& l_psi, const SoAComplex& l_dpsi_dx, const SoAComplex& l_dpsi_dy, const SoAComplex& l_dpsi_dz, const GridRangeMPI& l_grid, const Nucleus* nuclei, int num_nuclei, OneComplex* l_inner, int id_spin_kpoint, int num_all_nonlocal, double dV) {


        double* cut_psi_0 = new double[m_max_grid_size * 2 * 4];

        int index = 0;
        for (int ni = 0; ni < num_nuclei; ++ni) {
            const int max_l = m_nonlocal_projectors[ni].max_l;
            std::vector<SoAComplex> cut_psi(max_l + 1);


            auto* proj_p = m_nonlocal_projectors[ni].projector_RY_expikx[id_spin_kpoint];
            int64_t offset = 0;
            const auto& blocks = m_nonlocal_blocks[ni];

            for (int l = 0; l <= max_l; ++l) {

                const auto grid_size = blocks[l].grid_sizes;

                //size_t nl_size = blocks[l].grid_sizes;
                const int width = (m_nonlocal_projectors[ni].head_projector_l[l + 1] - m_nonlocal_projectors[ni].head_projector_l[l]);

                if (grid_size == 0) {

                    for (int j = 0; j < width * 4; ++j) {
                        l_inner[index + j].r = 0.0;
                        l_inner[index + j].i = 0.0;
                    }

                    index += 4 * width;
                    //offset += grid_size * width * 2; //because grid_size==0
                    continue;
                }

                cut_psi[l].re = cut_psi_0;
                cut_psi[l].im = cut_psi[l].re + grid_size;

                CutSubgridByRanges(blocks[l].range_blocks, cut_psi[l].re, l_grid, l_psi.re);
                CutSubgridByRanges(blocks[l].range_blocks, cut_psi[l].im, l_grid, l_psi.im);
                CutSubgridByRanges(blocks[l].range_blocks, cut_psi[l].re + grid_size * 2, l_grid, l_dpsi_dx.re);
                CutSubgridByRanges(blocks[l].range_blocks, cut_psi[l].im + grid_size * 2, l_grid, l_dpsi_dx.im);
                CutSubgridByRanges(blocks[l].range_blocks, cut_psi[l].re + grid_size * 4, l_grid, l_dpsi_dy.re);
                CutSubgridByRanges(blocks[l].range_blocks, cut_psi[l].im + grid_size * 4, l_grid, l_dpsi_dy.im);
                CutSubgridByRanges(blocks[l].range_blocks, cut_psi[l].re + grid_size * 6, l_grid, l_dpsi_dz.re);
                CutSubgridByRanges(blocks[l].range_blocks, cut_psi[l].im + grid_size * 6, l_grid, l_dpsi_dz.im);
                

                auto* psi_re = cut_psi[l].re;


                auto* proj_re = proj_p + offset;


                double* inner = new double[8 * width * 2];
#if 1
                blas_DGEMM_t(8, width * 2, grid_size, psi_re, proj_re, inner, 8, dV, 0.0);
                #pragma ivdep
                for (int j = 0; j < width; ++j) {
                    #pragma ivdep
                    for (int m = 0; m < 4; ++m) {
                        l_inner[index + j * 4 + m].r = inner[16 * j + m * 2 + 0] - inner[16 * j + m * 2 + 8 + 1];
                        l_inner[index + j * 4 + m].i = inner[16 * j + m * 2 + 1] + inner[16 * j + m * 2 + 8 + 0];

                    }
                }
#else
                blas_DGEMM_t(width * 2, 8, grid_size, proj_re, psi_re, inner, width * 2, dV, 0.0);
                #pragma ivdep
                for (int j = 0; j < width; ++j) {
                    #pragma ivdep
                    for (int m = 0; m < 4; ++m) {
                        l_inner[index + j * 4 + m].r = inner[2 * j + m * width * 4 + 0] - inner[2 * j + m * width * 4 + width * 2 + 1];
                        l_inner[index + j * 4 + m].i = inner[2 * j + m * width * 4 + 1] + inner[2 * j + m * width * 4 + width * 2 + 0];

                    }
                }
#endif
                index += 4 * width;
                offset += grid_size * width * 2;

                delete[] inner;
            }


        }

        delete[] cut_psi_0;

        if (num_all_nonlocal * 4 != index) {
            printf("ERROR: number of nonlocal is invalid"); fflush(stdout);
        }

    }

    void mInnerNonlocalPsi_b3_gamma(const SoAComplex& l_psi, const SoAComplex& l_dpsi_dx, const SoAComplex& l_dpsi_dy, const SoAComplex& l_dpsi_dz, const GridRangeMPI& l_grid, const Nucleus* nuclei, int num_nuclei, OneComplex* l_inner, int id_spin_kpoint, int num_all_nonlocal, double dV) {


        double* cut_psi_0 = new double[m_max_grid_size * 2 * 4];

        int index = 0;
        for (int ni = 0; ni < num_nuclei; ++ni) {
            const int max_l = m_nonlocal_projectors[ni].max_l;
            std::vector<SoAComplex> cut_psi(max_l + 1);


            auto* proj_p = m_nonlocal_projectors[ni].projector_RY_expikx[id_spin_kpoint];
            int64_t offset = 0;
            const auto& blocks = m_nonlocal_blocks[ni];

            for (int l = 0; l <= max_l; ++l) {

                const auto grid_size = blocks[l].grid_sizes;

                //size_t nl_size = blocks[l].grid_sizes;
                const int width = (m_nonlocal_projectors[ni].head_projector_l[l + 1] - m_nonlocal_projectors[ni].head_projector_l[l]);

                if (grid_size == 0) {

                    for (int j = 0; j < width * 4; ++j) {
                        l_inner[index + j].r = 0.0;
                        l_inner[index + j].i = 0.0;
                    }

                    index += 4 * width;
                    //offset += grid_size * width * 2; //because grid_size==0
                    continue;
                }

                cut_psi[l].re = cut_psi_0;
                cut_psi[l].im = cut_psi[l].re + grid_size;

                CutSubgridByRanges(blocks[l].range_blocks, cut_psi[l].re, l_grid, l_psi.re);
                CutSubgridByRanges(blocks[l].range_blocks, cut_psi[l].im, l_grid, l_psi.im);
                CutSubgridByRanges(blocks[l].range_blocks, cut_psi[l].re + grid_size * 2, l_grid, l_dpsi_dx.re);
                CutSubgridByRanges(blocks[l].range_blocks, cut_psi[l].im + grid_size * 2, l_grid, l_dpsi_dx.im);
                CutSubgridByRanges(blocks[l].range_blocks, cut_psi[l].re + grid_size * 4, l_grid, l_dpsi_dy.re);
                CutSubgridByRanges(blocks[l].range_blocks, cut_psi[l].im + grid_size * 4, l_grid, l_dpsi_dy.im);
                CutSubgridByRanges(blocks[l].range_blocks, cut_psi[l].re + grid_size * 6, l_grid, l_dpsi_dz.re);
                CutSubgridByRanges(blocks[l].range_blocks, cut_psi[l].im + grid_size * 6, l_grid, l_dpsi_dz.im);


                auto* psi_re = cut_psi[l].re;


                auto* proj_re = proj_p + offset;


                double* inner = new double[8 * width];
#if 1
                blas_DGEMM_t(8, width , grid_size, psi_re, proj_re, inner, 8, dV, 0.0);
#pragma ivdep
                for (int j = 0; j < width; ++j) {
#pragma ivdep
                    for (int m = 0; m < 4; ++m) {
                        l_inner[index + j * 4 + m].r = inner[8 * j + m * 2 + 0];
                        l_inner[index + j * 4 + m].i = inner[8 * j + m * 2 + 1];

                    }
                }
#else
                blas_DGEMM_t(width, 8, grid_size, proj_re, psi_re, inner, width, dV, 0.0);
#pragma ivdep
                for (int j = 0; j < width; ++j) {
#pragma ivdep
                    for (int m = 0; m < 4; ++m) {
                        l_inner[index + j * 4 + m].r = inner[2 * j + m * width * 2 + 0];
                        l_inner[index + j * 4 + m].i = inner[2 * j + m * width * 2 + width];

                    }
                }
#endif
                index += 4 * width;
                offset += grid_size * width;

                delete[] inner;
            }


        }

        delete[] cut_psi_0;

        if (num_all_nonlocal * 4 != index) {
            printf("ERROR: number of nonlocal is invalid"); fflush(stdout);
        }

    }
#endif

#ifndef USE_COMPLEX_PROJECTOR
    void mInnerNonlocalPsi_b4(int num_psi, const SoAComplex* l_psi, const double* l_dpsi_dxyz, double* work, const GridRangeMPI& l_grid, const Nucleus* nuclei, int num_nuclei, OneComplex* l_inner, int id_spin_kpoint, int num_all_nonlocal, double dV) {

        const size_t local_size = l_grid.Size3D();

        int index = 0;
        for (int ni = 0; ni < num_nuclei; ++ni) {
            const int max_l = m_nonlocal_projectors[ni].max_l;
            //std::vector<SoAComplex> cut_psi(max_l + 1);


            auto* proj_p = m_nonlocal_projectors[ni].projector_RY_expikx[id_spin_kpoint];
            int64_t offset = 0;
            const auto& blocks = m_nonlocal_blocks[ni];

            for (int l = 0; l <= max_l; ++l) {

                const auto grid_size = blocks[l].grid_sizes;

                //size_t nl_size = blocks[l].grid_sizes;
                const int width = (m_nonlocal_projectors[ni].head_projector_l[l + 1] - m_nonlocal_projectors[ni].head_projector_l[l]);

                if (grid_size == 0) {

                    for (int j = 0; j < width * 4* num_psi; ++j) {
                        l_inner[index + j].r = 0.0;
                        l_inner[index + j].i = 0.0;
                    }

                    index += 4 * width* num_psi;
                    //offset += grid_size * width * 2; //because grid_size==0
                    continue;
                }
                double* cut_psi_all = work;
                for (int n = 0; n < num_psi; ++n) {
                    auto* cut_psi = work + n * grid_size * 8;
                    

                    CutSubgridByRanges(blocks[l].range_blocks, cut_psi, l_grid, l_psi[n].re);
                    CutSubgridByRanges(blocks[l].range_blocks, cut_psi + grid_size * 1, l_grid, l_psi[n].im);
                    CutSubgridByRanges(blocks[l].range_blocks, cut_psi + grid_size * 2, l_grid, l_dpsi_dxyz + local_size * (6 * n));
                    CutSubgridByRanges(blocks[l].range_blocks, cut_psi + grid_size * 3, l_grid, l_dpsi_dxyz + local_size * (6 * n + 1));
                    CutSubgridByRanges(blocks[l].range_blocks, cut_psi + grid_size * 4, l_grid, l_dpsi_dxyz + local_size * (6 * n + 2));
                    CutSubgridByRanges(blocks[l].range_blocks, cut_psi + grid_size * 5, l_grid, l_dpsi_dxyz + local_size * (6 * n + 3));
                    CutSubgridByRanges(blocks[l].range_blocks, cut_psi + grid_size * 6, l_grid, l_dpsi_dxyz + local_size * (6 * n + 4));
                    CutSubgridByRanges(blocks[l].range_blocks, cut_psi + grid_size * 7, l_grid, l_dpsi_dxyz + local_size * (6 * n + 5));
                }

                //auto* psi_re = cut_psi[l].re;


                auto* proj_re = proj_p + offset;


                double* inner = new double[8 * width * 2 * num_psi];
                blas_DGEMM_t(8*num_psi, width * 2, grid_size, cut_psi_all, proj_re, inner, 8 * num_psi, dV, 0.0);


                #pragma ivdep
                for (int j = 0; j < width; ++j) {
                    #pragma ivdep
                    for (int m = 0; m < 4*num_psi; ++m) {
                        l_inner[index + j * 4 * num_psi + m].r = inner[16 * num_psi * j + m * 2 + 0] - inner[16 * num_psi * j + m * 2 + 8 * num_psi + 1];
                        l_inner[index + j * 4 * num_psi + m].i = inner[16 * num_psi * j + m * 2 + 1] + inner[16 * num_psi * j + m * 2 + 8 * num_psi + 0];

                    }
                }
                index += 4 * num_psi * width;
                offset += grid_size * width * 2;

                delete[] inner;
                //delete[] cut_psi[l].re;
            }


        }


        if (num_all_nonlocal * 4 * num_psi != index) {
            printf("ERROR: number of nonlocal is invalid"); fflush(stdout);
        }
        //printf("[%d] ProjectionPP 2: 0x%zx\n", GetProcessID(l_grid.mpi_comm), l_inner); fflush(stdout);
        //MPI_Barrier(l_grid.mpi_comm);
    }
#endif
#endif


    void ProjectionPP_gamma(SoAComplex& Hp, const SoAComplex& l_psi, const GridRangeMPI& l_grid, const Nucleus* nuclei, int num_nuclei,
        int id_spin_kpoint) {

        const int proc_id = GetProcessID(l_grid.mpi_comm);

        const double coef = (m_dx * m_dy * m_dz);
        const double dV = coef;
        const double Y_00 = 1.0 / sqrt(4.0 * M_PI);


        watch_pp.Restart();


        const int num_all_nonlocal = m_num_all_nonlocal;
        auto& info_range = m_num_list_nonlocal;

        watch_pp.Record(0);

        OneComplex* l_inner = new OneComplex[num_all_nonlocal * 2];
        OneComplex* sum_inner = l_inner + num_all_nonlocal;

        const int c_num_valid = mInnerNonlocalPsi_v2_gamma(l_psi, l_grid, nuclei, num_nuclei, l_inner, id_spin_kpoint, info_range, dV);
        //mInnerNonlocalPsi_v3_gamma(l_psi, l_grid, nuclei, num_nuclei, l_inner, id_spin_kpoint, info_range, dV);
        MPI_Waitall(c_num_valid, &m_request4atoms[0], &m_status4atoms[0]);
        watch_pp.Record(5);


        std::vector<SoAComplex> cut_Hp(MAX_L_SYSTEM);
        double* cut_psi_0 = new double[m_max_grid_size * MAX_L_SYSTEM * 2];
        for (int l = 0; l < MAX_L_SYSTEM; ++l) {
            cut_Hp[l].re = cut_psi_0 + m_max_grid_size * l * 2;
            cut_Hp[l].im = cut_psi_0 + m_max_grid_size * (l * 2 + 1);
        }



        for (int ni = 0; ni < num_nuclei; ++ni) {
            auto mycomm = m_comm4atoms.GetComm(ni);
            if (mycomm == MPI_COMM_NULL) continue;
            int index = info_range[ni];

            const int max_l = m_nonlocal_projectors[ni].max_l;
            //const auto& cut_range_block_l = m_nonlocal_block[ni].range_blocks;
            //const auto& nonlocal_size = m_nonlocal_block[ni].grid_sizes;
            const auto& blocks = m_nonlocal_blocks[ni];

            for (int l = 0; l <= max_l; ++l) {
                memset(cut_Hp[l].re, 0, sizeof(double) * blocks[l].grid_sizes);
                memset(cut_Hp[l].im, 0, sizeof(double) * blocks[l].grid_sizes);
            }

            watch_pp.Record(2);
            for (int k = 0; k < m_nonlocal_projectors[ni].num_projectors; ++k) {
                const int l = m_nonlocal_projectors[ni].quantum_l[k];
                const auto* proj_re = m_nonlocal_projectors[ni].projector_expikx[id_spin_kpoint][k].re;
                //const auto* proj_im = m_nonlocal_projectors[ni].projector_expikx[id_spin_kpoint][k].im;
                auto* cut_Hp_l_re = cut_Hp[l].re;
                auto* cut_Hp_l_im = cut_Hp[l].im;

                if (l == 0) {
                    OneComplex inner = sum_inner[index];
                    ++index;
                    inner.r *= Y_00 * m_nonlocal_projectors[ni].projector_energy[k];
                    inner.i *= Y_00 * m_nonlocal_projectors[ni].projector_energy[k];

                    const size_t nl_size = blocks[l].grid_sizes;
#pragma ivdep
                    for (int i = 0; i < nl_size; ++i) {
                        cut_Hp_l_re[i] += proj_re[i] * inner.r;
                        cut_Hp_l_im[i] += proj_re[i] * inner.i;
                    }


                } else {
                    for (int m1 = 0; m1 < 2 * l + 1; ++m1) {
                        const auto* Ylm = m_nonlocal_projectors[ni].Ylm[l * l + m1];
                        OneComplex inner = sum_inner[index];
                        ++index;
                        inner.r *= m_nonlocal_projectors[ni].projector_energy[k];
                        inner.i *= m_nonlocal_projectors[ni].projector_energy[k];

                        const size_t nl_size = blocks[l].grid_sizes;
#pragma ivdep
                        for (int i = 0; i < nl_size; ++i) {
                            const double YP = Ylm[i] * proj_re[i];
                            cut_Hp_l_re[i] += YP * inner.r;
                            cut_Hp_l_im[i] += YP * inner.i;
                        }
                    }
                }
            }
            watch_pp.Record(6);

            for (int l = 0; l <= max_l; ++l) {
                AddSubgridByRanges(l_grid, Hp.re, blocks[l].range_blocks, cut_Hp[l].re);
            }
            for (int l = 0; l <= max_l; ++l) {
                AddSubgridByRanges(l_grid, Hp.im, blocks[l].range_blocks, cut_Hp[l].im);
            }
            watch_pp.Record(7);

        }


        delete[] cut_psi_0;
        delete[] l_inner;
    }


    void ProjectionPP_gamma_s(SoAComplex& Hp, const SoAComplex& l_psi, const GridRangeMPI& l_grid, const Nucleus* nuclei, int num_nuclei,
        int id_spin_kpoint) {

        const int proc_id = GetProcessID(l_grid.mpi_comm);

        const double coef = (m_dx * m_dy * m_dz);
        const double dV = coef;
        const double Y_00 = 1.0 / sqrt(4.0 * M_PI);


        watch_pp.Restart();


        const int num_all_nonlocal = m_num_all_nonlocal;
        auto& info_range = m_num_list_nonlocal;

        watch_pp.Record(0);

        OneComplex* l_inner = new OneComplex[num_all_nonlocal * 2];
        OneComplex* sum_inner = l_inner + num_all_nonlocal;

        const int c_num_valid = mInnerNonlocalPsi_v2_gamma(l_psi, l_grid, nuclei, num_nuclei, l_inner, id_spin_kpoint, info_range, dV);
        //mInnerNonlocalPsi_v3_gamma(l_psi, l_grid, nuclei, num_nuclei, l_inner, id_spin_kpoint, info_range, dV);
        //MPI_Waitall(c_num_valid, &m_request4atoms[0], &m_status4atoms[0]);
        //watch_pp.Record(5);


        std::vector<SoAComplex> cut_Hp(MAX_L_SYSTEM);
        double* cut_psi_0 = new double[m_max_grid_size * MAX_L_SYSTEM * 2];
        for (int l = 0; l < MAX_L_SYSTEM; ++l) {
            cut_Hp[l].re = cut_psi_0 + m_max_grid_size * l * 2;
            cut_Hp[l].im = cut_psi_0 + m_max_grid_size * (l * 2 + 1);
        }

        int num_valid = 0;

        for (int ni = 0; ni < num_nuclei; ++ni) {
            auto mycomm = m_comm4atoms.GetComm(ni);
            if (mycomm == MPI_COMM_NULL) continue;
            int index = info_range[ni];

            const int max_l = m_nonlocal_projectors[ni].max_l;
            //const auto& cut_range_block_l = m_nonlocal_block[ni].range_blocks;
            //const auto& nonlocal_size = m_nonlocal_block[ni].grid_sizes;
            const auto& blocks = m_nonlocal_blocks[ni];

            for (int l = 0; l <= max_l; ++l) {
                memset(cut_Hp[l].re, 0, sizeof(double) * blocks[l].grid_sizes);
                memset(cut_Hp[l].im, 0, sizeof(double) * blocks[l].grid_sizes);
            }
            watch_pp.Record(2);
            
            MPI_Wait(&m_request4atoms[num_valid], &m_status4atoms[num_valid]);
            ++num_valid;
            watch_pp.Record(5);


            for (int k = 0; k < m_nonlocal_projectors[ni].num_projectors; ++k) {
                const int l = m_nonlocal_projectors[ni].quantum_l[k];
                const auto* proj_re = m_nonlocal_projectors[ni].projector_expikx[id_spin_kpoint][k].re;
                //const auto* proj_im = m_nonlocal_projectors[ni].projector_expikx[id_spin_kpoint][k].im;
                auto* cut_Hp_l_re = cut_Hp[l].re;
                auto* cut_Hp_l_im = cut_Hp[l].im;

                if (l == 0) {
                    OneComplex inner = sum_inner[index];
                    ++index;
                    inner.r *= Y_00 * m_nonlocal_projectors[ni].projector_energy[k];
                    inner.i *= Y_00 * m_nonlocal_projectors[ni].projector_energy[k];

                    const size_t nl_size = blocks[l].grid_sizes;
#pragma ivdep
                    for (int i = 0; i < nl_size; ++i) {
                        cut_Hp_l_re[i] += proj_re[i] * inner.r;
                        cut_Hp_l_im[i] += proj_re[i] * inner.i;
                    }


                } else {
                    for (int m1 = 0; m1 < 2 * l + 1; ++m1) {
                        const auto* Ylm = m_nonlocal_projectors[ni].Ylm[l * l + m1];
                        OneComplex inner = sum_inner[index];
                        ++index;
                        inner.r *= m_nonlocal_projectors[ni].projector_energy[k];
                        inner.i *= m_nonlocal_projectors[ni].projector_energy[k];

                        const size_t nl_size = blocks[l].grid_sizes;
#pragma ivdep
                        for (int i = 0; i < nl_size; ++i) {
                            const double YP = Ylm[i] * proj_re[i];
                            cut_Hp_l_re[i] += YP * inner.r;
                            cut_Hp_l_im[i] += YP * inner.i;
                        }
                    }
                }
            }
            watch_pp.Record(6);

            for (int l = 0; l <= max_l; ++l) {
                AddSubgridByRanges(l_grid, Hp.re, blocks[l].range_blocks, cut_Hp[l].re);
            }
            for (int l = 0; l <= max_l; ++l) {
                AddSubgridByRanges(l_grid, Hp.im, blocks[l].range_blocks, cut_Hp[l].im);
            }
            watch_pp.Record(7);

        }


        delete[] cut_psi_0;
        delete[] l_inner;
    }


    void ProjectionPP_v3_gamma(SoAComplex& Hp, const SoAComplex& l_psi, const GridRangeMPI& l_grid, const Nucleus* nuclei, int num_nuclei,
        int id_spin_kpoint) {

        const int proc_id = GetProcessID(l_grid.mpi_comm);

        const double coef = (m_dx * m_dy * m_dz);
        const double dV = coef;
        const double Y_00 = 1.0 / sqrt(4.0 * M_PI);


        watch_pp.Restart();


        const int num_all_nonlocal = m_num_all_nonlocal;
        auto& info_range = m_num_list_nonlocal;

        watch_pp.Record(0);

        OneComplex* l_inner = new OneComplex[num_all_nonlocal * 2];
        OneComplex* sum_inner = l_inner + num_all_nonlocal;

        //mInnerNonlocalPsi_v2_gamma(l_psi, l_grid, nuclei, num_nuclei, l_inner, id_spin_kpoint, info_range, dV);
        mInnerNonlocalPsi_v3_gamma(l_psi, l_grid, nuclei, num_nuclei, l_inner, id_spin_kpoint, info_range, dV);


        std::vector<SoAComplex> cut_Hp(MAX_L_SYSTEM);
        double* cut_psi_0 = new double[m_max_grid_size * MAX_L_SYSTEM * 2];
        for (int l = 0; l < MAX_L_SYSTEM; ++l) {
            cut_Hp[l].re = cut_psi_0 + m_max_grid_size * l * 2;
            cut_Hp[l].im = cut_psi_0 + m_max_grid_size * (l * 2 + 1);
        }



        for (int ni = 0; ni < num_nuclei; ++ni) {
            auto mycomm = m_comm4atoms.GetComm(ni);
            if (mycomm == MPI_COMM_NULL) continue;
            int index = info_range[ni];

            auto* proj_p = m_nonlocal_projectors[ni].projector_RY_expikx[id_spin_kpoint];
            int64_t offset = 0;

            const auto& blocks = m_nonlocal_blocks[ni];
            int k_offset = 0;
            const int max_l = m_nonlocal_projectors[ni].max_l;
            for (int l = 0; l <= max_l; ++l) {
                memset(cut_Hp[l].re, 0, sizeof(double) * blocks[l].grid_sizes);
                memset(cut_Hp[l].im, 0, sizeof(double) * blocks[l].grid_sizes);
            }
            watch_pp.Record(2);
            for (int l = 0; l <= max_l; ++l) {
                auto* cut_Hp_l_re = cut_Hp[l].re;
                auto* cut_Hp_l_im = cut_Hp[l].im;


                const size_t nl_size = blocks[l].grid_sizes;
                const int width = (m_nonlocal_projectors[ni].head_projector_l[l + 1] - m_nonlocal_projectors[ni].head_projector_l[l]);
                for (int w = 0; w < width; ++w) {
                    const int k = k_offset + w / (2*l + 1);

                    const auto* proj_re = proj_p + offset;
                    offset += nl_size;

                    OneComplex inner = sum_inner[index];
                    ++index;
                    inner.r *= m_nonlocal_projectors[ni].projector_energy[k];
                    inner.i *= m_nonlocal_projectors[ni].projector_energy[k];


                    #pragma ivdep
                    for (int i = 0; i < nl_size; ++i) {
                        cut_Hp_l_re[i] += proj_re[i] * inner.r;
                        cut_Hp_l_im[i] += proj_re[i] * inner.i;
                    }


                }

                k_offset += width / (2*l + 1);
            }
            watch_pp.Record(6);

            for (int l = 0; l <= max_l; ++l) {
                AddSubgridByRanges(l_grid, Hp.re, blocks[l].range_blocks, cut_Hp[l].re);
            }
            for (int l = 0; l <= max_l; ++l) {
                AddSubgridByRanges(l_grid, Hp.im, blocks[l].range_blocks, cut_Hp[l].im);
            }
            watch_pp.Record(7);
        }


        delete[] cut_psi_0;
        delete[] l_inner;
    }

    void ProjectionPP(SoAComplex& Hp, const SoAComplex& l_psi, const GridRangeMPI& l_grid, const Nucleus* nuclei, int num_nuclei,
        int id_spin_kpoint) {

#ifdef FAST_GAMMA_POINT
        if (m_is_gamma_point_list[id_spin_kpoint]) {
            //ProjectionPP_gamma(Hp, l_psi, l_grid, nuclei, num_nuclei, id_spin_kpoint);
            ProjectionPP_gamma_s(Hp, l_psi, l_grid, nuclei, num_nuclei, id_spin_kpoint);
            //ProjectionPP_v3_gamma(Hp, l_psi, l_grid, nuclei, num_nuclei, id_spin_kpoint);
            return;
        }
#endif

        const int proc_id = GetProcessID(l_grid.mpi_comm);

        const double coef = (m_dx * m_dy * m_dz);
        const double dV = coef;
        const double Y_00 = 1.0 / sqrt(4.0 * M_PI);


        watch_pp.Restart();


        const int num_all_nonlocal = m_num_all_nonlocal;
        auto& info_range = m_num_list_nonlocal;
        
        watch_pp.Record(0);

        OneComplex* l_inner = new OneComplex[num_all_nonlocal * 2];
        OneComplex* sum_inner = l_inner + num_all_nonlocal;

        mInnerNonlocalPsi_v2(l_psi, l_grid, nuclei, num_nuclei, l_inner, id_spin_kpoint, info_range, dV);
        //mInnerNonlocalPsi_v2_1(l_psi, l_grid, nuclei, num_nuclei, l_inner, id_spin_kpoint, info_range, dV);
        //mInnerNonlocalPsi_v2_2(l_psi, l_grid, nuclei, num_nuclei, l_inner, id_spin_kpoint, info_range, dV);
        //mInnerNonlocalPsi_v3(l_psi, l_grid, nuclei, num_nuclei, l_inner, id_spin_kpoint, info_range, dV);


        std::vector<SoAComplex> cut_Hp(MAX_L_SYSTEM);
        double* cut_psi_0 = new double[m_max_grid_size * MAX_L_SYSTEM * 2];
        for (int l = 0; l < MAX_L_SYSTEM; ++l) {
            cut_Hp[l].re = cut_psi_0 + m_max_grid_size * l * 2;
            cut_Hp[l].im = cut_psi_0 + m_max_grid_size * (l * 2 + 1);
        }



        for (int ni = 0; ni < num_nuclei; ++ni) {
            auto mycomm = m_comm4atoms.GetComm(ni);
            if (mycomm == MPI_COMM_NULL) continue;
            int index = info_range[ni];

            const int max_l = m_nonlocal_projectors[ni].max_l;
            //const auto& cut_range_block_l = m_nonlocal_block[ni].range_blocks;
            //const auto& nonlocal_size = m_nonlocal_block[ni].grid_sizes;
            const auto& blocks = m_nonlocal_blocks[ni];

            for (int l = 0; l <= max_l; ++l) {
                memset(cut_Hp[l].re, 0, sizeof(double) * blocks[l].grid_sizes);
                memset(cut_Hp[l].im, 0, sizeof(double) * blocks[l].grid_sizes);
            }

            watch_pp.Record(2);
            for (int k = 0; k < m_nonlocal_projectors[ni].num_projectors; ++k) {
                const int l = m_nonlocal_projectors[ni].quantum_l[k];
                const auto* proj_re = m_nonlocal_projectors[ni].projector_expikx[id_spin_kpoint][k].re;
                const auto* proj_im = m_nonlocal_projectors[ni].projector_expikx[id_spin_kpoint][k].im;
                auto* cut_Hp_l_re = cut_Hp[l].re;
                auto* cut_Hp_l_im = cut_Hp[l].im;

                if (l == 0) {
                    OneComplex inner = sum_inner[index];
                    ++index;
                    inner.r *= Y_00 * m_nonlocal_projectors[ni].projector_energy[k];
                    inner.i *= Y_00 * m_nonlocal_projectors[ni].projector_energy[k];

                    const size_t nl_size = blocks[l].grid_sizes;
                    //slow if use #pragma ivdep
                    for (int i = 0; i < nl_size; ++i) {
                        cut_Hp_l_re[i] += proj_re[i] * inner.r + proj_im[i] * inner.i;
                        cut_Hp_l_im[i] += -proj_im[i] * inner.r + proj_re[i] * inner.i;
                    }


                } else {
                    const size_t nl_size = blocks[l].grid_sizes;
                    for (int m1 = 0; m1 < 2 * l + 1; ++m1) {
                        const auto* Ylm = m_nonlocal_projectors[ni].Ylm[l * l + m1];
                        OneComplex inner = sum_inner[index];
                        ++index;
                        inner.r *= m_nonlocal_projectors[ni].projector_energy[k];
                        inner.i *= m_nonlocal_projectors[ni].projector_energy[k];

                        //slow if use #pragma ivdep
                        for (int i = 0; i < nl_size; ++i) {
                            cut_Hp_l_re[i] += Ylm[i] * (proj_re[i] * inner.r + proj_im[i] * inner.i);
                            cut_Hp_l_im[i] += Ylm[i] * (-proj_im[i] * inner.r + proj_re[i] * inner.i);
                        }
                    }
                }
            }
            watch_pp.Record(6);

            for (int l = 0; l <= max_l; ++l) {
                AddSubgridByRanges(l_grid, Hp.re, blocks[l].range_blocks, cut_Hp[l].re);
                AddSubgridByRanges(l_grid, Hp.im, blocks[l].range_blocks, cut_Hp[l].im);
            }
            watch_pp.Record(7);

        }

        delete[] cut_psi_0;
        delete[] l_inner;
    }



    double EnergyPPnonlocal(const SoAComplex& l_psi, const GridRangeMPI& l_grid, const Nucleus* nuclei, int num_nuclei, int id_spin_kpoint) {

        const double coef = (m_dx * m_dy * m_dz);
        const double dV = coef;


        const int num_all_nonlocal = m_num_all_nonlocal;
      

        OneComplex* l_inner = new OneComplex[num_all_nonlocal];


        mInnerNonlocalPsi(l_psi, l_grid, nuclei, num_nuclei, l_inner, id_spin_kpoint, num_all_nonlocal, dV);



        OneComplex* sum_inner = IsRoot(l_grid.mpi_comm) ? new OneComplex[num_all_nonlocal] : nullptr;
        MPI_Reduce(l_inner, sum_inner, num_all_nonlocal * 2, MPI_DOUBLE, MPI_SUM, 0, l_grid.mpi_comm);
        delete[]l_inner;


        double ene = 0.0;
        if (IsRoot(l_grid.mpi_comm)) {

            int index = 0;
            for (int ni = 0; ni < num_nuclei; ++ni) {

                for (int k = 0; k < m_nonlocal_projectors[ni].num_projectors; ++k) {
                    const int l = m_nonlocal_projectors[ni].quantum_l[k];

                    if (l == 0) {

                        const auto& inner = sum_inner[index];
                        ++index;
                        const double ee = (inner.r * inner.r + inner.i * inner.i) * m_nonlocal_projectors[ni].projector_energy[k];
                        ene += ee;
                        //#define DEBUG_PRINT1
#ifdef DEBUG_PRINT1
                        printf("E_nonlocal(%d,0) = %f, %f\n", l, ee, inner);
#endif
                    } else {
                        for (int m1 = 0; m1 < 2 * l + 1; ++m1) {
                            const auto& inner = sum_inner[index];
                            ++index;
                            const double ee = (inner.r * inner.r + inner.i * inner.i) * m_nonlocal_projectors[ni].projector_energy[k];
                            ene += ee;
#ifdef DEBUG_PRINT1
                            printf("E_nonlocal(%d,%d) = %f, %f\n", l, m1, ee, inner);
#endif
                            //#undef DEBUG_PRINT1
                        }
                    }
                }
            }

            delete[]sum_inner;
        }


        return ene;
    }



    double ForceNonlocal(double* forces, const SoAComplex& l_psi, const SoAComplex& l_dpsi_dx, const SoAComplex& l_dpsi_dy, const SoAComplex& l_dpsi_dz, const GridRangeMPI& l_grid, const Nucleus* nuclei, int num_nuclei,
        int id_spin_kpoint, double gx, double gy, double gz) {

        

        watch_pp.Restart();

        const int proc_id = GetProcessID(l_grid.mpi_comm);
        constexpr int stride = 3;
        const int local_size = l_grid.Size3D();

        const double dV = m_dx * m_dy * m_dz;

        const int num_all_nonlocal = m_num_all_nonlocal;
     
        OneComplex* l_inner = new OneComplex[num_all_nonlocal * 4];

#ifdef EXPAND_PROJECTOR_YLM_Z
        
#ifdef FAST_GAMMA_POINT
        if (m_is_gamma_point_list[id_spin_kpoint]) {
            mInnerNonlocalPsi_b3_gamma(l_psi, l_dpsi_dx, l_dpsi_dy, l_dpsi_dz, l_grid, nuclei, num_nuclei, l_inner, id_spin_kpoint, num_all_nonlocal, dV);
        } 
        else 
#endif

        {
            mInnerNonlocalPsi_b3(l_psi, l_dpsi_dx, l_dpsi_dy, l_dpsi_dz, l_grid, nuclei, num_nuclei, l_inner, id_spin_kpoint, num_all_nonlocal, dV);
            //mInnerNonlocalPsi_b2(l_psi, l_dpsi_dx, l_dpsi_dy, l_dpsi_dz, l_grid, nuclei, num_nuclei, l_inner, id_spin_kpoint, num_all_nonlocal, dV);
        }
#else
        mInnerNonlocalPsi(l_psi, l_grid, nuclei, num_nuclei, l_inner, id_spin_kpoint, num_all_nonlocal, dV);
        mInnerNonlocalPsi(l_dpsi_dx, l_grid, nuclei, num_nuclei, l_inner + num_all_nonlocal, id_spin_kpoint, num_all_nonlocal, dV);
        mInnerNonlocalPsi(l_dpsi_dy, l_grid, nuclei, num_nuclei, l_inner + num_all_nonlocal * 2, id_spin_kpoint, num_all_nonlocal, dV);
        mInnerNonlocalPsi(l_dpsi_dz, l_grid, nuclei, num_nuclei, l_inner + num_all_nonlocal * 3, id_spin_kpoint, num_all_nonlocal, dV);
#endif


        OneComplex* sum_inner = IsRoot(l_grid.mpi_comm) ? new OneComplex[num_all_nonlocal * 4] : nullptr;
        MPI_Reduce(l_inner, sum_inner, num_all_nonlocal * 8, MPI_DOUBLE, MPI_SUM, 0, l_grid.mpi_comm);
        delete[] l_inner;



        double ene = 0.0;
        if (IsRoot(l_grid.mpi_comm)) {
            watch_pp.Record(8);

#ifdef EXPAND_PROJECTOR_YLM_Z
#else
            const OneComplex* sum_diff_x = sum_inner + num_all_nonlocal;
            const OneComplex* sum_diff_y = sum_inner + num_all_nonlocal * 2;
            const OneComplex* sum_diff_z = sum_inner + num_all_nonlocal * 3;
#endif
            int index = 0;
            for (int ni = 0; ni < num_nuclei; ++ni) {

                double force_i[3] = { 0.0, 0.0, 0.0 };

                for (int k = 0; k < m_nonlocal_projectors[ni].num_projectors; ++k) {
                    const int l = m_nonlocal_projectors[ni].quantum_l[k];


                    for (int m1 = 0; m1 < 2 * l + 1; ++m1) {

                        const OneComplex& inner = sum_inner[index];

#ifdef EXPAND_PROJECTOR_YLM_Z
                        OneComplex diff_x = sum_inner[index + 1];
                        OneComplex diff_y = sum_inner[index + 2];
                        OneComplex diff_z = sum_inner[index + 3];
                        index += 4;

#else
                        OneComplex diff_x = sum_diff_x[index];
                        OneComplex diff_y = sum_diff_y[index];
                        OneComplex diff_z = sum_diff_z[index];
                        ++index;
#endif
                        diff_x.i += gx * inner.r;
                        diff_x.r += -gx * inner.i;
                        diff_y.i += gy * inner.r;
                        diff_y.r += -gy * inner.i;
                        diff_z.i += gz * inner.r;
                        diff_z.r += -gz * inner.i;

                        const double ee = (inner.r * inner.r + inner.i * inner.i) * m_nonlocal_projectors[ni].projector_energy[k];
                        ene += ee;
                        force_i[0] += (inner.r * diff_x.r + inner.i * diff_x.i) * m_nonlocal_projectors[ni].projector_energy[k];
                        force_i[1] += (inner.r * diff_y.r + inner.i * diff_y.i) * m_nonlocal_projectors[ni].projector_energy[k];
                        force_i[2] += (inner.r * diff_z.r + inner.i * diff_z.i) * m_nonlocal_projectors[ni].projector_energy[k];
                        //NOTE: imaginary part should be vanish to add c.c. part//
                    }
                }

                forces[ni * stride + 0] += force_i[0];
                forces[ni * stride + 1] += force_i[1];
                forces[ni * stride + 2] += force_i[2];
            }

            delete[]sum_inner;
            watch_pp.Record(9);
        }


        return ene;
    }

#if 0
    //slower than original ForceNonlocal that uses just one psi and its diferentials //
    double ForceNonlocal2(double* forces, int num_psi, const SoAComplex* l_psi, const double* l_dpsi_dxyz,
        double* work,
        const double* occupancy, const GridRangeMPI& l_grid, const Nucleus* nuclei, int num_nuclei,
        int id_spin_kpoint, double gx, double gy, double gz) {

        watch_pp.Restart();

        const int proc_id = GetProcessID(l_grid.mpi_comm);
        constexpr int stride = 3;
        const int local_size = l_grid.Size3D();

        const double dV = m_dx * m_dy * m_dz;

        const int num_all_nonlocal = m_num_all_nonlocal;

        OneComplex* l_inner = new OneComplex[num_all_nonlocal * 4 * num_psi];

        mInnerNonlocalPsi_b4(num_psi, l_psi, l_dpsi_dxyz, work, l_grid, nuclei, num_nuclei, l_inner, id_spin_kpoint, num_all_nonlocal, dV);
        



        OneComplex* sum_inner = IsRoot(l_grid.mpi_comm) ? new OneComplex[num_all_nonlocal * 4 * num_psi] : nullptr;
        MPI_Reduce(l_inner, sum_inner, num_all_nonlocal * 8 * num_psi, MPI_DOUBLE, MPI_SUM, 0, l_grid.mpi_comm);
        delete[] l_inner;



        double ene = 0.0;
        if (IsRoot(l_grid.mpi_comm)) {
            watch_pp.Record(8);


            int index = 0;
            for (int ni = 0; ni < num_nuclei; ++ni) {

                double force_i[3] = { 0.0, 0.0, 0.0 };

                for (int k = 0; k < m_nonlocal_projectors[ni].num_projectors; ++k) {
                    const int l = m_nonlocal_projectors[ni].quantum_l[k];


                    for (int m1 = 0; m1 < 2 * l + 1; ++m1) {

                        for (int n = 0; n < num_psi; ++n) {
                            const OneComplex& inner = sum_inner[index];

#ifdef EXPAND_PROJECTOR_YLM_Z
                            OneComplex diff_x = sum_inner[index + 1];
                            OneComplex diff_y = sum_inner[index + 2];
                            OneComplex diff_z = sum_inner[index + 3];
                            index += 4;

#else
                            OneComplex diff_x = sum_diff_x[index];
                            OneComplex diff_y = sum_diff_y[index];
                            OneComplex diff_z = sum_diff_z[index];
                            ++index;
#endif
                            diff_x.i += gx * occupancy[n] * inner.r;
                            diff_x.r += -gx * occupancy[n] * inner.i;
                            diff_y.i += gy * occupancy[n] * inner.r;
                            diff_y.r += -gy * occupancy[n] * inner.i;
                            diff_z.i += gz * occupancy[n] * inner.r;
                            diff_z.r += -gz * occupancy[n] * inner.i;

                            const double ee = (inner.r * inner.r + inner.i * inner.i) * m_nonlocal_projectors[ni].projector_energy[k];
                            ene += ee * occupancy[n];
                            force_i[0] += (inner.r * diff_x.r + inner.i * diff_x.i) * m_nonlocal_projectors[ni].projector_energy[k];
                            force_i[1] += (inner.r * diff_y.r + inner.i * diff_y.i) * m_nonlocal_projectors[ni].projector_energy[k];
                            force_i[2] += (inner.r * diff_z.r + inner.i * diff_z.i) * m_nonlocal_projectors[ni].projector_energy[k];
                            //NOTE: imaginary part should be vanish to add c.c. part//
                        }
                    }
                }

                forces[ni * stride + 0] += force_i[0];
                forces[ni * stride + 1] += force_i[1];
                forces[ni * stride + 2] += force_i[2];
            }

            delete[]sum_inner;
            watch_pp.Record(9);
        }


        return ene;
    }
#endif

    void mPrintTime() {
#ifdef TIME_PP_MPI
        //if (watch_pp.GetTime(4) > watch_pp.GetTime(5)) 
        {

            printf("PP calculation times====\n");
            watch_pp.Print("Count_num", 0);
            watch_pp.Print("Find_pp  ", 1);
            watch_pp.Print("Allocate ", 2);
            watch_pp.Print("Cut_psi  ", 3);
            watch_pp.Print("Inner    ", 4);
            watch_pp.Print("MPI      ", 5);
            watch_pp.Print("Add local", 6);
            watch_pp.Print("Paste    ", 7);

            watch_pp.Print("ForceNonlocal1", 8);
            watch_pp.Print("ForceNonlocal2", 9);
        }
#endif
    }



};

#endif
