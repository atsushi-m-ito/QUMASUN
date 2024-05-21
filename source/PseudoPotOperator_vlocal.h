#ifdef USE_MPI
#pragma once
#include <mpi.h>
#include "mpi_helper.h"
#include "PseudoPotOperator.h"
#include "vps_mpi_broadcast.h"
#include "GridSubgrid_arithmetic.h"
#include "soacomplex.h"
#include "comm4atom.h"
#include "SetBlockYlm.h"

//この定義を有効にすると遠方のポアソンを解いたVlocalの遠方での値はばらつかなくなるが
//格子定数を変化させたときのエネルギーカーブはガタつく
//#define CORE_CHARGE_SCALING



/********************************
* note: block化
* 擬ポテンシャルの存在する部分領域と、
* 系全体の領域分割によって本mpi-processが保持する波動関数の領域の
* オーバーラップ部分だけに相当する領域のデータをblock化されたデータと呼ぶことにする.
* 周期境界の為に、1つの擬ポテンシャルが複数のブロックで構成されることがあり得るが,
* 配列内では連続的につないでしまって一つの配列として保持する.
* 内積を取る演算時も一続きにアクセスして良い.
*****************************/

using NonlocalBlocks = std::vector< SubgridBlock>;

/*
* VlocalおよびPCC-charge関連の基本的な処理
* また、原子ごとのMPI通信も含む
* このクラスを継承してnonlocalを含む実数関数版(mpi版)および複素関数版(kpoint版)を作成
*/
class PseudoPotIntegrator_vlocal : public PseudoPotIntegrator {
protected:
    std::vector< SubgridBlock> m_vlocal_block;//block化されたVlocalと対応する電荷//
    std::vector< SubgridBlock> m_hr_vlocal_block;//block化されたVlocalと対応する電荷, high resolution//
    std::vector< NonlocalBlocks> m_nonlocal_blocks;
    std::vector< double*> m_rho_vlocal;
    std::vector< double*> m_hr_rho_vlocal;
    std::vector< double*> m_pcc_charge;
    std::vector< std::vector<GridRange>*> m_pcc_range_blocks; //just pointer of corresponding m_nonlocal_blocks or m_vlocal_block//
    

    static constexpr int MAX_L_SYSTEM = 4;
    size_t m_max_grid_size = 0;
    size_t m_hr_max_grid_size = 0;

    //High resolution grid for nuclear charge//
    int m_HR_ratio_x = 1;
    int m_HR_ratio_y = 1;
    int m_HR_ratio_z = 1;


    CommForAtoms m_comm4atoms;
    int m_num_valid_comm4atoms = 0;
    std::vector<MPI_Request> m_request4atoms;
    std::vector<MPI_Status> m_status4atoms;

    double m_E_nn_self[2]{ 0.0,0.0 };

    bool is_root = false;

public:
    virtual ~PseudoPotIntegrator_vlocal() {
        mDeleteBuffers(m_rho_vlocal);
        mDeleteBuffers(m_pcc_charge);
    }

protected:
    void mDeleteBuffers(std::vector<double*>& rho_vlocal) {
        for (auto& v : rho_vlocal) {
            delete[] v;
        }
    }

public:
    //mpi並列時に全データで呼ばれる関数
    void Load(int Z, const char* filepath, MPI_Comm& mpi_comm) {
        PseudoPot_MBK* pp = nullptr;
        if (IsRoot(mpi_comm)) {
            pp = LoadVPS(filepath);
        } else {
            pp = new PseudoPot_MBK;
        }

        //broadcast data//
        BroadcastPseudoPot(pp, mpi_comm, 0);

        //register//
        m_pp_list.emplace(Z, pp);

        is_root = IsRoot(mpi_comm);
    }

protected:
    size_t mSetSubspaceBlock(SubgridBlock& vlocal_block, SubgridBlock& hr_vlocal_block, NonlocalBlocks& nonlocal_blocks,
        const Nucleus nucleus, const GridRange& l_grid) {

        size_t max_grid_size = 0;

        const PseudoPot_MBK* pp = mFindPseudoPot(nucleus.Z);
        if (pp == nullptr) {
            printf("ERROR: VPS is not loaded: Z = %d\n", nucleus.Z);
            return 0;
        }

        //Vlocalのブロック化//
        {
            const double cutoff_r = pp->cutoff_vlocal;
            vlocal_block.cutoff = cutoff_r;
            GridRange subgrid = GridForSphere(nucleus.Rx, nucleus.Ry, nucleus.Rz, cutoff_r, m_dx, m_dy, m_dz);

            const int total_block_size = GetOverlapPeriodic(l_grid, subgrid, m_grid.size_x, m_grid.size_y, m_grid.size_z, vlocal_block.range_blocks, vlocal_block.shifted_grids);
            vlocal_block.grid_sizes = total_block_size;
            if (max_grid_size < total_block_size) {
                max_grid_size = total_block_size;
            }
#ifdef DEBUG_PRINT
            printf("block-size = %zd, %f\n", total_block_size, cutoff_r); fflush(stdout);
#endif

        }

        //HR版, Vlocalのブロック化//
        {
            const double hr_dx = m_dx / (double)m_HR_ratio_x;
            const double hr_dy = m_dy / (double)m_HR_ratio_y;
            const double hr_dz = m_dz / (double)m_HR_ratio_z;
            const double cutoff_r = pp->cutoff_vlocal;
            hr_vlocal_block.cutoff = cutoff_r;
            GridRange subgrid = GridForSphere(nucleus.Rx, nucleus.Ry, nucleus.Rz, cutoff_r, hr_dx, hr_dy, hr_dz);
            auto hr_grid = ScaleRange(l_grid, m_HR_ratio_x, m_HR_ratio_y, m_HR_ratio_z);
            const int total_block_size = GetOverlapPeriodic(hr_grid, subgrid, m_grid.size_x* m_HR_ratio_x, m_grid.size_y* m_HR_ratio_y, m_grid.size_z* m_HR_ratio_z, hr_vlocal_block.range_blocks, hr_vlocal_block.shifted_grids);
            hr_vlocal_block.grid_sizes = total_block_size;
            /*
            if (max_grid_size < total_block_size) {
                max_grid_size = total_block_size;
            }
#ifdef DEBUG_PRINT
            printf("block-size = %zd, %f\n", total_block_size, cutoff_r); fflush(stdout);
#endif
*/

        }

        const int num_projectors = pp->num_projectors;
        if (num_projectors == 0) {
            return 0;
        }

        const int MAX_L = pp->MaxProjectorL();
        nonlocal_blocks.resize(MAX_L + 1);

        for (int l = 0; l <= MAX_L; ++l) {
            const double cutoff_r = pp->cutoff_r[l];
            nonlocal_blocks[l].cutoff = cutoff_r;

            GridRange subgrid_Ylm = GridForSphere(nucleus.Rx, nucleus.Ry, nucleus.Rz, cutoff_r, m_dx, m_dy, m_dz);

            const int total_block_size = GetOverlapPeriodic(l_grid, subgrid_Ylm, m_grid.size_x, m_grid.size_y, m_grid.size_z, nonlocal_blocks[l].range_blocks, nonlocal_blocks[l].shifted_grids);
            nonlocal_blocks[l].grid_sizes = total_block_size;
            if (max_grid_size < total_block_size) {
                max_grid_size = total_block_size;
            }
        }

        return max_grid_size;
    }

    /*
    * memory確保しているので再度呼び出す前には一度deleteが必要
    * mpi版とkpoint版でほぼ一緒なので共通化すべし
    */
    void mCreateChargeVlocalOnBlock(SubgridBlock& vlocal_block, double** p_rho_vlocal,
        const Nucleus nucleus, const GridRange& l_grid) {


        if (vlocal_block.grid_sizes == 0) {
            *p_rho_vlocal = nullptr;
            return;
        }

        const PseudoPot_MBK* pp = mFindPseudoPot(nucleus.Z);

        {

            const size_t total_block_size = vlocal_block.grid_sizes;
#ifdef DEBUG_PRINT
            printf("block-size2[l=%d] = %zd\n", l, total_block_size); fflush(stdout);
#endif

            double* rho_block = new double[total_block_size];
            *p_rho_vlocal = rho_block;

            const double rr_min = (pp->radius[0]) * (pp->radius[0]);
            const double rr_max = (vlocal_block.cutoff) * (vlocal_block.cutoff);
            const auto& range_blocks = vlocal_block.range_blocks;
            const auto& shifted_grid = vlocal_block.shifted_grids;
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


                            if (rr_max < rr) {
                                rho_block[i] = 0.0;
                            } else {
                                double r1;
                                if (rr_min > rr) {
                                    r1 = sqrt(rr_min + 1.0e-14);
                                } else {
                                    r1 = sqrt(rr);
                                }

                                rho_block[i] = (1.0 / (2.0 * M_PI)) * RadialGrid2::GetLaplacian(r1, pp->V_local, 0, pp->num_radial_grids, pp->xi_min, pp->xi_delta);
                            }

                        }
                    }
                }
                offset_i += range_blocks[ib].Size3D();

            }

        }

    }


    /*
    * memory確保しているので再度呼び出す前には一度deleteが必要
    * mpi版とkpoint版でほぼ一緒なので共通化すべし
    * 二階微分を実空間直行グリッドに焼き直した後に行う
    */
    double mCreateChargeVlocalOnBlock_v2(int ni, SubgridBlock& vlocal_block, double** p_rho_vlocal,
                const double dx, const double dy, const double dz, const Nucleus& nucleus, const GridRange& l_grid) {

        const double coef = -1.0 / (4.0 * M_PI);
        const double c0 = -14350.0 / 5040.0 * coef * (1.0 / (dx * dx) + 1.0 / (dy * dy) + 1.0 / (dz * dz));
        const double c4x = -9.0 * coef / (5040.0 * dx * dx);
        const double c3x = 128.0 * coef / (5040.0 * dx * dx);
        const double c2x = -1008.0 * coef / (5040.0 * dx * dx);
        const double c1x = 8064.0 * coef / (5040.0 * dx * dx);
        const double c4y = -9.0 * coef / (5040.0 * dy * dy);
        const double c3y = 128.0 * coef / (5040.0 * dy * dy);
        const double c2y = -1008.0 * coef / (5040.0 * dy * dy);
        const double c1y = 8064.0 * coef / (5040.0 * dy * dy);
        const double c4z = -9.0 * coef / (5040.0 * dz * dz);
        const double c3z = 128.0 * coef / (5040.0 * dz * dz);
        const double c2z = -1008.0 * coef / (5040.0 * dz * dz);
        const double c1z = 8064.0 * coef / (5040.0 * dz * dz);

        auto&& ni_comm = m_comm4atoms.GetComm(ni);

        const PseudoPot_MBK* pp = mFindPseudoPot(nucleus.Z);
        double valence_electron = pp->valence_electron;

        if (vlocal_block.grid_sizes == 0) {
            *p_rho_vlocal = nullptr;

            if (ni_comm == MPI_COMM_NULL) {
                return 0.0;
            }

            double l_values[2]{ 0.0,0.0 };
            double sum_values[2]{ 0.0,0.0 };
            double& sum_int_rho = sum_values[0];
            double& sum_int_V_rho = sum_values[1];
#ifdef CORE_CHARGE_SCALING
            MPI_Allreduce(l_values, sum_values, 2, MPI_DOUBLE, MPI_SUM, ni_comm);
            sum_int_rho *= dx * dy * dz;
            sum_int_V_rho *= dx * dy * dz;
            const double scale_ratio = valence_electron / sum_int_rho;
            sum_int_V_rho *= scale_ratio * scale_ratio;
#else
            MPI_Reduce(l_values, sum_values, 2, MPI_DOUBLE, MPI_SUM, 0, ni_comm);
#endif
            if (IsRoot(ni_comm)) {
                return sum_int_V_rho;
            }
            return 0.0;

        }


        double l_values[2]{ 0.0,0.0 };
        double sum_values[2]{ 0.0,0.0 };
        double& l_int_rho = l_values[0];
        double& l_int_V_rho = l_values[1];
        double& sum_int_rho = sum_values[0];
        double& sum_int_V_rho = sum_values[1];
        {

            const size_t total_block_size = vlocal_block.grid_sizes;
#ifdef DEBUG_PRINT
            printf("block-size2[l=%d] = %zd\n", l, total_block_size); fflush(stdout);
#endif

            double* rho_block = new double[total_block_size];
            *p_rho_vlocal = rho_block;

            if (ni_comm == MPI_COMM_NULL) {
                for (size_t i = 0; i < total_block_size; ++i) {
                    rho_block[i] = 0.0;
                }
                return 0.0;
            }

            std::vector<double> buf_Vlocal(total_block_size * 4);

            const double rr_min = (pp->radius[0]) * (pp->radius[0]);
            const double rr_max = (vlocal_block.cutoff) * (vlocal_block.cutoff);
            const auto& range_blocks = vlocal_block.range_blocks;
            const auto& shifted_grid = vlocal_block.shifted_grids;
            size_t num_blocks = range_blocks.size();
            int offset_i = 0;
            for (size_t ib = 0; ib < num_blocks; ++ib) {
                const int ox_begin = range_blocks[ib].begin_x;
                const int oy_begin = range_blocks[ib].begin_y;
                const int oz_begin = range_blocks[ib].begin_z;
                const int ox_end = range_blocks[ib].end_x;
                const int oy_end = range_blocks[ib].end_y;
                const int oz_end = range_blocks[ib].end_z;

                const int ix_begin = ox_begin - 4;
                const int iy_begin = oy_begin - 4;
                const int iz_begin = oz_begin - 4;
                const int ix_end = ox_end + 4;
                const int iy_end = oy_end + 4;
                const int iz_end = oz_end + 4;

                const int size_x = ix_end - ix_begin;
                const int size_y = iy_end - iy_begin;

                const int size_ox = ox_end - ox_begin;
                const int size_oy = oy_end - oy_begin;
                buf_Vlocal.resize(size_x * size_y * (iz_end - iz_begin));

                for (int iz = iz_begin; iz < iz_end; ++iz) {
                    const double z = dz * (double)(iz - shifted_grid[ib].z) - nucleus.Rz;
                    const double zz = z * z;
                    for (int iy = iy_begin; iy < iy_end; ++iy) {
                        const double y = dy * (double)(iy - shifted_grid[ib].y) - nucleus.Ry;
                        const double yy_zz = y * y + zz;
                        for (int ix = ix_begin; ix < ix_end; ++ix) {
                            const double x = dx * (double)(ix - shifted_grid[ib].x) - nucleus.Rx;
                            const double rr = x * x + yy_zz;
                            const int i = (ix - ix_begin) + size_x * ((iy - iy_begin) + size_y * (iz - iz_begin));

                            if (rr_min > rr) {
                                buf_Vlocal[i] = pp->V_local[0];
                            } else {
                                buf_Vlocal[i] = RadialGrid2::GetValueBySquare(rr, pp->V_local, pp->num_radial_grids, pp->xi_min, pp->xi_delta);
                            }
                        }
                    }
                }

                auto& p = buf_Vlocal;


                for (int iz = oz_begin; iz < oz_end; ++iz) {
                    const double z = dz * (double)(iz - shifted_grid[ib].z) - nucleus.Rz;
                    const double zz = z * z;

                    const int idz_m4 = -4 * size_x * size_y;
                    const int idz_m3 = -3 * size_x * size_y;
                    const int idz_m2 = -2 * size_x * size_y;
                    const int idz_m = -1 * size_x * size_y;
                    const int idz_p = 1 * size_x * size_y;
                    const int idz_p2 = 2 * size_x * size_y;
                    const int idz_p3 = 3 * size_x * size_y;
                    const int idz_p4 = 4 * size_x * size_y;

                    for (int iy = oy_begin; iy < oy_end; ++iy) {
                        const double y = dy * (double)(iy - shifted_grid[ib].y) - nucleus.Ry;
                        const double yy_zz = y * y + zz;

                        const int idy_m4 = -4 * size_x;
                        const int idy_m3 = -3 * size_x;
                        const int idy_m2 = -2 * size_x;
                        const int idy_m = -1 * size_x;
                        const int idy_p = 1 * size_x;
                        const int idy_p2 = 2 * size_x;
                        const int idy_p3 = 3 * size_x;
                        const int idy_p4 = 4 * size_x;

                        for (int ix = ox_begin; ix < ox_end; ++ix) {
                            const double x = dx * (double)(ix - shifted_grid[ib].x) - nucleus.Rx;
                            const double rr = x * x + yy_zz;

                            const int idx_m4 = -4;
                            const int idx_m3 = -3;
                            const int idx_m2 = -2;
                            const int idx_m = -1;
                            const int idx_p = 1;
                            const int idx_p2 = 2;
                            const int idx_p3 = 3;
                            const int idx_p4 = 4;

                            const int i = (ix - ix_begin) + size_x * ((iy - iy_begin) + size_y * (iz - iz_begin));
                            const int oi = offset_i + (ix - ox_begin) + size_ox * ((iy - oy_begin) + size_oy * (iz - oz_begin));


                            if (rr_max < rr) {
                                rho_block[oi] = 0.0;
                            } else {
                                auto d2psidx2 = c0 * p[i];
                                d2psidx2 += c4x * (p[i + idx_p4] + p[i + idx_m4]) + c3x * (p[i + idx_p3] + p[i + idx_m3]) + c2x * (p[i + idx_p2] + p[i + idx_m2]) + c1x * (p[i + idx_p] + p[i + idx_m]);
                                d2psidx2 += c4y * (p[i + idy_p4] + p[i + idy_m4]) + c3y * (p[i + idy_p3] + p[i + idy_m3]) + c2y * (p[i + idy_p2] + p[i + idy_m2]) + c1y * (p[i + idy_p] + p[i + idy_m]);
                                d2psidx2 += c4z * (p[i + idz_p4] + p[i + idz_m4]) + c3z * (p[i + idz_p3] + p[i + idz_m3]) + c2z * (p[i + idz_p2] + p[i + idz_m2]) + c1z * (p[i + idz_p] + p[i + idz_m]);
                                rho_block[oi] = d2psidx2;
                                l_int_rho += d2psidx2;
                                l_int_V_rho += d2psidx2 * p[i];
                            }

                        }
                    }
                }
                offset_i += range_blocks[ib].Size3D();

            }//end of block loop//

#ifdef CORE_CHARGE_SCALING
            MPI_Allreduce(l_values, sum_values, 2, MPI_DOUBLE, MPI_SUM, ni_comm);
            sum_int_rho *= dx * dy * dz;
            sum_int_V_rho *= dx * dy * dz;
            double scale_ratio = -valence_electron / sum_int_rho;
            //scale_ratio = std::sqrt(scale_ratio);
            for (int i = 0; i < total_block_size; ++i) {
                rho_block[i] *= scale_ratio;
            }
            sum_int_V_rho *= scale_ratio * scale_ratio;
            //sum_int_V_rho *= scale_ratio ;
            /*
            if (IsRoot(ni_comm)) {
                printf("scale_ratio[%d] = %f\n", ni, scale_ratio);
            }
            */
#else
            l_int_rho *= dx * dy * dz;
            l_int_V_rho *= dx * dy * dz;
            MPI_Reduce(l_values, sum_values, 2, MPI_DOUBLE, MPI_SUM, 0, ni_comm);
            
#endif
        }

        if (IsRoot(ni_comm)) {
            return sum_int_V_rho;
        }
        return 0.0;
    }

    /*
    * memory確保しているので再度呼び出す前には一度deleteが必要
    * mpi版とkpoint版でほぼ一緒なので共通化すべし
    */
    void mCreatePccChargeOnBlock(int ni, double** p_pcc_charge, std::vector<GridRange>** p_pcc_range_block,
        const Nucleus nucleus, const GridRange& l_grid) {

        const PseudoPot_MBK* pp = mFindPseudoPot(nucleus.Z);

        if (pp->has_pcc_charge == 0) {
            *p_pcc_charge = nullptr;
            return;
        }

        //PCC-chargeのカットオフを選ぶ//
        double max_cutoff = pp->cutoff_vlocal;
        int adopt_cutoff_l = -1;
        for (int l = 0; l < 4; ++l) {
            if (max_cutoff < pp->cutoff_r[l]) {
                max_cutoff = pp->cutoff_r[l];
                adopt_cutoff_l = l;
            }
        }

        const size_t total_block_size = (adopt_cutoff_l < 0) ? m_vlocal_block[ni].grid_sizes : m_nonlocal_blocks[ni][adopt_cutoff_l].grid_sizes;
        if (total_block_size == 0) {
            *p_pcc_charge = nullptr;
            return;
        }


        {


#ifdef DEBUG_PRINT
            printf("block-size2[l=%d] = %zd\n", l, total_block_size); fflush(stdout);
#endif

            double* rho_block = new double[total_block_size];
            *p_pcc_charge = rho_block;

            const double rr_min = (pp->radius[0]) * (pp->radius[0]);
            const double rr_max = (max_cutoff) * (max_cutoff);
            const auto& range_blocks = (adopt_cutoff_l < 0) ? m_vlocal_block[ni].range_blocks : m_nonlocal_blocks[ni][adopt_cutoff_l].range_blocks;
            *p_pcc_range_block = (adopt_cutoff_l < 0) ? &(m_vlocal_block[ni].range_blocks) : &(m_nonlocal_blocks[ni][adopt_cutoff_l].range_blocks);
            const auto& shifted_grid = (adopt_cutoff_l < 0) ? m_vlocal_block[ni].shifted_grids : m_nonlocal_blocks[ni][adopt_cutoff_l].shifted_grids;
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


                            if (rr_max < rr) {
                                rho_block[i] = 0.0;
                            } else {
                                if (rr_min > rr) {
                                    rho_block[i] = pp->pcc_charge[0];
                                } else {
                                    rho_block[i] = RadialGrid2::GetValueBySquare(rr, pp->pcc_charge, pp->num_radial_grids, pp->xi_min, pp->xi_delta);
                                }


                            }

                        }
                    }
                }
                offset_i += range_blocks[ib].Size3D();

            }

        }

    }



    /*
    * Generate CommForAtoms class (m_comm4atoms)
    * Before this function is called at 2nd time,
    *
    */
protected:
    void CreateComm4Atoms(const Nucleus* nuclei, int num_nuclei, const GridRangeMPI& l_grid) {

        std::vector<CommForAtoms::AtomInfo> atoms;
        for (int ni = 0; ni < num_nuclei; ++ni) {
            const PseudoPot_MBK* pp = mFindPseudoPot(nuclei[ni].Z);
            double max_cutoff = pp->cutoff_vlocal;
            for (int l = 0; l < 4; ++l) {
                if (max_cutoff < pp->cutoff_r[l])max_cutoff = pp->cutoff_r[l];
            }
            atoms.push_back({ nuclei[ni].Rx,nuclei[ni].Ry,nuclei[ni].Rz, max_cutoff });
        }
        //これ以前の処理はPositionを動かしても変化しないので最初に一度初期化するだけにすべし//

        const double box_x = (double)m_grid.size_x * m_dx;
        const double box_y = (double)m_grid.size_y * m_dy;
        const double box_z = (double)m_grid.size_z * m_dz;

        m_comm4atoms.DeleteComms();
        m_comm4atoms.CreateComms(&atoms[0], num_nuclei, l_grid, box_x, box_y, box_z, m_dx, m_dy, m_dz);
        m_num_valid_comm4atoms = m_comm4atoms.CountValidComms();

        m_request4atoms.resize(m_num_valid_comm4atoms);
        m_status4atoms.resize(m_num_valid_comm4atoms);


    }

    void DeleteComm4Atoms() {
        m_comm4atoms.DeleteComms();
    }

    void DeleteBlocks() {
        //内部で保持するmemory解放(メンバのstd::vectorも)の為に,clear必須//
        m_vlocal_block.clear();
        m_hr_vlocal_block.clear();
        for (auto&& nb : m_nonlocal_blocks) {
            nb.clear();
        }
        //m_rho_vlocal.clear();

    }


    /*
    * 擬ポテンシャルを原子核位置に合わせて必要な情報を用意.
    */
public:
    void UpdatePositionLocal(const Nucleus* nuclei, int num_nuclei, const GridRangeMPI& l_grid) {


        //原子核ごとのMPI_Commを用意する(cutoff長の範囲内で担当グリッドとoverlapするプロセスだけのcomm)
        DeleteComm4Atoms();
        CreateComm4Atoms(nuclei, num_nuclei, l_grid);

        //reset memory size//
        mDeleteBuffers(m_rho_vlocal);
        mDeleteBuffers(m_hr_rho_vlocal);
        mDeleteBuffers(m_pcc_charge);
        DeleteBlocks();
        m_rho_vlocal.resize(num_nuclei);
        m_hr_rho_vlocal.resize(num_nuclei);
        m_pcc_charge.resize(num_nuclei);
        m_pcc_range_blocks.resize(num_nuclei);
        m_vlocal_block.resize(num_nuclei);
        m_hr_vlocal_block.resize(num_nuclei);
        m_nonlocal_blocks.resize(num_nuclei);
        m_max_grid_size = 0;
        m_hr_max_grid_size = 0;

        double l_E_nn_self[2] = { 0.0, 0.0 };

        for (int ni = 0; ni < num_nuclei; ++ni) {
            //set m_vlocal_block and m_nonlocal_block//
            //mSetNonlocalProjectorInfo(m_nonlocal_projectors[ni], nuclei[ni]);
            size_t grid_size = mSetSubspaceBlock(m_vlocal_block[ni], m_hr_vlocal_block[ni], m_nonlocal_blocks[ni], nuclei[ni], l_grid);
            if (m_max_grid_size < grid_size) {
                m_max_grid_size = grid_size;
            }

            l_E_nn_self[0] += mCreateChargeVlocalOnBlock_v2(ni, m_vlocal_block[ni], &m_rho_vlocal[ni], m_dx, m_dy, m_dz, nuclei[ni], l_grid);

            const double hr_dx = m_dx / (double)m_HR_ratio_x;
            const double hr_dy = m_dy / (double)m_HR_ratio_y;
            const double hr_dz = m_dz / (double)m_HR_ratio_z;
            l_E_nn_self[1] += mCreateChargeVlocalOnBlock_v2(ni, m_hr_vlocal_block[ni], &m_hr_rho_vlocal[ni], hr_dx, hr_dy, hr_dz, nuclei[ni], l_grid);
            if (m_hr_max_grid_size < m_hr_vlocal_block[ni].grid_sizes) {
                m_hr_max_grid_size = m_hr_vlocal_block[ni].grid_sizes;
            }


            mCreatePccChargeOnBlock(ni, &m_pcc_charge[ni], &m_pcc_range_blocks[ni], nuclei[ni], l_grid);
            //mCreateProjectorYlmOnBlock(m_nonlocal_blocks[ni], m_nonlocal_projectors[ni], nuclei[ni], l_grid);
        }

        m_E_nn_self[0] = 0.0;
        m_E_nn_self[1] = 0.0;
        MPI_Allreduce(&l_E_nn_self, &m_E_nn_self, 2, MPI_DOUBLE, MPI_SUM, l_grid.mpi_comm);

        
        
    }

    double GetSelfCoreEnergy() {
        return m_E_nn_self[0];
    }

    double GetSelfCoreEnergy_HR() {
        return m_E_nn_self[1];
    }


    void SetChargeVlocal_v2(double* l_rho, const GridRangeMPI& l_grid, const Nucleus* nuclei, int num_nuclei) {
        const int proc_id = GetProcessID(l_grid.mpi_comm);
        

        size_t local_size = l_grid.Size3D();
        for (size_t i = 0; i < local_size; ++i) {
            l_rho[i] = 0.0;
        }


        for (int ni = 0; ni < num_nuclei; ++ni) {
            AddSubgridByRanges(l_grid, l_rho, m_vlocal_block[ni].range_blocks, m_rho_vlocal[ni]);
        }
        

    }

    void SetChargeVlocal_HR(double* l_rho, const GridRangeMPI& l_grid, const Nucleus* nuclei, int num_nuclei) {
        const int proc_id = GetProcessID(l_grid.mpi_comm);


        size_t local_size = l_grid.Size3D() * m_HR_ratio_x* m_HR_ratio_y* m_HR_ratio_z;
        for (size_t i = 0; i < local_size; ++i) {
            l_rho[i] = 0.0;
        }

        auto hr_grid = ScaleRange(l_grid, m_HR_ratio_x , m_HR_ratio_y , m_HR_ratio_z);
        for (int ni = 0; ni < num_nuclei; ++ni) {
            AddSubgridByRanges(hr_grid, l_rho, m_hr_vlocal_block[ni].range_blocks, m_hr_rho_vlocal[ni]);
        }


    }


    void SetPccCharge_v2(double* l_rho, const GridRangeMPI& l_grid, const Nucleus* nuclei, int num_nuclei) {
        const int proc_id = GetProcessID(l_grid.mpi_comm);
        

        size_t local_size = l_grid.Size3D();
        for (size_t i = 0; i < local_size; ++i) {
            l_rho[i] = 0.0;
        }


        for (int ni = 0; ni < num_nuclei; ++ni) {
            if (m_pcc_charge[ni] == nullptr) continue;
            AddSubgridByRanges(l_grid, l_rho, *m_pcc_range_blocks[ni], m_pcc_charge[ni]);
        }
        

    }

protected:
    void mInnerChargeVlocal(const double* l_V, const GridRangeMPI& l_grid, std::vector< SubgridBlock>& vlocal_blocks, std::vector< double*>& rho_vlocal, const Nucleus* nuclei, int num_nuclei, double* l_inner, double dV, double* cut_psi) {


        const int num_inner = num_nuclei;
        double* sum_inner = l_inner + num_inner;


        

        //watch_pp.Record(2);

        //slower than Iallreduce with comm 4 atoms
#define TEST_WITH_REDUCE_TO_ROOT


        int num_valid = 0;
        for (int ni = 0; ni < num_nuclei; ++ni) {
            auto mycomm = m_comm4atoms.GetComm(ni);
            if (mycomm == MPI_COMM_NULL) {
#ifdef TEST_WITH_REDUCE_TO_ROOT
                l_inner[ni] = 0.0;
#endif
                continue;
            }
            const int& index = ni;

            const auto& cut_range_block_l = vlocal_blocks[ni].range_blocks;
            const auto grid_size = vlocal_blocks[ni].grid_sizes;

            CutSubgridByRanges(cut_range_block_l, cut_psi, l_grid, l_V);

            //watch_pp.Record(3);
            int index_PYlm = 0;

            const auto* rho = rho_vlocal[ni];

            double inner = 0.0;
            for (int i = 0; i < grid_size; ++i) {
                inner += rho[i] * cut_psi[i];
            }
            l_inner[index] = inner * dV;




            //watch_pp.Record(4);


#ifndef TEST_WITH_REDUCE_TO_ROOT


            MPI_Iallreduce(l_inner + info_range[ni], sum_inner + info_range[ni],
                info_range[ni + 1] - info_range[ni], MPI_DOUBLE, MPI_SUM, mycomm, &m_request4atoms[num_valid]);
            ++num_valid;

#endif

        }
        

#ifndef TEST_WITH_REDUCE_TO_ROOT

        MPI_Waitall(num_valid, &m_request4atoms[0], &m_status4atoms[0]);
        //watch_pp.Record(5);

#else

        MPI_Reduce(l_inner, sum_inner, num_inner, MPI_DOUBLE, MPI_SUM, 0, l_grid.mpi_comm);
        //watch_pp.Record(5);
#endif

    }

public:
    void InnerForChargeVlocal(double* inner, int stride, const double* l_V, const GridRangeMPI& l_grid, const Nucleus* nuclei, int num_nuclei, bool is_high_reso = false) {
        const int proc_id = GetProcessID(l_grid.mpi_comm);


        //watch_pp.Restart();

        const int& num_inner = num_nuclei;

        //watch_pp.Record(0);

        double* l_inner = new double[num_inner * 2];
        double* sum_inner = l_inner + num_inner;
        for (int i = 0; i < num_inner; ++i) {
            l_inner[i] = 0.0;
        }

        if (is_high_reso) {
            GridRangeMPI hr_grid = l_grid;
            hr_grid.begin_x *= m_HR_ratio_x;
            hr_grid.begin_y *= m_HR_ratio_y;
            hr_grid.begin_z *= m_HR_ratio_z;
            hr_grid.end_x *= m_HR_ratio_x;
            hr_grid.end_y *= m_HR_ratio_y;
            hr_grid.end_z *= m_HR_ratio_z;
            double* cut_psi = new double[m_hr_max_grid_size];
            const double dV = m_dx * m_dy * m_dz / (double)(m_HR_ratio_x* m_HR_ratio_y* m_HR_ratio_z);
            mInnerChargeVlocal(l_V, hr_grid, m_hr_vlocal_block, m_hr_rho_vlocal, nuclei, num_nuclei, l_inner, dV, cut_psi);
            delete[] cut_psi;
        }else{
            double* cut_psi = new double[m_max_grid_size];
            const double dV = m_dx * m_dy * m_dz;
            mInnerChargeVlocal(l_V, l_grid, m_vlocal_block, m_rho_vlocal, nuclei, num_nuclei, l_inner, dV, cut_psi);
            delete[] cut_psi;
        }

        for (int ni = 0; ni < num_nuclei; ++ni) {
            inner[ni * stride] = sum_inner[ni];
        }
        delete[]l_inner;
    }

    //Vlocalに相当する原子核の電荷密度を高解像度で解くグリッドの倍率(整数)
    void SetHighResolution(int ratio_x, int ratio_y, int ratio_z) {
        m_HR_ratio_x = ratio_x;
        m_HR_ratio_y = ratio_y;
        m_HR_ratio_z = ratio_z;
    }
};


/*
* memory確保しているので再度呼び出す前には一度deleteが必要
* nonlocal_blockのdestructorにメモリのdeleteの処理が入っているので、
* 実際にはvector<PP_Proj_Ylm_each_atom_kpoint>をclearすればよい。
*/
template<class PP_PROJ_Y>
void tCreateProjectorYlmOnBlock(NonlocalBlocks& nonlocal_blocks,
    PP_PROJ_Y& nonlocal_projectors,
    const Nucleus nucleus, double dx, double dy, double dz, const GridRange& l_grid) {


    int num_projectors = nonlocal_projectors.num_projectors;
    if (num_projectors == 0) {
        return;
    }

    const int MAX_L = nonlocal_projectors.max_l;

    //subspace data of Ylm/////////////////////////////
    nonlocal_projectors.Ylm.resize((MAX_L + 1) * (MAX_L + 1));

    for (int l = 0; l <= MAX_L; ++l) {
        for (int m = -l; m <= l; ++m) {
            auto& Ylm_block = nonlocal_projectors.Ylm[l * l + m + l];
            Ylm_block = new double[nonlocal_blocks[l].grid_sizes];
            SetBlockYlmAround(l, m, Ylm_block, nonlocal_blocks[l].range_blocks, nonlocal_blocks[l].shifted_grids, nucleus.Rx, nucleus.Ry, nucleus.Rz, dx, dy, dz);
            //CutSubgridByRanges(range_block_in_2[l], Ylm_block, subgrid_Ylm, yml_buf);
        }
    }


}

#endif
