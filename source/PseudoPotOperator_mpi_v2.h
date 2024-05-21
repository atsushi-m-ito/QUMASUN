#ifdef USE_MPI
#pragma once
#include <mpi.h>
#include <cmath>
#include "mpi_helper.h"
#include "PseudoPotOperator_vlocal.h"
#include "SetBlockYlm.h"
#include "StopWatch.h"

//#define USE_PYlm  //遅い




class PseudoPotIntegrator_mpi : public PseudoPotIntegrator_vlocal {
private:

    //Yml * exp(ikx), where k is k-sampling point//
    struct PP_Proj_Ylm_each_atom {
        int max_l;   //max of angular momentum l.
        int num_projectors;
        std::vector<double*> projector; //block化されたprojector, nonlocal項の数だけ存在//
        std::vector<double*> Ylm; //block化されたYlm, (max_l+1)^2の数だけ存在//

        int* quantum_l = nullptr;   //projectorごとのangular momentum l, nonlocal項の数だけ存在//
        double* projector_energy = nullptr;   //projectorごとのangular momentum l, nonlocal項の数だけ存在//

#ifdef USE_PYlm
        int num_PYlm;
        std::vector<double*> PYlm;
#endif


        ~PP_Proj_Ylm_each_atom() {
            for (auto&& p : projector) {
                delete[] p;
            }
            for (auto&& y : Ylm) {
                delete[] y;
            }
#ifdef USE_PYlm
            delete[] PYlm[0];
#endif
        }
    };

	std::vector< PP_Proj_Ylm_each_atom > m_nonlocal_projectors;//block化されたprojectorとYlm
    int m_num_all_nonlocal = 0;
    std::vector<int> m_num_list_nonlocal;


#ifdef TIME_PP_MPI
	StopWatch<10, true> watch_pp;
#else
	StopWatch<10, false> watch_pp;
#endif


public:
	~PseudoPotIntegrator_mpi() {
#ifdef TIME_PP_MPI
		if (is_root) {
			mPrintTime();
		}
#endif
	}


    void mSetNonlocalProjectorInfo(PP_Proj_Ylm_each_atom& nonlocal_projectors,
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



	void mCreateProjectorRadialOnBlock(NonlocalBlocks& nonlocal_blocks, 
        PP_Proj_Ylm_each_atom& nonlocal_projectors,
		const Nucleus nucleus, const GridRange& l_grid) {


		if (nonlocal_projectors.num_projectors == 0) {
			return;
		}

		const int MAX_L = nonlocal_projectors.max_l;
		//auto& range_block_in_1 = nonlocal_block.range_blocks;
		//auto& shifted_grids = nonlocal_block.shifted_grids;
		//auto& grid_sizes = nonlocal_block.grid_sizes;
        //auto& blocks = nonlocal_projectors.blocks;

		const int num_projectors = nonlocal_projectors.num_projectors;
        nonlocal_projectors.projector.resize(num_projectors);

		int gridsize_PYlm = 0;
		int num_PYlm = 0;
		

		const PseudoPot_MBK* pp = mFindPseudoPot(nucleus.Z);

		for (int k = 0; k < num_projectors; ++k) {
			const int l = nonlocal_projectors.quantum_l[k];

			const size_t total_block_size = nonlocal_blocks[l].grid_sizes;
			gridsize_PYlm += total_block_size * (2 * l + 1);
			num_PYlm += (2 * l + 1);

#ifdef DEBUG_PRINT
			printf("block-size2[l=%d] = %zd\n", l, total_block_size); fflush(stdout);
#endif
			auto& proj_block = nonlocal_projectors.projector[k];
			proj_block = new double[total_block_size];

			const double rr_min = (pp->radius[0]) * (pp->radius[0]);
			const double rr_max = (nonlocal_blocks[l].cutoff) * (nonlocal_blocks[l].cutoff);
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
								proj_block[i] = (proj_spin_orbit_up + proj_spin_orbit_dn) / 2.0;

							} else if (rr_max < rr) {
								proj_block[i] = 0.0;
							} else {
								//NOTE: if it is not relative DFT with spin-orbit interaction, projectors (j+1/2) and (j-1/2) should be averaged.
								const double proj_spin_orbit_up = RadialGrid2::GetValueBySquare(rr, pp->projector[k * 2], pp->num_radial_grids, pp->xi_min, pp->xi_delta);
								const double proj_spin_orbit_dn = RadialGrid2::GetValueBySquare(rr, pp->projector[k * 2 + 1], pp->num_radial_grids, pp->xi_min, pp->xi_delta);
								proj_block[i] = (proj_spin_orbit_up + proj_spin_orbit_dn) / 2.0;
							}

						}
					}
				}
				offset_i += range_blocks[ib].Size3D();

			}

		}

#ifdef USE_PYlm
		nonlocal_block.num_PYlm = num_PYlm;
		nonlocal_block.PYlm.resize(num_PYlm);
		double* head_PYlm = new double[gridsize_PYlm];
		int index = 0;
		for (int k = 0; k < num_projectors; ++k) {
			const int l = nonlocal_block.quantum_l[k];
			const int nl_gridsize = nonlocal_block.grid_sizes[l];
			const auto* proj = nonlocal_block.projector[k];
			for (int m1 = 0; m1 < 2 * l + 1; ++m1) {
				const auto* Ylm = nonlocal_block.Ylm[l * l + m1];
				nonlocal_block.PYlm[index] = head_PYlm;
				double* target = nonlocal_block.PYlm[index];
				for (int i = 0; i < nl_gridsize; ++i) {
					target[i] = proj[i] * Ylm[i];
				}
				head_PYlm += nl_gridsize;
				++index;
			}
			
		}
#endif

		
	}



    void CountNonlocalProjector(const Nucleus* nuclei, int num_nuclei) {

        int num_all_nonlocal = 0;
        //std::vector<int> m_info_range;
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
			mCreateProjectorRadialOnBlock(m_nonlocal_blocks[ni], m_nonlocal_projectors[ni], nuclei[ni], l_grid);
		}

	}


	void mInnerNonlocalPsi(const double* l_psi, const GridRangeMPI& l_grid, const Nucleus* nuclei, int num_nuclei, double* l_inner, double dV) {


		const double Y_00 = 1.0 / sqrt(4.0 * M_PI);
		
		double** cut_psi = new double* [MAX_L_SYSTEM];
		cut_psi[0] = new double[m_max_grid_size * MAX_L_SYSTEM];
		for (int l = 1; l < MAX_L_SYSTEM; ++l) {
			cut_psi[l] = cut_psi[0] + m_max_grid_size*l;
		}

		watch_pp.Record(2);

		int index = 0;
		for (int ni = 0; ni < num_nuclei; ++ni) {
			const int max_l = m_nonlocal_projectors[ni].max_l;
			

            const auto& blocks = m_nonlocal_blocks[ni];

			for (int l = 0; l <= max_l; ++l) {

				CutSubgridByRanges(blocks[l].range_blocks, cut_psi[l], l_grid, l_psi);
			}

			watch_pp.Record(3);
			int index_PYlm=0;
			for (int k = 0; k < m_nonlocal_projectors[ni].num_projectors; ++k) {
				const int l = m_nonlocal_projectors[ni].quantum_l[k];
				const auto* psi = cut_psi[l];


#ifdef USE_PYlm
				int nl_size = nonlocal_size[l];
				for (int m1 = 0; m1 < 2 * l + 1; ++m1) {
					const auto* PYlm = m_nonlocal_block[ni].PYlm[index_PYlm];
					double inner = 0.0;
					for (int i = 0; i < nl_size; ++i) {
						inner += PYlm[i] * psi[i];
					}
					l_inner[index] = inner * dV ;
					
					index_PYlm++;
					index++;
				}

#else
				const auto* proj = m_nonlocal_projectors[ni].projector[k];
				
				if (l == 0) {
					/*
					double inner = ForReduce2(cut_grid[l], 0.0,
					[&proj, &l_psi](int64_t i) { return proj[i] * l_psi[i]; });
					*/
					size_t nl_size = blocks[l].grid_sizes;
					double inner = 0.0;
					for (int i = 0; i < nl_size; ++i) {
						inner += proj[i] * psi[i];
					}
					l_inner[index] = inner * dV * Y_00;
					++index;

				} else {
					for (int m1 = 0; m1 < 2 * l + 1; ++m1) {
						const auto* Ylm = m_nonlocal_projectors[ni].Ylm[l * l + m1];
						//double inner = ForReduce2(cut_grid[l], 0.0,
						//	[&proj, &Ylm, &l_psi](int64_t i) { return proj[i] * Ylm[i] * l_psi[i]; });
						size_t nl_size = blocks[l].grid_sizes;
						double inner = 0.0;
						for (int i = 0; i < nl_size; ++i) {
							inner += proj[i] * Ylm[i] * psi[i];
						}
						l_inner[index] = inner * dV;
						++index;
					}
				}
#endif
			}
			watch_pp.Record(4);
		}
		delete[] cut_psi[0];
		delete[] cut_psi;
		
	}

	void mInnerNonlocalPsi_v2(const double* l_psi, const GridRangeMPI& l_grid, const Nucleus* nuclei, int num_nuclei, double* l_inner, const std::vector<int>& info_range, double dV) {


		const double Y_00 = 1.0 / sqrt(4.0 * M_PI);

		const int num_all_nonlocal = info_range[num_nuclei];
		double* sum_inner = l_inner + num_all_nonlocal;

		double** cut_psi = new double* [MAX_L_SYSTEM];
		cut_psi[0] = new double[m_max_grid_size * MAX_L_SYSTEM];
		for (int l = 1; l < MAX_L_SYSTEM; ++l) {
			cut_psi[l] = cut_psi[0] + m_max_grid_size * l;
		}

		watch_pp.Record(2);

//slower than Iallreduce with comm 4 atoms
//#define TEST_WITH_ALLREDUCE


		int num_valid = 0;
		for (int ni = 0; ni < num_nuclei; ++ni) {
			auto mycomm = m_comm4atoms.GetComm(ni);
			if (mycomm == MPI_COMM_NULL) {
#ifdef TEST_WITH_ALLREDUCE
				for (int idx = info_range[ni]; idx < info_range[ni + 1]; ++idx) {
					l_inner[idx] = 0.0;
				}
#endif
				continue;
			}
			int index = info_range[ni];
			const int max_l = m_nonlocal_projectors[ni].max_l;

			//const auto& cut_range_block_l = m_nonlocal_block[ni].range_blocks;
			//const auto& nonlocal_size = m_nonlocal_block[ni].grid_sizes;
            const auto& blocks = m_nonlocal_blocks[ni];

			for (int l = 0; l <= max_l; ++l) {

				CutSubgridByRanges(blocks[l].range_blocks, cut_psi[l], l_grid, l_psi);
			}

			watch_pp.Record(3);
			int index_PYlm = 0;
			for (int k = 0; k < m_nonlocal_projectors[ni].num_projectors; ++k) {
				const int l = m_nonlocal_projectors[ni].quantum_l[k];
				const auto* psi = cut_psi[l];


#ifdef USE_PYlm
				int nl_size = nonlocal_size[l];
				for (int m1 = 0; m1 < 2 * l + 1; ++m1) {
					const auto* PYlm = m_nonlocal_projectors[ni].PYlm[index_PYlm];
					double inner = 0.0;
					for (int i = 0; i < nl_size; ++i) {
						inner += PYlm[i] * psi[i];
					}
					l_inner[index] = inner * dV;

					index_PYlm++;
					index++;
				}

#else
				const auto* proj = m_nonlocal_projectors[ni].projector[k];

				if (l == 0) {
					/*
					double inner = ForReduce2(cut_grid[l], 0.0,
					[&proj, &l_psi](int64_t i) { return proj[i] * l_psi[i]; });
					*/
					size_t nl_size = blocks[l].grid_sizes;
					double inner = 0.0;
					for (int i = 0; i < nl_size; ++i) {
						inner += proj[i] * psi[i];
					}
					l_inner[index] = inner * dV * Y_00;
					++index;

				} else {
					for (int m1 = 0; m1 < 2 * l + 1; ++m1) {
						const auto* Ylm = m_nonlocal_projectors[ni].Ylm[l * l + m1];
						//double inner = ForReduce2(cut_grid[l], 0.0,
						//	[&proj, &Ylm, &l_psi](int64_t i) { return proj[i] * Ylm[i] * l_psi[i]; });
						size_t nl_size = blocks[l].grid_sizes;
						double inner = 0.0;
						for (int i = 0; i < nl_size; ++i) {
							inner += proj[i] * Ylm[i] * psi[i];
						}
						l_inner[index] = inner * dV;
						++index;
					}
				}
#endif
			}

			watch_pp.Record(4);

			
#ifndef TEST_WITH_ALLREDUCE
#if 1
				
			MPI_Iallreduce(l_inner + info_range[ni], sum_inner + info_range[ni],
				info_range[ni + 1] - info_range[ni], MPI_DOUBLE, MPI_SUM, mycomm, &m_request4atoms[num_valid]);
			++num_valid;
#else
			MPI_Allreduce(l_inner + info_range[ni], sum_inner + info_range[ni],
				info_range[ni + 1] - info_range[ni], MPI_DOUBLE, MPI_SUM, mycomm);
			watch_pp.Record(5);
#endif
#endif		
			

		}
		delete[] cut_psi[0];
		delete[] cut_psi;

#ifndef TEST_WITH_ALLREDUCE
		MPI_Waitall(num_valid, &m_request4atoms[0], &m_status4atoms[0]);
		watch_pp.Record(5);
#else		

		MPI_Allreduce(l_inner, sum_inner, num_all_nonlocal, MPI_DOUBLE, MPI_SUM, l_grid.mpi_comm);
		watch_pp.Record(5);
#endif

	}

	void ProjectionPP(double* Hp, const double* l_psi, const GridRangeMPI& l_grid, const Nucleus* nuclei, int num_nuclei) {
		const int proc_id = GetProcessID(l_grid.mpi_comm);

		const double dV = m_dx * m_dy * m_dz;
		const double Y_00 = 1.0/sqrt(4.0 * M_PI);

		watch_pp.Restart();
        const int num_all_nonlocal = m_num_all_nonlocal;
        auto& info_range = m_num_list_nonlocal;
        /*
		std::vector<int> info_range;
		info_range.push_back(0);
		int num_all_nonlocal = 0;
		for (int ni = 0; ni < num_nuclei; ++ni) {
			const PseudoPot_MBK* pp = mFindPseudoPot(nuclei[ni].Z);
			if (pp == nullptr) {
				printf("ERROR: VPS is not loaded: Z = %d\n", nuclei[ni].Z);
				return;
			}
			
			num_all_nonlocal += pp->TotalProjectorLM();
			info_range.push_back(num_all_nonlocal);
		}
        */
		watch_pp.Record(0);
		
		double* l_inner = new double[num_all_nonlocal*2];
		double* sum_inner = l_inner + num_all_nonlocal;
		for (int i = 0; i < num_all_nonlocal; ++i) {
			l_inner[i] = 0.0;
		}

		mInnerNonlocalPsi_v2(l_psi, l_grid, nuclei, num_nuclei, l_inner, info_range, dV);


		double** cut_Hp = new double* [MAX_L_SYSTEM];
		cut_Hp[0] = new double[m_max_grid_size * MAX_L_SYSTEM];
		for (int l = 1; l < MAX_L_SYSTEM; ++l) {
			cut_Hp[l] = cut_Hp[0] + m_max_grid_size * l;
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
				const int nonlocal_sz = blocks[l].grid_sizes;
				for (int i = 0; i < nonlocal_sz; ++i) {
					cut_Hp[l][i] = 0.0;
				}
			}

			watch_pp.Record(2);
			for (int k = 0; k < m_nonlocal_projectors[ni].num_projectors; ++k) {
				const int l = m_nonlocal_projectors[ni].quantum_l[k];
				const auto* proj = m_nonlocal_projectors[ni].projector[k];
				auto* cut_Hp_l = cut_Hp[l];

				if (l == 0) {
					double inner = sum_inner[index];
					++index;
					inner *= m_nonlocal_projectors[ni].projector_energy[k];
					inner *= Y_00;

					size_t nl_size = blocks[l].grid_sizes;
					#pragma ivdep
					for (int i = 0; i < nl_size; ++i) {
						cut_Hp_l[i] += proj[i] * inner;
					}


				} else {
					for (int m1 = 0; m1 < 2 * l + 1; ++m1) {
						const auto* Ylm = m_nonlocal_projectors[ni].Ylm[l * l + m1];
						double inner = sum_inner[index];
						++index;
						inner *= m_nonlocal_projectors[ni].projector_energy[k];

						size_t nl_size = blocks[l].grid_sizes;
						#pragma ivdep
						for (int i = 0; i < nl_size; ++i) {
							cut_Hp_l[i] += proj[i] * Ylm[i] * inner;
						}
					}
				}
			}
			watch_pp.Record(6);

			for (int l = 0; l <= max_l; ++l) {
				AddSubgridByRanges(l_grid, Hp, blocks[l].range_blocks, cut_Hp[l]);
			}
			watch_pp.Record(7);
		}

		delete[] cut_Hp[0];
		delete[] cut_Hp;
		delete[]l_inner;
	}


	double EnergyPPnonlocal(const double* l_psi, const GridRangeMPI& l_grid, const Nucleus* nuclei, int num_nuclei) {

		const double dV = (m_dx * m_dy * m_dz);

        const int num_all_nonlocal = m_num_all_nonlocal;
/*
		int num_all_nonlocal = 0;

		std::vector<int> info_range;
		info_range.push_back(0);
		for (int ni = 0; ni < num_nuclei; ++ni) {
			const PseudoPot_MBK* pp = mFindPseudoPot(nuclei[ni].Z);
			if (pp == nullptr) {
				printf("ERROR: VPS is not loaded: Z = %d\n", nuclei[ni].Z);
				return 0.0;
			}
			
			num_all_nonlocal += pp->TotalProjectorLM();
			info_range.push_back(num_all_nonlocal);
		}
        */

		
		double* l_inner = new double[num_all_nonlocal];
		mInnerNonlocalPsi(l_psi, l_grid, nuclei, num_nuclei, l_inner, dV);


		double* sum_inner = IsRoot(l_grid.mpi_comm) ? new double[num_all_nonlocal] : nullptr;		
		MPI_Reduce(l_inner, sum_inner, num_all_nonlocal, MPI_DOUBLE, MPI_SUM, 0, l_grid.mpi_comm);
		delete[]l_inner;

		double ene = 0.0;
		if (IsRoot(l_grid.mpi_comm) ){


			int index = 0;
			for (int ni = 0; ni < num_nuclei; ++ni) {
				
				for (int k = 0; k < m_nonlocal_projectors[ni].num_projectors; ++k) {
					const int l = m_nonlocal_projectors[ni].quantum_l[k];

					if (l == 0) {

						const double inner = sum_inner[index];
						++index;
						const double ee = inner * inner * m_nonlocal_projectors[ni].projector_energy[k];
						ene += ee;
						//#define DEBUG_PRINT1
#ifdef DEBUG_PRINT1
						printf("E_nonlocal(%d,0) = %f, %f\n", l, ee, inner);
#endif
					} else {
						for (int m1 = 0; m1 < 2 * l + 1; ++m1) {
							const double inner = sum_inner[index];
							++index;
							const double ee = inner * inner * m_nonlocal_projectors[ni].projector_energy[k];
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



	void mInnerNonlocalPsi_bundle(int num_bundle, const double* l_psi, const GridRangeMPI& l_grid, const Nucleus* nuclei, int num_nuclei, double* l_inner, std::vector<int>& info_range, double dV) {


		const double Y_00 = 1.0 / sqrt(4.0 * M_PI);

		const int local_size = l_grid.Size3D();

		const int num_all_nonlocal = info_range[num_nuclei];
		double* sum_inner = l_inner + num_all_nonlocal * num_bundle;



		double** cut_psi = new double* [MAX_L_SYSTEM];
		cut_psi[0] = new double[(size_t)m_max_grid_size * (size_t)MAX_L_SYSTEM* (size_t)num_bundle];
		for (int l = 1; l < MAX_L_SYSTEM; ++l) {
			cut_psi[l] = cut_psi[0] + (size_t)num_bundle * (size_t)m_max_grid_size * l;
		}

		watch_pp.Record(2);

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
				CutSubgridByRanges_bundle(num_bundle, blocks[l].range_blocks, cut_psi[l], l_grid, l_psi);				
			}

			watch_pp.Record(3);

			int num_PYlm=0;
			for (int k = 0; k < m_nonlocal_projectors[ni].num_projectors; ++k) {
				const int l = m_nonlocal_projectors[ni].quantum_l[k];
				num_PYlm += (l + 1) * (l + 1);
			}


			//int index_PYlm = 0;
			for (int k = 0; k < m_nonlocal_projectors[ni].num_projectors; ++k) {
				const int l = m_nonlocal_projectors[ni].quantum_l[k];
				size_t nl_size = blocks[l].grid_sizes;
				

#ifdef USE_PYlm
				int nl_size = nonlocal_size[l];
				for (int m1 = 0; m1 < 2 * l + 1; ++m1) {
					const auto* PYlm = m_nonlocal_block[ni].PYlm[index_PYlm];
					double inner = 0.0;
					for (int i = 0; i < nl_size; ++i) {
						inner += PYlm[i] * psi[i];
					}
					l_inner[index] = inner * dV;

					index_PYlm++;
					index++;
				}

#else
				const auto* proj = m_nonlocal_projectors[ni].projector[k];

				if (l == 0) {
					/*
					double inner = ForReduce2(cut_grid[l], 0.0,
					[&proj, &l_psi](int64_t i) { return proj[i] * l_psi[i]; });
					*/
					for (int n = 0; n < num_bundle; ++n) {
						l_inner[num_bundle * index + n] = 0.0;
					}
					const auto* psi = cut_psi[l];
					for (int i = 0; i < nl_size; ++i) {
						for (int n = 0; n < num_bundle; ++n) {
							l_inner[num_bundle * index + n] += proj[i] * psi[i * num_bundle + n];
						}
						
					}
					for (int n = 0; n < num_bundle; ++n) {
						l_inner[num_bundle * index + n] *= dV* Y_00;
					}
					++index;

				} else {
					for (int m1 = 0; m1 < 2 * l + 1; ++m1) {
						const auto* Ylm = m_nonlocal_projectors[ni].Ylm[l * l + m1];

						for (int n = 0; n < num_bundle; ++n) {
							l_inner[num_bundle * index + n] = 0.0;
						}
						const auto* psi = cut_psi[l];
						for (int i = 0; i < nl_size; ++i) {
							for (int n = 0; n < num_bundle; ++n) {
								l_inner[num_bundle * index + n] += proj[i] * Ylm[i] * psi[i * num_bundle + n];
							}
						}
						for (int n = 0; n < num_bundle; ++n) {
							l_inner[num_bundle * index + n] *= dV;
						}
						++index;
					}
				}
#endif
			}
			
			watch_pp.Record(4);

			if (mycomm != MPI_COMM_NULL) {
#ifndef TEST_WITH_ALLREDUCE

				MPI_Iallreduce(l_inner + info_range[ni]*num_bundle, sum_inner + info_range[ni] * num_bundle,
					(info_range[ni + 1] - info_range[ni]) * num_bundle, MPI_DOUBLE, MPI_SUM, mycomm, &m_request4atoms[num_valid]);
				++num_valid;

#endif		
			}

		}
		delete[] cut_psi[0];
		delete[] cut_psi;




#ifndef TEST_WITH_ALLREDUCE
		MPI_Waitall(num_valid, &m_request4atoms[0], &m_status4atoms[0]);
		watch_pp.Record(5);
#else		

		MPI_Allreduce(l_inner, sum_inner, num_all_nonlocal* num_bundle, MPI_DOUBLE, MPI_SUM, l_grid.mpi_comm);
		watch_pp.Record(5);
#endif

	}

	void ProjectionPP_bundle(int num_bundle, double* Hp, const double* l_psi, const GridRangeMPI& l_grid, const Nucleus* nuclei, int num_nuclei) {
		const int proc_id = GetProcessID(l_grid.mpi_comm);

		const double dV = m_dx * m_dy * m_dz;
		const double Y_00 = 1.0 / sqrt(4.0 * M_PI);
		const size_t local_size = l_grid.Size3D();

		watch_pp.Restart();

        const int num_all_nonlocal = m_num_all_nonlocal;
        auto& info_range = m_num_list_nonlocal;
        /*
		std::vector<int> info_range;
		info_range.push_back(0);
		int num_all_nonlocal = 0;
		for (int ni = 0; ni < num_nuclei; ++ni) {
			const PseudoPot_MBK* pp = mFindPseudoPot(nuclei[ni].Z);
			if (pp == nullptr) {
				printf("ERROR: VPS is not loaded: Z = %d\n", nuclei[ni].Z);
				return;
			}

			num_all_nonlocal += pp->TotalProjectorLM();
			info_range.push_back(num_all_nonlocal);
		}
        */
		watch_pp.Record(0);

		double* l_inner = new double[num_all_nonlocal * num_bundle*2];
		double* sum_inner = l_inner + num_all_nonlocal * num_bundle;
		mInnerNonlocalPsi_bundle(num_bundle, l_psi, l_grid, nuclei, num_nuclei, l_inner, info_range, dV);



		double** cut_Hp = new double* [MAX_L_SYSTEM ];
		cut_Hp[0] = new double[m_max_grid_size * MAX_L_SYSTEM * num_bundle];
		for (int l = 1; l < MAX_L_SYSTEM; ++l) {
			cut_Hp[l] = cut_Hp[0] + m_max_grid_size * num_bundle * l;
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
				const int nonlocal_sz = blocks[l].grid_sizes;
				for (int64_t i = 0; i < (int64_t)num_bundle * (int64_t)nonlocal_sz; ++i) {
					cut_Hp[l][i] = 0.0;					
				}
			}

			watch_pp.Record(2);
			int index_PYlm = 0;
			for (int k = 0; k < m_nonlocal_projectors[ni].num_projectors; ++k) {
				const int l = m_nonlocal_projectors[ni].quantum_l[k];
				const auto* proj = m_nonlocal_projectors[ni].projector[k];
				size_t nl_size = blocks[l].grid_sizes;

				if (l == 0) {
					double factor = m_nonlocal_projectors[ni].projector_energy[k] * Y_00;
					auto* cut_Hp_l = cut_Hp[l];

					#pragma ivdep
					for (int i = 0; i < nl_size; ++i) {
						for (int n = 0; n < num_bundle; ++n) {
							const double inner = sum_inner[num_bundle * index + n];
							cut_Hp_l[i * num_bundle + n] += proj[i] * inner * factor;
						}
					}
					++index;

				} else {
					double factor = m_nonlocal_projectors[ni].projector_energy[k];
					auto* cut_Hp_l = cut_Hp[l];
					for (int m1 = 0; m1 < 2 * l + 1; ++m1) {
						const auto* Ylm = m_nonlocal_projectors[ni].Ylm[l * l + m1];
						#pragma ivdep
						for (int i = 0; i < nl_size; ++i) {
							for (int n = 0; n < num_bundle; ++n) {
								const double inner = sum_inner[num_bundle * index + n];							
								cut_Hp_l[i * num_bundle + n] += proj[i] * Ylm[i] * inner * factor;
							}
						}
						++index;
					}
				}
			}
			watch_pp.Record(6);

			for (int l = 0; l <= max_l; ++l) {
				AddSubgridByRanges_bundle(num_bundle, l_grid, Hp, blocks[l].range_blocks, cut_Hp[l]);
			}
			
			watch_pp.Record(7);
		}
		delete[] cut_Hp[0];
		delete[] cut_Hp;





		delete[] l_inner;
	}


	double ForceNonlocal(double* forces, const double* l_psi, const double* l_dpsi_dx, const GridRangeMPI& l_grid, const Nucleus* nuclei, int num_nuclei) {
		const int proc_id = GetProcessID(l_grid.mpi_comm);
		constexpr int stride = 3;
		const int local_size = l_grid.Size3D();
		const double* l_dpsi_dy = l_dpsi_dx + local_size;
		const double* l_dpsi_dz = l_dpsi_dx + local_size * 2;


		//watch_pp.Restart();

		const double dV = (m_dx * m_dy * m_dz);


        const int num_all_nonlocal = m_num_all_nonlocal;
        /*
		int num_all_nonlocal = 0;

		std::vector<int> info_range;
		info_range.push_back(0);
		for (int ni = 0; ni < num_nuclei; ++ni) {
			const PseudoPot_MBK* pp = mFindPseudoPot(nuclei[ni].Z);
			if (pp == nullptr) {
				printf("ERROR: VPS is not loaded: Z = %d\n", nuclei[ni].Z);
				return 0.0;
			}

			num_all_nonlocal += pp->TotalProjectorLM();
			info_range.push_back(num_all_nonlocal);
		}
        */
		//watch_pp.Record(0);


		double* l_inner = new double[num_all_nonlocal*4];
		mInnerNonlocalPsi(l_psi, l_grid, nuclei, num_nuclei, l_inner, dV);
		mInnerNonlocalPsi(l_dpsi_dx, l_grid, nuclei, num_nuclei, l_inner + num_all_nonlocal, dV);
		mInnerNonlocalPsi(l_dpsi_dy, l_grid, nuclei, num_nuclei, l_inner + num_all_nonlocal*2, dV);
		mInnerNonlocalPsi(l_dpsi_dz, l_grid, nuclei, num_nuclei, l_inner + num_all_nonlocal*3, dV);


		double* sum_inner = IsRoot(l_grid.mpi_comm) ? new double[num_all_nonlocal*4] : nullptr;
		MPI_Reduce(l_inner, sum_inner, num_all_nonlocal*4, MPI_DOUBLE, MPI_SUM, 0, l_grid.mpi_comm);
		delete[]l_inner;


		double ene = 0.0;
		if (IsRoot(l_grid.mpi_comm)) {


			const double* sum_diff_x = sum_inner + num_all_nonlocal;
			const double* sum_diff_y = sum_inner + num_all_nonlocal * 2;
			const double* sum_diff_z = sum_inner + num_all_nonlocal * 3;

			int index = 0;
			for (int ni = 0; ni < num_nuclei; ++ni) {

				double force_i[3] = { 0.0, 0.0, 0.0 };

				for (int k = 0; k < m_nonlocal_projectors[ni].num_projectors; ++k) {
					const int l = m_nonlocal_projectors[ni].quantum_l[k];

					if (l == 0) {

						const double inner = sum_inner[index];
						const double diff_x = sum_diff_x[index];
						const double diff_y = sum_diff_y[index];
						const double diff_z = sum_diff_z[index];
						++index;
						const double ee = inner * inner * m_nonlocal_projectors[ni].projector_energy[k];
						ene += ee;
						force_i[0] += diff_x * inner * m_nonlocal_projectors[ni].projector_energy[k];
						force_i[1] += diff_y * inner * m_nonlocal_projectors[ni].projector_energy[k];
						force_i[2] += diff_z * inner * m_nonlocal_projectors[ni].projector_energy[k];
						
						
#ifdef DEBUG_PRINT1
						printf("E_nonlocal(%d,0) = %f, %f\n", l, ee, inner);
#endif
					} else {
						for (int m1 = 0; m1 < 2 * l + 1; ++m1) {
							const double inner = sum_inner[index];
							const double diff_x = sum_diff_x[index];
							const double diff_y = sum_diff_y[index];
							const double diff_z = sum_diff_z[index];
							++index;
							const double ee = inner * inner * m_nonlocal_projectors[ni].projector_energy[k];
							force_i[0] += diff_x * inner * m_nonlocal_projectors[ni].projector_energy[k];
							force_i[1] += diff_y * inner * m_nonlocal_projectors[ni].projector_energy[k];
							force_i[2] += diff_z * inner * m_nonlocal_projectors[ni].projector_energy[k];
							ene += ee;
#ifdef DEBUG_PRINT1
							printf("E_nonlocal(%d,%d) = %f, %f\n", l, m1, ee, inner);
#endif
							//#undef DEBUG_PRINT1
						}
					}
				}

				forces[ni * stride + 0] += force_i[0];
				forces[ni * stride + 1] += force_i[1];
				forces[ni * stride + 2] += force_i[2];
			}
		
			delete[]sum_inner;
		}


		return ene;
	}

	void mPrintTime() {
#ifdef TIME_PP_MPI
        printf("PP calculation times====\n");
		watch_pp.Print("Count_num", 0);
		watch_pp.Print("Find_pp  ", 1);
		watch_pp.Print("Allocate ", 2);
		watch_pp.Print("Cut_psi  ", 3);
		watch_pp.Print("Inner    ", 4);
		watch_pp.Print("MPI      ", 5);
		watch_pp.Print("Add local", 6);
		watch_pp.Print("Paste    ", 7);
#endif
	}



};

#endif
