#pragma once
#include <map>
#include "vps_loader.h"
#include "wave_function.h"
#include "nucleus.h"
#include "SubspaceField.h"
#include "RealSphericalHarmonics.h"
#include "vecmath.h"
#include "ellipsoidal_integral.h"

//#define COMMON_CUTOFF_L     //slower than org.

inline double pow2(double x) { return x * x; };
inline double fold(double x, double box_w) { return (x *2.0> box_w) ? x - box_w : (x*2.0 < -box_w) ? x + box_w : x; };

/*
* Pseudo Potentialをsubspace(範囲の限定された)のグリッドデータで保持する
* subspaceの中心は原子座標であり、このクラスも原子毎に保持する
* 同じ原子番号であっても位置が異なればこのクラスのデータは異なるため
*/
struct PP_nonlocal_each_atom {
	std::vector<SubspaceField> projector;
	std::vector<SubspaceField> Ylm;
};


template<class YLM>
void SetYlmAroundAny(YLM& Y, double value_org, double* buf, int ix_begin, int iy_begin, int iz_begin,
		int subsize_x, int subsize_y, int subsize_z, 
		double R0_x, double R0_y, double R0_z,
		double dx, double dy, double dz	)
{

	
	constexpr double rr_min = 1.0e-14;

	for (int iz = 0; iz < subsize_z; ++iz) {
		const double z = dz * (double)(iz + iz_begin) - R0_z;
		const double rz2 = z * z;
		for (int iy = 0; iy < subsize_y; ++iy) {
			const double y = dy * (double)(iy + iy_begin) - R0_y;
			const double ry2 = y * y;
			for (int ix = 0; ix < subsize_x; ++ix) {
				const double x = dx * (double)(ix + ix_begin) - R0_x;
				const double rx2 = x * x;
				const size_t i = (size_t)ix + (size_t)subsize_x * ((size_t)iy + ((size_t)subsize_y * (size_t)iz));

				const double rr = rx2 + ry2 + rz2;
				if (rr_min > rr) {
					buf[i] = value_org;
				} else {
					const double r = sqrt(rr);
					buf[i] = Y(x / r, y / r, z / r);
				}
			}
		}
	}
}

template<int L, int M>
void SetYlmAround(double* buf, int ix_begin, int iy_begin, int iz_begin,
	int subsize_x, int subsize_y, int subsize_z,
	double R0_x, double R0_y, double R0_z,
	double dx, double dy, double dz)
{
	const double Y00 = 1.0 / sqrt(4.0 * M_PI);
	Ylm<L, M> Y;
	SetYlmAroundAny(Y, (L == 0) ? Y00 : 0.0, buf, ix_begin, iy_begin, iz_begin,
		subsize_x, subsize_y, subsize_z,
		R0_x, R0_y, R0_z, dx, dy, dz);
}

class PseudoPotIntegrator {
protected:
    GridInfo m_grid;
    std::map<int, PseudoPot_MBK*> m_pp_list;

    std::vector<PP_nonlocal_each_atom> m_pp_by_atom;

    double m_dx;
    double m_dy;
    double m_dz;

public:

    virtual ~PseudoPotIntegrator() {
        //printf("destruct PseudoPotIntegrator\n"); fflush(stdout);
        for (auto&& a : m_pp_list) {
            a.second->Release();
        }
    }

    void InitializeGrid(const GridInfo& grid_, double dx, double dy, double dz) {
        m_grid = grid_;
        m_dx = dx;
        m_dy = dy;
        m_dz = dz;
    }

    void Load(int Z, const char* filepath) {
        m_pp_list.emplace(Z, LoadVPS(filepath));
    }



    PseudoPot_MBK* mFindPseudoPot(int Z)
    {
        auto it = m_pp_list.find(Z);
        if (it == m_pp_list.end()) {
            return nullptr;
        }
        return it->second;
    }

    /*
    * 原子番号Zの原子のVlocal(動径関数)において、距離rの値を取得
    */
    double GetCutoffVlocal(int Z) {
        const PseudoPot_MBK* pp = mFindPseudoPot(Z);
        return pp->cutoff_vlocal;
    }

    /*
    * 原子番号Zの原子のVlocal(動径関数)において、距離rの値を取得
    */
    double GetVlocalRadial(int Z, double r) {
        const PseudoPot_MBK* pp = mFindPseudoPot(Z);

        const double r_min = pp->radius[0];
        //const double r_max = pp->radius[pp->num_radial_grids - 2];
        const double r_max = pp->cutoff_vlocal;
        if (r_min > r) {
            return pp->V_local[0];
        } else if (r_max <= r) {
            return -(pp->valence_electron) / r;
        } else {
            return RadialGrid2::GetValue(r, pp->V_local, pp->num_radial_grids, pp->xi_min, pp->xi_delta);
        }
    }

    double NumValenceElectron(int Z) {
        const PseudoPot_MBK* pp = mFindPseudoPot(Z);
        return pp->valence_electron;
    }


    /*
    * 原子核間がvlocalのカットオフ長の和よりも短い距離になった場合の補正
    */
    double CoreCoreEnergyCorrection(double R, int iZ_a, int iZ_b, double* force) {

        const PseudoPot_MBK* pp_a = mFindPseudoPot(iZ_a);
        const PseudoPot_MBK* pp_b = mFindPseudoPot(iZ_b);
        const double cutoff_a = pp_a->cutoff_vlocal;
        const double cutoff_b = pp_b->cutoff_vlocal;
        const double Q_a = pp_a->valence_electron;
        const double Q_b = pp_b->valence_electron;
        const double Z_a = (double)iZ_a;
        const double Z_b = (double)iZ_b;


        struct ENN_SELF {
            double E_Vlocal_rho = 0.0;
            double Force_Vlocal_rho = 0.0;

            ENN_SELF& operator+=(const ENN_SELF& b) {
                E_Vlocal_rho += b.E_Vlocal_rho;
                Force_Vlocal_rho += b.Force_Vlocal_rho;
                return *this;
            }

            ENN_SELF operator*(const double b) const {
                ENN_SELF r;
                r.E_Vlocal_rho = E_Vlocal_rho * b;
                r.Force_Vlocal_rho = Force_Vlocal_rho * b;
                return r;
            }
        };

        ENN_SELF Enn = ReGZ::IntegrateEllipsoidal<ENN_SELF, 96, 64>(R, cutoff_a, cutoff_b,
            //double Enn = ReGZ::IntegrateEllipsoidal<96, 64>(R, 1.0, 1.0,
            [&R, &pp_a, &pp_b, &cutoff_a, &cutoff_b](double r_a, double r_b) {
                if (r_b < cutoff_b) {
                    double V = RadialGrid2::GetValue(r_a, pp_a->V_local, pp_a->num_radial_grids, pp_a->xi_min, pp_a->xi_delta);
                    double dV_dr = RadialGrid2::GetDifferential1(r_a, pp_a->V_local, pp_a->num_radial_grids, pp_a->xi_min, pp_a->xi_delta);
                    double rho = (1.0 / (2.0 * M_PI)) * RadialGrid2::GetLaplacian(r_b, pp_b->V_local, 0, pp_b->num_radial_grids, pp_b->xi_min, pp_b->xi_delta);

                    double cos_theta_a = (R * R + r_a * r_a - r_b * r_b) / (2.0 * R * r_a);
                    return ENN_SELF{ V * rho,  -dV_dr * rho * cos_theta_a };
                } else {
                    return ENN_SELF{ 0.0, 0.0 };
                }

            });
        *force = Enn.Force_Vlocal_rho;
        return Enn.E_Vlocal_rho;
    }
};


class PseudoPotIntegrator_single : public PseudoPotIntegrator{
private:
    
    std::vector<PP_nonlocal_each_atom> m_pp_by_atom;

public:

	/*
	* 動径関数として読み込まれたPPのVlocalを、
	* 直交グリッドに焼き直して、Vに足しこむ
	* 
	*/
	void SetPotentialPPVlocal( RspaceFunc<double>& V, const Nucleus* nuclei, int num_nuclei) {
		const int size_3d = m_grid.size_3d;
		vecmath::SetZero<double>(V, size_3d);
		
		const int size_x = m_grid.size_x;
		const int size_y = m_grid.size_y;
		const int size_z = m_grid.size_z;
		const double dx = m_dx;
		const double dy = m_dy;
		const double dz = m_dz;
		const double box_x = dx * (double)size_x;
		const double box_y = dy * (double)size_y;
		const double box_z = dz * (double)size_z;

//#define DEBUG_PRINT2
#ifdef DEBUG_PRINT2
		FILE* fp = fopen("pp.local.txt", "w");
		fprintf(fp, "#rr, Vlocal\n");
#endif
		for (int ni = 0; ni < num_nuclei; ++ni) {
			const PseudoPot_MBK* pp = mFindPseudoPot(nuclei[ni].Z);
			if (pp == nullptr) {
				printf("ERROR: VPS is not loaded: Z = %d\n", nuclei[ni].Z);
				return;
			}

			const double rr_min = pow2(pp->radius[0]);
			const double rr_max = pow2(pp->radius[pp->num_radial_grids-2]);
			
			const double R0_x = nuclei[ni].Rx;
			const double R0_y = nuclei[ni].Ry;
			const double R0_z = nuclei[ni].Rz;
			for (int iz = 0; iz < size_z; ++iz) {
				const double rz2 = pow2(fold(dz * (double)iz - R0_z, box_z));
				for (int iy = 0; iy < size_y; ++iy) {
					const double ry2 = pow2(fold(dy * (double)iy - R0_y, box_y));
					for (int ix = 0; ix < size_x; ++ix) {
						const double rx2 = pow2(fold(dx * (double)ix - R0_x, box_x));
						const size_t i = (size_t)ix + (size_t)size_x * ((size_t)iy + ((size_t)size_y * (size_t)iz));

						const double rr = rx2 + ry2 + rz2;
						if (rr_min > rr) {
							V[i] += pp->V_local[0];
						}else if(rr_max <= rr){
							V[i] += -(pp->valence_electron) / sqrt(rr);
						} else {
							V[i] += RadialGrid2::GetValueBySquare(rr, pp->V_local, pp->num_radial_grids, pp->xi_min, pp->xi_delta);
						}

#ifdef DEBUG_PRINT2
						fprintf(fp, "%f\t%f\n", sqrt(rr), V[i]);
#endif

					}
				}
			}
		}

#ifdef DEBUG_PRINT2
		fclose(fp);
#endif

	}
	

	/*
	* Pcc Charge: Partial Core Correction Charge
	*/
	void SetPccCharge(RspaceFunc<double>& rho, const Nucleus* nuclei, int num_nuclei) {
		//const int size_3d = m_grid.size_3d;

		const int size_x = m_grid.size_x;
		const int size_y = m_grid.size_y;
		const int size_z = m_grid.size_z;
		const double dx = m_dx;
		const double dy = m_dy;
		const double dz = m_dz;
		const double box_x = dx * (double)size_x;
		const double box_y = dy * (double)size_y;
		const double box_z = dz * (double)size_z;

		for (int ni = 0; ni < num_nuclei; ++ni) {
			const PseudoPot_MBK* pp = mFindPseudoPot(nuclei[ni].Z);
			if (pp == nullptr) {
				printf("ERROR: VPS is not loaded: Z = %d\n", nuclei[ni].Z);
				return;
			}

			if (pp->has_pcc_charge) {

                //max cutoff//
                double max_cutoff = pp->cutoff_vlocal;
                
                for (int l = 0; l < 4; ++l) {
                    if (max_cutoff < pp->cutoff_r[l])max_cutoff = pp->cutoff_r[l];
                }
                

				const double rr_min = pp->radius[0] * pp->radius[0];
				//const double rr_max = pow2(pp->radius[pp->num_radial_grids - 2]);
                const double rr_max = pow2(max_cutoff);

				const double R0_x = nuclei[ni].Rx;
				const double R0_y = nuclei[ni].Ry;
				const double R0_z = nuclei[ni].Rz;
				for (int iz = 0; iz < size_z; ++iz) {
					const double rz2 = pow2(fold(dz * (double)iz - R0_z, box_z));
					for (int iy = 0; iy < size_y; ++iy) {
						const double ry2 = pow2(fold(dy * (double)iy - R0_y, box_y));
						for (int ix = 0; ix < size_x; ++ix) {
							const double rx2 = pow2(fold(dx * (double)ix - R0_x, box_x));
							const size_t i = (size_t)ix + (size_t)size_x * ((size_t)iy + ((size_t)size_y * (size_t)iz));

							const double rr = rx2 + ry2 + rz2;
							if (rr_min > rr) {
								rho[i] += pp->pcc_charge[0];
							}else if(rr_max<= rr){
								//rho[i] += 0.0;
							} else {
								rho[i] += RadialGrid2::GetValueBySquare(rr, pp->pcc_charge, pp->num_radial_grids, pp->xi_min, pp->xi_delta);
							}
						}
					}
				}
			}
		}

	}

	double TotalCoreCharge(const Nucleus* nuclei, int num_nuclei) {
		double total = 0.0;
		for (int ni = 0; ni < num_nuclei; ++ni) {
			const PseudoPot_MBK* pp = mFindPseudoPot(nuclei[ni].Z);
			if (pp == nullptr) {
				printf("ERROR: VPS is not loaded: Z = %d\n", nuclei[ni].Z);
				return 0.0;
			}
			total += pp->valence_electron;
		}
		return total;
	}

	/*
	* PPとして読み込まれたVlocalを数値的に二階微分を計算することで
	* VlocalとPoisson方程式で結ばれるChargeを計算する。
	* さらに、そのChargeを動径座標から直交グリッドに焼き直し、rhoに格納する
	* Charge corresponding to Vlocal
	*/
	void SetChargeVlocal(RspaceFunc<double>& rho, RspaceFunc<double>& work, const Nucleus* nuclei, int num_nuclei) {
		//const int size_3d = m_grid.size_3d;

		const int size_x = m_grid.size_x;
		const int size_y = m_grid.size_y;
		const int size_z = m_grid.size_z;
		const size_t size_xyz = size_x * size_y * size_z;
		const double dx = m_dx;
		const double dy = m_dy;
		const double dz = m_dz;
		const double box_x = dx * (double)size_x;
		const double box_y = dy * (double)size_y;
		const double box_z = dz * (double)size_z;

		double sum = 0.0;

		for (int ni = 0; ni < num_nuclei; ++ni) {
			const PseudoPot_MBK* pp = mFindPseudoPot(nuclei[ni].Z);			
			if (pp == nullptr) {
				printf("ERROR: VPS is not loaded: Z = %d\n", nuclei[ni].Z);
				return;
			}
			const double Qi = pp->valence_electron;

			vecmath::SetZero(work.Pointer(), size_xyz);
			double sum = 0.0;

			
			
			const double rr_min = pp->radius[0] * pp->radius[0];
			//const double rr_max = pow2(pp->radius[pp->num_radial_grids - 2]);
			const double rr_max = pow2(pp->cutoff_vlocal);

#if 0
			double max_inside = 0.0;
			double max_outside = 0.0;
#endif

			const double R0_x = nuclei[ni].Rx;
			const double R0_y = nuclei[ni].Ry;
			const double R0_z = nuclei[ni].Rz;
			for (int iz = 0; iz < size_z; ++iz) {
				const double rz2 = pow2(fold(dz * (double)iz - R0_z, box_z));
				for (int iy = 0; iy < size_y; ++iy) {
					const double ry2 = pow2(fold(dy * (double)iy - R0_y, box_y));
					for (int ix = 0; ix < size_x; ++ix) {
						const double rx2 = pow2(fold(dx * (double)ix - R0_x, box_x));
						const size_t i = (size_t)ix + (size_t)size_x * ((size_t)iy + ((size_t)size_y * (size_t)iz));

						const double rr = rx2 + ry2 + rz2;
						double r1;
						if (rr_max < rr) {
							//rho[i] = 0.0;
						}else{

							if (rr_min > rr) {
								r1 = sqrt(rr_min + 1.0e-14);
							} else {
								r1 = sqrt(rr);
							}

							work[i] = (1.0 / (2.0 * M_PI)) * RadialGrid2::GetLaplacian(r1, pp->V_local, 0, pp->num_radial_grids, pp->xi_min, pp->xi_delta);
							sum += work[i];
#if 0
							if (rr_max_vc <= rr) {
								if (max_outside < work[i]) max_outside = work[i];
							} else {
								if (max_inside < work[i]) max_inside = work[i];
							}
#endif
						}
					}
				}
			}
			
			sum *= - dx * dy * dz;

#ifdef DEBUG_PRINT
			printf("Charge corresponding to Vlocal: %f / %f\n", sum, Qi);

#endif
			double factor = Qi / sum;
			vecmath::AddV(rho.Pointer(), work.Pointer(), factor, size_xyz);
		}

	}


	/*
	* 原子核の位置を登録
	* PPの動径形式のデータから、subspaceのxyzグリッドデータに変換
	*/
	void UpdatePosition(const Nucleus* nuclei, int num_nuclei) {
		
		m_pp_by_atom.resize(num_nuclei);
		printf("Size of PP_by_atom = %zd\n", m_pp_by_atom.size()); fflush(stdout);

		for (int ni = 0; ni < num_nuclei; ++ni) {			
			mMakeSubspacePPnonlocal(m_pp_by_atom[ni], nuclei[ni]);
		}


		for (int ni = 0; ni < num_nuclei; ++ni) {
			mCheckOrthogonal(m_pp_by_atom[ni], nuclei[ni]);
		}
	}


	/*
	* 原子核の位置を登録
	* PPの動径形式のデータから、subspaceのxyzグリッドデータに変換
	*/
	void mMakeSubspacePPnonlocal_Radial(std::vector<SubspaceField>& projectors, const Nucleus nucleus) {
		const double dx = m_dx;
		const double dy = m_dy;
		const double dz = m_dz;


		const PseudoPot_MBK* pp = mFindPseudoPot(nucleus.Z);
		if (pp == nullptr) {
			printf("ERROR: VPS is not loaded: Z = %d\n", nucleus.Z);
			return;
		}
		if (pp->num_projectors == 0) {
			return;
		}

#ifdef COMMON_CUTOFF_L
		auto MaxCutoff = [&]() {
			double max_cutoff = 0.0;
			for (int l = 0; l < 4; ++l) {
				if (max_cutoff < pp->cutoff_r[l]) max_cutoff = pp->cutoff_r[l];
			}
			return max_cutoff;
			};
		const double max_cutoff = MaxCutoff();
#endif

		//subspace data of projector w(r)/////////////////////////////
		for (int k = 0; k < pp->num_projectors; ++k) {
			const int l = pp->projector_quantum_l[k];

#ifdef COMMON_CUTOFF_L
			const double cutoff_r = max_cutoff;
#else
			const double cutoff_r = pp->cutoff_r[l];
#endif
			const double rr_max = cutoff_r * cutoff_r;

			GridRange subgrid = GridForSphere(nucleus.Rx, nucleus.Ry, nucleus.Rz, cutoff_r, m_dx, m_dy, m_dz);

			const int ix_begin = subgrid.begin_x;
			const int ix_end = subgrid.end_x;
			const int iy_begin = subgrid.begin_y;
			const int iy_end = subgrid.end_y;
			const int iz_begin = subgrid.begin_z;
			const int iz_end = subgrid.end_z;


			projectors.emplace_back(ix_begin, iy_begin, iz_begin,
				ix_end - ix_begin, iy_end - iy_begin, iz_end - iz_begin);
			SubspaceField& pp_rsubspace = projectors.back();
			double* buf = pp_rsubspace.Reserve();

			//pp_rsubspace.push_back(pp_rsubspace);

			const int subsize_x = ix_end - ix_begin;
			const int subsize_y = iy_end - iy_begin;
			const int subsize_z = iz_end - iz_begin;



			const double rr_min = pp->radius[0] * pp->radius[0];

			const double R0_x = nucleus.Rx;
			const double R0_y = nucleus.Ry;
			const double R0_z = nucleus.Rz;
			for (int iz = 0; iz < subsize_z; ++iz) {
				const double rz2 = pow2(dz * (double)(iz + iz_begin) - R0_z);
				for (int iy = 0; iy < subsize_y; ++iy) {
					const double ry2 = pow2(dy * (double)(iy + iy_begin) - R0_y);
					for (int ix = 0; ix < subsize_x; ++ix) {
						const double rx2 = pow2(dx * (double)(ix + ix_begin) - R0_x);
						const size_t i = (size_t)ix + (size_t)subsize_x * ((size_t)iy + ((size_t)subsize_y * (size_t)iz));

						const double rr = rx2 + ry2 + rz2;
						if (rr_min > rr) {
							const double proj_spin_orbit_up = pp->projector[k * 2][0];
							const double proj_spin_orbit_dn = pp->projector[k * 2 + 1][0];
							buf[i] = (proj_spin_orbit_up + proj_spin_orbit_dn) / 2.0;

						} else if(rr_max < rr){
							buf[i] = 0.0;
						} else {
							//NOTE: if it is not relative DFT with spin-orbit interaction, projectors (j+1/2) and (j-1/2) should be averaged.
							const double proj_spin_orbit_up = RadialGrid2::GetValueBySquare(rr, pp->projector[k * 2], pp->num_radial_grids, pp->xi_min, pp->xi_delta);
							const double proj_spin_orbit_dn = RadialGrid2::GetValueBySquare(rr, pp->projector[k * 2 + 1], pp->num_radial_grids, pp->xi_min, pp->xi_delta);
							buf[i] = (proj_spin_orbit_up + proj_spin_orbit_dn) / 2.0;
						}
					}
				}
			}

		}
	}

	void mMakeSubspacePPnonlocal_Ylm(std::vector<SubspaceField>& Ylm, const Nucleus nucleus) {
		const double dx = m_dx;
		const double dy = m_dy;
		const double dz = m_dz;

		//const double Y00 = 1.0 / sqrt(4.0 * M_PI);

		const PseudoPot_MBK* pp = mFindPseudoPot(nucleus.Z);
		if (pp == nullptr) {
			printf("ERROR: VPS is not loaded: Z = %d\n", nucleus.Z);
			return;
		}
		if (pp->num_projectors == 0) {
			return;
		}

#ifdef COMMON_CUTOFF_L
		auto MaxCutoff = [&]() {
			double max_cutoff = 0.0;
			for (int l = 0; l < 4; ++l) {
				if (max_cutoff < pp->cutoff_r[l]) max_cutoff = pp->cutoff_r[l];
			}
			return max_cutoff;
			};
		const double max_cutoff = MaxCutoff();
#endif

		//subspace data of Ylm/////////////////////////////
		const int MAX_L = pp->MaxProjectorL();
		//constexpr const int MAX_L = 3;
		//double* buf[2 * MAX_L + 1 + 1];
		for (int l = 0; l <= MAX_L; ++l) {

#ifdef COMMON_CUTOFF_L
			const double cutoff_r = max_cutoff;
#else
			const double cutoff_r = pp->cutoff_r[l];
			
#endif
			


			GridRange subgrid = GridForSphere(nucleus.Rx, nucleus.Ry, nucleus.Rz, cutoff_r, m_dx, m_dy, m_dz);

			const int ix_begin = subgrid.begin_x;
			const int ix_end = subgrid.end_x;
			const int iy_begin = subgrid.begin_y;
			const int iy_end = subgrid.end_y;
			const int iz_begin = subgrid.begin_z;
			const int iz_end = subgrid.end_z;


			const int subsize_x = ix_end - ix_begin;
			const int subsize_y = iy_end - iy_begin;
			const int subsize_z = iz_end - iz_begin;

			const double rr_min = 1.0e-14;

			const double R0_x = nucleus.Rx;
			const double R0_y = nucleus.Ry;
			const double R0_z = nucleus.Rz;


			for (int m = -l; m <= l; ++m) {
				Ylm.emplace_back(ix_begin, iy_begin, iz_begin,
					ix_end - ix_begin, iy_end - iy_begin, iz_end - iz_begin);
			}

			switch(l){
			case 0:
			{
				double* buf = Ylm[0].Reserve();
				SetYlmAround<0, 0>(buf, ix_begin, iy_begin, iz_begin, subsize_x, subsize_y, subsize_z,
					R0_x, R0_y, R0_z, dx, dy, dz);
				break;
			}
			case 1:
			{
				const int LL = 1 * 1;
				double* buf = Ylm[LL].Reserve();
				SetYlmAround<1, -1>(buf, ix_begin, iy_begin, iz_begin, subsize_x, subsize_y, subsize_z,
					R0_x, R0_y, R0_z, dx, dy, dz);

				buf = Ylm[LL + 1].Reserve();
				SetYlmAround<1, 0>(buf, ix_begin, iy_begin, iz_begin, subsize_x, subsize_y, subsize_z,
					R0_x, R0_y, R0_z, dx, dy, dz);

				buf = Ylm[LL + 2].Reserve();
				SetYlmAround<1, 1>(buf, ix_begin, iy_begin, iz_begin, subsize_x, subsize_y, subsize_z,
					R0_x, R0_y, R0_z, dx, dy, dz);
				break;
			}
			case 2:
			{
				const int LL = 2 * 2;
				double* buf = Ylm[LL].Reserve();
				SetYlmAround<2, -2>(buf, ix_begin, iy_begin, iz_begin, subsize_x, subsize_y, subsize_z,
					R0_x, R0_y, R0_z, dx, dy, dz);

				buf = Ylm[LL + 1].Reserve();
				SetYlmAround<2, -1>(buf, ix_begin, iy_begin, iz_begin, subsize_x, subsize_y, subsize_z,
					R0_x, R0_y, R0_z, dx, dy, dz);

				buf = Ylm[LL + 2].Reserve();
				SetYlmAround<2, 0>(buf, ix_begin, iy_begin, iz_begin, subsize_x, subsize_y, subsize_z,
					R0_x, R0_y, R0_z, dx, dy, dz);
				
				buf = Ylm[LL + 3].Reserve();
				SetYlmAround<2, 1>(buf, ix_begin, iy_begin, iz_begin, subsize_x, subsize_y, subsize_z,
					R0_x, R0_y, R0_z, dx, dy, dz);
				
				buf = Ylm[LL + 4].Reserve();
				SetYlmAround<2, 2>(buf, ix_begin, iy_begin, iz_begin, subsize_x, subsize_y, subsize_z,
					R0_x, R0_y, R0_z, dx, dy, dz);
				break;
			}
			case 3:
			{
				const int LL = 3 * 3;
				double* buf = Ylm[LL].Reserve();
				SetYlmAround<3, -3>(buf, ix_begin, iy_begin, iz_begin, subsize_x, subsize_y, subsize_z,
					R0_x, R0_y, R0_z, dx, dy, dz);
				
				buf = Ylm[LL+1].Reserve();
				SetYlmAround<3, -2>(buf, ix_begin, iy_begin, iz_begin, subsize_x, subsize_y, subsize_z,
					R0_x, R0_y, R0_z, dx, dy, dz);

				buf = Ylm[LL + 2].Reserve();
				SetYlmAround<3, -1>(buf, ix_begin, iy_begin, iz_begin, subsize_x, subsize_y, subsize_z,
					R0_x, R0_y, R0_z, dx, dy, dz);

				buf = Ylm[LL + 3].Reserve();
				SetYlmAround<3, 0>(buf, ix_begin, iy_begin, iz_begin, subsize_x, subsize_y, subsize_z,
					R0_x, R0_y, R0_z, dx, dy, dz);

				buf = Ylm[LL + 4].Reserve();
				SetYlmAround<3, 1>(buf, ix_begin, iy_begin, iz_begin, subsize_x, subsize_y, subsize_z,
					R0_x, R0_y, R0_z, dx, dy, dz);

				buf = Ylm[LL + 5].Reserve();
				SetYlmAround<3, 2>(buf, ix_begin, iy_begin, iz_begin, subsize_x, subsize_y, subsize_z,
					R0_x, R0_y, R0_z, dx, dy, dz);
				
				buf = Ylm[LL + 6].Reserve();
				SetYlmAround<3, 3>(buf, ix_begin, iy_begin, iz_begin, subsize_x, subsize_y, subsize_z,
					R0_x, R0_y, R0_z, dx, dy, dz);
				break;
			}
			default:

				for (int m = 0; m <= 2*l; ++m) {
					double* buf = Ylm[l*l + m].Reserve(); 						
					YlmSimple Y(l, m);
					SetYlmAroundAny(Y, 0.0, buf, ix_begin, iy_begin, iz_begin, subsize_x, subsize_y, subsize_z,
						R0_x, R0_y, R0_z, dx, dy, dz);
				}
				
			}
		}

	}

	void mMakeSubspacePPnonlocal(PP_nonlocal_each_atom& pp_by_atom, const Nucleus nucleus) {
		const double dx = m_dx;
		const double dy = m_dy;
		const double dz = m_dz;

		mMakeSubspacePPnonlocal_Radial(pp_by_atom.projector, nucleus);
		mMakeSubspacePPnonlocal_Ylm(pp_by_atom.Ylm, nucleus);
	}



	void ProjectionPP(RspaceFunc<double>& Hp, const RspaceFunc<double>& psi, const Nucleus* nuclei, int num_nuclei) {

		const double coef = (m_dx * m_dy * m_dz);
		const double Y_00 = 1.0 / sqrt(4.0 * M_PI);

		for (int ni = 0; ni < num_nuclei; ++ni) {
			const PseudoPot_MBK* pp = mFindPseudoPot(nuclei[ni].Z);
			if (pp == nullptr) {
				printf("ERROR: VPS is not loaded: Z = %d\n", nuclei[ni].Z);
				return;
			}
			if (pp->num_projectors == 0) continue;

			auto pp_by_atom = m_pp_by_atom[ni];

			for (int k = 0; k < pp->num_projectors; ++k) {			
				const int l = pp->projector_quantum_l[k];
				
				if (l == 0) {				
					double inner = coef * Y_00 * InnerProd_1S_1R(pp_by_atom.projector[k], psi, m_grid.size_x, m_grid.size_y, m_grid.size_z);
					inner *= pp->projector_energy_up[k];
					inner *= Y_00; //Add_1Sの演算用のY_00//
					Add_1S(Hp, pp_by_atom.projector[k], inner, m_grid.size_x, m_grid.size_y, m_grid.size_z);
				} else {
					for (int m1 = 0; m1 < 2 * l + 1; ++m1) {
						double inner = coef * InnerProd_2S_1R(pp_by_atom.projector[k], pp_by_atom.Ylm[l * l + m1], psi, m_grid.size_x, m_grid.size_y, m_grid.size_z);
						inner *= pp->projector_energy_up[k];
						Add_2S(Hp, pp_by_atom.projector[k], pp_by_atom.Ylm[l * l + m1], inner, m_grid.size_x, m_grid.size_y, m_grid.size_z);
					}
				}
			}
		}
	}


	double EnergyPPnonlocal(const RspaceFunc<double>& psi, const Nucleus* nuclei, int num_nuclei) {

		const double coef = (m_dx * m_dy * m_dz);
		const double Y_00 = 1.0 / sqrt(4.0 * M_PI);
		double ene = 0.0;

		for (int ni = 0; ni < num_nuclei; ++ni) {
			const PseudoPot_MBK* pp = mFindPseudoPot(nuclei[ni].Z);
			if (pp == nullptr) {
				printf("ERROR: VPS is not loaded: Z = %d\n", nuclei[ni].Z);
				return 0.0;
			}
			if (pp->num_projectors == 0) continue;

			auto pp_by_atom = m_pp_by_atom[ni];

			for (int k = 0; k < pp->num_projectors; ++k) {
				const int l = pp->projector_quantum_l[k];

				if (l == 0) {

					double inner = coef * Y_00 * InnerProd_1S_1R(pp_by_atom.projector[k], psi, m_grid.size_x, m_grid.size_y, m_grid.size_z);
					double ee= inner* inner *pp->projector_energy_up[k];
					ene += ee;

#ifdef DEBUG_PRINT1
					printf("E_nonlocal(%d,0) = %f, %f\n", l, ee, inner);
#endif
				} else {
					for (int m1 = 0; m1 < 2 * l + 1; ++m1) {
						double inner = coef * InnerProd_2S_1R(pp_by_atom.projector[k], pp_by_atom.Ylm[l * l + m1], psi, m_grid.size_x, m_grid.size_y, m_grid.size_z);
						double ee = inner * inner * pp->projector_energy_up[k]; 
						ene += ee;
#ifdef DEBUG_PRINT1
						printf("E_nonlocal(%d,%d) = %f, %f\n", l, m1, ee, inner);
#endif
					}
				}
			}
		}
		return ene;
	}




	/*
	* 原子核の位置を登録
	* PPの動径形式のデータから、subspaceのxyzグリッドデータに変換
	*/
	void mCheckOrthogonal(PP_nonlocal_each_atom& pp_by_atom, const Nucleus nucleus) {
		
		const PseudoPot_MBK* pp = mFindPseudoPot(nucleus.Z);
		if (pp == nullptr) {
			printf("ERROR: VPS is not loaded: Z = %d\n", nucleus.Z);
			return;
		}
		if (pp->num_projectors == 0) {
			return;
		}

		const double xi_delta = pp->xi_delta;
		const double xi_min = pp->xi_min;
		const int ixi_end = pp->num_radial_grids;


		//printf("xi_min, delta, ixi_end = %.10f, %.10f, %d\n", xi_min, xi_delta, ixi_end);

		//subspace data of projector w(r)/////////////////////////////
		for (int k = 0; k < pp->num_projectors; ++k) {
			const int l = pp->projector_quantum_l[k];


			//const double curoff_r = pp->cutoff_r[l];

			{
				double norm = 0.0;
				for (int ixi = 0; ixi < ixi_end; ++ixi) {
					const double r = pp->radius[ixi];
					const double pro1 = pp->projector[k][ixi];
					norm += pro1 * pro1 * r * r * r;
				}

				
#ifdef _DEBUG
				printf("norm_%d = %.10f\n", k, norm);
#endif
			}

			
			//有限サイズグリッド上での直交かどうか検証//
			//MBKの場合はprojectorが直交のはず --> 結果：直交ではなかった//
			
			
			for (int j = 0; j < k; ++j) {
				const int bro_l = pp->projector_quantum_l[j];
				if (bro_l != l)continue;

				double norm = 0.0;
				double inner = 0.0;

				for (int ixi = 0; ixi < ixi_end; ++ixi) {
					const double r = pp->radius[ixi];
					const double pro1 = pp->projector[k][ixi];
					const double pro2 = pp->projector[j][ixi];
					norm += pro1 * pro1 * r * r * r;
					inner += pro1 * pro2 * r * r * r;
				}

				inner *= (xi_delta / (4.0 * M_PI)) /sqrt(norm * xi_delta / (4.0 * M_PI));
#ifdef _DEBUG
				printf("<p_%d|p_%d> = %.10f\n", k,j,inner);
#endif
			}


		}

#ifdef _DEBUG
		fflush(stdout);
#endif

	}

};
