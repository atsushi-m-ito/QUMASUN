#pragma once
#include "qumasun_input.h"

//k点の折りたたみルール//
//空空間をK倍した場合と結果が一致する//
//double kpx = (kpoint * 2 > kpoint_limit) ? (double)(kpoint - kpoint_limit) : (double)kpoint;

//これを有効にするとガンマ点を通らなくなり、空間をK倍した場合と結果が一致しなくなる//
//逆空間の平均値としてどちらが適切かは議論の余地あり//
//#define KPOINT_SHIFT_IF_EVEN
//#ifdef KPOINT_SHIFT_IF_EVEN
//if ((kpoint_limit & 0x1) == 0) kpx -= 0.5;//even number
//#endif

/*
inline
double actual_kpoint(int kpoint, int kpoint_limit) {
#if 1
	//k点の折りたたみ//
	//空空間をK倍した場合と結果が一致する//
	double kpx = (kpoint * 2 > kpoint_limit) ? (double)(kpoint - kpoint_limit) : (double)kpoint;
#else
	//k点を折りたたまない//
	//空間をK倍した場合と結果が一致しない
	double kpx = (double)kpoint_x;
#endif


	//これを有効にするとガンマ点を通らなくなり、空間をK倍した場合と結果が一致しなくなる//
	//逆空間の平均値としてどちらが適切かは議論の余地あり//
#ifdef KPOINT_SHIFT_IF_EVEN
	if ((kpoint_limit & 0x1)==0) kpx -= 0.5;//even number
#endif

	return kpx;
}
*/

struct Kpoint3D {
	int kx;
	int ky;
	int kz;
	int weight;
};

inline
int ListupKpoints(int k_num_x, int k_num_y, int k_num_z, uint32_t symmetry, std::vector<Kpoint3D>& kpoint_list) {
	using namespace QUMASUN;
	/*
	for (int kz = 0; kz < k_num_z; ++kz) {
		int kp_z = kz;
		if (kz * 2 > k_num_z) {
			if (symmetry & KPOINT_SYMMETRY::NEGAPOSI_Z){
				continue;
			} else {
				kp_z = kp_z - k_num_z;
			}
		}
		for (int ky = 0; ky < k_num_y; ++ky) {
			if (symmetry & KPOINT_SYMMETRY::MIRROR_YZ) {
				if (ky < kz) {
					continue;
				}
			}
			int kp_y = ky;
			if (ky * 2 > k_num_y) {
				if (symmetry & KPOINT_SYMMETRY::NEGAPOSI_Y) {
					continue;
				} else {
					kp_y = kp_y - k_num_y;
				}
			}
			for (int kx = 0; kx < k_num_x; ++kx) {
				if (symmetry & KPOINT_SYMMETRY::MIRROR_XY) {
					if (kx < ky) {
						continue;
					}
				}
				if (symmetry & KPOINT_SYMMETRY::MIRROR_ZX) {
					if (kx < kz) {
						continue;
					}
				}
				int kp_x = kx;
				if (kx * 2 > k_num_x) {
					if (symmetry & KPOINT_SYMMETRY::NEGAPOSI_X) {
						continue;
					} else {
						kp_x = kp_x - k_num_x;
					}
				}

				kpoint_list.push_back(Kpoint3D{ kp_x, kp_y, kp_z });
			}
		}
	}
	*/


	for (int kz = 0; kz < k_num_z; ++kz) {
		for (int ky = 0; ky < k_num_y; ++ky) {
			for (int kx = 0; kx < k_num_x; ++kx) {
				int kp_x = (kx * 2 > k_num_x) ? kx - k_num_x : kx;
				int kp_y = (ky * 2 > k_num_y) ? ky - k_num_y : ky;
				int kp_z = (kz * 2 > k_num_z) ? kz - k_num_z : kz;
				kpoint_list.push_back(Kpoint3D{ kp_x, kp_y, kp_z ,1});
			}
		}
	}

	auto Find = [&kpoint_list](int kx, int ky, int kz) -> Kpoint3D*{
		for (auto& a : kpoint_list) {
			if ((a.kx == kx) && (a.ky == ky) && (a.kz == kz)) {
				return &a;
			}
		}
		return nullptr;
		};

	if (symmetry & KPOINT_SYMMETRY::NEGAPOSI_Z) {
		for (auto& a : kpoint_list) {
			if (a.weight == 0) continue;
			if (a.kz < 0) {
				auto* m = Find(a.kx, a.ky,  - a.kz);
				if (m == nullptr) {
					return 0; //error//
				}
				m->weight += a.weight;
				a.weight = 0;
			}
		}
	}
	if (symmetry & KPOINT_SYMMETRY::NEGAPOSI_Y) {
		for (auto& a : kpoint_list) {
			if (a.weight == 0) continue;
			if (a.ky <0) {
				auto* m = Find(a.kx,  - a.ky, a.kz);
				if (m == nullptr) {
					return 0; //error//
				}
				m->weight += a.weight;
				a.weight = 0;
			}
		}
	}
	if (symmetry & KPOINT_SYMMETRY::NEGAPOSI_X) {
		for (auto& a : kpoint_list) {
			if (a.weight == 0) continue;
			if (a.kx <0) {
				auto* m = Find( - a.kx, a.ky, a.kz);
				if (m == nullptr) {
					return 0; //error//
				}
				m->weight += a.weight;
				a.weight = 0;
			}
		}
	}
	if (symmetry & KPOINT_SYMMETRY::MIRROR_YZ) {
		for (auto& a : kpoint_list) {
			if (a.weight == 0) continue;
			if (a.kz > a.ky) {
				auto* m = Find(a.kx, a.kz, a.ky);
				if (m == nullptr) {
					return 0; //error//
				}
				m->weight += a.weight;
				a.weight = 0;
			}
		}
	}
	if (symmetry & KPOINT_SYMMETRY::MIRROR_ZX) {
		for (auto& a : kpoint_list) {
			if (a.weight == 0) continue;
			if (a.kz > a.kx) {
				auto* m = Find(a.kz, a.ky, a.kx);
				if (m == nullptr) {
					return 0; //error//
				}
				m->weight += a.weight;
				a.weight = 0;
			}
		}
	}
	if (symmetry & KPOINT_SYMMETRY::MIRROR_XY) {
		for (auto& a : kpoint_list) {
			if (a.weight == 0) continue;
			if (a.ky > a.kx) {
				auto* m = Find(a.ky, a.kx,a.kz);
				if (m == nullptr) {
					return 0; //error//
				}
				m->weight += a.weight;
				a.weight = 0;
			}
		}
	}
	auto res = std::remove_if(kpoint_list.begin(), kpoint_list.end(), [](auto& a) {return a.weight == 0; });
	kpoint_list.erase(res, kpoint_list.end());
	return (int)kpoint_list.size();
}



