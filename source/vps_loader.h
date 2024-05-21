#pragma once
#include <cstring>
#include <string>
#include <cmath>
#include <vector>
#include "JobParam.h"
#include "RadialGrid2.h"

class PseudoPot_MBK {
public:
	int valence_electron = 0;
	int num_projectors = 0;
	int num_radial_grids = 0;
	double cutoff_vlocal = 0.0;
	//cutoff length in radial coordinate for l=0,1,2,3(angler q number)//
	double cutoff_r[4]{ 0.0 ,0.0 ,0.0 ,0.0 };  

	//std::vector<int> projector_quantum_l;
	int* projector_quantum_l = nullptr;
	double* projector_energy_up = nullptr;
	double* projector_energy_down = nullptr;

	double* buffer_all = nullptr;
	double* radius = nullptr;
	double* V_local = nullptr;
	std::vector<double*> projector;

	double* pcc_charge = nullptr;
    int has_pcc_charge = 0;
	
	double xi_min;   //required to interpolate values on xi grid//
	double xi_delta; //required to interpolate values on xi grid//

public:
	/*
	* Release all memory, and this is not called in destructor
	*/
	void Release() {
		delete[] buffer_all;
		delete[] projector_quantum_l;
		delete[] projector_energy_up; 
	}

	int TotalProjectorLM() const{
		int num = 0;
		for (int k = 0; k < num_projectors; ++k) {
			const int l = projector_quantum_l[k];
			num += 2 * l + 1;
		}
		return num;
	}

	int MaxProjectorL() const {
		int max_l = 0;
		for (int k = 0; k < num_projectors; ++k) {
			const int l = projector_quantum_l[k];
			if (max_l < l)max_l = l;
		}
		return max_l;
	}

	template<class T>
	void mSwap(T& a, T& b) {
		T c = a;
		b = a;
		a = c;
	}

	void SortByL() {
		for (int k = 1; k < num_projectors; ++k) {
			for (int j = 2; j < num_projectors - j; ++j) {
				const int l0 = projector_quantum_l[j - 1];
				const int l = projector_quantum_l[j];
				if (l < l0) {
					mSwap(projector_quantum_l[j - 1], projector_quantum_l[j]);
					mSwap(projector_energy_up[j - 1], projector_energy_up[j]);
					mSwap(projector_energy_down[j - 1], projector_energy_down[j]);
					mSwap(projector[(j - 1) * 2], projector[j * 2]);
					mSwap(projector[(j - 1) * 2 + 1], projector[j * 2 + 1]);
				}
			}
		}
	}
};

inline
PseudoPot_MBK* LoadVPS(const char* filepath) {
	JobParam job_param(filepath, "<*", "*>", "*#\r\n", "<VPS>");
	
	if (!job_param.IsEqualString("vps.type", "MBK")){
		printf("ERROR: VPS of MBK is only supported.");
		return nullptr; //error//
	}

	PseudoPot_MBK* p_pseudo_pot = new PseudoPot_MBK;
	PseudoPot_MBK& pp = *p_pseudo_pot;
	{
		//double max_cutoff = 0.0;
		//int max_l = 0;
		const char* BLOCK_KEY = "pseudo.NandL";
		const int num = job_param.GetBlockNumLines(BLOCK_KEY);
		for (int i = 0; i < num; ++i) {
			std::string line (job_param.GetBlockString(BLOCK_KEY,i));
			const char* pbegin = line.c_str();
			char* pend = nullptr;
			const int idx = strtol(pbegin, &pend, 10);
			pbegin = pend;
			const int quantum_n = strtol(pbegin, &pend, 10);
			pbegin = pend;
			const int quantum_l = strtol(pbegin, &pend, 10);
			pbegin = pend;
			const double cutoff = strtod(pbegin, &pend);
			if (quantum_l < 4) {
				if (pp.cutoff_r[quantum_l] < cutoff) {
					pp.cutoff_r[quantum_l] = cutoff;
					//if (max_cutoff < cutoff) max_cutoff = cutoff;
					//if (max_l < quantum_l) max_l = quantum_l;
				}
			}
		}	


	}
	
	{
		const char* BLOCK_KEY = "project.energies";
		const int i_end = job_param.GetBlockNumLines(BLOCK_KEY);
		const int num = strtol(job_param.GetBlockString(BLOCK_KEY, 0),nullptr,10);
		pp.num_projectors = num;
		pp.projector_quantum_l = new int[num];
		pp.projector_energy_up = new double[num*2];
		pp.projector_energy_down = pp.projector_energy_up + num;
		for (int i = 1; i < i_end; ++i) {
			std::string line(job_param.GetBlockString(BLOCK_KEY, i));
			const char* pbegin = line.c_str();
			char* pend = nullptr;
			const int quantum_l = strtol(pbegin, &pend, 10);
			pbegin = pend;
			const double E_up = strtod(pbegin, &pend);
			pbegin = pend;
			const double E_down = strtod(pbegin, &pend);
			pp.projector_quantum_l[i - 1] = quantum_l;
			pp.projector_energy_up[i - 1] = E_up;
			pp.projector_energy_down[i - 1] = E_down;
		}
	}




	{
		const char* BLOCK_KEY = "Pseudo.Potentials";
		const int num_grids = job_param.GetBlockNumLines(BLOCK_KEY);
		pp.num_radial_grids = num_grids;
		pp.buffer_all = new double[num_grids * (pp.num_projectors * 2  + 2 + 1 )];//+3 means, radius, v_local, pcc_charge//
		pp.radius = pp.buffer_all;
		pp.V_local = pp.buffer_all + num_grids;
		for (int i = 0; i < pp.num_projectors * 2; ++i){
			pp.projector.push_back(pp.buffer_all + (2+i)* num_grids);
		}

		double xi_min;
		double xi;
		for (int i = 0; i < num_grids; ++i) {
			std::string line(job_param.GetBlockString(BLOCK_KEY, i));
			const char* pbegin = line.c_str();
			char* pend = nullptr;
			xi = strtod(pbegin, &pend);
			if (i == 0)xi_min = xi;
			pbegin = pend;
			pp.radius[i] = strtod(pbegin, &pend);
			pbegin = pend;
			pp.V_local[i] = strtod(pbegin, &pend);
			pbegin = pend;
			for (int k = 0; k < pp.num_projectors * 2; ++k) {
				pp.projector[k][i] = strtod(pbegin, &pend);
				pbegin = pend;
			}

		}
		pp.xi_min = xi_min;
		pp.xi_delta = (xi - xi_min) / (double)(num_grids - 1);

	}


	{
        
        const char* BLOCK_KEY = "density.PCC";
        const int num_grids = job_param.GetBlockNumLines(BLOCK_KEY);
		if (num_grids > 0) {
			if (pp.num_radial_grids != num_grids) {
				printf("ERROR: %s has incorrect grid size.", BLOCK_KEY);
			}
            
            pp.pcc_charge = pp.buffer_all + (2 + (pp.num_projectors * 2)) * num_grids;
            pp.has_pcc_charge = 1;
			
			double sum_rho=0.0;
			double prev_pcc = 0.0;
			double prev_r = 0.0;
			for (int i = 0; i < num_grids; ++i) {
				std::string line(job_param.GetBlockString(BLOCK_KEY, i));
				const char* pbegin = line.c_str();
				char* pend = nullptr;
				double xi = strtod(pbegin, &pend);
				//if (i == 0)xi_min = xi;
				pbegin = pend;
				const double r = strtod(pbegin, &pend);
				pbegin = pend;
				const double pcc_c = strtod(pbegin, &pend);				
				pbegin = pend;
				pp.pcc_charge[i] = pcc_c;
				//printf("%d, %.12f, %.12f\n", i, r, pcc_c);

#if 1
				sum_rho += pcc_c*r*r*r;
#elif 1
				sum_rho += prev_pcc * prev_r * prev_r * (r - prev_r);
				prev_pcc = pcc_c;
				prev_r = r;
#else
				if (i < num_grids - 1) {
					sum_rho += pcc_c * pp.radius[i] * pp.radius[i] * (pp.radius[i + 1] - pp.radius[i]);
				}
#endif
			}
			const double xi_delta = pp.xi_delta;
			sum_rho *= xi_delta;
			sum_rho *= 4.0 * M_PI;
			//sum_rho *= 4.0 * std::acos(-1.0);

			printf("PCC charge in Pseudo Potential = %f\n", sum_rho);
		}
	}
	
	
	pp.valence_electron = (int)job_param.GetDouble("valence.electron");
	pp.cutoff_vlocal = job_param.GetDouble("local.cutoff");
	
	pp.SortByL();

	return p_pseudo_pot;
			
	
}

