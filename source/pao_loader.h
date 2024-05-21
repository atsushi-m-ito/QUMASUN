#pragma once
#include <cstring>
#include <string>
#include <vector>
#include "JobParam.h"
#include "RadialGrid2.h"

/*
* 現状ではL=0のorbitalだけを読む
*/
class Orbital_PAO {
public:
	int num_orbitals = 0;
	
	
	int num_radial_grids;
	double* buffer_all = nullptr;
	double* radius = nullptr;
	double* valence_charge = nullptr;
	std::vector<double*> orbitals;

	
	double xi_min;   //required to interpolate values on xi grid//
	double xi_delta; //required to interpolate values on xi grid//

public:
	/*
	* Release all memory, and this is not called in destructor
	*/
	void Release() {
		delete[] radius;
		delete[] buffer_all;
		
	}
};

inline
Orbital_PAO LoadPAO(const char* filepath) {
	JobParam job_param(filepath, "<*", "*>", "*#\r\n");
	Orbital_PAO pao;
	
	
	pao.num_orbitals = job_param.GetInt("PAO.Mul");

	
	{
		const char* BLOCK_KEY = "valence.charge.density";
		const int num_grids = job_param.GetBlockNumLines(BLOCK_KEY);
		pao.num_radial_grids = num_grids;
		//pao.buffer_all = new double[num_grids * (pao.num_orbitals + 2)];
		pao.radius = new double[num_grids * 2];
		pao.valence_charge = pao.radius + num_grids;

		double sum_rho = 0.0;
		double xi_min;
		double xi;
		for (int i = 0; i < num_grids; ++i) {
			std::string line(job_param.GetBlockString(BLOCK_KEY, i));
			const char* pbegin = line.c_str();
			char* pend = nullptr;
			xi = strtod(pbegin, &pend);
			if (i == 0)xi_min = xi;
			pbegin = pend;
			const double r = strtod(pbegin, &pend);
			pao.radius[i] = r;
			pbegin = pend;
			pao.valence_charge[i] = strtod(pbegin, &pend);
			pbegin = pend;

			sum_rho += pao.valence_charge[i] * r * r * r;
		}
		
		pao.xi_min = xi_min;
		pao.xi_delta = (xi - xi_min) / (double)(num_grids - 1);

		sum_rho *= pao.xi_delta;
		sum_rho *= 4.0 * M_PI;
		printf("File read: valence charge = %f\n", sum_rho);

	}


	{
		const char* BLOCK_KEY = "pseudo.atomic.orbitals.L=0";
		const int num_grids = job_param.GetBlockNumLines(BLOCK_KEY);
		pao.num_radial_grids = num_grids;
		pao.buffer_all = new double[num_grids * (pao.num_orbitals)];
		for (int i = 0; i < pao.num_orbitals; ++i) {
			pao.orbitals.push_back(pao.buffer_all + i * num_grids);
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
			const double r = strtod(pbegin, &pend);
			//pao.radius[i] = r;
			pbegin = pend;
			for (int k = 0; k < pao.num_orbitals; ++k) {
				pao.orbitals[k][i] = strtod(pbegin, &pend);
				pbegin = pend;
			}

		}
		//pao.xi_min = xi_min;
		//pao.xi_delta = (xi - xi_min) / (double)(num_grids - 1);

	}


		
	return pao;
			
	
}
