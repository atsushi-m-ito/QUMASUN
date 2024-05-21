#if 0
#pragma once
#ifdef USE_MPI
#include <mpi.h>
#include "mpi_helper.h"
#endif
#include "GridRange.h"


//The loop for bandled multi-waves and x,y,z
//Ns: number of states
template <class FUNC>
void ForN(int Ns, const GridRange& grid, FUNC&& func) {
	const size_t size_xyz = grid.Size3D();
	for (size_t i = 0; i < size_xyz; ++i) {
		for (size_t n = 0; n < Ns; ++n) {
			func(n, i);
		}
	}
}

template <class T, class FUNC>
void ForNReduce(int Ns, const GridRange& grid, T* result, FUNC&& func) {

	const size_t size_xyz = grid.Size3D();
	
	for (size_t i = 0; i < size_xyz; ++i) {
		for (size_t n = 0; n < Ns; ++n) {
			result[n] += func(i);
		}
	}	
}

#endif
