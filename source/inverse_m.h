#pragma once
#include <float.h>
#include <cmath>
#if 0
#ifdef _WIN32
#define finite(x) _finite(x)
#endif
#endif

namespace {
	/*
	逆行列の計算
	m: 行列
	size:	行列のサイズ(size x sizeの二次元行列となる)
	inv_m:	逆行列の格納場所
	*/
	inline void InverseM(double* m, int size, double* inv_m) {

		//#pragma omp parallel
		{

			//単位行列を作る
#pragma omp for
			for (int j = 0; j < size; j++) {
				for (int i = 0; i < size; i++) {
					inv_m[i + j*size] = (i == j) ? 1.0 : 0.0;
				}
			}
		}
		//掃き出し法

		for (int i = 0; i < size; i++) {

			double buf = 1 / m[i + i*size];

#if 1
			if (std::isfinite(buf)) {
#else
			if (finite(buf)) {
#endif
				for (int j = 0; j < size; j++) {
					m[i + j*size] *= buf;
					inv_m[i + j*size] *= buf;
				}
			}

			for (int j = 0; j < size; j++) {
				if (i != j) {
					double buf2 = m[j + i*size];


					for (int k = 0; k < size; k++) {
						m[j + k*size] -= m[i + k*size] * buf2;
						inv_m[j + k*size] -= inv_m[i + k*size] * buf2;
					}
				}
			}
		}


	}
}
