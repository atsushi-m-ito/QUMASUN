#pragma once
#include <cfloat>
#include <cmath>
#include <cstdio>
#include <cstdlib>
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
	inline void OutputMatrix(double* m, int size1, int size2, const char* filepath) {

		FILE* fp = fopen(filepath, "w");

		for (int j = 0; j < size2; j++) {
			for (int i = 0; i < size1 - 1; i++) {
				fprintf(fp, "%f\t", m[i + size1* j]);
			}
			fprintf(fp, "%f\n", m[(size1 - 1) + size1 * j]);
		}

		fclose(fp);
	}
}
