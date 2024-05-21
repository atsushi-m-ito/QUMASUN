/*******************************

QUantum MAterial Simulation UNraveler (QUMASUN)

QUMASUN is numerical simulation code for Density Functional Theory(DFT) and Time-dependent DFT based on the real space grid.

********************************/

#include <cstdlib>
#include <cstdio>
#include <cstdint>
#include <cmath>
#include "field_interpolation.h"



//for Test
//全体のfieldから部分領域(cut)を切り取るサンプルプログラム
inline
void TestFieldCut() {
	int src_org_x = 100;
	int src_org_y = 100;
	int src_org_z = 0;
	int src_size_x = 100;
	int src_size_y = 100;
	int src_size_z = 100;
	int global_size_x = 200;
	int global_size_y = 200;
	int global_size_z = 200;

	double* field = new double[src_size_x * src_size_y * src_size_z];

	for (int iz = 0; iz < src_size_z; ++iz) {
		for (int iy = 0; iy < src_size_y; ++iy) {
			for (int ix = 0; ix < src_size_x; ++ix) {
				int i = ix + src_size_x * (iy + src_size_y * iz);
				int val = (ix + src_org_x) + global_size_x * ((iy + src_org_y) + global_size_y * (iz + src_org_z));
				field[i] = (double)val;
			}
		}
	}

	const int cut_begin_x = -5;
	const int cut_begin_y = 115;
	const int cut_begin_z = 1;
	const int cut_size_x = 4;
	const int cut_size_y = 4;
	const int cut_size_z = 4;
	double dest[cut_size_x * cut_size_y * cut_size_z] = { 0.0 };

	CutoutField_periodic(dest, { cut_begin_x, cut_begin_y, cut_begin_z,
						cut_begin_x + cut_size_x, cut_begin_y + cut_size_y, cut_begin_z + cut_size_z },
		field,
		{ src_org_x, src_org_y, src_org_z, src_org_x + src_size_x, src_org_y + src_size_y, src_org_z + src_size_z },
		{ global_size_x, global_size_y, global_size_z });


	auto imodpositive = [](int x, int N) {
		return (x >= 0) ? x % N : (N - ((-x - 1) % N + 1));
		};

	for (int kz = 0; kz < cut_size_z; ++kz) {
		for (int ky = 0; ky < cut_size_y; ++ky) {
			for (int kx = 0; kx < cut_size_x; ++kx) {
				const int k = kx + cut_size_x * (ky + cut_size_y * kz);
				const int val = imodpositive(kx + cut_begin_x, global_size_x)
					+ global_size_x * (imodpositive(ky + cut_begin_y, global_size_y)
						+ global_size_y * imodpositive(kz + cut_begin_z, global_size_z));

				if (fabs(dest[k] - (double)val) > 1.0e-7) {
					printf("field != val: (%d,%d,%d) %f, %f\n", kx, ky, kz, dest[k], (double)val);
				}
			}
		}
	}

	delete[] field;
}



//for Test
//全体のfieldから部分領域(cut)を切り取るサンプルプログラム
inline
void TestInterpolation() {
	int src_org_x = 0;
	int src_org_y = 0;
	int src_org_z = 0;
	int src_size_x = 100;
	int src_size_y = 100;
	int src_size_z = 100;
	int global_size_x = 100;
	int global_size_y = 100;
	int global_size_z = 100;

	double* field = new double[src_size_x * src_size_y * src_size_z];

	const double dx = 1.0 / (double)src_size_x;
	const double dy = 1.0 / (double)src_size_y;
	const double dz = 1.0 / (double)src_size_z;

	auto Value = [](double x, double y, double z) {
		return cos(2.0 * M_PI * (double)x) *
			cos(2.0 * M_PI * (double)y) *
			cos(2.0 * M_PI * (double)z);
		};

	for (int iz = 0; iz < src_size_z; ++iz) {
		for (int iy = 0; iy < src_size_y; ++iy) {
			for (int ix = 0; ix < src_size_x; ++ix) {
				const int i = ix + src_size_x * (iy + src_size_y * iz);

				const double val = Value((double)ix / (double)src_size_x,
					(double)iy / (double)src_size_y,
					(double)iz / (double)src_size_z);

				field[i] = val;
			}
		}
	}


	std::mt19937 engine(123456789);
	std::uniform_real_distribution dist(-0.5, 1.5);
	for (int istep = 0; istep < 100; ++istep)
	{
		double x = dist(engine);
		double y = dist(engine);
		double z = dist(engine);
		if (istep == 0) {
			x = -0.0011534;
			y = 0.9944256;
			z = 0.0089780;
		}

		const double res = Interpolation(x / dx, y / dy, z / dz, field, { src_org_x, src_org_y, src_org_z, src_org_x + src_size_x, src_org_y + src_size_y, src_org_z + src_size_z },
			{ global_size_x, global_size_y, global_size_z });

		const double val = Value(x, y, z);
		if (fabs(res - val) > 1.0e-5) {
			printf("ERROR: interpol != result: (%f,%f,%f) %f, %f\n", x, y, z, res, val);
		} else {
			printf("OK: interpol == result: (%f,%f,%f) %f, %f\n", x, y, z, res, val);
		}

		/*
		auto imodpositive = [](int x, int N) {
			return (x >= 0) ? x % N : (N - ((-x - 1) % N + 1));
			};

		const int cut_begin_x = (int)floor(x / dx) - 1;
		const int cut_begin_y = (int)floor(y / dy) - 1;
		const int cut_begin_z = (int)floor(z / dz) - 1;
		const int cut_size_x = 4;
		const int cut_size_y = 4;
		const int cut_size_z = 4;
		double dest[cut_size_x * cut_size_y * cut_size_z] = { 0.0 };

		CutoutField_periodic(dest, { cut_begin_x, cut_begin_y, cut_begin_z,
							cut_begin_x + cut_size_x, cut_begin_y + cut_size_y, cut_begin_z + cut_size_z },
			field,
			{ src_org_x, src_org_y, src_org_z, src_org_x + src_size_x, src_org_y + src_size_y, src_org_z + src_size_z },
			{ global_size_x, global_size_y, global_size_z });


		for (int kz = 0; kz < cut_size_z; ++kz) {
			for (int ky = 0; ky < cut_size_y; ++ky) {
				for (int kx = 0; kx < cut_size_x; ++kx) {
					const int k = kx + cut_size_x * (ky + cut_size_y * kz);


					const double val = Value((double)(kx + cut_begin_x)/ (double)global_size_x,
						(double)(ky + cut_begin_y)/ (double)global_size_y,
						(double)(kz + cut_begin_z)/ (double)global_size_z);

					if (fabs(dest[k] - (double)val) > 1.0e-7) {
						printf("ERROR: field != val: (%d,%d,%d) %f, %f\n", kx, ky, kz, dest[k], (double)val);
					} else {
					//	printf("OK: field == val: (%d,%d,%d) %f, %f\n", kx, ky, kz, dest[k], (double)val);
					}
				}
			}
		}
		*/
	}


	delete[] field;
}


//test to cut sub firld and to interporate from field
int main(int argc, char* argv[]) {
	TestFieldCut();
	TestInterpolation();
}



