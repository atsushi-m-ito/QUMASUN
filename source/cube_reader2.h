#pragma once
//#include <ctchar>
#include <cstdio>
#include <cstring>
//#include <vec3.h>
//#include <mat33.h>
#include <charconv>
//#include "ATOMS_DATA.h"
#include "transpose.h"

class CubeReader2 {
public:
	CubeReader2() = default;
	~CubeReader2() = default;

	struct Frame {
		int grid_x;//グリッド数
		int grid_y;
		int grid_z;
		double boxaxis[9];
		double boxorg[3];
        char comment[1024];
	};

#if 1

    static
        double* LoadCube(const char* filename, Frame* frame) {

        FILE* fp = fopen(filename, "r");

        const size_t SIZE = 1024;
        char line[SIZE];

        if (fgets(line, 1024, fp) == NULL) { fclose(fp); return NULL; }	//タイトル
        if (fgets(line, 1024, fp) == NULL) { fclose(fp); return NULL; }	//コメント
        strcpy(frame->comment, line);

        if (fgets(line, 1024, fp) == NULL) { fclose(fp); return NULL; }	//原子数とvolumeデータの原点座標
        int pcnt;
        sscanf(line, "%d %lf %lf %lf", &pcnt, &(frame->boxorg[0]), &(frame->boxorg[1]), &(frame->boxorg[2]));

        if (fgets(line, 1024, fp) == NULL) { fclose(fp); return NULL; }	//x方向のメッシュ数とx軸
        sscanf(line, "%d %lf %lf %lf", &(frame->grid_x), &(frame->boxaxis[0]), &(frame->boxaxis[1]), &(frame->boxaxis[2]));

        if (fgets(line, 1024, fp) == NULL) { fclose(fp); return NULL; }	//y方向のメッシュ数とx軸
        sscanf(line, "%d %lf %lf %lf", &(frame->grid_y), &(frame->boxaxis[3]), &(frame->boxaxis[4]), &(frame->boxaxis[5]));

        if (fgets(line, 1024, fp) == NULL) { fclose(fp); return NULL; }	//z方向のメッシュ数とx軸
        sscanf(line, "%d %lf %lf %lf", &(frame->grid_z), &(frame->boxaxis[6]), &(frame->boxaxis[7]), &(frame->boxaxis[8]));

        //frame->boxorg[0] -= (frame->boxaxis[0] + frame->boxaxis[3] + frame->boxaxis[6]) / 2.0;
        //frame->boxorg[1] -= (frame->boxaxis[1] + frame->boxaxis[4] + frame->boxaxis[7]) / 2.0;
        //frame->boxorg[2] -= (frame->boxaxis[2] + frame->boxaxis[5] + frame->boxaxis[8]) / 2.0;
        frame->boxaxis[0] *= (double)frame->grid_x;
        frame->boxaxis[1] *= (double)frame->grid_x;
        frame->boxaxis[2] *= (double)frame->grid_x;
        frame->boxaxis[3] *= (double)frame->grid_y;
        frame->boxaxis[4] *= (double)frame->grid_y;
        frame->boxaxis[5] *= (double)frame->grid_y;
        frame->boxaxis[6] *= (double)frame->grid_z;
        frame->boxaxis[7] *= (double)frame->grid_z;
        frame->boxaxis[8] *= (double)frame->grid_z;


        {
            for (int i = 0; i < pcnt; i++) {
                if (fgets(line, 1024, fp) == NULL) { fclose(fp); return nullptr; }	//原子位置読み込み				
            }
        }

        //ボリュームデータの読み込み
        //const int64_t width = frame->grid_x;
        size_t mesh_sz = (frame->grid_x) * (frame->grid_y) * (frame->grid_z);
        double* values = new double[mesh_sz];
        //char* tp;
        double maxvalue = 0.0;
        double minvalue = 1e10;
        int64_t ix = 0;
        int64_t iy = 0;
        int64_t iz = 0;
        int64_t ioffset = 0;
        const int64_t width = (frame->grid_z);
        const char* const endptr = line + SIZE - 1;
        while (fgets(line, SIZE, fp) != NULL) {
            char* tp = line;
            const int64_t i_end = std::min(iz + 6, width);
            for (; iz < i_end; ++iz) {
                char* ep = nullptr;
                values[iz + ioffset] = strtod(tp, &ep);
                tp = ep+1;
            }
            if (iz == width) {
                iz = 0; 
                ++iy;
                if (iy == frame->grid_y) {
                    iy = 0;
                    ++ix;
                }
                ioffset = width * (ix + frame->grid_x * iy);
            }            
        }

        fclose(fp);

        double* trans = new double[mesh_sz];
        transpose(width, (frame->grid_x) * (frame->grid_y), values, width, trans, (frame->grid_x) * (frame->grid_y));
        delete[] values;

        return trans;

    }
#else
	static 
	double* LoadCube(const char* filename, Frame* frame) {

		FILE* fp = fopen(filename, "r");

		const size_t SIZE = 1024;
		char line[SIZE];

		if (fgets(line, 1024, fp) == NULL) { fclose(fp); return NULL; }	//タイトル
		if (fgets(line, 1024, fp) == NULL) { fclose(fp); return NULL; }	//コメント
        strcpy(frame->comment, line);

		if (fgets(line, 1024, fp) == NULL) { fclose(fp); return NULL; }	//原子数とvolumeデータの原点座標
		int pcnt;
		sscanf(line, "%d %lf %lf %lf", &pcnt, &(frame->boxorg[0]), &(frame->boxorg[1]), &(frame->boxorg[2]));

		if (fgets(line, 1024, fp) == NULL) { fclose(fp); return NULL; }	//x方向のメッシュ数とx軸
		sscanf(line, "%d %lf %lf %lf", &(frame->grid_x), &(frame->boxaxis[0]), &(frame->boxaxis[1]), &(frame->boxaxis[2]));

		if (fgets(line, 1024, fp) == NULL) { fclose(fp); return NULL; }	//y方向のメッシュ数とx軸
		sscanf(line, "%d %lf %lf %lf", &(frame->grid_y), &(frame->boxaxis[3]), &(frame->boxaxis[4]), &(frame->boxaxis[5]));

		if (fgets(line, 1024, fp) == NULL) { fclose(fp); return NULL; }	//z方向のメッシュ数とx軸
		sscanf(line, "%d %lf %lf %lf", &(frame->grid_z), &(frame->boxaxis[6]), &(frame->boxaxis[7]), &(frame->boxaxis[8]));

		//frame->boxorg[0] -= (frame->boxaxis[0] + frame->boxaxis[3] + frame->boxaxis[6]) / 2.0;
		//frame->boxorg[1] -= (frame->boxaxis[1] + frame->boxaxis[4] + frame->boxaxis[7]) / 2.0;
		//frame->boxorg[2] -= (frame->boxaxis[2] + frame->boxaxis[5] + frame->boxaxis[8]) / 2.0;
		frame->boxaxis[0] *= (double)frame->grid_x;
		frame->boxaxis[1] *= (double)frame->grid_x;
		frame->boxaxis[2] *= (double)frame->grid_x;
		frame->boxaxis[3] *= (double)frame->grid_y;
		frame->boxaxis[4] *= (double)frame->grid_y;
		frame->boxaxis[5] *= (double)frame->grid_y;
		frame->boxaxis[6] *= (double)frame->grid_z;
		frame->boxaxis[7] *= (double)frame->grid_z;
		frame->boxaxis[8] *= (double)frame->grid_z;


		{
			for (int i = 0; i < pcnt; i++) {
				if (fgets(line, 1024, fp) == NULL) { fclose(fp); return nullptr; }	//原子位置読み込み				
			}
		}

		//ボリュームデータの読み込み
		size_t mesh_sz = (frame->grid_x) * (frame->grid_y) * (frame->grid_z);
		double* values = new double[mesh_sz];
		char* tp;
		double maxvalue = 0.0;
		double minvalue = 1e10;
		size_t ix = 0;
		size_t iy = 0;
		size_t iz = 0;
		const size_t iz_end = (frame->grid_z);
		const char* const endptr = line + SIZE - 1;
		while (fgets(line, SIZE, fp) != NULL) {
			
			const char* tp = strtok(line, " \t\r\n");
			while (tp){
				double vp = strtod(tp, nullptr);
				
				size_t index = ix + (frame->grid_x) * (iy + (frame->grid_y) * iz);
				values[index] = vp;


				//最大値の取得
				if (maxvalue < (vp)) { maxvalue = vp; }
				if (minvalue > (vp) && (vp > 0.0)) { minvalue = vp; }
				tp = strtok(nullptr, " \t\r\n");
				
				++iz;
				if (iz == iz_end) {
					iz = 0;
					++iy;
					break;
				}

				
				//res = std::from_chars(tp, endptr, vp);
			}

			if (iy == (frame->grid_y)) {
				iy = 0;
				++ix;
				if (ix == (frame->grid_x)) {
					break;
				}

			}
		}

		fclose(fp);


		return values;

	}
#endif

};
