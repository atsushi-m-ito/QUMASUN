#pragma once
#ifndef md3_writer2_h     	//version 0001//
#define md3_writer2_h

#include <stdio.h>
#include <stdlib.h>
#include "vec3.h"
#include "atomic_number.h"

namespace msz {

	template <typename T>
	class MD3_Writer2 {
	private:
		FILE* m_fp;	//ファイル
					
		int m_num_atoms;	//粒子数.
		int m_frameno;

	public:


		MD3_Writer2() : m_fp(nullptr), m_frameno(0){}
		virtual ~MD3_Writer2() {
			Close();
		}

		///////////////////////////////////////////////////////
		//ファイルをオープンする.
		//	正しくオープンできたときは true を返す.
		/////////////////////////////////////////////////////////////
		bool Open(const char* filepath) {
			if (m_fp) {
				return false;
			}

			m_fp = fopen(filepath, "w");

			if (!m_fp) {
				return false;
			}

			return true;
		}


		///////////////////////////////////////////////////////
		//ファイルを閉じる.
		//	destructerからもこの関数が自動で呼ばれる.
		///////////////////////////////////////////////////////////
		void Close() {
			if (m_fp) { fclose(m_fp); }
			m_fp = nullptr;
		}

	private:

		const int* m_Z;
		//const int* m_id;
		const vec3<T>* m_r;
		const T* m_rx;
		const T* m_ry;
		const T* m_rz;
		const vec3<T>* m_p;
		const T* m_px;
		const T* m_py;
		const T* m_pz;
		const vec3<T>* m_f;
		const T* m_fx;
		const T* m_fy;
		const T* m_fz;
		const T* m_box_axis_origin;
		const char* m_text;
		const T* m_stress;

	public:
		////////////////////////////////////////////////////////////////////////////////
		//必要なデータをファイルに書き出す.
		//
		//  書き出したいデータの種類は,BeginFrame()とEndFrame()の間で,WriteR/P/F等を呼ぶことで指定する
		//  呼ばなかったデータは書き出さない。
		//　Userの指定したバッファに実際にデータが格納されるのはEndFrame()をコールした時点.
		//
		///////////////////////////////////////////////////////////////////////////////

		void BeginFrame(int count) {
			m_num_atoms = count;
			m_Z = nullptr;
			//m_id = nullptr;
			m_r = nullptr;
			m_rx = nullptr;
			m_ry = nullptr;
			m_rz = nullptr;
			m_p = nullptr;
			m_px = nullptr;
			m_py = nullptr;
			m_pz = nullptr;
			m_f = nullptr;
			m_fx = nullptr;
			m_fy = nullptr;
			m_fz = nullptr;
			m_box_axis_origin = nullptr;
			m_text = nullptr;
			m_stress = nullptr;
		}

		void WriteZ(const int* knd) {
			m_Z = knd;
		}

		void WriteR(const vec3<T>* r) {
			m_r = r;
		}

		void WriteR_SoA(const T* rx, const T* ry, const T* rz) {
			m_rx = rx;
			m_ry = ry;
			m_rz = rz;
		}

		void WriteP(const vec3<T>* p) {
			m_p = p;
		}

		void WriteP_SoA(const T* px, const T* py, const T* pz) {
			m_px = px;
			m_py = py;
			m_pz = pz;
		}

		void WriteF(const vec3<T>* f) {
			m_f = f;
		}

		void WriteF_SoA(const T* fx, const T* fy, const T* fz) {
			m_fx = fx;
			m_fy = fy;
			m_fz = fz;
		}

		void WriteBox(const T* box_axis_origin) {
			m_box_axis_origin = box_axis_origin;
		}

		void WriteStress(const T* stress) {
			m_stress = stress;
		}

		void WriteComment(const char* comment) {
			m_text = comment;
		}

		//////////////////////////////////////////////
		//
		//  実際の読み込みを行う
		//
		//  戻り値として読み込んだ粒子数を返す
		//  エラーの場合は負の値を返す
		//////////////////////////////////////////////
		int EndFrame() {

			if (!m_fp) { return -1; }

			//粒子数の出力//
			fprintf(m_fp, "%d\n", m_num_atoms);

			//コメント行の出力//
			if (m_text) {
				fprintf(m_fp, "%s\n", m_text);
			} else {
				fprintf(m_fp, "frame = %d\n", m_frameno);
			}

			//BOXの出力//
			if (m_box_axis_origin) {
				fprintf(m_fp, "BOX   %.15f  %.15f  %.15f  %.15f  %.15f  %.15f  %.15f  %.15f  %.15f\n",
					m_box_axis_origin[0], m_box_axis_origin[1], m_box_axis_origin[2],
					m_box_axis_origin[3], m_box_axis_origin[4], m_box_axis_origin[5],
					m_box_axis_origin[6], m_box_axis_origin[7], m_box_axis_origin[8]);
			} else {
				fprintf(m_fp, "BOX   %.15f  %.15f  %.15f  %.15f  %.15f  %.15f  %.15f  %.15f  %.15f\n",
					10.0, 0.0, 0.0,
					0.0, 10.0, 0.0,
					0.0, 0.0, 10.0);
			}

			for (int i = 0; i < m_num_atoms; i++) {
				if (m_r) {
					fprintf(m_fp, "%4s   %.15f  %.15f  %.15f",
						msz::GetAtomicSymbol(m_Z[i]), m_r[i].x, m_r[i].y, m_r[i].z);
				} else {
					fprintf(m_fp, "%4s   %.15f  %.15f  %.15f",
						msz::GetAtomicSymbol(m_Z[i]), m_rx[i], m_ry[i], m_rz[i]);
				}
				if(m_p){
					fprintf(m_fp, "  %.15f  %.15f  %.15f",
						m_p[i].x, m_p[i].y, m_p[i].z);
				} else if (m_px) {
					fprintf(m_fp, "  %.15f  %.15f  %.15f",
						m_px[i], m_py[i], m_pz[i]);
				} else if (m_f || m_fx) {
					fprintf(m_fp, "  0.0  0.0  0.0");
				}
				if (m_f) {
					fprintf(m_fp, "  %.15f  %.15f  %.15f\n",
						m_f[i].x, m_f[i].y, m_f[i].z);
				} else if (m_fx) {
					fprintf(m_fp, "  %.15f  %.15f  %.15f\n",
						m_fx[i], m_fy[i], m_fz[i]);
				} else {
					fprintf(m_fp, "\n");
				}
			} 

			m_frameno++;
			return m_frameno;
		}

	};



}//namespace msz//
#endif //!md3_writer2_h//
