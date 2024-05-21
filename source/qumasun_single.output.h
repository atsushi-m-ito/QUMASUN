#pragma once
#include "qumasun_single.h"
#include <cstdio>
#include "cube_writer.h"

inline
void QUMASUN_SINGLE::OutputDensity(QUMASUN::OUTPUT_TARGET target, const char* filepath) {
	using namespace QUMASUN; 
	
	CubeWriter::Frame frame;
	frame.grid_x = m_size_x;
	frame.grid_y = m_size_y;
	frame.grid_z = m_size_z;
	frame.boxaxis[0] = m_box_x;
	frame.boxaxis[1] = 0.0;
	frame.boxaxis[2] = 0.0;
	frame.boxaxis[3] = 0.0;
	frame.boxaxis[4] = m_box_y;
	frame.boxaxis[5] = 0.0;
	frame.boxaxis[6] = 0.0;
	frame.boxaxis[7] = 0.0;
	frame.boxaxis[8] = m_box_z;
	frame.boxorg[0] = 0.0;
	frame.boxorg[1] = 0.0;
	frame.boxorg[2] = 0.0;

	CubeWriter writer;

	
	switch (target) {
	case OUTPUT_TARGET::Density :
		writer.SaveCube(filepath, frame, m_num_nuclei, m_nuclei, m_rho);
		break;
	case OUTPUT_TARGET::Vhart:
		writer.SaveCube(filepath, frame, m_num_nuclei, m_nuclei, m_Vhart);
		break;
	}

}
