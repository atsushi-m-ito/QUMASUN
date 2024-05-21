/*******************************

QUantum MAterial Simulation UNraveler (QUMASUN)

QUMASUN is numerical simulation code for Density Functional Theory(DFT) and Time-dependent DFT based on the real space grid.

********************************/

#include <cstdlib>
#include <cstdio>
#include "qumasun_single.h"
#include "atomic_number.h"
#include "field_interpolation.h"
#include "qumasun_make_input.h"



//simple DFT

int main(int argc, char* argv[]) {

	if (argc < 2) {
		printf("ERROR: no input file.\n");
		return -1;
	}

	printf("Begin calculation on QUMASUN\n"); fflush(stdout);
	
	QUMASUN::Input input = QUMASUN::MakeInput(argv[1]);

	//begin calculation///////////////////////////////
	QUMASUN_SINGLE qumasun(input);
	qumasun.Execute();
	qumasun.OutputDensity(QUMASUN::OUTPUT_TARGET::Density, "test.cube");
	qumasun.OutputDensity(QUMASUN::OUTPUT_TARGET::Vhart, "test_Vhart.cube");
}





