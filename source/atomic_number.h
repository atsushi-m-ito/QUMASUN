#pragma once
#ifndef atomic_number_h   	//version 0001//
#define atomic_number_h

#include <string.h>


namespace msz {

	
	//元素記号の文字列から原子番号を返す//
	inline const char* GetAtomicSymbol(int atomic_number) {
		//元素記号//
		static const char glips_species_list[119][3]{
			"E",		//empty軌道//
			"H", "He",	//第1周期
			"Li", "Be", "B", "C", "N", "O", "F", "Ne",		//第2周期
			"Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar",	//第3周期
			"K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr",	//第4周期
			"Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe",	//第5周期
			"Cs", "Ba",	//第6周期(Cs, Ba)
			"La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu",	////第6周期(ランタノイド)
			"Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn",//第6周期(Hf-Rn)
			"Fr", "Ra",	//第7周期(Fr, Ra)
			"Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr",	////第7周期(アクチノイド)
																										//104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,//第7周期(Rf-Uuo)
		};

		return glips_species_list[atomic_number];
	}

	//元素記号の文字列から原子番号を返す//
	inline int GetAtomicNumber(const char* atomicsymbol) {
		constexpr int SPECIES_COUNT = 118 + 1;	

		for (int i = 0; i < SPECIES_COUNT; i++) {
			//if (strcmp(atomicsymbol, glips_species_list[i]) == 0) {
			if (strcmp(atomicsymbol, GetAtomicSymbol(i)) == 0) {			
				return i;
			}
		}
		return -1;
	}
}//namespae msz//

#endif //!atomic_number_h//
