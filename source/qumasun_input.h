#pragma once
#include <complex>
#include <string>
#include <map>
#include <vector>
#include "wave_function.h"
#include "PseudoPotOperator.h"
#include "vps_loader.h"
#include "nucleus.h"



namespace QUMASUN {

	using PseudoPotSet = std::map<int, std::string>;
	using AtomicWaveSet = std::map<int, std::string>;
	enum class SOLVER {
		Lanczos, LOBPCG
	};

	enum class HAMILTONIAN {
		Unsupported, Schrodinger, KohnSham_AE, KohnSham_PP, KohnSham_PAW,
	};


	static inline
	HAMILTONIAN ToHamiltonianType(const std::string& word) {
		if (word == "Schrodinger") return HAMILTONIAN::Schrodinger;
		if (word == "KohnSham_AE") return HAMILTONIAN::KohnSham_AE;
		if (word == "KohnSham_PP") return HAMILTONIAN::KohnSham_PP;
		if (word == "KohnSham_PAW") return HAMILTONIAN::KohnSham_PAW;
		return HAMILTONIAN::Unsupported;
	}

	namespace KPOINT_SYMMETRY{
		static constexpr uint32_t NONE = 0;
		static constexpr uint32_t NEGAPOSI_X = 0x1;
		static constexpr uint32_t NEGAPOSI_Y = 0x2;
		static constexpr uint32_t NEGAPOSI_Z = 0x4;
		static constexpr uint32_t MIRROR_XY = 0x8;
		static constexpr uint32_t MIRROR_YZ = 0x10;
		static constexpr uint32_t MIRROR_ZX = 0x20;
		static constexpr uint32_t FULL_XYZ = 0x3F;
	};


    enum DynamicsMode : int {
        None = 0,
        Test = 99,
        Relaxation = 1,
        //note: time dependent Khon-Sham is sed when value is 1000 or higher//
        TDDFT = 1000,
        EhrenfestMDv1 = 2001,
        EhrenfestMDv2 = 2002,
    };

    struct AddedVelocityForWave {//原子と共に動く場合など、波動関数に初速度を与えるためのもの//
        double center_x;
        double center_y;
        double center_z;
        double velocity_x;
        double velocity_y;
        double velocity_z;
    };
    
    
	struct Input {
		//condition of SCF solver///////////////
		int scf_step;
		int eigen_step_per_scf;
		int eigen_step_initial;
		HAMILTONIAN hamiltonian_type;
		PseudoPotSet pseudo_pot_set;
		AtomicWaveSet atomic_wave_set;
		//condition of system size///////////////
		int grid_size[3];
		double box_axis[9];
		//condition of electrons and nuclei///////////////
		int num_solutions;
		int num_electrons;
		int num_nuclei;
		std::vector<Nucleus> nuclei;
		int spin_polarization = 0;
        int initial_spin_difference = 0;
        int HR_ratio = 1;
		int kpoint_sample[3]{ 1,1,1 };
		uint32_t kpoint_symmetry= KPOINT_SYMMETRY::NONE;
		double temperature_K = 300.0;
		double scf_mixing_ratio = 0.5;

        std::string initial_density;
        std::string initial_density_difference;
        std::vector < std::string> initial_state;
        int initial_expand_from_kpoint=0;

        DynamicsMode dynamics_mode = DynamicsMode::None;
        int dynamics_step = 1;
        int dynamics_output_step = 1;
        int dynamics_start_step = 0;
        double dynamics_timestep_dt = 0.01;
        double dynamics_force_threshold = 0.01;
        std::string output_state_vector;
        std::vector<double> mass;
        std::vector<double> velocity;
        std::map<std::string, AddedVelocityForWave> velocity_for_wave;
	};


	enum OUTPUT_TARGET {
		Density, Vhart, DensityHR
	};

}

