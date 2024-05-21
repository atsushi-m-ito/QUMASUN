#pragma once
#include "qumasun_input.h"
#include "JobParam.h"
#include "folding.h"
#include "atomic_number.h"

namespace QUMASUN {

	void SetSupercell(JobParam& job_params, Input& input)
	{
		//super-cellの読み込み//	
		int super_cell[3];
		job_params.GetIntArray("Material.SuperCell", super_cell, 3);
		//generator.SetSuperCell(super_cell[0], super_cell[1], super_cell[2]);
		const int num_total_cells = super_cell[0] * super_cell[1] * super_cell[2];

		//unit-lattice-vectorの読み込み//
		constexpr int LINE_SIZE = 1024;
		const char DELIMITER[] = " \t\r\n";
		char line[LINE_SIZE];
		struct V3 {
			double x, y, z;
		};
		V3 unit_vector[3];
		for (int i = 0; i < 3; ++i) {
			strcpy(line, job_params.GetBlockString("Material.UnitLatticeVector", i));

			char* p = strtok(line, DELIMITER);
			unit_vector[i].x = strtod(p, nullptr);

			p = strtok(nullptr, DELIMITER);
			unit_vector[i].y = strtod(p, nullptr);

			p = strtok(nullptr, DELIMITER);
			unit_vector[i].z = strtod(p, nullptr);
		}

		
		if (job_params.IsEqualString("Material.UnitLatticeVector.Unit", "Ang")) {
			//default is angstrome//
			for (int i = 0; i < 3; ++i) {
				unit_vector[i].x *= Ang_as_Bohr;
				unit_vector[i].y *= Ang_as_Bohr;
				unit_vector[i].z *= Ang_as_Bohr;
			}
		}

        double vacuum_cell_left[3] = { 0.0 };
        double vacuum_cell_right[3] = { 0.0 };
        if (job_params.GetDoubleArray("Material.VacuumCell.Leftside", vacuum_cell_left, 3)) {
            //generator.SetVacuumCellLeft(vacuum_cell[0], vacuum_cell[1], vacuum_cell[2]);
        }
        if (job_params.GetDoubleArray("Material.VacuumCell.Rightside", vacuum_cell_right, 3)) {
            //generator.SetVacuumCellRight(vacuum_cell[0], vacuum_cell[1], vacuum_cell[2]);
        }

        
		{
            double sx = (double)super_cell[0] + vacuum_cell_left[0] + vacuum_cell_right[0];
			input.box_axis[0] = sx * unit_vector[0].x;
			input.box_axis[1] = sx * unit_vector[0].y;
			input.box_axis[2] = sx * unit_vector[0].z;

            double sy = (double)super_cell[1] + vacuum_cell_left[1] + vacuum_cell_right[1];
			input.box_axis[3] = sy * unit_vector[1].x;
			input.box_axis[4] = sy * unit_vector[1].y;
			input.box_axis[5] = sy * unit_vector[1].z;

            double sz = (double)super_cell[2] + vacuum_cell_left[2] + vacuum_cell_right[2];
			input.box_axis[6] = sz * unit_vector[2].x;
			input.box_axis[7] = sz * unit_vector[2].y;
			input.box_axis[8] = sz * unit_vector[2].z;
		}




		//粒子数の読み込み//
		const int num_unit_atoms = job_params.GetBlockNumLines("Material.UnitCell");
		std::vector<Nucleus> atom_u_pos;
		for (int i = 0; i < num_unit_atoms; i++) {
			strcpy(line, job_params.GetBlockString("Material.UnitCell", i));

			char* p = strtok(line, DELIMITER);
			const int Z = msz::GetAtomicNumber(p);

			char* sx = strtok(nullptr, DELIMITER);
			char* sy = strtok(nullptr, DELIMITER);
			char* sz = strtok(nullptr, DELIMITER);

			//double unit_mass = ATOMIC_MASS[unit_element];

			//generator.AddIntoUnitcell(unit_element, unit_mass,
			atom_u_pos.emplace_back(Nucleus{ Z, strtod(sx, nullptr), strtod(sy, nullptr), strtod(sz, nullptr) });
		}

		if (job_params.IsEqualString("Material.UnitCell.Unit", "Ang")) {
			for (auto&& nuc : atom_u_pos) {
				nuc.Rx *= Ang_as_Bohr;
				nuc.Ry *= Ang_as_Bohr;
				nuc.Rz *= Ang_as_Bohr;
			}
		} else if (job_params.IsEqualString("Material.UnitCell.Unit", "FRAC")) {
			for (auto&& nuc : atom_u_pos) {
				const double rx = nuc.Rx * unit_vector[0].x + nuc.Ry * unit_vector[1].x + nuc.Rz * unit_vector[2].x;
				const double ry = nuc.Rx * unit_vector[0].y + nuc.Ry * unit_vector[1].y + nuc.Rz * unit_vector[2].y;
				const double rz = nuc.Rx * unit_vector[0].z + nuc.Ry * unit_vector[1].z + nuc.Rz * unit_vector[2].z;
				nuc.Rx = rx;
				nuc.Ry = ry;
				nuc.Rz = rz;
			}
		}


		input.num_nuclei = num_unit_atoms * num_total_cells;
		for (int cz = 0; cz < super_cell[2]; ++cz) {
			for (int cy = 0; cy < super_cell[1]; ++cy) {
				for (int cx = 0; cx < super_cell[0]; ++cx) {
					const double rx = (vacuum_cell_left[0] + (double)cx) * unit_vector[0].x + (vacuum_cell_left[1] + (double)cy) * unit_vector[1].x + (vacuum_cell_left[2] + (double)cz) * unit_vector[2].x;
					const double ry = (vacuum_cell_left[0] + (double)cx) * unit_vector[0].y + (vacuum_cell_left[1] + (double)cy) * unit_vector[1].y + (vacuum_cell_left[2] + (double)cz) * unit_vector[2].y;
					const double rz = (vacuum_cell_left[0] + (double)cx) * unit_vector[0].z + (vacuum_cell_left[1] + (double)cy) * unit_vector[1].z + (vacuum_cell_left[2] + (double)cz) * unit_vector[2].z;

					for (int k = 0; k < num_unit_atoms; ++k) {
						
						input.nuclei.emplace_back(Nucleus{ atom_u_pos[k].Z, 
							atom_u_pos[k].Rx + rx, +atom_u_pos[k].Ry + ry, +atom_u_pos[k].Rz + rz});
					}
				}
			}
		}
		
        for (auto&& nuc : input.nuclei) {
            Folding(nuc.Rx, nuc.Ry, nuc.Rz, input.box_axis);
        }
	}


	Input MakeInput(const char* input_file) {
		JobParam jobparam(input_file, "*.Begin", "*.End");

        //仕様変更したkeywordをreplace-keyとして登録(内部でキーを置換する)//
        jobparam.ReplaceKey("SCF.Initial.Density", "Initial.Density");
        jobparam.ReplaceKey("SCF.Initial.State", "Initial.State");
        jobparam.ReplaceKey("SCF.NumSolusions", "System.NumSolusions");


        constexpr int LINE_SIZE = 1024;
        const char DELIMITER[] = " \t\r\n";
        char line[LINE_SIZE];

		//condition of SCF solver///////////////
		Input input;
		input.scf_step = jobparam.GetInt("SCF.Step", 1);
		input.eigen_step_per_scf = jobparam.GetInt("SCF.EigenSolver.Step", 3);
		input.eigen_step_initial = jobparam.GetInt("SCF.EigenSolver.InitialStep", input.eigen_step_per_scf);
		const char* hamiltonian_name = jobparam.GetString("System.Hamiltonian");
		input.hamiltonian_type = ToHamiltonianType(hamiltonian_name  ? hamiltonian_name : "KohnSham_PP");

		//condition of system size///////////////		
		jobparam.GetIntArray("System.SpaceGrid", input.grid_size, 3);
        input.HR_ratio = jobparam.GetInt("System.HighResolution.Ratio", 1);

		//SetArray(input.box_axis, 20.0, 0.0, 0.0, 0.0, 20.0, 0.0, 0.0, 0.0, 20.0);


		const int num_elements = jobparam.GetBlockNumLines("System.AtomicElement");
		for (int i = 0; i < num_elements; ++i) {
			strcpy(line, jobparam.GetBlockString("System.AtomicElement", i));

			char* p = strtok(line, DELIMITER);
			const int Z = msz::GetAtomicNumber(p);

			p = strtok(nullptr, DELIMITER);
			input.pseudo_pot_set.emplace(std::make_pair(Z, p));

			p = strtok(nullptr, DELIMITER);
			if (p) {
				input.atomic_wave_set.emplace(std::make_pair(Z, p));
			} else {
				input.atomic_wave_set.emplace(std::make_pair(Z, ""));
			}
		}

		//condition of electrons and nuclei///////////////
		input.num_electrons = jobparam.GetInt("System.NumElectrons", 1);
        input.num_solutions = jobparam.GetInt("System.NumSolusions", 1);//alias
        SetSupercell(jobparam, input);
		
        if (jobparam.Find("Initial.Density")) {
            input.initial_density = jobparam.GetString("Initial.Density");
        } else {
			input.initial_density = "none";
		}
        if (jobparam.Find("Initial.DensityDifference")) {
            input.initial_density_difference = jobparam.GetString("Initial.DensityDifference");
        } else {
            input.initial_density_difference = "none";
        }

		input.spin_polarization = jobparam.GetInt("System.SpinPolarization.Mode", 0);
        input.initial_spin_difference = jobparam.GetInt("Initial.SpinDifference", 0);
		{
			auto res = jobparam.GetIntArray("System.KpointSample", input.kpoint_sample, 3);
			if (!res) {
				input.kpoint_sample[0] = 1; input.kpoint_sample[1] = 1; input.kpoint_sample[2] = 1;
			};

			if (jobparam.IsEqualString("System.KpointSymmetry", "FULL_XYZ")) {
				input.kpoint_symmetry = QUMASUN::KPOINT_SYMMETRY::FULL_XYZ;
			} else {
				input.kpoint_symmetry = QUMASUN::KPOINT_SYMMETRY::NONE;
			}
		}

        
		int num_init_state = jobparam.GetBlockNumLines("Initial.State");
		for (int i = 0; i < num_init_state; ++i) {
			input.initial_state.push_back(jobparam.GetBlockString("Initial.State", i));
		}

        input.initial_expand_from_kpoint = jobparam.GetInt("Initial.ExpandFromKpoint", 0);

		input.temperature_K = jobparam.GetDouble("SCF.Temperature", 300.0);
		input.scf_mixing_ratio = jobparam.GetDouble("SCF.MixingRatio", 0.25);
        
        //Dynamics////////////////////////////////////////////////////////
        if (jobparam.IsEqualString("Dynamics.Mode", "Test")) {
            input.dynamics_mode = DynamicsMode::Test;
        } else if (jobparam.IsEqualString("Dynamics.Mode", "Relaxation")) {
            input.dynamics_mode = DynamicsMode::Relaxation;
        } else if (jobparam.IsEqualString("Dynamics.Mode", "TDDFT")) {
            input.dynamics_mode = DynamicsMode::TDDFT;
        } else if (jobparam.IsEqualString("Dynamics.Mode", "EhrenfestMD")) {
            input.dynamics_mode = DynamicsMode::EhrenfestMDv1;
        }

        input.dynamics_step = jobparam.GetInt("Dynamics.Step", 1);
        input.dynamics_start_step = jobparam.GetInt("Dynamics.Start.Step", 0);
        input.dynamics_output_step = jobparam.GetInt("Dynamics.Output.Step", 100);
        if (input.dynamics_output_step < 1) {
            input.dynamics_output_step = 1;
        }
        input.dynamics_force_threshold = jobparam.GetDouble("Dynamics.Force.Threshold", 0.01);
        input.dynamics_timestep_dt = jobparam.GetDouble("Dynamics.TimeStep", 0.01);

        //mass and velocity///////////////////////////////////////////////////
        double mass_as_au = 1.0;
        if (jobparam.IsEqualString("Dynamics.Mass.Unit", "me") || jobparam.IsEqualString("Dynamics.Mass.Unit", "a.u.")) {
            mass_as_au = 1.0;
        } else if (jobparam.IsEqualString("Dynamics.Mass.Unit", "u") || jobparam.IsEqualString("Dynamics.Mass.Unit", "Da")) {
            mass_as_au = mass_u_as_au;
        }

        double velocity_as_au = 1.0;
        if (jobparam.IsEqualString("Dynamics.Velocity.Unit", "a.u.")) {
            velocity_as_au = 1.0;
        } else if (jobparam.IsEqualString("Dynamics.Velocity.Unit", "m/s") ) {
            velocity_as_au = velocity_ms_as_au;
        }


        //粒子数の読み込み//
        const int num_line_velocity = jobparam.GetBlockNumLines("Dynamics.Velocity");
        for (int i = 0; i < num_line_velocity; i++) {
            strcpy(line, jobparam.GetBlockString("Dynamics.Velocity", i));

            char* p = strtok(line, DELIMITER);
            input.mass.emplace_back(strtod(p, nullptr)* mass_as_au);

            char* sx = strtok(nullptr, DELIMITER);
            char* sy = strtok(nullptr, DELIMITER);
            char* sz = strtok(nullptr, DELIMITER);

            //double unit_mass = ATOMIC_MASS[unit_element];

            //generator.AddIntoUnitcell(unit_element, unit_mass,
            input.velocity.emplace_back(strtod(sx, nullptr)* velocity_as_au);
            input.velocity.emplace_back(strtod(sy, nullptr)* velocity_as_au);
            input.velocity.emplace_back(strtod(sz, nullptr)* velocity_as_au);
                    
        }


        double position_as_au = 1.0;
        if (jobparam.IsEqualString("Material.UnitLatticeVector.Unit", "Ang")) {
            position_as_au = Ang_as_Bohr;
        }

        //波動関数に初速を与えるルール//
        const int num_line_added_velocity = jobparam.GetBlockNumLines("Initial.State.Option");
        for (int i = 0; i < num_line_added_velocity; i++) {
            strcpy(line, jobparam.GetBlockString("Initial.State.Option", i));

            char* key = strtok(line, DELIMITER);
            //input.mass.emplace_back(strtod(p, nullptr) * mass_as_au);
            char* opecode = strtok(nullptr, DELIMITER);
            if (std::strncmp(opecode, "AddV", 4) != 0) continue;

            char* sx = strtok(nullptr, DELIMITER);
            char* sy = strtok(nullptr, DELIMITER);
            char* sz = strtok(nullptr, DELIMITER);

            //generator.AddIntoUnitcell(unit_element, unit_mass,
            AddedVelocityForWave add_to_wave;
            add_to_wave.velocity_x = (strtod(sx, nullptr) * velocity_as_au);
            add_to_wave.velocity_y = (strtod(sy, nullptr) * velocity_as_au);
            add_to_wave.velocity_z = (strtod(sz, nullptr) * velocity_as_au);

            sx = strtok(nullptr, DELIMITER);
            sy = strtok(nullptr, DELIMITER);
            sz = strtok(nullptr, DELIMITER);

            //generator.AddIntoUnitcell(unit_element, unit_mass,
            add_to_wave.center_x = (strtod(sx, nullptr) * position_as_au);
            add_to_wave.center_y = (strtod(sy, nullptr) * position_as_au);
            add_to_wave.center_z = (strtod(sz, nullptr) * position_as_au);
            input.velocity_for_wave.emplace(key, add_to_wave);
        }


        if (jobparam.Find("Output.StateVector")) {
            input.output_state_vector = jobparam.GetString("Output.StateVector");
        } else {
            input.output_state_vector = "none";
        }





		return input;
	}
}

