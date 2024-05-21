#pragma once


constexpr const double eV_as_J = (1.602176462e-19);		//(J / eV)
constexpr const double Hartree_as_J = 4.3597447222071e-18;   //
constexpr const double Hartree_as_eV = Hartree_as_J / eV_as_J;// 27.21138;
//constexpr const double JHartree = (JeV * eV_Hartree);		//(J / Hartree)
constexpr const double Kb = (1.3806503e-23);			//Boltzmann constant (J / K)
constexpr const double KbeV = (Kb / eV_as_J);			//coefficient (eV / K)
constexpr const double KbHartree = (Kb / Hartree_as_J);			//coefficient (Hartree / K)
constexpr const double Bohr_as_Ang = 0.529177210903;
constexpr const double Ang_as_Bohr = 1.0 / Bohr_as_Ang;
constexpr const double time_au  = 2.4188843265857e-17; // \hbar / Hartree, where \hbar is action in atomic unit
constexpr const double mass_u_as_kg = 1.66053906660e-27;
constexpr const double mass_au_as_kg = 9.1093837015e-31;
constexpr const double mass_u_as_au = mass_u_as_kg / mass_au_as_kg;
constexpr const double velocity_ms_as_au = time_au /Bohr_as_Ang * 1.0e+10;

//macro//
//note: 0.0063497943974046[a.u.] = std::sqrt(2.0 / (Hartree_as_eV * mass_u_as_au))
constexpr const double velocity_1eV_1u_as_au = 0.0063497933278141;
//note: 13891.38738921143[m/s] = std::sqrt(2.0 * JeV / mass_u_per_kg)
constexpr const double velocity_1eV_1u_as_ms = 13891.38738921143;


