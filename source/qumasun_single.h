#pragma once
#include <complex>
#include "wave_function.h"
#include "nucleus.h"
#include "vps_loader.h"
#include "SubspaceField.h"
#include "PseudoPotOperator.h"
#include "qumasun_input.h"
#include "physical_param.h"
#include "StopWatch.h"



class QUMASUN_SINGLE
{
	using HAMILTONIAN = QUMASUN::HAMILTONIAN;
	using Input = QUMASUN::Input;
	using PseudoPotSet = QUMASUN::PseudoPotSet;
	using AtomicWaveSet = QUMASUN::AtomicWaveSet;
	using SOLVER = QUMASUN::SOLVER;
	using dcomplex = std::complex<double>;
	


public:
	QUMASUN_SINGLE(const Input & input);
	~QUMASUN_SINGLE();

	void Execute();

	void OutputDensity(QUMASUN::OUTPUT_TARGET target, const char* filepath);

private:
	
	const SOLVER m_solver;
	const HAMILTONIAN m_hamiltonian_type;

	const int num_solution;//number of solution//
	const int LIMIT_SCF_STEP;//
	const int LOBPCG_STEP_PER_SCF;
	const int LOBPCG_STEP_FIRST;
	//LOBPCG_Solver lobpcg;
	//LOBPCG_Solver lobpcg_down;

	//parameters for system size and resolution//
	const int m_size_x;//
	const int m_size_y;//
	const int m_size_z;//
	const double m_box_x;
	const double m_box_y;
	const double m_box_z;

	//the following definitions should be written after m_size_x and m_box_x//
	const int m_size_3d;
	const double m_dx;
	const double m_dy;
	const double m_dz;
	const double m_volume;

	//parameters for atoms and electrons/////////
	const bool is_spin_on;
	const int m_num_spin = 1;  //if spin polarization is calculated, m_num_spin becomes 2//
	const int num_electrons;
	const int m_num_nuclei;
	Nucleus* m_nuclei;
	double* m_nuclei_valence_elecron;
	double m_diff_Coulomb_Vlocal; //Zi*Zk/r - Zi*Vlocal_k(r)を事前に計算して格納する//
	double m_Ecore_self;          //自己相互作用 Zi*Vlocal_i(r=0)を事前に計算して格納する//

	//parameter for DFT//////////////////////////
	const double m_kbT = 1000.0 * KbHartree;// 0.03;

	//data for DFT //////////////////////////////
	RspaceFunc<double> m_Vtot;
	RspaceFunc<double> m_Vtot_down;
	RspaceFunc<double> m_Vext;
	RspaceFunc<double> m_rho;
	RspaceFunc<double> m_rho_diff;
	RspaceFunc<double> m_rho_prev;
	RspaceFunc<double> m_Vhart;
	RspaceFunc<double> m_Vx;
	RspaceFunc<double> m_Vc;
	RspaceFunc<double> m_Vx_down;
	RspaceFunc<double> m_Vc_down;
	RspaceFunc<double>* m_psi_set;
	RspaceFunc<double> m_psi_buffer;
	
	double* m_eigen_values = nullptr;
	double* m_eigen_values_down = nullptr;   //the second buffer for spin down//
	double* m_occupancy = nullptr;
	double* m_occupancy_down = nullptr;      //the second buffer for spin down//
	//data for the case with Pseudo Potential///
	RspaceFunc<double> m_pcc_rho;
	RspaceFunc<double> m_nucl_rho;
	////////////////////////////////////////////

	double* m_lobpcg_keep_p = nullptr;
	double* m_work = nullptr;

	//for Pseudo Potential//
	
	PseudoPotSet m_pseudo_pot_set;
	AtomicWaveSet m_atomic_wave_set;

	PseudoPotIntegrator_single m_pp_integrator;
	
	//way to set initial wave function//
	
	std::string m_initial_density;

private:
	bool mCheckConditions();
	void mInitializeState();
	void mInitializeDensity();
	void mFinalize();

	void mPrepareCore();

	double mSetDensity(bool is_mixing);
	void mSetDensityOne(size_t size_3d, double* rho, RspaceFunc<double> const* psi_set, const double* occupancy, int num_solution);
	void mSetPotential();
	void mSetPotentialVext(RspaceFunc<double>& V, const Nucleus* nuclei, int num_nuclei); //common between PW & RS//
	void mSetPotentialVhart(RspaceFunc<double>& Vhart, RspaceFunc<double>& rho);          //common between PW & RS//
	void mSetPotentialVxc(RspaceFunc<double>& Vx, RspaceFunc<double>& Vc, RspaceFunc<double>& rho, RspaceFunc<double>& pcc_rho); //common between PW & RS//
	void mSetPotentialVxcSpin(double* Vx_up, double* Vc_up, double* Vx_down, double* Vc_down, const double* rho, const double* pcc_rho, const double* rho_diff);
	void mCheckPotentialVhart(RspaceFunc<double>& Vhart, RspaceFunc<double>& rho);  //common between PW & RS//
	void mSetLaplasianFFT(RspaceFunc<double>& out, RspaceFunc<double>& in, double coef);          //common between PW & RS//
	void mSetPotentialVlocalFromCorrespondingRho(RspaceFunc<double>& Vext, RspaceFunc<double>& nucl_rho, const Nucleus* nuclei, int num_nuclei);

	void mSolveLanczos();
	void mSolveLOBPCG(int steps, bool is_first_step);
	void mCheckEigenVector();

	void mHamiltonianMatrix(RspaceFunc<double >& Hp, const RspaceFunc<double>& V, const RspaceFunc<double>& p);
	void mKineticMatrixAdd(RspaceFunc<double>& Kp, const RspaceFunc<double>& p);
	//void mKineticMatrixAdd(KspaceFunc<dcomplex >& p);
	void mPotentialMatrix(RspaceFunc<double>& Vp, const RspaceFunc<double>& V, const RspaceFunc<double >& p);

	void mSetOccupancy();//common between PW & RS//
	void mSetOccupancySpinOn();//common between PW & RS//
	void mSetOccupancySpinOff();//common between PW & RS//

	double mGetTotalEnergy();
	double mGetEnergyKinetic();
	double mGetOneEnergyKinetic(RspaceFunc<double>& p);

	double mGetEnergy_V_rho(const double* V, const double* rho);
	double mGetEnergyExtByVextRho();
	double mGetEnergyHartree();
	double mGetEnergyExtByVhartNuclRho();
	double mGetEnergyExtByVhartAtPoint();
	//double mGetEnergyXC(double* pEx, double* pEc);
	double mGetEnergyXC(double* pEx, double* pEc, double* pVxRho, double* pVcRho);
	double mGetEnergyXCSpin(double* pEx, double* pEc, double* pVxRho_up, double* pVcRho_up, double* pVxRho_down, double* pVcRho_down);
	
	double mGetEnergyPseudoNonlocal();

	double mGetEnergyCoreCore_direct();
	double mGetEnergyPotentialAtNucl(const double* V);
	double mGetEnergyCoreCoreByVextAtPoint();
	double mGetEnergyCoreCoreNuclDensity();


	inline double Length(double x, double box_w) { return (x > box_w*0.5) ? x - box_w : (x < -box_w*0.5) ? x + box_w : x; };
	inline double SQ(double x) { return x * x; };


private:
	//for benchmark///////////////////////////////////////
	StopWatch<20, true> watch;
	void mPrintTime();
};

#include "qumasun_single.constructor.h"
#include "qumasun_single.core.h"
#include "qumasun_single.density.h"
#include "qumasun_single.energy.h"
#include "qumasun_single.initialize.h"
#include "qumasun_single.lanczos.h"
#include "qumasun_single.lobpcg.h"
#include "qumasun_single.matrix.h"
#include "qumasun_single.occupancy.h"
#include "qumasun_single.potential.h"
#include "qumasun_single.scfloop.h"
#include "qumasun_single.output.h"

