#ifdef USE_MPI
#pragma once
#include <complex>
#include <mpi.h>
#include "wave_function.h"
#include "nucleus.h"
#include "vps_loader.h"
#include "SubspaceField.h"
#include "PseudoPotOperator_mpi_v2.h"
#include "qumasun_input.h"
#include "GridRange.h"
#include "physical_param.h"
#include "vec3.h"
#include "StopWatch.h"



class QUMASUN_MPI
{
	using HAMILTONIAN = QUMASUN::HAMILTONIAN;
	using Input = QUMASUN::Input;
	using PseudoPotSet = QUMASUN::PseudoPotSet;
	using AtomicWaveSet = QUMASUN::AtomicWaveSet;
	using SOLVER = QUMASUN::SOLVER;
	using dcomplex = std::complex<double>;
	
	using Field = double*;

public:
	QUMASUN_MPI(const Input& input, const MPI_Comm& mpi_comm_, const int split_num[3]);
	~QUMASUN_MPI();

	void Execute(int dynamics_step = 0);



private:
	bool is_construction_successful = false;
	
	const SOLVER m_solver;
	const HAMILTONIAN m_hamiltonian_type;

	const int num_solution;//number of solution//
	const int LIMIT_SCF_STEP;//
	const int LOBPCG_STEP_PER_SCF;
	const int LOBPCG_STEP_FIRST;

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
    //High resolution grid for nuclear charge//
    const int m_HR_ratio_x = 2;
    const int m_HR_ratio_y = 2;
    const int m_HR_ratio_z = 2;


	//parameters for atoms and electrons/////////
	const bool is_spin_on;
	bool has_up_spin = true;
	int m_num_spin = 1;  //if spin polarization is calculated, m_num_spin becomes 2//
	const int num_electrons;
	const int m_num_nuclei;
	Nucleus* m_nuclei = nullptr;
	double* m_nuclei_valence_elecron = nullptr;
	double m_diff_Coulomb_Vlocal; //Zi*Zk/r - Zi*Vlocal_k(r)を事前に計算して格納する//
    double m_Enn_close_correction;//Zi*Zk/r - \int rho*Vlocal_k(r)を事前に計算して格納する//
    double m_Ecore_self;          //自己相互作用 Zi*Vlocal_i(r=0)を事前に計算して格納する//
    double m_Ecore_self_v2;       //自己相互作用 Zi*Vlocal_i(r=0)を事前に計算して格納する//
    double m_Ecore_self_HR;       //HR版:自己相互作用 Zi*Vlocal_i(r=0)を事前に計算して格納する//

	//parameter for DFT//////////////////////////
	double m_kbT = 300.0 * KbHartree;// 0.03;
	double m_mixing_ratio = 0.25;

	//Data common to all spins and k-points//////////////////////////////
	//RspaceFunc<double> m_Vtot;
	//RspaceFunc<double> m_Vtot_down;
	RspaceFunc<double> m_Vext;
	RspaceFunc<double> m_rho;
	RspaceFunc<double> m_rho_diff;
	RspaceFunc<double> m_rho_prev;
	RspaceFunc<double> m_Vhart;
	//RspaceFunc<double> m_Vx;
	//RspaceFunc<double> m_Vc;
	//RspaceFunc<double> m_Vx_down;
	//RspaceFunc<double> m_Vc_down;
    RspaceFunc<double> m_hr_Vext;
    RspaceFunc<double> m_hr_rho;
    RspaceFunc<double> m_hr_Vhart;

	//slave nodes have data for each spin and k-point. 
	//master node has data for all spins and k-points.
	double* m_eigen_values = nullptr;
	double* m_eigen_values_down = nullptr;   //the second buffer for spin down//
	double* m_occupancy = nullptr;
	double* m_occupancy_down = nullptr;      //the second buffer for spin down//

	//data for the case with Pseudo Potential///
	RspaceFunc<double> m_pcc_rho;
	////////////////////////////////////////////

	size_t m_work_size = 0;
	double* m_work = nullptr;
	
	//for Pseudo Potential//
	
	PseudoPotSet m_pseudo_pot_set;
	AtomicWaveSet m_atomic_wave_set;

	PseudoPotIntegrator_mpi m_pp_integrator;

	//way to set initial wave function//
	std::string m_initial_density;

    //force///////////////
    std::vector<vec3d> m_nucl_forces;
    std::vector<vec3d> m_force_nn_correction;
    std::vector<vec3d> m_force_nn_TF;

	
private:
	bool mCheckConditions();
	void mInitializeBuffer(const Input& input);
	void mInitializeState();
	void mInitializeDensity();
	void mFinalize();
	void mInitStateArnoldi();

	void mPrepareCore();

	double mSetDensity(bool is_mixing);
	void mSetDensityOne(GridRange& l_grid, double* l_rho, const double* const* l_psi, const double* occupancy, int num_solution);
	void mSetPotential();
	void mSetPotentialVext(RspaceFunc<double>& V, const Nucleus* nuclei, int num_nuclei); //common between PW & RS//
	//void mSetPotentialVhart(RspaceFunc<double>& Vhart, RspaceFunc<double>& rho);          //common between PW & RS//
	void mSetPotentialVxc(RspaceFunc<double>& Vx, RspaceFunc<double>& Vc, RspaceFunc<double>& rho, RspaceFunc<double>& pcc_rho); //common between PW & RS//
	void mSetPotentialVxcSpin(double* Vx_up, double* Vc_up, double* Vx_down, double* Vc_down, const double* rho, const double* pcc_rho, const double* rho_diff);
    void mSetPotentialVxc_local(double* Vx, double* Vc, const double* rho, const double* pcc_rho);
    void mSetPotentialVxcSpin_local(double* Vx_up, double* Vc_up, double* Vx_down, double* Vc_down, const double* rho, const double* pcc_rho, const double* rho_diff);
    void mCheckPotentialVhart(RspaceFunc<double>& Vhart, RspaceFunc<double>& rho);  //common between PW & RS//
	void mSetLaplasianFFT(RspaceFunc<double>& out, RspaceFunc<double>& in, double coef, double* work);          //common between PW & RS//
	//void mSetPotentialVlocalFromCorrespondingRho(RspaceFunc<double>& Vext, RspaceFunc<double>& nucl_rho);

	void mSolveLanczos();
	void mSolveLOBPCG(int steps, int SCF_current_step);
	void mCheckEigenVector();

	//void mHamiltonianMatrix(RspaceFunc<double >& Hp, const RspaceFunc<double>& V, const RspaceFunc<double>& p);
	void mHamiltonianMatrix_ddm(double* l_Hp, const double* l_Vtot, const double* l_phi);
	void mHamiltonianMatrix_bundle_ddm(int num_states, double* l_Hp, const double* l_Vtot, const double* l_phi, double* l_tmp);
	//void mKineticMatrixAdd(RspaceFunc<double>& Kp, const RspaceFunc<double>& p);
	void mKineticMatrixAdd_ddm(double* Kp, const double* p);
	//void mKineticMatrixAdd(KspaceFunc<dcomplex >& p);
	void mPotentialMatrix(RspaceFunc<double>& Vp, const RspaceFunc<double>& V, const RspaceFunc<double >& p);
	void mPotentialMatrix_ddm(double* Vp, const double* V, const double* p);

	void mSetOccupancy();//common between PW & RS//
	void mSetOccupancySpinOn();//common between PW & RS//
	void mSetOccupancySpinOff();//common between PW & RS//

	double mGetTotalEnergy();
	double mGetEnergyKinetic();

	double mGetEnergy_V_rho(const double* V, const double* rho);
    double mGetEnergy_V_rho_HR(const double* V, const double* rho);
	double mGetEnergyExtByVextRho();
    double mGetEnergyExtByVextRho_HR(); 
	double mGetEnergyHartree(); 
	double mGetEnergyExtByVhartNuclRho();
    double mGetEnergyExtByVhartNuclRho_HR();
	double mGetEnergyExtByVhartAtPoint();
	//double mGetEnergyXC(double* pEx, double* pEc);
	double mGetEnergyXC(double* pEx, double* pEc, double* pVxRho, double* pVcRho);
	double mGetEnergyXCSpin(double* pEx, double* pEc, double* pVxRho_up, double* pVcRho_up, double* pVxRho_down, double* pVcRho_down);
	double mGetEnergyPseudoNonlocal();

	double mGetEnergyCoreCore_direct();
	double mGetEnergyPotentialAtNucl(const double* V);
	double mGetEnergyCoreCoreByVextAtPoint();
	double mGetEnergyCoreCoreNuclDensity();
    double mGetEnergyCoreCoreNuclDensity_HR();

	void mGetForce();
	void mGetForceOnNuclRho(double* force, const double* l_Vpot);
    void mGetForceOnNuclRho_HR(double* force, const double* l_hr_Vpot);
	void mGetForceHartreeNuclRho(double* force);
    void mGetForceHartreeNuclRho_HR(double* force);
    void mGetForceCoreCoreByVextNuclRho(double* force);
    void mGetForceCoreCoreByVextNuclRho_HR(double* force);
	double mGetForcePseudoNonlocal(double* force);



	inline double Length(double x, double box_w) { return (x > box_w*0.5) ? x - box_w : (x < -box_w*0.5) ? x + box_w : x; };
	


	//for mpi////////////////////
private:
	// DDMにおいて分散したローカルデータはml_で始まる変数とする//
	const int m_root_id = 0;
	MPI_Comm m_mpi_comm;
	MPI_Comm m_ddm_comm;
	GridRange m_global_grid;      //common to all spins and k-points
	GridRangeMPI ml_grid;         //common to all spins and k-points
	int m_mpi_split_color;        //color of MPI_Comm_split, which corresponds to process group and is used ScaLAPACK multi-grid//
	
	
	Field* ml_psi = nullptr;      //data for each spin and k-point//
	Field ml_lobpcg_keep_p = nullptr;      //data for each spin and k-point//
	double* m_lobpcg_keep_S_matrix = nullptr;//data for each spin and k-point//
	Field ml_rho = nullptr;       //common to all spins and k-points, used by calculation of density and energy .
	Field ml_Vtot = nullptr;      //depends on spin//
	Field ml_Vtot_down = nullptr; //depends on spin//
	Field ml_Vx = nullptr;        //depends on spin//
	Field ml_Vc = nullptr;        //depends on spin//	
    Field ml_Vx_down = nullptr;        //depends on spin//
    Field ml_Vc_down = nullptr;        //depends on spin//	
	Field ml_Vext = nullptr;      //common to all spins and k-points
	Field ml_Vhart = nullptr;     //common to all spins and k-points
	Field ml_pcc_rho = nullptr;   //common to all spins and k-points, used by calculation of energy.
	Field ml_rho_diff = nullptr;  //common to all spins and k-points, used by calculation of energy.
	Field ml_nucl_rho = nullptr;  //common to all spins and k-points, used by calculation of energy.
    Field ml_hr_nucl_rho = nullptr;  //common to all spins and k-points, used by calculation of energy.
    Field ml_hr_Vext = nullptr;      //common to all spins and k-points
    Field ml_hr_rho = nullptr;      //common to all spins and k-points
    Field ml_hr_Vhart = nullptr;      //common to all spins and k-points
	
	void mScatterField(double* local_dest, const double* global_src);
	void mGatherField(double* global_dest, const double* local_src);
    void mScatterField_HR(double* local_dest, const double* global_src);
    void mGatherField_HR(double* global_dest, const double* local_src);
    void mHierarchyScatterField(double* local_dest, const double* global_src);

	//for benchmark/////////////////////////////////////////////
	StopWatch<80, true> watch;
	void mPrintTime();

public:

	void OutputDensity(QUMASUN::OUTPUT_TARGET target, const char* filepath);
	void OutputEigenValue(const char* filepath);
    void OutputEigenVector(const char* filepath) {};

	void PrintCondition();
    void MoveNuclei(int num_nuclei, const Nucleus* next_nucleis);
    void GetForce(vec3d* forces);
};

#include "qumasun_mpi.constructor.h"
#include "qumasun_mpi.core.h"
#include "qumasun_mpi.density.h"
#include "qumasun_mpi.energy.h"
#include "qumasun_mpi.initialize.h"
#include "qumasun_mpi.lanczos.h"
#include "qumasun_mpi.lobpcg.h"
#include "qumasun_mpi.matrix.h"
#include "qumasun_mpi.occupancy.h"
#include "qumasun_mpi.potential.h"
#include "qumasun_mpi.scfloop.h"
#include "qumasun_mpi.scatter.h"
#include "qumasun_mpi.output.h"
#include "qumasun_mpi.force.h"
#include "qumasun_mpi.move.h"
#include "qumasun_mpi.print.h"

#endif
