#ifdef USE_MPI
#pragma once
#include <complex>
#include <mpi.h>
#include "wave_function.h"
#include "nucleus.h"
#include "vps_loader.h"
#include "SubspaceField.h"
#include "PseudoPotOperator_kpoint_v2.h"
#include "qumasun_input.h"
#include "GridRange.h"
#include "soacomplex.h"
#include "physical_param.h"
#include "vec3.h"
#include "StopWatch.h"



class QUMASUN_KPOINT
{
	using HAMILTONIAN = QUMASUN::HAMILTONIAN;
	using Input = QUMASUN::Input;
	using PseudoPotSet = QUMASUN::PseudoPotSet;
	using AtomicWaveSet = QUMASUN::AtomicWaveSet;
	using SOLVER = QUMASUN::SOLVER;
	using dcomplex = std::complex<double>;
	
	using Field = double*;


public:
	QUMASUN_KPOINT(const Input& input, const MPI_Comm& mpi_comm_, const int ddm_num[3]);
	~QUMASUN_KPOINT();

	void Execute(int dynamics_step = 0);
	
private:
	bool is_construction_successful = false;

	const SOLVER m_solver;
	const HAMILTONIAN m_hamiltonian_type;

	const int num_solution;//number of solution//
    int m_initial_spin_differnce = 0;
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

	const double m_dkx;
	const double m_dky;
	const double m_dkz;

    //High resolution grid for nuclear charge//
    const int m_HR_ratio_x = 2;
    const int m_HR_ratio_y = 2;
    const int m_HR_ratio_z = 2;


	//parameters for atoms and electrons/////////
	const bool is_spin_on;
	//int m_num_spin = 1;  //if spin polarization is calculated, m_num_spin becomes 2//
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

	//data for the case with Pseudo Potential///
	RspaceFunc<double> m_pcc_rho;
	////////////////////////////////////////////

	size_t m_work_size = 0;
	double* m_work = nullptr;
	
	//for Pseudo Potential//
	
	PseudoPotSet m_pseudo_pot_set;
	AtomicWaveSet m_atomic_wave_set;

	PseudoPotIntegrator_kpoint m_pp_integrator;

	//way to set initial wave function//
	std::string m_initial_density;
    std::string m_initial_density_difference;
	std::vector<std::string> m_initial_state;
    bool m_initial_expand_from_kpoint = false;

	int m_all_kinds_spin_kpoint = 1;  //all spin and k sampling points(Born and von Karman numbers)//
	int m_kpoint_sampling[3]{ 1,1,1 }; //(Born and von Karman numbers)
	uint32_t kpoint_symmetry;
	
    //force///////////////
    std::vector<vec3d> m_nucl_forces;
    std::vector<vec3d> m_force_nn_correction;
    std::vector<vec3d> m_force_nn_TF;

private:
	bool mCheckConditions();
	void mInitializeBuffer(const Input& input);
	void mInitializeState();
	void mInitializeStateRandom();
	void mInitializeStateFile();
	void mInitializeDensity();
	void mFinalize();

	void mPrepareCore();

	double mSetDensity(bool is_mixing);
	void mSetDensityOne(GridRange& l_grid, double* l_rho, const SoAComplex* l_psi, const double* occupancy, int num_solution);
	void mSetPotential();
	void mSetPotentialVext(RspaceFunc<double>& V, const Nucleus* nuclei, int num_nuclei); //common between PW & RS//
	//void mSetPotentialVhart(RspaceFunc<double>& Vhart, RspaceFunc<double>& rho);          //common between PW & RS//
	void mSetPotentialVxc(RspaceFunc<double>& Vx, RspaceFunc<double>& Vc, RspaceFunc<double>& rho, RspaceFunc<double>& pcc_rho); //common between PW & RS//
	void mSetPotentialVxcSpin(double* Vx_up, double* Vc_up, double* Vx_down, double* Vc_down, const double* rho, const double* pcc_rho, const double* rho_diff);
    void mSetPotentialVxc_local(double* Vx, double* Vc, const double* rho, const double* pcc_rho);
    void mSetPotentialVxcSpin_local(double* Vx_up, double* Vc_up, double* Vx_down, double* Vc_down, const double* rho, const double* pcc_rho, const double* rho_diff);
    void mCheckPotentialVhart(RspaceFunc<double>& Vhart, RspaceFunc<double>& rho);  //common between PW & RS//
	void mSetLaplasianFFT(RspaceFunc<double>& out, RspaceFunc<double>& in, double coef, double* work);          //common between PW & RS//
	//void mSetPotentialVlocalFromCorrespondingRho(double* Vext, const double* nucl_rho, double* work);

	//void mSolveLanczos();
	void mSolveLOBPCG(int steps, int SCF_current_step);
	//void mCheckEigenVector();

	//void mHamiltonianMatrix(RspaceFunc<double >& Hp, const RspaceFunc<double>& V, const RspaceFunc<double>& p);
	void mHamiltonianMatrix_ddm(SoAComplex& l_Hp, const double* l_Vtot, const SoAComplex& l_phi, int kpoint_x, int kpoint_y, int kpoint_z, int id_spin_kpoint);
	//void mKineticMatrixAdd(RspaceFunc<double>& Kp, const RspaceFunc<double>& p);
	void mKineticMatrixAdd_ddm(double* Kp, const double* p);
	void mKineticMatrixAdd_ddm_kpint(SoAComplex& Kp, const SoAComplex& p, int kpoint_x, int kpoint_y, int kpoint_z);
	//void mKineticMatrixAdd(KspaceFunc<dcomplex >& p);
	//void mPotentialMatrix(RspaceFunc<double>& Vp, const RspaceFunc<double>& V, const RspaceFunc<double >& p);
	void mPotentialMatrix_ddm(double* Vp, const double* V, const double* p);

	void mSetOccupancy();
    void mSetOccupancyZero();
    void mSetOccupancyTemperature();
	
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
	inline double SQ(double x) { return x * x; };


	//for mpi////////////////////
private:
	// DDMにおいて分散したローカルデータはml_で始まる変数とする//
	const int m_root_id = 0;
	MPI_Comm m_mpi_comm;
	MPI_Comm m_ddm_comm;
	MPI_Comm m_same_ddm_place_comm;
	GridRange m_global_grid;      //common to all spins and k-points
	GridRangeMPI ml_grid;         //common to all spins and k-points
	int m_num_procs_ddm = 1;
	int m_having_spin_kpoint_begin = 0;
	int m_num_having_spin_kpoint = 0;
	int m_mpi_split_color = 0;

	enum class SPIN : int {
		UP, DOWN
	};

#if 1
	
	//kpoinおよびspinごとに保持するデータ構造
	struct WaveSet {		//data for each spin and k-point//
		SoAComplex* l_psi_set = nullptr;		//local grid data split by ddm//
		double* eigen_values = nullptr;
		double* occupancy = nullptr;
		int kpoint_x = 0;
		int kpoint_y = 0;
		int kpoint_z = 0;
		SPIN spin = SPIN::UP;			//1 or -1//
		int kpoint_weight = 1;

		void Allocate(size_t local_size, size_t num_solution, int kx, int ky, int kz, SPIN spin_, int k_weight) {
			l_psi_set = new SoAComplex[num_solution];
			//BLAS利用の為に連続領域としてallocateするのが必須//
			double* buffer = new double[num_solution * local_size * 2]; //2 means complex numebr//
			///l_psi_set[0] = new W[num_solution * local_size];
			for (size_t s = 0; s < num_solution; ++s) {
				l_psi_set[s].re = buffer + (s * 2) * local_size;
				l_psi_set[s].im = buffer + (s * 2 + 1) * local_size;
			}

			eigen_values = new double[num_solution * 2];
			occupancy = eigen_values + num_solution;

			kpoint_x = kx;
			kpoint_y = ky;
			kpoint_z = kz;
			spin = spin_;
			kpoint_weight = k_weight;

			//printf("alloc: k-spin, %d, %d, %d, %d", kpoint_x, kpoint_y, kpoint_z, spin==SPIN::UP?0:1);
			//fflush(stdout);
		}

		~WaveSet() {
			delete[] l_psi_set[0].re;
			delete[] l_psi_set;
			delete[] eigen_values;
		}
	};

	WaveSet* ml_wave_set = nullptr;

	//LOBPCG方を利用するための保存領域
	//WaveSetと対応関係があり、同じ感ず存在する。
	struct KeepLOBPCGSet {
		double* keep_p = nullptr;
		OneComplex* keep_S_matrix = nullptr;
		

		void Allocate(size_t local_size, size_t num_solution) {
			keep_p = new double[num_solution * local_size * 2]; //2 means complex numebr//
			keep_S_matrix = new OneComplex[num_solution * num_solution * 3];
		}

		~KeepLOBPCGSet() {
			delete[] keep_p;
			delete[] keep_S_matrix;
		}
	};

	KeepLOBPCGSet* ml_keep_lobpcg = nullptr;


#else
	template<class W>
	struct WaveSet {		//data for each spin and k-point//
		W** l_psi_set = nullptr;		//local grid data split by ddm//
		double* eigen_values = nullptr;	
		double* occupancy = nullptr;
		int kpoint_x = 0;
		int kpoint_y = 0;
		int kpoint_z = 0;
		SPIN spin = SPIN::UP;			//1 or -1//

		void Allocate(size_t local_size, size_t num_solution, int kx, int ky, int kz, SPIN spin_) {
			l_psi_set = new W*[num_solution];
			l_psi_set[0] = new W[num_solution * local_size];
			for (size_t s = 1; s < num_solution; ++s) {
				l_psi_set[s] = l_psi_set[0] + s * local_size;
			}

			eigen_values = new double[num_solution * 2];
			occupancy = eigen_values + num_solution;

			kpoint_x = kx;
			kpoint_y = ky;
			kpoint_z = kz;
			spin = spin_;
		}

		~WaveSet() {
			delete[] l_psi_set[0];
			delete[] l_psi_set;
			delete[] eigen_values;
		}
	};

	WaveSet<double>* ml_wave_set=nullptr;
#endif

	//double* ml_lobpcg_keep_p = nullptr;

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
	void mHierarchyScatterField(double* local_dest, const double* global_src);
	void mHierarchyGatherField(double* global_dest, const double* local_src);
	
    void mScatterField_HR(double* local_dest, const double* global_src);
    void mGatherField_HR(double* global_dest, const double* local_src);

	//for benchmark/////////////////////////////////////////////
	StopWatch<80, true> watch;
	void mPrintTime();

public:
	void OutputDensity(QUMASUN::OUTPUT_TARGET target, const char* filepath);
	void OutputEigenValue(const char* filepath);
	void OutputEigenVector(const char* filepath);
	void PrintCondition();

	void MoveNuclei(int num_nuclei, const Nucleus* next_nucleis);
    void GetForce(vec3d* forces);
};

#include "qumasun_kpoint.constructor.h"
#include "qumasun_kpoint.core.h"
#include "qumasun_kpoint.density.h"
#include "qumasun_kpoint.energy.h"
#include "qumasun_kpoint.initialize.h"
#include "qumasun_kpoint.lobpcg.h"
#include "qumasun_kpoint.matrix.h"
#include "qumasun_kpoint.occupancy.h"
#include "qumasun_kpoint.potential.h"
#include "qumasun_kpoint.scfloop.h"
#include "qumasun_kpoint.scatter.h"
#include "qumasun_kpoint.output.h"
#include "qumasun_kpoint.force.h"
#include "qumasun_kpoint.move.h"
#include "qumasun_kpoint.print.h"


#endif
