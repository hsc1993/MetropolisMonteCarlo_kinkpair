#ifndef __PARAMETERS_H
#define __PARAMETERS_H

/// length of burgers vector

//#define blength 2.7223 // for tungsten

// #define blength 2.4829 // for iron
// #define blength 2.598 // for V
// #define blength 2.875 // for Nb
// #define blength 2.867 // for Ta
#define blength 2.78 // for gray

/// check parallel
#define epsilon_parallel 1e-1



struct SimulationParameters
{
	int nPhi;
	int nQ;
	int flag_iso_an;
	int n_h;
	double h_initial[1]; 
	double x_initial[1];
	
	int ind_collapse[1][1];							// first index is height, second index is stress						
	double tau[1];
	double tau_p[1];
	double enthalpy_saddle[1];
	double enthalpy_anisotropic_saddle[1];
	double tension_saddle[1];
	double tension_anisotropic_saddle[1];
	double core_saddle[1];
	double U_saddle[1];
	double W_saddle[1];
	double h_saddle[1];
	double a1_saddle[1];
	double a2_saddle[1];
	double w_saddle[1];
	double h_max[1];
	double h_saddle_predict[1];
	double x_equi[1];
	double U_equi[1];
	double x_max[1];
	int ind_saddle[1];
	double x_saddle_predict[1];
	int n_tau;										// Number of interval stresses between stress_min and stress_max
	Eigen::Matrix<double,3,3> T; 					// Trasformation matrix 
	Eigen::Matrix<double,3,3> Tinv; 				// Inverse of Trasformation matrix 
	double tolerance;								// Tolerance for stoping
	double tolerance_W;								// Tolerance for stoping
	double fracOsci;								// The amount of oscilation for monte carlo 
	int N_faliure;									// Maximum faliure attempts allowed for monte carlo
	int N_success;
    double a;						   		    	// Lattice constant (Angstrom units).
	int n_faliure;
	int n_success;									// number of successful attempts that needs to be made
	int iteration;
	Eigen::Matrix<double,3,1> lineDirection;		// [1 1 1]
	Eigen::Matrix<double,3,1> glideDirection;		// [-1 2 -1]
	Eigen::Matrix<double,3,1> lineDirection_ref;	// [0 0 1]
	Eigen::Matrix<double,3,1> glideDirection_ref;	// [1 0 0]

	Tensor<4,3> Cref;				// This is the tensor of elastic constants C. Components are in a reference system parallel to the cubic axes
	double Cref11,Cref12,Cref44;
	Tensor<4,3> C; 				// C in rotated reference system
	double C11,C12,C44;
	Eigen::Matrix<double,3,1> bvector3d;								// The Burgers vector in spatial coordinates.

	// for non-singular theory, L tensor is needed
	Eigen::Matrix<double,3,3> Lref;
	Eigen::Matrix<double,3,3> L;

	
	

	/// parameters for Line Tension Model
	double f1;
	double fa1;
	double f2;
	double fa2;
	double fw;
	int n_segments; 							// Number of segments
	double U_step;
	double W_step;
	double enthalpy_step;
	double x_saddle_predict_ele;
	double h_saddle_predict_ele;
	double h_max_ele;
	double x_max_ele;
	double x_equi_ele;
	double U_equi_ele;


	/// for Elastic Interaction Model
	int lineLength;									// Length of (straight) screw dislocation in units of the Burgers vector. This determines the simulation size.
						
	double h0;
	SymmetricTensor2 externalStress;				// The external stress that drives the dislocation (units: Pa).
	double mu;										// The elastic shear modulus (mu) of the material (units: Pa).
	double pr;										// The elastic Poisson ratio (nu) of the material.
	double a0;										// The core width parameter of the non-singular dislocation theory (units: b).
	double temperature;								// The temperature parameter of the simulation (units: Kelvin).
	double unitt;									// 1 eV/Å3 = 1/unitt GPa

	double tau_min_eva;								// Minimun stress in eva
	double tau_max_eva;								// Maximum stress in eva
	double tau_min_gpa;								// Minimun stress in gva
	double tau_max_gpa;								// Maximum stress in gva

	double w;										// width of kink, measured from the two top points, used for simulation mainly.
	double a1,a2;
	double term1;

	double currentTau;
	double currentTaup;

	double alpha;
	double U0;
};


	
/// see /Users/siconghe/jaime_projects/los_model/eam/elastic/eam_Up_UW  for explanation
struct PotentialFitParams_W
{
	double U0_low1 = 1.4871;    
	double U0_low2 = 0.0230;   
	double U0_low3 = 0.0018;
	double U0_high1 = 0.0046;  
	double U0_high2 = -0.0500;    
	double U0_high3 = 0.0755;    
	double U0_high4 = 0.7676;
	double alpha1 = 0.0769;    
	double alpha2 = 0.2388;    
	double alpha3 = 0.1707;
	int N_land = 5000;
};

struct PotentialFitParams
{
	double U0 = 0.018; // eV/A  
	// double U0 = 0.0112; // eV/A  
	int N_land = 5000;
};



struct EcoreParameters_w
{
	// fitting parameters for p
    double c0pCore = 1.1040;
    double c1pCore = 0.0192;
    double c2pCore = -0.7785;
    double c3pCore = 0.0087;
    double c4pCore = -0.0672;
    double c5pCore = -0.0122;
    double c6pCore = -0.0157;

    // fitting parameters for q
    double d0qCore = 0.7067;
    double d1qCore = -0.1141;
};

struct EcoreParameters
{
	//ecore = c1+c2*cos(2*angle)
    double c1 = 0.51;
    double c2 = -0.29;

};






#endif // __PARAMETERS_H