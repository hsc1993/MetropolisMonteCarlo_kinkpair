#include "../StandardIncludes.h"
#include "../Settings.h"
#include <tuple>
#include "../Parameters.h"

#include "Simulation.h"
#include "../dislocation/Dislocation.h"
#include "../Anisotropic_headers.h"
#include "../SegmentStressCal/SegmentStressCal.h"


using namespace std;



void Simulation::initialize(SimulationParameters& simParams)
{
	simParams.nPhi = 60; // aniso tension: 15:0.00596044   30:0.00550908
	simParams.nQ = 10;    // 5
	simParams.iteration = 500000; //aniso 10000
	simParams.flag_iso_an = 0; // 0 is iso, 1 is anisotropic
	simParams.n_h = 1;
	simParams.n_tau = 1;

	simParams.fracOsci = 0.5;
	simParams.tolerance = 1e-10;
	simParams.tolerance_W = 1e-10;
	simParams.N_faliure = 400;
	simParams.N_success = 3;//3 for aniso
	
	// initializing parameters
	simParams.a = 4e-6*blength;
	simParams.a0 = blength*2/sqrt(3);
	simParams.h0 = simParams.a0*sqrt(6)/3;
	simParams.bvector3d << blength/pow(3,0.5),blength/pow(3,0.5),blength/pow(3,0.5); 				
	simParams.lineLength = 100;				
	simParams.externalStress = SymmetricTensor2(0);				



	simParams.lineDirection << 1/pow(3,0.5),1/pow(3,0.5),1/pow(3,0.5);			// [1 1 1]
	simParams.glideDirection << -1/pow(6,0.5),2/pow(6,0.5),-1/pow(6,0.5);		// [-1 2 -1]
	simParams.lineDirection_ref << 0,0,1;										// [0 0 1]
	simParams.glideDirection_ref << 1,0,0;										// [1 0 0]
	simParams.f1 = 0.4;
	simParams.fa1 = 0.04;
	simParams.f2 = simParams.f1;
	simParams.fa2 = simParams.fa1;
	simParams.fw = 1-simParams.f1-simParams.f2-simParams.fa1-simParams.fa2;
	simParams.n_segments = 100; 

	simParams.Cref11 = simParams.Cref(0,0,0,0);
	simParams.Cref12 = simParams.Cref(0,0,1,1);
	simParams.Cref44 = simParams.Cref(0,1,1,0);

	simParams.Lref << pow((simParams.a),2), 0.0,  0.0,  
					0.0,  pow((simParams.a),2),  0.0,  
					0.0,  0.0,  pow((simParams.a),2);
	
	simParams.mu = simParams.Cref44;				
	simParams.pr = simParams.Cref12/2/(simParams.Cref12+simParams.Cref44);
	// cout << "simParams.mu=" << simParams.mu << endl;

	simParams.mu = 51;
	simParams.pr = 0.37;
	

	
	simParams.temperature = 0;				
	simParams.unitt = 1e-21/(1.6*1e-19);
	simParams.term1 = simParams.mu*pow(blength,2)/(4*M_PI*(1-simParams.pr));	
}

void Simulation::calStressDependents(SimulationParameters& simParams, PotentialFitParams& potFitParams, DislocationNode* pnode)
{
		double U_land[potFitParams.N_land];
		double x_landscape[potFitParams.N_land];
		double U_land_gradiant[potFitParams.N_land];
		double UW_land[potFitParams.N_land];

		simParams.U0 = potFitParams.U0;

		for (int i = 0; i < potFitParams.N_land; i++){
			x_landscape[i] = simParams.h0/potFitParams.N_land*i;
			U_land[i] = pnode->Up_vs_x(simParams,x_landscape[i]);
			U_land_gradiant[i] = M_PI*simParams.U0/simParams.h0*sin(2*M_PI*x_landscape[i]/simParams.h0);
			UW_land[i] = U_land[i]-simParams.currentTau*blength*x_landscape[i];
		}		
		
		double smallest = abs(U_land_gradiant[0]-simParams.currentTau*blength);
		int ind_equi = 0;
		for (int i = 0; i < potFitParams.N_land/4; i++){
			double temp = abs(U_land_gradiant[i]-simParams.currentTau*blength);
			if(temp < smallest){
				smallest = temp;
				ind_equi = i;
				}
		}
		simParams.x_equi_ele = x_landscape[ind_equi];
		simParams.U_equi_ele = pnode->Up_vs_x(simParams, simParams.x_equi_ele);
		cout << "simParams.x_equi_ele  = " << simParams.x_equi_ele  << endl;
		

		double U_equi = U_land[ind_equi];
		double biggest = UW_land[ind_equi];
		int ind_max = 0;
		double temp_max;
		for (int i = ind_equi; i < potFitParams.N_land; i++){
			temp_max = UW_land[i];
			if(temp_max > biggest){
				biggest = temp_max;
				ind_max = i;
				}
		}

		

		double temp_saddle = 1;
		int ind_saddle;
		for (int i = ind_max; i < potFitParams.N_land; i++){
			if (temp_saddle > abs(UW_land[i]-UW_land[ind_equi]))
			{
				temp_saddle = abs(UW_land[i]-UW_land[ind_equi]);
				simParams.x_saddle_predict_ele = x_landscape[i];
				ind_saddle = i;
			}
		}
		
		simParams.x_max_ele = x_landscape[ind_max];
		simParams.h_max_ele = simParams.x_max_ele-simParams.x_equi_ele;
		simParams.h_saddle_predict_ele = simParams.x_saddle_predict_ele-simParams.x_equi_ele;
		
		cout << "simParams.h_saddle_predict_ele = " << simParams.h_saddle_predict_ele << endl;
		cout << "simParams.h0 = " << simParams.h0 << endl;
		for (int k = 0; k < simParams.n_h; k++){
			simParams.h_initial[k] = simParams.h0-(k*1.0)/simParams.n_h*simParams.h0;
			simParams.h_initial[k] = simParams.h_saddle_predict_ele;
			simParams.x_initial[k] = simParams.h_initial[k]+simParams.x_equi_ele;

			cout << "simParams.h_initial[k] = " << simParams.h_initial[k] << endl;
			cout << "simParams.x_initial[k] = " << simParams.x_initial[k] << endl;
		}	
		cout << "simParams.x_equi_ele = " << simParams.x_equi_ele << endl;
		
}



void Simulation::calStressDependents_W(SimulationParameters& simParams, PotentialFitParams_W& potFitParams)
{
		double U_land[potFitParams.N_land];
		double x_landscape[potFitParams.N_land];
		double U_land_gradiant[potFitParams.N_land];
		double UW_land[potFitParams.N_land];
		if (simParams.currentTaup < 0.8){
			simParams.U0 = potFitParams.U0_low3*pow(simParams.currentTaup,potFitParams.U0_low1)+potFitParams.U0_low2;
		}else{
			simParams.U0 = potFitParams.U0_low3*pow(simParams.currentTaup,potFitParams.U0_low1)+potFitParams.U0_low2;
			//simParams.U0 = potFitParams.U0_high3*pow((simParams.currentTaup-potFitParams.U0_high4),potFitParams.U0_high1)+potFitParams.U0_high2;
		}
		simParams.alpha = potFitParams.alpha1*(simParams.currentTaup-potFitParams.alpha2)+potFitParams.alpha3;

		for (int i = 0; i < potFitParams.N_land; i++){
			x_landscape[i] = simParams.h0/potFitParams.N_land*i;
			U_land[i] = simParams.U0/(2-2*simParams.alpha)*(1-cos(2*M_PI*x_landscape[i]/simParams.h0)-simParams.alpha/2*pow(1-cos(2*M_PI*x_landscape[i]/simParams.h0),2));
			U_land_gradiant[i] = simParams.U0/(1-simParams.alpha)*M_PI/simParams.h0*sin(2*M_PI*x_landscape[i]/simParams.h0)*(1-simParams.alpha*(1-cos(2*M_PI*x_landscape[i]/simParams.h0)));
			UW_land[i] = U_land[i]-simParams.currentTau*blength*x_landscape[i];
			
		}		
		
		double smallest = abs(U_land_gradiant[0]-simParams.currentTau*blength);
		int ind_equi = 0;
		for (int i = 0; i < potFitParams.N_land/4; i++){
			double temp = abs(U_land_gradiant[i]-simParams.currentTau*blength);
			if(temp < smallest){
				smallest = temp;
				ind_equi = i;
				}
		}

		double U_equi = U_land[ind_equi];

		double biggest = UW_land[ind_equi];
		int ind_max = 0;
		double temp_max;
		for (int i = ind_equi; i < potFitParams.N_land; i++){
			temp_max = UW_land[i];
			if(temp_max > biggest){
				biggest = temp_max;
				ind_max = i;
				}
		}

		simParams.x_equi_ele = x_landscape[ind_equi];

		double temp_saddle = 1;
		int ind_saddle;
		for (int i = ind_max; i < potFitParams.N_land; i++){
			if (temp_saddle > abs(UW_land[i]-UW_land[ind_equi]))
			{
				temp_saddle = abs(UW_land[i]-UW_land[ind_equi]);
				simParams.x_saddle_predict_ele = x_landscape[i];
				ind_saddle = i;
			}
		}
		
		simParams.x_max_ele = x_landscape[ind_max];
		simParams.h_max_ele = simParams.x_max_ele-simParams.x_equi_ele;
		simParams.h_saddle_predict_ele = simParams.x_saddle_predict_ele-simParams.x_equi_ele;

		for (int k = 0; k < simParams.n_h; k++){
			simParams.h_initial[k] = simParams.h0-(k*1.0)/simParams.n_h*simParams.h0;
		}	
}



void Simulation::updateStressDependents(SimulationParameters& simParams, int indTau)
{
	simParams.tau[indTau] = simParams.currentTau;
	simParams.tau_p[indTau] = simParams.currentTaup;
	simParams.x_equi[indTau] = simParams.x_equi_ele;
	simParams.x_max[indTau] = simParams.x_max_ele;
	simParams.h_max[indTau] = simParams.h_max_ele;
	simParams.h_saddle_predict[indTau] = simParams.h_saddle_predict_ele;
	simParams.x_saddle_predict[indTau] = simParams.h_saddle_predict[indTau]+simParams.x_equi[indTau];
}



/*******************************************************************************
* creating trasformation matrix T  x[1 0 0],y[0 1 0],z[0 0 1]
********************************************************************************/
Eigen::Matrix<double,3,3> Simulation::createTransformationMatrix()
{
    Eigen::Matrix<double,3,3> T;
    Eigen::Matrix<double,3,1> position;
    
    Eigen::Matrix<double,3,1> x_coor;
    Eigen::Matrix<double,3,1> y_coor;
    Eigen::Matrix<double,3,1> z_coor;

    x_coor << 1, 0, 0;
    y_coor << 0, 1, 0;
    z_coor << 0, 0, 1;

    Eigen::Matrix<double,3,1> x_coor_trans;
    Eigen::Matrix<double,3,1> y_coor_trans;
    Eigen::Matrix<double,3,1> z_coor_trans;

    x_coor_trans << -1/pow(6,0.5), 2/pow(6,0.5), -1/pow(6,0.5);
    y_coor_trans << -1/pow(2,0.5), 0, 1/pow(2,0.5);
    z_coor_trans << 1/pow(3,0.5), 1/pow(3,0.5), 1/pow(3,0.5);

    T(0,0) = x_coor_trans.transpose()*x_coor;
    T(0,1) = x_coor_trans.transpose()*y_coor;
    T(0,2) = x_coor_trans.transpose()*z_coor;

    T(1,0) = y_coor_trans.transpose()*x_coor;
    T(1,1) = y_coor_trans.transpose()*y_coor;
    T(1,2) = y_coor_trans.transpose()*z_coor;

    T(2,0) = z_coor_trans.transpose()*x_coor;
    T(2,1) = z_coor_trans.transpose()*y_coor;
    T(2,2) = z_coor_trans.transpose()*z_coor;

    return T;    
}

void Simulation::rotate_tensor_4th(Tensor<4,3> Cref,Eigen::Matrix<double,3,3> T, SimulationParameters& simParams)
{
	Tensor<4,3> Cp;

	for (int i = 0; i < 3; i++){
		for (int j = 0; j < 3; j++){
			for (int k = 0; k < 3; k++){
				for (int l = 0; l < 3; l++){
					for (int m = 0; m < 3; m++){
						for (int n = 0; n < 3; n++){
							for (int s = 0; s < 3; s++){
								for (int t = 0; t < 3; t++){
									Cp(i,j,k,l) += T(i,m)*T(j,n)*T(k,s)*T(l,t)*Cref(m,n,s,t);
								}
							}
						}
					}
				}
			}
		}
	}
	
	simParams.C = Cp;
}

void Simulation::rotate_tensor_3th(Eigen::Matrix<double,3,3> Lref, Eigen::Matrix<double,3,3> T, SimulationParameters& simParams)
{
	Eigen::Matrix<double,3,3> Lp;

	for (int i = 0; i < 3; i++){
		for (int j = 0; j < 3; j++){
			for (int m = 0; m < 3; m++){
				for (int n = 0; n < 3; n++){
					Lp(i,j) += T(i,m)*T(j,n)*Lref(m,n);
				}
			}
		}
	}
	simParams.L = Lp;
}

















