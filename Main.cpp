#include "StandardIncludes.h"
#include "Settings.h"
#include "simulation/Simulation.h"
#include <tuple>
#include "Anisotropic_headers.h"
#include "SegmentStressCal/SegmentStressCal.h"
#include<cstdlib>
#include<ctime>
//#include "../context/Context.h"

// to compile: g++ Main.cpp simulation/Simulation.cpp dislocation/Dislocation.cpp SegmentStressCal/SegmentStressCal.cpp -o main -std=c++11

using namespace std;


int main(int argc, char* argv[])
{	
	srand(time(0));
	Simulation simulation;
	SimulationParameters simParams;
	PotentialFitParams potParams;
	EcoreParameters ecoreParams_w;
	EcoreParameters ecoreParams;


	simParams.tau_min_gpa = 0.8;
	simParams.tau_max_gpa = simParams.tau_min_gpa;

	ofstream stressDependents("stressDependents.txt");
	ofstream initial_config("initial_config.txt");
	ofstream saddle_config("saddle_config.txt");
	ofstream initial_config_ref("initial_config_ref.txt");
	ofstream saddle_config_ref("saddle_config_ref.txt");

	/*
	ofstream all_config("all_config.txt");
	ofstream ind_collapse_txt("ind_collapse.txt");
	ofstream final_report("final_report.txt");
	ofstream all_data("all_data.txt");
	*/

/*
	ofstream Wself_an_file("Wself_an.txt");
	ofstream Wself_iso_file("Wself_iso.txt");
	ofstream Wself_ratio_file("Wself_ratio.txt");
	ofstream r1coor_file("r1coor.txt");

	ofstream Wself_file("Wself.txt");
	ofstream Wself0_file("Wself0.txt");
	ofstream lt_file("lt.txt");
	ofstream core_file("core.txt");
*/
	ofstream enthalpy_anisotropic_evolution("enthalpy_anisotropic_evolution.txt");
	ofstream enthalpy_evolution("enthalpy_evolution.txt");
	ofstream tension_anisotropic_evolution("tension_anisotropic_evolution.txt");
	ofstream tension_evolution("tension_evolution.txt");
	ofstream core_evolution("core_evolution.txt");
	ofstream U_evolution("U_evolution.txt");
	ofstream W_evolution("W_evolution.txt");

	

	/****************  Anisotropic set up ****************/

	//std::string potentialFolder="Fe_bcc/BO_Mueller_2007";
	//double bs = 2.477001910395397; // Burgers vector magnitude
	//std::string potentialFolder="Fe_bcc/MEAM_Lee_2001";
	std::string potentialFolder="Fe_bcc/Jaime_2020";
	//std::string potentialFolder="W_bcc";
	
	readCtensor("./ElasticConstants/"+potentialFolder, simParams.Cref);
	simulation.initialize(simParams);	
	cout << "simParams.mu = " << simParams.mu << endl;
	cout << "simParams.pr = " << simParams.pr << endl;
	cout << "Cref [GPa]= " << endl 
	<< "Cref11 = " << simParams.Cref11 << endl
	<< "Cref12 = " << simParams.Cref12 << endl
	<< "refC44 = " << simParams.Cref44 << endl << endl;
	cout << "Lref = [A^2]" << endl 
	<< simParams.Lref << endl << endl;

	simParams.T = simulation.createTransformationMatrix();
	simParams.Tinv = simParams.T.inverse();
	//simulation.rotate_tensor_4th(simParams.Cref, simParams.T, simParams);
/*
	simParams.C11 = simParams.C(0,0,0,0);
	simParams.C12 = simParams.C(0,0,1,1);
	simParams.C44 = simParams.C(0,1,1,0);

	cout << "C [GPa]= " << endl 
	<< "C11 = " << simParams.C11 << endl
	<< "C12 = " << simParams.C12 << endl
	<< "C44 = " << simParams.C44 << endl << endl;

	
	simulation.rotate_tensor_3th(simParams.Lref, simParams.T, simParams);
	cout << "L = [A^2]" << endl 
	<< simParams.L << endl << endl;
*/

	// Alan 
	// screw: glide = [-1 2 -1], plane normal = [-1 0 1], line direction = [1 1 1]. burgers = [1 1 1]

	// giacomo 
	// edge:  glide = [1 1 1], plane normal = [-1 1 0], line direction[-1 -1 2] 	burgers = [1 1 1] 
	DislocationNode nodeSaddle[simParams.n_segments][simParams.n_tau];
	DislocationNode* pnodeSaddle[simParams.n_segments][simParams.n_tau];

	DislocationLoop disLoop[simParams.n_h][simParams.n_tau];
	DislocationLoop* pdisLoop[simParams.n_h][simParams.n_tau];
	DislocationNode node[simParams.n_segments][simParams.n_h][simParams.n_tau];
	DislocationNode* pnode[simParams.n_segments][simParams.n_h][simParams.n_tau];


	/****************  Anisotropic set up ****************/
	
	/// for each given stress
	simParams.tau_min_eva = simParams.tau_min_gpa*simParams.unitt;   		    //   1 eV/Å3 = 160.2176621 GPa.  unitt = 0.0063, 1GPa = 1 unitt * 1 eV/A^3 
	simParams.tau_max_eva = simParams.tau_max_gpa*simParams.unitt;				//   1 eV/Å3 = 1/unitt GPa
	for (int indTau = 0; indTau < simParams.n_tau; indTau++){  
		// calculate current stress
		if(simParams.n_tau!=1){
			simParams.currentTau = simParams.tau_min_eva+indTau*(simParams.tau_max_eva-simParams.tau_min_eva)/(simParams.n_tau-1);
		}else{
			simParams.currentTau = simParams.tau_min_eva;
		}
		simParams.currentTaup = simParams.currentTau/simParams.unitt;    
		cout << endl << "currentTaup = " << simParams.currentTaup << "GPa" << endl;
		cout << "currentTau = " << simParams.currentTau << "eV/A3" << endl;
		cout << "currentTau*b = " << simParams.currentTau*blength << "eV/A2" << endl << endl;

		simulation.calStressDependents(simParams,potParams,pnode[0][0][0]);
		simulation.updateStressDependents(simParams,indTau);
		stressDependents << simParams.x_equi_ele << " " << simParams.x_saddle_predict_ele << endl;
		simParams.enthalpy_saddle[indTau] = 0;

		for (int indh = 0; indh < simParams.n_h; indh++){
			for (int inode = 0; inode < simParams.n_segments; inode++){
				pnode[inode][indh][indTau] = &node[inode][indh][indTau];
				pnode[inode][indh][indTau]->set_inode(inode);
				pnode[inode][indh][indTau]->set_indh(indh);
				pnode[inode][indh][indTau]->set_indTau(indTau);
				pnode[inode][indh][indTau]->initialize(simParams);
				pnode[inode][indh][indTau]->setPositionTransformed(pnode[inode][indh][indTau]->transform(pnode[inode][indh][indTau]->position(),simParams.T)); 
				initial_config << pnode[inode][indh][indTau]->posXTransformed() << " "
								<< pnode[inode][indh][indTau]->posYTransformed() << " "
								<< pnode[inode][indh][indTau]->posZTransformed() << " " 
								<< endl;
				initial_config_ref << pnode[inode][indh][indTau]->posX() << " "
								<< pnode[inode][indh][indTau]->posY() << " "
								<< pnode[inode][indh][indTau]->posZ() << " " 
								<< endl;
				}
			for (int inode = 0; inode < simParams.n_segments-1; inode++){
					pnode[inode][indh][indTau]->setNext(pnode[inode+1][indh][indTau]);
					pnode[inode][indh][indTau]->setxF(0.5*(pnode[inode][indh][indTau]->positionTransformed()
													 +pnode[inode][indh][indTau]->nextNode()->positionTransformed())
													 );
				}
				pnode[simParams.n_segments-1][indh][indTau]->setxF(pnode[simParams.n_segments-1][indh][indTau]->positionTransformed());
				pnode[simParams.n_segments-1][indh][indTau]->setNext(NULL);
		}

		for (int indh = 0; indh < simParams.n_h; indh++){
			pdisLoop[indh][indTau] = &disLoop[indh][indTau];
			pdisLoop[indh][indTau]->set_firstNode(pnode[0][indh][indTau]);
		}


		/****** Calculating Screw Segment Wself *******/
		if (simParams.flag_iso_an == 1){
		Eigen::Matrix<double,3,1> originPoint, screwBurgers, fieldPoint;
		originPoint << 0,0,0;
		screwBurgers = simParams.bvector3d;
		fieldPoint = simParams.bvector3d/2;

		SegmentStressClassical screwSegEnergy(simParams.Cref,simParams.Lref,simParams.nPhi,simParams.nQ,simParams.bvector3d,originPoint,screwBurgers,1); // the segment stress evaluator object
		SegmentStressClassical* pscrewSegEnergy;
		pscrewSegEnergy = &screwSegEnergy;
		for (int indh = 0; indh < simParams.n_h; indh++){
			for (int inode = 0; inode < simParams.n_segments; inode++){
				pnode[inode][indh][indTau]->ScrewSegWself(pscrewSegEnergy,fieldPoint);
			}
		}
		}

		for (int ii = 0;  ii < simParams.n_segments; ii++){
			pnodeSaddle[ii][indTau] = &nodeSaddle[ii][indTau];
		}
		simParams.ind_saddle[indTau] = 0;

		
		int inode;
		/// for each tentative saddle height
		for (int indh = 0; indh < simParams.n_h; indh++){
			/// as stress increases, saddle point height decreases, so if one height leads to collapse for lower stress, 
			/// higher configuration than "that height" should be ignored for higher stress.
			if (indTau > 0 && indh < simParams.ind_saddle[indTau-1]){
				//simParams.ind_collapse[indh][indTau] = 3;
				//continue;
			}
			if (indh > simParams.ind_saddle[indTau]){
				//simParams.ind_collapse[indh][indTau] = 3;
				//continue;
			}
			
			/// calculating each energy term
			for (inode = 0; inode < simParams.n_segments-1; inode++){
				pnode[inode][indh][indTau]->calTheta(simParams);
				pnode[inode][indh][indTau]->calGlideLength(simParams);
			}
			for (inode = 0; inode < simParams.n_segments-1; inode++){
				if (simParams.flag_iso_an == 0){
					pnode[inode][indh][indTau]->calNodeEnergy(ecoreParams, simParams);
					// cout << "pnode[inode][indh][indTau] = " << pnode[inode][indh][indTau]->enthalpy() << endl;
				}

				if (simParams.flag_iso_an == 1){
					SegmentStressClassical ss(simParams.Cref,simParams.Lref,simParams.nPhi,simParams.nQ,simParams.bvector3d,
					pnode[inode][indh][indTau]->positionTransformed(),pnode[inode][indh][indTau]->nextNode()->positionTransformed(),1);
					SegmentStressClassical* pss;
					pss = &ss;
					pnode[inode][indh][indTau]->calNodeEnergy_anisotropic(pss, ecoreParams, simParams);
				}

			}
			if (simParams.flag_iso_an == 0){
				pdisLoop[indh][indTau]->calLoopEnergy(simParams);
			}
			if (simParams.flag_iso_an == 1){
				pdisLoop[indh][indTau]->calLoopEnergy_anisotropic(simParams);
			}

			if (simParams.flag_iso_an == 0){
			cout << "initial enthalpy = " << pdisLoop[indh][indTau]->enthalpy() << endl
			<< "initial tension() = " << pdisLoop[indh][indTau]->tension() << endl
			<< "initial core = " << pdisLoop[indh][indTau]->core() << endl
			<< "initial U = " << pdisLoop[indh][indTau]->U() << endl
			<< "initial W = " << pdisLoop[indh][indTau]->W() << endl << endl;
			}

			if (simParams.flag_iso_an == 1){
			cout << "initial enthalpy anisotropic = " << pdisLoop[indh][indTau]->enthalpy_anisotropic() << endl
			<< "initial tension() anisotropic = " << pdisLoop[indh][indTau]->tension_anisotropic() << endl
			<< "initial core = " << pdisLoop[indh][indTau]->core() << endl
			<< "initial U = " << pdisLoop[indh][indTau]->U() << endl
			<< "initial W = " << pdisLoop[indh][indTau]->W() << endl << endl;
			}


			/// creating two trial nodes for calculation of enthalpy_trial
			DislocationNode nodeTrial1;
			DislocationNode* pnodeTrial1;
			pnodeTrial1 = &nodeTrial1;
			pnodeTrial1->set_screwSegWself(pnode[0][indh][indTau]->screwSegWself());

			DislocationNode nodeTrial2;
			DislocationNode* pnodeTrial2;
			pnodeTrial2 = &nodeTrial2;
			pnodeTrial2->set_screwSegWself(pnode[0][indh][indTau]->screwSegWself());
			double enthalpy_trial,enthalpy_compare;

			/// Monte Carlo
			simParams.n_success = 0;
			simParams.n_faliure = 0;
			while(0 < 1)
			// for (int itemp = 0; itemp<simParams.iteration; itemp++)
			{
				// cout << "itemp = " << itemp << endl;
				enthalpy_trial = 0;
				enthalpy_compare = 0;
				int ind_trial = rand() % (simParams.n_segments-2) + 1;
				double random_pre = rand() % simParams.n_segments + 1;
				double random = random_pre/(simParams.n_segments*1.0);
				double random_fluc = (0.5-random)*simParams.fracOsci;

				pnodeTrial2->setPositionTransformed(pnode[ind_trial][indh][indTau]->positionTransformed()+0.2*random_fluc*simParams.glideDirection);//for isotropic
				//pnodeTrial2->setPositionTransformed(pnode[ind_trial][indh][indTau]->positionTransformed()+random_fluc*simParams.glideDirection);//for anisotropic
				pnodeTrial1->setPositionTransformed(pnode[ind_trial-1][indh][indTau]->positionTransformed());
				pnodeTrial2->setNext(pnode[ind_trial+1][indh][indTau]);
				pnodeTrial1->setNext(pnodeTrial2);
				

				pnodeTrial2->calTheta(simParams);
				pnodeTrial2->calGlideLength(simParams);
				pnodeTrial2->setxF(0.5*(pnodeTrial2->positionTransformed()+pnodeTrial2->nextNode()->positionTransformed()));
				pnodeTrial2->set_x_equi(pnode[ind_trial][indh][indTau]->x_equi());
				

				pnodeTrial1->calTheta(simParams);
				pnodeTrial1->calGlideLength(simParams);
				pnodeTrial1->setxF(0.5*(pnodeTrial1->positionTransformed()+pnodeTrial1->nextNode()->positionTransformed()));
				pnodeTrial1->set_x_equi(pnode[ind_trial-1][indh][indTau]->x_equi());

				pnode[ind_trial-1][indh][indTau]->calGlideLength(simParams);
				pnode[ind_trial][indh][indTau]->calGlideLength(simParams);
					
				/// calculate enthalpy corresponding to a step
				if (simParams.flag_iso_an == 0){
					pdisLoop[indh][indTau]->calEnthalpyStep(ecoreParams, simParams, pnodeTrial1, pnodeTrial2,
						pnode[ind_trial-1][indh][indTau], pnode[ind_trial][indh][indTau]);
					enthalpy_trial = pnodeTrial1->enthalpy() + pnodeTrial2->enthalpy();
					enthalpy_compare = pnode[ind_trial-1][indh][indTau]->enthalpy() + pnode[ind_trial][indh][indTau]->enthalpy();
					// cout << "enthalpy_trial = " << enthalpy_trial << endl;
					// cout << "enthalpy_compare = " << enthalpy_compare << endl;
				}
				
				if (simParams.flag_iso_an == 1){
					SegmentStressClassical ss_step_trial1(simParams.Cref,simParams.Lref,simParams.nPhi,simParams.nQ,simParams.bvector3d,
					pnodeTrial1->positionTransformed(),pnodeTrial1->nextNode()->positionTransformed(),1);
					SegmentStressClassical* pss_step_trial1;
					pss_step_trial1 = &ss_step_trial1;
					pnodeTrial1->calNodeEnergy_anisotropic(pss_step_trial1, ecoreParams, simParams);

					SegmentStressClassical ss_step_trial2(simParams.Cref,simParams.Lref,simParams.nPhi,simParams.nQ,simParams.bvector3d,
					pnodeTrial2->positionTransformed(),pnodeTrial2->nextNode()->positionTransformed(),1);
					SegmentStressClassical* pss_step_trial2;
					pss_step_trial2 = &ss_step_trial2;
					pnodeTrial2->calNodeEnergy_anisotropic(pss_step_trial2, ecoreParams, simParams);
					
/*
					SegmentStressClassical ss_step_1(simParams.Cref,simParams.Lref,simParams.nPhi,simParams.nQ,simParams.bvector3d,
					pnode[ind_trial-1][indh][indTau]->positionTransformed(),pnode[ind_trial-1][indh][indTau]->nextNode()->positionTransformed(),1);
					SegmentStressClassical* pss_step_1;
					pss_step_1 = &ss_step_1;
					pnode[ind_trial-1][indh][indTau]->calNodeEnergy_anisotropic(pss_step_1, ecoreParams, simParams);

					SegmentStressClassical ss_step_2(simParams.Cref,simParams.Lref,simParams.nPhi,simParams.nQ,simParams.bvector3d,
					pnode[ind_trial][indh][indTau]->positionTransformed(),pnode[ind_trial][indh][indTau]->nextNode()->positionTransformed(),1);
					SegmentStressClassical* pss_step_2;
					pss_step_2 = &ss_step_2;
					pnode[ind_trial][indh][indTau]->calNodeEnergy_anisotropic(pss_step_2, ecoreParams, simParams);
*/

					enthalpy_trial = pnodeTrial1->enthalpy_anisotropic() + pnodeTrial2->enthalpy_anisotropic();
					enthalpy_compare = pnode[ind_trial-1][indh][indTau]->enthalpy_anisotropic() + pnode[ind_trial][indh][indTau]->enthalpy_anisotropic();


					if (pdisLoop[indh][indTau]->enthalpy_anisotropic() < 0){
						break;
					}
				}

				/// final convergence, exit current loop and update all values
				if (simParams.n_success > simParams.N_success){
					simParams.ind_collapse[indh][indTau] = 0;
					break;
				}

				/// too much faliure means that tolerance is too small for current scheme
				if (simParams.n_faliure > simParams.N_faliure){
					simParams.n_success += 1;
					simParams.n_faliure = 0;
					continue;
				}

				/// converged value
				if ((abs(enthalpy_compare-enthalpy_trial) < abs(simParams.tolerance*enthalpy_compare) && enthalpy_trial>enthalpy_compare)){
					cout << "n_success = " << simParams.n_success << endl;	
					simParams.n_success += 1;
					simParams.n_faliure = 0;
					continue;
				}
				

				/// accepted lower enthalpy
			 	if (enthalpy_trial < enthalpy_compare){	
					 //cout << "difference = " << enthalpy_trial-enthalpy_compare << endl;
					//  cout << "pdisLoop[indh][indTau]->enthalpy_anisotropic() = " << pdisLoop[indh][indTau]->enthalpy_anisotropic() << endl;


					pnode[ind_trial][indh][indTau]->setPositionTransformed(pnodeTrial2->positionTransformed());

					if (simParams.flag_iso_an == 0){
						for (inode = 0; inode < simParams.n_segments-1; inode++){
							pnode[inode][indh][indTau]->calTheta(simParams);
							pnode[inode][indh][indTau]->calGlideLength(simParams);
						}
						for (inode = 0; inode < simParams.n_segments-1; inode++){
							pnode[inode][indh][indTau]->calNodeEnergy(ecoreParams, simParams);
						}
						pdisLoop[indh][indTau]->calLoopEnergy(simParams);
						
					}

					if (simParams.flag_iso_an == 1){
					//for (inode = ind_trial-1; inode < ind_trial+1; inode++){
					for (inode = 0; inode < simParams.n_segments-1; inode++){
						pnode[inode][indh][indTau]->calTheta(simParams);
						pnode[inode][indh][indTau]->calGlideLength(simParams);
						pnode[inode][indh][indTau]->setxF(0.5*(pnode[inode][indh][indTau]->positionTransformed()
													 +pnode[inode+1][indh][indTau]->positionTransformed())
													 );
					}
					//for (inode = ind_trial-1; inode < ind_trial+1; inode++){
					for (inode = 0; inode < simParams.n_segments-1; inode++){
						SegmentStressClassical ss(simParams.Cref,simParams.Lref,simParams.nPhi,simParams.nQ,simParams.bvector3d,
						pnode[inode][indh][indTau]->positionTransformed(),pnode[inode][indh][indTau]->nextNode()->positionTransformed(),1);
						SegmentStressClassical* pss;
						pss = &ss;
						pnode[inode][indh][indTau]->calNodeEnergy_anisotropic(pss, ecoreParams, simParams);
					}	
					pdisLoop[indh][indTau]->calLoopEnergy_anisotropic(simParams);
					}
					
					enthalpy_anisotropic_evolution << pdisLoop[indh][indTau]->enthalpy_anisotropic() << endl;
					enthalpy_evolution << pdisLoop[indh][indTau]->enthalpy() << endl;
					tension_anisotropic_evolution << pdisLoop[indh][indTau]->tension_anisotropic() << endl;
					tension_evolution << pdisLoop[indh][indTau]->tension() << endl;
					core_evolution << pdisLoop[indh][indTau]->core() << endl;
					U_evolution << pdisLoop[indh][indTau]->U() << endl;
					W_evolution << pdisLoop[indh][indTau]->W() << endl;

					simParams.n_success = 0;
					simParams.n_faliure = 0;
					continue;
				}

				/// deneyed higher enthalpy
				if (enthalpy_trial > enthalpy_compare){
					simParams.n_faliure += 1;
					continue;
				}


				/**************************************************
				 *  check collapse 
				 *************************************************/      
				/// expand up or down
				/*
				if (abs(pnode[simParams.n_segments/2][indh][indTau]->posX()-simParams.x_initial[indh])>5*simParams.h0/simParams.n_h){
					simParams.ind_collapse[indh][indTau] = 1;
					break;
				}


				/// expand to left
				int ind_left = ceil(simParams.n_segments*(simParams.f1+simParams.fa1))-simParams.n_segments/10.0;
				if (abs(pnode[ind_left][indh][indTau]->posX()-simParams.x_initial[indh])<simParams.h0/simParams.n_h){
					simParams.ind_collapse[indh][indTau] = 2;
					break;
				}
				
				///expand to right
				int ind_right = ceil(simParams.n_segments*(simParams.f1+simParams.fa1+simParams.fw))+simParams.n_segments/10.0;
				if (abs(pnode[ind_right][indh][indTau]->posX()-simParams.x_initial[indh])<simParams.h0/simParams.n_h){
					simParams.ind_collapse[indh][indTau] = 2;
					break;
				}
				*/
				
			}// end of monte carlo loop
			cout << "height No." << indh << " = " << simParams.h_initial[indh] 
			<< "   enthalpy = " << pdisLoop[indh][indTau]->enthalpy()
			<< "   ind_collapse[" << indTau << "][" << indh << "] = " << simParams.ind_collapse[indh][indTau] 
			<< endl;

			for (int ii = 0; ii < simParams.n_segments; ii++){
				//all_config << pnode[ii][indh][indTau]->posX() << endl;
			}
		}// end of looping for height

	} // end of looping tau
		/// finding index of saddle point config and calculate all values
		
		for (int jj = 0; jj < simParams.n_tau; jj++){
			double temp_saddle_enthalpy = 0;
			for (int hh = 0; hh < simParams.n_h; hh++){
				if (simParams.ind_collapse[hh][jj] == 0 && (pdisLoop[hh][jj]->enthalpy() > temp_saddle_enthalpy)){
					temp_saddle_enthalpy = pdisLoop[hh][jj]->enthalpy();
					simParams.ind_saddle[jj] = hh;
				}
			}
		}

		for (int jj = 0; jj < simParams.n_tau; jj++){
			simParams.enthalpy_saddle[jj] = pdisLoop[simParams.ind_saddle[jj]][jj]->enthalpy();
			simParams.enthalpy_anisotropic_saddle[jj] = pdisLoop[simParams.ind_saddle[jj]][jj]->enthalpy_anisotropic();
			simParams.tension_saddle[jj] = pdisLoop[simParams.ind_saddle[jj]][jj]->tension();
			simParams.tension_anisotropic_saddle[jj] = pdisLoop[simParams.ind_saddle[jj]][jj]->tension_anisotropic();
			simParams.core_saddle[jj] = pdisLoop[simParams.ind_saddle[jj]][jj]->core();
			simParams.U_saddle[jj] = pdisLoop[simParams.ind_saddle[jj]][jj]->U();
			simParams.W_saddle[jj] = pdisLoop[simParams.ind_saddle[jj]][jj]->W();
			for (int ii = 0; ii < simParams.n_segments; ii++){
				pnodeSaddle[ii][jj]->setPositionTransformed(pnode[ii][simParams.ind_saddle[jj]][jj]->positionTransformed());
				pnodeSaddle[ii][jj]->setPosition(pnodeSaddle[ii][jj]->transform(pnodeSaddle[ii][jj]->positionTransformed(),simParams.Tinv));
				pnodeSaddle[ii][jj]->set_tension(pnode[ii][simParams.ind_saddle[jj]][jj]->tension());
				pnodeSaddle[ii][jj]->set_tension_anisotropic(pnode[ii][simParams.ind_saddle[jj]][jj]->tension_anisotropic());
				pnodeSaddle[ii][jj]->set_core(pnode[ii][simParams.ind_saddle[jj]][jj]->core());
				pnodeSaddle[ii][jj]->set_U(pnode[ii][simParams.ind_saddle[jj]][jj]->U());
				pnodeSaddle[ii][jj]->set_W(pnode[ii][simParams.ind_saddle[jj]][jj]->W());
			}
		}





		for (int jj = 0; jj < simParams.n_tau; jj++){
			for (int ii = 0; ii < simParams.n_segments; ii++){
				saddle_config << pnodeSaddle[ii][jj]->posXTransformed() << " "
							<< pnodeSaddle[ii][jj]->posYTransformed() << " "
							<< pnodeSaddle[ii][jj]->posZTransformed() << " " 
							<< endl;
				saddle_config_ref << pnodeSaddle[ii][jj]->posX() << " "
							<< pnodeSaddle[ii][jj]->posY() << " "
							<< pnodeSaddle[ii][jj]->posZ() << " " 
							<< endl;
			}
			cout << "stress = " << simParams.tau_p[jj] << endl
			<< "tension_saddle = " << simParams.tension_saddle[jj] << endl
			<< "tension_anisotropic_saddle = " << simParams.tension_anisotropic_saddle[jj] << endl
			<< "core_saddle = " << simParams.core_saddle[jj] << endl
			<< "U_saddle = " << simParams.U_saddle[jj] << endl
			<< "W_saddle = " << simParams.W_saddle[jj] << endl
			<< "height = " << simParams.h_initial[simParams.ind_saddle[jj]] << endl
			<< "enthalpy_saddle = " << simParams.enthalpy_saddle[jj] << endl
			<< "enthalpy_anisotropic_saddle = " << simParams.enthalpy_anisotropic_saddle[jj] << endl
			<< "ind_saddle = " << simParams.ind_saddle[jj] << endl
			<< endl;
		}

		stressDependents.close();
		saddle_config.close();
		initial_config.close();
		saddle_config_ref.close();
		initial_config_ref.close();
		/*
		all_config.close();
		ind_collapse_txt.close();
		final_report.close();
		
		Wself_an_file.close();
		Wself_iso_file.close();
		Wself_ratio_file.close();
		r1coor_file.close();

		lt_file.close();
		core_file.close();
		Wself_file.close();
		Wself0_file.close();
		*/
	 	enthalpy_anisotropic_evolution.close();
	 	enthalpy_evolution.close();
		tension_anisotropic_evolution.close();
		tension_evolution.close();
		core_evolution.close();
		U_evolution.close();
		W_evolution.close();

	return 0;
}//end of main
