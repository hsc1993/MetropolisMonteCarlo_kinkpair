#include "../StandardIncludes.h"
#include "../Settings.h"
#include <tuple>
#include "../Parameters.h"

#include "../simulation/Simulation.h"
#include "Dislocation.h"

// ordering matters: anisotropic_hearders.h must be put in front of SegmentStressCal.h
#include "../Anisotropic_headers.h"
#include "../SegmentStressCal/SegmentStressCal.h"


using namespace std;



/*******************************************************************************
* initialize all segments positions  x[1 0 0],y[0 1 0],z[0 0 1]
********************************************************************************/
/*
void DislocationNode::initialize(SimulationParameters& simParams, int _inode, double h_initial){
    
    if (_inode < simParams.n_segments*simParams.f1){
        _pos = {simParams.x_equi_ele,0,_inode*blength};
    }
    if (_inode >= simParams.n_segments*simParams.f1 && _inode < simParams.n_segments*(simParams.f1+simParams.fa1)){
        _pos = {simParams.x_equi_ele+h_initial/(simParams.n_segments*simParams.fa1)*_inode-h_initial*simParams.f1/simParams.fa1,0,_inode*blength};
    }
    if (_inode >= simParams.n_segments*(simParams.f1+simParams.fa1) && _inode < simParams.n_segments*(simParams.f1+simParams.fa1+simParams.fw)){
        _pos = {simParams.x_equi_ele+h_initial,0,_inode*blength};
    }
    if (_inode >= simParams.n_segments*(simParams.f1+simParams.fa1+simParams.fw) && _inode < simParams.n_segments*(simParams.f1+simParams.fa1+simParams.fw+simParams.fa2)){
        _pos = {simParams.x_equi_ele-h_initial/(simParams.n_segments*simParams.fa2)*_inode+h_initial*(simParams.f1+simParams.fa1+simParams.fw+simParams.fa2)/simParams.fa2,0,_inode*blength};
    }
    if (_inode >= simParams.n_segments*(simParams.f1+simParams.fa1+simParams.fw+simParams.fa2)){
        _pos = {simParams.x_equi_ele,0,_inode*blength};
    }
}
*/

void DislocationNode::initialize(SimulationParameters& simParams){
	
    _h_initial = simParams.h_initial[_indh];
	_x_equi = simParams.x_equi_ele;


    if (_inode < simParams.n_segments*simParams.f1){
        _position << _x_equi,0,_inode*blength;
    }
    if (_inode >= simParams.n_segments*simParams.f1 && _inode < simParams.n_segments*(simParams.f1+simParams.fa1)){
        _position << _x_equi+_h_initial/(simParams.n_segments*simParams.fa1)*_inode-_h_initial*simParams.f1/simParams.fa1,0,_inode*blength;
    }
    if (_inode >= simParams.n_segments*(simParams.f1+simParams.fa1) && _inode < simParams.n_segments*(simParams.f1+simParams.fa1+simParams.fw)){
        _position << _x_equi+_h_initial,0,_inode*blength;
    }
    if (_inode >= simParams.n_segments*(simParams.f1+simParams.fa1+simParams.fw) && _inode < simParams.n_segments*(simParams.f1+simParams.fa1+simParams.fw+simParams.fa2)){
        _position << _x_equi-_h_initial/(simParams.n_segments*simParams.fa2)*_inode+_h_initial*(simParams.f1+simParams.fa1+simParams.fw+simParams.fa2)/simParams.fa2,0,_inode*blength;
    }
    if (_inode >= simParams.n_segments*(simParams.f1+simParams.fa1+simParams.fw+simParams.fa2)){
        _position << _x_equi,0,_inode*blength;
    }
    

}


void DislocationNode::calTheta(SimulationParameters& simParams)
{
    Eigen::Matrix<double,3,1> r = _positionTransformed;
    Eigen::Matrix<double,3,1> r_next = this->nextNode()->positionTransformed();
    Eigen::Matrix<double,3,1> dr = r_next-r;
    Eigen::Matrix<double,1,1> costheta_matrix = dr.transpose()*simParams.lineDirection/dr.norm();
    double costheta = costheta_matrix(0);
    if (costheta>1){
        costheta = 1;
    }
    if (costheta<-1){
        costheta = -1;
    }
    this->set_theta(acos(costheta));
}

void DislocationNode::calGlideLength(SimulationParameters& simParams)
{
    _distanceToOrigin = _positionTransformed.norm();

    if (_distanceToOrigin == 0 || _distanceToOrigin<0){ //tau=0 case, first node's distance to origin is zero, leading to nan bug
        _cos_thetaWrtOrigin = 1; 
    }
    else{
        Eigen::Matrix<double,1,1> cos_thetaWrtOrigin_matrix = (this->positionTransformed()).transpose()*simParams.lineDirection/_distanceToOrigin;
	    _cos_thetaWrtOrigin = cos_thetaWrtOrigin_matrix(0);
    }

    if (_cos_thetaWrtOrigin>1){
        _cos_thetaWrtOrigin = 1;
    }
    if (_cos_thetaWrtOrigin<-1){
        _cos_thetaWrtOrigin = -1;
    }
	_thetaWrtOrigin = acos(_cos_thetaWrtOrigin);

	if(_positionTransformed.transpose()*simParams.glideDirection>0){
		this->set_glideLength(_distanceToOrigin*sin(_thetaWrtOrigin));
	}
	else{
		this->set_glideLength(-(_distanceToOrigin*sin(_thetaWrtOrigin)));
	}
}


/*******************************************************************************
* transforming trapezoid into x[-1 2 -1],y[-1 0 1],z[1 1 1]
********************************************************************************/

Eigen::Matrix<double,3,1> DislocationNode::transform(Eigen::Matrix<double,3,1> vector, Eigen::Matrix<double,3,3> T)
{
   return vector.transpose()*T;
}



void DislocationNode::ScrewSegWself(SegmentStressClassical* segStressCla, Eigen::Matrix<double,3,1> fieldPoint)
{
	this->set_screwSegWself(segStressCla->Wself(fieldPoint));
}


void DislocationNode::Tension(SimulationParameters& simParams)
{
    /// straight line tension
    double cs_straight = cos(0);
    double term2_straight = (1-simParams.pr*pow(cs_straight,2))*log((blength+sqrt(pow(blength,2)+pow(simParams.a,2)))/simParams.a);
    double term3_straight = (3-simParams.pr)/2*(sqrt(pow(blength,2)+pow(simParams.a,2))-simParams.a)/blength*(pow(cs_straight,2));
    double T_straight = simParams.term1*(term2_straight-term3_straight)*simParams.unitt;

    double L = blength*100/simParams.n_segments;
    double La = sqrt(pow(L,2)+pow(simParams.a,2)); 

    double cs = cos(this->theta());
    double z_calcu = abs(blength/cs);
    double term2 = (1-simParams.pr*pow(cs,2))*log((blength+sqrt(pow(blength,2)+pow(simParams.a,2)))/simParams.a);
    double term3 = (3-simParams.pr)/2.0*(sqrt(pow(blength,2)+pow(simParams.a,2))-simParams.a)/blength*(pow(cs,2));
    double T = simParams.term1*(term2-term3)*simParams.unitt;
    double H1_T = T/2.0*z_calcu-T_straight/2.0*blength;
    this->set_tension(H1_T);
}


void DislocationNode::Tension_anisotropic(SegmentStressClassical* segStressCla)
{
	this->set_tension_anisotropic(segStressCla->Wself(this->xF())-this->screwSegWself());
}

void DislocationNode::Core(EcoreParameters& core, SimulationParameters& simParams)
{
	// for Iron, ecore = core.c1+core.c2*cos(2*theta)
		/// straight line core
		double thetaStrai_left = 0.0;
		double ECoreStrai_left = core.c1+core.c2*cos(2*thetaStrai_left);
		double thetaStrai_right = M_PI;
		double ECoreStrai_right = core.c1+core.c2*cos(2*thetaStrai_right);

		bool flag_right = 0;
		if (this->theta()<0){
			this->set_theta(this->theta()+M_PI);
			bool flag_right = 1;
		}
		
		double ECore = core.c1+core.c2*cos(2*this->theta());
		double H1_C;
		if (flag_right == 0){
			H1_C = (ECore-ECoreStrai_left)*blength;
		}else{
			H1_C = (ECore-ECoreStrai_right)*blength;
		}
		this->set_core(H1_C);
}

void DislocationNode::UPeierls(SimulationParameters& simParams)
{
	double x = this->glideLength();
	//double U = Up_vs_x(simParams, this->glideLength());
	double U = simParams.U0/2.0*(1-cos(2.0*M_PI*x/simParams.h0));
	double U_peierls= U*blength-simParams.U_equi_ele*blength;
	this->set_U(U_peierls);
}

double DislocationNode::Up_vs_x(SimulationParameters& simParams, double x)
{
	return simParams.U0/2.0*(1-cos(2.0*M_PI*x/simParams.h0));
}

void DislocationNode::Work(SimulationParameters& simParams)
{
	double x = this->glideLength();
	double x_next = this->nextNode()->glideLength();
	double W_return;

	if (x > simParams.x_initial[_indh]){
		W_return = -simParams.currentTau*blength*blength/2.0*(-(x_next+x-2.0*simParams.x_equi_ele));
		//W_return = -simParams.currentTau*blength*blength/2.0*(x_next+x-2.0*simParams.x_equi_ele);
	}
	else{
		W_return = -simParams.currentTau*blength*blength/2.0*(x_next+x-2.0*simParams.x_equi_ele);
	}
	this->set_W(W_return);
}

void DislocationNode::calNodeEnergy(EcoreParameters& ecoreParams, SimulationParameters& simParams)
{
	this->Tension(simParams);
	this->Core(ecoreParams,simParams);
	this->UPeierls(simParams);
	this->Work(simParams);
	// this->set_enthalpy(_tension + _core + _U + _W);
	this->set_enthalpy(_tension + _U + _W);
}


void DislocationNode::calNodeEnergy_anisotropic(SegmentStressClassical* segStressCla, EcoreParameters& ecoreParams, SimulationParameters& simParams)
{
	this->Tension_anisotropic(segStressCla);
	this->Core(ecoreParams,simParams);
	this->UPeierls(simParams);
	this->Work(simParams);
	this->set_enthalpy_anisotropic(this->tension_anisotropic() + this->core() + this->U() + this->W());

}



//***************** Dislocation Loop calculation functions *****************//

void DislocationLoop::calEnthalpyStep(EcoreParameters& ecoreParams, SimulationParameters& simParams, 
					DislocationNode* pnodeTrial1, DislocationNode* pnodeTrial2,
					DislocationNode* pnode1, DislocationNode* pnode2)
{
	/// calculate change in tension and incorporate them into nodeTrial object
	pnodeTrial1->Tension(simParams);
	pnodeTrial2->Tension(simParams);

	/// calculate change in core and incorporate them into nodeTrial object
	pnodeTrial1->Core(ecoreParams, simParams);
	pnodeTrial2->Core(ecoreParams, simParams);

	/// calculate change in U and incorporate them into nodeTrial object
	pnodeTrial1->UPeierls(simParams);
	pnodeTrial2->UPeierls(simParams);

	/// calculate change in W and incorporate them into nodeTrial object
	pnodeTrial1->Work(simParams);
	pnodeTrial2->Work(simParams);

	// pnodeTrial1->set_enthalpy(pnodeTrial1->tension() + pnodeTrial1->core() + pnodeTrial1->U() + pnodeTrial1->W());
	// pnodeTrial2->set_enthalpy(pnodeTrial2->tension() + pnodeTrial2->core() + pnodeTrial2->U() + pnodeTrial2->W());
	pnodeTrial1->set_enthalpy(pnodeTrial1->tension() + pnodeTrial1->U() + pnodeTrial1->W());
	pnodeTrial2->set_enthalpy(pnodeTrial2->tension() + pnodeTrial2->U() + pnodeTrial2->W());

	/// calculate each term for previous configuration's
	pnode1->Tension(simParams);
	pnode2->Tension(simParams);

	pnode1->Core(ecoreParams, simParams);
	pnode2->Core(ecoreParams, simParams);

	pnode1->UPeierls(simParams);
	pnode2->UPeierls(simParams);

	pnode1->Work(simParams);
	pnode2->Work(simParams);

	// pnode1->set_enthalpy(pnode1->tension() + pnode1->core() + pnode1->U() + pnode1->W());
	// pnode2->set_enthalpy(pnode2->tension() + pnode2->core() + pnode2->U() + pnode2->W());
	pnode1->set_enthalpy(pnode1->tension() + pnode1->U() + pnode1->W());
	pnode2->set_enthalpy(pnode2->tension() + pnode2->U() + pnode2->W());


 	_tensionStep = pnodeTrial1->tension()+pnodeTrial2->tension()-pnode1->tension()-pnode2->tension();
	_coreStep = pnodeTrial1->core()+pnodeTrial2->core()-pnode1->core()-pnode2->core();
	_UStep = pnodeTrial1->U()+pnodeTrial2->U()-pnode1->U()-pnode2->U();
	_WStep = pnodeTrial1->W()+pnodeTrial2->W()-pnode1->W()-pnode2->W();
	_enthalpyStep = pnodeTrial1->enthalpy()+pnodeTrial2->enthalpy()-pnode1->enthalpy()-pnode2->enthalpy();
}


void DislocationLoop::calEnthalpyStep_anisotropic(SegmentStressClassical* segStressCla_trail1, SegmentStressClassical* segStressCla_trail2, 
											SegmentStressClassical* segStressCla1, SegmentStressClassical* segStressCla2, 
											EcoreParameters& ecoreParams, SimulationParameters& simParams, 
											DislocationNode* pnodeTrial1, DislocationNode* pnodeTrial2,
											DislocationNode* pnode1, DislocationNode* pnode2)
{
	/// calculate change in tension and incorporate them into nodeTrial object
	pnodeTrial1->Tension_anisotropic(segStressCla_trail1);
	pnodeTrial2->Tension_anisotropic(segStressCla_trail2);

	/// calculate change in core and incorporate them into nodeTrial object
	pnodeTrial1->Core(ecoreParams, simParams);
	pnodeTrial2->Core(ecoreParams, simParams);

	/// calculate change in U and incorporate them into nodeTrial object
	pnodeTrial1->UPeierls(simParams);
	pnodeTrial2->UPeierls(simParams);

	/// calculate change in W and incorporate them into nodeTrial object
	pnodeTrial1->Work(simParams);
	pnodeTrial2->Work(simParams);

	pnodeTrial1->set_enthalpy_anisotropic(pnodeTrial1->tension_anisotropic() + pnodeTrial1->core() + pnodeTrial1->U() + pnodeTrial1->W());
	pnodeTrial2->set_enthalpy_anisotropic(pnodeTrial2->tension_anisotropic() + pnodeTrial2->core() + pnodeTrial2->U() + pnodeTrial2->W());


	/// calculate each term for previous configuration's
	
	pnode1->Tension_anisotropic(segStressCla1);
	pnode2->Tension_anisotropic(segStressCla2);

	pnode1->Core(ecoreParams, simParams);
	pnode2->Core(ecoreParams, simParams);

	pnode1->UPeierls(simParams);
	pnode2->UPeierls(simParams);

	pnode1->Work(simParams);
	pnode2->Work(simParams);
	

	pnode1->set_enthalpy_anisotropic(pnode1->tension_anisotropic() + pnode1->core() + pnode1->U() + pnode1->W());
	pnode2->set_enthalpy_anisotropic(pnode2->tension_anisotropic() + pnode2->core() + pnode2->U() + pnode2->W());

/*
	_tension_anisotropicStep;
	_coreStep;
	_UStep;
	_WStep;
	_enthalpy_anisotropicStep;
	*/
}



void DislocationLoop::calLoopEnergy(SimulationParameters& simParams)
{
	this->set_enthalpy(0);
	this->set_tension(0);
	this->set_core(0);
	this->set_U(0);
	this->set_W(0);
	DislocationNode* pnode;
	pnode = _firstNode;
	for (int i = 0; i < simParams.n_segments-1; i++){
		this->set_enthalpy(pnode->enthalpy() + this->enthalpy());
		this->set_tension(pnode->tension() + this->tension());
		this->set_core(pnode->core() + this->core());
		this->set_U(pnode->U() + this->U());
		this->set_W(pnode->W() + this->W());
		pnode = pnode->nextNode();
	}
}

void DislocationLoop::calLoopEnergy_anisotropic(SimulationParameters& simParams)
{
	DislocationNode* pnode;
	pnode = _firstNode;	
	this->set_enthalpy_anisotropic(0);
	this->set_tension_anisotropic(0);
	this->set_core(0);
	this->set_U(0);
	this->set_W(0);

	for (int i = 0; i < simParams.n_segments-1; i++){
		this->set_enthalpy_anisotropic(pnode->enthalpy_anisotropic() + this->enthalpy_anisotropic());
		this->set_tension_anisotropic(pnode->tension_anisotropic() + this->tension_anisotropic());
		this->set_core(pnode->core() + this->core());
		this->set_U(pnode->U() + this->U());
		this->set_W(pnode->W() + this->W());
		pnode = pnode->nextNode();
	}
}


			