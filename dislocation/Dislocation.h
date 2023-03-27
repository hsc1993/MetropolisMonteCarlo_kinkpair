#ifndef __DISLOCATIONS_H
#define __DISLOCATIONS_H


#include "../StandardIncludes.h"
#include "../Settings.h"
#include <tuple>
#include "../Parameters.h"


/******************************************************************************
* originally from dislocation.h.
*****************************************************************************/

/// Stores the position of a node. X and Y components are always integer lattice positions, while Z can vary continuously.
typedef Point_3<int,int,double> NodalPosition;

/// Stores a vector between two nodes. X and Y components are always integer lattice positions, while Z can vary continuously.
typedef Vector_3<int,int,double> NodalVector;

/// Epsilon value that defines when two positions are considered equal.
#define NODAL_POS_EPSILON		1e-6

/// Epsilon value that defines when two nodal velocities are considered equal.
#define NODAL_VELOCITY_EPSILON	1e-6

typedef Eigen::Matrix<double, 3, 1> Vector3d;


class Simulation;
class SegmentStressClassical;
class SegmentStressClassical0;

class DislocationNode
{
	public:
		/// initializaing position of nodes
		void initialize(SimulationParameters& simParameters);

		Eigen::Matrix<double,3,3> createTransformationMatrix();

		Eigen::Matrix<double,3,1> transform(Eigen::Matrix<double,3,1> vector, Eigen::Matrix<double,3,3> T);

		void translate(double x, double y, double z);

		void calTheta(SimulationParameters& simParams);

		void calGlideLength(SimulationParameters& simParams);

		void ScrewSegWself(SegmentStressClassical* segStressCla, Eigen::Matrix<double,3,1> fieldPoint);
		
		void set_inode (const int inode_in) { _inode = inode_in; }
		void set_indh (const int indh_in) { _indh = indh_in; }
		void set_indTau (const int indTau_in) { _indTau = indTau_in; }
		void set_theta (const double theta_in) { _theta = theta_in; }
		void set_screwSegWself(const double screwSegWself_in) { _screwSegWself = screwSegWself_in; }
		void set_glideLength (const double glideLength_in) { _glideLength = glideLength_in; }
		void set_tension (const double tension_in) { _tension = tension_in; }
		void set_tension_anisotropic (const double tension_anisotropic_in) { _tension_anisotropic = tension_anisotropic_in; }
		void set_core (const double core_in) { _core = core_in; }
		void set_U (const double U_in) { _U = U_in; }
		void set_W (const double W_in) { _W = W_in; }
		void set_enthalpy  (const double enthalpy_in) { _enthalpy = enthalpy_in; }
		void set_enthalpy_anisotropic (const double enthalpy_anisotropic_in) { _enthalpy_anisotropic = enthalpy_anisotropic_in; }
		void set_x_equi (const double x_equi_in) { _x_equi = x_equi_in; }
		
		int inode() const { return _inode; }
		int indh() const { return _indh; }
		int indTau() const { return _indTau; }
		double theta() const { return _theta; }
		double screwSegWself() const { return _screwSegWself; }
		double glideLength() const { return _glideLength; }
		double tension() const { return _tension; }
		double tension_anisotropic() const { return _tension_anisotropic; }
		double core() const { return _core; }
		double U() const { return _U; }
		double W() const { return _W; }
		double enthalpy() const { return _enthalpy; }
		double enthalpy_anisotropic() const { return _enthalpy_anisotropic; }
		double x_equi() const { return _x_equi; }
		double h_initial_node() const { return _h_initial; }

		double distanceToOrigin() const { return _distanceToOrigin; }
		double cos_thetaWrtOrigin() const { return _cos_thetaWrtOrigin; }
		double thetaWrtOrigin() const { return _thetaWrtOrigin; }


		void Tension(SimulationParameters& simParams);

		void Tension_anisotropic(SegmentStressClassical* segStressCla);

		void Core(EcoreParameters& core, SimulationParameters& simParams);

		void UPeierls(SimulationParameters& simParams);

		double Up_vs_x(SimulationParameters& simParams, double x);

		void Work(SimulationParameters& simParams);

		void calNodeEnergy(EcoreParameters& ecoreParams, SimulationParameters& simParams);

		void calNodeEnergy_anisotropic(SegmentStressClassical* segStressCla, EcoreParameters& ecoreParams, SimulationParameters& simParams);
												

		/// Returns the next node in the dislocation loop.
		DislocationNode* nextNode() const { return _next; }

		/// Returns the previous node in the dislocation loop.
		DislocationNode* prevNode() const { return _prev; }

		/// Changes the node's position on the lattice grid.
		void setPosition(const Vector3d p) { _position = p; }

		void setPositionTransformed(const Vector3d p) { _positionTransformed = p; }

		void setxF(const Vector3d xF) { _xF = xF; }

		void setNext(DislocationNode* next) { _next = next;}

		/// print position
		void printPosition() { std::cout << _position << std::endl; }

		void printPosTransformed() { std::cout << _positionTransformed << std::endl; }

		Vector3d xF() { return _xF;}


		Vector3d position() { return _position;}
		double posX(){ return _position(0);}
		double posY(){ return _position(1);}
		double posZ(){ return _position(2);}

		Vector3d positionTransformed() { return _positionTransformed;}
		double posXTransformed(){ return _positionTransformed(0);}
		double posYTransformed(){ return _positionTransformed(1);}
		double posZTransformed(){ return _positionTransformed(2);}


	private:
		// field point of a dislocation node object
		Vector3d _xF;

		// coordinates in [100],[010],[001]
		Vector3d _position;

		// coordinates for [-1 2 -1],[-1 0 1],[1 1 1] kink pair
		Vector3d _positionTransformed;

		DislocationNode* _prev;

		DislocationNode* _next;
		int _inode;
		int _indh;
		int _indTau;

		double _x_equi;
		double _h_initial;

		double _distanceToOrigin;
		double _cos_thetaWrtOrigin;
		double _thetaWrtOrigin;
		
		double _theta;
		double _screwSegWself;
		double _glideLength;
		double _tension;
		double _tension_anisotropic;
		double _core;
		double _U;
		double _W;
		double _enthalpy;
		double _enthalpy_anisotropic;

};


/// Define the linear segment type. This is the same as the node type because they are stored
/// in the same structure.
typedef DislocationNode DislocationSegment;

/// Define handle types for nodes.
typedef DislocationNode* NodeHandle;

/// Define handle types for nodes.
typedef DislocationSegment* SegmentHandle;


/*
 * A dislocation line.
 */
class DislocationLoop
{
	public:
		void set_firstNode (DislocationNode* pnode) { _firstNode = pnode; }
		void set_tension (const double tension_in) { _tension = tension_in; }
		void set_tension_anisotropic (const double tension_anisotropic_in) { _tension_anisotropic = tension_anisotropic_in; }
		void set_core (const double core_in) { _core = core_in; }
		void set_U (const double U_in) { _U = U_in; }
		void set_W (const double W_in) { _W = W_in; }
		void set_enthalpy  (const double enthalpy_in) { _enthalpy = enthalpy_in; }
		void set_enthalpy_anisotropic (const double enthalpy_anisotropic_in) { _enthalpy_anisotropic = enthalpy_anisotropic_in; }


		double tension() const { return _tension; }
		double tension_anisotropic() const { return _tension_anisotropic; }
		double core() const { return _core; }
		double U() const { return _U; }
		double W() const { return _W; }
		double enthalpy() const { return _enthalpy; }
		double enthalpy_anisotropic() const { return _enthalpy_anisotropic; }

		void calEnthalpyStep(EcoreParameters& ecoreParams, SimulationParameters& simParams, 
							DislocationNode* pnodeTrial1, DislocationNode* pnodeTrial2,
							DislocationNode* pnode1, DislocationNode* node2);

		void calEnthalpyStep_anisotropic(SegmentStressClassical* segStressCla_trail1, SegmentStressClassical* segStressCla_trail2, 
											SegmentStressClassical* segStressCla1, SegmentStressClassical* segStressCla2, 
											EcoreParameters& ecoreParams, SimulationParameters& simParams, 
											DislocationNode* pnodeTrial1, DislocationNode* pnodeTrial2,
											DislocationNode* pnode1, DislocationNode* pnode2);

		void calLoopEnergy(SimulationParameters& simParams);

		void calLoopEnergy_anisotropic(SimulationParameters& simParams);

		double tensionStep() const { return _tensionStep; };
		double tension_anisotropicStep() const { return _tension_anisotropicStep; };
		double coreStep() const { return _coreStep; };
		double UStep() const { return _UStep; };
		double WStep() const { return _WStep; };
		double enthalpyStep() const { return _enthalpyStep; };
		double enthalpy_anisotropicStep() const { return _enthalpy_anisotropicStep; };

	private:

		/// Pointer to the head of the linked-list of nodes/segments.
		DislocationNode* _firstNode;

		double _tension;

		double _tension_anisotropic;

		double _core;

		double _U;
		
		double _W;
		
		double _enthalpy;
		
		double _enthalpy_anisotropic;

		double _tensionStep;
		double _tension_anisotropicStep;
		double _coreStep;
		double _UStep;
		double _WStep;
		double _enthalpyStep;
		double _enthalpy_anisotropicStep;
};




#endif // __DISLOCATIONS_H
