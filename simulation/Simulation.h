#ifndef __SIMULATION_H
#define __SIMULATION_H

#include "../StandardIncludes.h"
#include "../Settings.h"
#include <tuple>
#include "../dislocation/Dislocation.h"


//#include "../SegmentStressCal/SegmentStressCal.h"
/*
 * Stores the precomputed information
 * on each possible kink direction.
 */


class SegmentStressClassical;
class SegmentStressClassical0;




class Simulation 
{
	public:
	/// Constructor
	//Simulation();

	/// Returns a reference to the parameter set.
	SimulationParameters& params() { return _params; }

	void initialize(SimulationParameters& simParams);
	
	void calStressDependents_W(SimulationParameters& simParams, PotentialFitParams_W& potFitParams);

	void calStressDependents(SimulationParameters& simParams, PotentialFitParams& potFitParams, DislocationNode* pnode);

	void updateStressDependents(SimulationParameters& simParams, int indTau);

	Eigen::Matrix<double,3,3> createTransformationMatrix();

	void rotate_tensor_4th(Tensor<4,3> Cref,Eigen::Matrix<double,3,3> T, SimulationParameters& simParams);

	void rotate_tensor_3th(Eigen::Matrix<double,3,3> Lref, Eigen::Matrix<double,3,3> T, SimulationParameters& simParams);

	private:
	/// The input parameter set.
	SimulationParameters _params;

};


#endif // __SIMULATION_H
