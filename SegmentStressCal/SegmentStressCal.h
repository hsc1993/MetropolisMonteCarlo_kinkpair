#ifndef _SEG_STRESS_CAL_H
#define _SEG_STRESS_CAL_H

#include "../StandardIncludes.h"
#include "../Settings.h"
#include "../simulation/Simulation.h"

#include <tuple>


/************************************************************
 * 		non-singular anisotropic stress functions
 * *********************************************************/

class SegmentStressClassical
{
    
    const Tensor<4,3>& C;               // Elastic constants
    GreenTensor green;                  // Classical anisotropic Green tensor. Input is the number of integration point
    const Eigen::Matrix<double,3,1> b;  // Burgers vector
    const Eigen::Matrix<double,3,1> x0; // starting point of the segment
    const Eigen::Matrix<double,3,1> x1; // end point of the segment
    const Eigen::Matrix<double,3,1> dL; // infinitesimal length along segment
    const Tensor<3,3> CbedL;            // product C*b*e*dL
    std::vector<Eigen::Matrix<double,3,1>> segmentSubdivisionCenters; // a set of equally-spaced points dicretizing the segment
        
	public:
    SegmentStressClassical(const Tensor<4,3>& Cin,
			      const Eigen::Matrix<double,3,3>& Lin,  //Eigen::Matrix<double,3,3>
                  const size_t& greenQuadPoints,
				  const size_t& integrationSteps,
                  const Eigen::Matrix<double,3,1>& bin,
                  const Eigen::Matrix<double,3,1>& x0in,
                  const Eigen::Matrix<double,3,1>& x1in,
                  const size_t& segmentSubDivisions ) :
        /* init */ C(Cin)
        /* init */,green(C,Lin,greenQuadPoints,integrationSteps)
        /* init */,b(bin)
        /* init */,x0(x0in)
        /* init */,x1(x1in)
        /* init */,dL((x1-x0)/segmentSubDivisions)
        /* init */,CbedL(cbedl(C,b,dL))
    {
        segmentSubdivisionCenters.reserve(segmentSubDivisions);
        segmentSubdivisionCenters.push_back(x0+0.5*dL);
        for(size_t k=0;k<segmentSubDivisions-1;++k)
        {
            segmentSubdivisionCenters.push_back(segmentSubdivisionCenters.back()+dL);
        }
        
    }

    Eigen::Matrix<double,3,3> stress(const Eigen::Matrix<double,3,1>& x) const;

    Eigen::Matrix<double,3,3> distortion(const Eigen::Matrix<double,3,1>& xF) const;
    
    Eigen::Matrix<double,3,3> strain(const Eigen::Matrix<double,3,1>& x) const;

    //Eigen::Matrix<double,4,3> FTensor(const Eigen::Matrix<double,3,1>& x) const;
    Tensor<4,3> FTensor(const Eigen::Matrix<double,3,1>& x) const;

    double Wself(const Eigen::Matrix<double,3,1>& x) const;

    Eigen::Matrix<double,3,1> getb() const {return b;};

    Eigen::Matrix<double,3,1> getLine() const {return x1-x0;};


};

/************************************************************
 * 		singular anisotropic stress functions
 * *********************************************************/

class SegmentStressClassical0
{
    
    const Tensor<4,3>& C;               // Elastic constants
    GreenTensor0 green0;                // Classical anisotropic Green tensor. Input is the number of integration points
    const Eigen::Matrix<double,3,1> b;  // Burgers vector
    const Eigen::Matrix<double,3,1> x0; // starting point of the segment
    const Eigen::Matrix<double,3,1> x1; // end point of the segment
    const Eigen::Matrix<double,3,1> dL; // infinitesimal length along segment
    const Tensor<3,3> CbedL;            // product C*b*e*dL
    std::vector<Eigen::Matrix<double,3,1>> segmentSubdivisionCenters; // a set of equally-spaced points dicretizing the segment
    
public:
    
    SegmentStressClassical0(const Tensor<4,3>& Cin,
                  const size_t& greenQuadPoints,
                  const Eigen::Matrix<double,3,1>& bin,
                  const Eigen::Matrix<double,3,1>& x0in,
                  const Eigen::Matrix<double,3,1>& x1in,
                  const size_t& segmentSubDivisions ) :
        /* init */ C(Cin)
        /* init */,green0(C,greenQuadPoints)
        /* init */,b(bin)
        /* init */,x0(x0in)
        /* init */,x1(x1in)
        /* init */,dL((x1-x0)/segmentSubDivisions)
        /* init */,CbedL(cbedl(C,b,dL))
    {
        segmentSubdivisionCenters.reserve(segmentSubDivisions);
        segmentSubdivisionCenters.push_back(x0+0.5*dL);
        for(size_t k=0;k<segmentSubDivisions-1;++k)
        {
            segmentSubdivisionCenters.push_back(segmentSubdivisionCenters.back()+dL);
        }
        
    }
    
    Eigen::Matrix<double,3,3> stress(const Eigen::Matrix<double,3,1>& x) const;
    
    Eigen::Matrix<double,3,3> distortion(const Eigen::Matrix<double,3,1>& xF) const;
    
    Eigen::Matrix<double,3,3> strain(const Eigen::Matrix<double,3,1>& x) const;

    Tensor<4,3> FTensor0(const Eigen::Matrix<double,3,1>& xF) const;

    double Wself0(const Eigen::Matrix<double,3,1>& xF) const;

    Eigen::Matrix<double,3,1> getb() const {return b;};

    Eigen::Matrix<double,3,1> getLine() const {return x1-x0;};

};



#endif  // _SEG_STRESS_CAL_H