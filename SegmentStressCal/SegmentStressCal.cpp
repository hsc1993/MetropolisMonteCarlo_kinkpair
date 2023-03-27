#include "../StandardIncludes.h"
#include "../Settings.h"
#include "../simulation/Simulation.h"
#include "../Anisotropic_headers.h"
#include "SegmentStressCal.h"



/************************************************************
 * 		non-singular anisotropic stress functions
 * *********************************************************/


Eigen::Matrix<double,3,3> SegmentStressClassical::stress(const Eigen::Matrix<double,3,1>& x) const
{
    Eigen::Matrix<double,3,3> beta(distortion(x));
    Eigen::Matrix<double,3,3> sigma=Eigen::Matrix<double,3,3>::Zero();
    for (size_t i=0;i<3;++i){
        for (size_t j=0;j<3;++j){
            for (size_t k=0;k<3;++k){
                for (size_t l=0;l<3;++l){
                    sigma(i,j)+=C(i,j,k,l)*beta(k,l);
                }
            }
        }
    }
    return sigma;
}

Eigen::Matrix<double,3,3> SegmentStressClassical::distortion(const Eigen::Matrix<double,3,1>& xF) const
{
    Eigen::Matrix<double,3,3> beta0(Eigen::Matrix<double,3,3>::Zero());
    for (const auto& xS :  segmentSubdivisionCenters){
        const auto R=xF-xS;
        /**************** Alan commentted ***********/
        const Tensor<3,3> dG=green.grad(R);
        for (size_t k=0;k<3;++k){
            for (size_t l=0;l<3;++l){
                for (size_t m=0;m<3;++m){
                    for (size_t n=0;n<3;++n){
                        beta0(k,l)+=dG(k,m,n)*CbedL(m,n,l);
                    }
                }
            }
        }
    }
    return beta0;
}

Eigen::Matrix<double,3,3> SegmentStressClassical::strain(const Eigen::Matrix<double,3,1>& x) const
{
    Eigen::Matrix<double,3,3> beta(distortion(x));
    return 0.5*(beta+beta.transpose());
}

//Eigen::Matrix<double,4,3> SegmentStressClassical::FTensor(const Eigen::Matrix<double,3,1>& x) const
Tensor<4,3> SegmentStressClassical::FTensor(const Eigen::Matrix<double,3,1>& xF) const
{
    Tensor<4,3> Ftensor(green.greenFtensor(xF));
    return Ftensor;
}

double SegmentStressClassical::Wself(const Eigen::Matrix<double,3,1>& xF) const
{
    /*
    - Pick a xF on a dislocation line
    - initialize Wself=0
    - loop over xS (ok to use one point per segment)
    - compute R=xF-xS
    -evaluate F(R)
    - loop over indices to accumulate to Wself
    */

    double Wself_gpa = 0; // units of GPa
    double Wself_gpa_temp = 0;
    double Wself_eva = 0;

    for (const auto& xS :  segmentSubdivisionCenters){
        const auto R=xF-xS;
        //const auto R=xS-xS;
        //std::cout << "R = " << R.transpose() << std::endl;
        

        Tensor<4,3> Ftensor(green.greenFtensor(R)); //F(s,k,m,r)
        //std::cout << "Ftensor = " << Ftensor << std::endl;
        for (size_t s=0;s<3;++s){
            for (size_t k=0;k<3;++k){
                for (size_t m=0;m<3;++m){
                    for (size_t r=0;r<3;++r){
                        for (size_t n=0;n<3;++n){
                            Wself_gpa_temp += 0.5*CbedL(m,n,k)*CbedL(r,s,n)*Ftensor(s,k,m,r);
                        }
                    }
                }
            }
        }
        //std::cout << "Ftensor = " << Ftensor << std::endl;
    }
    
    Wself_gpa = Wself_gpa_temp/dL.norm();
    double unitt = 1e-21/(1.6*1e-19);
    Wself_eva = Wself_gpa*unitt; // units of ev

    return Wself_eva;
}



/************************************************************
 * 		singular anisotropic stress functions
 * *********************************************************/

Eigen::Matrix<double,3,3> SegmentStressClassical0::stress(const Eigen::Matrix<double,3,1>& x) const
{
    Eigen::Matrix<double,3,3> beta(distortion(x));
    Eigen::Matrix<double,3,3> sigma=Eigen::Matrix<double,3,3>::Zero();
    for (size_t i=0;i<3;++i){
        for (size_t j=0;j<3;++j){
            for (size_t k=0;k<3;++k){
                for (size_t l=0;l<3;++l){
                    sigma(i,j)+=C(i,j,k,l)*beta(k,l);
                }
            }
        }
    }
    return sigma;
}


Eigen::Matrix<double,3,3> SegmentStressClassical0::distortion(const Eigen::Matrix<double,3,1>& xF) const
{
    Eigen::Matrix<double,3,3> beta0(Eigen::Matrix<double,3,3>::Zero());
    for (const auto& xS :  segmentSubdivisionCenters){
        const auto R=xF-xS;
        /**************** Alan commentted ***********/
        const Tensor<3,3> dG0=green0.grad(R);

        //const Tensor<3,3> dG=green.grad(R);
        for (size_t k=0;k<3;++k){
            for (size_t l=0;l<3;++l){
                for (size_t m=0;m<3;++m){
                    for (size_t n=0;n<3;++n){
                        beta0(k,l)+=dG0(k,m,n)*CbedL(m,n,l);
                    }
                }
            }
        }
    }
    return beta0;   // it has no unit 
}

Eigen::Matrix<double,3,3> SegmentStressClassical0::strain(const Eigen::Matrix<double,3,1>& x) const
{
    Eigen::Matrix<double,3,3> beta(distortion(x));
    return 0.5*(beta+beta.transpose());
}

Tensor<4,3> SegmentStressClassical0::FTensor0(const Eigen::Matrix<double,3,1>& xF) const
{
    Tensor<4,3> Ftensor0(green0.greenFtensor0(xF));
    return Ftensor0;
}

double SegmentStressClassical0::Wself0(const Eigen::Matrix<double,3,1>& xF) const
{
    /*
    - Pick a xF on a dislocation line
    - initialize Wself=0
    - loop over xS (ok to use one point per segment)
    - compute R=xF-xS
    -evaluate F(R)
    - loop over indices to accumulate to Wself
    */

    double Wself0_gpa = 0;
    double Wself0_eva = 0;

    for (const auto& xS :  segmentSubdivisionCenters){
        //const auto R=xF-xS;
        const auto R=xF-xF;
        

        Tensor<4,3> Ftensor0(green0.greenFtensor0(R)); //F(s,k,m,r)
        //std::cout << "Ftensor0 = "<<Ftensor0<< std::endl;
        for (size_t s=0;s<3;++s){
            for (size_t k=0;k<3;++k){
                for (size_t m=0;m<3;++m){
                    for (size_t r=0;r<3;++r){
                        for (size_t n=0;n<3;++n){
                            Wself0_gpa += 0.5*CbedL(m,n,k)*CbedL(r,s,n)*Ftensor0(s,k,m,r);
                        }
                    }
                }
            }
        }
    }
    double unitt = 1e-21/(1.6*1e-19);
    Wself0_eva = Wself0_gpa*unitt; // units of ev

    return Wself0_eva;
}
