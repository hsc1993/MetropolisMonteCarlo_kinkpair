#ifndef _FTENSOR_H0
#define _FTENSOR_H0


#include <iostream>
#include <cfloat>
#include <vector>
#include <Voigt2Tensor.h>
#include <Tensor.h>


class Ftensor0
{
    
    /**************************************************************************/
    Eigen::Matrix<double,3,1> getRandomNormal(const Eigen::Matrix<double,3,1>& x) const
    {
        Eigen::Matrix<double,3,1> temp(Eigen::Matrix<double,3,1>::Random());
        Eigen::Matrix<double,3,1> temp1(temp.cross(x));
        
        while (temp1.squaredNorm()<FLT_EPSILON)
        {
            temp=Eigen::Matrix<double,3,1>::Random();
            temp1=temp.cross(x);
        }
        
        return temp1.normalized();
    }
    
    
    
    
    /**************************************************************************/
    Eigen::Matrix<double,3,3> getLocalTriad(const Eigen::Matrix<double,3,1>& x) const
    {
        Eigen::Matrix<double,3,3> temp;
        
        if (x.squaredNorm()>0.0)
        {
            temp.col(2)=x.normalized();
        }
        else
        {
            //temp.col(2)=Eigen::Matrix<double,3,1>::Random().normalized();
            temp.col(2)<<1.0, 0.0, 0.0;
        }
        
        temp.col(0)=getRandomNormal(temp.col(2));
        temp.col(1)=temp.col(2).cross(temp.col(0));
        
        return temp;
    }
    
    
    /**************************************************************************/
    Eigen::Matrix<double,3,3> Lik(const Eigen::Matrix<double,3,1>& K) const
    {/*!Computes Lik=C_ijkl*K_j*K_l
      */
        Eigen::Matrix<double,3,3> L(Eigen::Matrix<double,3,3>::Zero());
        
        for (int i=0;i<3;++i)
        {
            for (int k=0;k<3;++k)
            {
                for (int l=0;l<3;++l)
                {
                    for (int j=0;j<3;++j)
                    {
                        L(i,k)+= cT(i,j,k,l)*K(l)*K(j);
                    }
                }
                
            }
        }
        
        return L;
    }
    
//    /**************************************************************************/
//    double ell(const Eigen::Matrix<double,3,1>& K) const
//    {/*!Computes ell=L_{ij}*K_i*K_j
//      */
//        return sqrt(K.transpose()*L*K);
//    }
    
    
public:
    
//    Eigen::Matrix<double,6,6> C; // Matrix of material constants
    const Tensor<4,3> cT;
//    const Voigt2Tensor<4,3> cT; // tensor of material constants
    const int nPhi; // number of quadrature points, phi coordinate
    const double   dPhi; // integration step, phi coordinate
    const double f;
    
    /**************************************************************************/
    Ftensor0(const Tensor<4,3>& C_in,
                const int& nPhi_in) :
    /* init list */ cT(C_in),
//    /* init list */ cT(C),
    /* init list */ nPhi(nPhi_in),
    /* init list */ dPhi(2.0*M_PI/nPhi),
    /* init list */ f(-dPhi/(8.0*M_PI*M_PI))
    {
        std::cout<<"Creating F-tensor0"<<std::endl;
//        std::cout<<"C=\n"<<C<<std::endl;

    }
    
    /**************************************************************************/
    Tensor<4,3> operator()(const Eigen::Matrix<double,3,1>& x) const
    {
        
        const Eigen::Matrix<double,3,3> E(getLocalTriad(x));
        
        Tensor<4,3> F;
        

            for (int p=0;p<nPhi;++p) // integrate in phi
            {
                //          const double theta(0.5*dTheta+t*dTheta);
                const double phi(0.5*dPhi+p*dPhi); // current values of phi
                
                const Eigen::Matrix<double,3,1> K(cos(phi)*E.col(0)+sin(phi)*E.col(1));
                assert(std::fabs(K.squaredNorm()-1.0)<FLT_EPSILON);
                const Eigen::Matrix<double,3,3> temp(f*Lik(K).inverse()/x.norm());
                
                for(int i=0;i<3;++i)
                {
                    for(int j=0;j<3;++j)
                    {
                        for(int k=0;k<3;++k)
                        {
                            for(int m=0;m<3;++m)
                            {
                                F(i,j,k,m) += temp(k,m)*K(i)*K(j);
                            }
                        }
                    }
                }
            }
    
        
        return F;
    }
    
    
    
};





#endif   // _FTENSOR_H0
