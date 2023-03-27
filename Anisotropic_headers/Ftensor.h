#ifndef _FTENSOR_H
#define _FTENSOR_H


#include <iostream>
#include <cfloat>
#include <vector>
#include "Voigt2Tensor.h"
#include "Tensor.h"


class Ftensor
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
    
    /**************************************************************************/
    double ell(const Eigen::Matrix<double,3,1>& K) const
    {/*!Computes ell=L_{ij}*K_i*K_j
      */
        return sqrt(K.transpose()*L*K);
    }
    
    
public:
    
//    Eigen::Matrix<double,6,6> C; // Matrix of material constants
    const Tensor<4,3> cT;
    const Eigen::Matrix<double,3,3> L;               // gradient parameter
//    const Voigt2Tensor<4,3> cT; // tensor of material constants
    const int nPhi; // number of quadrature points, phi coordinate
    const int nQ; // number of quadrature points, q coordinate
    const double   dPhi; // integration step, phi coordinate
    const double   dQ; // integration step, q coordinate
    const double f;
    
    /**************************************************************************/
    Ftensor(const Tensor<4,3>& C_in, const Eigen::Matrix<double,3,3>& L_in,
                const int& nPhi_in, const int& nQ_in) :
    /* init list */ cT(C_in),
    /* init list */ L(L_in),
//    /* init list */ cT(C),
    /* init list */ nPhi(nPhi_in),
    /* init list */ nQ(nQ_in),
    /* init list */ dPhi(2.0*M_PI/nPhi),
    /* init list */ dQ(1.0/nQ),
    /* init list */ f(-dQ*dPhi/(8.0*M_PI*M_PI))
    {
        std::cout<<"Creating F-tensor"<<std::endl;
//        std::cout<<"C=\n"<<C<<std::endl;
        std::cout<<"L=\n"<<L<<std::endl;

    }
    
    /**************************************************************************/
    Tensor<4,3> operator()(const Eigen::Matrix<double,3,1>& x) const
    {
        
        const Eigen::Matrix<double,3,3> E(getLocalTriad(x));
        
        Tensor<4,3> F;
        
        for (int q=0;q<nQ;++q)       // integrate in q
        {
            for (int p=0;p<nPhi;++p) // integrate in phi
            {
                //          const double theta(0.5*dTheta+t*dTheta);
                const double phi(0.5*dPhi+p*dPhi); // current values of phi
                const double   Q(0.5*dQ+q*dQ);     // current values of q
                
                const double sinT(sqrt(1.0-Q*Q));
                const Eigen::Matrix<double,3,1> K(sinT*cos(phi)*E.col(0)+sinT*sin(phi)*E.col(1)+Q*E.col(2));
                assert(std::fabs(K.squaredNorm()-1.0)<FLT_EPSILON);
                const double eLL(ell(K));
                const Eigen::Matrix<double,3,3> temp(f*Lik(K).inverse()*exp(-x.norm()*Q/eLL)/eLL);
                
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
        }
        
        return F;
    }
    
    
    
};





#endif //_FTENSOR_H