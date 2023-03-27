
#ifndef _GREENTENSOR_H0
#define _GREENTENSOR_H0


#include <iostream>
#include <vector>
#include "Voigt2Tensor.h"
#include "Tensor.h"



class GreenTensor0
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
    
    
    /**************************************************************************
    * Alan: equation (10)
    * ************************************************************************/
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
    
    const Tensor<4,3> cT;
    //    const Eigen::Matrix<double,3,3> L;               // gradient parameter
    //    const Voigt2Tensor<4,3> cT; // tensor of material constants
    const int nPhi; // number of quadrature points, phi coordinate
    //    const int nQ; // number of quadrature points, q coordinate
    const double   dPhi; // integration step, phi coordinate
    //    const double   dQ; // integration step, q coordinate
    const double f;
    
    /**************************************************************************/
    GreenTensor0(const Tensor<4,3>& C_in,
                 const int& nPhi_in) :
    /* init list */ cT(C_in),
    /* init list */ nPhi(nPhi_in),
    //  /* init list */ dPhi(M_PI/nPhi),
    //  /* init list */ f(dPhi/(4.0*M_PI*M_PI))  // different from Ftensor0.h, need to make changes in greenFtensor0 function

    /* init list */ dPhi(2.0*M_PI/nPhi),
    /* init list */ f(-dPhi/(8.0*M_PI*M_PI))

    {
        //std::cout<<"Creating GreenTensor0"<<std::endl;
        //std::cout<<"C=\n"<<<<std::endl;
        
    }
    
    //    /**************************************************************************/
    
    Tensor<4,3> greenFtensor0(const Eigen::Matrix<double,3,1>& x) const
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

    
    Tensor<3,3> grad(const Eigen::Matrix<double,3,1>& x) const
    {

        const Eigen::Matrix<double,3,3> E(getLocalTriad(x));

        //std::cout 
        //<< "E = " << std::endl
        //<< E << std::endl
        //<< "getLocalTriad(x) = " << std::endl
        //<< getLocalTriad(x) << std::endl;

        
        Tensor<3,3> dG;
        const double x2=x.squaredNorm();
        
        
        for (int p=0;p<nPhi;++p) // integrate in phi
        {
            //          const double theta(0.5*dTheta+t*dTheta);
            const double phi(0.5*dPhi+p*dPhi); // current values of phi
            const Eigen::Matrix<double,3,1> K(cos(phi)*E.col(0)+sin(phi)*E.col(1));
            assert(std::fabs(K.squaredNorm()-1.0)<FLT_EPSILON);
            
            // Compute inverse L
            Eigen::Matrix<double,3,3> Linv=Lik(K).inverse();

            
            // Compute F_ir
            Eigen::Matrix<double,3,3> c(Eigen::Matrix<double,3,3>::Zero());
            for(int j=0;j<3;++j){
                for(int n=0;n<3;++n){
                    for(int p=0;p<3;++p){
                        for(int w=0;w<3;++w){
                            c(j,n) += cT(j,p,n,w)*(K(p)*E(w,2)+K(w)*E(p,2));
                        }
                    }
                }
            }
            Eigen::Matrix<double,3,3> F(Linv*c*Linv);

            
            
            for(int i=0;i<3;++i){
                for(int j=0;j<3;++j){
                    for(int m=0;m<3;++m){
                        dG(i,j,m) += (-E(m,2)*Linv(i,j)+K(m)*F(i,j))*f/x2;
                    }
                }
            }
        }
        
        return dG;
    }


    
    
};





#endif  // _GREENTENSOR_H0