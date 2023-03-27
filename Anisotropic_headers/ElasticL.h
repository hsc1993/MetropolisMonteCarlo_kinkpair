
#ifndef _ElasticL_h_
#define _ElasticL_h_

#include <array>
#include <math.h>       /* pow */
// #include <eigen3/Eigen/Dense>
// #include <Eigen/Dense>
#include <Eigen/Dense>


//#include <ElasticC.h>
//#include <ElasticD.h>
#include "Tensor.h"

struct Cubic{};
struct Triclinic{};
struct Orthorombic{};
struct Hexagonal{};

template <typename Symmetry>
class ElasticL;


/******************************************************************************/
template<>
class ElasticL<Cubic> : public Tensor<2,3>
{
    
    
public:
    ElasticL(const Tensor<4,3>& C, const Tensor<6,3>& D)
    {
        
        Tensor<4,3> Dmm;
        for(int i=0;i<3;++i)
        {
            for(int j=0;j<3;++j)
            {
                for(int k=0;k<3;++k)
                {
                    for(int l=0;l<3;++l)
                    {
                        for(int m=0;m<3;++m)
                        {
                            Dmm(i,j,k,l)+=D(i,j,m,k,l,m);
                        }
                    }
                }
            }
        }
        
//        for(int m=0;m<3;++m)
//        {
//            for(int n=0;n<3;++n)
//            {
                double num=0.0;
                double den=0.0;
                
                for(int i=0;i<3;++i)
                {
                    for(int j=0;j<3;++j)
                    {
                        for(int k=0;k<3;++k)
                        {
                            for(int l=0;l<3;++l)
                            {
                                num+=1.0/3.0*Dmm(i,j,k,l)*C(i,j,k,l);
                                den+=C(i,j,k,l)*C(i,j,k,l);
                            }
                        }
                    }
                }
        
        const double L2=num/den;
//                this->operator()(m,n)=num/den;
//            }
//        }
        
        for(int i=0;i<3;++i)
        {

                this->operator()(i,i)=L2;
        }
    }
    
    Eigen::Matrix<double,3,3> matrix() const
    {
        Eigen::Matrix<double,3,3> temp;
        for(int i=0;i<3;++i)
        {
            for(int j=0;j<3;++j)
            {
                temp(i,j)=this->operator()(i,j);
            }
        }
        return temp;
    }
    
};


/******************************************************************************/
template<>
class ElasticL<Orthorombic> : public Tensor<2,3>
{
    

public:
    ElasticL(const Tensor<4,3>& C, const Tensor<6,3>& D)
    {
        for(int m=0;m<3;++m)
        {
//            for(int n=0;n<3;++n)
//            {
                double num=0.0;
                double den=0.0;
                
                for(int i=0;i<3;++i)
                {
                    for(int j=0;j<3;++j)
                    {
                        for(int k=0;k<3;++k)
                        {
                            for(int l=0;l<3;++l)
                            {
                                num+=D(i,j,m,k,l,m)*C(i,j,k,l);
                                den+=C(i,j,k,l)*C(i,j,k,l);
                            }
                        }
                    }
                }
                
                this->operator()(m,m)=num/den;
            //}
        }
    }
    
    Eigen::Matrix<double,3,3> matrix() const
    {
        Eigen::Matrix<double,3,3> temp;
        for(int i=0;i<3;++i)
        {
            for(int j=0;j<3;++j)
            {
                temp(i,j)=this->operator()(i,j);
            }
        }
        return temp;
    }
    
};

/******************************************************************************/
template<>
class ElasticL<Hexagonal> : public Tensor<2,3>
{
    
    
public:
    ElasticL(const Tensor<4,3>& C, const Tensor<6,3>& D)
    {
        std::array<double,3> temp;
        for(int m=0;m<3;++m)
        {
            //            for(int n=0;n<3;++n)
            //            {
            double num=0.0;
            double den=0.0;
            
            for(int i=0;i<3;++i)
            {
                for(int j=0;j<3;++j)
                {
                    for(int k=0;k<3;++k)
                    {
                        for(int l=0;l<3;++l)
                        {
                            num+=D(i,j,m,k,l,m)*C(i,j,k,l);
                            den+=C(i,j,k,l)*C(i,j,k,l);
                        }
                    }
                }
            }
            
//            this->operator()(m,m)
            temp[m]=num/den;
            //}
        }
        
        this->operator()(0,0)=0.5*(temp[0]+temp[1]);
        this->operator()(1,1)=0.5*(temp[0]+temp[1]);
        this->operator()(2,2)=temp[2];
    }
    
    Eigen::Matrix<double,3,3> matrix() const
    {
        Eigen::Matrix<double,3,3> temp;
        for(int i=0;i<3;++i)
        {
            for(int j=0;j<3;++j)
            {
                temp(i,j)=this->operator()(i,j);
            }
        }
        return temp;
    }
    
};


/******************************************************************************/
template<>
class ElasticL<Triclinic> : public Tensor<2,3>
{
    
    
public:
    ElasticL(const Tensor<4,3>& C, const Tensor<6,3>& D)
    {
        for(int m=0;m<3;++m)
        {
            for(int n=0;n<3;++n)
            {
                double num=0.0;
                double den=0.0;
                
                for(int i=0;i<3;++i)
                {
                    for(int j=0;j<3;++j)
                    {
                        for(int k=0;k<3;++k)
                        {
                            for(int l=0;l<3;++l)
                            {
                                num+=0.5*(D(i,j,m,k,l,n)+D(i,k,n,k,l,m))*C(i,j,k,l);
                                den+=C(i,j,k,l)*C(i,j,k,l);
                            }
                        }
                    }
                }
                
                this->operator()(m,n)=num/den;
            }
        }
    }
    
    Eigen::Matrix<double,3,3> matrix() const
    {
        Eigen::Matrix<double,3,3> temp;
        for(int i=0;i<3;++i)
        {
            for(int j=0;j<3;++j)
            {
                temp(i,j)=this->operator()(i,j);
            }
        }
        return temp;
    }
    
};




//    /***********************************/
//    template <class T>
//    friend T& operator << (T& os, const ElasticL<>& L)
//    {
//        for(int m=0;m<3;++m)
//        {
//            for(int n=0;n<3;++n)
//            {
//                os<<std::setprecision(15)<<std::scientific<<L(m,n)<<" ";
//            }
//        }
//        os<<"\n";
//        return os;
//    }


#endif