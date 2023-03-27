#include <Eigen/Dense>
#include <Tensor.h>
#include <leviCivita.h>

Tensor<4,3> cbe(const Tensor<4,3>& C,
                  const Eigen::Matrix<double,3,1>& b)
{
    /* computes the tensor
     * cbe(m,n,l,r)=C(m,n,p,q)*b(p)*leviCivita(q,r,l)
     */
    
    Tensor<3,3> Cb;
    for (size_t m=0;m<3;++m)
    {
        for (size_t n=0;n<3;++n)
        {
            for (size_t q=0;q<3;++q)
            {
                for (size_t p=0;p<3;++p)
                {
                    Cb(m,n,q)+=C(m,n,p,q)*b(p);
                }
            }
        }
    }
    
//    Tensor<2,3> edL;
//    for (size_t q=0;q<3;++q)
//    {
//        for (size_t l=0;l<3;++l)
//        {
//            for (size_t r=0;r<3;++r)
//            {
//                edL(q,l)+=leviCivita(q,r,l)*dL(r);
//            }
//        }
//    }
    
    Tensor<4,3> Cbe;
    for (size_t m=0;m<3;++m)
    {
        for (size_t n=0;n<3;++n)
        {
            for (size_t l=0;l<3;++l)
            {
                for (size_t q=0;q<3;++q)
                {
                    for (size_t r=0;r<3;++r)
                    {
                    Cbe(m,n,l,r)+=Cb(m,n,q)*leviCivita(q,r,l);
                    }
                }
            }
        }
    }
    return Cbe;
}
