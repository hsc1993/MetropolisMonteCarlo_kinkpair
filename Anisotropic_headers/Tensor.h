
#ifndef _Tensor_h_
#define _Tensor_h_

#include <array>       /* pow */
#include <string>
template<int N,int dim>
struct ComponentLabels
{
    
    static std::string label(const std::string& s)
    {
        std::string temp;
        for(int i=0;i<dim;++i)
        {
            temp+=ComponentLabels<N-1,dim>::label(s+std::to_string(i+1));
        }
        return temp;
    }
    
    /***********************************/
    template <class A>
    friend A& operator << (A& os, const ComponentLabels<N,dim>& )
    {
        os<<label("");
        return os;
    }
    
};

template<int dim>
struct ComponentLabels<1,dim>
{
    
    static std::string label(const std::string& s)
    {
        std::string temp;
        for(int i=0;i<dim;++i)
        {
            temp+=s+std::to_string(i+1)+" ";
        }
        return temp;
    }
    
};

template<int N>
struct TensorRotation
{};


template<int N,int dim,typename Scalar=double>
struct Tensor : std::array<Tensor<N-1,dim,Scalar>,dim>
{
    typedef Tensor<N,dim,Scalar> TensorType;
    
    template<typename...Ints>
    Scalar& operator()(const int&i,const Ints&...ints)
    {
        return this->operator[](i)(ints...);
    }
    
    template<typename...Ints>
    const Scalar& operator()(const int&i,const Ints&...ints) const
    {
        return this->operator[](i)(ints...);
    }
    
    /***********************************/
    Tensor<N,dim,Scalar> rotate(const Eigen::Matrix<double,dim,dim>& Q) const
    {
        return TensorRotation<N>::rotate(*this,Q);
    }
    
    /***********************************/
    template <class T>
    friend T& operator << (T& os, const Tensor<N,dim,Scalar>& t)
    {
        for(int i=0;i<dim;++i)
        {
            os<<t[i];
        }
        return os;
    }
    
};

template<int dim,typename Scalar>
class Tensor<1,dim,Scalar> : public std::array<Scalar,dim>
{
    
    Scalar data[dim];
    
public:
    
    Tensor()
    {
        this->fill(0.0);
    }
    
    Scalar& operator()(const int&i)
    {
        return this->operator[](i);
    }
    
    const Scalar& operator()(const int&i) const
    {
        return this->operator[](i);
    }
    
    /***********************************/
    template <class T>
    friend T& operator << (T& os, const Tensor<1,dim,Scalar>& t)
    {
        for(int i=0;i<dim;++i)
        {
            os<<t[i]<<" ";
        }
        return os;
    }
    
    
};

template<>
struct TensorRotation<2>
{
    
    template<int dim,typename Scalar>
    static Tensor<2,dim> rotate(const Tensor<2,dim,Scalar>& T,
                                const Eigen::Matrix<double,dim,dim>& Q)
    {
        
        Tensor<2,dim,Scalar> temp;
        for (int i=0;i<dim;++i)
        {
            for (int j=0;j<dim;++j)
            {
                for (int a=0;a<dim;++a)
                {
                    for (int b=0;b<dim;++b)
                    {
                        temp(i,j)+=T(a,b)*Q(a,i)*Q(b,j);
                    }
                }
            }
        }
        return temp;
    }
    
};

template<>
struct TensorRotation<3>
{
    template<int dim,typename Scalar>
    static Tensor<3,dim> rotate(const Tensor<3,dim,Scalar>& T,
                                const Eigen::Matrix<double,dim,dim>& Q)
    {
        Tensor<3,dim> temp;
        for (int i=0;i<dim;++i)
        {
            for (int j=0;j<dim;++j)
            {
                for (int k=0;k<dim;++k)
                {
                    for (int a=0;a<dim;++a)
                    {
                        for (int b=0;b<dim;++b)
                        {
                            for (int c=0;c<dim;++c)
                            {
                                temp(i,j,k)+=T(a,b,c)*Q(a,i)*Q(b,j)*Q(c,k);
                            }
                        }
                    }
                }
            }
        }
        return temp;
    }
    
};

template<>
struct TensorRotation<4>
{
    template<int dim,typename Scalar>
    static Tensor<4,dim> rotate(const Tensor<4,dim,Scalar>& T,
                                const Eigen::Matrix<double,dim,dim>& Q)
    {
        Tensor<4,dim> temp;
        for (int i=0;i<dim;++i)
        {
            for (int j=0;j<dim;++j)
            {
                for (int k=0;k<dim;++k)
                {
                    for (int l=0;l<dim;++l)
                    {
                        for (int a=0;a<dim;++a)
                        {
                            for (int b=0;b<dim;++b)
                            {
                                for (int c=0;c<dim;++c)
                                {
                                    for (int d=0;d<dim;++d)
                                    {
                                        temp(i,j,k,l)+=T(a,b,c,d)*Q(a,i)*Q(b,j)*Q(c,k)*Q(d,l);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        return temp;
    }
    
};


#endif
