#ifndef _ANISOTROPIC_HEADERS_H
#define _ANISOTROPIC_HEADERS_H


inline void readCtensor(const std::string& dirName, Tensor<4,3>& C){
    int i,j,k,l;
    double val;
    std::string filename=dirName+"/Cijkl_ext.out";
    std::ifstream ifs ( filename.c_str() , std::ifstream::in );
    
    int numRead=0;
    if (ifs.is_open()){
        std::cout<<"Reading: "<<filename<<std::endl;
        std::string line;
        while (std::getline(ifs, line)){
            std::stringstream ss(line);
            ss>>i;
            ss>>j;
            ss>>k;
            ss>>l;
            ss>>val;
            C(i-1,j-1,k-1,l-1)=val;
            //std::cout<<"C_"<<i<<j<<k<<l<<"="<<C(i-1,j-1,k-1,l-1)<<std::endl;
            numRead++;
        }
        assert(numRead==81 && "WRONG NUMBER OF LINES IN C-file");
    }
    else{
        std::cout<<"FILE "<<filename<<" cannot be opened"<<std::endl;
        assert(0);
    }   
}

inline void readDtensor(const std::string& dirName,
                 Tensor<6,3>& D)
{
    
    int i,j,k,l,m,n;
    double val;
	int numRead=0;
    std::string filename=dirName+"/Dijmkln_ext.out";
    std::ifstream ifs2 ( filename.c_str() , std::ifstream::in );
    
    if (ifs2.is_open()){
        std::cout<<"Reading: "<<filename<<std::endl;
        std::string line;
        while (std::getline(ifs2, line)){
            std::stringstream ss(line);
            
            ss>>i;
            ss>>j;
            ss>>m;
            ss>>k;
            ss>>l;
            ss>>n;
            ss>>val;
            
            D(i-1,j-1,m-1,k-1,l-1,n-1)=val;
            //std::cout<<"D_"<<i<<j<<m<<k<<l<<n<<"="<<D(i-1,j-1,m-1,k-1,l-1,n-1)<<std::endl;
            numRead++;
        }
        assert(numRead==729 && "WRONG NUMBER OF LINES IN D-file"); 
    }
    else{
        std::cout<<"FILE "<<filename<<" cannot be opened"<<std::endl;
        assert(0);
    }
}

inline void readBurgers(const std::string& dirName,double& bs)
{

    
    std::string filename=dirName+"/Burgers.out";
    std::ifstream ifs3 ( filename.c_str() , std::ifstream::in );

    if (ifs3.is_open())
    {
        std::cout<<"Reading: "<<filename<<std::endl;
        std::string line;
        while (std::getline(ifs3, line))
        {
            std::stringstream ss(line);

            ss>>bs;
            std::cout<<"Burgers="<<bs<<std::endl;
        }
    }
    else
    {
        std::cout<<"FILE "<<filename<<" cannot be opened"<<std::endl;
        assert(0);
    }
    
    

    
    
}


inline double leviCivita(const double& i,
                  const double& j,
                  const double& k)
{
    return 0.5*(i-j)*(j-k)*(k-i);
}


inline Tensor<3,3> cbedl(const Tensor<4,3>& C,
                  const Eigen::Matrix<double,3,1>& b,
                  const Eigen::Matrix<double,3,1>& dL)
{
    /* computes the tensor
     * cbedl(m,n,l)=C(m,n,p,q)*b(p)*leviCivita(q,r,l)*dL(r)
     */
    
    Tensor<3,3> Cb;
    for (size_t m=0;m<3;++m){
        for (size_t n=0;n<3;++n){
            for (size_t q=0;q<3;++q){
                for (size_t p=0;p<3;++p){
                    Cb(m,n,q)+=C(m,n,p,q)*b(p);
                }
            }
        }
    }
    
    Tensor<2,3> edL;
    for (size_t q=0;q<3;++q){
        for (size_t l=0;l<3;++l){
            for (size_t r=0;r<3;++r){
                edL(q,l)+=leviCivita(q,r,l)*dL(r);
            }
        }
    }
    
    Tensor<3,3> CbedL;
    for (size_t m=0;m<3;++m){
        for (size_t n=0;n<3;++n){
            for (size_t l=0;l<3;++l){
                for (size_t q=0;q<3;++q){
                    CbedL(m,n,l)+=Cb(m,n,q)*edL(q,l);
                }
            }
        }
    }
    return CbedL;
}


#endif  // _ANISOTROPIC_HEADERS_H