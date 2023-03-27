#ifndef _READ_D_TENSOR_H
#define _READ_D_TENSOR_H



void readDtensor(const std::string& dirName,
                 Tensor<6,3>& D)
{
    
    int i,j,k,l,m,n;
    double val;
    
    
    int numRead=0;
    
    
    
    
    
    std::string filename=dirName+"/Dijmkln_ext.out";
    std::ifstream ifs2 ( filename.c_str() , std::ifstream::in );
    
    if (ifs2.is_open())
    {
        std::cout<<"Reading: "<<filename<<std::endl;
        std::string line;
        while (std::getline(ifs2, line))
        {
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
    else
    {
        std::cout<<"FILE "<<filename<<" cannot be opened"<<std::endl;
        assert(0);
    }
    
    //    filename=dirName+"/Burgers.out";
    //    std::ifstream ifs3 ( filename.c_str() , std::ifstream::in );
    //
    //    if (ifs3.is_open())
    //    {
    //        std::cout<<"Reading: "<<filename<<std::endl;
    //        std::string line;
    //        while (std::getline(ifs3, line))
    //        {
    //            std::stringstream ss(line);
    //
    //            ss>>bs;
    //            std::cout<<"Burgers="<<bs<<std::endl;
    //        }
    //    }
    //    else
    //    {
    //        std::cout<<"FILE "<<filename<<" cannot be opened"<<std::endl;
    //        assert(0);
    //    }
    
    
    //    std::ofstream dcl(dirName + "/DvsCL.txt");
    //
    //    ElasticL<Cubic> L(C,D);
    //    for(int i=0;i<3;++i)
    //        for(int j=0;j<3;++j)
    //            for(int k=0;k<3;++k)
    //                for(int l=0;l<3;++l)
    //                    for(int m=0;m<3;++m)
    //                        for(int n=0;n<3;++n)
    //                {
    //                    dcl<<i<<j<<m<<k<<l<<n<<" "<<D(i,j,m,k,l,n)<<" "<<C(i,j,k,l)*L(m,n)<<"\n";
    //                }
    
    
}


#endif  // _READ_D_TENSOR_H