#ifndef _READ_C_TENSOR
#define _READ_C_TENSOR


void readCtensor(const std::string& dirName,
                 Tensor<4,3>& C)
{
    int i,j,k,l;
    double val;
    
    std::string filename=dirName+"/Cijkl_ext.out";
    std::ifstream ifs ( filename.c_str() , std::ifstream::in );
    
    int numRead=0;
    
    if (ifs.is_open())
    {
        std::cout<<"Reading: "<<filename<<std::endl;
        
        std::string line;
        while (std::getline(ifs, line))
        {
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
    else
    {
        std::cout<<"FILE "<<filename<<" cannot be opened"<<std::endl;
        assert(0);
    }
    
}


#endif // _READ_C_TENSOR