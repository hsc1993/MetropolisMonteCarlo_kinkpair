#ifndef _READ_CL_H
#define _READ_CL_H 


void readCL(const std::string& filename, Eigen::Matrix<double,6,6>& C, Eigen::Matrix<double,3,3>& L, int& nPhi, int& nQ)
{
    
    double temp;
    
    std::string line;
    std::ifstream ifs ( filename.c_str() , std::ifstream::in );
    
    if (ifs.is_open())
    {
        std::cout<<"Reading: "<<filename<<std::endl;
        
        
        for (int i=0;i<6;++i)
        {
            for (int j=0;j<6;++j)
            {
                std::getline(ifs, line);
                std::stringstream ss(line);
                ss >> temp;
                C(i,j)=temp;
            }
        }
        
        for (int i=0;i<3;++i)
        {
            for (int j=0;j<3;++j)
            {
                std::getline(ifs, line);
                std::stringstream ss(line);
                ss >> temp;
                L(i,j)=temp;
            }
        }

        {
        std::getline(ifs, line);
        std::stringstream ss(line);
        ss >> temp;
        nPhi=temp;
        }
        
        {
        std::getline(ifs, line);
        std::stringstream ss(line);
        ss >> temp;
        nQ=temp;
        }
        
        std::cout<<"Matrix C=\n"<<C<<std::endl;
        std::cout<<"L=\n"<<L<<std::endl;
        std::cout<<"nPhi="<<nPhi<<std::endl;
        std::cout<<"nQ="<<nQ<<std::endl;
        
    }
    else
    {
        std::cout<<"FILE "<<filename<<" cannot be opened"<<std::endl;
        assert(0);
    }
    
}


#endif  // _READ_CL_H