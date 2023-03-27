#ifndef _readBurgers_h_
#define _readBurgers_h_

void readBurgers(const std::string& dirName,double& bs)
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

#endif  // _readBurgers_h_
