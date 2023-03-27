#ifndef _READ_XY_H
#define _READ_XY_H


#include <deque>

template<int cols>
std::deque<Eigen::Matrix<double,cols,1> > readXY(const std::string& filename, const Eigen::Matrix<double,3,1>& center=Eigen::Matrix<double,3,1>::Zero())
{
    
    double temp;
//    std::string filename="XY.txt";
    
    std::deque<Eigen::Matrix<double,cols,1> > deq;
    
    std::string line;
    std::ifstream ifs ( filename.c_str() , std::ifstream::in );
    std::cout<<"filename= "<<filename<<std::endl;

    
    if (ifs.is_open())
    {
        std::cout<<"Reading: "<<filename<<std::endl;
        
        while (std::getline(ifs, line)) {
            std::stringstream ss(line);
            
            int col=0;
            
            Eigen::Matrix<double,cols,1> tempV(Eigen::Matrix<double,cols,1>::Zero());
            
            while (ss >> temp) {
                
                tempV(col)=temp;
                
                
                col++;
                
                
            }
            
//            std::cout<<tempV.transpose()<<std::endl;
            tempV.template block<3,1>(0,0)-=center;
//            tempV(0)-=center(0);
//            tempV(1)-=center(1);
//            tempV(2)-=center(2);
            
            deq.push_back(tempV);
        }
    }
    else
    {
        assert(0 && "UNABLE TO OPEN FILE.");
    }
    
    return deq;
}


#endif  // _READ_XY_H


