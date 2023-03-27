
#ifndef  mmdl_Voigt2Tensor_H_
#define mmdl_Voigt2Tensor_H_

template <short unsigned int rank, short unsigned int dim=3>
class Voigt2Tensor{
	
	
	};
	
	
	
	
template <>
class Voigt2Tensor<4,3>{
	enum{dim=3};
	  Eigen::Matrix<double, dim*(dim+1)/2,dim*(dim+1)/2>& V	;
	
	
	int voigthIndex(const int& i, const int&j) const {
		
		if (j<i){
			return voigthIndex(j,i);
		}
		else if (j==i){
			return i;
		}
		else if (i==1 && j==2){
			return 3;
		}
		else if (i==0 && j==2){
			return 4;
		}
		else if (i==0 && j==1){
			return 5;
		}		
		else{
			std::cout << "Voigt2Tensor.h, voigthIndex function having trouble" << std::endl;
			return 10000;
		}
	}
	
	public:
		
	Voigt2Tensor( Eigen::Matrix<double, dim*(dim+1)/2,dim*(dim+1)/2>& V_in) : V(V_in) {}
		
	
	double& operator()(const int& i,const int& j,const int& k,const int& l) {
		return V(voigthIndex(i,j),voigthIndex(k,l));
		}
		
	const double& operator()(const int& i,const int& j,const int& k,const int& l) const {
		return V(voigthIndex(i,j),voigthIndex(k,l));
		}
		
	};


#endif