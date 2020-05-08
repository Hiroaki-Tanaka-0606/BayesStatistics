#include "normalDistribution.hpp"
#include <iostream>

// Test program for normal distribution

int main(int argc, char** argv){
	int i;
	double sigma, mu;
	int imax;
	if(argc>=4){
		sigma=atof(argv[1]);
		mu=atof(argv[2]);
		imax=atoi(argv[3]);
	}else{
		sigma=1.0;
		mu=0.0;
		if(argc>=2){
			imax=atoi(argv[1]);
		}else{
			imax=100000;
		}
	}
	cout << "# sigma = " << sigma << endl;
	cout << "# mu = " << mu << endl;
	cout << "# Number of data = " << imax << endl;
	ND* nd=new ND(sigma,mu);
	for(i=0;i<imax;i++){
		cout << nd->rand() << endl;
	}
}

