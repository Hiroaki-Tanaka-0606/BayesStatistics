#include "LorentzDistribution.cpp"
#include <iostream>

// Test program for lorentz distribution

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
	LD* ld=new LD(sigma,mu);
	for(i=0;i<imax;i++){
		cout << ld->rand() << endl;
	}
}

