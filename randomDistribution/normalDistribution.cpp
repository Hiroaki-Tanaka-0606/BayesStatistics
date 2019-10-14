#define _USE_MATH_DEFINES
#include "mersenne_twister.cpp"
#include <random>
#include <cmath>
#include <iostream>

class ND
{
public:
	ND(){
		ND::mt1=new MT();
		ND::mt2=new MT();
		ND::sigma=1;
		ND::mu=0;
	}
	ND(double sigma, double mu){
		ND::mt1=new MT();
		ND::mt2=new MT();
		ND::sigma=sigma;
		ND::mu=mu;
	}
	ND(int seed1, int seed2, double sigma, double mu){
		ND::mt1=new MT(seed1);
		ND::mt2=new MT(seed2);
		ND::sigma=sigma;
		ND::mu=mu;
	}
	double rand(){
		double rand1=mt1->rand();
		double rand2=mt2->rand();
		double r=sqrt(-2*log(1-rand1));
		double theta=2*M_PI*rand2;
		return sigma*(r*cos(theta))+mu;
		
	}
private:
	MT* mt1;
	MT* mt2;
	double sigma;
	double mu;
};
	
/*
int main(){
	ND nd;
	int i;
	int imax=100000;
	for(i=0;i<imax;i++){
		cout << nd.rand() <<endl;
	}
}
*/
