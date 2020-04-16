#define _USE_MATH_DEFINES
#include "mersenne_twister.cpp"
#include <random>
#include <cmath>

class LD
{
public:
	LD(){
		LD::mt=new MT();
		LD::sigma=1;
		LD::mu=0;
	}
	LD(double sigma, double mu){
		LD::mt=new MT();
		LD::sigma=sigma;
		LD::mu=mu;
	}
	LD(int seed, double sigma, double mu){
		LD::mt=new MT(seed);
		LD::sigma=sigma;
		LD::mu=mu;
	}
	double rand(){
		double rand=mt->rand();
		return sigma*(tan(M_PI*(rand-0.5)))+mu;
		
	}
private:
	MT* mt;
	double sigma;
	double mu;
};