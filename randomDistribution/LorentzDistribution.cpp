#include "LorentzDistribution.hpp"
#define _USE_MATH_DEFINES
#include <cmath>


//class LD: generate Lorentz distribution(sigma, mu)
LD::LD(){
  LD::mt=new MT();
  LD::sigma=1;
  LD::mu=0;
}
LD::LD(double sigma, double mu){
  LD::mt=new MT();
  LD::sigma=sigma;
  LD::mu=mu;
}
LD::LD(int seed, double sigma, double mu){
  LD::mt=new MT(seed);
  LD::sigma=sigma;
  LD::mu=mu;
}

double LD::rand(){
  double rand=mt->rand();
  return sigma*(tan(M_PI*(rand-0.5)))+mu;
		
}
