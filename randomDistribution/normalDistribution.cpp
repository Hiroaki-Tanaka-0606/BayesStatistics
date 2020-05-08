#define _USE_MATH_DEFINES
#include "normalDistribution.hpp"
#include <cmath>

// class ND: generate normal distribution(sigma,mu) by Box-Muller method

ND::ND(){
  ND::mt1=new MT();
  ND::mt2=new MT();
  ND::sigma=1;
  ND::mu=0;
}
ND::ND(double sigma, double mu){
  ND::mt1=new MT();
  ND::mt2=new MT();
  ND::sigma=sigma;
  ND::mu=mu;
}
ND::ND(int seed1, int seed2, double sigma, double mu){
  ND::mt1=new MT(seed1);
  ND::mt2=new MT(seed2);
  ND::sigma=sigma;
  ND::mu=mu;
}

double ND::rand(){
  double rand1=mt1->rand();
  double rand2=mt2->rand();
  double r=sqrt(-2*log(1-rand1));
  double theta=2*M_PI*rand2;
  return sigma*(r*cos(theta))+mu;
		
}
