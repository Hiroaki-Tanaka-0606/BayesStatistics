//bisection.cpp: optimize parameter of probability distribution by bisection method
//probability variable: return 1 by probability f(x), 0 by 1-f(x)
//f(0)=1, f(infty)=0, x=[0, infty)
//target=[0,1] (input)
//considering error sqrt(N), determine whether f(x(t)) is larger than target or not
//F(x(t))+n*sigma < target -> f(x(t)) is smaller than target
//F(x(t))-n*sigma > target -> f(x(t)) is larger than target
//F(x(t))=ratio of 1=(number of 1)/(total number: N)
//sigma=1/sqrt(total number)=1/sqrt(N)
//n: coefficient (input)
//x(0) (input), developed by bisection method
//stop when sigma < sigma_stop (input)

#include <iostream>
#include "../randomDistribution/mersenne_twister.hpp"
#include <cmath>

//probability f(x)
//f(0)=1, x>=0, f(infty)=0
double f(double x){
  return exp(-x*x);
}

//probability trial until F(x(t)) is determined to be larger or smaller than target
int prob_trial(double x, double target, double n, double sigma_stop, MT* mt1, int* N1, int* N, double* ratio, double* sigma);

//argument: target, n, x(0), sigma_stop
int main(int argc, char** argv){
  if(argc<5){
    printf("Usage: %s target n x(0) sigma_stop\n", argv[0]);
    return 0;
  }
  double target=atof(argv[1]);
  double n=atof(argv[2]);
  double x0=atof(argv[3]);
  double sigma_stop=atof(argv[4]);

  //probability parameter
  double x=x0;
  //current step
  int t=0;
  //total number of trial (in the current step)
  int N=0;
  //number of 1 (in the current step)
  int N1=0;
  //ratio of 1 (negative value means undefined because N=0)
  double ratio=-1;
  //sigma (negative value means undefined because N=0)
  double sigma=-1;  

  cout << "# Optimize paramter of probability by bisection method" << endl;
  cout << "# target: " << target << endl;
  cout << "# sigma_stop: " << sigma_stop << endl;
  cout << "# sigma coefficient: " << n << endl;
  cout << "# t x N1 N ratio sigma" << endl;
  MT* mt1=new MT();

  //first step: find x(t) s.t. x(t) < target < x(0)=1
  int flag=1; // 1 when F(x(t)) > target, -1 when F(x(t)) < target
  while(flag>0){
    flag=prob_trial(x, target, n, sigma_stop, mt1, &N1, &N, &ratio, &sigma);
    printf("%d %f %d %d %f %f\n", t, x, N1, N, ratio, sigma);
    if(flag<=0){
      break;
    }
    x*=2;
  }
  double xLarge=x/2; // F(xLarge) > target
  double xSmall=x; // F(xSmall) < target
  x=(xLarge+xSmall)/2; //bisection
  while(sigma>sigma_stop){
    flag=prob_trial(x, target, n, sigma_stop, mt1, &N1, &N, &ratio, &sigma);
    printf("%d %f %d %d %f %f\n", t, x, N1, N, ratio, sigma);
    if(flag>0){
      //F(xSmall) < target < F(x), F(xLarge);
      xLarge=x;
    }else if(flag<0){
      //F(xSmall), F(x) < target < F(xLarge);
      xSmall=x;
    }else{
      //F(x) ~ target
      break;
    }
    x=(xLarge+xSmall)/2;
  }
}

//probability trial until F(x(t)) is determined to be larger or smaller than target
int prob_trial(double x, double target, double n, double sigma_stop, MT* mt1, int* N1, int* N, double* ratio, double* sigma){
  //initial value
  *N=0;
  *N1=0;
  *ratio=-1;
  *sigma=-1;
  //reference by which probability variable (1 or 0) is determined
  double ref=f(x);
  while((*ratio)<0 || !((*ratio)+(*sigma)*n < target || (*ratio)-(*sigma)*n > target || (*sigma) < sigma_stop)){
    double rand=mt1->rand();
    (*N)++;
    if(rand<ref){
      (*N1)++;
    }
    *ratio=(*N1)*1.0/(*N);
    *sigma=1.0/sqrt(*N);
  }
  if((*ratio)+(*sigma)*n < target){
    //F(x(t)) < target
    return -1;
  }else if((*ratio)-(*sigma)*n > target){
    //F(x(t)) > target
    return 1;
  }else{
    //F(x(t)) ~ target
    return 0;
  }
}
