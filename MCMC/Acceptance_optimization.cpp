#include <iostream>
#include "Acceptance_optimization.hpp"
#include <cmath>

using namespace std;

// f(sigma): acceptance ratio
// F(sigma)=accept/trial: simulated acceptance ratio
// considering error sqrt(trial), determine whether f(x(t)) is larger than target or not
// F(x(t))+error_coeff*error < target -> f(x(t)) is smaller than target
// F(x(t))-error_coeff*error > target -> f(x(t)) is larger than target
// error=1/sqrt(trial)
// error_coeff: coefficient (input)
// development of arguemnt sigma by bisection
// stop when error < error_stop (input)
int Acceptance_optimization(char* index_print, double* sigma_small, double* sigma, double* sigma_large, int* stage, int* accept, int* trial, double target, double error_coeff, double error_stop){

  double Fsigma=(*accept)*1.0/(*trial);
  double error=1.0/sqrt(*trial);
  switch(*stage){
  case 0:
    // find sigma s.t. F(sigma) < target < x(0)=1
    if(Fsigma + error_coeff*error < target){
      // sigma is found, start bisection
      cout << "Bisection optimizaztion start for parameter " << index_print << endl;
      cout << *accept << " / " << *trial << " +- " << error << " < " << target << endl;
      *sigma_small=*sigma;
      // sigma_large is determined in the previous step ( or initial value 0)
      *sigma=(*sigma_small+*sigma_large)/2.0;
      *stage=1;
      *accept=0;
      *trial=0;
      return 1;
    }
    if(Fsigma - error_coeff*error > target){
      // sigma is found to be too small to satisfy F(sigma) < target
      cout << "Update parameter " << index_print << endl;
      cout << *accept << " / " << *trial << " +- " << error << " > " << target << endl;
      *sigma_large=*sigma;
      *sigma*=2;
      *accept=0;
      *trial=0;
      return 1;
    }
    break;
  case 1:
    // bisection optimization
    if(Fsigma + error_coeff*error < target){
      // F(sigma_small), F(sigma) < target < F(sigma_large)
      cout << "Bisection update for parameter " << index_print << endl;
      cout << *accept << " / " << *trial << " +- " << error << " < " << target << endl;
      *sigma_small=*sigma;
      *sigma=(*sigma_small+*sigma_large)/2.0;
      *accept=0;
      *trial=0;
      return 1;
    }
    if(Fsigma - error_coeff*error > target){
      // F(sigma_small) < target < F(sigma), F(sigma_large)
      cout << "Bisection update for parameter " << index_print << endl;
      cout << *accept << " / " << *trial << " +- " << error << " > " << target << endl;
      *sigma_large=*sigma;
      *sigma=(*sigma_small+*sigma_large)/2.0;
      *accept=0;
      *trial=0;
      return 1;
    }
    if(error<error_stop){
      // F(sigma) ~ target
      cout << "Bisection optimization ended for parameter " << index_print << endl;
      cout << *accept << " / " << *trial << " +- " << error << " ~ " << target << endl;
      *stage=2;
      return 0;
    }
    break;
  case 2:
    // optimized
    break;
  }
  return 0;
}
