#include "MC_Metropolis.hpp"
#include <cmath>

// update s[index] by Metropolis method
int MC_Metropolis(State* s, State* s_next, Hamiltonian* h, int index, double sigma, double beta, MT* mt){
  
  s->next(index, sigma, s_next);
  double bh1=h->betaH(s,beta);
  double bh2=h->betaH(s_next,beta);
  if(!isfinite(bh2)){
    return 0;
  }

  double threshold=exp(-(bh2-bh1));
  double rand=mt->rand();
  if(rand<threshold){
    // accept
    return 1;
  }else{
    // reject
    return 0;
  }
}
