#include "MC_Metropolis.hpp"
#include <cmath>

// update s[index] by Metropolis method
int MC_Metropolis(State* s, State* s_next, Hamiltonian* h, int index, double sigma, double beta, MT* mt){
  
  s->next(index, sigma, s_next);
  double e1=h->energy(s);
  double e2=h->energy(s_next);

  double threshold=exp(-beta*(e2-e1));
  double rand=mt->rand();
  if(rand<threshold){
    // accept
    return 1;
  }else{
    // reject
    return 0;
  }
}
