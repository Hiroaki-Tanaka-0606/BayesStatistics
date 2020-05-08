#ifndef INCLUDED_MC_METROPOLIS
#define INCLUDED_MC_METROPOLIS

#include "Model_selection.hpp"
#include "../randomDistribution/mersenne_twister.hpp"

// update s[index] by Metropolis method
int MC_Metropolis(State* s, State* s_next, Hamiltonian* h, int index, double sigma, double beta, MT* mt);

#endif
