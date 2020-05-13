#ifndef INCLUDED_LOAD_2D_ISING
#define INCLUDED_LOAD_2D_ISING

#include <iostream>

int loadConfig(FILE* config, double* beta, int* Nstep, char* state_file, char* accept_file, double* J, char* initial_state);

#endif

