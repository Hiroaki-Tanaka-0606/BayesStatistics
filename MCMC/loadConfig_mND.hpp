#ifndef INCLUDED_LOAD_2D_ISING
#define INCLUDED_LOAD_2D_ISING

#include <iostream>

int loadConfig(FILE* config, double* beta, int* Nstep, char* state_file, char* accept_file, char* sample_file, int* Nbin, char* initial_state);

#endif

