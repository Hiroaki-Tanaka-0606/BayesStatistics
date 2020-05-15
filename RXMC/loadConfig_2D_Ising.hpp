#ifndef INCLUDED_LOAD_2D_ISING
#define INCLUDED_LOAD_2D_ISING

#include <iostream>

int loadConfig_RXMC(FILE* config, int* Nreplica, double** beta_pointer, int* Nstep, char* state_file_format, char* accept_file_format, char* exchange_file_name, double* J, int state_buffer_length, char*** initial_state_pointer);

#endif

