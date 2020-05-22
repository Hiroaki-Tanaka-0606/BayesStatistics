#ifndef INCLUDED_LOAD_MND
#define INCLUDED_LOAD_MND

#include <iostream>
int loadConfig_RXMC(FILE* config, int* Nreplica, double** beta_pointer, int* Nstep, char* state_file_format, char* accept_file_format, char* exchange_file_name, char* sample_file, int* Nbin, double*** sigma_pointer, int num_params, int state_buffer_length, char*** initial_state_pointer);

#endif
