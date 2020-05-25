#ifndef INCLUDE_LOAD_MND
#define INCLUDE_LOAD_mND

#include <iostream>
int loadConfig(FILE* config, int* Nchains, int* Nreplica, double** beta_pointer, int* Nstep, char* input_file_format, char* state_file_format, char* accept_file_format, char* exchange_file_format, char* log_file_format, char* sample_file, int* Nbin, double*** sigma_pointer, int num_params, char* shell_script, char* input2_file_format, int* Burnin, char* average_file_format, char* variance_file_format, char* log2_file_format, char* value_format, char* input3_file_format, char* psrf_file_format, char* log3_file_format);


#endif
