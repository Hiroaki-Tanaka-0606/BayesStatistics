#include <iostream>
#include "../MCMC/Model_selection_mND.hpp"
#include <string.h>
#include "loadConfig_mND.hpp"

using namespace std;

int main(int argc, char** argv){
  cout << "Generate multiple RXMC configuration files" << endl;
  if(argc<2){
    printf("Usage: %s config_file\n", argv[0]);
    return 0;
  }

  
  State* s=new State();
  char* initial_state=s->print();
  int num_params=s->num_params;

  // load configuration
  char* config_file=argv[1];
  cout << "Configuration file: " << config_file << endl;
  FILE* config=fopen(config_file, "r");
  if(config==NULL){
    cout << "Error in opening the config file" << endl;
    return -1;
  }
  // configuration for RXMC.o
  int Nchains; // number of Marcov chains
  int Nreplica;
  double* beta; // inverse temperature
  int Nstep; // Number of MC steps
  int buffer_length=256;
  char* input_file_format=new char[buffer_length]; // name of the input file, including %d (chain index)
  char* state_file_format=new char[buffer_length]; // name of the output file for development of the state, including %d (chain index), %%d (replica index)
  char* accept_file_format=new char[buffer_length]; // name of the output file for log of acceptance, including %d (chain index), %%d (replica index)
  char* exchange_file_format=new char[buffer_length]; // name of the output file for log of exchange, including %d (chain index)
  char* log_file_format=new char[buffer_length]; // name of the log (STDOUT), including %d (chain index)
  char* sample_file=new char[buffer_length]; // sample file
  int Nbin; // number of bins
  double** sigma; // development parameter
  char* shell_script=new char[buffer_length]; // name of the shell script file

  // configuration for Average_analysis.o
  // number of parameters, buffer size for fgets, state file[i], Nstep, are required
  char* input2_file_format=new char[buffer_length]; // name of the input, including %%d (chain index) and %d (replica index)
  int Burnin; // Burn-in steps
  char* average_file_format=new char[buffer_length]; // name of the output file for average, including %%d (chain index) and %d (replica index)
  char* variance_file_format=new char[buffer_length]; // name of the output file for variance, including  %%d (chain index) and %d (replica index)
  char* log2_file_format=new char[buffer_length]; // name of the log, including  %%d (chain index) and %d (replica index)
  char* value_format=new char[buffer_length]; // output format for a value
  
  // configuration for PSRF.o
  // number of parameters, buffer size for fgets, number of MCs, , average file format, variance file format, Nstep-Burnin, are required
  char* input3_file_format=new char[buffer_length]; // name of the input, including %d (replica index)
  char* psrf_file_format=new char[buffer_length]; // name of the output for psrf, including %d (replica index) and %%d (parameter index, not chain index)
  char* log3_file_format=new char[buffer_length]; // name of the log, including %d (replica index)

  // pointer for loadConfig arguments
  double** beta_pointer=new double*[1];
  double*** sigma_pointer=new double**[1];
  
  int loadConfig_status=loadConfig(config, &Nchains, &Nreplica, beta_pointer, &Nstep, input_file_format, state_file_format, accept_file_format, exchange_file_format, log_file_format, sample_file, &Nbin, sigma_pointer, num_params, shell_script, input2_file_format, &Burnin, average_file_format, variance_file_format, log2_file_format, value_format, input3_file_format, psrf_file_format, log3_file_format);
  if(loadConfig_status!=1){
    cout << "Failed in loading configuration" << endl;
    return -1;
  }

  beta=beta_pointer[0];
  sigma=sigma_pointer[0];
  int i,j;
  
  cout << "Number of Marcov chains: " << Nchains << endl;
  cout << "Number of replicas: " << Nreplica << endl;
  cout << "Inverse temperatures: ";
  for(i=0;i<Nreplica;i++){
    cout << beta[i] << ", ";
  }
  cout << endl;
  cout << "Number of MC steps: " << Nstep << endl;
  cout << "Sample file: " << sample_file << endl;
  cout << "Shell script file: " << shell_script << endl;

  // generate shell script file
  FILE* shell=fopen(shell_script, "w");
  if(shell==NULL){
    cout << "Error in opening the shell script file" << endl;
    return -1;
  }
  fprintf(shell, "echo 'Start RXMC.o and Average_analysis.o (Total %d)'\n", Nchains);
  
  FILE* input; // for RXMC.o
  char* input_file=new char[buffer_length];
  char* input_buffer=new char[buffer_length];

  // i -> chain index, j -> replica index
  for(i=0;i<Nchains;i++){
    // RXMC.o
    // open input file
    sprintf(input_file, input_file_format, i+1);
    input=fopen(input_file, "w");
    if(input==NULL){
      cout << "Error in opening the input file" << endl;
      return -1;
    }
    // line 1: Nreplica
    fprintf(input, "%d # Nreplica\n", Nreplica);
    // line 2-(Nreplica-1): beta
    for(j=0;j<Nreplica;j++){
      fprintf(input, "%.4e # beta for replica %d\n", beta[j], j+1);
    }
    // line 2+Nreplica: Nstep
    fprintf(input, "%d # Nstep\n", Nstep);
    // line 3+Nreplica: state file
    fprintf(input, state_file_format, i+1);
    fprintf(input, " # state file\n");
    // line 4+Nreplica: accept file
    fprintf(input, accept_file_format, i+1);
    fprintf(input, " # accept file\n");
    // line 5+Nreplica: exchange file
    fprintf(input, exchange_file_format, i+1);
    fprintf(input, " # exchange file\n");
    // line 6+Nreplica: sample_file
    fprintf(input, "%s # sample file\n", sample_file);
    // line 7+Nreplica: Nbin
    fprintf(input, "%d # Nbin\n", Nbin);
    // lin (8+Nreplica)-(7+Nreplica): development parameters
    for(j=0;j<Nreplica;j++){
      int k;
      for(k=0;k<num_params;k++){
	fprintf(input, "%12.4e ", sigma[j][k]);
      }
      fprintf(input, "# development parameters for replica %d\n",j+1);
    }
    // line (8+2*Nreplica)-(7+3*Nreplica): initial state
    for(j=0;j<Nreplica; j++){
      s->init_rand();
      delete initial_state;
      initial_state=s->print();
      fprintf(input, "%s # initial state for replica %d\n", initial_state, i+1);
    }
    
    fclose(input);

    // add command to shell file
    fprintf(shell, "./RXMC_mND.o %s", input_file);
    fprintf(shell, " > ");
    fprintf(shell, log_file_format, i+1);
    fprintf(shell, "\n");
    fprintf(shell, "echo 'Done RXMC_mND.o for chain %d'\n", i+1);
    
    // Average_analysis.o
    for(j=0;j<Nreplica;j++){
      sprintf(input_buffer, input2_file_format, j+1);
      sprintf(input_file, input_buffer, i+1);
      input=fopen(input_file, "w");
      if(input==NULL){
	cout << "Error in opening the input2 file" << endl;
	return -1;
      }
      // line 1: Number of parameters
      fprintf(input, "%d # Number of parameters\n", s->num_params);
      // line 2: buffer size (length of print() or StateHeader())
      char* header=StateHeader();
      int data_length=max(strlen(header), strlen(initial_state))+10; // +10 for \r, \n, \0
      fprintf(input, "%d # buffer size\n", data_length);
      // line 3: input (state file)
      sprintf(input_buffer, state_file_format, i+1);
      fprintf(input, input_buffer, j+1);
      fprintf(input, " # input file (state)\n");
      // line 4: Number of MC steps
      fprintf(input, "%d # Number of MC steps\n", Nstep);
      // line 5: Burnin
      fprintf(input, "%d # Burnin\n", Burnin);
      // line 6: Average file
      sprintf(input_buffer, average_file_format, j+1);
      fprintf(input, input_buffer, i+1);
      fprintf(input, " # Average file\n");
      // line 7: Variance file
      sprintf(input_buffer, variance_file_format, j+1);
      fprintf(input, input_buffer, i+1);
      fprintf(input, " # Variance file\n");
      // line 8: Value format
      fprintf(input, "%s # Value format\n", value_format);
    
      fclose(input);

      // add command
      fprintf(shell, "./Average_analysis.o %s", input_file);
      fprintf(shell, " > ");
      sprintf(input_buffer, log2_file_format, j+1);
      fprintf(shell, input_buffer, i+1);
      fprintf(shell, "\n");
      fprintf(shell, "echo 'Done Average_analysis.o for chain %d replica %d'\n", i+1, j+1);
    }
  }

  for(j=0;j<Nreplica;j++){
    // PSRF.o
    sprintf(input_file, input3_file_format, j+1);
    input=fopen(input_file,"w");
    if(input==NULL){
      cout << "Error in opening the input3 file" << endl;
      return -1;
    }
    // line 1: Number of parameters
    fprintf(input, "%d # Number of parameters\n", s->num_params);
    // line 2: Buffer size
    char* value_print=new char[buffer_length];
    sprintf(value_print, value_format, 0.0);
    fprintf(input, "%d # Buffer size\n", (s->num_params)*(1+strlen(value_print)));
    // line 3: Number of Marcov chains
    fprintf(input, "%d # Number of Marcov Chains\n", Nchains);
    // line 4: Average file format
    sprintf(input_buffer, average_file_format, j+1);
    fprintf(input, "%s # average file format\n", input_buffer);
    // line 5: variance file format
    sprintf(input_buffer, variance_file_format, j+1);
    fprintf(input, "%s # variance file format\n", input_buffer);
    // line 6: Nstep-Burnin
    fprintf(input, "%d # Nstep - Burnin\n", Nstep-Burnin);
    // line 7: output file format
    sprintf(input_buffer, psrf_file_format, j+1);
    fprintf(input, "%s # output file format\n", input_buffer);

    fclose(input);

    // add command

    fprintf(shell, "./PSRF.o %s", input_file);
    fprintf(shell, " > ");
    fprintf(shell, log3_file_format, j+1);
    fprintf(shell, "\n");
    fprintf(shell, "echo 'Done PSRF.o for replica %d'\n", j+1);

  }
  fclose(shell);
  return 1;
}
