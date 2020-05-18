#include <iostream>
#include "../MCMC/Model_selection.hpp"
#include <string.h>

using namespace std;

int loadConfig(FILE* config, int* Nchains, int* Nreplica, double** beta_pointer, int* Nstep, char* input_file_format, char* state_file_format, char* accept_file_format, char* exchange_file_format, char* log_file_format, double* J, char* shell_script, char* input2_file_format, int* Burnin, char* average_file_format, char* variance_file_format, char* log2_file_format, char* value_format, char* input3_file_format, char* psrf_file_format, char* log3_file_format);


int main(int argc, char** argv){
  cout << "Generate multiple RXMC configuration files" << endl;
  if(argc<2){
    printf("Usage: %s config_file\n", argv[0]);
    return 0;
  }

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
  double J; //bond strength
  char* shell_script=new char[buffer_length]; // name of the shell script file

  // configuration for Average_analysis.o
  // number of parameters, buffer size for fgets, state file[i], Nstep, are required
  char* input2_file_format=new char[buffer_length]; // name of the input 
  int Burnin; // Burn-in steps
  char* average_file_format=new char[buffer_length]; // name of the output file for average
  char* variance_file_format=new char[buffer_length]; // name of the output file for variance
  char* log2_file_format=new char[buffer_length]; // name of the log
  char* value_format=new char[buffer_length]; // output format for a value
  
  // configuration for PSRF.o
  // number of parameters, buffer size for fgets, number of MCs, , average file format, variance file format, Nstep-Burnin, are required
  char* input3_file_format=new char[buffer_length]; // name of the input
  char* psrf_file_format=new char[buffer_length]; // name of the output for psrf, including %d (parameter index, not chain index)
  char* log3_file_format=new char[buffer_length]; // name of the log

  // pointer for loadConfig arguments
  double** beta_pointer=new double*[1];
  
  int loadConfig_status=loadConfig(config, &Nchains, &Nreplica, beta_pointer, &Nstep, input_file_format, state_file_format, accept_file_format, exchange_file_format, log_file_format, &J, shell_script, input2_file_format, &Burnin, average_file_format, variance_file_format, log2_file_format, value_format, input3_file_format, psrf_file_format, log3_file_format);
  if(loadConfig_status!=1){
    cout << "Failed in loading configuration" << endl;
    return -1;
  }

  beta=beta_pointer[0];
  int i,j;
  
  cout << "Number of Marcov chains: " << Nchains << endl;
  cout << "Number of replicas: " << Nreplica << endl;
  cout << "Inverse temperatures: ";
  for(i=0;i<Nreplica;i++){
    cout << beta[i] << ", ";
  }
  cout << endl;
  cout << "Number of MC steps: " << Nstep << endl;
  cout << "J: " << J << endl;
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
  State* s=new State();
  char* initial_state=s->print();

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
    // line 6+Nreplica: J
    fprintf(input, "%.4e # J\n", J);
    // line (7+Nreplica)-(6+2*Nreplica): initial state
    for(j=0;j<Nreplica; j++){
      s->init_rand();
      delete initial_state;
      initial_state=s->print();
      fprintf(input, "%s # initial state for replica %d\n", initial_state, i+1);
    }
    
    fclose(input);

    // add command to shell file
    fprintf(shell, "./RXMC.o %s", input_file);
    fprintf(shell, " > ");
    fprintf(shell, log_file_format, i+1);
    fprintf(shell, "\n");
    fprintf(shell, "echo 'Done RXMC.o for chain %d'\n", i+1);
    
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
      // line 2: buffer size (length of print())
      fprintf(input, "%d # buffer size\n", strlen(initial_state));
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

int loadConfig(FILE* config, int* Nchains, int* Nreplica, double** beta_pointer, int* Nstep, char* input_file_format, char* state_file_format, char* accept_file_format, char* exchange_file_format, char* log_file_format, double* J, char* shell_script, char* input2_file_format, int* Burnin, char* average_file_format, char* variance_file_format, char* log2_file_format, char* value_format, char* input3_file_format, char* psrf_file_format, char* log3_file_format){
  int buffer_length=1024;
  char line[buffer_length];
  char* fgets_status;
  int sscanf_status;

  // line 1: Number of chains
  fgets_status=fgets(line, buffer_length, config);
  if(fgets_status==NULL){
    cout << "Error in loading line 1 of the configuration file" << endl;
    return 0;
  }
  sscanf_status=sscanf(line, "%d", Nchains);
  if(sscanf_status!=1){
    cout << "Error in parsing line 1 of the configuration file" << endl;
    return 0;
  }

  // line 2: Nreplica
  fgets_status=fgets(line, buffer_length, config);
  if(fgets_status==NULL){
    cout << "Error in loading line 2 of the configuration file" << endl;
    return 0;
  }
  sscanf_status=sscanf(line, "%d", Nreplica);
  if(sscanf_status!=1){
    cout << "Error in parsing line 2 of the configuration file" << endl;
    return 0;
  }

  // line 3-(Nreplica+2): inverse temperature beta
  int i;
  double* beta=new double[*Nreplica];
  beta_pointer[0]=beta;
  for(i=0;i<*Nreplica;i++){
    fgets_status=fgets(line, buffer_length, config);
    if(fgets_status==NULL){
      cout << "Error in loading line " << 3+i << " of the configuration file" << endl;
      return 0;
    }
    sscanf_status=sscanf(line, "%lf", &beta[i]);
    if(sscanf_status!=1){
      cout << "Error in parsing line " << 3+i << " of the configuration file" << endl;
      return 0;
    }
  }


  
  // line 3+Nreplica: Number of MC steps
  fgets_status=fgets(line, buffer_length, config);
  if(fgets_status==NULL){
    cout << "Error in loading line " << 3+*Nreplica << " of the configuration file" << endl;
    return 0;
  }
  sscanf_status=sscanf(line, "%d", Nstep);
  if(sscanf_status!=1){
    cout << "Error in parsing line " << 3+*Nreplica << " of the configuration file" << endl;
    return 0;
  }

  // line 4+Nreplica: input_file_format
  fgets_status=fgets(line, buffer_length, config);
  if(fgets_status==NULL){
    cout << "Error in loading line " << 4+*Nreplica << " of the configuration file" << endl;
    return 0;
  }
  sscanf_status=sscanf(line, "%s", input_file_format);
  if(sscanf_status!=1){
    cout << "Error in parsing line " << 4+*Nreplica << " of the configuration file" << endl;
    return 0;
  }
  
  // line 5+Nreplica: state_file_format
  fgets_status=fgets(line, buffer_length, config);
  if(fgets_status==NULL){
    cout << "Error in loading line " << 5+*Nreplica << " of the configuration file" << endl;
    return 0;
  }
  sscanf_status=sscanf(line, "%s", state_file_format);
  if(sscanf_status!=1){
    cout << "Error in parsing line " << 5+*Nreplica << " of the configuration file" << endl;
    return 0;
  }

  // line 6+Nreplica: accept_file_format
  fgets_status=fgets(line, buffer_length, config);
  if(fgets_status==NULL){
    cout << "Error in loading line " << 6+*Nreplica << " of the configuration file" << endl;
    return 0;
  }
  sscanf_status=sscanf(line, "%s", accept_file_format);
  if(sscanf_status!=1){
    cout << "Error in parsing line " << 6+*Nreplica << " of the configuration file" << endl;
    return 0;
  }

  // line 7+Nreplica: exchange_file_format
  fgets_status=fgets(line, buffer_length, config);
  if(fgets_status==NULL){
    cout << "Error in loading line " << 7+*Nreplica << " of the configuration file" << endl;
    return 0;
  }
  sscanf_status=sscanf(line, "%s", exchange_file_format);
  if(sscanf_status!=1){
    cout << "Error in parsing line " << 7+*Nreplica << " of the configuration file" << endl;
    return 0;
  }
  
  // line 8+Nreplica: log_file_format
  fgets_status=fgets(line, buffer_length, config);
  if(fgets_status==NULL){
    cout << "Error in loading line " << 8+*Nreplica << " of the configuration file" << endl;
    return 0;
  }
  sscanf_status=sscanf(line, "%s", log_file_format);
  if(sscanf_status!=1){
    cout << "Error in parsing line " << 8+*Nreplica << " of the configuration file" << endl;
    return 0;
  }

  
  // line 9+Nreplica: bond length J
  fgets_status=fgets(line, buffer_length, config);
  if(fgets_status==NULL){
    cout << "Error in loading line " << 9+*Nreplica << " of the configuration file" << endl;
    return 0;
  }
  sscanf_status=sscanf(line, "%lf", J);
  if(sscanf_status!=1){
    cout << "Error in parsing line " << 9+*Nreplica << " of the configuration file" << endl;
    return 0;
  }

  // line 10+Nreplica: shell_script
  fgets_status=fgets(line, buffer_length, config);
  if(fgets_status==NULL){
    cout << "Error in loading line " << 10+*Nreplica << " of the configuration file" << endl;
    return 0;
  }
  sscanf_status=sscanf(line, "%s", shell_script);
  if(sscanf_status!=1){
    cout << "Error in parsing line " << 10+*Nreplica << " of the configuration file" << endl;
    return 0;
  }

  // line 11+Nreplica: input2_file_format
  fgets_status=fgets(line, buffer_length, config);
  if(fgets_status==NULL){
    cout << "Error in loading line " << 11+*Nreplica << " of the configuration file" << endl;
    return 0;
  }
  sscanf_status=sscanf(line, "%s", input2_file_format);
  if(sscanf_status!=1){
    cout << "Error in parsing line " << 11+*Nreplica << " of the configuration file" << endl;
    return 0;
  }

  // line 12+Nreplica: Burnin
  fgets_status=fgets(line, buffer_length, config);
  if(fgets_status==NULL){
    cout << "Error in loading line " << 12+*Nreplica << " of the configuration file" << endl;
    return 0;
  }
  sscanf_status=sscanf(line, "%d", Burnin);
  if(sscanf_status!=1){
    cout << "Error in parsing line " << 12+*Nreplica << " of the configuration file" << endl;
    return 0;
  }

  // line 13+Nreplica: average file format
  fgets_status=fgets(line, buffer_length, config);
  if(fgets_status==NULL){
    cout << "Error in loading line " << 13+*Nreplica << " of the configuration file" << endl;
    return 0;
  }
  sscanf_status=sscanf(line, "%s", average_file_format);
  if(sscanf_status!=1){
    cout << "Error in parsing line " << 13+*Nreplica << " of the configuration file" << endl;
    return 0;
  }

  // line 14+Nreplica: variance file format
  fgets_status=fgets(line, buffer_length, config);
  if(fgets_status==NULL){
    cout << "Error in loading line " << 14+*Nreplica << " of the configuration file" << endl;
    return 0;
  }
  sscanf_status=sscanf(line, "%s", variance_file_format);
  if(sscanf_status!=1){
    cout << "Error in parsing line " << 14+*Nreplica << " of the configuration file" << endl;
    return 0;
  }

  // line 15+Nreplica: log file format
  fgets_status=fgets(line, buffer_length, config);
  if(fgets_status==NULL){
    cout << "Error in loading line " << 15+*Nreplica << " of the configuration file" << endl;
    return 0;
  }
  sscanf_status=sscanf(line, "%s", log2_file_format);
  if(sscanf_status!=1){
    cout << "Error in parsing line " << 15+*Nreplica << " of the configuration file" << endl;
    return 0;
  }

  // line 16+Nreplica: value format
  fgets_status=fgets(line, buffer_length, config);
  if(fgets_status==NULL){
    cout << "Error in loading line " << 16+*Nreplica << " of the configuration file" << endl;
    return 0;
  }
  sscanf_status=sscanf(line, "%s", value_format);
  if(sscanf_status!=1){
    cout << "Error in parsing line " << 16+*Nreplica << " of the configuration file" << endl;
    return 0;
  }

  
  // line 17+Nreplica: input3 file name
  fgets_status=fgets(line, buffer_length, config);
  if(fgets_status==NULL){
    cout << "Error in loading line " << 17+*Nreplica << " of the configuration file" << endl;
    return 0;
  }
  sscanf_status=sscanf(line, "%s", input3_file_format);
  if(sscanf_status!=1){
    cout << "Error in parsing line " << 17+*Nreplica << " of the configuration file" << endl;
    return 0;
  }

  // line 18+Nreplica: output file format
  fgets_status=fgets(line, buffer_length, config);
  if(fgets_status==NULL){
    cout << "Error in loading line " << 18+*Nreplica << " of the configuration file" << endl;
    return 0;
  }
  sscanf_status=sscanf(line, "%s", psrf_file_format);
  if(sscanf_status!=1){
    cout << "Error in parsing line " << 18+*Nreplica << " of the configuration file" << endl;
    return 0;
  }

  // line 19+Nreplica: log file
  fgets_status=fgets(line, buffer_length, config);
  if(fgets_status==NULL){
    cout << "Error in loading line " << 19+*Nreplica << " of the configuration file" << endl;
    return 0;
  }
  sscanf_status=sscanf(line, "%s", log3_file_format);
  if(sscanf_status!=1){
    cout << "Error in parsing line " << 19+*Nreplica << " of the configuration file" << endl;
    return 0;
  }
  
  return 1;
}
