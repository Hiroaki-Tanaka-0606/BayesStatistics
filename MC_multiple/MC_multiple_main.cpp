#include <iostream>
#include "../MCMC/Model_selection.hpp"

using namespace std;

int loadConfig(FILE* config, int* Nchains, double* beta, int* Nstep, char* input_file_format, char* state_file_format, char* accept_file_format, char* log_file_format, double* J, char* shell_script);


int main(int argc, char** argv){
  cout << "Generate multiple MCMC configuration files" << endl;
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

  int Nchains; // number of Marcov chains
  double beta; // inverse temperature
  int Nstep; // Number of MC steps
  int buffer_length=256;
  char* input_file_format=new char[buffer_length]; // name of the input file, including %d (chain index)
  char* state_file_format=new char[buffer_length]; // name of the output file for development of the state, including %d (chain index)
  char* accept_file_format=new char[buffer_length]; // name of the output file for log of acceptance, including %d (chain index)
  char* log_file_format=new char[buffer_length]; // name of the log (STDOUT), including %d (chain index)
  double J; //bond strength
  char* shell_script=new char[buffer_length]; // name of the shell script file
  int loadConfig_status=loadConfig(config, &Nchains, &beta, &Nstep, input_file_format, state_file_format, accept_file_format, log_file_format, &J, shell_script);
  if(loadConfig_status!=1){
    cout << "Failed in loading configuration" << endl;
    return -1;
  }

  cout << "Number of Marcov chains: " << Nchains << endl;
  cout << "Inverse temperature: " << beta << endl;
  cout << "Number of MC steps: " << Nstep << endl;
  cout << "J: " << J << endl;
  cout << "Shell script file: " << shell_script << endl;

  // generate shell script file
  FILE* shell=fopen(shell_script, "w");
  if(shell==NULL){
    cout << "Error in opening the shell script file" << endl;
    return -1;
  }
  
  FILE* input;
  int i;
  char* input_file=new char[buffer_length];
  State* s=new State();
  char* initial_state;
  for(i=0;i<Nchains;i++){
    // open input file
    sprintf(input_file, input_file_format, i+1);
    input=fopen(input_file, "w");
    if(input==NULL){
      cout << "Error in opening the input file" << endl;
      return -1;
    }
    // line 1: beta
    fprintf(input, "%.4e # beta\n", beta);
    // line 2: Nstep
    fprintf(input, "%d # Nstep\n", Nstep);
    // line 3: state file
    fprintf(input, state_file_format, i+1);
    fprintf(input, " # state file\n");
    // line 4: accept file
    fprintf(input, accept_file_format, i+1);
    fprintf(input, " # accept file\n");
    // line 5: J
    fprintf(input, "%.4e # J\n", J);
    // line 6: initial state
    s->init_rand();
    initial_state=s->print();
    fprintf(input, "%s # initial state\n", initial_state);
    
    fclose(input);

    // add MCMC command to shell file
    fprintf(shell, "./MCMC.o %s", input_file);
    fprintf(shell, " > ");
    fprintf(shell, log_file_format, i+1);
    fprintf(shell, "\n");
    fprintf(shell, "echo 'Done %d / %d'\n", i+1, Nchains);
  }

  
  fclose(shell);
  return 1;
}

int loadConfig(FILE* config, int* Nchains, double* beta, int* Nstep, char* input_file_format, char* state_file_format, char* accept_file_format, char* log_file_format, double* J, char* shell_script){
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

  // line 2: beta
  fgets_status=fgets(line, buffer_length, config);
  if(fgets_status==NULL){
    cout << "Error in loading line 2 of the configuration file" << endl;
    return 0;
  }
  sscanf_status=sscanf(line, "%lf", beta);
  if(sscanf_status!=1){
    cout << "Error in parsing line 2 of the configuration file" << endl;
    return 0;
  }

  
  // line 3: Number of MC steps
  fgets_status=fgets(line, buffer_length, config);
  if(fgets_status==NULL){
    cout << "Error in loading line 3 of the configuration file" << endl;
    return 0;
  }
  sscanf_status=sscanf(line, "%d", Nstep);
  if(sscanf_status!=1){
    cout << "Error in parsing line 3 of the configuration file" << endl;
    return 0;
  }

  // line 4: input_file_format
  fgets_status=fgets(line, buffer_length, config);
  if(fgets_status==NULL){
    cout << "Error in loading line 4 of the configuration file" << endl;
    return 0;
  }
  sscanf_status=sscanf(line, "%s", input_file_format);
  if(sscanf_status!=1){
    cout << "Error in parsing line 4 of the configuration file" << endl;
    return 0;
  }
  
  // line 5: state_file_format
  fgets_status=fgets(line, buffer_length, config);
  if(fgets_status==NULL){
    cout << "Error in loading line 5 of the configuration file" << endl;
    return 0;
  }
  sscanf_status=sscanf(line, "%s", state_file_format);
  if(sscanf_status!=1){
    cout << "Error in parsing line 5 of the configuration file" << endl;
    return 0;
  }

  // line 6: accept_file_format
  fgets_status=fgets(line, buffer_length, config);
  if(fgets_status==NULL){
    cout << "Error in loading line 6 of the configuration file" << endl;
    return 0;
  }
  sscanf_status=sscanf(line, "%s", accept_file_format);
  if(sscanf_status!=1){
    cout << "Error in parsing line 6 of the configuration file" << endl;
    return 0;
  }

  // line 7: log_file_format
  fgets_status=fgets(line, buffer_length, config);
  if(fgets_status==NULL){
    cout << "Error in loading line 7 of the configuration file" << endl;
    return 0;
  }
  sscanf_status=sscanf(line, "%s", log_file_format);
  if(sscanf_status!=1){
    cout << "Error in parsing line 7 of the configuration file" << endl;
    return 0;
  }

  
  // line 8: bond length J
  fgets_status=fgets(line, buffer_length, config);
  if(fgets_status==NULL){
    cout << "Error in loading line 8 of the configuration file" << endl;
    return 0;
  }
  sscanf_status=sscanf(line, "%lf", J);
  if(sscanf_status!=1){
    cout << "Error in parsing line 8 of the configuration file" << endl;
    return 0;
  }

  // line 9: shell_script
  fgets_status=fgets(line, buffer_length, config);
  if(fgets_status==NULL){
    cout << "Error in loading line 9 of the configuration file" << endl;
    return 0;
  }
  sscanf_status=sscanf(line, "%s", shell_script);
  if(sscanf_status!=1){
    cout << "Error in parsing line 9 of the configuration file" << endl;
    return 0;
  }

  return 1;
}
