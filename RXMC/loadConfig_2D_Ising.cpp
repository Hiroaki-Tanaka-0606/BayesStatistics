#include "loadConfig_2D_ising.hpp"
#include <iostream>
#include <string.h>

using namespace std;
int loadConfig_RXMC(FILE* config, int* Nreplica, double** beta_pointer, int* Nstep, char* state_file_format, char* accept_file_format, char* exchange_file_name, double* J, int state_buffer_length, char*** initial_state_pointer){
  int buffer_length=256;
  char line[buffer_length];
  char* fgets_status;
  int sscanf_status;


  // line 1: Number of replicas
  fgets_status=fgets(line, buffer_length, config);
  if(fgets_status==NULL){
    cout << "Error in loading line 1 of the configuration file" << endl;
    return 0;
  }
  sscanf_status=sscanf(line, "%d", Nreplica);
  if(sscanf_status!=1){
    cout << "Error in parsing line 1 of the configuration file" << endl;
    return 0;
  }

  // line 2-(Nreplica+1): inverse temperature beta
  int i;
  double* beta=new double[*Nreplica];
  beta_pointer[0]=beta;
  for(i=0;i<*Nreplica;i++){
    fgets_status=fgets(line, buffer_length, config);
    if(fgets_status==NULL){
      cout << "Error in loading line " << 2+i << " of the configuration file" << endl;
      return 0;
    }
    sscanf_status=sscanf(line, "%lf", &beta[i]);
    if(sscanf_status!=1){
      cout << "Error in parsing line " << 2+i << " of the configuration file" << endl;
      return 0;
    }
  }

  // line 2+Nreplica: number of MC steps
  fgets_status=fgets(line, buffer_length, config);
  if(fgets_status==NULL){
    cout << "Error in loading line " << 2+*Nreplica << " of the configuration file" << endl;
    return 0;
  }
  sscanf_status=sscanf(line, "%d", Nstep);
  if(sscanf_status!=1){
    cout << "Error in parsing line " << 2+*Nreplica << " of the configuration file" << endl;
    return 0;
  }

  // line 3+Nreplica: state_file_format
  fgets_status=fgets(line, buffer_length, config);
  if(fgets_status==NULL){
    cout << "Error in loading line " << 3+*Nreplica << " of the configuration file" << endl;
    return 0;
  }
  sscanf_status=sscanf(line, "%s", &state_file_format[0]);
  if(sscanf_status!=1){
    cout << "Error in parsing line " << 3+*Nreplica << " of the configuration file" << endl;
    return 0;
  }

  // line 4+Nreplica: accept_file_format
  fgets_status=fgets(line, buffer_length, config);
  if(fgets_status==NULL){
    cout << "Error in loading line " << 4+*Nreplica << " of the configuration file" << endl;
    return 0;
  }
  sscanf_status=sscanf(line, "%s", &accept_file_format[0]);
  if(sscanf_status!=1){
    cout << "Error in parsing line " << 4+*Nreplica << " of the configuration file" << endl;
    return 0;
  }

  // line 5+Nreplica: exchange_file_name
  fgets_status=fgets(line, buffer_length, config);
  if(fgets_status==NULL){
    cout << "Error in loading line " << 5+*Nreplica << " of the configuration file" << endl;
    return 0;
  }
  sscanf_status=sscanf(line, "%s", &exchange_file_name[0]);
  if(sscanf_status!=1){
    cout << "Error in parsing line " << 5+*Nreplica << " of the configuration file" << endl;
    return 0;
  }
  
  // line 6+Nreplica: bond strength
  fgets_status=fgets(line, buffer_length, config);
  if(fgets_status==NULL){
    cout << "Error in loading line " << 6+*Nreplica << " of the configuration file" << endl;
    return 0;
  }
  sscanf_status=sscanf(line, "%lf", J);
  if(sscanf_status!=1){
    cout << "Error in parsing line " << 6+*Nreplica << " of the configuration file" << endl;
    return 0;
  }

  // line (7+Nreplica)-:(6+2*Nreplica) (optional) initial state
  char** initial_state=new char*[*Nreplica];
  initial_state_pointer[0]=initial_state;
  for(i=0;i<*Nreplica;i++){
    initial_state[i]=new char[state_buffer_length+1];
    fgets_status=fgets(&initial_state[i][0], state_buffer_length, config);
    if(fgets_status!=NULL){
      cout << "Read initial state for replica " << i+1 << endl;
    }else{
      cout << "Not read initial state for replica " << i+1 << endl;
    }
  }
  
  return 1;
}
