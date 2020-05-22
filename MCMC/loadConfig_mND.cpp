#include "loadConfig_mND.hpp"
#include <iostream>
#include <string.h>

using namespace std;
int loadConfig(FILE* config, double* beta, int* Nstep, char* state_file, char* accept_file, char* sample_file, int* Nbin, double* sigma, int num_params, char* initial_state){
  int buffer_length=256;
  char line[buffer_length];
  char* fgets_status;
  int sscanf_status;
  
  // line 1: inverse temperature beta
  fgets_status=fgets(line, buffer_length, config);
  if(fgets_status==NULL){
    cout << "Error in loading line 1 of the configuration file" << endl;
    return 0;
  }
  sscanf_status=sscanf(line, "%lf", beta);
  if(sscanf_status!=1){
    cout << "Error in parsing line 1 of the configuration file" << endl;
    return 0;
  }

  // line 2: number of MC steps
  fgets_status=fgets(line, buffer_length, config);
  if(fgets_status==NULL){
    cout << "Error in loading line 2 of the configuration file" << endl;
    return 0;
  }
  sscanf_status=sscanf(line, "%d", Nstep);
  if(sscanf_status!=1){
    cout << "Error in parsing line 2 of the configuration file" << endl;
    return 0;
  }

  // line 3: state_file
  fgets_status=fgets(line, buffer_length, config);
  if(fgets_status==NULL){
    cout << "Error in loading line 3 of the configuration file" << endl;
    return 0;
  }
  sscanf_status=sscanf(line, "%s", &state_file[0]);
  if(sscanf_status!=1){
    cout << "Error in parsing line 3 of the configuration file" << endl;
    return 0;
  }

  // line 4: accept_file
  fgets_status=fgets(line, buffer_length, config);
  if(fgets_status==NULL){
    cout << "Error in loading line 4 of the configuration file" << endl;
    return 0;
  }
  sscanf_status=sscanf(line, "%s", &accept_file[0]);
  if(sscanf_status!=1){
    cout << "Error in parsing line 4 of the configuration file" << endl;
    return 0;
  }
  // line 5: sample_file
  fgets_status=fgets(line, buffer_length, config);
  if(fgets_status==NULL){
    cout << "Error in loading line 5 of the configuration file" << endl;
    return 0;
  }
  sscanf_status=sscanf(line, "%s", &sample_file[0]);
  if(sscanf_status!=1){
    cout << "Error in parsing line 5 of the configuration file" << endl;
    return 0;
  }

  // line 6: Number of bins
  fgets_status=fgets(line, buffer_length, config);
  if(fgets_status==NULL){
    cout << "Error in loading line 6 of the configuration file" << endl;
    return 0;
  }
  sscanf_status=sscanf(line, "%d", Nbin);
  if(sscanf_status!=1){
    cout << "Error in parsing line 6 of the configuration file" << endl;
    return 0;
  }  

  // line 7: development parameter
  int i;
  for(i=0;i<num_params;i++){
    sscanf_status=fscanf(config, "%lf", &sigma[i]);
    if(sscanf_status!=1){
      cout << "Error in parsing line 7 of the configuration file" << endl;
      return 0;
    }
  }
  // pass to the end of line 7
  fgets_status=fgets(line, buffer_length, config);
  if(fgets_status==NULL){
    cout << "Error in loading line 7 of the configuration file" << endl;
    return 0;
  }
  
  // line 8: (optional) initial state
  int state_length=strlen(initial_state);
  fgets_status=fgets(initial_state, state_length, config);
  if(fgets_status!=NULL){
    cout << "Read initial state" << endl;
  }else{
    cout << "Not read initial state" << endl;
  }
  
  return 1;
}
