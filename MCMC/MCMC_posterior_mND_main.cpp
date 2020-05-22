#include "Model_selection_mND.hpp"
#include <iostream>
#include <string.h>

using namespace std;

int loadConfig2(FILE* config2, double* min_x, double* max_x, int* num_points, int* Burnin, char* output_file);

int main(int argc, char** argv){
  // initial message
  cout << "# MCMC Posterior distribution calculation" << endl;
  
  if(argc<3){
    printf("Usage: %s config_file config_file2\n", argv[0]);
    return 0;
  }

  // load configuration
  char* config_file=argv[1];
  char* config_file2=argv[2];
  
  cout << "Configuration file: " << config_file << ", " << config_file2 << endl;
  FILE* config=fopen(config_file, "r");
  if(config==NULL){
    cout << "Error in opening the config file" << endl;
    return -1;
  }

  int i,j;
  
  State* s=new State();
  char* s_print=s->print();
  int state_length=strlen(s_print);
  int num_params=s->num_params;
  
  // <---- Parameters loaded from config_file
  
  double beta; // inverse temperature
  int Nstep; // Number of MC steps
  int buffer_length=256;
  char* state_file=new char[buffer_length]; // development of the state by MC steps
  char* accept_file=new char[buffer_length]; // log of acceptance
  char* sample_file=new char[buffer_length]; // sample data
  int Nbin; // number of bins
  double* sigma=new double[num_params]; // development parameter
  char* initial_state=new char[state_length+10]; // +10 for \r, \n, \0

  // ---->
  
  for(i=0;i<state_length+9;i++){
    initial_state[i]='a';
  }
  initial_state[state_length+9]='\0';

  // <---- Load configuration from config_file
  
  int loadConfig_status=loadConfig(config, &beta, &Nstep, state_file, accept_file, sample_file, &Nbin, sigma, num_params, initial_state);
  
  // ---->
  
  if(loadConfig_status!=1){
    cout << "Failed in loading configuration" << endl;
    return -1;
  }
  fclose(config);

  FILE* config2=fopen(config_file2, "r");
  if(config2==NULL){
    cout << "Error in opening the config2 file" << endl;
  }

  double min_x, max_x;
  int num_points;
  int Burnin;
  char* output_file=new char[buffer_length];
  int loadConfig2_status=loadConfig2(config2, &min_x, &max_x, &num_points, &Burnin, output_file);

  if(loadConfig2_status!=1){
    cout << "Failed in loading configuration2" << endl;
    return -1;
  }
  
  
  cout << "Inverse temperature: " << beta << endl;
  cout << "Number of MC steps: " << Nstep << endl;
  cout << "Posterior distribution range: " << min_x << " to " << max_x << endl;
  cout << "Output file: " << output_file << endl;

  char* header=StateHeader();
  int data_length=max(strlen(header), strlen(s_print))+10; // +10 for \r, \n, \0
  char* buffer=new char[data_length];
  FILE* state_input=fopen(state_file, "r");
  int state_count=0;
  int line_count=0;

  FILE* output=fopen(output_file, "w");
  fprintf(output, "# Posterior distribution\n");
  fprintf(output, "# x p(x)\n");

  double* x=new double[num_points];
  double* posterior=new double[num_points];
  for(i=0;i<num_points;i++){
    posterior[i]=0;
    x[i]=min_x+(max_x-min_x)*i/(num_points-1);
  }
  while(fgets(buffer, data_length, state_input)!=NULL){
    line_count++;
    if(buffer[0]=='#'){
      continue;
    }
    
    int load_status=s->load(buffer);
    if(load_status==1){
      // load succeeded
      state_count++;
      // we do not need to update obs, because the State is passed to obs by the pointer
      if(state_count>Burnin){
	for(j=0;j<num_points;j++){
	  posterior[j]+=s->probability(x[j]);
	}
      }
    }else{
      // load fail
      cout << "Failed in loading line " << line_count << ", " << load_status << endl;
    }
  } 
  fclose(state_input);

  if(state_count==Nstep){
    cout << "Calculation succeeded" << endl;
    for(i=0;i<num_points;i++){
      fprintf(output, "%12.4e %12.4e\n", x[i], posterior[i]/(Nstep-Burnin));
    }
  }else{
    cout << "Calculation failed" <<endl;
  }

  fclose(output);
}


int loadConfig2(FILE* config2, double* min_x, double* max_x, int* num_points, int* Burnin, char* output_file){
  int buffer_length=256;
  char line[buffer_length];
  char* fgets_status;
  int sscanf_status;

  // line 1: min_x, max_x, num_points
  fgets_status=fgets(line, buffer_length, config2);
  if(fgets_status==NULL){
    cout << "Error in loading line 1 of the configuration file" << endl;
    return 0;
  }
  sscanf_status=sscanf(line, "%lf%lf%d", min_x, max_x, num_points);
  if(sscanf_status!=3){
    cout << "Error in parsing line 1 of the configuration file" << endl;
    return 0;
  }

  // line 2: Burnin
  fgets_status=fgets(line, buffer_length, config2);
  if(fgets_status==NULL){
    cout << "Error in loading line 2 of the configuration file" << endl;
    return 0;
  }
  sscanf_status=sscanf(line, "%d", Burnin);
  if(sscanf_status!=1){
    cout << "Error in parsing line 2 of the configuration file" << endl;
    return 0;
  }

  // line 3: output_file
  fgets_status=fgets(line, buffer_length, config2);
  if(fgets_status==NULL){
    cout << "Error in loading line 3 of the configuration file" << endl;
    return 0;
  }
  sscanf_status=sscanf(line, "%s", output_file);
  if(sscanf_status!=1){
    cout << "Error in parsing line 3 of the configuration file" << endl;
    return 0;
  }
  return 1;
}
