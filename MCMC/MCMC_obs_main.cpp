#include "Model_selection.hpp"
#include <iostream>
#include <string.h>

using namespace std;

int main(int argc, char** argv){
  // initial message
  cout << "MCMC observable calculation" << endl;
  
  if(argc<3){
    printf("Usage: %s config_file\n", argv[0]);
    return 0;
  }

  // load configuration
  char* config_file=argv[1];
  char* output_file=argv[2];
  cout << "Configuration file: " << config_file << endl;
  cout << "Output file: " << output_file << endl;
  FILE* config=fopen(config_file, "r");
  if(config==NULL){
    cout << "Error in opening the config file" << endl;
    return -1;
  }

  double beta; // inverse temperature
  int Nstep; // Number of MC steps
  int buffer_length=256;
  char* state_file=new char[buffer_length]; // development of the state by MC steps
  char* accept_file=new char[buffer_length]; // log of acceptance
  double J; // bond strength
  char* initial_state=new char[buffer_length]; // not used
  int loadConfig_status=loadConfig(config, &beta, &Nstep, state_file, accept_file, &J, initial_state);
  if(loadConfig_status!=1){
    cout << "Failed in loading configuration" << endl;
    return -1;
  }
  fclose(config);
  
  cout << "Inverse temperature: " << beta << endl;
  cout << "Number of MC steps: " << Nstep << endl;
  cout << "Input file: " << state_file << endl;
  cout << "J: " << J << endl;

  Hamiltonian* h=new Hamiltonian(J);
  State* s=new State();
  Observables* obs=new Observables(s);
  char* header=StateHeader();
  char* s_print=s->print();
  int data_length=max(strlen(header), strlen(s_print))+10; // +10 for \r, \n, \0
  char* buffer=new char[data_length];
  FILE* state_input=fopen(state_file, "r");
  int state_count=0;
  int line_count=0;

  FILE* obs_output=fopen(output_file, "w");
  fprintf(obs_output, "# Energy Mag Mag^2\n");
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
      double e=h->energy(s);
      int mag=obs->Mag();
      int mag2=obs->Mag2();
      fprintf(obs_output, "%12.4e %4d %4d\n", e, mag, mag2);
    }else{
      // load fail
      cout << "Failed in loading line " << line_count << ", " << load_status << endl;
    }
  } 
  fclose(state_input);
  fclose(obs_output);

  if(state_count==Nstep){
    cout << "Calculation succeeded" << endl;
  }else{
    cout << "Calculation failed" <<endl;
  }
}
