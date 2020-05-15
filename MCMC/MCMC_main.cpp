#include "Model_selection.hpp"
#include "MC_Metropolis.hpp"
#include "../randomDistribution/mersenne_twister.hpp"
#include <iostream>
#include <string.h>

using namespace std;

int main(int argc, char** argv){
  // initial message
  cout << "MCMC calculation" << endl;
  
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

  // State used in the simulation
  State* s=new State();
  State* s_next=new State();
  State* s_buff;
  int num_params=s->num_params;

  int i,j;

  char* s_print=s->print();
  int state_length=strlen(s_print);
  
  double beta; // inverse temperature
  int Nstep; // Number of MC steps
  int buffer_length=256;
  char* state_file=new char[buffer_length]; // development of the state by MC steps
  char* accept_file=new char[buffer_length]; // log of acceptance
  double J; // bond strength
  char* initial_state=new char[state_length+10]; // +10 for \r, \n, \0
  for(i=0;i<state_length+9;i++){
    initial_state[i]='a';
  }
  initial_state[state_length+9]='\0';

  int loadConfig_status=loadConfig(config, &beta, &Nstep, state_file, accept_file, &J, initial_state);
  if(loadConfig_status!=1){
    cout << "Failed in loading configuration" << endl;
    return -1;
  }
  fclose(config);
  
  cout << "Inverse temperature: " << beta << endl;
  cout << "Number of MC steps: " << Nstep << endl;
  cout << "Output files: " << state_file << ", " << accept_file << endl;
  cout << "J: " << J << endl;

  // load initial state
  int loadInit_status=s->load(initial_state);
  if(loadInit_status==1){
    cout << "Initial state is loaded" << endl;
  }else{
    cout << "Initial state is not loaded: error " << loadInit_status << endl;
    s=new State();
  }
  cout << "Initial state: " << s->print() << endl;


  Hamiltonian* h=new Hamiltonian(J);
  double sigma=0;
  MT* mt=new MT();

  FILE* state_output=fopen(state_file, "w");
  FILE* accept_output=fopen(accept_file, "w");

  char* header=StateHeader();
  fprintf(state_output, "# MCMC calculation\n");
  fprintf(state_output, "# Development of states\n");
  fprintf(state_output, header);
  fprintf(state_output, "\n");

  
  fprintf(accept_output, "# MCMC calculation\n");
  fprintf(accept_output, "# Acceptance\n");
  fprintf(accept_output, header);
  fprintf(accept_output, "\n");
  
  int* accept_count=new int[num_params];
  for(i=0;i<num_params;i++){
    accept_count[i]=0;
  }
  
  char* state_string=s->print();
  int accept;
  for(i=0;i<Nstep;i++){
    for(j=0;j<num_params;j++){
      accept=MC_Metropolis(s,s_next,h,j,sigma,beta,mt);
      if(accept==1){
	//replace s and s_next
	s_buff=s_next;
	s_next=s;
	s=s_buff;
	accept_count[j]++;
      }
      fprintf(accept_output, "%1d ", accept);
    }
    fprintf(accept_output, "\n");
    delete state_string;
    state_string=s->print();
    fprintf(state_output, state_string);
    fprintf(state_output, "\n");
  }

  // output accept count
  fprintf(accept_output, "# Number of acceptance at each lattice\n");
  fprintf(accept_output, "# ");
  for(i=0;i<num_params;i++){
    fprintf(accept_output, "%d ", accept_count[i]);
  }
  fprintf(accept_output, "\n");
  
  fprintf(accept_output, "# Acceptance rate at each lattice\n");
  fprintf(accept_output, "# ");
  for(i=0;i<num_params;i++){
    fprintf(accept_output, "%f ", accept_count[i]*1.0/Nstep);
  }
  fprintf(accept_output, "\n");
  

  fclose(state_output);
  fclose(accept_output);
}
