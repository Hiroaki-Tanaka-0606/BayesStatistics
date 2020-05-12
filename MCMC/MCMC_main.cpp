#include "Model_selection.hpp"
#include "MC_Metropolis.hpp"
#include "../randomDistribution/mersenne_twister.hpp"
#include <iostream>

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

  double beta; // inverse temperature
  int Nstep; // Number of MC steps
  int buffer_length=256;
  char* state_file=new char[buffer_length]; // development of the state by MC steps
  char* accept_file=new char[buffer_length]; // log of acceptance
  double J; // bond strength
  int loadConfig_status=loadConfig(config, &beta, &Nstep, state_file, accept_file, &J);
  if(loadConfig_status!=1){
    cout << "Failed in loading configuration" << endl;
    return -1;
  }
  fclose(config);
  
  cout << "Inverse temperature: " << beta << endl;
  cout << "Number of MC steps: " << Nstep << endl;
  cout << "Output files: " << state_file << ", " << accept_file << endl;
  cout << "J: " << J << endl;

  Hamiltonian* h=new Hamiltonian(J);
  State* s=new State();
  State* s_next=new State();
  State* s_buff;
  int i,j;
  int num_params=s->num_params;
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
  
  for(i=0;i<Nstep;i++){
    for(j=0;j<num_params;j++){
      int accept=MC_Metropolis(s,s_next,h,j,sigma,beta,mt);
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
    
    char* state_string=s->print();
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
