
#include "loadConfig_2D_Ising.hpp"
#include "../MCMC/Model_selection.hpp"
#include "../MCMC/MC_Metropolis.hpp"
#include "../randomDistribution/mersenne_twister.hpp"
#include <iostream>
#include <string.h>

using namespace std;

int main(int argc, char** argv){
  // initial message
  cout << "RXMC observable calculation" << endl;
  
  if(argc<3){
    printf("Usage: %s config_file output_file_format\n", argv[0]);
    return 0;
  }

  // load configuration
  char* config_file=argv[1];
  char* output_file_format=argv[2];
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

  int i,j,k;

  char* s_print=s->print();
  int state_length=strlen(s_print);
  
  int Nreplica; // number of replicas
  double* beta; // inverse temperatures of replicas
  int Nstep; // Number of MC steps
  int buffer_length=256;
  char* state_file_format=new char[buffer_length]; // development of the state by MC steps
  char* accept_file_format=new char[buffer_length]; // log of acceptance
  char* exchange_file_name=new char[buffer_length]; // log of exchange
  double J; // bond strength
  int state_buffer_length=state_length+10; // +10 for \r, \n, \0
  char** initial_state; // initial states for each replica

  // pointer for loadConfig arguments
  double** beta_pointer=new double*[1];
  char*** initial_state_pointer=new char**[1];

  int loadConfig_status=loadConfig_RXMC(config, &Nreplica, beta_pointer, &Nstep, state_file_format, accept_file_format, exchange_file_name, &J, state_buffer_length, initial_state_pointer);
  if(loadConfig_status!=1){
    cout << "Failed in loading configuration" << endl;
    return -1;
  }
  fclose(config);

  beta=beta_pointer[0];
  initial_state=initial_state_pointer[0];
  
  cout << "Number of replicas: " << Nreplica << endl;
  cout << "Inverse temperatures: ";
  for(i=0;i<Nreplica;i++){
    cout << beta[i] << ", ";
  }
  cout << endl;
  cout << "Number of MC steps: " << Nstep << endl;
  cout << "J: " << J << endl;


  Hamiltonian* h=new Hamiltonian(J);
  Observables* obs=new Observables(s);
  char* header=StateHeader();
  int data_length=max(strlen(header), strlen(s_print))+10; // +10 for \r, \n, \0
  char* buffer=new char[data_length];

  char* state_file=new char[buffer_length];
  char* output_file=new char[buffer_length];

  for(i=0;i<Nreplica;i++){
    int state_count=0;
    int line_count=0;
    sprintf(state_file, state_file_format, i+1);
    sprintf(output_file, output_file_format, i+1);
    FILE* state_input=fopen(state_file, "r");
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

}
