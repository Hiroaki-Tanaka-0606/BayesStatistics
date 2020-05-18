
#include "loadConfig_2D_Ising.hpp"
#include "../MCMC/Model_selection.hpp"
#include "../MCMC/MC_Metropolis.hpp"
#include "../randomDistribution/mersenne_twister.hpp"
#include <iostream>
#include <string.h>

using namespace std;

int main(int argc, char** argv){
  // initial message
  cout << "RXMC calculation" << endl;
  
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
  
  // load initial state
  State** ss=new State*[Nreplica];
  for(i=0;i<Nreplica;i++){
    ss[i]=new State();
    int loadInit_status=ss[i]->load(initial_state[i]);
    if(loadInit_status==1){
      cout << "Initial state for replica " << i+1 << " is loaded" << endl;
    }else{
      cout << "Initial state for replica " << i+1 << " is not loaded: error " << loadInit_status << endl;
      ss[i]=new State();
    }
    cout << "Initial state for replica " << i+1 << ": " << ss[i]->print() << endl;
  }

  
  Hamiltonian* h=new Hamiltonian(J);
  double sigma=0;
  MT* mt=new MT();

  FILE** state_output=new FILE*[Nreplica];
  FILE** accept_output=new FILE*[Nreplica];
  char* file_name=new char[buffer_length];
  char* header=StateHeader();
  for(i=0;i<Nreplica;i++){
    // state_file
    sprintf(file_name, state_file_format, i+1);
    state_output[i]=fopen(file_name, "w");
    fprintf(state_output[i], "# MCMC calculation for replica %d (beta = %f)\n", i+1, beta[i]);
    fprintf(state_output[i], "# Development of states\n");
    fprintf(state_output[i], header);
    fprintf(state_output[i], "\n");

    // accept_file
    sprintf(file_name, accept_file_format, i+1);
    accept_output[i]=fopen(file_name, "w");
    fprintf(accept_output[i], "# MCMC calculation for replica %d (beta = %f)\n", i+1, beta[i]);
    fprintf(accept_output[i], "# Acceptance\n");
    fprintf(accept_output[i], header);
    fprintf(accept_output[i], "\n");
  }
  // exchange file
  FILE* exchange_output=fopen(exchange_file_name, "w");
  fprintf(exchange_output, "# RXMC calculation\n");
  fprintf(exchange_output, "# Exchange\n");
  fprintf(exchange_output, "# ");
  for(i=0;i<Nreplica-1;i++){
    fprintf(exchange_output, "%d<->%d ", i+1, i+2);
  }
  fprintf(exchange_output, "\n");

  int** accept_count=new int*[Nreplica];
  for(i=0;i<Nreplica;i++){
    accept_count[i]=new int[num_params];
    for(j=0;j<num_params;j++){
      accept_count[i][j]=0;
    }
  }
  int* exchange_count=new int[Nreplica-1];
  for(i=0;i<Nreplica-1;i++){
    exchange_count[i]=0;
  }

  // MC step
  char* state_string=s->print();
  int accept;
  for(i=0;i<Nstep;i++){
    // MC step
    for(j=0;j<Nreplica;j++){
      for(k=0;k<num_params;k++){
	accept=MC_Metropolis(ss[j],s_next,h,j,sigma,beta[j],mt);
	if(accept==1){
	  //replace s and s_next
	  s_buff=s_next;
	  s_next=ss[j];
	  ss[j]=s_buff;
	  accept_count[j][k]++;
	}
	fprintf(accept_output[j], "%1d ", accept);
      }
      fprintf(accept_output[j], "\n");
    }

    // RX step
    int index_low;
    int index_high;
    if(i%2==0){
      index_low=0;
      index_high=1;
    }else{
      index_low=1;
      index_high=2;
    }

    while(index_high<Nreplica){
      double e1=h->energy(ss[index_low]);
      double e2=h->energy(ss[index_high]);
      double beta1=beta[index_low];
      double beta2=beta[index_high];

      double threshold=exp(-(beta1-beta2)*(e2-e1));
      double rand=mt->rand();

      // space-filling
      if(i%2==1){
	fprintf(exchange_output, " 0 ");
      }
      int accept;
      if(rand<threshold){
	// accept -> exchange
	s_buff=ss[index_low];
	ss[index_low]=ss[index_high];
	ss[index_high]=s_buff;
	accept=1;
	exchange_count[index_low]++;
      }else{
	// reject
	accept=-1;
      }
      fprintf(exchange_output, "%+2d ", accept);
      // space-filling
      if(i%2==0 && !((Nreplica-1)%2==1 && index_high==Nreplica-1)){
	fprintf(exchange_output, " 0 ");
      }
      index_low+=2;
      index_high+=2;
    }
    if((Nreplica-1)%2==1 && index_low == Nreplica-1){
      fprintf(exchange_output, " 0 ");
    }
    fprintf(exchange_output, "\n");

    // output state
    for(j=0;j<Nreplica;j++){
      delete state_string;
      state_string=ss[j]->print();
      fprintf(state_output[j], state_string);
      fprintf(state_output[j], "\n");
    }
  }

  // output accept count
  for(j=0;j<Nreplica;j++){
    fprintf(accept_output[j], "# Number of acceptance at each lattice\n");
    fprintf(accept_output[j], "# ");
    for(k=0;k<num_params;k++){
      fprintf(accept_output[j], "%d ", accept_count[j][k]);
    }
    fprintf(accept_output[j], "\n");
  
    fprintf(accept_output[j], "# Acceptance rate at each lattice\n");
    fprintf(accept_output[j], "# ");
    for(k=0;k<num_params;k++){
      fprintf(accept_output[j], "%f ", accept_count[j][k]*1.0/Nstep);
    }
    fprintf(accept_output[j], "\n");
  }

  // output exchange count
  int* exchange_trial=new int[Nreplica-1];
  fprintf(exchange_output, "# Number of exchange trials\n# ");
  for(j=0;j<Nreplica-1;j++){
    if(Nstep%2==0){
      // Nstep/2
      fprintf(exchange_output, "%d ", Nstep/2);
      exchange_trial[j]=Nstep/2;
    }else{
      if(j%2==0){
	// ceil(Nstep/2)
	fprintf(exchange_output, "%d ", Nstep/2+1);
	exchange_trial[j]=Nstep/2+1;
      }else{
	// floor(Nstep/2)
	fprintf(exchange_output, "%d ", Nstep/2);
	exchange_trial[j]=Nstep/2;
      }
    }
  }
  fprintf(exchange_output, "\n");
  fprintf(exchange_output, "# Number of exchange\n# ");
  for(j=0;j<Nreplica-1;j++){
    fprintf(exchange_output, "%d ", exchange_count[j]);
  }
  fprintf(exchange_output, "\n");
  fprintf(exchange_output, "# exchange rate\n# ");
  for(j=0;j<Nreplica-1;j++){
    fprintf(exchange_output, "%lf ", exchange_count[j]*1.0/exchange_trial[j]);
  }
  fprintf(exchange_output, "\n");
  
    
  // close files
  for(j=0;j<Nreplica;j++){
    fclose(state_output[j]);
    fclose(accept_output[j]);
  }
  fclose(exchange_output);  
}
