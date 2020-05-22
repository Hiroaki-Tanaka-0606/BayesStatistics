#include "loadConfig_mND.hpp"
#include "../MCMC/Model_selection_mND.hpp"
#include "../MCMC/MC_Metropolis.hpp"
#include "../randomDistribution/mersenne_twister.hpp"
#include <iostream>
#include <string.h>
#include "../MCMC/Acceptance_optimization.hpP"

using namespace std;

int loadConfig2(FILE* config2, double** acceptance_target, int Nreplica, int num_params, double* error_coeff, double* error_stop);

int main(int argc, char** argv){
  // initial message
  cout << "RXMC calculation" << endl;
  
  if(argc<3){
    printf("Usage: %s config_file\n", argv[0]);
    return 0;
  }

  // load configuration
  char* config_file=argv[1];
  char* config2_file=argv[2];
  cout << "Configuration files: " << config_file << ", " << config2_file << endl;
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
  char* sample_file=new char[buffer_length]; //sample data
  int Nbin; // number of bins
  double** sigma; // development parameter
  int state_buffer_length=state_length+10; // +10 for \r, \n, \0
  char** initial_state; // initial states for each replica

  // pointer for loadConfig arguments
  double** beta_pointer=new double*[1];
  char*** initial_state_pointer=new char**[1];
  double*** sigma_pointer=new double**[1];

  int loadConfig_status=loadConfig_RXMC(config, &Nreplica, beta_pointer, &Nstep, state_file_format, accept_file_format, exchange_file_name, sample_file, &Nbin, sigma_pointer, num_params,  state_buffer_length, initial_state_pointer);
  if(loadConfig_status!=1){
    cout << "Failed in loading configuration" << endl;
    return -1;
  }
  fclose(config);

  beta=beta_pointer[0];
  initial_state=initial_state_pointer[0];
  sigma=sigma_pointer[0];

  FILE* config2=fopen(config2_file, "r");
  if(config2==NULL){
    cout << "Error in opening config2" << endl;
    return -1;
  }
  
  // parameters loaded from config2
  double** acceptance_target=new double*[Nreplica]; // goal of optimized acceptance
  for(i=0;i<Nreplica;i++){
    acceptance_target[i]=new double[num_params];
  }
  // parameters for bisection optimization
  double error_coeff;
  double error_stop; 

  int loadConfig2_status=loadConfig2(config2, acceptance_target, Nreplica, num_params, &error_coeff, &error_stop);
  if(loadConfig2_status!=1){
    cout << "Failed in loading configuration2" << endl;
    return -1;
  }
  fclose(config2);
  
  cout << "Number of replicas: " << Nreplica << endl;
  cout << "Inverse temperatures: ";
  for(i=0;i<Nreplica;i++){
    cout << beta[i] << ", ";
  }
  cout << endl;
  cout << "Number of MC steps: " << Nstep << endl;
  cout << "Sample data: " << sample_file << endl;
  cout << "Number of bins: " << Nbin << endl;
  
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

  cout << "Initial development parameters: " << endl;
  for(i=0;i<Nreplica;i++){
    cout << "Replica " << i+1 << ": ";
    for(j=0;j<num_params;j++){
      cout << sigma[i][j] << ", ";
    }
    cout << endl;
  }
  
  Hamiltonian* h=new Hamiltonian(sample_file, Nbin);
  MT* mt=new MT();

  // see MCMC/MCMC_opt_mND_main.cpp for details
  int** opt_stage=new int*[Nreplica];
  int** accept_count=new int*[Nreplica];
  int** trial_count=new int*[Nreplica];
  double** sigma_large=new double*[Nreplica];
  double** sigma_small=new double*[Nreplica];
  for(i=0;i<Nreplica;i++){
    opt_stage[i]=new int[num_params];
    accept_count[i]=new int[num_params];
    trial_count[i]=new int[num_params];
    sigma_large[i]=new double[num_params];
    sigma_small[i]=new double[num_params];
    for(j=0;j<num_params;j++){
      accept_count[i][j]=0;
      trial_count[i][j]=0;
      opt_stage[i][j]=0;
      sigma_large[i][j]=0.0;
    }
  }

  // MC step
  char* state_string=s->print();
  int accept;
  int sigma_update;
  char index_print[32];
  int opt_complete;
  for(i=0;i<Nstep;i++){
    // MC step
    sigma_update=0;
    opt_complete=1;
    for(j=0;j<Nreplica;j++){
      for(k=0;k<num_params;k++){
	trial_count[j][k]++;
	accept=MC_Metropolis(ss[j],s_next,h,k,sigma[j][k],beta[j],mt);
	if(accept==1){
	  //replace s and s_next
	  s_buff=s_next;
	  s_next=ss[j];
	  ss[j]=s_buff;
	  accept_count[j][k]++;
	}
	sprintf(index_print, "%d-%d", j+1, k+1);
	sigma_update+=Acceptance_optimization(index_print, &sigma_small[j][k], &sigma[j][k], &sigma_large[j][k], &opt_stage[j][k], &accept_count[j][k], &trial_count[j][k], acceptance_target[j][k], error_coeff, error_stop);
      }
    }
    // if all parameters are optimized, break the while loop
    for(j=0;j<Nreplica;j++){
      for(k=0;k<num_params;k++){
	if(opt_stage[j][k]!=2){
	  opt_complete=0;
	  break;
	}
      }
      if(opt_complete==0){
	break;
      }
    }
    if(opt_complete==1){
      break;
    }
    /*
    if(sigma_update>0){
      cout << "Current development parameters: " << endl;
      for(j=0;j<Nreplica;j++){
	cout << "Replica " << j+1 << ": ";
	for(k=0;k<num_params;k++){
	  cout << sigma[j][k] << ", ";
	}
	cout << endl;
      }
      cout << endl;
      }*/
    
    
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

      if(rand<threshold){
	// accept -> exchange
	s_buff=ss[index_low];
	ss[index_low]=ss[index_high];
	ss[index_high]=s_buff;
      }else{
	// reject
      }
      index_low+=2;
      index_high+=2;
    }
  }

  cout << "--------" << endl;
  if(opt_complete==1){
    cout << "Optimization completed in " << i << " steps" << endl;
  }else{
    cout << "Optimization not completed" << endl;
  }
  
  cout << "Development parameters: " << endl;
  for(i=0;i<Nreplica;i++){
    for(j=0;j<num_params;j++){
      printf("%12.4e ", sigma[i][j]);
    }
    cout << endl;
  }
  cout << endl;

  cout << "Optimization stage: " << endl;
  for(i=0;i<Nreplica;i++){
    for(j=0;j<num_params;j++){
      printf("%1d ", opt_stage[i][j]);
    }
    cout << endl;
  }
  cout << endl;

  cout << "Number of acceptance: " << endl;
  for(i=0;i<Nreplica;i++){
    for(j=0;j<num_params;j++){
      printf("%8d ", accept_count[i][j]);
    }
    cout << endl;
  }
  cout << endl;

  
  cout << "Number of trial: " << endl;
  for(i=0;i<Nreplica;i++){
    for(j=0;j<num_params;j++){
      printf("%8d ", trial_count[i][j]);
    }
    cout << endl;
  }
  cout << endl;

  
  cout << "Acceptance rate: " << endl;
  for(i=0;i<Nreplica;i++){
    for(j=0;j<num_params;j++){
      if(trial_count[i][j]!=0){
	printf("%6.3f ", accept_count[i][j]*1.0/trial_count[i][j]);
      }else{
	printf("%6.3f ", 0);
      }
    }
    cout << endl;
  }
  cout << endl;
  
}



int loadConfig2(FILE* config, double** acceptance_target, int Nreplica, int num_params, double* error_coeff, double* error_stop){

  int buffer_length=256;
  char line[buffer_length];
  char* fgets_status;
  int sscanf_status;
  
  // line 1-Nreplica: target acceptance ratio
  int i,j;
  
  for(i=0;i<Nreplica;i++){
    for(j=0;j<num_params;j++){
      sscanf_status=fscanf(config, "%lf", &acceptance_target[i][j]);
      if(sscanf_status!=1){
	cout << "Error in parsing line " << i+1 << " of the configuration file" << endl;
	return 0;
      }
    }
    // pass to the end of the line
    fgets_status=fgets(line, buffer_length, config);
    if(fgets_status==NULL){
      cout << "Error in loading line " << i+1 << " of the configuration file" << endl;
      return 0;
    }
  }
  
  // line 1+Nreplica: error_coeff
  fgets_status=fgets(line, buffer_length, config);
  if(fgets_status==NULL){
    cout << "Error in loading line " << 1+Nreplica << " of the configuration file" << endl;
    return 0;
  }
  sscanf_status=sscanf(line, "%lf", error_coeff);
  if(sscanf_status!=1){
    cout << "Error in parsing line " << 1+Nreplica << " of the configuration file" << endl;
    return 0;
  }

  // line 3: error_stop
  fgets_status=fgets(line, buffer_length, config);
  if(fgets_status==NULL){
    cout << "Error in loading line " << 1+Nreplica << " of the configuration file" << endl;
    return 0;
  }
  sscanf_status=sscanf(line, "%lf", error_stop);
  if(sscanf_status!=1){
    cout << "Error in parsing line " << 1+Nreplica << " of the configuration file" << endl;
    return 0;
  }

  return 1;
    
}
