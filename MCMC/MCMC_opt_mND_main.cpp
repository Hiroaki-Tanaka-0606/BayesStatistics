#include "Model_selection_mND.hpp"
#include "MC_Metropolis.hpp"
#include "../randomDistribution/mersenne_twister.hpp"
#include <iostream>
#include <string.h>

// <---- ----> means an area to be changed according to the model used in the MCMC calculation

using namespace std;

int loadConfig2(FILE* config2, double* acceptance_target, int num_params, double* sigma_coeff, double* sigma_stop);

int Acceptance_optimization(int parm_index, double* sigma_small, double* sigma, double* sigma_large, int* stage, int* accept, int* trial, double target, double sigma_coeff, double sigma_stop);

int main(int argc, char** argv){
  // initial message
  cout << "MCMC acceptance optimization calculation" << endl;
  
  if(argc<3){
    printf("Usage: %s config_file config2_file\n", argv[0]);
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
  char* s_print=s->print();
  int state_length=strlen(s_print);
  int num_params=s->num_params;

  int i,j;

  // <---- Parameters loaded from config_file
  
  double beta; // inverse temperature
  int Nstep; // Number of MC steps
  int buffer_length=256;
  char* state_file=new char[buffer_length]; // development of the state by MC steps
  char* accept_file=new char[buffer_length]; // log of acceptance
  char* sample_file=new char[buffer_length]; // sample data
  int Nbin; // number of bins
  double* sigma=new double[num_params]; // development parameter for each parameter
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

  FILE* config2=fopen(config2_file, "r");
  if(config2==NULL){
    cout << "Error in opening config2" << endl;
    return -1;
  }

  // parameters loaded from config2
  double* acceptance_target=new double[num_params]; // goal of optimized acceptance
  // parameters for bisection optimization
  double error_coeff;
  double error_stop; 

  int loadConfig2_status=loadConfig2(config2, acceptance_target, num_params, &error_coeff, &error_stop);
  if(loadConfig2_status!=1){
    cout << "Failed in loading configuration2" << endl;
    return -1;
  }
  fclose(config2);
  
  
  cout << "Inverse temperature: " << beta << endl;
  cout << "Number of MC steps: " << Nstep << endl;
  cout << "Sample data: " << sample_file << endl;
  cout << "Number of bins: " << Nbin << endl;

  // load initial state
  int loadInit_status=s->load(initial_state);
  if(loadInit_status==1){
    cout << "Initial state is loaded" << endl;
  }else{
    cout << "Initial state is not loaded: error " << loadInit_status << endl;
    s=new State();
  }
  cout << "Initial state: " << s->print() << endl;

  cout << "Initial development parameters: ";
  for(i=0;i<num_params;i++){
    cout << sigma[i] << ", ";
  }
  cout << endl;

  cout << "Acceptance target: " ;
  for(i=0;i<num_params;i++){
    cout << acceptance_target[i] << ", ";
  }
  cout << endl;

  State* s_next=new State();
  State* s_buff;

  Hamiltonian* h=new Hamiltonian(sample_file, Nbin);
  MT* mt=new MT();

  // opt_stage[i]
  // 0 -> find x(t) < target < x(0)=1
  // 1 -> bisection
  // 2 -> optimization succeeded
  int* opt_stage=new int[num_params];
  int* accept_count=new int[num_params];
  int* trial_count=new int[num_params];
  // F(sigma_large) > target
  // F(sigma_small) < target
  // sigma=average(sigma_large, sigma_small);
  double* sigma_large=new double[num_params];
  double* sigma_small=new double[num_params];
  for(i=0;i<num_params;i++){
    accept_count[i]=0;
    trial_count[i]=0;
    opt_stage[i]=0;
    sigma_large[i]=0.0;
  }
  
  char* state_string=s->print();
  int accept;
  int sigma_update;
  for(i=0;i<Nstep;i++){
    sigma_update=0;
    for(j=0;j<num_params;j++){
      trial_count[j]++;
      accept=MC_Metropolis(s,s_next,h,j,sigma[j],beta,mt);
      if(accept==1){
	//replace s and s_next
	s_buff=s_next;
	s_next=s;
	s=s_buff;
	accept_count[j]++;
      }
      sigma_update+=Acceptance_optimization(j, &sigma_small[j], &sigma[j], &sigma_large[j], &opt_stage[j], &accept_count[j], &trial_count[j], acceptance_target[j], error_coeff, error_stop);
    }
    if(sigma_update>0){
      cout << "Current development parameters: ";
      for(j=0;j<num_params;j++){
	cout << sigma[j] << ", ";
      }
      cout << endl << endl;
    }
  }

  // output accept count at final condition
  cout << "--------" << endl;
  cout << "Development parameters:" << endl;
  for(i=0;i<num_params;i++){
    cout << sigma[i] << " ";
  }
  cout << endl;
  
  cout << "Optimization stage: ";
  for(i=0;i<num_params;i++){
    cout << opt_stage[i] << ", ";
  }
  cout << endl;
  
  cout << "Number of acceptance: ";
  for(i=0;i<num_params;i++){
    cout << accept_count[i] << ", ";
  }
  cout << endl;
  
  cout << "Number of trial: ";
  for(i=0;i<num_params;i++){
    cout << trial_count[i] << ", ";
  }
  cout << endl;
  
  cout << "Acceptance rate: ";
  for(i=0;i<num_params;i++){
    cout << accept_count[i]*1.0/trial_count[i] << ", ";
  }
  cout << endl;

}

// f(sigma): acceptance ratio
// F(sigma)=accept/trial: simulated acceptance ratio
// considering error sqrt(trial), determine whether f(x(t)) is larger than target or not
// F(x(t))+error_coeff*error < target -> f(x(t)) is smaller than target
// F(x(t))-error_coeff*error > target -> f(x(t)) is larger than target
// error=1/sqrt(trial)
// error_coeff: coefficient (input)
// development of arguemnt sigma by bisection
// stop when error < error_stop (input)
int Acceptance_optimization(int parm_index, double* sigma_small, double* sigma, double* sigma_large, int* stage, int* accept, int* trial, double target, double error_coeff, double error_stop){

  double Fsigma=(*accept)*1.0/(*trial);
  double error=1.0/sqrt(*trial);
  switch(*stage){
  case 0:
    // find sigma s.t. F(sigma) < target < x(0)=1
    if(Fsigma + error_coeff*error < target){
      // sigma is found, start bisection
      cout << "Bisection optimizaztion start for parameter " << parm_index+1 << endl;
      cout << *accept << " / " << *trial << " +- " << error << " < " << target << endl;
      *sigma_small=*sigma;
      // sigma_large is determined in the previous step ( or initial value 0)
      *sigma=(*sigma_small+*sigma_large)/2.0;
      *stage=1;
      *accept=0;
      *trial=0;
      return 1;
    }
    if(Fsigma - error_coeff*error > target){
      // sigma is found to be too small to satisfy F(sigma) < target
      cout << "Update parameter " << parm_index+1 << endl;
      cout << *accept << " / " << *trial << " +- " << error << " > " << target << endl;
      *sigma_large=*sigma;
      *sigma*=2;
      *accept=0;
      *trial=0;
      return 1;
    }
    break;
  case 1:
    // bisection optimization
    if(Fsigma + error_coeff*error < target){
      // F(sigma_small), F(sigma) < target < F(sigma_large)
      cout << "Bisection update for parameter " << parm_index+1 << endl;
      cout << *accept << " / " << *trial << " +- " << error << " < " << target << endl;
      *sigma_small=*sigma;
      *sigma=(*sigma_small+*sigma_large)/2.0;
      *accept=0;
      *trial=0;
      return 1;
    }
    if(Fsigma - error_coeff*error > target){
      // F(sigma_small) < target < F(sigma), F(sigma_large)
      cout << "Bisection update for parameter " << parm_index+1 << endl;
      cout << *accept << " / " << *trial << " +- " << error << " > " << target << endl;
      *sigma_large=*sigma;
      *sigma=(*sigma_small+*sigma_large)/2.0;
      *accept=0;
      *trial=0;
      return 1;
    }
    if(error<error_stop){
      // F(sigma) ~ target
      cout << "Bisection optimization ended for parameter " << parm_index+1 << endl;
      cout << *accept << " / " << *trial << " +- " << error << " ~ " << target << endl;
      *stage=2;
      return 0;
    }
    break;
  case 2:
    // optimized
    break;
  }
  return 0;
}
int loadConfig2(FILE* config, double* acceptance_target, int num_params, double* error_coeff, double* error_stop){

  int buffer_length=256;
  char line[buffer_length];
  char* fgets_status;
  int sscanf_status;
  
  // line 1: target acceptance ratio
  int i;
  for(i=0;i<num_params;i++){
    sscanf_status=fscanf(config, "%lf", &acceptance_target[i]);
    if(sscanf_status!=1){
      cout << "Error in parsing line 1 of the configuration file" << endl;
      return 0;
    }
  }
  // pass to the end of line 1
  fgets_status=fgets(line, buffer_length, config);
  if(fgets_status==NULL){
    cout << "Error in loading line 1 of the configuration file" << endl;
    return 0;
  }
  
  // line 2: error_coeff
  fgets_status=fgets(line, buffer_length, config);
  if(fgets_status==NULL){
    cout << "Error in loading line 2 of the configuration file" << endl;
    return 0;
  }
  sscanf_status=sscanf(line, "%lf", error_coeff);
  if(sscanf_status!=1){
    cout << "Error in parsing line 2 of the configuration file" << endl;
    return 0;
  }

  // line 3: error_stop
  fgets_status=fgets(line, buffer_length, config);
  if(fgets_status==NULL){
    cout << "Error in loading line 3 of the configuration file" << endl;
    return 0;
  }
  sscanf_status=sscanf(line, "%lf", error_stop);
  if(sscanf_status!=1){
    cout << "Error in parsing line 3 of the configuration file" << endl;
    return 0;
  }

  return 1;
    
}
