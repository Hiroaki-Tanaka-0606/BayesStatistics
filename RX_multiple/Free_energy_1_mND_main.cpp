#include <iostream>
#include "../MCMC/Model_selection_mND.hpp"
#include <string.h>
#include <cmath>
#include "loadConfig_mND.hpp"

using namespace std;

double log_sum_exp(double a, double b);

int main(int argc, char** argv){
  cout << "# Calculate Free energy (log(Z)) by method 1" << endl;
  if(argc<2){
    printf("Usage: %s config_file\n", argv[0]);
    return 0;
  }

  State* s=new State();
  int num_params=s->num_params;
  
  
  // load configuration
  char* config_file=argv[1];
  cout << "# Configuration file: " << config_file << endl;
  FILE* config=fopen(config_file, "r");
  if(config==NULL){
    cout << "Error in opening the config file" << endl;
    return -1;
  }
  // input file is the same as RX_multiple.o, so some parameters are not necessary for this calculation
  
  // configuration for RXMC.o
  int Nchains; // number of Marcov chains
  int Nreplica;
  double* beta; // inverse temperature
  int Nstep; // Number of MC steps
  int buffer_length=256;
  char* input_file_format=new char[buffer_length]; // name of the input file, including %d (chain index)
  char* state_file_format=new char[buffer_length]; // name of the output file for development of the state, including %d (chain index), %%d (replica index)
  char* accept_file_format=new char[buffer_length]; // name of the output file for log of acceptance, including %d (chain index), %%d (replica index)
  char* exchange_file_format=new char[buffer_length]; // name of the output file for log of exchange, including %d (chain index)
  char* log_file_format=new char[buffer_length]; // name of the log (STDOUT), including %d (chain index)
  char* sample_file=new char[buffer_length]; // sample file
  int Nbin; // number of bins
  double** sigma; // development parameter
  char* shell_script=new char[buffer_length]; // name of the shell script file

  // configuration for Average_analysis.o
  // number of parameters, buffer size for fgets, state file[i], Nstep, are required
  char* input2_file_format=new char[buffer_length]; // name of the input 
  int Burnin; // Burn-in steps
  char* average_file_format=new char[buffer_length]; // name of the output file for average
  char* variance_file_format=new char[buffer_length]; // name of the output file for variance
  char* log2_file_format=new char[buffer_length]; // name of the log
  char* value_format=new char[buffer_length]; // output format for a value
  
  // configuration for PSRF.o
  // number of parameters, buffer size for fgets, number of MCs, , average file format, variance file format, Nstep-Burnin, are required
  char* input3_file_format=new char[buffer_length]; // name of the input
  char* psrf_file_format=new char[buffer_length]; // name of the output for psrf, including %d (parameter index, not chain index)
  char* log3_file_format=new char[buffer_length]; // name of the log

  // pointer for loadConfig arguments
  double** beta_pointer=new double*[1];
  double*** sigma_pointer=new double**[1];
  int loadConfig_status=loadConfig(config, &Nchains, &Nreplica, beta_pointer, &Nstep, input_file_format, state_file_format, accept_file_format, exchange_file_format, log_file_format, sample_file, &Nbin, sigma_pointer, num_params, shell_script, input2_file_format, &Burnin, average_file_format, variance_file_format, log2_file_format, value_format, input3_file_format, psrf_file_format, log3_file_format);
  if(loadConfig_status!=1){
    cout << "Failed in loading configuration" << endl;
    return -1;
  }

  beta=beta_pointer[0];
  int i,j;
  
  cout << "# Number of Marcov chains: " << Nchains << endl;
  cout << "# Number of replicas: " << Nreplica << endl;
  cout << "# Inverse temperatures: ";
  for(i=0;i<Nreplica;i++){
    cout << beta[i] << ", ";
  }
  cout << endl;
  cout << "# Number of MC steps: " << Nstep << endl;
  cout << "# Burn-in steps: " << Burnin << endl;

  double* log_sum_ebH=new double[Nreplica];
  double* log_sum_e2bH=new double[Nreplica];
  int log_sum_exp_init;
  Hamiltonian* h=new Hamiltonian(sample_file, Nbin);

  // i -> chain index, j -> replica index
  char* state_file=new char[buffer_length];
  char* state_file_buffer=new char[buffer_length];
  FILE* input;

  char* header=StateHeader();
  char* state_print=s->print();  
  int state_buffer_length=max(strlen(header), strlen(state_print))+10; // +10 for \r, \n, \0
  char* state_buffer=new char[state_buffer_length+10]; // +10 for \r, \n, \0

  int state_load_status;
  int state_count;
  double bH;
  for(j=0;j<Nreplica;j++){
    log_sum_exp_init=0;
    for(i=0;i<Nchains;i++){
      sprintf(state_file_buffer, state_file_format, i+1);
      sprintf(state_file, state_file_buffer, j+1);
      input=fopen(state_file, "r");
      if(input==NULL){
	cout << "Error in opening state file" << endl;
	return -1;
      }
      state_count=0;
      while(fgets(state_buffer, state_buffer_length+9, input)!=NULL){
	if(state_buffer[0]=='#'){
	  continue;
	}
	state_load_status=s->load(state_buffer);
	if(state_load_status==1){
	  // load succeeded
	  state_count++;
	  if(state_count>Burnin){
	    // calculate beta*H(s)
	    bH=beta[j]*h->energy(s);
	    // add exp(beta*H(s)) by log-sum-exp method
	    if(log_sum_exp_init==0){
	      log_sum_exp_init=1;
	      log_sum_ebH[j]=bH;
	      log_sum_e2bH[j]=2*bH;
	    }else{
	      log_sum_ebH[j]=log_sum_exp(log_sum_ebH[j], bH);
	      log_sum_e2bH[j]=log_sum_exp(log_sum_e2bH[j], 2*bH);
	    }
	    if(!isfinite(bH)){
	      cout << "!!NAN value appear!!" << endl;
	      cout << s->print() << endl;
	    }
	  }
	}else{
	  // load failed
	  cout << "Failed in loading a state" << endl;
	  cout << state_buffer << endl;
	}
      }
      
      fclose(input);
      if(state_count==Nstep){
	printf("# Succeeded in loading states from chain %d replica %d\n", i+1, j+1);
      }else{
	printf("Failed in loading states from chain %d replica %d\n", i+1, j+1);
	return -1;
      }	
    }
  }

  // calculate log(Z)
  printf("# beta log(Z) Var(log(Z)) Delta(log(Z))\n");
  double N=(Nstep-Burnin)*Nchains*1.0;
  for(j=0;j<Nreplica;j++){
    // Zp=1/Z
    // logZp=log(1/N sum(exp(bH(s[i]))))
    double logZp=log_sum_ebH[j]-log(N);
    // logZp2=log(1/N sum(exp(2bH(s[i]))))
    double logZp2=log_sum_e2bH[j]-log(N);
    // log(Var(Zp))=log(N/N-1 (sum(exp(2bH(s[i])))/N-(sum(exp(bH(s[i])))/N)^2))
    //       =log(N/N-1)+log(exp(logZp2)-exp(2*logpZ))
    //       =log(N/N-1)+logZp2+log(1-exp(2*logZp-logZp2))
    double logVarZp=log(N/(N-1))+logZp2+log(1-exp(2.0*logZp-logZp2));
    // logZ=-logZp
    double logZ=-logZp;
    // log(VarZ)=log(VarZp/Zp^4)=logVarZp-4logZp
    double logVarZ=logVarZp-4.0*logZp;
    // Var(logZ)=VarZ/Z^2=exp(log(VarZ)-2*logZ)
    double VarlogZ=exp(logVarZ-2.0*logZ);
    printf("%12.4e %12.4e %12.4e %12.4e\n", beta[j], logZ, VarlogZ, sqrt(VarlogZ));
  }

  
  return 1;
}

double log_sum_exp(double a, double b){
  // calculate log(exp(a)+exp(b)) = max(a,b)+log(1+exp(-|a-b|);
  return max(a,b)+log(1+exp(-abs(a-b)));
}
