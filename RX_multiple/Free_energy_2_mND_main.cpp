#include <iostream>
#include "../MCMC/Model_selection_mND.hpp"
#include <string.h>
#include <cmath>
#include "loadConfig_mND.hpp"

using namespace std;

double log_sum_exp(double a, double b);


int main(int argc, char** argv){
  cout << "# Calculate Free energy (log(Z)) by method 2" << endl;
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

  // edbH means exp(-(beta_{i+1}-beta_{i})*H(s))
  double* log_sum_edbH=new double[Nreplica-1];
  double* log_sum_e2dbH=new double[Nreplica-1];
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
  double dbH;
  for(j=0;j<Nreplica-1;j++){
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
	    dbH=-(beta[j+1]-beta[j])*h->energy(s);
	    // add exp(beta*H(s)) by log-sum-exp method
	    if(log_sum_exp_init==0){
	      log_sum_exp_init=1;
	      log_sum_edbH[j]=dbH;
	      log_sum_e2dbH[j]=2*dbH;
	    }else{
	      log_sum_edbH[j]=log_sum_exp(log_sum_edbH[j], dbH);
	      log_sum_e2dbH[j]=log_sum_exp(log_sum_e2dbH[j], 2*dbH);
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
  printf("# beta log(dZ) Var(log(dZ)) log(Z) Var(log(Z)) Delta(log(Z))\n");
  double logdZ_sum=0;
  double VarlogdZ_sum;
  double N=(Nstep-Burnin)*Nchains*1.0;
  for(j=0;j<Nreplica-1;j++){
    // logdZ=log(Z[j+1]/Z[j])=log(1/N sum(exp(-(b[j+1]-b[j])H(s[i]))))
    double logdZ=log_sum_edbH[j]-log(N);
    // logdZ2=log(1/N sum(exp(2dbH(s[i]))))
    double logdZ2=log_sum_e2dbH[j]-log(N);
    // log(Var(dZ))=log(Var(Z[j+1]/Z[j]))=log(N/N-1 (sum(exp(2bH(s[i])))/N-(sum(exp(bH(s[i])))/N)^2))
    //       =log(N/N-1)+log(exp(logZ2)-exp(2*logZ))
    //       =log(N/N-1)+logZ2+log(1-exp(2*logZ-logZ2))
    double logVardZ=log(N/(N-1))+logdZ2+log(1-exp(2.0*logdZ-logdZ2));
    // Var(log(dZ))=Var(dZ)/dZ^2=exp(log(Var(dZ))-2*log(dZ))
    double VarlogdZ=exp(logVardZ-2.0*logdZ);

    // logdZ_sum=sum_j logdZ[j]=Z[j]
    logdZ_sum+=logdZ;
    // Var(Z[j])=sum_j VarlogdZ[j]
    VarlogdZ_sum+=VarlogdZ;
    printf("%12.4e %12.4e %12.4e %12.4e %12.4e %12.4e\n", beta[j+1], logdZ, VarlogdZ, logdZ_sum, VarlogdZ_sum, sqrt(VarlogdZ_sum));
  }

  
  return 1;
}

double log_sum_exp(double a, double b){
  // calculate log(exp(a)+exp(b)) = max(a,b)+log(1+exp(-|a-b|);
  return max(a,b)+log(1+exp(-abs(a-b)));
}
