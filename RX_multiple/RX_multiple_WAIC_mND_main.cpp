#include "../MCMC/Model_selection_mND.hpp"
#include <iostream>
#include <string.h>

using namespace std;

int loadConfig(FILE* config, int* Nchains, char* state_file_format, int* Nstep, int* Burnin, char* sample_file, int* Nbin);

int main(int argc, char** argv){
  // initial message
  cout << "# RX_multiple WAIC calculation" << endl;
  
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

  int i,j;
  
  State* s=new State();
  char* s_print=s->print();
  int state_length=strlen(s_print);
  int num_params=s->num_params;
  
  // <---- Parameters loaded from config_file

  int Nchains; // number of Marcov chains
  int buffer_length=256;
  char* state_file_format=new char[buffer_length]; // development of the state
  int Nstep; // Number of MC steps
  int Burnin; // number of burn-in steps
  char* sample_file=new char[buffer_length]; // sample data
  int Nbin; // number of bins

  int loadConfig_status=loadConfig(config, &Nchains, state_file_format, &Nstep, &Burnin, sample_file, &Nbin);
  // ---->
  
  if(loadConfig_status!=1){
    cout << "Failed in loading configuration" << endl;
    return -1;
  }
  fclose(config);

  cout << "Number of Marcov chains: " << Nchains << endl;
  cout << "Number of MC steps: " << Nstep << endl;
  cout << "Sample data: " << sample_file << endl;
  cout << "Number of bins: " << Nbin << endl;

  // load sample data (see Hamiltonian)
  double* Xvalue=new double[Nbin];
  int* Count=new int[Nbin];
  FILE* sample=fopen(sample_file, "r");
  if(sample==NULL){
    cout << "Error in opening sample file" << endl;
    std::exit(0);
  }
  int data_count=0;
  int sscanf_status;
  char* buffer=new char[buffer_length];
  int Nsample=0;
  while(fgets(buffer, buffer_length, sample)!=NULL){
    if(buffer[0]=='#'){
      continue;
    }
    sscanf_status=sscanf(buffer, "%lf%d", &Xvalue[data_count], &Count[data_count]);
    if(sscanf_status==2){
      Nsample+=Count[data_count];
      data_count++;
    }else{
      cout << "Error in parsing sample data" << endl;
    }
  }
  if(data_count==Nbin){
    cout << "Sample data is loaded successfully" << endl;
    cout << "Number of samples: " << Nsample << endl;
  }else{
    cout << "Failed in loading sample data" << endl;
    return -1;
  }
  

  char* header=StateHeader();
  int data_length=max(strlen(header), strlen(s_print))+10; // +10 for \r, \n, \0
  FILE* state_input;
  char* file_name=new char[data_length];
  int state_count=0;
  int line_count=0;


  double* posterior_xi=new double[Nbin];
  double* log_p_xi=new double[Nbin];
  double* log_p_xi2=new double[Nbin];
  double p_xi;
  for(i=0;i<Nbin;i++){
    posterior_xi[i]=0.0;
    log_p_xi[i]=0.0;
    log_p_xi2[i]=0.0;
  }
  
  
  for(i=0;i<Nchains;i++){
    state_count=0;
    sprintf(file_name, state_file_format, i+1);
    state_input=fopen(file_name, "r");
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
	  for(j=0;j<Nbin;j++){
	    p_xi=s->probability(Xvalue[j]);
	    posterior_xi[j]+=p_xi;
	    log_p_xi[j]+=log(p_xi);
	    log_p_xi2[j]+=pow(log(p_xi),2);
	  }
	}
      }else{
	// load fail
	cout << "Failed in loading line " << line_count << ", " << load_status << endl;
      }
    } 
    fclose(state_input);

    if(state_count==Nstep){
      cout << "Load succeeded from chain " << i+1 << endl;
    }else{
      cout << "Load failed from chain " << i+1 <<endl;
      return -1;
    }
  }
  
  cout << "All loading completed" << endl;
  cout << "bin index, Xvalue, Count, posterior, variance" << endl;
  // T_n: experience loss = -1/Nsample \sum_{sample} log[posterior(sample)]
  // sample data has the histogram form,
  // T_n=-1/Nsample \sum_{i: bin index} Count[i]*log[posterior(Xvalue[i])]
  double T_n=0.0;
  double V_n=0.0;
  double T_ni, V_ni;
  for(i=0;i<Nbin;i++){
    // divide by (Nstep-Burnin)*Nchains to normalize(average)
    posterior_xi[i]/=((Nstep-Burnin)*Nchains*1.0);
    log_p_xi[i]/=((Nstep-Burnin)*Nchains*1.0);
    log_p_xi2[i]/=((Nstep-Burnin)*Nchains*1.0);
    T_n+=-Count[i]*1.0*log(posterior_xi[i])/Nsample;
    V_ni=log_p_xi2[i]-pow(log_p_xi[i],2);
    V_n+=Count[i]*V_ni;
    printf("%d, %12.4e, %d, %12.4e, %12.4e\n", i+1, Xvalue[i], Count[i], posterior_xi[i], V_ni);
  }

  printf("T_n: %12.4e\n", T_n);
  printf("V_n: %12.4e\n", V_n);
  printf("WAIC: %12.4e\n", T_n+V_n/Nsample);
  
}
int loadConfig(FILE* config, int* Nchains, char* state_file_format, int* Nstep, int* Burnin, char* sample_file, int* Nbin){
  int buffer_length=256;
  char line[buffer_length];
  char* fgets_status;
  int sscanf_status;

  // line 1: Number of chains
  fgets_status=fgets(line, buffer_length, config);
  if(fgets_status==NULL){
    cout << "Error in loading line 1 of the configuration file" << endl;
    return 0;
  }
  sscanf_status=sscanf(line, "%d", Nchains);
  if(sscanf_status!=1){
    cout << "Error in parsing line 1 of the configuration file" << endl;
    return 0;
  }

  // line 2: state_file_format
  fgets_status=fgets(line, buffer_length, config);
  if(fgets_status==NULL){
    cout << "Error in loading line 2 of the configuration file" << endl;
    return 0;
  }
  sscanf_status=sscanf(line, "%s", state_file_format);
  if(sscanf_status!=1){
    cout << "Error in parsing line 2 of the configuration file" << endl;
    return 0;
  }
  
  // line 3: Number of steps
  fgets_status=fgets(line, buffer_length, config);
  if(fgets_status==NULL){
    cout << "Error in loading line 3 of the configuration file" << endl;
    return 0;
  }
  sscanf_status=sscanf(line, "%d", Nstep);
  if(sscanf_status!=1){
    cout << "Error in parsing line 3 of the configuration file" << endl;
    return 0;
  }
  
  // line 4: Burnin steps
  fgets_status=fgets(line, buffer_length, config);
  if(fgets_status==NULL){
    cout << "Error in loading line 4 of the configuration file" << endl;
    return 0;
  }
  sscanf_status=sscanf(line, "%d", Burnin);
  if(sscanf_status!=1){
    cout << "Error in parsing line 4 of the configuration file" << endl;
    return 0;
  }
  
  // line 5: sample_file
  fgets_status=fgets(line, buffer_length, config);
  if(fgets_status==NULL){
    cout << "Error in loading line 5 of the configuration file" << endl;
    return 0;
  }
  sscanf_status=sscanf(line, "%s", sample_file);
  if(sscanf_status!=1){
    cout << "Error in parsing line 5 of the configuration file" << endl;
    return 0;
  }

  // line 6: Nbin
  fgets_status=fgets(line, buffer_length, config);
  if(fgets_status==NULL){
    cout << "Error in loading line 6 of the configuration file" << endl;
    return 0;
  }
  sscanf_status=sscanf(line, "%d", Nbin);
  if(sscanf_status!=1){
    cout << "Error in parsing line 6 of the configuration file" << endl;
    return 0;
  }
  
  return 1;
}
