#include "../MCMC/Model_selection_mND.hpp"
#include <iostream>
#include <string.h>

using namespace std;

int loadConfig(FILE* config, int* Nchains, char* state_file_format, int* Nstep, double* min_x, double* max_x, int* num_points, int* Burnin, char* output_file);

int main(int argc, char** argv){
  // initial message
  cout << "# RX_multiple Posterior distribution calculation" << endl;
  
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
  double min_x, max_x; // x range for the posterior distribution
  int num_points; // number of points for which the distribution is calculated
  int Burnin; // number of burn-in steps
  char* output_file=new char[buffer_length]; // output file

  int loadConfig_status=loadConfig(config, &Nchains, state_file_format, &Nstep, &min_x, &max_x, &num_points, &Burnin, output_file);
  
  // ---->
  
  if(loadConfig_status!=1){
    cout << "Failed in loading configuration" << endl;
    return -1;
  }
  fclose(config);

  cout << "Number of Marcov chains: " << Nchains << endl;
  cout << "Number of MC steps: " << Nstep << endl;
  cout << "Posterior distribution range: " << min_x << " to " << max_x << endl;
  cout << "Output file: " << output_file << endl;

  char* header=StateHeader();
  int data_length=max(strlen(header), strlen(s_print))+10; // +10 for \r, \n, \0
  char* buffer=new char[data_length];
  FILE* state_input;
  char* file_name=new char[data_length];
  int state_count=0;
  int line_count=0;

  FILE* output=fopen(output_file, "w");
  fprintf(output, "# Posterior distribution\n");
  fprintf(output, "# x p(x)\n");

  double* x=new double[num_points];
  double* posterior=new double[num_points];
  for(i=0;i<num_points;i++){
    posterior[i]=0;
    x[i]=min_x+(max_x-min_x)*i/(num_points-1);
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
	  for(j=0;j<num_points;j++){
	    posterior[j]+=s->probability(x[j]);
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
  for(i=0;i<num_points;i++){
    fprintf(output, "%12.4e %12.4e\n", x[i], posterior[i]/((Nstep-Burnin)*Nchains));
  }
  
  fclose(output);
}

int loadConfig(FILE* config, int* Nchains, char* state_file_format, int* Nstep, double* min_x, double* max_x, int* num_points, int* Burnin, char* output_file){
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
  
  // line 4: min_x, max_x, num_points
  fgets_status=fgets(line, buffer_length, config);
  if(fgets_status==NULL){
    cout << "Error in loading line 4 of the configuration file" << endl;
    return 0;
  }
  sscanf_status=sscanf(line, "%lf%lf%d", min_x, max_x, num_points);
  if(sscanf_status!=3){
    cout << "Error in parsing line 4 of the configuration file" << endl;
    return 0;
  }

  // line 5: Burnin
  fgets_status=fgets(line, buffer_length, config);
  if(fgets_status==NULL){
    cout << "Error in loading line 5 of the configuration file" << endl;
    return 0;
  }
  sscanf_status=sscanf(line, "%d", Burnin);
  if(sscanf_status!=1){
    cout << "Error in parsing line 5 of the configuration file" << endl;
    return 0;
  }

  // line 6: output_file
  fgets_status=fgets(line, buffer_length, config);
  if(fgets_status==NULL){
    cout << "Error in loading line 6 of the configuration file" << endl;
    return 0;
  }
  sscanf_status=sscanf(line, "%s", output_file);
  if(sscanf_status!=1){
    cout << "Error in parsing line 6 of the configuration file" << endl;
    return 0;
  }
  return 1;
}
