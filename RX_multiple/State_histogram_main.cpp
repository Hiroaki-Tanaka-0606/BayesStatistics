#include "../MCMC/Model_selection_mND.hpp"
#include <iostream>
#include <string.h>

using namespace std;

int loadConfig(FILE* config, int* Nchains, char* state_file_format, int* Nstep, double* min_mu, double* max_mu, int* num_points_mu, double* min_a, double* max_a, int* num_points_a, int* Burnin, char* output_file);

int main(int argc, char** argv){
  // initial message
  cout << "# RX_multiple state histogram calculation" << endl;
  
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
	int NND=s->NND;
  
  // <---- Parameters loaded from config_file

  int Nchains; // number of Marcov chains
  int buffer_length=256;
  char* state_file_format=new char[buffer_length]; // development of the state
  int Nstep; // Number of MC steps
  double min_mu, max_mu; // mu range of histogram
  int num_points_mu; // number of points for mu direction
	double min_a, max_a; // a range of histogram
	int num_points_a; // number of points for a direction
  int Burnin; // number of burn-in steps
  char* output_file=new char[buffer_length]; // output file

  int loadConfig_status=loadConfig(config, &Nchains, state_file_format, &Nstep, &min_mu, &max_mu, &num_points_mu, &min_a, &max_a, &num_points_a, &Burnin, output_file);
  
  // ---->
  
  if(loadConfig_status!=1){
    cout << "Failed in loading configuration" << endl;
    return -1;
  }
  fclose(config);

  cout << "Number of Marcov chains: " << Nchains << endl;
  cout << "Number of MC steps: " << Nstep << endl;
  cout << "mu range: " << min_mu << " to " << max_mu << endl;
	cout << "a range: " << min_a << " to " << max_a << endl;
  cout << "Output file: " << output_file << endl;

  char* header=StateHeader();
  int data_length=max(strlen(header), strlen(s_print))+10; // +10 for \r, \n, \0
  char* buffer=new char[data_length];
  FILE* state_input;
  char* file_name=new char[data_length];
  int state_count=0;
  int line_count=0;

  FILE* output=fopen(output_file, "w");
  fprintf(output, "# State histogram\n");
  fprintf(output, "# mu_center sigma_center count\n");

	// histogram
	// [min_mu, min_mu+dmu], [min_mu+dmu, min_mu+2*dmu], ... , [max_mu-dmu, max_mu]
	// dmu=(max_mu-min_mu)/num_points_mu
	double dmu=(max_mu-min_mu)/num_points_mu;
	double da=(max_a-min_a)/num_points_a;
	int** histogram=new int*[num_points_mu];
	for(i=0;i<num_points_mu;i++){
		histogram[i]=new int[num_points_a];
		for(j=0;j<num_points_a;j++){
			histogram[i][j]=0;
		}
	}
  
  for(i=0;i<Nchains;i++){
    state_count=0;
    line_count=0;
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
					for(j=0;j<NND;j++){
						double mu_j=s->mu[j];
						double a_j=(s->a[j+1])-(s->a[j]);
						int mu_index=floor((mu_j-min_mu)/dmu);
						int a_index=floor((a_j-min_a)/da);
						if(0<=mu_index && mu_index<num_points_mu &&
							 0<=a_index && a_index<num_points_a){
							histogram[mu_index][a_index]++;
						}
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
	int total_histogram=0;
  for(i=0;i<num_points_mu;i++){
		for(j=0;j<num_points_a;j++){
			double mu_center=min_mu+(i+0.5)*dmu;
			double a_center=min_a+(j+0.5)*da;
			fprintf(output, "%12.4e %12.4e %d\n", mu_center, a_center, histogram[i][j]);
			total_histogram+=histogram[i][j];
		}
		fprintf(output, "\n");
  }
  cout << total_histogram << " data are in the histogram." << endl;
  fclose(output);
}

int loadConfig(FILE* config, int* Nchains, char* state_file_format, int* Nstep, double* min_mu, double* max_mu, int* num_points_mu, double* min_a, double* max_a, int* num_points_a, int* Burnin, char* output_file){
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
  
  // line 4: min_mu, max_mu, num_points_mu
  fgets_status=fgets(line, buffer_length, config);
  if(fgets_status==NULL){
    cout << "Error in loading line 4 of the configuration file" << endl;
    return 0;
  }
  sscanf_status=sscanf(line, "%lf%lf%d", min_mu, max_mu, num_points_mu);
  if(sscanf_status!=3){
    cout << "Error in parsing line 4 of the configuration file" << endl;
    return 0;
  }
  
  // line 5: min_a, max_a, num_points_a
  fgets_status=fgets(line, buffer_length, config);
  if(fgets_status==NULL){
    cout << "Error in loading line 5 of the configuration file" << endl;
    return 0;
  }
  sscanf_status=sscanf(line, "%lf%lf%d", min_a, max_a, num_points_a);
  if(sscanf_status!=3){
    cout << "Error in parsing line 5 of the configuration file" << endl;
    return 0;
  }
	
  // line 6: Burnin
  fgets_status=fgets(line, buffer_length, config);
  if(fgets_status==NULL){
    cout << "Error in loading line 6 of the configuration file" << endl;
    return 0;
  }
  sscanf_status=sscanf(line, "%d", Burnin);
  if(sscanf_status!=1){
    cout << "Error in parsing line 6 of the configuration file" << endl;
    return 0;
  }

  // line 7: output_file
  fgets_status=fgets(line, buffer_length, config);
  if(fgets_status==NULL){
    cout << "Error in loading line 7 of the configuration file" << endl;
    return 0;
  }
  sscanf_status=sscanf(line, "%s", output_file);
  if(sscanf_status!=1){
    cout << "Error in parsing line 7 of the configuration file" << endl;
    return 0;
  }
  return 1;
}
