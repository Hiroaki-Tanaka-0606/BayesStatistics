#include "loadConfig_mND.hpp"

#include <iostream>

using namespace std;

int loadConfig(FILE* config, int* Nchains, int* Nreplica, double** beta_pointer, int* Nstep, char* input_file_format, char* state_file_format, char* accept_file_format, char* exchange_file_format, char* log_file_format, char* sample_file, int* Nbin, double*** sigma_pointer, int num_params, char* shell_script, char* input2_file_format, int* Burnin, char* average_file_format, char* variance_file_format, char* log2_file_format, char* value_format, char* input3_file_format, char* psrf_file_format, char* log3_file_format){
  int buffer_length=1024;
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

  // line 2: Nreplica
  fgets_status=fgets(line, buffer_length, config);
  if(fgets_status==NULL){
    cout << "Error in loading line 2 of the configuration file" << endl;
    return 0;
  }
  sscanf_status=sscanf(line, "%d", Nreplica);
  if(sscanf_status!=1){
    cout << "Error in parsing line 2 of the configuration file" << endl;
    return 0;
  }

  // line 3-(Nreplica+2): inverse temperature beta
  int i;
  double* beta=new double[*Nreplica];
  beta_pointer[0]=beta;
  for(i=0;i<*Nreplica;i++){
    fgets_status=fgets(line, buffer_length, config);
    if(fgets_status==NULL){
      cout << "Error in loading line " << 3+i << " of the configuration file" << endl;
      return 0;
    }
    sscanf_status=sscanf(line, "%lf", &beta[i]);
    if(sscanf_status!=1){
      cout << "Error in parsing line " << 3+i << " of the configuration file" << endl;
      return 0;
    }
  }

  // line 3+Nreplica: Number of MC steps
  fgets_status=fgets(line, buffer_length, config);
  if(fgets_status==NULL){
    cout << "Error in loading line " << 3+*Nreplica << " of the configuration file" << endl;
    return 0;
  }
  sscanf_status=sscanf(line, "%d", Nstep);
  if(sscanf_status!=1){
    cout << "Error in parsing line " << 3+*Nreplica << " of the configuration file" << endl;
    return 0;
  }

  // line 4+Nreplica: input_file_format
  fgets_status=fgets(line, buffer_length, config);
  if(fgets_status==NULL){
    cout << "Error in loading line " << 4+*Nreplica << " of the configuration file" << endl;
    return 0;
  }
  sscanf_status=sscanf(line, "%s", input_file_format);
  if(sscanf_status!=1){
    cout << "Error in parsing line " << 4+*Nreplica << " of the configuration file" << endl;
    return 0;
  }
  
  // line 5+Nreplica: state_file_format
  fgets_status=fgets(line, buffer_length, config);
  if(fgets_status==NULL){
    cout << "Error in loading line " << 5+*Nreplica << " of the configuration file" << endl;
    return 0;
  }
  sscanf_status=sscanf(line, "%s", state_file_format);
  if(sscanf_status!=1){
    cout << "Error in parsing line " << 5+*Nreplica << " of the configuration file" << endl;
    return 0;
  }

  // line 6+Nreplica: accept_file_format
  fgets_status=fgets(line, buffer_length, config);
  if(fgets_status==NULL){
    cout << "Error in loading line " << 6+*Nreplica << " of the configuration file" << endl;
    return 0;
  }
  sscanf_status=sscanf(line, "%s", accept_file_format);
  if(sscanf_status!=1){
    cout << "Error in parsing line " << 6+*Nreplica << " of the configuration file" << endl;
    return 0;
  }

  // line 7+Nreplica: exchange_file_format
  fgets_status=fgets(line, buffer_length, config);
  if(fgets_status==NULL){
    cout << "Error in loading line " << 7+*Nreplica << " of the configuration file" << endl;
    return 0;
  }
  sscanf_status=sscanf(line, "%s", exchange_file_format);
  if(sscanf_status!=1){
    cout << "Error in parsing line " << 7+*Nreplica << " of the configuration file" << endl;
    return 0;
  }
  
  // line 8+Nreplica: log_file_format
  fgets_status=fgets(line, buffer_length, config);
  if(fgets_status==NULL){
    cout << "Error in loading line " << 8+*Nreplica << " of the configuration file" << endl;
    return 0;
  }
  sscanf_status=sscanf(line, "%s", log_file_format);
  if(sscanf_status!=1){
    cout << "Error in parsing line " << 8+*Nreplica << " of the configuration file" << endl;
    return 0;
  }

  
  // line 9+Nreplica: sample_file
  fgets_status=fgets(line, buffer_length, config);
  if(fgets_status==NULL){
    cout << "Error in loading line " << 9+*Nreplica << " of the configuration file" << endl;
    return 0;
  }
  sscanf_status=sscanf(line, "%s", sample_file);
  if(sscanf_status!=1){
    cout << "Error in parsing line " << 9+*Nreplica << " of the configuration file" << endl;
    return 0;
  }

  
  // line 10+Nreplica: number of bins
  fgets_status=fgets(line, buffer_length, config);
  if(fgets_status==NULL){
    cout << "Error in loading line " << 10+*Nreplica << " of the configuration file" << endl;
    return 0;
  }
  sscanf_status=sscanf(line, "%d", Nbin);
  if(sscanf_status!=1){
    cout << "Error in parsing line " << 10+*Nreplica << " of the configuration file" << endl;
    return 0;
  }
 
  // line (11+Nreplica)-(10+2*Nreplica): development parameters
  double** sigma=new double*[*Nreplica];
  sigma_pointer[0]=sigma;
  for(i=0;i<*Nreplica;i++){
    int j;
    sigma[i]=new double[num_params];
    for(j=0;j<num_params;j++){
      sscanf_status=fscanf(config, "%lf", &sigma[i][j]);
      if(sscanf_status!=1){
	cout << "Error in parsing line " << 11+*Nreplica << " of the configuration file" << endl;
	return 0;
      }
    }
    // pass to the end of the line
    fgets_status=fgets(line, buffer_length, config);
    if(fgets_status==NULL){
      cout << "Error in loading line " <<  11+*Nreplica << " of the configuration file" << endl;
	return 0;
    }
  }

  // line 11+2*Nreplica: shell_script
  fgets_status=fgets(line, buffer_length, config);
  if(fgets_status==NULL){
    cout << "Error in loading line " << 11+2*(*Nreplica) << " of the configuration file" << endl;
    return 0;
  }
  sscanf_status=sscanf(line, "%s", shell_script);
  if(sscanf_status!=1){
    cout << "Error in parsing line " << 11+2*(*Nreplica) << " of the configuration file" << endl;
    return 0;
  }

  // line 12+2*Nreplica: input2_file_format
  fgets_status=fgets(line, buffer_length, config);
  if(fgets_status==NULL){
    cout << "Error in loading line " << 12+2*(*Nreplica) << " of the configuration file" << endl;
    return 0;
  }
  sscanf_status=sscanf(line, "%s", input2_file_format);
  if(sscanf_status!=1){
    cout << "Error in parsing line " << 12+2*(*Nreplica) << " of the configuration file" << endl;
    return 0;
  }

  // line 13+2*Nreplica: Burnin
  fgets_status=fgets(line, buffer_length, config);
  if(fgets_status==NULL){
    cout << "Error in loading line " << 13+2*(*Nreplica) << " of the configuration file" << endl;
    return 0;
  }
  sscanf_status=sscanf(line, "%d", Burnin);
  if(sscanf_status!=1){
    cout << "Error in parsing line " << 13+2*(*Nreplica) << " of the configuration file" << endl;
    return 0;
  }

  // line 14+2*Nreplica: average file format
  fgets_status=fgets(line, buffer_length, config);
  if(fgets_status==NULL){
    cout << "Error in loading line " << 14+2*(*Nreplica) << " of the configuration file" << endl;
    return 0;
  }
  sscanf_status=sscanf(line, "%s", average_file_format);
  if(sscanf_status!=1){
    cout << "Error in parsing line " << 14+2*(*Nreplica) << " of the configuration file" << endl;
    return 0;
  }

  // line 15+2*Nreplica: variance file format
  fgets_status=fgets(line, buffer_length, config);
  if(fgets_status==NULL){
    cout << "Error in loading line " << 15+2*(*Nreplica) << " of the configuration file" << endl;
    return 0;
  }
  sscanf_status=sscanf(line, "%s", variance_file_format);
  if(sscanf_status!=1){
    cout << "Error in parsing line " << 15+2*(*Nreplica) << " of the configuration file" << endl;
    return 0;
  }

  // line 16+2*Nreplica: log file format
  fgets_status=fgets(line, buffer_length, config);
  if(fgets_status==NULL){
    cout << "Error in loading line " << 16+2*(*Nreplica) << " of the configuration file" << endl;
    return 0;
  }
  sscanf_status=sscanf(line, "%s", log2_file_format);
  if(sscanf_status!=1){
    cout << "Error in parsing line " << 16+2*(*Nreplica) << " of the configuration file" << endl;
    return 0;
  }

  // line 17+2*Nreplica: value format
  fgets_status=fgets(line, buffer_length, config);
  if(fgets_status==NULL){
    cout << "Error in loading line " << 17+2*(*Nreplica) << " of the configuration file" << endl;
    return 0;
  }
  sscanf_status=sscanf(line, "%s", value_format);
  if(sscanf_status!=1){
    cout << "Error in parsing line " << 17+2*(*Nreplica) << " of the configuration file" << endl;
    return 0;
  }

  
  // line 18+2*Nreplica: input3 file name
  fgets_status=fgets(line, buffer_length, config);
  if(fgets_status==NULL){
    cout << "Error in loading line " << 18+2*(*Nreplica) << " of the configuration file" << endl;
    return 0;
  }
  sscanf_status=sscanf(line, "%s", input3_file_format);
  if(sscanf_status!=1){
    cout << "Error in parsing line " << 18+2*(*Nreplica) << " of the configuration file" << endl;
    return 0;
  }

  // line 19+2*Nreplica: output file format
  fgets_status=fgets(line, buffer_length, config);
  if(fgets_status==NULL){
    cout << "Error in loading line " << 19+2*(*Nreplica) << " of the configuration file" << endl;
    return 0;
  }
  sscanf_status=sscanf(line, "%s", psrf_file_format);
  if(sscanf_status!=1){
    cout << "Error in parsing line " << 19+2*(*Nreplica) << " of the configuration file" << endl;
    return 0;
  }

  // line 20+2*Nreplica: log file
  fgets_status=fgets(line, buffer_length, config);
  if(fgets_status==NULL){
    cout << "Error in loading line " << 20+2*(*Nreplica) << " of the configuration file" << endl;
    return 0;
  }
  sscanf_status=sscanf(line, "%s", log3_file_format);
  if(sscanf_status!=1){
    cout << "Error in parsing line " << 20+2*(*Nreplica) << " of the configuration file" << endl;
    return 0;
  }
  
  return 1;
}
