#include "mixedNormalDistribution.hpp"
#include <iostream>
#include <string.h>

// generate mixed Normal Distribution
int main(int argc, char** argv){
  if(argc<2){
    printf("Usage: %s configuration_file\n", argv[0]);
    return 0;
  }
  cout << "# mixed Normal Distribution" << endl;

  int fileName_length=256;
  char config_file[fileName_length];
  strcpy(config_file, argv[1]);
  cout << "# Configuration file: " << config_file << endl;

  //load configuration
  int num_samples, N;
  FILE* config=fopen(config_file, "r");

  if(config==NULL){
    cout << "Error in opening the configuration file" << endl;
    return 0;
  }

  int buffer_length=256;
  char line[buffer_length];
  char* fgets_status;
  int sscanf_status;

  //line 1: number of samples
  fgets_status=fgets(line, buffer_length, config);
  if(fgets_status==NULL){
    cout << "Error in loading line 1 of the configuration file" << endl;
    return 0;
  }
  sscanf_status=sscanf(line, "%d", &num_samples);
  if(sscanf_status!=1){
    cout << "Error in parsing line 1 of the configuration file" << endl;
    return 0;
  }

  //line 2: number of NDs
  fgets_status=fgets(line, buffer_length, config);
  if(fgets_status==NULL){
    cout << "Error in loading line 2 of the configuration file" << endl;
    return 0;
  }
  sscanf_status=sscanf(line, "%d", &N);
  if(sscanf_status!=1){
    cout << "Error in parsing line 2 of the configuration file" << endl;
    return 0;
  }

  cout << "# Number of samples: " << num_samples << endl;
  cout << "# Number of normal distributions: " << N << endl;
  int i;
  double** conf=new double*[N];
  for(i=0; i<N; i++){
    //lines 3 to N+2: ND configuration
    conf[i]=new double[3];

    fgets_status=fgets(line, buffer_length, config);
    if(fgets_status==NULL){
      cout << "Error in loading line " << (3+i) << " of the configuration file" << endl;
      return 0;
    }
    sscanf_status=sscanf(line, "%lf%lf%lf", &conf[i][0], &conf[i][1], &conf[i][2]);
    if(sscanf_status!=3){
      cout << "Error in parsing line " << (3+i) << " of the configuration file" << endl;
      return 0;
    }
    printf("# ND[%d]: sigma=%f, mu=%f, weight=%f\n", i, conf[i][0], conf[i][1], conf[i][2]);
  }

  mND* mnd=new mND(N, conf);
  for(i=0; i<num_samples; i++){
    cout << mnd->rand() << endl;
  }
  
}

