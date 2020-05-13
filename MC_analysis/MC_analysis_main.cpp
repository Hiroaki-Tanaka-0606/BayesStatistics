#include <iostream>

using namespace std;

int loadConfig(FILE* config, char* data_file, int* Nstep, int* Burnin, int* Nobs, char*** Obs_names_pointer, int* cor_max, char* output_format);

int main(int argc, char** argv){
  cout << "MC observable analysis" << endl;
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

  int fileName_length=256;
  char* data_file=new char[fileName_length]; // data file
  int Nstep; // number of MC steps
  int Burnin; // burn-in steps
  int Nobs; // number of observables
  char** Obs_names; // names of observables (used in the output file name)
  char*** Obs_names_pointer=new char**[1]; // to pass to loadConfig
  int cor_max; // max value of the correlation time (used in the correlation function calculation)

  char* output_format=new char[fileName_length]; // format of the output file name, including "%s" (replaced by the Obs_names[i])

  int loadConfig_status=loadConfig(config, data_file, &Nstep, &Burnin, &Nobs, Obs_names_pointer, &cor_max, output_format);

  if(loadConfig_status!=1){
    cout << "Failed in loading configuration" << endl;
    return -1;
  }
  fclose(config);

  Obs_names=Obs_names_pointer[0];

  cout << "Data file: " << data_file << endl;
  cout << "Number of MC steps: " << Nstep << endl;
  cout << "Burn-in steps: " << Burnin << endl;
  cout << "Observables: ";
  int i;
  for(i=0;i<Nobs;i++){
    cout << Obs_names[i] << ", ";
  }
  cout << endl;
  cout << "Max value of the correlation time: " << cor_max << endl << endl;

  // load data
  FILE* data_input=fopen(data_file,"r");
  int data_length=1024;
  char* buffer=new char[data_length];
  if(data_input==NULL){
    cout << "Error in opening the data file" << endl;
    return -1;
  }

  double** data=new double*[Nobs];
  for(i=0;i<Nobs;i++){
    data[i]=new double[Nstep-Burnin];
  }
  int line_count=0;
  int state_count=0;
  char* format_buff1=new char[Nobs*4+1];
  char* format_buff2=new char[Nobs*4+1];
  char* format_buff3;
  const char* format_iter="%%*lf%s";

  int loadOk=0;
  int sscanf_status;
  while(fgets(buffer, data_length, data_input)!=NULL){
    line_count++;
    if(buffer[0]=='#'){
      continue;
    }
    format_buff1[0]='%';
    format_buff1[1]='l';
    format_buff1[2]='f';
    format_buff1[3]='\0';
    int sscanf_status_all=1;
    for(i=0;i<Nobs;i++){
      double loaded_value;
      sscanf_status=sscanf(buffer, format_buff1, &loaded_value);
      if(sscanf_status!=1){
	sscanf_status_all=0;
	break;
      }
      if(state_count>=Burnin){
	data[i][state_count-Burnin]=loaded_value;
      }
      sprintf(format_buff2, format_iter, format_buff1); // buff2="%*lf"+buff1
      // replace buff1 and buff2
      format_buff3=format_buff2;
      format_buff2=format_buff1;
      format_buff1=format_buff3;
    }
    if(sscanf_status_all==1){
      state_count++;
      if(state_count==Nstep){
	loadOk=1;
	break;
      }
    }
  }
  if(loadOk!=1){
    cout << "Failed in loading the data" << endl;
    return -1;
  }else{
    cout << "Loaded data" << endl;
  }

  // data analysis
  char* output_file=new char[fileName_length];
  FILE* output;
  int j,k;
  int Ndata=Nstep-Burnin;
  double* correlation=new double[cor_max+1];
  for(i=0;i<Nobs;i++){
    sprintf(output_file, output_format, Obs_names[i]);
    output=fopen(output_file, "w");
    if(output==NULL){
      cout << "Error in opening output file" << endl;
      return -1;
    }

    // calculate average
    double data_sum=0;
    for(j=0;j<Ndata;j++){
      data_sum+=data[i][j];
    }
    double average=data_sum/Ndata;
    fprintf(output, "# Average: %.4e\n", average);

    // calculate correlation
    for(j=0;j<=cor_max;j++){
      double corr_sum=0;
      int Ndata_correlation=Ndata-j;
      for(k=0;k<Ndata_correlation;k++){
	corr_sum+=data[i][k]*data[i][k+j];
      }
      correlation[j]=corr_sum/Ndata_correlation - average*average;
    }

    // output correlation & normalized correlation
    double corr_accum=0;
    for(j=0;j<=cor_max;j++){
      double corr_norm=correlation[j]/correlation[0];
      corr_accum+=corr_norm;
      fprintf(output, "%6d %12.4e %12.4e %12.4e\n", j, correlation[j], corr_norm, corr_accum);
    }
    fclose(output);
  }
  return 1;
}

int loadConfig(FILE* config, char* data_file, int* Nstep, int* Burnin, int* Nobs, char*** Obs_names_pointer, int* cor_max, char* output_format){
  int buffer_length=1024;
  char line[buffer_length];
  char* fgets_status;
  int sscanf_status;

  // line 1: datafile
  fgets_status=fgets(line, buffer_length, config);
  if(fgets_status==NULL){
    cout << "Error in loading line 1 of the configuration file" << endl;
    return 0;
  }
  sscanf_status=sscanf(line, "%s", data_file);
  if(sscanf_status!=1){
    cout << "Error in parsing line 1 of the configuration file" << endl;
    return 0;
  }
  
  // line 2: Number of MC steps
  fgets_status=fgets(line, buffer_length, config);
  if(fgets_status==NULL){
    cout << "Error in loading line 2 of the configuration file" << endl;
    return 0;
  }
  sscanf_status=sscanf(line, "%d", Nstep);
  if(sscanf_status!=1){
    cout << "Error in parsing line 2 of the configuration file" << endl;
    return 0;
  }

  
  // line 3: Burn-in steps
  fgets_status=fgets(line, buffer_length, config);
  if(fgets_status==NULL){
    cout << "Error in loading line 3 of the configuration file" << endl;
    return 0;
  }
  sscanf_status=sscanf(line, "%d", Burnin);
  if(sscanf_status!=1){
    cout << "Error in parsing line 3 of the configuration file" << endl;
    return 0;
  }

  // line 4: Number of Observables
  fgets_status=fgets(line, buffer_length, config);
  if(fgets_status==NULL){
    cout << "Error in loading line 4 of the configuration file" << endl;
    return 0;
  }
  sscanf_status=sscanf(line, "%d", Nobs);
  if(sscanf_status!=1){
    cout << "Error in parsing line 4 of the configuration file" << endl;
    return 0;
  }

  int i;
  char** Obs_names=new char*[*Nobs];
  Obs_names_pointer[0]=Obs_names;
  // lines 5-(4+Nobs): Names of Observables
  int obs_length=64;
  for(i=0;i<*Nobs;i++){
    Obs_names[i]=new char[obs_length];

    fgets_status=fgets(line, buffer_length, config);
    if(fgets_status==NULL){
      cout << "Error in loading line " << 5+i << " of the configuration file" << endl;
      return 0;
    }
    sscanf_status=sscanf(line, "%s", Obs_names[i]);
    if(sscanf_status!=1){
      cout << "Error in parsing line " << 5+i << " of the configuration file" << endl;
      return 0;
    }    
  }

  // line 5+Nobs: max value of correlation time
  fgets_status=fgets(line, buffer_length, config);
  if(fgets_status==NULL){
    cout << "Error in loading line " << 5+*Nobs << " of the configuration file" << endl;
    return 0;
  }
  sscanf_status=sscanf(line, "%d", cor_max);
  if(sscanf_status!=1){
    cout << "Error in parsing line " << 5+*Nobs << " of the configuration file" << endl;
    return 0;
  }

  // line 6+Nobs: output file format
  fgets_status=fgets(line, buffer_length, config);
  if(fgets_status==NULL){
    cout << "Error in loading line " << 6+*Nobs << " of the configuration file" << endl;
    return 0;
  }
  sscanf_status=sscanf(line, "%s", output_format);
  if(sscanf_status!=1){
    cout << "Error in parsing line " << 6+*Nobs << " of the configuration file" << endl;
    return 0;
  }
  return 1;
}
