#include <iostream>
#include <string.h>

using namespace std;

int loadConfig(FILE* config, int* Nparams, int* buffer_length_data, char* input_file_name, int* Nstep,  int* Burnin, char* average_file_name, char* variance_file_name, char* value_format);

int main(int argc, char** argv){
  cout << "Calculate average and variance for each parameter" << endl;
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
  // configuration for MCMC.o
  int Nparams; // number of parameters
  int buffer_length_data; // buffer size for data loading
  int buffer_length=256;
  char* input_file_name=new char[buffer_length]; // name of the input file
  int Nstep; // Number of MC steps
  int Burnin; // burn-in steps
  char* average_file_name=new char[buffer_length]; // name of the output (average) file
  char* variance_file_name=new char[buffer_length]; // name of the output (variance) file
  char* value_format=new char[buffer_length]; // value format
  
  int loadConfig_status=loadConfig(config, &Nparams, &buffer_length_data, input_file_name, &Nstep, &Burnin, average_file_name, variance_file_name, value_format);
  if(loadConfig_status!=1){
    cout << "Failed in loading configuration" << endl;
    return -1;
  }

  cout << "Number of parameters: " << Nparams << endl;
  cout << "Input file: " << input_file_name << endl;
  cout << "MC steps: " << Nstep << endl;
  cout << "Burn-in steps: " << Burnin << endl;
  cout << "Output files: " << average_file_name << ", " << variance_file_name << endl;

  buffer_length_data+=10; // for \r, \n, \0
  FILE* average=fopen(average_file_name,"w");
  if(average==NULL){
    cout << "Error in opening average file" << endl;
    return -1;
  }
  FILE* variance=fopen(variance_file_name,"w");
  if(variance==NULL){
    cout << "Error in opening variance file" << endl;
    return -1;
  }

  double* sum_array=new double[Nparams];
  double* sum2_array=new double[Nparams];
  int i;
  for(i=0;i<Nparams;i++){
    sum_array[i]=0;
    sum2_array[i]=0;
  }
      
  FILE* input=fopen(input_file_name,"r");
  if(input==NULL){
    cout << "Error in opening input file" << endl;
    return -1;
  }
  int state_count=0;
  int line_count=0;

  fprintf(average, "# average calculation\n");
  fprintf(average, "# input file: %s\n", input_file_name);

  fprintf(variance, "# variance calculation\n");
  fprintf(variance, "# input file: %s\n", input_file_name);

  char* buffer=new char[buffer_length_data];
  char* format_buff1=new char[Nparams*4+1];
  char* format_buff2=new char[Nparams*4+1];
  char* format_buff3;
  const char* format_iter="%%*lf%s";
  

  int sscanf_status;
  int sscanf_status_all;
  double* loaded_value=new double[Nparams];
  while(fgets(buffer, buffer_length_data, input)!=NULL){
    line_count++;
    if(buffer[0]=='#'){
      continue;
    }
    // data load 
    format_buff1[0]='%';
    format_buff1[1]='l';
    format_buff1[2]='f';
    format_buff1[3]='\0';
    sscanf_status_all=1;
    for(i=0;i<Nparams;i++){
      sscanf_status=sscanf(buffer, format_buff1, &loaded_value[i]);
      if(sscanf_status!=1){
	sscanf_status_all=0;
	break;
      }
      sprintf(format_buff2, format_iter, format_buff1); // buff2="%*lf"+buff1
      // replace buff1 and buff2
      format_buff3=format_buff2;
      format_buff2=format_buff1;
      format_buff1=format_buff3;
    }
    if(sscanf_status_all==1){
      state_count++;
      if(state_count>Burnin){ // state_count is already increased
	for(i=0;i<Nparams;i++){
	  sum_array[i]+=loaded_value[i];
	  sum2_array[i]+=loaded_value[i]*loaded_value[i];
	}
	int state_used=state_count-Burnin;
	// output average and variance
	for(i=0;i<Nparams;i++){
	  double average_so_far=sum_array[i]/state_used;
	  double average2_so_far=sum2_array[i]/state_used;
	  double variance_so_far=0;
	  if(state_used>1){
	    variance_so_far=(average2_so_far-average_so_far*average_so_far)*(state_used)/(state_used-1);
	  }
	  fprintf(average, value_format, average_so_far);
	  fprintf(average, " ");
	  fprintf(variance, value_format, variance_so_far);
	  fprintf(variance, " ");
	}
	fprintf(average, "\n");
	fprintf(variance, "\n");
      }
      if(state_count==Nstep){
	// break;
      }
    }else{
      cout << "sscanf failed" << endl;
    }
  }

  fclose(average);
  fclose(variance);

  if(state_count==Nstep){
    cout << "Calculation succeeded" << endl;
  }else{
    cout << "Calculation failed" << endl;
  }
  
  return 1;
}

int loadConfig(FILE* config, int* Nparams, int* buffer_length_data, char* input_file_name, int* Nstep, int* Burnin, char* average_file_name, char* variance_file_name, char* value_format){
  int buffer_length=1024;
  char line[buffer_length];
  char* fgets_status;
  int sscanf_status;

  // line 1: Number of parameters
  fgets_status=fgets(line, buffer_length, config);
  if(fgets_status==NULL){
    cout << "Error in loading line 1 of the configuration file" << endl;
    return 0;
  }
  sscanf_status=sscanf(line, "%d", Nparams);
  if(sscanf_status!=1){
    cout << "Error in parsing line 1 of the configuration file" << endl;
    return 0;
  }

  // line 2: Buffer size for loading data
  fgets_status=fgets(line, buffer_length, config);
  if(fgets_status==NULL){
    cout << "Error in loading line 2 of the configuration file" << endl;
    return 0;
  }
  sscanf_status=sscanf(line, "%d", buffer_length_data);
  if(sscanf_status!=1){
    cout << "Error in parsing line 2 of the configuration file" << endl;
    return 0;
  }

  // line 3: input_file_name
  fgets_status=fgets(line, buffer_length, config);
  if(fgets_status==NULL){
    cout << "Error in loading line 3 of the configuration file" << endl;
    return 0;
  }
  sscanf_status=sscanf(line, "%s", input_file_name);
  if(sscanf_status!=1){
    cout << "Error in parsing line 3 of the configuration file" << endl;
    return 0;
  }

  // line 4: Nstep
  fgets_status=fgets(line, buffer_length, config);
  if(fgets_status==NULL){
    cout << "Error in loading line 4 of the configuration file" << endl;
    return 0;
  }
  sscanf_status=sscanf(line, "%d", Nstep);
  if(sscanf_status!=1){
    cout << "Error in parsing line 4 of the configuration file" << endl;
    return 0;
  }
  
  // line 5: Burn-in
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

  // line 6: average_file_name
  fgets_status=fgets(line, buffer_length, config);
  if(fgets_status==NULL){
    cout << "Error in loading line 6 of the configuration file" << endl;
    return 0;
  }
  sscanf_status=sscanf(line, "%s", average_file_name);
  if(sscanf_status!=1){
    cout << "Error in parsing line 6 of the configuration file" << endl;
    return 0;
  }

  // line 7: variance_file_name
  fgets_status=fgets(line, buffer_length, config);
  if(fgets_status==NULL){
    cout << "Error in loading line 7 of the configuration file" << endl;
    return 0;
  }
  sscanf_status=sscanf(line, "%s", variance_file_name);
  if(sscanf_status!=1){
    cout << "Error in parsing line 7 of the configuration file" << endl;
    return 0;
  }

  // line 8: value_format
  fgets_status=fgets(line, buffer_length, config);
  if(fgets_status==NULL){
    cout << "Error in loading line 8 of the configuration file" << endl;
    return 0;
  }
  sscanf_status=sscanf(line, "%s", value_format);
  if(sscanf_status!=1){
    cout << "Error in parsing line 8 of the configuration file" << endl;
    return 0;
  }
  return 1;
}
