#include <iostream>
#include <string.h>
#include <cmath>

using namespace std;

int loadConfig(FILE* config, int* Nparams, int* buffer_length_data, int* Nchains, char* average_file_format, char* variance_file_format, int* Nstep, char* psrf_file_format);

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
  int Nchains; // number of Marcov chains
  int buffer_length=256;
  char* average_file_format=new char[buffer_length]; // name of the input file (average)
  char* variance_file_format=new char[buffer_length]; // name of the input file (variance)
  int Nstep; // Number of MC steps (- Burnin)
  char* psrf_file_format=new char[buffer_length]; // name of the output file, including %d (parameter index)
  
  int loadConfig_status=loadConfig(config, &Nparams, &buffer_length_data, &Nchains, average_file_format, variance_file_format, &Nstep, psrf_file_format);
  if(loadConfig_status!=1){
    cout << "Failed in loading configuration" << endl;
    return -1;
  }

  cout << "Number of parameters: " << Nparams << endl;
  cout << "Number of Marcov chains: " << Nchains << endl;
  cout << "MC steps: " << Nstep << endl;

  buffer_length_data+=10; // for \r, \n, \0
  
  // open input files
  FILE** average_array=new FILE*[Nchains];
  FILE** variance_array=new FILE*[Nchains];
  // i=parameter index, j=chain index
  int i,j;
  char* file_name=new char[buffer_length];
  for(j=0;j<Nchains;j++){
    // average file
    sprintf(file_name, average_file_format, j+1);
    average_array[j]=fopen(file_name, "r");
    if(average_array[j]==NULL){
      cout << "Error in opening an average file" << endl;
      cout << file_name << endl;
      return -1;
    }
    
    // variance file
    sprintf(file_name, variance_file_format, j+1);
    variance_array[j]=fopen(file_name, "r");
    if(variance_array[j]==NULL){
      cout << "Error in opening a variance file" << endl;
      return -1;
    }
  }

  // open output files
  FILE** psrf_array=new FILE*[Nparams];
  for(i=0;i<Nparams;i++){
    sprintf(file_name, psrf_file_format, i+1);
    psrf_array[i]=fopen(file_name, "w");
    if(psrf_array[i]==NULL){
      cout << "Error in opening a psrf file" << endl;
      return -1;
    }
    fprintf(psrf_array[i], "# PSRF calculation for parameter %d\n", i+1);
    fprintf(psrf_array[i], "# Number of Marcov chains: %d\n", Nchains);
    fprintf(psrf_array[i], "# n B/n W sigma^2 V Var(V) df sqrt(R) sqrt(R df/df-2) sqrt(R df+3/df+1)\n");
  }

  int state_count=0;
  char* buffer=new char[buffer_length_data];
  char* format_buff1=new char[Nparams*4+1];
  char* format_buff2=new char[Nparams*4+1];
  char* format_buff3;
  const char* format_iter="%%*lf%s";

  // i=parameter index, j=chain index
  double* ave_sum=new double[Nparams]; // [i]= \sum_j average
  double* ave2_sum=new double[Nparams]; // [i]= \sum_j average^2
  double* var_sum=new double[Nparams]; // [i]= \sum_j variance
  double* var2_sum=new double[Nparams]; // [i]= \sum\J variance^2
  double* cov_va_sum=new double[Nparams]; // [i]= \sum_j variance*average
  double* cov_va2_sum=new double[Nparams]; // [i]= \sum_j variance*average^2


  char* fgets_status;
  int load_status=1;
  int sscanf_status;
  int sscanf_status_all;
  double* loaded_average=new double[Nparams];
  double* loaded_variance=new double[Nparams];
  while(true){
    // calculate psrf factor
    // initialize arrays
    for(i=0;i<Nparams;i++){
      ave_sum[i]=0;
      ave2_sum[i]=0;
      var_sum[i]=0;
      var2_sum[i]=0;
      cov_va_sum[i]=0;
      cov_va2_sum[i]=0;
    }
    // read one line
    for(j=0;j<Nchains;j++){
      // average
      sscanf_status_all=0;
      while(sscanf_status_all!=1){
	fgets_status=fgets(buffer, buffer_length_data, average_array[j]);
	// skip comment rows
	while(fgets_status!=NULL && buffer[0]=='#'){
	  fgets_status=fgets(buffer, buffer_length_data, average_array[j]);
	}
	if(fgets_status==NULL){
	  // reached to EOF
	  load_status=0;
	  break;
	}
	// try to get data
	format_buff1[0]='%';
	format_buff1[1]='l';
	format_buff1[2]='f';
	format_buff1[3]='\0';
	for(i=0;i<Nparams;i++){
	  sscanf_status=sscanf(buffer, format_buff1, &loaded_average[i]);
	  if(sscanf_status!=1){
	    // scan failed
	    break;
	  }
	  if(i==Nparams-1){
	    // scan succeeded
	    sscanf_status_all=1;
	    break;
	  }
	  sprintf(format_buff2, format_iter, format_buff1); // buff2="%*lf"+buff1
	  // replace buff1 and buff2
	  format_buff3=format_buff2;
	  format_buff2=format_buff1;
	  format_buff1=format_buff3;
	}
      }
      // reached to EOF before scan succeeds
      if(load_status==0){
	break;
      }
      // variance
      sscanf_status_all=0;
      while(sscanf_status_all!=1){
	fgets_status=fgets(buffer, buffer_length_data, variance_array[j]);
	// skip comment rows
	while(fgets_status!=NULL && buffer[0]=='#'){
	  fgets_status=fgets(buffer, buffer_length_data, variance_array[j]);
	}
	if(fgets_status==NULL){
	  // reached to EOF
	  load_status=0;
	  break;
	}
	// try to get data
	format_buff1[0]='%';
	format_buff1[1]='l';
	format_buff1[2]='f';
	format_buff1[3]='\0';
	for(i=0;i<Nparams;i++){
	  sscanf_status=sscanf(buffer, format_buff1, &loaded_variance[i]);
	  if(sscanf_status!=1){
	    // scan failed
	    break;
	  }
	  if(i==Nparams-1){
	    // scan succeeded
	    sscanf_status_all=1;
	    break;
	  }
	  sprintf(format_buff2, format_iter, format_buff1); // buff2="%*lf"+buff1
	  // replace buff1 and buff2
	  format_buff3=format_buff2;
	  format_buff2=format_buff1;
	  format_buff1=format_buff3;
	}
      }
      // reached to EOF before scan succeeds
      if(load_status==0){
	break;
      }

      // add value to ave_sum etc.
      for(i=0;i<Nparams;i++){
	ave_sum[i]+=loaded_average[i];
	ave2_sum[i]+=loaded_average[i]*loaded_average[i];
	var_sum[i]+=loaded_variance[i];
	var2_sum[i]+=loaded_variance[i]*loaded_variance[i];
	cov_va_sum[i]+=loaded_variance[i]*loaded_average[i];
	cov_va2_sum[i]+=loaded_variance[i]*loaded_average[i]*loaded_average[i];
      }
    }
    if(load_status==0){
      break;
    }
    state_count++;
    // calculate PSRF factor
    for(i=0;i<Nparams;i++){
      // number of MC steps
      double N=state_count*1.0;
      // number of Marcov chains
      double M=Nchains*1.0;
      
      // average (w.r.t. chain index) of average (w.r.t. MC step)
      double ave_ave=ave_sum[i]/M;
      // average of average^2
      double ave_ave2=ave2_sum[i]/M;
      // variance of average
      double var_ave=(ave_ave2-pow(ave_ave,2))*M/(M-1);
      // average of variance
      double ave_var=var_sum[i]/M;
      // variance of variance
      double var_var=(var2_sum[i]/M-pow(ave_var,2))*M/(M-1);
      // covariance of variance and average
      double cov_va=(cov_va_sum[i]/M-ave_var*ave_ave)*M/(M-1);
      // covariance of variance and average^2
      double cov_va2=(cov_va2_sum[i]/M-ave_var*ave_ave2)*M/(M-1);

      double Bn=var_ave;
      double W=ave_var;
      double sigma2=(N-1)/N*W + Bn;
      double V=sigma2+Bn/M;
      double VarV=pow((N-1)/N,2)/M*var_var
	+pow((M+1)/M,2)*2*pow(Bn,2)/(M-1)
	+2*(M+1)*(N-1)/(pow(M,2)*N)*(cov_va2-2*ave_ave*cov_va);
      double df=2*pow(V,2)/VarV;

      double PSRF1=sqrt(V/W);
      double PSRF2=sqrt(V/W*df/(df-2));
      double PSRF3=sqrt(V/W*(df+3)/(df+1));
      
      fprintf(psrf_array[i], "%8d ", state_count);
      fprintf(psrf_array[i], "%12.4e ", Bn);
      fprintf(psrf_array[i], "%12.4e ", W);
      fprintf(psrf_array[i], "%12.4e ", sigma2);
      fprintf(psrf_array[i], "%12.4e ", V);
      fprintf(psrf_array[i], "%12.4e ", VarV);
      fprintf(psrf_array[i], "%12.4e ", df);
      fprintf(psrf_array[i], "%12.4e ", PSRF1);
      fprintf(psrf_array[i], "%12.4e ", PSRF2);
      fprintf(psrf_array[i], "%12.4e \n", PSRF3);
    }
    if(state_count==Nstep){
      break;
    }
  }

  if(state_count==Nstep){
    cout << "Calculation succeeded" << endl;
  }else{
    cout << "Calculation failed" << endl;
  }


  // close files
  for(i=0;i<Nchains;i++){
    fclose(average_array[i]);
    fclose(variance_array[i]);
  }
  for(i=0;i<Nparams;i++){
    fclose(psrf_array[i]);
  }

  return 1;
}

int loadConfig(FILE* config, int* Nparams, int* buffer_length_data, int* Nchains, char* average_file_format, char* variance_file_format, int* Nstep, char* psrf_file_format){
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

  // line 3: Number of Marcov chains
  fgets_status=fgets(line, buffer_length, config);
  if(fgets_status==NULL){
    cout << "Error in loading line 3 of the configuration file" << endl;
    return 0;
  }
  sscanf_status=sscanf(line, "%d", Nchains);
  if(sscanf_status!=1){
    cout << "Error in parsing line 3 of the configuration file" << endl;
    return 0;
  }
  
  // line 4: average_file_format
  fgets_status=fgets(line, buffer_length, config);
  if(fgets_status==NULL){
    cout << "Error in loading line 4 of the configuration file" << endl;
    return 0;
  }
  sscanf_status=sscanf(line, "%s", average_file_format);
  if(sscanf_status!=1){
    cout << "Error in parsing line 4 of the configuration file" << endl;
    return 0;
  }
  
  // line 5: variance_file_format
  fgets_status=fgets(line, buffer_length, config);
  if(fgets_status==NULL){
    cout << "Error in loading line 5 of the configuration file" << endl;
    return 0;
  }
  sscanf_status=sscanf(line, "%s", variance_file_format);
  if(sscanf_status!=1){
    cout << "Error in parsing line 5 of the configuration file" << endl;
    return 0;
  }

  // line 6: Nstep (- Burnin)
  fgets_status=fgets(line, buffer_length, config);
  if(fgets_status==NULL){
    cout << "Error in loading line 6 of the configuration file" << endl;
    return 0;
  }
  sscanf_status=sscanf(line, "%d", Nstep);
  if(sscanf_status!=1){
    cout << "Error in parsing line 6 of the configuration file" << endl;
    return 0;
  }
  
  // line 7: psrf_file_format
  fgets_status=fgets(line, buffer_length, config);
  if(fgets_status==NULL){
    cout << "Error in loading line 7 of the configuration file" << endl;
    return 0;
  }
  sscanf_status=sscanf(line, "%s", psrf_file_format);
  if(sscanf_status!=1){
    cout << "Error in parsing line 7 of the configuration file" << endl;
    return 0;
  }
  
  return 1;
}
