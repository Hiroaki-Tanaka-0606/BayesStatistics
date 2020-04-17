#define _USE_MATH_DEFINES
#include <iostream>
#include <string.h>
#include <cmath>
using namespace std;

//calculate posterior distribution function (model: one Gaussian) from histogram data

int main(int argc, char** argv){
	if(argc<2){
		printf("Usage: %s configuration_file\n",argv[0]);
		return 0;
	}
	cout << "# Bayes Estimation calculation" << endl;
	cout << "# Model: one Gaussian" << endl;

	//get file names
	int fileName_length=256;
	char config_file[fileName_length];
	strcpy(config_file,argv[1]);
	cout << "# Configuration file: " << config_file << endl;

	//load configuration
	double sigma_start, sigma_end;
	int sigma_split;
	double mu_start, mu_end;
	int mu_split;
	char data_file[fileName_length];
	int data_start_row;
	int num_data;
	double x_start, x_end;
	int x_split;
	FILE* config=fopen(config_file,"r");
	int buffer_length=256;
	char line[buffer_length];
	char* fgets_status;
	int sscanf_status;

	if(config==NULL){
		cout << "Error in opening the configuration file" << endl;
		return 0;
	}

	//first line: sigma range
	fgets_status=fgets(line, buffer_length, config);
	if(fgets_status==NULL){
		cout << "Error in loading 1st line of the configuration file" << endl;
		return 0;
	}
	sscanf_status=sscanf(line, "%lf%lf%d", &sigma_start, &sigma_end, &sigma_split);
	if(sscanf_status!=3){
		cout << "Error in parsing 1st line of the configuration file" << endl;
		return 0;
	}

	//second line: mu range
	fgets_status=fgets(line, buffer_length, config);
	if(fgets_status==NULL){
		cout << "Error in loading 2nd line of the configuration file" << endl;
		return 0;
	}
	sscanf_status=sscanf(line, "%lf%lf%d", &mu_start, &mu_end, &mu_split);
	if(sscanf_status!=3){
		cout << "Error in parsing 2nd line of the configuration file" << endl;
		return 0;
	}

	//third line: data file name
	fgets_status=fgets(line, buffer_length, config);
	if(fgets_status==NULL){
		cout << "Error in loading 3rd line of the configuration file" << endl;
		return 0;
	}
	sscanf_status=sscanf(line, "%s", &data_file[0]);
	if(sscanf_status!=1){
		cout << "Error in parsing 3rd line of the configuration file" << endl;
		return 0;
	}

	//fourth line: start row of the data file
	fgets_status=fgets(line, buffer_length, config);
	if(fgets_status==NULL){
		cout << "Error in loading 4th line of the configuration file" << endl;
		return 0;
	}
	sscanf_status=sscanf(line, "%d", &data_start_row);
	if(sscanf_status!=1){
		cout << "Error in parsing 4th line of the configuration file" << endl;
		return 0;
	}

	//fifth line: number of the data
	fgets_status=fgets(line, buffer_length, config);
	if(fgets_status==NULL){
		cout << "Error in loading 5th line of the configuration file" << endl;
		return 0;
	}
	sscanf_status=sscanf(line, "%d", &num_data);
	if(sscanf_status!=1){
		cout << "Error in parsing 5th line of the configuration file" << endl;
		return 0;
	}

	//sixth line: mu range
	fgets_status=fgets(line, buffer_length, config);
	if(fgets_status==NULL){
		cout << "Error in loading 6th line of the configuration file" << endl;
		return 0;
	}
	sscanf_status=sscanf(line, "%lf%lf%d", &x_start, &x_end, &x_split);
	if(sscanf_status!=3){
		cout << "Error in parsing 6th line of the configuration file" << endl;
		return 0;
	}

	//validation of the configuration
	if(sigma_start<1e-5 || sigma_end<1e-5){
		cout << "Configuration error: sigma must be positive" << endl;
		return 0;
	}
	if(sigma_split<0){
		cout << "Configuration error: sigma_split must be 0 or positive" << endl;
		return 0;
	}
	if(mu_split<0){
		cout << "Configuration error: mu_split must be 0 or positive" << endl;
		return 0;
	}
	if(data_start_row<0){
		cout << "Configuration error: start row of the data file must be 0 or positive" << endl;
		return 0;
	}
	if(num_data<0){
		cout << "Configuration error: number of data must be 0 or positive" << endl;
	}
	if(x_split<=0){
		cout << "Configuration error: x_split must be positive" << endl;
		return 0;
	}

	//output configuration
	
	if(sigma_split==0){
		cout << "# Sigma: " << sigma_start << endl;
	}else{
		double sigma_step=(sigma_end-sigma_start)/sigma_split;
		cout << "# Sigma: " << sigma_start << " to " << sigma_end << ", step " << sigma_step << endl;
	}
	if(mu_split==0){
		cout << "# Mu: " << mu_start << endl;
	}else{
		double mu_step=(mu_end-mu_start)/mu_split;
		cout << "# Mu: " << mu_start << " to " << mu_end << ", step " << mu_step << endl;
	}
	cout << "# Data file: " << data_file << endl;
	cout << "# Data row: " << data_start_row << " to " << (data_start_row+num_data-1) << endl;
	double x_step=(x_end-x_start)/x_split;
	cout << "# x: " << x_start << " to " << x_end << ", step " << x_step << endl;

	//generate sigma_list, mu_list, x_list
	double* sigma_list=new double[sigma_split+1];
	double* mu_list=new double[mu_split+1];
	double* x_list=new double[x_split+1];
	int i;
	sigma_list[0]=sigma_start;
	for(i=1; i<sigma_split+1; i++){
		sigma_list[i]=sigma_start+(sigma_end-sigma_start)/sigma_split*i;
	}
	for(i=0; i<sigma_split+1; i++){
		cout << "# sigma_list[" << i << "] = " << sigma_list[i] << endl;
	}
	mu_list[0]=mu_start;
	for(i=1; i<mu_split+1; i++){
		mu_list[i]=mu_start+(mu_end-mu_start)/mu_split*i;
	}
	for(i=0; i<mu_split+1; i++){
		cout << "# mu_list[" << i << "] = " << mu_list[i] << endl;
	}
	for(i=0; i<x_split+1; i++){
		x_list[i]=x_start+(x_end-x_start)/x_split*i;
	}
	fclose(config);
	//load data
	double* value_list=new double[num_data];
	int* count_list=new int[num_data];
	FILE* data=fopen(data_file,"r");
	if(data==NULL){
		cout << "Error in opening the data file" << endl;
		return 0;
	}
	//skip rows
	for(i=1;i<data_start_row;i++){
		fgets_status=fgets(line, buffer_length, data);
		if(fgets_status==NULL){
			cout << "Error in skipping line(s) of the data file" << endl;
			return 0;
		}
	}

	//load data
	for(i=0; i<num_data; i++){
		fgets_status=fgets(line, buffer_length, data);
		if(fgets_status==NULL){
			cout << "Error in loading data[" << i << "]" << endl;
			return 0;
		}
		sscanf_status=sscanf(line, "%lf%d", &value_list[i], &count_list[i]);
		if(sscanf_status!=2){
			cout << "Error in parsing data[" << i << "]" << endl;
			return 0;
		}
	}
	fclose(data);

	//calculate posterior distribution
	//model probability function p(x; sigma, mu) = \frac{\sigma}{\pi*((x-\mu)^2+\sigma^2)}
	double prior_distribution=1.0/((sigma_split+1)*(mu_split+1));
	double log_prior_distribution=-log(sigma_split+1)-log(mu_split+1);
	//posterioir distribution (before normalization) is stored after taking logarithm
	double* log_posterior_distribution[sigma_split+1];
	double log_distribution_function;
	int j, k;
	double sigma, mu;
	for(i=0; i<sigma_split+1; i++){
		sigma=sigma_list[i];
		double log_sigma=log(sigma);
		log_posterior_distribution[i]= new double[mu_split+1];
		for(j=0; j<mu_split+1; j++){
			mu=mu_list[j];
			log_posterior_distribution[i][j]=log_prior_distribution;
			for(k=0; k<num_data; k++){
				//add log(p(value_list[k]; sigma, mu))^count_list[k]=count_list[k]*(log_sigma-log(pi*((value_list[k]-mu)^2+sigma^2)))
				log_posterior_distribution[i][j]+=count_list[k]*(log_sigma-log(M_PI*(pow(value_list[k]-mu,2)+pow(sigma,2))));
			}
			if(i==0 && j==0){
				//first term
				log_distribution_function=log_posterior_distribution[i][j];
			}else{
				//use log-sum-exp method
				double a=log_distribution_function;
				double b=log_posterior_distribution[i][j];
				log_distribution_function=max(a,b)+log(1+exp(-abs(a-b)));
			}
		}
	}
	cout << "# Free energy: " << -log_distribution_function << endl;
	//normalization by distribution function
	double* posterior_distribution[sigma_split+1];
	for(i=0; i<sigma_split+1; i++){
		posterior_distribution[i]=new double[mu_split+1];
		for(j=0; j<mu_split+1; j++){
			posterior_distribution[i][j]=exp(log_posterior_distribution[i][j]-log_distribution_function);
		}
	}

	//calculation of estimated distribution
	cout << "# x p_x" << endl;
	double x;
	double p_x;
	for(i=0; i<x_split+1; i++){
		x=x_list[i];
		p_x=0.0;
		for(j=0; j<sigma_split+1; j++){
			sigma=sigma_list[j];
			double twopisigma=1.0/(sqrt(2*M_PI)*sigma);
			for(k=0; k<mu_split+1; k++){
				mu=mu_list[k];
				p_x+=posterior_distribution[j][k]*twopisigma*exp(-pow(x-mu,2)/(2*pow(sigma,2)));
			}
		}
		cout << x << "\t" << p_x << endl;
	}
	return 1;
}
