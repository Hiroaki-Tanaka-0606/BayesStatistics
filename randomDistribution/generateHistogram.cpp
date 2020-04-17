#include <iostream>
#include <string.h>
#include <cmath>
#include <algorithm>
#include <vector>
using namespace std;

//generate histogram from the list of sample data

int main(int argc, char** argv){
	if(argc<5){
		printf("Usage: %s data_file start_row num_data bin_width",argv[0]);
		return 0;
	}
	//get data file name
	int fileName_length=256;
	char data_file[fileName_length];
	strcpy(data_file,argv[1]);
	//get data row information
	int start_row=atoi(argv[2]);
	int num_data=atoi(argv[3]);
	//get bin_width
	double bin_width=atof(argv[4]);

	if(start_row<0){
		cout << "Configuration error: start row of the data file must be 0 or positive" << endl;
		return 0;
	}
	if(num_data<=0){
		cout << "Configuration error: number of data must be positive" << endl;
	}
	if(bin_width<1e-5){
		cout << "Configuration error: width of the bin must be positive" <<endl;
	}

	cout << "# Data file: " << data_file << endl;
	cout << "# Data row: " << start_row << " to " << (start_row+num_data-1) << endl;
	cout << "# Bin width: " << bin_width << endl;

	//load data
	char* fgets_status;
	int sscanf_status;
	int buffer_length=256;
	char line[buffer_length];
	vector<int> data_list(num_data); //already quantized to bin_width*n (n: integer) and n is stored
	double data_i; //temporal data (not yet quantized)
	FILE* data=fopen(data_file,"r");
	int i;
	//skip rows
	for(i=1; i<start_row; i++){
		fgets_status=fgets(line, buffer_length, data);
		if(fgets_status==NULL){
			cout << "Error in skipping line(s) of the data file" << endl;
			return 0;
		}
	}
	//load data
	double min, max;
	for(i=0; i<num_data; i++){
		fgets_status=fgets(line,buffer_length,data);
		if(fgets_status==NULL){
			cout << "Error in loading data[" << i << "]" << endl;
			return 0;
		}
		sscanf_status=sscanf(line, "%lf", &data_i);
		if(sscanf_status!=1){
			cout << "Error in parsing data[" << i << "]" << endl;
			return 0;
		}
		if(i==0){
			min=data_i;
			max=data_i;
		}else{
			if(data_i<min){
				min=data_i;
			}
			if(data_i>max){
				max=data_i;
			}
		}
		data_list[i]=round(data_i/bin_width);
	}
	fclose(data);
	cout << "# Min: " << min << ", Max: " << max << endl;
	
	//sort
	sort(data_list.begin(), data_list.end());
	
	//output
	int num_bins=0;
	for(i=0; i<num_data-1; i++){
		if(data_list[i]!=data_list[i+1]){
			// data_list[x,i] composes a bin (x (< i) is unknown)
			num_bins++;
		}else{
			// duplicate
		}
	}
	num_bins++;
	cout << "# Number of bins: " << num_bins << endl;

	cout << "# value count" << endl;
	int count=1;
	for(i=0; i<num_data-1; i++){
		if(data_list[i]!=data_list[i+1]){
			cout << (data_list[i]*bin_width) << "\t" << count << endl;
			count=1;
		}else{
			count++;
		}
	}
	cout << (data_list[num_data-1]*bin_width) << "\t" << count << endl;
}