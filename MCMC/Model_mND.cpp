#include "Model_mND.hpp"
#include <stdio.h>
#include <string.h>
#include <iostream>
#define _USE_MATH_DEFINES
#include <cmath>
#include "../randomDistribution/mersenne_twister.hpp"

// constructor
State::State(){
  // initialization of paramters[]
  int i;
  for(i=0;i<NND;i++){
    State::mu[i]=0.0;
  }
  for(i=0;i<NND+1;i++){
    State::a[i]=i*1.0/NND;
  }

  // initialization of ND(1,0)
  nd=new ND(1.0, 0.0);
}

// return data header
char* StateHeader(){
  // generate header
  // "mu[i]" (i=0, ..., NND-1)  and "a[%d]" (i=1, ..., NND-1)
  
  // get the number of digits of NND
  int NND=State::NND;
  double NNDd=(double)NND;
  int digit=0;
  while(NNDd>1){
    digit++;
    NNDd/=10;
  }
  // header length (a little longer than the actual length);
  int headerLength=NND*(4+digit+1)+(NND-1)*(3+digit+1)+3;
  char* header_buff1=new char[headerLength];
  char* header_buff2=new char[headerLength];
  char* header_buff3;
  header_buff1[0]='\0';
  header_buff2[0]='#';
  header_buff2[1]=' ';
  header_buff2[2]='\0';
  int i;
  for(i=0;i<NND;i++){
    // write to header_buff1
    sprintf(header_buff1, "%smu[%d] ", header_buff2, i);
    // replace header_buff1 and buff2
    header_buff3=header_buff2;
    header_buff2=header_buff1;
    header_buff1=header_buff3;
  }
  for(i=1;i<NND;i++){
    sprintf(header_buff1, "%sa[%d] ", header_buff2, i);
    header_buff3=header_buff2;
    header_buff2=header_buff1;
    header_buff1=header_buff3;
  }
  delete header_buff1;
  return header_buff2;
}

// output the state
char* State::print(){
  int i,j;
  // format=%16.8e + " "
  char* ret_buff1=new char[17*NND+17*(NND-1)+1];
  char* ret_buff2=new char[17*NND+17*(NND-1)+1];
  char* ret_buff3;
  ret_buff1[0]='\0';
  ret_buff2[0]='\0';
  for(i=0;i<NND;i++){
    // write to ret_buff1
    sprintf(ret_buff1, format_iter, ret_buff2, mu[i]);
    // replace ret_buff1 and buff2
    ret_buff3=ret_buff2;
    ret_buff2=ret_buff1;
    ret_buff1=ret_buff3;
  }
  for(i=1;i<NND;i++){
    // write to ret_buff1
    sprintf(ret_buff1, format_iter, ret_buff2, a[i]);
    // replace ret_buff1 and buff2
    ret_buff3=ret_buff2;
    ret_buff2=ret_buff1;
    ret_buff1=ret_buff3;
  }
  delete ret_buff1;
  return ret_buff2;
}

// generate a next state with update of parameters[index]
void State::next(int index, double sigma, State* s_next){
  // copy state
  int i;
  for(i=0;i<NND;i++){
    s_next->mu[i]=mu[i];
  }
  for(i=1;i<NND;i++){
    s_next->a[i]=a[i];
  }
  
  if(index<NND){
    // update of mu[index]
    double dmu=sigma*nd->rand();
    if(abs(dmu)>DMU_MAX){
      // will be rejected
      s_next->mu[index]=NAN;
    }else{
      s_next->mu[index]+=dmu;
    }
  }else{
    // update of a[index-NND+1];
    // if a[i-1] < a[i] < a[i+1] does not hold, Hamiltonian::energy returns inf -> MC_Metropolis returns rejection
    s_next->a[index-NND+1]+=sigma*nd->rand();
  }
}

// load parameters from text
int State::load(char* dataRow){
  int i;
  int len=strlen(dataRow);
  char* buffer=new char[len+1];
  strcpy(buffer, dataRow);
  int format_l_len=strlen(format_l);
  char* format_buff1=new char[(format_l_len+1)*num_params+1]; // format_l_len"+1" for "*"
  char* format_buff2=new char[(format_l_len+1)*num_params+1];
  char* format_buff3; // for replace
  strcpy(format_buff1, format_l);
  format_buff2[0]='\0';

  int sscanf_status;
  int j;
  for(i=0;i<NND;i++){
    sscanf_status=sscanf(buffer, format_buff1, &mu[i]);
    if(sscanf_status!=1){
      return -i;
    }
    sprintf(format_buff2, format_l_iter, format_buff1); //buff2="%*d"+buff1
    //replace buff1 and buff2
    format_buff3=format_buff2;
    format_buff2=format_buff1;
    format_buff1=format_buff3;
  }
  for(i=1;i<NND;i++){
    sscanf_status=sscanf(buffer, format_buff1, &a[i]);
    if(sscanf_status!=1){
      return -i;
    }
    sprintf(format_buff2, format_l_iter, format_buff1); //buff2="%*d"+buff1
    //replace buff1 and buff2
    format_buff3=format_buff2;
    format_buff2=format_buff1;
    format_buff1=format_buff3;
  }
  delete buffer;
  delete format_buff1;
  delete format_buff2;
  return 1;;
}

void State::init_rand(){
  MT* mt=new MT();
  // initialization of mu
  // random between MU_MIN and MU_MAX
  int i;
  double rand;
  for(i=0;i<NND;i++){
    rand=mt->rand();
    State::mu[i]=MU_MIN*1.0*rand+MU_MAX*1.0*(1.0-rand);
  }
  // initialization of a
  // sort (NND-1) random values in (0,1)
  double* rand_array=new double[NND-1];
  for(i=0;i<NND-1;i++){
    rand_array[i]=mt->rand();
  }
  a[0]=0.0;
  a[NND]=1.0;
  int j;
  int index;
  for(i=0;i<NND-1;i++){
    index=1;
    for(j=0;j<NND-1;j++){
      if(i!=j && rand_array[i]>rand_array[j]){
	index++;
      }
    }
    a[index]=rand_array[i];
  }
}

double State::probability(double x){
  double p=0;
  int i;
  double mu_i;
  double weight_i;
  double sigma=1.0;
  for(i=0;i<NND;i++){
    weight_i=a[i+1]-a[i];
    mu_i=mu[i];
    p+=weight_i/(sqrt(2*M_PI)*sigma)*exp(-pow((x-mu_i)/sigma,2)/2.0);
  }
  return p;
}

// constructor
Hamiltonian::Hamiltonian(char* sample_file, int Nbin){
  Hamiltonian::Nbin=Nbin;
  Hamiltonian::Xvalue=new double[Nbin];
  Hamiltonian::Count=new int[Nbin];

  // load sample data
  FILE* sample=fopen(sample_file, "r");
  if(sample==NULL){
    cout << "Error in opening sample file" << endl;
    std::exit(0);
  }
  int data_count=0;
  int buffer_length=256;
  int sscanf_status;
  char* buffer=new char[buffer_length];
  while(fgets(buffer, buffer_length, sample)!=NULL){
    if(buffer[0]=='#'){
      continue;
    }
    sscanf_status=sscanf(buffer, "%lf%d", &Xvalue[data_count], &Count[data_count]);
    if(sscanf_status==2){
      data_count++;
    }else{
      cout << "Error in parsing sample data" << endl;
    }
  }
  if(data_count==Nbin){
    cout << "Sample data is loaded successfully" << endl;
    int i;
    for(i=0;i<Nbin;i++){
      cout << "X = " << Xvalue[i] << ", Count = " << Count[i] << endl;
    }
  }else{
    cout << "Failed in loading sample data" << endl;
    std::exit(0);
  }
}

// calculate energy
// H = - \sum_i log p_s(x_i)
// p_s = mixed normal distribution with parameters s
// x_i: sample data
// Now that the data is provided by histrogram,
// H = - \sum_i Count[i] * log p_s(Xvalue[i])
double Hamiltonian::energy(State* s){
  // check if a[i] is correctly ordered
  int NND=State::NND;
  int i;
  for(i=0;i<NND;i++){
    if((s->a[i+1])-(s->a[i]) <= 0){
      return NAN;
    }
    if(!isfinite(s->mu[i])){
      return NAN;
    }
  }
  
  double e=0;
  double psxi; // probability with parameters s and Xvalue x_i
  int j;
  double sigma=1.0;
  double mu;
  for(i=0;i<Nbin;i++){
    psxi=s->probability(Xvalue[i]);
    e+=-Count[i]*1.0*log(psxi);
  }
  
  return e;
}

// calculate beta*(H(s)-1/beta log(psi(s)))=beta*H(s)-log(psi(s))
// log(psi(s)) term is neglected
double Hamiltonian::betaH(State* s, double beta){
  return beta*Hamiltonian::energy(s);
}
