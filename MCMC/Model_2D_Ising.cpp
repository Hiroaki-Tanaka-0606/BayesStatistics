#include "Model_2D_Ising.hpp"
#include <stdio.h>
#include <string.h>
#include "../randomDistribution/mersenne_twister.hpp"

// constructor
State::State(){
  // initialization of paramters[]
  int i;
  for(i=0;i<N;i++){
    State::parameters[i]=1;
  }
}

// return data header
char* StateHeader(){
  // generate header
  // sequence of "(x,y)"
  int L=State::L;
  int N=State::N;
  // get the number of digits of L
  double Ld=(double)L;
  int digit=0;
  while(Ld>1){
    digit++;
    Ld/=10;
  }
  // header length (a little longer than the actual length);
  int headerLength=N*(3+2*digit+1)+3;
  char* header_buff1=new char[headerLength];
  char* header_buff2=new char[headerLength];
  char* header_buff3;
  header_buff1[0]='\0';
  header_buff2[0]='#';
  header_buff2[1]=' ';
  header_buff2[2]='\0';
  int i,j;
  // i -> y, j -> x
  for(i=0;i<L;i++){
    for(j=0;j<L;j++){
      // write to header_buff1
      sprintf(header_buff1, "%s(%d,%d) ", header_buff2, j, i);
      // replace header_buff1 and buff2
      header_buff3=header_buff2;
      header_buff2=header_buff1;
      header_buff1=header_buff3;
    }
  }
  delete header_buff1;
  return header_buff2;
}

// output the state
char* State::print(){
  int i,j;
  char* ret_buff1=new char[3*N];
  char* ret_buff2=new char[3*N];
  char* ret_buff3;
  ret_buff1[0]='\0';
  ret_buff2[0]='\0';
  // i -> y, j -> x
  for(i=0;i<L;i++){
    for(j=0;j<L;j++){
      // write to ret_buff1
      sprintf(ret_buff1, format_iter, ret_buff2, parameters[i*L+j]);
      // replace ret_buff1 and buff2
      ret_buff3=ret_buff2;
      ret_buff2=ret_buff1;
      ret_buff1=ret_buff3;
    }
  }
  delete ret_buff1;
  return ret_buff2;
}

// generate a next state with update of parameters[index]
void State::next(int index, double sigma, State* s_next){
  int i;
  for(i=0;i<N;i++){
    if(i==index){
      s_next->parameters[i]=(-1)*parameters[i];
    }else{
      s_next->parameters[i]=parameters[i];
    }
  }
}

// load parameters from text
int State::load(char* dataRow){
  int i;
  int len=strlen(dataRow);
  char* buffer=new char[len+1];
  strcpy(buffer, dataRow);
  int format_l_len=strlen(format_l);
  char* format_buff1=new char[(format_l_len+1)*N+1]; // format_l_len"+1" for "*"
  char* format_buff2=new char[(format_l_len+1)*N+1];
  char* format_buff3; // for replace
  strcpy(format_buff1, format_l);
  format_buff2[0]='\0';

  int sscanf_status;
  int j;
  for(i=0;i<N;i++){
    sscanf_status=sscanf(buffer, format_buff1, &parameters[i]);
    if(sscanf_status!=1){
      return -i;
    }
    sprintf(format_buff2, format_l_iter, format_buff1); //buff2="%*d"+buff1
    //replace buff1 and buff2
    format_buff3=format_buff2;
    format_buff2=format_buff1;
    format_buff1=format_buff3;
  }
  return 1;
}

void State::init_rand(){
  MT* mt=new MT();
  int i;
  for(i=0;i<N;i++){
    double rand=mt->rand();
    if(rand<0.5){
      parameters[i]=-1;
    }else{
      parameters[i]=1;
    }
  }
}

// constructor
Hamiltonian::Hamiltonian(double J){
  Hamiltonian::J=J;
}

// calculate energy
// H=-J*\sum_{<i,j>}S_iS_j <i,j> = nearest-heighbor
double Hamiltonian::energy(State* s){
  int L=State::L;
  int i,j;
  double e=0;
  // i -> y, j -> x
  for(i=0;i<L;i++){
    for(j=0;j<L;j++){
      int index=i*L+j;
      // left lattice point
      int indexL=i*L+(j-1+L)%L;
      // below lattice point
      int indexB=((i-1+L)%L)*L+j;
      // do not use right and above to avoid double-count

      int Si=s->parameters[index];
      int SjL=s->parameters[indexL];
      int SjB=s->parameters[indexB];
      // add -J*\sum_{j=left or bottom} S_iS_j
      e+=(-J)*Si*(SjL+SjB);
    }
  }
  return e;
}

// constructor
Observables::Observables(State* s){
  Observables::s=s;
}

// total magnetization
int Observables::Mag(){
  int N=State::N;
  int i;
  int magtot=0;
  for(i=0;i<N;i++){
    magtot+=s->parameters[i];
  }
  return magtot;
}

// (mag)^2
int Observables::Mag2(){
  int magtot=Observables::Mag();
  return magtot*magtot;
}
