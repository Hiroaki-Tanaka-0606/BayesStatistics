#include "Model_selection.hpp"
#include "MC_Metropolis.hpp"
#include "../randomDistribution/mersenne_twister.hpp"
#include <iostream>

using namespace std;

int main(){

  Hamiltonian* h=new Hamiltonian(1);
  State* s=new State();
  State* s_next=new State();
  State* s_buff;
  int Nstep=100000;
  int i,j;
  int num_params=s->num_params;
  double beta=0.1;
  double sigma=0;

  MT* mt=new MT();
  
  for(i=0;i<Nstep;i++){
    for(j=0;j<num_params;j++){
      int accept=MC_Metropolis(s,s_next,h,j,sigma,beta,mt);
      if(accept==1){
	//replace s and s_next
	s_buff=s_next;
	s_next=s;
	s=s_buff;
      }
      cout << accept << " ";
    }
    // calculate m
    int m=0;
    for(j=0;j<s->num_params;j++){
      m+=s->parameters[j];
    }
    cout << m << endl;
  }
    
}
