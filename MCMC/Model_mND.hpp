#ifndef INCLUDED_MODEL_MND
#define INCLUDED_MODEL_MND

// mixed normal distribution used in MCMC simulation
// number of normal distributions
#define NUM_NORM_DISTRIBS 2
// min, max of mu
// only used in init_rand
#define MU_MIN -5
#define MU_MAX 5

#include "../randomDistribution/normalDistribution.hpp"

// class State
class State{
public:
  // constructor (no-argument function), in which all paramters are initialized
  State();
  
  // number of parameters in the state (index 0..num_params-1)
  const int num_params=NND*2-1;
    
  // output the state (one-line text)
  char* print();
  
  // load parameters from one-line text, with the format same as the output of print(), return 1 when succeeded, negative (error code) otherwise
  int load(char* dataRow); 

  // compose a state with a parameter (index) updated according to the parameter sigma (do not always depend on the parameter)
  void next(int index, double sigma, State* s_next);

  // randomly initialize
  void init_rand();
  
  //---- added variables ----//
public:
  // number of normal distributions
  static const int NND=NUM_NORM_DISTRIBS;
  
  // mu[NND]
  double mu[NND];

  // a[NND+1]: intensity separator
  // a[0]=0<a[1]< ...<a[NND-1]<a[NND]=1
  // a[0] and a[NND] are not paramters, but constants
  // intensity of ND[i] = a[i+1]-a[i]
  double a[NND+1];

  // calculate probability at x
  double probability(double x);

private:
  // format for the load of one parameter
  // mu[0], ..., mu[NND-1], a[1], ..., a[NND-1]
  const char* format_l="%lf";
  
  // iterative format for print();
  const char* format_iter="%s%16.8e ";

  // iterative format for load();
  const char* format_l_iter="%%*lf%s";
  
  ND* nd;
};

// return data header, which appears at the first row of output
char* StateHeader();

// class Hamiltonian
class Hamiltonian{
public:
  // constructor (adequately modify MCMC_main.cpp according to the arguments)
  // sample_file: data file of samples (histogram)
  // Nbin: number of bins in the file
  Hamiltonian(char* sample_file, int Nbin);

  // calculate the energy of the state s H(s)
  // used in RX prcess
  double energy(State* s);

  // calculate beta*(H(s)-1/beta log(psi(s)))=beta*H(s)-log(psi(s))
  // if psi(s) is an uniform prior distribution, log(psi(s)) term can be neglected
  // used in MC process
  double betaH(State* s, double beta);

  //---- added variables ----//
  // Nbin: Number of bins
  int Nbin;
  // Xvalue[Nbin]: values of bins
  double* Xvalue;

  // Count[Nbin]: heights of bins
  int* Count;
};

#endif

  
