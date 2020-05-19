#ifndef INCLUDED_MODEL_TEMPLATE
#define INCLUDED_MODEL_TEMPLATE

// template for the Model used in MCMC simulation

// class State

class State{
public:
  // constructor (no-argument function), in which all paramters are initialized
  State();
  
  // number of parameters in the state (index 0..num_params-1)
  static const int num_params;
  
  // data header, which appears at the first row of output
  static const char* header;
  
  // output the state (one-line text)
  char* print();
  
  // load parameters from one-line text, with the format same as the output of print(), return 1 when succeeded, negative (error code) otherwise
  int load(char* dataRow); 

  // compose a state with a parameter (index) updated according to the parameter sigma (do not always depend on the parameter)
  void next(int index, double sigma, State* s_next);

  // randomly initialize
  void init_rand();
};

// Class Hamiltonian
class Hamiltonian{
  // constructor (adequately modify MCMC_main.cpp according to the arguments)
  Hamiltonian();

  // calculate the energy of the state s H(s)
  // used in RX prcess
  double energy(State* s);

  // calculate beta*(H(s)-1/beta log(psi(s)))=beta*H(s)-log(psi(s))
  // if psi(s) is an uniform prior distribution, log(psi(s)) term can be neglected
  // used in MC process
  double betaH(State* s, double beta);
  
}

#endif

  
