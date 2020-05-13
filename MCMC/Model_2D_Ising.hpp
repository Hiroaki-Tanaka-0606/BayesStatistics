#ifndef INCLUDED_MODEL_2D_ISING
#define INCLUDED_MODEL_2D_ISING

// 2D Ising Model (periodic boundary condition) used in MCMC simulation
// lattice size
#define LATTICE_SIZE 4

// class State
class State{
public:
  // constructor (no-argument function), in which all paramters are initialized
  State();
  
  // number of parameters in the state (index 0..num_params-1)
  const int num_params=N;
    
  // output the state (one-line text)
  char* print();
  
  // load parameters from one-line text, with the format same as the output of print(), return 1 when succeeded, negative (error code) otherwise
  int load(char* dataRow); 

  // compose a state with a parameter (index) updated according to the parameter sigma (do not always depend on the parameter)
  void next(int index, double sigma, State* s_next);

  // randomly initialize
  void init_rand();
  
  //---- added variables ----//
  // lattice size
  static const int L=LATTICE_SIZE;

  // lattice volume
  static const int N=L*L;
  
  // lattice spins, spin of (x, y) = parameters[y*L+x] = -1 or 1 (x,y=0..L-1)
  int parameters[N];

  // format for the output of one spin (not used)
  // const char* format="%+2d";

  // format for the load of one spin
  const char* format_l="%d";
  
  // iterative format for print();
  const char* format_iter="%s%+2d ";

  // iterative format for load();
  const char* format_l_iter="%%*d%s";
};

// return data header, which appears at the first row of output
char* StateHeader();

// class Hamiltonian
class Hamiltonian{
public:
  // constructor (adequately modify MCMC_main.cpp according to the arguments)
  Hamiltonian(double J);

  // calculate the energy of the state s
  double energy(State* s);

  //---- added variables ----//
  double J;
};

// class Observables
class Observables{
public:
  // constructor
  Observables(State* s);

  // calculate the total magnetization
  int Mag();

  // calculate (total magnetization)^2
  int Mag2();
private:
  State* s;
};
  

#endif

  
