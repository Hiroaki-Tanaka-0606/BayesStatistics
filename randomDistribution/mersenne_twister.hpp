#ifndef INCLUDED_MT_HPP
#define INCLUDED_MT_HPP

#include <random>
using namespace std;

class MT{
public:
  MT();
  MT(int seed);

  double rand();
  
private:
  void init(int seed);
  mt19937* mt;
  uniform_real_distribution<double>* rand_double;
};

#endif
