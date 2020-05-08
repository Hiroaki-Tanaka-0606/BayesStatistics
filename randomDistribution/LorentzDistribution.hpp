#ifndef INCLUDED_LD_HPP
#define INCLUDED_LD_HPP

#include "mersenne_twister.hpp"

class LD{
public:
  LD();
  LD(double sigma, double mu);
  LD(int seed, double sigma, double mu);

  double rand();

private:
  MT* mt;
  double sigma;
  double mu;
};

#endif
