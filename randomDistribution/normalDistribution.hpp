#ifndef INCLUDED_ND_HPP
#define INCLUDED_ND_HPP

// class ND: generate normal distribution(sigma,mu) by Box-Muller method

#include "mersenne_twister.hpp"

class ND{
public:
  ND();
  ND(double sigma, double mu);
  ND(int seed, int seed2, double sigma, double mu);

  double rand();

private:
  MT* mt1;
  MT* mt2;
  double sigma;
  double mu;
};

#endif
