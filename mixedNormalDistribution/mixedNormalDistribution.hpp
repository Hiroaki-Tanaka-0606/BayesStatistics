#ifndef INCLUDED_MND_HPP
#define INCLUDED_MND_HPP

//class mND: generate mixed normal distribution
#include "../randomDistribution/normalDistribution.hpp"
#include "../randomDistribution/mersenne_twister.hpp"

class mND{
public:
  mND(int N, double** conf);
  void print_weight_accum();
  double rand();

private:
  int N;
  double** conf;
  int valid=1;
  ND** nds;
  double* weight_accum;
  MT* mt;
};

#endif
