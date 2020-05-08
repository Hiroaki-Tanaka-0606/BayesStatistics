#include "mersenne_twister.hpp"
#include <random>

using namespace std;

MT::MT(){
  random_device rd;
  MT::init(rd());
}
MT::MT(int seed){
  MT::init(seed);
}

double MT::rand(){
  return rand_double->operator()(*mt);
}

void MT::init(int seed){
  mt=new mt19937(seed);
  rand_double=new uniform_real_distribution<double>(0.0,1.0);
}
