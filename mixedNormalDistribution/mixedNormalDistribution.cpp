#include <iostream>

// class mND: generate mixed normal distribution
#include "mixedNormalDistribution.hpp"

mND::mND(int N, double** conf){
  // conf[N][3], conf[i]={sigma, mu, weight}
  mND::N=N;
  mND::conf=conf;
  // validation conf[i][2]>0
  double weightSum=0;
  int i;
  for(i=0; i<N; i++){
    if(conf[i][2]>0){
      weightSum+=conf[i][2];
    }else{
      cout << "Error: weight must be positive" << endl;
      mND::valid=0;
      return;
    }
  }
  nds=new ND*[N];
  weight_accum=new double[N+1];
  double weightSum_norm=0;
  for(i=0; i<N; i++){
    conf[i][2]/=weightSum;
    weight_accum[i]=weightSum_norm;
    weightSum_norm+=conf[i][2];
    nds[i]= new ND(conf[i][0], conf[i][1]);
  }
  weight_accum[N]=1.0;
  mt=new MT();
}

void mND::print_weight_accum(){
  int i;
  for(i=0; i<N+1; i++){
    printf("# weight_accum[%d] = %f\n", i, weight_accum[i]);
  }
}

double mND::rand(){
  if(valid!=1){
    cout << "Error: initialization failure" << endl;
    return (double)NULL;
  }
  double rand=mt->rand();
  int index=-1; // which ND should be used
  int i;
  for(i=0; i<N; i++){
    if(weight_accum[i] < rand && rand <= weight_accum[i+1]){
      index=i;
      break;
    }
  }
  if(index<0 || index>N){
    cout << "# Error: invalid random value by MT" << endl;
    return (double)NULL;
  }
  return nds[index]->rand();
}

