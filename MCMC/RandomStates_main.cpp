#include <iostream>
#include "Model_selection.hpp"

using namespace std;
int main(int argc, char** argv){
  if(argc<2){
    printf("Usage: %s Nstate\n", argv[0]);
  }
  int Nstate=atoi(argv[1]);
  int i;
  State* s=new State();
  for(i=0;i<Nstate;i++){
    s->init_rand();
    cout << s->print() << endl;
  }
}
