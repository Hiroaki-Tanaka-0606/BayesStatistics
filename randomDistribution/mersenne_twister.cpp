#include <random>
#include <iostream>

#ifndef INCLUDED_MT
#define INCLUDED_MT

using namespace std;

class MT
{
public:
	MT(){
		random_device rd;
		MT::init(rd());
	}
	MT(int seed){
		MT::init(seed);
	}
	double rand(){
		return rand_double->operator()(*mt);
	}
private:
	void init(int seed){
		mt=new mt19937(seed);
		rand_double=new uniform_real_distribution<double>(0.0,1.0);
	}
	mt19937* mt;
	uniform_real_distribution<double>* rand_double;
};

/*
int main(){
	MT mt;
	cout << mt.rand() << endl;
}
*/

#endif
