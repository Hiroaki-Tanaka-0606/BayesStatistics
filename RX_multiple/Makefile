all: RX_multiple.o

model = ../MCMC/Model_2D_Ising.o
mt = ../randomDistribution/mersenne_twister.o

RX_multiple.o: RX_multiple_main.o
	$(CXX) $^ $(model) $(mt) -o $@

.cpp.o:
	$(CXX) -c $<

.hpp.p:
	$(CXX) -c $<
