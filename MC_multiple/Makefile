all: MC_multiple.o

model = ../MCMC/Model_2D_Ising.o
mt = ../randomDistribution/mersenne_twister.o

MC_multiple.o: MC_multiple_main.o
	$(CXX) $^ $(model) $(mt) -o $@

.cpp.o:
	$(CXX) -c $<

.hpp.p:
	$(CXX) -c $<
