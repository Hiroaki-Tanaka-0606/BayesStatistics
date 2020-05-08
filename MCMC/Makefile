all: MCMC.o

MCMC.o: MCMC_main.o Model_2D_Ising.o MC_Metropolis.o
	$(CXX) $^ ../randomDistribution/mersenne_twister.o -o $@

.cpp.o:
	$(CXX) -c $<

.hpp.o:
	$(CXX) -c $<

MCMC_main.o: Model_selection.hpp MC_Metropolis.hpp
MC_Metropolis.o: Model_selection.hpp
