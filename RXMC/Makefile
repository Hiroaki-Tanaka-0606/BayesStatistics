all: RXMC.o

model = ../MCMC/Model_2D_Ising.o ../MCMC/MC_Metropolis.o ../randomDistribution/mersenne_twister.o

RXMC.o: RXMC_main.o loadConfig_2D_Ising.o
	$(CXX) $^ $(model) -o $@

.cpp.o:
	$(CXX) -c $<
