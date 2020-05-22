all: RXMC.o RXMC_obs.o RXMC_mND.o RXMC_opt_mND.o

model = ../MCMC/Model_2D_Ising.o ../MCMC/MC_Metropolis.o ../randomDistribution/mersenne_twister.o

model_mND = ../MCMC/Model_mND.o ../MCMC/MC_Metropolis.o ../randomDistribution/normalDistribution.a

RXMC.o: RXMC_main.o loadConfig_2D_Ising.o
	$(CXX) $^ $(model) -o $@

RXMC_mND.o: RXMC_mND_main.o loadConfig_mND.o
	$(CXX) $^ $(model_mND) -o $@


RXMC_opt_mND.o: RXMC_opt_mND_main.o loadConfig_mND.o
	$(CXX) $^ $(model_mND) ../MCMC/Acceptance_optimization.o -o $@

RXMC_obs.o: RXMC_obs_main.o loadConfig_2D_Ising.o
	$(CXX) $^ $(model) -o $@

.cpp.o:
	$(CXX) -c $<
