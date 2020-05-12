all: MCMC.o MCMC_obs.o

model = Model_2D_Ising.o loadConfig_2D_Ising.o

MCMC.o: MCMC_main.o $(model) MC_Metropolis.o
	$(CXX) $^ ../randomDistribution/mersenne_twister.o -o $@

MCMC_obs.o: MCMC_obs_main.o $(model)
	$(CXX) $^ -o $@

.cpp.o:
	$(CXX) -c $<

.hpp.o:
	$(CXX) -c $<

MCMC_main.o: Model_selection.hpp MC_Metropolis.hpp
MCMC_obs_main.o: Model_selection.hpp
MC_Metropolis.o: Model_selection.hpp
