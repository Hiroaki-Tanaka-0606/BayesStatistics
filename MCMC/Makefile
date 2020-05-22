all: MCMC.o MCMC_obs.o RandomStates.o MCMC_mND.o MCMC_posterior_mND.o MCMC_opt_mND.o

model = Model_2D_Ising.o loadConfig_2D_Ising.o
model_hpp = Model_2D_Ising.hpp loadConfig_2D_Ising.hpp

MCMC.o: MCMC_main.o $(model) MC_Metropolis.o
	$(CXX) $^ ../randomDistribution/mersenne_twister.o -o $@

model_mND = Model_mND.o loadConfig_mND.o
model_hpp_mND = Model_mND.hpp loadConfig_mND.hpp

MCMC_mND.o: MCMC_mND_main.o $(model_mND) MC_Metropolis.o
	$(CXX) $^ ../randomDistribution/normalDistribution.a -o $@

MCMC_posterior_mND.o: MCMC_posterior_mND_main.o $(model_mND)
	$(CXX) $^ ../randomDistribution/normalDistribution.a -o $@


MCMC_opt_mND.o: MCMC_opt_mND_main.o $(model_mND) MC_Metropolis.o
	$(CXX) $^ ../randomDistribution/normalDistribution.a -o $@


MCMC_obs.o: MCMC_obs_main.o $(model)
	$(CXX) $^ ../randomDistribution/mersenne_twister.o -o $@

RandomStates.o: RandomStates_main.o $(model)
	$(CXX) $^ ../randomDistribution/mersenne_twister.o -o $@

.cpp.o:
	$(CXX) -c $<

MCMC_main.o: $(model_hpp) MC_Metropolis.hpp
MCMC_obs_main.o: $(model_hpp)
MC_Metropolis.o: $(model_hpp)
Model_2D_Ising.o: Model_2D_Ising.hpp

MCMC_main_mND.o: $(model_hpp_mND) MC_Metropolis.hpp
Model_mND.o: Model_mND.hpp
