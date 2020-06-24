all: RX_multiple.o Free_energy_1.o Free_energy_2.o RX_multiple_mND.o Free_energy_1_mND.o Free_energy_2_mND.o RX_multiple_posterior_mND.o RX_multiple_WAIC_mND.o State_histogram.o

model = ../MCMC/Model_2D_Ising.o
model_mND = ../MCMC/Model_mND.o
mt = ../randomDistribution/mersenne_twister.o
ND = ../randomDistribution/normalDistribution.a

RX_multiple.o: RX_multiple_main.o
	$(CXX) $^ $(model) $(mt) -o $@

RX_multiple_mND.o: RX_multiple_mND_main.o loadConfig_mND.o
	$(CXX) $^ $(model_mND) $(ND) -o $@

Free_energy_1.o: Free_energy_1_main.o
	$(CXX) $^ $(model) $(mt) -o $@

Free_energy_1_mND.o: Free_energy_1_mND_main.o loadConfig_mND.o
	$(CXX) $^ $(model_mND) $(ND) -o $@

Free_energy_2.o: Free_energy_2_main.o
	$(CXX) $^ $(model) $(mt) -o $@

Free_energy_2_mND.o: Free_energy_2_mND_main.o loadConfig_mND.o
	$(CXX) $^ $(model_mND) $(ND) -o $@

RX_multiple_posterior_mND.o: RX_multiple_posterior_mND_main.o
	$(CXX) $^ $(model_mND) $(ND) -o $@

RX_multiple_WAIC_mND.o: RX_multiple_WAIC_mND_main.o
	$(CXX) $^ $(model_mND) $(ND) -o $@

State_histogram.o: State_histogram_main.o
	$(CXX) $^ $(model_mND) $(ND) -o $@

.cpp.o:
	$(CXX) -c $<

.hpp.p:
	$(CXX) -c $<
