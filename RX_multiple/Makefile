all: RX_multiple.o Free_energy_1.o Free_energy_2.o

model = ../MCMC/Model_2D_Ising.o
mt = ../randomDistribution/mersenne_twister.o

RX_multiple.o: RX_multiple_main.o
	$(CXX) $^ $(model) $(mt) -o $@

Free_energy_1.o: Free_energy_1_main.o
	$(CXX) $^ $(model) $(mt) -o $@


Free_energy_2.o: Free_energy_2_main.o
	$(CXX) $^ $(model) $(mt) -o $@

.cpp.o:
	$(CXX) -c $<

.hpp.p:
	$(CXX) -c $<
