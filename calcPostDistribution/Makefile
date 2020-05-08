all: ND.o LD.o genHisto.o

ND.o: ND_test.o normalDistribution.a
	$(CXX) $^ -o $@

LD.o: LD_test.o LorentzDistribution.a
	$(CXX) $^ -o $@

genHisto.o: generateHistogram.o
	$(CXX) $^ -o $@

normalDistribution.a: normalDistribution.o mersenne_twister.o
	$(AR) rcs $@ $^

LorentzDistribution.a: LorentzDistribution.o mersenne_twister.o
	$(AR) rcs $@ $^

.cpp.o:
	$(CXX) -c $<

normalDistribution.o: mersenne_twister.hpp
LorentzDistribution.o: mersenne_twister.hpp
