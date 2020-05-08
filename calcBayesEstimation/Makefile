all: Bayes_oneGauss.o Bayes_oneLorentz.o

Bayes_oneGauss.o: BayesEstimation_oneGauss.o
	$(CXX) $^ -o $@
Bayes_oneLorentz.o: BayesEstimation_oneLorentz.o
	$(CXX) $^ -o $@

.cpp.o:
	$(CXX) -c $<

