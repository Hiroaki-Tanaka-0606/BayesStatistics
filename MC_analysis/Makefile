all: MC_analysis.o

MC_analysis.o: MC_analysis_main.o
	$(CXX) $^ -o $@

.cpp.o:
	$(CXX) -c $<

.hpp.o:
	$(CXX) -c $<
