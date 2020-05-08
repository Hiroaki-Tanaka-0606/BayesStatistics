all: mND.o

mND.o: mND_test.o mixedNormalDistribution.a
	$(CXX) $^ -o $@

mixedNormalDistribution.a: mixedNormalDistribution.o
	cp ../randomDistribution/normalDistribution.a $@
	$(AR) rcs $@ $^

.cpp.o:
	$(CXX) -c $<
