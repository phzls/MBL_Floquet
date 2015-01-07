evolutions = flo_evol_model.o mtrand.o flo_level.o
functions = generic_func.o update.o
results = level_stats.o
outputs = output_func.o

main_objects = main.o $(evolutions) $(functions) $(results) $(outputs)
test_objects = test.o $(evolutions) $(functions) $(results) $(outputs)

CXXFLAGS = -O3 -fopenmp
CXX = g++

all: mbl

mbl: $(main_objects)
	$(CXX) -o $@ $(CXXFLAGS) $^

test: $(test_objects)
	$(CXX) -o $@ $(CXXFLAGS) $^

clean:
	$(RM) *.o
	$(RM) .depend

depend:
	$(CXX) $(CXXFLAGS) *.cpp > .depend

-include .depend
