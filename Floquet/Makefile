evolutions = flo_evol_model.o mtrand.o
results = level_stats.o
outputs = output_func.o

main_objects = random_flo_level.o $(evolutions) $(results) $(outputs)
test_objects = test.o $(evolutions) $(results) $(outputs)

CXXFLAGS = -fopenmp
CXX = g++

all: main

main: $(main_objects)
	$(CXX) -o $@ $(CXXFLAGS) $^

test: $(test_objects)
	$(CXX) -o $@ $(CXXFLAGS) $^

clean:
	$(RM) *.o
	$(RM) .depend

depend:
	$(CXX) $(CXXFLAGS) *.cpp > .depend

-include .depend