models = flo_evol_random.o flo_evol_random_rotation.o flo_evol_xxz.o flo_evol.o
evolutions = mtrand.o flo_level.o flo_rightmost_sigma_z.o evol_data.o entropy_per_model.o
functions = generic_func.o update.o tasks_models.o model_func.o sort.o reduced_density_left_2.o
transitions = transition.o basic_full.o basic_parity.o parity_full.o
results = level_stats.o
outputs = output_func.o eigen_output.o
initials = initial_obj.o random_product.o product_random.o random_pure_amplitude.o norm_check.o

main_objects = main.o $(models) $(evolutions) $(functions) $(results) $(outputs) $(transitions) $(initials)
test_objects = test.o $(models) $(evolutions) $(functions) $(results) $(outputs) $(transitions) $(initials)

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
