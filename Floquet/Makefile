//INC = ~/

models = flo_evol_random.o flo_evol_random_rotation.o flo_evol_xxz.o flo_evol.o flo_evol_inter_random.o flo_evol_markov_inter_random.o flo_evol_markov_inter_random_both.o flo_evol_markov_inter_random_both_x.o flo_evol_xxz_random.o flo_evol_markov_xxz_random_both_x.o

evolutions = mersenne.o flo_level.o flo_rightmost_sigma_z.o evol_data.o entropy_per_model.o flo_evolution.o flo_evolution_simple_markov.o leftmost_spin_z_per_model.o leftmost_spin_z_one_run.o flo_leftmost_sigma_z.o flo_chain_end_sigma_z.o flo_evolution_simple_markov_one_model.o single_model_func.o full_leftmost_spin_z_per_model.o evol_data_total.o flo_evolution_density.o

functions = generic_func.o update.o tasks_models.o model_func.o sort.o reduced_density_left_2.o

transitions = transition.o basic_full.o basic_parity.o parity_full.o

results = level_stats.o

outputs = output_func.o eigen_output.o

initials = initial_obj.o random_product.o product_random.o random_pure_amplitude.o norm_check.o state_to_density.o random_pure.o largest_leftmost_spin_z_eigenstate.o leftmost_spin_z_eigenstate.o leftmost_spin_random_state.o

single_models = flo_chain_end_sigma_z_under_one_model.o flo_evolution_simple_markov_under_one_model.o

main_objects = main.o $(models) $(evolutions) $(functions) $(results) $(outputs) $(transitions) $(initials) $(single_models)

test_objects = test.o $(models) $(evolutions) $(functions) $(results) $(outputs) $(transitions) $(initials) $(single_models)

CXXFLAGS = -O3 -fopenmp $(INC)
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
