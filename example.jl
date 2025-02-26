using GridOptiPlan

SCENARIO_NAME = "1_conservative"

# directories
root_dir = pwd()
datapath = joinpath(root_dir, "data", SCENARIO_NAME)
output_dir = joinpath(root_dir, "output")

# simulation main settings 
solver_threads = 1  # Number of cores that the solver can use
iterations = 50
population = 100  
parallel_evaluation = true
simulation_settings = OptiSettings(joinpath(datapath,"settings.json"), current_simulation_output_folder, solver_threads, iterations, population)

# global data dictionary
ref = load_data_from_files(datapath, simulation_settings)

# metaheuristic algorithm search space bounds 
bounds = create_search_space(ref[:settings].upgrades, ref[:upgrades_list], ref[:upgrades_nitems])
GridOptiPlan.print_search_space(bounds, ref[:settings].upgrades, ref[:upgrades_nitems])

ga = initialize_search_algorithm(bounds, type="GA", population_size=population, iter=iterations, isparallel=parallel_evaluation)
GridOptiPlan.print_algorithm_info(ga)

# Run Bi-Level Optimization
results, best_solutions = optimize(ref, bounds, ga)