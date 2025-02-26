SIMULATION_TYPENAME = "run_complete_all_years"
SCENARIO_NAME = "3_worst_case"
time_overall_start = time()

# directories
root_dir = pwd()
datapath = joinpath(root_dir,"data",SCENARIO_NAME)

# General modules
using JuMP
using JSON3
using Dates
using Logging
using LoggingExtras
# Local modules
include(joinpath(root_dir,"src/GridOptiPlan.jl"))
using .GridOptiPlan

# output directories
output_dir = joinpath(root_dir, "output")
isdir(output_dir) || mkdir(output_dir)
# create parent folder of specific simulation type
simulation_type_output_folder = joinpath(output_dir, SIMULATION_TYPENAME)
isdir(simulation_type_output_folder) || mkdir(simulation_type_output_folder)
# create folder to store results for current simulation
current_simulation_output_folder = GridOptiPlan.mk_output_dir(simulation_type_output_folder)

# Enable logging
log_file = joinpath(current_simulation_output_folder,"log.txt")
logger = FileLogger(log_file)
global_logger(logger)

@info "Started Julia using $(Base.Threads.nthreads()) threads"
@info "Scenario name: $SCENARIO_NAME"

# simulation main settings 
solver_threads = 1                          # Number of cores for each gurobi instance
iterations = parse(Int, get(ENV,"GA_ITERATIONS", "1"))  # Default: 0 (=500)
population = parse(Int, get(ENV,"GA_POPULATION", "1"))  # Default: 100
parallel_evaluation = true
time_limit = GridOptiPlan.get_time_limit()              # Time limit for the optimization process in seconds
simulation_settings = OptiSettings(joinpath(datapath,"settings.json"), current_simulation_output_folder, solver_threads, iterations, population)

# global data dictionary
@info "Loading data from given file directory"
ref = load_data_from_files(datapath, simulation_settings)
@info "Finished"

# metaheuristic algorithm search space bounds 
bounds = create_search_space(ref[:settings].upgrades, ref[:upgrades_list], ref[:upgrades_nitems])
GridOptiPlan.print_search_space(bounds, ref[:settings].upgrades, ref[:upgrades_nitems])

ga = initialize_search_algorithm(bounds, type="GA", population_size=population, iter=iterations, isparallel=parallel_evaluation, time_limit_s=time_limit)
GridOptiPlan.print_algorithm_info(ga)

# save general info
general_info = Dict(
    :scenario_name => SCENARIO_NAME,
    :tot_threads => Base.Threads.nthreads(),
    :output_folder => current_simulation_output_folder,
    :initial_iterations => ga.options.iterations, 
    :population => ga.parameters.initializer.N,
    :parallel_evaluation => ga.options.parallel_evaluation, 
    :time_limit => Dates.canonicalize(Dates.Second(round(ga.options.time_limit))), 
)
GridOptiPlan.save_dict_to_json(general_info, current_simulation_output_folder, "info_general")

# Run Bi-Level Optimization
results, best_solutions = optimize(ref, bounds, ga)

# save info
solution_info = Dict(
    :total_simulation_time => canonicalize(Dates.Second(round(time() - time_overall_start))),
    :best_objective => results.best_sol.f,
    :best_network => results.best_sol.x,
    :optimization_time => canonicalize(Dates.Second(round(results.overall_time))),
    :iteration_performed => results.iteration,
    :heuristic_termination_status => string(results.termination_status_code)
)
GridOptiPlan.save_dict_to_json(solution_info, current_simulation_output_folder, "info_solution")