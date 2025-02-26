SIMULATION_TYPENAME = "1y_repd"
SCENARIO_NAME = "3_worst_case"

@info "Started Julia using $(Base.Threads.nthreads()) threads"
@info "Scenario name: $SCENARIO_NAME"
# directories
root_dir = pwd()
datapath = joinpath(root_dir,"data",SCENARIO_NAME)

@info "Loading necessary modules"
include(joinpath(root_dir,"src/GridOptiPlan.jl"))
using .GridOptiPlan
using JuMP
@info "Finished"

# output directories
output_dir = joinpath(root_dir, "output")
# create parent folder of specific simulation type
simulation_type_output_folder = joinpath(output_dir,SIMULATION_TYPENAME)
isdir(simulation_type_output_folder) || mkdir(simulation_type_output_folder)
# create folder to store results for current simulation
current_simulation_output_folder = GridOptiPlan.mk_output_dir(simulation_type_output_folder)

# Year to run the optimization for
year = 11

# simulation main settings 
solver_threads = 1        # Number of cores for each gurobi instance
simulation_settings = OptiSettings(joinpath(datapath,"settings.json"), output_dir, solver_threads)
simulation_settings.first_flexibility_year = year
simulation_settings.final_year = year
simulation_settings.upgrades = [:oltc, :svr1, :lines, :batteries]

# global data dictionary
@info "Loading data from given file directory"
ref = load_data_from_files(datapath, simulation_settings)
@info "Finished"

# metaheuristic algorithm search space bounds 
# just for printing the search space
bounds = create_search_space(ref[:settings].upgrades, ref[:upgrades_list], ref[:upgrades_nitems])
GridOptiPlan.print_search_space(bounds, ref[:settings].upgrades, ref[:upgrades_nitems])

# initialize gurobi
GridOptiPlan.init_gurobi_envs!(1)

# create a state vector
x_oltc = [0]
x_svr = [1]
x_line = [0,0,0,0,0,0,0,0]
x_bat = [50,50,0]
xvec = vcat(x_oltc,x_svr,x_line,x_bat)
xvec = [1, 0, 0, 0, 0, 0, 2, 2, 2, 0, 1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 52, 0]

states = GridOptiPlan.decompose_x_vector(xvec, ref[:settings].upgrades, ref[:upgrades_nitems]);
local_ref = GridOptiPlan.active_network_configuration_unique_data(states, ref);
# Initialize model
model = GridOptiPlan.init_model(1, ref[:settings].threads_milp)
GridOptiPlan.initialize_daily_optimization_problem!(model, ref, local_ref);

# days to run
days = sort(collect(keys(ref[:days][year])))

# Initialize arrays to store the results
# Dict to story operational costs of each day
cost_of_day = Dict(day => 0.0 for day in days)
# Battery soc curves
batteries_daily_soc = Dict(day => Dict{Int,Vector{Float64}}() for day in days)
# Solution status
termination_status_code = Dict(day => 0 for day in days)
# Count progress
days_solved = 0

@info "Start of optimization"
t_start = time()
# Solve all days sequentially
for i in eachindex(days)
    day = days[i]
    daily_cost, day_termination_status_code, violations = GridOptiPlan.solve_daily_operational_problem!(year, day, model, 1, ref, local_ref)

    # update results
    cost_of_day[day] = daily_cost
    batteries_daily_soc[day] = GridOptiPlan.get_batteries_day_soc(model,collect(values(local_ref[:batteries])),ref[:params]["timesteps"])
    termination_status_code[day] = GridOptiPlan.update_termination_status(day_termination_status_code, violations, 1)

    # Print and increment solved days
    global days_solved += 1
    println("Solved $(days_solved) out of $(length(days))")
end

total_time = round(time() - t_start; digits=2)
@info "Optimziation finished"
@info "Total time to solve: $(total_time)s"

# calculate outputs
year_termination_status = maximum(values(termination_status_code))
# calculate total cost using regression
total_cost = GridOptiPlan.total_year_oper_cost(ref[:clusters][year], cost_of_day)
# battery deg - Dict
batteries_vec = collect(values(local_ref[:batteries]))
battery_deg_vec = GridOptiPlan.yearly_battery_degredation(batteries_daily_soc, ref[:clusters][year], batteries_vec, ref)
total_battery_deg = Dict(bat.id => battery_deg_vec[i]*100 for (i,bat) in enumerate(batteries_vec))

# save important optimization results
GridOptiPlan.save_dict_to_json(cost_of_day, current_simulation_output_folder, "daily_costs")
GridOptiPlan.save_dict_to_json(batteries_daily_soc, current_simulation_output_folder, "daily_bat_soc_curve")
GridOptiPlan.save_dict_to_json(termination_status_code, current_simulation_output_folder, "daily_termination_status")

# create dict with general info
info_dict = Dict(
    :year => year,
    :days => days,
    :tot_threads => Base.Threads.nthreads(),
    :total_operational_cost => total_cost,
    :total_battery_degradation => total_battery_deg,
    :solve_time => total_time,
    :termination_status => year_termination_status,
    :state_vector => xvec,
    :costs_data => ref[:costs],
    :battery_costs_data => [bat.w for bat in values(local_ref[:batteries])]
)
GridOptiPlan.save_dict_to_json(info_dict, current_simulation_output_folder, "info")