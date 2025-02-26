simulation_typename = "1y_alld"
println("Started Julia using $(Base.Threads.nthreads()) threads")
# directories
root_dir = pwd()
datapath = joinpath(root_dir,"data","initial")

println("Loading necessary modules")
include(joinpath(root_dir,"src/GridOptiPlan.jl"))
using .GridOptiPlan
using JuMP
println("Finished")

# output directories
output_dir = joinpath(root_dir, "output")
# create parent folder of specific simulation type
simulation_type_output_folder = joinpath(output_dir,simulation_typename)
isdir(simulation_type_output_folder) || mkdir(simulation_type_output_folder)
# create folder to store results for current simulation
current_simulation_output_folder = GridOptiPlan.mk_output_dir(simulation_type_output_folder)

# Run Optimization for One Year
year = 15
days = collect(1:365)

# simulation main settings 
#TODO: Those ones should be inputs by the user with a display
solver_threads = 1        # Number of cores for each gurobi instance
simulation_settings = OptiSettings(joinpath(datapath,"settings.json"), output_dir, solver_threads)
simulation_settings.first_flexibility_year = year
simulation_settings.final_year = year

# global data dictionary
println("Loading data from given file directory")
ref = load_data_from_files(datapath, simulation_settings)
println("Finished")

# initialize gurobi
GridOptiPlan.init_gurobi_envs!(length(days))

# create a state vector
x_oltc = [0]
x_svr = [3,3]
x_line = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]
x_bat = [50,50,50,50,50]
xvec = vcat(x_oltc,x_svr,x_line,x_bat)
xvec = [0, 2, 0, 1, 1, 0, 0, 2, 0, 1, 5, 4, 5, 1, 4, 0, 1, 1, 10, 8, 87, 145, 82]
#xvec = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

states = GridOptiPlan.decompose_x_vector(xvec, ref[:settings].upgrades, ref[:upgrades_nitems]);
local_ref = GridOptiPlan.active_network_configuration_unique_data(states, ref);

# Initialize arrays to store the results
# Dict to story operational costs of each day
cost_of_day = Dict(day => 0.0 for day in days)
# Battery soc curves
batteries_daily_soc = Dict(day => Dict{Int,Vector{Float64}}() for day in days)
# Solution status
termination_status_code = Dict(day => 0 for day in days)
# Count progress
days_solved = 0
lk = ReentrantLock()

println("Start of optimization")
t_start = time()
# Solve all days using multithreading
Threads.@threads for i in eachindex(days)
    day = days[i]
    # Initialize model
    model = GridOptiPlan.init_model(i, ref[:settings].threads_milp)
    GridOptiPlan.initialize_daily_optimization_problem!(model, ref, local_ref);
    daily_cost, day_termination_status_code, violations = GridOptiPlan.solve_daily_operational_problem!(year, day, model, ref, local_ref)

    # update results
    cost_of_day[day] = daily_cost
    batteries_daily_soc[day] = GridOptiPlan.get_batteries_day_soc(model,collect(values(local_ref[:batteries])),ref[:params]["timesteps"])
    termination_status_code[day] = GridOptiPlan.update_termination_status(day_termination_status_code, violations, 1)

    # Print and increment solved days
    lock(lk)
    try
        global days_solved
        days_solved += 1
        println("Solved $(days_solved) out of $(length(days))")
    finally
        unlock(lk)
    end
end

total_time = round(time() - t_start; digits=2)
println("Optimziation finished")
println("Total time to solve: $(total_time)s")

# calculate outputs
year_termination_status = maximum(values(termination_status_code))
total_cost = sum(values(cost_of_day))
# battery deg - Dict
total_battery_deg = GridOptiPlan.yearly_battery_degredation(batteries_daily_soc, collect(values(local_ref[:batteries])), ref)

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