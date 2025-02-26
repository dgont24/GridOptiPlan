SIMULATION_TYPENAME = "1d"
SCENARIO_NAME = "3_worst_case"
year = 12

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
x_svr = [3,3]
x_line = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]
x_bat = [50,50,50,50,50]
xvec = vcat(x_oltc,x_svr,x_line,x_bat)
# scenario 3
xvec = [1, 0, 0, 0, 0, 0, 2, 2, 2, 0, 1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 52, 0]

states = GridOptiPlan.decompose_x_vector(xvec, ref[:settings].upgrades, ref[:upgrades_nitems]);
local_ref = GridOptiPlan.active_network_configuration_unique_data(states, ref);
# Initialize model
model = GridOptiPlan.init_model(1, ref[:settings].threads_milp, verbose = true, direct=false)
GridOptiPlan.initialize_daily_optimization_problem!(model, ref, local_ref);

# initialize results DataFrames
result_vars = GridOptiPlan._results_variables(local_ref[:sets])
result_obj = GridOptiPlan._results_objective(local_ref[:sets]) 

# Run Optimization for One Day
# day = 244 # Test either 143 or 244 (they have the higher costs: 27 and 37 respectively)
for day in sort(collect(keys(ref[:days][year])))
    # Update the model constraints with the network consumption/generation values for the specific day
    # Active power injections at each node. Size of each array:(n_bus,96)
    loads_pu = GridOptiPlan.node_daily_timeseries(ref[:load_profiles], local_ref[:sets][:N_L], day)
    evs_pu = GridOptiPlan.node_daily_timeseries(ref[:ev_profiles], local_ref[:sets][:N_EV], day)
    gens_pu = GridOptiPlan.node_daily_timeseries(ref[:gen_profiles], local_ref[:sets][:N_G], day)
    GridOptiPlan.update_const_terms(loads_pu, evs_pu, gens_pu, model, local_ref[:sets])
    t_start = time()
    optimize!(model)
    total_time = round(time() - t_start; digits=2)

    # save results in DataFrames
    is_solved_and_feasible(model) && GridOptiPlan.save_daily_operational_problem_results(year, day, model, ref, local_ref, result_vars, result_obj, evs_pu, gens_pu)
end

# Save result dataframes to folder
GridOptiPlan.save_results_to_folder(result_vars, result_obj, current_simulation_output_folder)
info_dict = Dict(
    :year => year,
    :day => day,
    :solve_time => total_time,
    :termination_status => JuMP.termination_status(model),
    :state_vector => xvec,
    :costs_data => ref[:costs],
    :battery_costs_data => [bat.w for bat in values(local_ref[:batteries])]
)
GridOptiPlan.save_dict_to_json(info_dict, current_simulation_output_folder, "info")

# Plot
gdf_ev = GridOptiPlan.DataFrames.groupby(result_vars[:evs], [:year, :day])
gdf_pv = GridOptiPlan.DataFrames.groupby(result_vars[:gens], [:year, :day])
gdf_oltc = GridOptiPlan.DataFrames.groupby(result_vars[:oltc], [:year, :day])
gdf_svr = GridOptiPlan.DataFrames.groupby(result_vars[:svr], [:year, :day])
gdf_bat = GridOptiPlan.DataFrames.groupby(result_vars[:batteries], [:year, :day])
gdf_bus = GridOptiPlan.DataFrames.groupby(result_vars[:buses], [:year, :day])
gdf_dict = Dict(:ev => gdf_ev, :pv => gdf_pv, :oltc => gdf_oltc, :svr => gdf_svr, :bat => gdf_bat)

GridOptiPlan.ev_heatmap_daily(gdf_ev, local_ref[:sets], year, day)
GridOptiPlan.pv_heatmap_daily(gdf_pv, :pcurt, local_ref[:sets], year, day)
GridOptiPlan.pv_heatmap_daily(gdf_pv, :qgen, local_ref[:sets], year, day) 
GridOptiPlan.battery_plot(1, gdf_bat, local_ref[:batteries][43], local_ref[:sets], year, day)
GridOptiPlan.taps_plot(gdf_svr, gdf_oltc, local_ref[:sets], year, day)
GridOptiPlan.daily_overall_area_plot(gdf_dict, local_ref[:sets], year, day)
GridOptiPlan.daily_overall_area_plot2(gdf_dict, local_ref[:sets], year, day)
GridOptiPlan.all_days_bar_plot(gdf_dict, local_ref[:sets], params, year, days)
GridOptiPlan.daily_costs_bar_plot(result_obj, days)

# save image
using Plots
Plots.pdf(joinpath(pwd(),"report","results","3_worst_case","battery_52kwh_node43_yr$year" * "_day$day"))
Plots.svg(joinpath(current_simulation_output_folder,"images","all_daily_yr$year" * "_day$day"))


# Relax model to find constraints causing infeasibilities
# map = relax_with_penalty!(model)
# optimize!(model)
# termination_status(model)
# for (con, penalty) in map
#     violation = value(penalty)
#     if violation > 0
#         println("Constraint `$(name(con))` is violated by $violation")
#     end
# end

# Function to find conflicts if the optimization is infeasible
# compute_conflict!(model)
# if get_attribute(model, MOI.ConflictStatus()) == MOI.CONFLICT_FOUND
#     iis_model, _ = copy_conflict(model)
#     print(iis_model)
# end