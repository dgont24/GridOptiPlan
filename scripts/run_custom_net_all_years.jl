SIMULATION_TYPENAME = "custom_net_all_years"
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

# simulation main settings 
solver_threads = 1        # Number of cores for each gurobi instance
simulation_settings = OptiSettings(joinpath(datapath,"settings.json"), output_dir, solver_threads)
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
x_oltc = Int[]
x_svr = Int[]
x_line = [1,1,1,1,1,1,1,1]
x_bat = [50,50,0]
xvec = vcat(x_oltc,x_svr,x_line,x_bat)
xvec = zeros(11) # scenario 2
# scenario 3
xvec = [1, 0, 0, 0, 0, 0, 2, 2, 2, 0, 1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 52, 0]

states = GridOptiPlan.decompose_x_vector(xvec, ref[:settings].upgrades, ref[:upgrades_nitems]);
local_ref = GridOptiPlan.active_network_configuration_unique_data(states, ref);
# Initialize model
model = GridOptiPlan.init_model(1, ref[:settings].threads_milp)
GridOptiPlan.initialize_daily_optimization_problem!(model, ref, local_ref);

# ==========================
# Capital Costs Calculation
# ==========================
# OLTC, SVR
tap_changers = union(local_ref[:oltc], local_ref[:svr])
tap_changer_cost = [GridOptiPlan.tap_changer_total_cost(element, ref[:settings]) for element in tap_changers]
tap_changer_cost_total_cost = sum(tap_changer_cost)

# Lines
line_ids = sort(collect(keys(ref[:upgrades_list][:lines])))
filtered_lines = [id for (i, id) in enumerate(line_ids) if states[:lines][i] != 0]
line_cost = Float64[]
for id in filtered_lines
        length = local_ref[:lines][id, :length]
        linecode = local_ref[:lines][id, :linecode]
        object_cost = GridOptiPlan.line_total_cost(linecode, length, ref[:linecodes], ref[:settings])
        push!(line_cost, object_cost)
end
line_cost_total_cost = sum(line_cost)

# Batteries
batteries = local_ref[:batteries]
battery_cost = [GridOptiPlan.battery_total_cost(element, ref[:settings]) for element in values(batteries)]
battery_cost_total_cost = sum(battery_cost)

# Total capital cost
total_capital_cost = tap_changer_cost_total_cost + line_cost_total_cost + battery_cost_total_cost

# ==========================
# Operational Costs Calculation
# ==========================
@info "Start of optimization"
t_start = time()
operational_cost, termination_status = GridOptiPlan.solve_multiyear_operational_problem!(model, 1, ref, local_ref, print_output=true)
total_time = round(time() - t_start; digits=2)
@info "Optimziation finished"
@info "Total time to solve: $(total_time)s"

# calculate total cost
total_cost = total_capital_cost + operational_cost
println("""
    Total capital cost: $total_capital_cost
    Total operational cost: $operational_cost
    Total cost: $total_cost
""")
