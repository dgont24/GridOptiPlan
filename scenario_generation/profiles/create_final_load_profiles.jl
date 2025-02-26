# allocate loads to nodes
# allocate EVs, PVs to loads
# create aggregated profiles
# scale aggregated profiles

import CSV
import Statistics
using DataFrames
import Random
import StatsBase

# define paths
ROOT_DIR = pwd()
EV_DATA_PATH = joinpath(ROOT_DIR, raw"data\EV_data") 
LOAD_DATA_PATH = joinpath(ROOT_DIR, raw"data\Load_Profiles_CER") 
PV_DATA_PATH = joinpath(ROOT_DIR, raw"data\PV_data") 
OUTPUT_PATH = joinpath(ROOT_DIR, raw"data\Aggregated_load_profiles")

# read dataframes with info
loads_info = CSV.read(joinpath(LOAD_DATA_PATH,raw"load_info.csv"), DataFrame)
ev_info = CSV.read(joinpath(EV_DATA_PATH,"ev_info.csv"), DataFrame)
pv10 = CSV.read(joinpath(PV_DATA_PATH,"pv_10kWp.csv"), DataFrame)
pv24 = CSV.read(joinpath(PV_DATA_PATH,"pv_24kWp.csv"), DataFrame)

total_loads = nrow(loads_info)

# Create an array to store dataframes of loadprofiles
loadprofiles_container = Vector{DataFrame}(undef, total_loads)
for (i,row) in enumerate(eachrow(loads_info))
    file = row.profile
    loadprofiles_container[i] = CSV.read(joinpath(LOAD_DATA_PATH, file), DataFrame)
end
# Create an array to store dataframes of ev_profiles
evprofiles_container = Vector{DataFrame}(undef, nrow(ev_info))
for (i,row) in enumerate(eachrow(ev_info))
    file = row.profile_name
    evprofiles_container[i] = CSV.read(joinpath(EV_DATA_PATH, file), DataFrame)
end

# Find available nodes for Loads
lvfeeder_path = raw"data\European_LV_Test_Feeder_v2\European_LV_CSV"
initialLoads_df = CSV.read(joinpath(lvfeeder_path, "Loads.csv"), DataFrame)
load_nodes = initialLoads_df[:, :Bus]
# (nope) add manually 25 more nodes so that there are 80 nodes in total
# 241,651, 544, 530, 690, 891, 739, 763, 884, 868, 336, 484, 261, 283, 188, 27, 36, 66, 226

# Assign loads to nodes
# Shuffle the loads randomly
load_names = loads_info.id
shuffled_nodes = Random.shuffle(load_nodes)
# Initialize an empty dictionary to store the load assignments
load_assignments = Dict()
# Assign at least one load to each node
for node_index in shuffled_nodes
    load_assignments[node_index] = [pop!(load_names)]
end
# Assign at most one additional load to each node
for node_index in shuffled_nodes
    if length(load_assignments[node_index]) < 2 && !isempty(load_names)
        push!(load_assignments[node_index], pop!(load_names))
    end
end
# Display the load assignments
for (node_index, assigned_loads) in sort(load_assignments)
    println("Node $node_index: ", join(assigned_loads, ", "))
end
select!(loads_info, Not(:LoadNames, :BusName, :Load_Names))
CSV.write(joinpath(LOAD_DATA_PATH, "load_info.csv"), loads_info)

# Assign PVs to Loads
pv_percentage = 2/3
totalPVs = ceil(Int, pv_percentage*80)
countPV24 = 10
countPV10 = totalPVs-countPV24
group_under18 = loads_info[loads_info.group .== "<=18", :id]
pv10_loads = sort(StatsBase.sample(group_under18, countPV10, replace=false))
pv24_loads = loads_info[loads_info.group .== ">18", :id]
# Add list to load_info df
pv_node_list = Vector{Int64}(undef,80)
for i in 1:80
    if i in pv10_loads
        pv_node_list[i] = 10
    elseif i in pv24_loads
        pv_node_list[i] = 24
    else
        pv_node_list[i] = 0
    end
end
loads_info.pv_size = pv_node_list

# Assign EVs to Loads
ev_load_id = Vector{Int64}(undef,54)
for (i,row) in enumerate(eachrow(ev_info))
    ev_load_id[i] = loads_info[loads_info.meter_id .== row.meter_id, :id][1]
end
ev_info.load_id = ev_load_id
ev_profile_names = Vector{String}(undef,54)
for i in 1:54
    ev_profile_names[i] = "ev_profile_$i.csv"
end
ev_info.profile_name = ev_profile_names


# Compute Aggregate profiles and save them
# BE CAREFUL!! Don't add the PV. The are generation not loads. For pandapower it's wrong to add load as negative
# You should use the static generator instead.
for (i,row) in enumerate(eachrow(loads_info))
    load_profile = loadprofiles_container[i][:, :consumption_kW]
    #check if it has an EV
    ev_index = findfirst(x -> x == i, ev_info.load_id)
    if !isnothing(ev_index)
        ev_profile = evprofiles_container[ev_index][:, :charging_power_kW]
    else
        ev_profile = zeros(Int, 35040)
    end
    # create profile
    aggregate_profile = load_profile .+ ev_profile
    
    #create df
    df = select(loadprofiles_container[i], Not(:consumption_kW))
    df.consumption_kW = aggregate_profile

    filename = "aggregate_profile_$i.csv"
    CSV.write(joinpath(OUTPUT_PATH, filename), df)
end

# make info CSV
total_loads_info = select(loads_info, Not(:daily_avg_kWh, :yearly_avg_MWh, :profile))
# add file names
aggregate_profile_names = Vector{String}(undef,80)
for i in 1:80
    aggregate_profile_names[i] = "aggregate_profile_$i.csv"
end
total_loads_info.profile_name = aggregate_profile_names
# find inital power consumption
aggregate_profiles_container = Vector{DataFrame}(undef, nrow(total_loads_info))
for (i,row) in enumerate(eachrow(total_loads_info))
    file = row.profile_name
    aggregate_profiles_container[i] = CSV.read(joinpath(OUTPUT_PATH, file), DataFrame)
end
initial_consumption = Vector{Float64}(undef,nrow(total_loads_info))
for (i,row) in enumerate(eachrow(total_loads_info))
    load_profile = aggregate_profiles_container[i][:,:consumption_kW]
    initial_consumption[i] = load_profile[1] 
end
total_loads_info.initial_consumption_kW = initial_consumption

# add kV
value_kV = 0.416
total_loads_info.kV = fill(value_kV, nrow(total_loads_info))
total_loads_info

# add loadNames L1 ....
aggregate_load_names = Vector{String}(undef, nrow(total_loads_info))
for i in 1:nrow(total_loads_info)
    aggregate_load_names[i] = "L$i"
end
total_loads_info.Name = aggregate_load_names
CSV.write(joinpath(OUTPUT_PATH, "aggregate_load_info.csv"), total_loads_info)


# make info CSV for Generators
generators_info = select(loads_info, [:id, :bus, :pv_size, :group, :code])
rename!(generators_info,Dict(:id => "load_id",))
# create inital generation
generators_info.initial_generation_kW = zeros(Int, nrow(generators_info))
combine(groupby(generators_info,:pv_size),nrow)
generators_info_filtered = filter(row -> row.pv_size != 0, generators_info)
CSV.write(joinpath(ROOT_DIR,raw"data\PV_data", "generators_info.csv"), generators_info_filtered)
