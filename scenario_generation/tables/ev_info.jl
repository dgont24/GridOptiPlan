import CSV
using DataFrames

current_scenario = "1_"

# define paths
rootdir = pwd()
scenariofolder = (joinpath(rootdir, "data", current_scenario))
evpath = joinpath(rootdir, "data", current_scenario, "ev_data") 
loadpath = joinpath(rootdir, "data", current_scenario, "load_data")
pvpath = joinpath(rootdir, "data", current_scenario, "pv_data")

scenario_info = CSV.read(joinpath(scenariofolder, "scenario_info.csv"), DataFrame)

# Create pv_info DataFrame
# Filter rows where pv_size_kW is not missing and select only the required columns
ev_info = scenario_info[.!ismissing.(scenario_info.ev_max_ac_charging_kW), 
    [:bus, :ev_max_ac_charging_kW ]]
# Rename the column pv_size_kW to pv_size
rename!(ev_info, :ev_max_ac_charging_kW => :max_ac_charging_kW)
ev_info.id = 1:size(ev_info, 1)
ev_info.name = "EV" .* string.(ev_info.id)

# Load all available profiles
ev_profiles_data = CSV.read(joinpath(evpath,"ev_initial_data.csv"), DataFrame)
# Remove rows where the 'code' column is equal to 2
ev_profiles_data = filter(row -> row.code != 2, ev_profiles_data)
ev_profiles_11 = ev_profiles_data[ev_profiles_data.max_ac_charging_power_kW .== 11, :profile_name]
ev_profiles_22 = ev_profiles_data[ev_profiles_data.max_ac_charging_power_kW .== 22, :profile_name]

ev_profile_names = String[]
for row in eachrow(ev_info)
    if row.max_ac_charging_kW .== 11
        index = rand(1:length(ev_profiles_11))
        push!(ev_profile_names, ev_profiles_11[index])
        splice!(ev_profiles_11, index)
    elseif row.max_ac_charging_kW .== 22
        index = rand(1:length(ev_profiles_22))
        push!(ev_profile_names, ev_profiles_22[index])
        splice!(ev_profiles_22, index)
    else
        throw(ArgumentError("EV not found"))
    end
end

ev_initial_consumption = Float64[]
# Iterate over each file
for filename in ev_profile_names
    # Read the CSV file into a DataFrame
    df = CSV.read(joinpath(evpath, "15min", filename), DataFrame)
    first_value = df[1, :charging_power_kW]
    push!(ev_initial_consumption, first_value)
end

ev_info.profile_name = ev_profile_names
ev_info.initial_consumption_kW = ev_initial_consumption

ev_info = ev_info[:, 
        [:bus, :id, :name, :profile_name, :initial_consumption_kW, :max_ac_charging_kW]
    ]

# save
CSV.write(joinpath(evpath, "ev_info.csv"), ev_info)