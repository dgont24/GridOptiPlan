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
pv_info = scenario_info[.!ismissing.(scenario_info.pv_size_kW), 
    [:bus, :pv_size_kW, :load_meter_id, :load_code,	:load_yearly_avg_MWh]]
# Rename the column pv_size_kW to pv_size
rename!(pv_info, :pv_size_kW => :pv_size)
# Sort the dataframe by the bus column
sorted_df = sort(pv_info, :bus)
# Convert the pv_size column to Int64, ensuring no missing values are present
pv_info.pv_size = convert(Vector{Int64}, pv_info.pv_size)
pv_info.id = 1:size(pv_info, 1)
pv_info.profile_name = ["pv_$(i)kWp.csv" for i in pv_info.pv_size]
pv_info.initial_generation_kW = [0 for i in 1:size(pv_info, 1)]
pv_info = pv_info[:, 
        [:bus, :id, :pv_size, :profile_name, :load_meter_id, :load_code, :load_yearly_avg_MWh, :initial_generation_kW]
    ]

CSV.write(joinpath(pvpath, "pv_info.csv"), pv_info)