# First manually add the integration year to the ev_info and pv_info tables. Then run the following code:
# Update aggregate info file with the integration year of each new load

import CSV
using DataFrames

current_scenario = "2_bau"

# define paths
rootdir = pwd()
scenariofolder = (joinpath(rootdir, "data", current_scenario))
evpath = joinpath(rootdir, "data", current_scenario, "ev_data") 
loadpath = joinpath(rootdir, "data", current_scenario, "load_data")
pvpath = joinpath(rootdir, "data", current_scenario, "pv_data")

loads_info = CSV.read(joinpath(loadpath, "load_info.csv"), DataFrame)
ev_info = CSV.read(joinpath(evpath, "ev_info.csv"), DataFrame)
pv_info = CSV.read(joinpath(pvpath, "pv_info.csv"), DataFrame)
aggregate_info = CSV.read(joinpath(scenariofolder, "aggregate_node_info.csv"), DataFrame)

temp_df = select(pv_info, :bus, :integration_year => :pv_integration_year)
aggregate_info = leftjoin(aggregate_info, temp_df, on = :bus, order = :left)
temp_df = select(ev_info, :bus, :integration_year => :ev_integration_year)
aggregate_info = leftjoin(aggregate_info, temp_df, on = :bus, order = :left)

CSV.write(joinpath(scenariofolder, "aggregate_node_info.csv"), aggregate_info)