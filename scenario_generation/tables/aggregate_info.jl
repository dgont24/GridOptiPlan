# allocate loads to nodes
# allocate EVs, PVs to loads
# create aggregated profiles
# scale aggregated profiles

import CSV
import Statistics
using DataFrames
import Random
import StatsBase

current_scenario = "1_"

# define paths
ROOT_DIR = pwd()
EV_DATA_PATH = joinpath(ROOT_DIR, "data", current_scenario, "ev_data") 
LOAD_DATA_PATH = joinpath(ROOT_DIR, "data", current_scenario, "load_data")
PV_DATA_PATH = joinpath(ROOT_DIR, "data", current_scenario, "pv_data")
OUTPUT_PATH = joinpath(ROOT_DIR, "data", current_scenario)

# read dataframes with info
loads_info = CSV.read(joinpath(LOAD_DATA_PATH, "load_info.csv"), DataFrame)
ev_info = CSV.read(joinpath(EV_DATA_PATH, "ev_info.csv"), DataFrame)
pv_info = CSV.read(joinpath(PV_DATA_PATH, "pv_info.csv"), DataFrame)

aggregate_node_info = CSV.read(joinpath(LOAD_DATA_PATH, "load_info.csv"), DataFrame)
aggregate_node_info = select(
    rename(
        aggregate_node_info,
        :id => :load_id,
        :name => :load_name,
        :profile_name => :load_profile_name,
        :daily_avg_kWh => :load_daily_avg_kWh,
        :yearly_avg_MWh => :load_yearly_avg_MWh,
        :meter_id => :load_meter_id,
        :code => :load_code
    ),
    :bus,
    :load_meter_id,
    :load_code,
    :load_name,
    :load_profile_name,
    :load_daily_avg_kWh,
    :load_yearly_avg_MWh
)


pv_size = [in(current_bus, pv_info.bus) ? pv_info[pv_info.bus .== current_bus, :pv_size][1] : missing for current_bus in aggregate_node_info.bus]
ev_name = [in(current_bus, ev_info.bus) ? ev_info[ev_info.bus .== current_bus, :name][1] : missing for current_bus in aggregate_node_info.bus]
ev_profile_name = [in(current_bus, ev_info.bus) ? ev_info[ev_info.bus .== current_bus, :profile_name][1] : missing for current_bus in aggregate_node_info.bus]
ev_max_ac_charging_kW = [in(current_bus, ev_info.bus) ? ev_info[ev_info.bus .== current_bus, :max_ac_charging_kW][1] : missing for current_bus in aggregate_node_info.bus]
ev_initial_consumption = [in(current_bus, ev_info.bus) ? ev_info[ev_info.bus .== current_bus, :initial_consumption_kW][1] : 0.0 for current_bus in aggregate_node_info.bus]


aggregate_node_info.init_kW = loads_info.initial_consumption_kW .+ ev_initial_consumption
rename!(aggregate_node_info, :init_kW => Symbol("load-ev_init_kW"))
aggregate_node_info.pv_size	= pv_size
aggregate_node_info.ev_name	= ev_name
aggregate_node_info.ev_profile_name	= ev_profile_name
aggregate_node_info.ev_max_ac_charging_kW = ev_max_ac_charging_kW

CSV.write(joinpath(OUTPUT_PATH, "aggregate_node_info.csv"), aggregate_node_info)