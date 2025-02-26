# choose which loads to keep
# and allocate PV profiles to each load

import CSV
import Statistics
using DataFrames
import StatsBase

ROOT_DIR = pwd()
DATA_PATH = joinpath(ROOT_DIR, raw"data\Load_Profiles_CER") 
CER_DATA_PATH = joinpath(raw"C:\Users\dgont\GradProj\CER_Electricity\CER Electricity Revised March 2012")

PV_threshold = 18
TOTAL_LOADS = 80
SME_LOADS = 8
RESIDENTIAL_LOADS = 72 #daily_avg >= 4kWh

# read all csv load profiles
filepaths = String[]
for i in 1:120
    push!(filepaths, joinpath(DATA_PATH,"load_profile_$i.csv"))
end 
load_data_df = CSV.read(filepaths, DataFrame)
consumer_info = CSV.read(joinpath(CER_DATA_PATH,"CER_Electricity_Documentation","SME and Residential allocations.csv"), DataFrame)
# read csv pv profiles
pv10_df = CSV.read(joinpath(ROOT_DIR, raw"data\PV_data", "pv_10kWp.csv"), DataFrame)
pv24_df = CSV.read(joinpath(ROOT_DIR, raw"data\PV_data", "pv_24kWp.csv"), DataFrame)


# prepare load profiles data
group_loads = groupby(load_data_df, :meter_id)
loads_info  = combine(group_loads, :consumption_kW => (x -> Statistics.mean(x) .* (96/4)) => :daily_avg_kWh, :consumption_kW => (x -> Statistics.mean(x) .* (96/4*365/1000)) => :yearly_avg_MWh)
loads_info.code = [consumer_info[consumer_info.ID .== loads_info[i,:meter_id], :Code][1] for i in 1:nrow(loads_info)]

# find threshold to allocate 10kWp and 20kWp pv installations
yearly_avg_pv10 = Statistics.mean(pv10_df[:, :p_w]) * (96/4) * (365/1000000) # in MWh
yearly_avg_pv24 = Statistics.mean(pv24_df[:, :p_w]) * (96/4) * (365/1000000) # in MWh
# allocate threshold somewhere in between so that there are some ZeroEnergyHomes for both categories.

# Apply the custom grouping function to create a new column for grouping
loads_info[!, :group] = custom_group.(loads_info.yearly_avg_MWh, PV_threshold)
# Group by the new categorical column
grouped_df = groupby(loads_info, [:code, :group])
combine(grouped_df, nrow)

# Choose which loads to keep
group1_big = findall((loads_info.code .== 1) .&& (loads_info.group .== ">18"))
group2 = findall(loads_info.code .== 2)
group1_small = findall((loads_info.code .== 1) .&& (loads_info.group .== "<=18"))
# define the remaining number of loads needed from group1
n_group1small = TOTAL_LOADS - length(group2) - length(group1_big)
# choose n loads from the group1_small
group1_small_random = StatsBase.sample(group1_small, n_group1small, replace=false)
# combine the three vector above the a single one
loads_indexes = copy(group1_big)
append!(loads_indexes, group1_small_random)
append!(loads_indexes, group2)
sort!(loads_indexes)
# create final load dataframe
loadsFiltered_info = loads_info[loads_indexes, :]
# add to the dataframe the load profile filename
loadsFiltered_info.profile = ["loads_profile_$i.csv" for i in loads_indexes]
# save dataframe to CSV file
CSV.write(joinpath(DATA_PATH, "load_info.csv"), loadsFiltered_info)

# Define a custom grouping function
function custom_group(value, threshold)
    return value <= threshold ? "<=$threshold" : ">$threshold"
end




