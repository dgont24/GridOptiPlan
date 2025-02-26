using DataFrames
import CSV

ROOT_DIR = pwd()
DATA_DIR = joinpath(ROOT_DIR, raw"data\Load_Profiles_CER")

filepaths = String[]
for i in 1:120
    push!(filepaths, joinpath(DATA_DIR,"load_profile_$i.csv"))
end 

df = CSV.read(filepaths, DataFrame)
gd_consumers = groupby(df, :meter_id)
gd_consumers_info = combine(gd_consumers, nrow)
@assert all(gd_consumers_info.nrow .== 35040)
@assert nrow(gd_consumers_info) == 120
@assert allunique(gd_consumers_info[:, :meter_id])

# test each profle separately
# test days are unique
# test start and end dates are correct
# test number of days is correct
# test each day has 96 datapoints


@assert