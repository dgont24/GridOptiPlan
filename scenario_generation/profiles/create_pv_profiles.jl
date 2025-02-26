import CSV
using DataFrames

ROUTE_DIR = pwd()
DATA_PATH = joinpath(ROUTE_DIR, raw"data\PV_data")

pv10_df = CSV.read(joinpath(DATA_PATH, "pv_10kWp.csv"), DataFrame)
pv24_df = CSV.read(joinpath(DATA_PATH, "pv_24kWp.csv"), DataFrame)

sanitize_pv_data(pv10_df)
sanitize_pv_data(pv24_df)

CSV.write(joinpath(DATA_PATH, "pv_10kWp.csv"), pv10_df)
CSV.write(joinpath(DATA_PATH, "pv_24kWp.csv"), pv24_df)

function sanitize_pv_data(df)
    transform!(df, :time => ByRow(x -> split.(x, ":")) => [:date, :time])
    
    # repeat to make 15min intervals instead of hourly
    select!(df, [:date, :p_w])
    repeat!(df, inner=4, outer=1) 
    
    #add time
    hours = [lpad(string(i), 2, '0') for i in 0:23]
    minutes = ["00", "15", "30", "45"]
    time_df = DataFrame(hours=repeat(hours, inner=4), minutes=repeat(minutes, outer=24))
    time_df.time = time_df.hours .* time_df.minutes

    df.time = repeat(time_df.time, outer=365)
    select!(df, [:date, :time, :p_w])
end



