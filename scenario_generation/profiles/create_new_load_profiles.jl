# function to exctract yearly load profiles 
# creates n yearly load profiles from the first n consumers
# 15 min resolution

import CSV
using DataFrames

ROOT_DIR = pwd()
DATA_PATH = joinpath(raw"C:\Users\dgont\GradProj\CER_Electricity\CER Electricity Revised March 2012") 
OUTPUT_PATH = joinpath(ROOT_DIR, raw"data\Load_Profiles_CER") 

N_CONSUMERS = 120
N_DAYS = 365
RESOLUTION = 15 # in min
START_DATE = 366 # Jan 1, 2010 (Tue)
END_DATE = START_DATE + N_DAYS - 1

# consumer_info = CSV.read(joinpath(DATA_PATH,"CER_Electricity_Documentation","SME and Residential allocations.csv"), DataFrame)

# load dataframes
input_data = CSV.read(joinpath(DATA_PATH,"File1.txt"), DataFrame; header=false)
# pre-processing
# check for missing data
missing_data = names(input_data, any.(ismissing, eachcol(input_data))) # pick columns that contain missing values
@assert isempty(missing_data)
# rename columns
rename!(input_data, [:meter_id, :five_digit_code, :consumption_kW])
# sort 
sort!(input_data, [:meter_id, :five_digit_code])
# split five digit code
transform!(input_data, :five_digit_code => ByRow(x -> split_code(x)) => [:date, :time])
#change consumption from kWh to kW. Thus mulitply by 2 as there are 30min intervals
transform!(input_data, :consumption_kW => ByRow(x -> x * 2) => :consumption_kW)

# group dataframes by each consumer
g_consumers = groupby(input_data, :meter_id)
g_consumers_info = combine(g_consumers, :five_digit_code .=> [minimum maximum], nrow => :total_measurements)

index=1;
total_profiles=1;
while total_profiles < N_CONSUMERS + 1
    df_consumer = select(g_consumers[index], :)

    # sort dataframe according to time
    sort!(df_consumer, :five_digit_code)
    
    # check if the data is suitable to make a profile
    if !suitable_profile(df_consumer, START_DATE, END_DATE, N_DAYS)
        index += 1
        continue
    end

    load_profile_df = make_load_profile(df_consumer, START_DATE, END_DATE)
    CSV.write(joinpath(OUTPUT_PATH,"load_profile_$total_profiles.csv"), load_profile_df)
    index += 1
    total_profiles +=1
end

using Plots
df_test = CSV.read(joinpath(OUTPUT_PATH, "load_profile_116.csv"), DataFrame)

plot_range = 672
x = range(0, 7, length=plot_range)
y_data = df_test[673:1344, :consumption_kW]
plot(x, y_data)

function make_load_profile(df_consumer::DataFrame, start_date, end_date)
    # create empty dataframe to append the final data
    final_yearly_data = DataFrame(meter_id=Int64[], five_digit_code=Int64[], consumption_kW=Float64[], date=Int64[], time=Int64[])

    # find starting and ending rows of data
    gd = groupby(df_consumer, :date)
    start_loc = find_day_index(gd, start_date)
    end_loc = find_day_index(gd, end_date)
    # Group consumer data day by day
    g_daily_data = gd[start_loc:end_loc]
    daily_data_info = combine(g_daily_data, nrow => :total_timestamps)

    # append daily profiles one by one
    for (key, subdf) in pairs(g_daily_data)
        df = select(subdf,:)
        date = key.date
        # fix daily timestamps
        if daily_data_info[daily_data_info.date .== date,  :total_timestamps][1] > 48
            rows_to_remove = nrow(df) - 48
            for _ in 1:rows_to_remove
                pop!(df)
            end
        end
        if daily_data_info[daily_data_info.date .== date,  :total_timestamps][1] < 48
            rows_to_add = 48 - nrow(df)
            for _ in 1:rows_to_add
                push!(df, df[end, :])
            end
        end
        extended_df = change_resolution(df)
        append!(final_yearly_data, extended_df);
    end
    
    return final_yearly_data
end

function change_resolution(df::DataFrame)
    modified_df = df[:, [:meter_id,:consumption_kW,:date]]
    # duplicate every entry to convert from 30min to 15min resolution
    repeat!(modified_df, inner=2, outer=1);
    modified_df.time = 1:nrow(modified_df)

    codes_list = Int64[]
    for i in 1:nrow(modified_df)
        code = combine_integers(modified_df[i,:date][1], modified_df[i,:time][1])
        push!(codes_list, code)
    end
    modified_df.five_digit_code = codes_list
    reordered_df = modified_df[:, [:meter_id, :five_digit_code, :consumption_kW, :date, :time]]
    return reordered_df
end

function suitable_profile(consumer_data::DataFrame, start_date, end_date, n_days)
    # Group consumer data day by day
    g_daily_data = groupby(consumer_data, :date)

    # check that first day exists
    start_loc = find_day_index(g_daily_data, start_date)
    if isnothing(start_loc)
        return false
    end
    # check that last day exists
    end_loc = find_day_index(g_daily_data, end_date)
    if isnothing(end_loc)
        return false
    end

    g_daily_data_clipped = g_daily_data[start_loc:end_loc]
    daily_data_info = combine(g_daily_data_clipped, nrow => :total_timestamps)

    # Check that all the days are consecutive
    if !allunique(daily_data_info, :date) || nrow(daily_data_info) != n_days
        return false
    end
    # Check that all total_timestamps are between 46-50
    if !all(45 .<  daily_data_info[:, :total_timestamps] .< 51)
        return false
    end
    return true
end

function find_day_index(gd::GroupedDataFrame, date)
    combined_gd = combine(gd, nrow)
    col = combined_gd[:, :date]
    day_index = findfirst(==(date), col)
    return day_index
end

function split_code(code)
    code_str = string(code)
    first_part = parse(Int64, code_str[1:3])
    second_part = parse(Int64, code_str[4:end])
    return first_part, second_part
end

function combine_integers(first_part::Int64, second_part::Int64)::Int64
    # Convert to strings for concatenation
    first_part_str = string(first_part)
    second_part_str = lpad(string(second_part), 2, '0')

    # Concatenate and parse as an integer
    combined_integer = parse(Int64, first_part_str * second_part_str)
    return combined_integer
end

