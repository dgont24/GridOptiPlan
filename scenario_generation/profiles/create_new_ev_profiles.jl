# uses distributions from ElaadNL
# assume that up to 2/3 of the end-customers can own an EV in each building (Residential/SME). *Efkarpides ERA-NET SG+

import CSV
import Statistics
using DataFrames
import Dates
import StatsBase
import Interpolations
using Plots
import Random

ROOT_DIR = pwd()
DATA_PATH = joinpath(ROOT_DIR, raw"data\EV_data") 
loads_info = CSV.read(joinpath(ROOT_DIR, raw"data\Load_Profiles_CER\load_info.csv"), DataFrame)

N_LOADS = nrow(loads_info)
EV_PERCENT = 2/3
N_EVs = ceil(Int,N_LOADS * EV_PERCENT)
YEAR = 2010

# load files
arrival_weekdays = load_arrival_time(joinpath(DATA_PATH,"distribution-of-arrival-weekdays.csv"))
arrival_weekends = load_arrival_time(joinpath(DATA_PATH,"distribution-of-arrival-weekend.csv"))
connection_time_df = load_connection_time(joinpath(DATA_PATH,"distribution-of-connection-time-per-event.csv"))
energy_demand_df = load_energy_demand(joinpath(DATA_PATH,"distribution-of-energy-demand.csv"))
cars_data = CSV.read(joinpath(DATA_PATH,"car_models.csv"), DataFrame)

#determine weekend dates
weekends_2010 = weekends_of_year(YEAR)
#evprofile_container = [DataFrame(datetime = Dates.DateTime[], power_kW = Float64[]) for _ in 1:N_EVs]

# determine profile characteristics
# consider all SME loads
evs_workplace = findall((loads_info.code .== 2))
# consider the remaing home loads
evs_home = findall((loads_info.code .== 1) .&& (loads_info.daily_avg_kWh .>= 12))
remaining_evs = N_EVs - length(evs_workplace)
evs_home = StatsBase.sample(evs_home, remaining_evs, replace=false)
# combine to one list
evs = copy(evs_workplace)
append!(evs, evs_home)
sort!(evs)
meterID = [loads_info[i, :meter_id] for i in evs]
codeConsumer = [loads_info[i, :code] for i in evs]
evs_df = DataFrame([:meter_id => meterID, :code => codeConsumer])
# fix evs Percentage to sum to 100
evs_percentages = cars_data[:, :share_percentage]
evs_percentages ./= sum(evs_percentages)
cars_data[!, :share_percentage] = evs_percentages
# find number of cars from each brand
brands_count = Int64[]
for row in eachrow(cars_data)
    car_count = round(Int, row.share_percentage * N_EVs)
    append!(brands_count, car_count)
end
sum(brands_count)
cars_data.ncount = brands_count
# Generate a vector with repeated letters based on the count vector
car_assignments = [string(car_type) for (car_type, count) in zip(cars_data[:,:car_model], cars_data[:,:ncount]) for _ in 1:count]
# Shuffle the car assignments
shuffled_assignments = Random.randperm(length(car_assignments))
# Assign the shuffled cars to the DataFrame
evs_df.car_model = car_assignments[shuffled_assignments[1:nrow(evs_df)]]
# Add capacity column
evs_df.usable_energy_capacity_kWh = [cars_data[cars_data.car_model .== evs_df[i, :car_model], :usable_energy_capacity_kWh][1] for i in 1:nrow(evs_df)]
evs_df.usable_energy_capacity_kWh = round.(Int64, evs_df.usable_energy_capacity_kWh)
# Add charging power column
evs_df.max_ac_charging_power_kW = [cars_data[cars_data.car_model .== evs_df[i, :car_model], :max_ac_charging_power_kW][1] for i in 1:nrow(evs_df)]
evs_df
#CSV.write(joinpath(DATA_PATH, "ev_info.csv"), evs_df)

# produce profiles
index=1;
for row in eachrow(evs_df)
    if row.code == 1
        sector = "private"
    else
        sector = "workplace"
    end
    df = create_ev_load(sector, row)
    df.meter_id = fill(row.meter_id , nrow(df))
    output_name = "ev_profile_$index.csv"
    CSV.write(joinpath(ROOT_DIR, DATA_PATH, output_name), df)
    index += 1;
end

function create_ev_load(sector::String, car_info::DataFrameRow, year::Int64=YEAR, weekends::Vector{Dates.Date}=weekends_2010)
# sector = "private", "workplace"
    # create the dataframe
    dates = repeat(1:1:365, inner = 96)
    times = repeat(1:1:96, outer = 365)
    df = DataFrame(date = dates, time = times)
    # car specs
    car_energy_capacity = car_info["usable_energy_capacity_kWh"] * 4 # in kW(15min)
    car_charging_power = car_info["max_ac_charging_power_kW"]
    
    residual_connection_time = 0
    residual_energy_demand = 0
    chargingprofile = Vector{Int64}(undef, nrow(df))
    for i in 1:365
        date = Dates.Date(year, 1, 1) + Dates.Day(i - 1)
        # Step 1
        # determine arrival time, connection duration
        if (date in weekends) && sector == "workplace"
            arrival_time = 0
            connection_duration = 0
        elseif residual_connection_time == 0
            # find arrival time (in 15 min multiples)
            arrival_time = find_arrival_time(date, sector)
            # find connection duration (in 15 min multiples)
            connection_duration = find_connection_duration_h(sector) * 4
        else
            arrival_time = 0
            connection_duration = residual_connection_time
        end

        # Step 2
        # determine energy demand in kW(15min) = 4 kWh
        if (date in weekends) && sector == "workplace"
            energy_demand = 0
        elseif residual_energy_demand == 0
            energy_demand = find_energy_demand(sector) * 4 # convert from kWh to kW(15min)
            # make sure eneryg demand is not more than vehicles capacity
            if energy_demand > car_energy_capacity
                energy_demand = car_energy_capacity
            end           
        else
            energy_demand = residual_energy_demand
        end #

        #= Step 2 - Modified - Because if there is some residual_energy_demand the person can still use this car today
        # determine energy demand in kW(15min) = 4 kWh
        today_energy_use = find_energy_demand(sector) * 4 # convert from kWh to kW(15min)
        energy_demand = today_energy_use + residual_energy_demand
        # make sure energy demand is not more than vehicles capacity. Only then she cannot use it and has to charge the vehicle first
        if energy_demand > car_energy_capacity
            energy_demand = car_energy_capacity
        end =#        
        
        # Step 1.1
        # find connection duration for the specific day only
        connection_duration_today = 0
        if (96 - arrival_time - connection_duration) >= 0
            connection_duration_today = connection_duration
        else
            connection_duration_today  = 96 - arrival_time
        end

        # Step 2.1
        # calculate charging time needed to reach SoC=1 (in 15min multiples)
        charging_time_SoC1 = ceil(Int, (energy_demand / car_charging_power))
        # stop the charging if SoC reaches value 1
        if charging_time_SoC1 <= connection_duration_today
            connection_duration_today = charging_time_SoC1
            connection_duration =  charging_time_SoC1
        end 

        # Step 2.2 
        # calculate future residual energy demand
        residual_charging_time_SoC1 = charging_time_SoC1 - connection_duration_today
        residual_energy_demand = residual_charging_time_SoC1 * car_charging_power

        # Step 1.2
        #check for future residual connection time
        if (arrival_time + connection_duration) > 96
            residual_connection_time = arrival_time + connection_duration - 96
        else
            residual_connection_time = 0
        end

        # Step 3 - create daily charging profiles
        # create list profile_container
        j = 96*(i-1) + 1
        if arrival_time == 0 && connection_duration_today == 0
            chargingprofile[j:(j+95)] .= 0
        elseif arrival_time == 0 && connection_duration_today > 0 && connection_duration_today < 96
            chargingprofile[j:(j+connection_duration_today-1)] .= car_charging_power 
            chargingprofile[(j+connection_duration_today):(j+95)] .= 0 
        elseif arrival_time == 0 && connection_duration_today == 96
            chargingprofile[j:(j+95)] .= car_charging_power 

        elseif arrival_time > 0 && connection_duration_today == 0
            chargingprofile[j:(j+95)] .= 0
        elseif arrival_time > 0 && connection_duration_today > 0 && connection_duration_today < (96 - arrival_time)
            chargingprofile[j:(j+arrival_time-1)] .= 0
            chargingprofile[(j+arrival_time):(j+arrival_time+connection_duration_today-1)] .= car_charging_power
            chargingprofile[(j+arrival_time+connection_duration_today):(j+95)] .= 0 
        elseif arrival_time > 0 && connection_duration_today == (96 - arrival_time)
            chargingprofile[j:(j+arrival_time-1)] .= 0
            chargingprofile[(j+arrival_time):(j+arrival_time+connection_duration_today-1)] .= car_charging_power
        end
    end 

    df.charging_power_kW = chargingprofile
    return df
end

function find_energy_demand(sector::String, demand_distr::DataFrame=energy_demand_df) :: Int64
    values = demand_distr[:, :energy_demand_kWh]
    cdf_col = "cdf_$sector"
    cdf = demand_distr[:, cdf_col]
    energy_demand_kWh = generate_random_number(values, cdf) 
    return Int64(energy_demand_kWh)
end

function find_connection_duration_h(sector::String, connections_distr::DataFrame=connection_time_df) :: Int64
    values = connections_distr[:, :charging_time]
    cdf_col = "cdf_$sector"
    cdf = connections_distr[:, cdf_col]
    connection_duration = generate_random_number(values, cdf) 
    return Int64(connection_duration)
end

function find_arrival_time(date::Dates.Date, sector::String, weekends::Vector{Dates.Date}=weekends_2010, 
                            weekdays_distr::DataFrame=arrival_weekdays, weekends_distr::DataFrame=arrival_weekends) :: Int64
    arrival_time = 0
    cdf_col = "cdf_$sector"
    if date in weekends
        arrival_time = generate_random_number(weekends_distr[:, :time], weekends_distr[:, cdf_col])
    else
        arrival_time = generate_random_number(weekdays_distr[:, :time], weekdays_distr[:, cdf_col])
    end
    return Int(arrival_time)
end       
        

function weekends_of_year(year)
    # Get the first and last day of the year
    year=YEAR
    start_date = Dates.Date(year, 1, 1)
    end_date = Dates.Date(year, 12, 31)

    # Generate a range of dates from the first day to the last day
    all_dates = start_date:Dates.Day(1):end_date

    # Filter out only the weekends (Saturday and Sunday)
    weekends = filter(date -> Dates.dayofweek(date) in [6, 7], all_dates)

    return collect(weekends)
end


function load_arrival_time(path)
    df = CSV.read(path, DataFrame)
    # rename the private,public, workplace columns
    rename!(df, Dict(:private => :pdf_private, :public => :pdf_public))
    if "workplace" in names(df)
        rename!(df, Dict(:workplace => :pdf_workplace))
    end
    # Normalize each PDF to make sure it sums to 1
    normalize_pdf!(df)
    # Compute the cumulative distribution function (CDF) for each column
    cdf_from_pdf!(df)
    # add an extra column for counting total daily 15min intervals
    df.time = 1:1:96
    return df
end

##############
##############
#=
The data looks a bit confusing. What is given exactly in the data. If the complementary CDF is given then
the commented functions below are the correct ones. So you just have to normalize and do 1-y
But if you plot the pdf this way it is very abrapt and so many short connection times are not logically
jsutified.

function load_connection_time(path)
    df = CSV.read(path, DataFrame)
    # normalize the cdf function
    normalize_cdf!(df)
    # convert data from 1 - F(x) to F(x)
    complementary_to_cdf!(df)
    # rename the private,public, workplace columns
    rename!(df, Dict(:private => :cdf_private, :public => :cdf_public, :workplace => :cdf_workplace))
    # change the first column name and scale the values
    rename!(df, Dict(:"Percentage of charging events" => :"charging_time"))
    df.charging_time = 0:0.72:72
    return df
end

function load_energy_demand(path)
    df = CSV.read(path, DataFrame)
    # normalize the cdf function
    normalize_cdf!(df)
    # convert data from 1 - F(x) to F(x)
    complementary_to_cdf!(df)
    # rename the private,public, workplace columns
    rename!(df, Dict(:private => :cdf_private, :public => :cdf_public, :workplace => :cdf_workplace))
    # change the first column name and scale the values
    rename!(df, Dict(:"Percentage of charging events" => :"energy_demand_kWh"))
    df.charging_time = 0:1:100
    return df
end

function complementary_to_cdf!(df)
    col_names = [:private, :public, :workplace]
    for col in col_names
        y = df[:, col]
        df[!, col] = 1 .- y
    end
end
=#     

function load_connection_time(path)
    df = CSV.read(path, DataFrame)
    # Reverse the order of the columns private, public, workplace
    reverse_column_order!(df)
    # Compute PDF for each column
    compute_pdf!(df, "connection_time")
    # Compute the cumulative distribution function (CDF) for each column
    cdf_from_pdf!(df)
    # Drop initial columns
    select!(df, Not(:private, :public, :workplace))
    df = df[1:73, :]
    # Change data types from Float64? to Float64 for specified columns
    columns_to_convert = [:pdf_private, :pdf_public, :pdf_workplace, :cdf_private, :cdf_public, :cdf_workplace]
    for col in columns_to_convert
        df[!, col] = convert(Vector{Float64}, df[!, col])
    end
    # change the first column name and scale the values
    rename!(df, Dict(:"Percentage of charging events" => :"charging_time"))
    df.charging_time = 0:1:72

    return df
end


function load_energy_demand(path)
    df = CSV.read(path, DataFrame)
    # Reverse the order of the columns private, public, workplace
    reverse_column_order!(df)
    # Compute PDF for each column
    compute_pdf!(df, "energy_demand")
    # Compute the cumulative distribution function (CDF) for each column
    cdf_from_pdf!(df)
    # Drop initial columns
    select!(df, Not(:private, :public, :workplace))
    # change the first column name and scale the values
    rename!(df, Dict(:"Percentage of charging events" => :"energy_demand_kWh"))
    df.energy_demand_kWh = 0:1:100

    return df
end


function normalize_pdf!(df::DataFrame)
# Function to normalize each column

    df[!, :pdf_private] ./= sum(df.pdf_private)
    df[!, :pdf_public] ./= sum(df.pdf_public)
    if "pdf_workplace" in names(df)
        df[!, :pdf_workplace] ./= sum(df.pdf_workplace)
    end
end

function normalize_cdf!(df::DataFrame)
# Function to normalize each column

    df[!, :private] ./= maximum(df.private)
    df[!, :public] ./= maximum(df.public)
    df[!, :workplace] ./= maximum(df.workplace)
end

function cdf_from_pdf!(df::DataFrame)
# Function to compute the cumulative distribution functions (CDF)
    
    df[!, :cdf_private] = cumsum(df.pdf_private)
    df[!, :cdf_public] = cumsum(df.pdf_public)
    if "pdf_workplace" in names(df)
        df[!, :cdf_workplace] = cumsum(df.pdf_workplace)
    end
end

function compute_pdf!(df::DataFrame, mode::String = "connection_time")
# Compute PDF from CDF using discrete differences
# mode = "connection_time" or "energy_demand"

    col_names = [:private, :public, :workplace]
    new_cols_names = [:pdf_private, :pdf_public, :pdf_workplace]
    for (i, col) in enumerate(col_names)
        # Interpolate
        # Define the new number of points
        new_num_points = 100001
        x = range(0, 1, nrow(df))
        y = df[!, col]
        itp_cubic = Interpolations.cubic_spline_interpolation(x, y)
        x_new = range(0, 1, new_num_points)
        y_new = itp_cubic(x_new)
        y_new[1]=0

        # Find PDF
        # compute histogram
        left_edge = 72.0 # hours
        if mode == "energy_demand"
            left_edge = 100.0 # kWh
        end
        hist = StatsBase.fit(StatsBase.Histogram, y_new, 0.0:1.0:left_edge)
        # Calculate the empirical density function
        interval_width = hist.edges[1][2] - hist.edges[1][1]
        empirical_density = hist.weights / sum(hist.weights) / interval_width
        # Plot the PDF
        #bar(hist.edges[1][1:end], empirical_density)
        #plot(range(0,72,73), empirical_density)

        # Add an extra element with value 0 at the first position to adjust the total length
        prepend!(empirical_density, 0)
        # Pad the list with missing values to match the number of rows in the DataFrame
        if mode == "connection_time"
            empirical_density = vcat(empirical_density, fill(missing, nrow(df) - length(empirical_density)))
        end
        
        df[!, new_cols_names[i]] .= empirical_density
    end
end

function reverse_column_order!(df)
    sort!(df.private)
    sort!(df.public)
    sort!(df.workplace)  
end


function generate_random_number(values::Union{Vector{Int64}, Vector{Float64}}, cdf::Vector{Float64})
    # Step 1: Generate a random number between 0 and 1
    r = rand()

    # Step 2: Find the corresponding value in the CDF column
    index = searchsortedfirst(cdf, r)

    # Step 3: Use the corresponding value in the original range
    sampled_value = values[index]

    return sampled_value
end

function moving_average(x, window_size)
# Function to calculate moving average
    result = similar(x, Float64)
    for i in 1:length(x)
        lower = max(1, i - window_size ÷ 2)
        upper = min(length(x), i + window_size ÷ 2)
        result[i] = mean(x[lower:upper])
    end
    return result
end