function load_data_from_files(path::String, settings::OptiSettings)
    # simulation global parameters
    params = create_global_parameters()
    clusters = load_clusters(joinpath(path,"clusters.csv"), joinpath(path,"clusters_representative.json"))
    days = create_representative_days_info(clusters)

    # grid data
    source_df = load_df(joinpath(path,"network","source.csv"), _source_columns)
    transformer_df = load_trafo(joinpath(path,"network","transformer.csv"), params["S_base"])
    buses_df = load_buses(joinpath(path,"network","buses.csv"))
    linecodes_df = load_df(joinpath(path,"equipment","linecodes.csv"), _linecode_columns)
    lines_df = load_lines(joinpath(path,"network","lines.csv"), linecodes_df, params)

    # create network directed graph
    node_parentsDict, node_childsDict = create_graphDict(lines_df.from_bus, lines_df.to_bus)

    # node injection data
    loads_df, load_profiles = load_node_data(joinpath(path, "load_data"), "load", params["Resolution"], _load_columns, _load_profile_columns, params["S_base"])  
    evs_df, ev_profiles = load_node_data(joinpath(path, "ev_data"), "ev", params["Resolution"], _ev_columns, _ev_profile_columns, params["S_base"])  
    gens_df, gen_profiles = load_node_data(joinpath(path, "pv_data"), "pv", params["Resolution"], _pv_columns, _pv_profile_columns, params["S_base"])  

    # load equipment types
    oltctypes = load_tap_changer_types(joinpath(path,"equipment","oltccodes.csv"))
    svrtypes = load_tap_changer_types(joinpath(path,"equipment","svrcodes.csv"))
    battypes = load_battery_types(joinpath(path,"equipment","batterycodes.csv"))

    # costs
    costs = load_cost_data(params["S_base"])

    # upgrades
    network_upgrades = load_upgrades(joinpath(path, "upgrades", "upgrades.json"), settings.upgrades, loads_df.bus)
    # Named tuple with the number of possible upgrades for each category
    upgrades_nitems = items_per_category(network_upgrades)

    # simulation outputs
    output_arr = Array{Any}(undef, settings.population_size, settings.iterations, 6)

    # create ref dictionary
    ref = Dict{Symbol, Any}(
        :settings => settings,
        :params => params,
        :clusters => clusters,
        :days => days,

        :source => source_df,
        :trafo => transformer_df,
        :buses => buses_df,
        :linecodes => linecodes_df,
        :lines_init => lines_df,

        :node_parents => node_parentsDict,
        :node_childs => node_childsDict,

        :loads => loads_df,
        :load_profiles => load_profiles,
        :evs => evs_df,
        :ev_profiles => ev_profiles,
        :gens => gens_df,
        :gen_profiles => gen_profiles,

        :oltctypes => oltctypes,
        :svrtypes => svrtypes,
        :battypes => battypes,

        :costs => costs,

        :upgrades_list => network_upgrades,
        :upgrades_nitems => upgrades_nitems,

        :output => output_arr
    )
    ref[:pshift] = ev_pshift(evs_df, params["S_base"])
    ref[:p_g_nom] = pvs_nominal_power(gens_df, params["S_base"])

    return ref
end


function load_df(path::String)
    return CSV.read(path, DataFrames.DataFrame)
end

function load_df(path::String, columns_to_keep::Vector{Tuple{String, DataType}})
    df = CSV.read(path, DataFrames.DataFrame)
    
    columns = [Symbol(i[1]) for i in columns_to_keep]
    DataFrames.select!(df, columns)

    return df
end

# buses
function load_buses(path::String)
# Slack bus: id ->0. Rest of the buses: id -> 1,2.....
    buses_df = load_df(path, _bus_columns)

    # skip the slack bus
    first_id = buses_df[1, :id]::Int
    if first_id == 0
        buses_df = buses_df[2:end,:]
    end

    return buses_df
end

# lines
function load_lines(path::String, linecodes::DataFrames.DataFrame, params::Dict{String, Real})
    lines_df = load_df(path, _line_columns)
    lines_df.id = 1:DataFrames.nrow(lines_df)

    #impedance values and thermal limits of each line in pu
    lineparams = line_params_pu(lines_df, linecodes, params["S_base"], params["Z_base"], params["V_NOM"])
    lines_df.r_pu = Float64[line[1] for line in lineparams]
    lines_df.x_pu = Float64[line[2] for line in lineparams]
    lines_df.smax_pu = Float64[line[3] for line in lineparams]

    return lines_df
end

# line parameters
function line_params_pu(data::T, linecodes::DataFrames.DataFrame, sbase::Float64, zbase::Float64, vnom::Float64) where {T<:DataFrames.DataFrameRow}
    r_per_km::Float64 = linecodes[linecodes.name .== data.linecode, :r1][1] # ohm/km
    x_per_km::Float64 = linecodes[linecodes.name .== data.linecode, :x1][1] # ohm/km
    i_max_A::Float64  = linecodes[linecodes.name .== data.linecode, :ampacity_a][1] # A 
    length_km::Float64 = 4 * data.length * 1e-3   # km

    r_pu::Float64 = r_per_km * length_km / zbase # pu
    x_pu::Float64 = x_per_km * length_km / zbase # pu
    s_max_pu::Float64 = floor(sqrt(3) * vnom * i_max_A / sbase) # pu

    return (r_pu, x_pu, s_max_pu)
end

function line_params_pu(df::DataFrames.DataFrame, linecodes::DataFrames.DataFrame, sbase::Float64, zbase::Float64, vnom::Float64)
    lineparams = Vector{Tuple{Float64, Float64, Float64}}([line_params_pu(dfr, linecodes, sbase, zbase, vnom) for dfr in eachrow(df)])
    return lineparams
end

# transformer
function load_trafo(path::String, s_base::Float64)
    trafo_df = load_df(path, _trafo_columns)

    # rated power
    trafo_df.s_pu = trafo_df.mva * 1e6 / s_base

    # series reactance and resistance
    s_btr_old = trafo_df.mva * 1e6  # VA
    # the values below seem to be percentage value sto I have to devide by 100 to get the actual ones.
    trafo_df.r_pu = trafo_df.r_100 .* (s_base ./ s_btr_old) / 100
    trafo_df.x_pu = trafo_df.x_100 .* (s_base ./ s_btr_old) / 100

    return trafo_df
end

# loads, evs, generators(pv)
function load_node_data(
    path::String,
    type::String,
    resolution::Int, 
    info_df_columns::Vector{Tuple{String, DataType}},
    profile_df_columns::Vector{Tuple{String, DataType}},
    sbase::Float64
)
    info_df = load_df(joinpath(path, type * "_info.csv"), info_df_columns)
    profiles = Dict{Int, DataFrames.DataFrame}()

    for row in eachrow(info_df)
        node_profile = load_df(joinpath(path, "$(resolution)min", row[:profile_name]), profile_df_columns)
        if type == "load"
            node_profile.consumption_kW .*= 1e3 / sbase # read consumption data in pu
            DataFrames.rename!(node_profile, Dict(:consumption_kW => "p_pu"))
        elseif type == "ev"
            node_profile.charging_power_kW .*= 1e3 / sbase # read consumption data in pu
            DataFrames.rename!(node_profile, Dict(:charging_power_kW => "p_pu"))
        elseif type == "pv"
            node_profile.p_w .*= 1 / sbase # read consumption data in pu
            DataFrames.rename!(node_profile, Dict(:p_w => "p_pu"))
        end

        profiles[row.bus] = node_profile 
    end

    return info_df, profiles
end

function load_tap_changer_types(path::String)
    df = load_df(path, _tap_changer_columns)

    # Array of Structs
    types_dict = Dict{Int, TapChangerType}(row.code => TapChangerType(row) for row in eachrow(df))
    
    return types_dict
end

function load_battery_types(path::String)
    df = load_df(path, _batterycode_columns)

    # Array of Structs
    types_dict = Dict{Int, BatteryType}(row.code => BatteryType(row) for row in eachrow(df))
    
    return types_dict
end

function load_clusters(filename_df::String, filename_json::String)
    clusters_df = load_df(filename_df, _cluster_columns)
    clusters_represent_days = JSON3.read(filename_json, Dict{Int, Dict{String, Vector{Int}}})

    # for the year find the min year in the df
    min_year = minimum(clusters_df[:,:year])
    max_year = maximum(clusters_df[:,:year])

    # initialize dictionary to store the data for all clusters for each year
    all_clusters_data = Dict{Int, Dict{String, Cluster}}()
    for year in min_year:max_year
        active_clusters = unique(clusters_df[clusters_df.year .== year, :type])
        setdiff!(active_clusters, _normal_cluster_types)
        @assert all(t -> t in _violation_cluster_types, active_clusters)

        if !isempty(active_clusters)
            all_clusters_data[year] = Dict(type => cluster_data_from_df(clusters_df, clusters_represent_days, year, type) for type in active_clusters)
        end
    end

    return all_clusters_data
end

function cluster_data_from_df(
    clusters_df::DataFrames.AbstractDataFrame, 
    all_representative_days::AbstractDict, 
    year::Int, 
    name::AbstractString
)
    # exctract rows of the current cluster
    mask  = (clusters_df.year .== year) .&& (clusters_df.type .== name)
    current_cluster_df = clusters_df[mask, [:id, :f_value]]

    # create data dict for each day
    days_dict = Dict(row.id => DayData(year, row.id, name, row.f_value) for row in eachrow(current_cluster_df))
    
    # find representative days
    representative_days = all_representative_days[year][name]
    len_rep_days = length(representative_days)
    rep_days_dict = Dict{Int, String}()
    for (idx,day) in enumerate(representative_days)
        rep_days_dict[day] = "$(idx)/$(len_rep_days)"
    end

    cluster_data  = Cluster(year, name, days_dict, rep_days_dict, representative_days)

    return cluster_data
end

function create_representative_days_info(clusters::Dict{Int, Dict{String, Cluster}})
    min_year::Int = minimum(keys(clusters))
    max_year::Int = maximum(keys(clusters))
    
    days = Dict(i => Dict{Int, RepDayData}() for i in min_year:max_year)
    for (year, year_cluster_data) in clusters
        for (cluster_name, cluster) in year_cluster_data
            for (day,type) in cluster.repdays_dict
                fval = cluster.days[day].fvalue
                days[year][day] = RepDayData(year, day, cluster_name, type, fval)
            end
        end
    end

    return days
end

# Read upgrades data from JSON file
# Define custom types that reflect the JSON structure
struct OltcJSON
    oltc::Dict{Int, Vector{Int}}
end

struct Svr1JSON
    svr1::Dict{Int, Vector{Int}}
end

struct Svr2JSON
    svr2::Dict{Int, Vector{Int}}
end

struct LineJSON
    lines::Dict{Int, Vector{String}}
end

struct BatteryJSON
    batteries::Dict{Int, Vector{Int}}
end

function load_upgrades(filename::String, upgrade_types::Vector{Symbol}, buses_with_loads::Vector{Int64})
    # Read JSON data from file
    json_data = read(filename, String)

    dict = Dict{Symbol, Dict{Int}}()

    if :oltc in upgrade_types
        oltc = JSON3.read(json_data, OltcJSON)
        dict[:oltc] = oltc.oltc
    else
        dict[:oltc] = Dict{Int,Vector{Int}}()
    end
    if :svr1 in upgrade_types
        svr_1 = JSON3.read(json_data, Svr1JSON)
        dict[:svr1] = svr_1.svr1
    else
        dict[:svr1] = Dict{Int,Vector{Int}}()
    end
    if :svr2 in upgrade_types
        svr_2 = JSON3.read(json_data, Svr2JSON)
        dict[:svr2] = svr_2.svr2
    else
        dict[:svr2] = Dict{Int,Vector{Int}}()
    end
    if :lines in upgrade_types
        ln = JSON3.read(json_data, LineJSON)
        dict[:lines] = ln.lines
    end
    if :batteries in upgrade_types
        bat = JSON3.read(json_data, BatteryJSON)
        dict[:batteries] = bat.batteries
    end

    # Check if any of the battery nodes are also load nodes. This is not allowed due to the optimization formulation.
    if haskey(dict, :batteries)
        common_buses = intersect(collect(keys(dict[:batteries])), buses_with_loads)
        if !isempty(common_buses)
            error("Error: Battery placement node, cannot be a load node")
        end
    end

    return dict
end


"""
# It calculates the nominal power of pv, which is required to compute the reactive power limits
# However this is the maximum power output. The nominal power is already known from the installed capacity.
function pvs_nominal_power()
    p_g_nom = Dict{Int64, Float64}()
    for (key, df) in gen_profiles
        max_p_w = maximum(df[:, :p_w])
        p_g_nom[key] = max_p_w / params["S_base"]
    end 

    return p_g_nom
end
"""
function pvs_nominal_power(info_df::DataFrames.DataFrame, s_base::Float64)
    p_g_nom = Dict{Int64, Float64}(info_df.bus .=> info_df.pv_size * 1e3 / s_base)

    return p_g_nom
end

# EV load shifting value
function ev_pshift(info_df::DataFrames.DataFrame, s_base::Float64)
    pshift = Dict{Int,Float64}(info_df.bus .=> (info_df.max_ac_charging_kW * 1e3 / s_base))

    return pshift
end