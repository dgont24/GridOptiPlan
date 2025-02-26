function active_network_configuration_unique_data(states::Dict{Symbol, Vector{Int64}}, ref::Dict{Symbol,<:Any})
    
    oltc_config = determine_config(states[:oltc], IncludeOLTCConfig, NoOLTCConfig)
    svr_config = determine_config(vcat(states[:svr1],states[:svr2]), IncludeSVRConfig, NoSVRConfig)
    bat_config = determine_config(states[:batteries], IncludeBatteryConfig, NoBatteryConfig)

    # modify initial lines data with current network configuration
    lines_df = copy(ref[:lines_init])
    modify_current_lines!(states[:lines], lines_df, ref[:linecodes], ref[:upgrades_list][:lines], ref[:params]["S_base"], ref[:params]["Z_base"], ref[:params]["V_NOM"])

    oltc_data = create_oltc(states[:oltc], ref[:upgrades_list][:oltc], ref[:oltctypes], oltc_config)
    svr_1, id_end1 = create_svr(states[:svr1], ref[:upgrades_list][:svr1], ref[:svrtypes], 1)
    svr_2, _ = create_svr(states[:svr2], ref[:upgrades_list][:svr2], ref[:svrtypes], id_end1)
    svr_data = vcat(svr_1, svr_2)
    batteries_data = create_battery(states[:batteries], ref[:upgrades_list][:batteries], ref[:battypes], bat_config, ref[:settings], ref[:params]["S_base"])

    # modify network lines df after adding the voltage regulators
    lines_df.newf_bus = add_VR_connections_to_lines(lines_df, svr_data)
    sourcebus = isa(oltc_config, IncludeOLTCConfig) ? 1001 : 0
    
    local_ref = Dict{Symbol, Any}(
        :oltc_config => oltc_config,
        :svr_config => svr_config,
        :bat_config => bat_config,

        :lines => lines_df,
        :oltc => oltc_data,
        :svr => svr_data,
        :batteries => batteries_data,
        :sourcebus => sourcebus
    )

    # Define sets
    local_ref[:sets] = create_optimization_sets(ref, local_ref)

    return local_ref
end

function create_oltc(
    states::Vector{Int}, 
    upgrades::Dict{Int,Vector{Int}},
    types_data::Dict{Int, TapChangerType},
    config::T
) where {T<:AbstractOLTCConfiguration}
    
    arr = TapChanger[]

    if !isa(config, NoOLTCConfig)
        id=1
        for (idx, value) in enumerate(states)
            if value > 0
                tapnode = 1000 + id
                hvnode = upgrades[idx][1]
                lvnode = upgrades[idx][2]
                typecode = upgrades[idx][3]
                data = types_data[typecode]

                push!(arr, TapChanger(id, tapnode, hvnode, lvnode, data))

                id += 1
            end
        end
    end

    return arr
end

function create_svr(
    states::Vector{Int}, 
    upgrades::Dict{Int,Vector{Int}},
    types_data::Dict{Int, TapChangerType},
    id_start::Int
)
    
    arr = TapChanger[]
    id_current = id_start
    
    if sum(states, init=0) > 0
        for value in states
            if value > 0
                tapnode = 2000 + id_current
                hvnode = upgrades[value][1]
                lvnode = upgrades[value][2]
                typecode = upgrades[value][3]
                data = types_data[typecode]

                push!(arr, TapChanger(id_current, tapnode, hvnode, lvnode, data))

                id_current += 1
            end
        end
    end
    
    return arr, id_current
end

function create_battery(
    states::Vector{Int}, 
    upgrades::Dict{Int,Vector{Int}},
    types_data::Dict{Int, BatteryType},
    ::NoBatteryConfig,
    settings::OptiSettings,
    sbase::Float64;
    sensitivity::Float64=0.1
)

    return Dict{Int,Battery}()
end

function create_battery(
    states::Vector{Int}, 
    upgrades::Dict{Int,Vector{Int}},
    types_data::Dict{Int, BatteryType},
    ::IncludeBatteryConfig,
    settings::OptiSettings,
    sbase::Float64;
    sensitivity::Float64=0.1
)
    arr = Dict{Int,Battery}()
    sorted_keys::Vector{Int} = sort(collect(keys(upgrades)))

    id = 1
    for (idx,capacity) in enumerate(states)
        node = sorted_keys[idx]
        upper_bound = upgrades[node][2]

        # Define a boolean variable to check if the battery should be included or not.
        # Because the battery search space is continuous in a range eg[0-150 kWh] it is very unlikely that the capacity will be 0.
        # So for a range of values near zero, it is considered that the battery is not included.
        # This range is [0, upper_bound * sensitivity] eg [0, 150*0,1] = [0,15]
        includebat = capacity > (upper_bound * sensitivity)
        
        if includebat
            typecode = upgrades[node][3]
            soc_start = upgrades[node][4] / 100 # value is a percentage so that it is an integer.
            type = types_data[typecode]
            capacity_pu = capacity * 1000 / sbase # convert capacity (kWh) to pu/h.
            capital_cost = type.costcapacity_eur_per_kwh * capacity
            operational_cost = type.costmaint_eur_per_year_per_kWh * capacity
            # TODO: operational cost should be adapted here because it is based on the capacity in kWh and not in pu
            w = operational_cost / (type.cycles * (type.soc_max - type.soc_min) * capacity_pu)
            arr[node] = Battery(id, node, capacity, capacity_pu, Ref(Float64(capacity_pu)), Ref(Float64(type.efficiency)), soc_start, w, type, capital_cost, operational_cost)

            id += 1
        end
    end

    return arr
end   

# Define helper function to determine configuration of network type
function determine_config(
    states::Vector{T}, 
    include_config::Type{<:AbstractTypeConfiguration}, 
    no_config::Type{<:AbstractTypeConfiguration}
) where {T<:Real}

    if sum(states, init=0) > 0
        return include_config()
    else
        return no_config()
    end
end

# TODO: Change lines_df to be a struct instead of an df, as well as linecodes
# modify network lines based on network configuration given from genetic algorithm
function modify_current_lines!(
    line_states::Vector{Int},
    lines_df::DataFrames.DataFrame,
    linecodes::DataFrames.DataFrame,
    line_upgrades_list::Dict{Int,Vector{String}},
    sbase::Float64, 
    zbase::Float64, 
    vnom::Float64  
)
    sorted_ids = sort(collect(keys(line_upgrades_list)))
    for (i,id) in enumerate(sorted_ids)
        if line_states[i] != 0
            new_linecode = line_upgrades_list[id][line_states[i]]
            lines_df[id, :linecode] = new_linecode
            new_params = line_params_pu(lines_df[id,:], linecodes, sbase, zbase, vnom)
            lines_df[id, :r_pu] = new_params[1]
            lines_df[id, :x_pu] = new_params[2]
            lines_df[id, :smax_pu] = new_params[3]
        end
    end
end

# modify network lines after adding the voltage regulators
function add_VR_connections_to_lines(lines_df::DataFrames.DataFrame, svr_obj::Vector{TapChanger})
    newf_bus::Vector{Int} = lines_df[:, :from_bus]

    if !isempty(svr_obj)
        for obj in svr_obj
            mask = (lines_df.from_bus .== obj.hv_node) .& (lines_df.to_bus .== obj.lv_node)
            lineidx = findfirst(mask)
            newf_bus[lineidx] = obj.tap_node
        end
    end

    return newf_bus
end

function node_daily_timeseries(
    profiles::Dict{Int, DataFrames.DataFrame}, 
    N_X::Vector{Int}, 
    day::Int
)
    # Active power injection at each node
    p_dict = Dict{Int, Vector{Float64}}()

    for i in N_X
        df = profiles[i]
        # Store the profile in the dictionary with the node index as the key
        p_dict[i] = df[df.date .== day, :p_pu]
    end

    return p_dict
end