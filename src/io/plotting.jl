# Function to create a DataFrame based on id columns from a given DataFrame
# and a list of suffixes for each id.
function create_custom_df(ids::Vector{Int64}, suffixes::Vector{String})
    # Initial columns of the DataFrame
    columns = [:year, :cluster, :type, :day, :timestep]
    types = [Int64, String, String, Int64, Int64]
     
    # Dynamically add columns for each id with specified suffixes
    for suffix in suffixes
        for id in ids
            push!(columns, Symbol("$(suffix)_$(id)"))
            push!(types, Float64)  # Assuming all custom columns are Float64
        end
    end
    
    # Create an empty DataFrame with the specified columns
    df = DataFrames.DataFrame([Vector{T}(undef, 0) for T in types], columns)
    
    return df
end

# Function to create DataFrames to store the final values of all model variables
function _results_variables(
    s::Dict{Symbol, Vector{Int64}}, #sets
)
    # buses
    buses_suffixes = ["v_pu", "p_kw", "q_kw"]
    res_buses_df = create_custom_df(s[:N], buses_suffixes)

    # lines
    lines_suffixes = ["p_kw", "q_kw"]
    res_lines_df = create_custom_df(s[:LINES], lines_suffixes)

    # generators
    gens_suffixes = ["pinit_kW", "pcurt", "pgen", "qgen", "qctrl"]
    res_gens_df = create_custom_df(s[:G], gens_suffixes)

    # evs
    evs_suffixes = ["pinit_kW", "flex_up", "flex_dn", "pev_kW"]
    res_evs_df = create_custom_df(s[:EV], evs_suffixes)

    # batteries
    batteries_suffixes = ["energy_kWh", "pch_kW", "pdis_kW", "throughput_kWh"]
    res_batteries_df = create_custom_df(s[:B], batteries_suffixes)

    # oltc
    oltc_suffixes = ["tap", "tap_diff"]
    res_oltc_df = create_custom_df(s[:OLTC], oltc_suffixes)

    #svr
    svr_suffixes = ["tap", "tap_diff"]
    res_svr_df = create_custom_df(s[:SVR], svr_suffixes)

    res_dict = Dict{Symbol,DataFrames.DataFrame}(
        :buses => res_buses_df, :lines => res_lines_df, 
        :gens => res_gens_df, :evs => res_evs_df, :batteries => res_batteries_df, 
        :oltc => res_oltc_df, :svr => res_svr_df)

    return res_dict
end

# function to update the variables DataFrame with the results of the new optimization
function results_variables!(
    year_id::Int64, 
    day_id::Int64,
    model::JuMP.Model,
    ref::Dict{Symbol,<:Any},
    lref::Dict{Symbol,<:Any}, 
    res_df::Dict{Symbol, DataFrames.DataFrame}, 
    p_ev_pu::Dict{Int64, Vector{Float64}}, 
    p_gen_pu::Dict{Int64, Vector{Float64}}
)
    cluster_id::String = ref[:days][year_id][day_id].cluster
    typename::String = if day_id in ref[:clusters][year_id][cluster_id].repdays
        ref[:days][year_id][day_id].type
    else
        "other"
    end
    s::Dict{Symbol, Vector{Int64}} = lref[:sets]
    params::Dict{String, Real} = ref[:params]
    lines_df::DataFrames.DataFrame = lref[:lines]

    #buses
    bus_v = transpose(sqrt.([JuMP.value(model[:vVOLT2][i,t]) for i in s[:N], t in s[:T]])) #[pu]
    bus_p = transpose([JuMP.value(model[:PNODE][i,t]) for i in s[:N], t in s[:T]]) * params["S_base"]/1000  #[kW]
    bus_q = transpose([JuMP.value(model[:QNODE][i,t]) for i in s[:N], t in s[:T]]) * params["S_base"]/1000  #[kVA]
    for t in s[:T]
        push!(res_df[:buses], vcat(year_id, cluster_id, typename, day_id, t, bus_v[t,:], bus_p[t,:], bus_q[t,:]))
    end

    #lines
    ln_p = Array{Float64}(undef, length(s[:LINES]), length(s[:T]))
    ln_q = Array{Float64}(undef, length(s[:LINES]), length(s[:T]))
    for l in s[:LINES]
        for t in s[:T]
            ln_p[l,t] = JuMP.value(model[:PFLOW][lines_df[l,:from_bus],lines_df[l,:to_bus], t]) * params["S_base"]/1000  #[kW]
            ln_q[l,t] = JuMP.value(model[:QFLOW][lines_df[l,:from_bus],lines_df[l,:to_bus], t]) * params["S_base"]/1000  #[kVAr]
        end
    end
    ln_p = transpose(ln_p)
    ln_p = transpose(ln_q)
    for t in s[:T]
        push!(res_df[:lines], vcat(year_id, cluster_id, typename, day_id, t, ln_p[t,:], ln_p[t,:]))
    end

    #generators
    gen_pinit = reduce(hcat, [p_gen_pu[j] for j in s[:N_G]]) * params["S_base"]/1000  #[kW]
    gen_pcurt = transpose([JuMP.value(model[:PCURT][j,t]) for j in s[:N_G], t in s[:T]]) * params["S_base"]/1000  #[kW]
    gen_pgen  = transpose([JuMP.value(model[:PGEN][j,t])  for j in s[:N_G], t in s[:T]]) * params["S_base"]/1000  #[kW]
    gen_qgen  = transpose([JuMP.value(model[:QGEN][j,t])  for j in s[:N_G], t in s[:T]]) * params["S_base"]/1000  #[kVAr]
    gen_qctrl = transpose([JuMP.value(model[:QCTRL][j,t]) for j in s[:N_G], t in s[:T]]) * params["S_base"]/1000  #[kVAr]
    for t in s[:T]
        push!(res_df[:gens], vcat(year_id, cluster_id, typename, day_id, t, 
                                    gen_pinit[t,:], gen_pcurt[t,:], gen_pgen[t,:],
                                    gen_qgen[t,:], gen_qctrl[t,:]))
    end

    #evs
    ev_pinit  = reduce(hcat, [p_ev_pu[j] for j in s[:N_EV]]) * params["S_base"]/1000  #[kW]
    ev_flexup = transpose([JuMP.value(model[:flex_evUp][j,t]) for j in s[:N_EV], t in s[:T]])
    ev_flexdn = transpose([JuMP.value(model[:flex_evDn][j,t]) for j in s[:N_EV], t in s[:T]])
    ev_pev    = transpose([JuMP.value(model[:PEV][j,t]) for j in s[:N_EV], t in s[:T]]) * params["S_base"]/1000  #[kW]
    for t in s[:T]
        push!(res_df[:evs], vcat(year_id, cluster_id, typename, day_id, t, ev_pinit[t,:], ev_flexup[t,:], ev_flexdn[t,:], ev_pev[t,:]))
    end

    #batteries
    if isa(lref[:bat_config], IncludeBatteryConfig)
        bat_e =     transpose([JuMP.value(model[:EBAT][j,t])       for j in s[:N_B], t in s[:T]]) * params["S_base"]/1000  #[kWh]
        bat_pch =   transpose([JuMP.value(model[:PBAT_ch][j,t])    for j in s[:N_B], t in s[:T]]) * params["S_base"]/1000  #[kW]
        bat_pdis =  transpose([JuMP.value(model[:PBAT_dis][j,t])   for j in s[:N_B], t in s[:T]]) * params["S_base"]/1000  #[kW]
        bat_thr = transpose([JuMP.value(model[:EBAT_THR][j]) for j in s[:N_B]]) * params["S_base"]/1000  #[kWh]
        bat_thr = repeat(bat_thr, outer = length(s[:T]))
        for t in s[:T]
            push!(res_df[:batteries], vcat(year_id, cluster_id, typename, day_id, t, bat_e[t,:], bat_pch[t,:], bat_pdis[t,:], bat_thr[t,:]))
        end
    end

    #oltc
    if isa(lref[:oltc_config], IncludeOLTCConfig)
        oltc_tap = transpose([JuMP.value(model[:vNOLTC][j,t]) for j in s[:OLTC], t in s[:T]])
        oltc_tapDiff = transpose([JuMP.value(model[:vOTAPdiff][j,t]) for j in s[:OLTC], t in s[:T][2:end]])
        # add 0 in the row befoe oltc_tapDiff to account for the zero tap difference at the beginning
        oltc_tapDiff = vcat(zeros(1, length(s[:OLTC])), oltc_tapDiff)
        for t in s[:T]
            push!(res_df[:oltc], vcat(year_id, cluster_id, typename, day_id, t, oltc_tap[t,:], oltc_tapDiff[t,:]))
        end
    end

    #svr
    if isa(lref[:svr_config], IncludeSVRConfig)
        svr_tap = transpose([JuMP.value(model[:vNSVR][j,t]) for j in s[:SVR], t in s[:T]])
        svr_tapDiff = transpose([JuMP.value(model[:vSTAPdiff][j,t]) for j in s[:SVR], t in s[:T][2:end]])
        svr_tapDiff = vcat(zeros(1, length(s[:SVR])), svr_tapDiff)
        for t in s[:T]
            push!(res_df[:svr], vcat(year_id, cluster_id, typename, day_id, t, svr_tap[t,:], svr_tapDiff[t,:]))
        end
    end

    return
end

# Function to create a DataFrame to store the values of each term in the objective function
function _results_objective(
    s::Dict{Symbol, Vector{Int64}} #sets
)
    # Basic columns
    columns = [:year, :cluster, :type, :day, :tot]
    types = [Int64, String, String, Int64, Float64]
    
    # Add oltc columns
    for id in s[:OLTC]
        push!(columns, Symbol("oltc_$(id)"))
        push!(types, Float64)
    end

    # Add svr_n columns
    for id in s[:SVR]
        push!(columns, Symbol("svr_$(id)"))
        push!(types, Float64)
    end
    
    # Add bat_n columns
    for id in s[:B]
        push!(columns, Symbol("bat_$(id)"))
        push!(types, Float64)
    end

    # Add total flexibility costs
    push!(columns, Symbol("apc_tot_cost"))
    push!(columns, Symbol("rpc_tot_cost"))
    push!(columns, Symbol("evshift_tot_cost"))
    push!(types, Float64, Float64, Float64)
    
    # Add pcurt_n columns for gens_df
    for id in s[:G]
        push!(columns, Symbol("pcurt_$(id)"))
        push!(types, Float64)
    end
    # Add qctrl_n columns for gens_df
    for id in s[:G]
        push!(columns, Symbol("qctrl_$(id)"))
        push!(types, Float64)
    end
    
    # Add ev_shift_n columns for evs_df
    for id in s[:EV]
        push!(columns, Symbol("ev_shift_$(id)"))
        push!(types, Float64)
    end
    
    # Create an empty DataFrame with the specified columns and types
    df = DataFrames.DataFrame([Vector{T}(undef, 0) for T in types], columns)
    
    return df
end

# function to update the DataFrame which stores the values of the terms in the objective function
# function to update the variables DataFrame with the results of the new optimization
function results_objective!(
    year_id::Int64, 
    day_id::Int64,
    model::JuMP.Model,
    ref::Dict{Symbol,<:Any},
    lref::Dict{Symbol,<:Any},
    res_df::DataFrames.DataFrame
)
    cluster_id::String = ref[:days][year_id][day_id].cluster
    typename::String = if day_id in ref[:clusters][year_id][cluster_id].repdays
        ref[:days][year_id][day_id].type
    else
        "other"
    end
    s::Dict{Symbol, Vector{Int64}} = lref[:sets]
    c_p::Float64 = ref[:costs]["c_p"]
    c_q::Float64 = ref[:costs]["c_q"]
    c_ev::Float64 = ref[:costs]["c_ev"]
    p_shift::Dict{Int64, Float64} = ref[:pshift]
    ev_config::AbstractEVConfiguration = ref[:settings].ev_config

    total_cost = JuMP.objective_value(model)
    oltc_cost::Vector{Float64} = isa(lref[:oltc_config], IncludeOLTCConfig) ? [JuMP.value(model[:costOLTC][j]) for j in s[:OLTC]] : Float64[]
    svr_cost::Vector{Float64} = isa(lref[:svr_config], IncludeSVRConfig) ? [JuMP.value(model[:costSVR][j]) for j in s[:SVR]] : Float64[]
    bat_cost::Vector{Float64} = isa(lref[:bat_config], IncludeBatteryConfig) ? [JuMP.value(model[:bat_penalty][j]) for j in s[:N_B]] : Float64[]

    # PVs flex
    apc_cost = Vector{Float64}(undef, length(s[:G]))
    rpc_cost = Vector{Float64}(undef, length(s[:G]))
    for (idx, bus) in enumerate(s[:N_G])
        apc_cost[idx] = c_p * sum(JuMP.value(model[:PCURT][bus,t]) for t in s[:T]) * ref[:params]["Dt_hr_adjust"]
        rpc_cost[idx] = c_q * sum(JuMP.value(model[:QCTRL][bus,t]) for t in s[:T]) * ref[:params]["Dt_hr_adjust"]
    end
    apc_cost_tot = sum(apc_cost)
    rpc_cost_tot = sum(rpc_cost)

    # EVs flex
    evshift_cost = Vector{Float64}(undef, length(s[:EV]))
    for (idx, bus) in enumerate(s[:N_EV])
        if isa(ev_config, LinearEVConfig)
            evshift_cost[idx] = c_ev * p_shift[bus] * sum(JuMP.value(model[:flex_evDn][bus,t]) for t in s[:T]) * ref[:params]["Dt_hr_adjust"]
        elseif isa(ev_config, BinaryEVConfig)
            evshift_cost[idx] = c_ev * p_shift[bus] * sum(JuMP.value(model[:y_ev][bus,t]) for t in s[:T]) * ref[:params]["Dt_hr_adjust"]
        else
            error("Error: Unknown ev formulation. Cannot create result df of objective values")
        end
    end
    evshift_cost_tot = sum(evshift_cost)    

    push!(res_df, vcat(year_id, cluster_id, typename, day_id, total_cost, oltc_cost, svr_cost, bat_cost, 
                        apc_cost_tot, rpc_cost_tot, evshift_cost_tot, apc_cost, rpc_cost, evshift_cost))

    return
end

function save_results_to_folder(res_vars::Dict{Symbol, DataFrames.DataFrame}, res_obj::DataFrames.DataFrame, folder_path::String)
    # Ensure the folder exists
    if !isdir(folder_path)
        mkdir(folder_path)
    end

    # Save each DataFrame from the dictionary
    for (key, df) in res_vars
        file_path = joinpath(folder_path, string(key) * ".csv")
        CSV.write(file_path, df)
    end

    # Save the extra DataFrame
    objective_file_path = joinpath(folder_path, "objective" * ".csv")
    CSV.write(objective_file_path, res_obj)
end


#############################################
#################   EV   ####################
#############################################

#=
function value_to_color(val::Float64, min::Float64, max::Float64, palette, n_colors::Int64, showLimit::Float64, backgroundClr)
    if val < showLimit && val > -showLimit
        # Dont assign any color to arbitrary small values
        return backgroundClr
    else
        val_position = (val - min) / (max - min)  # position of value in the input range, relative to the length of the input range
        ind = round(Int, val_position * (n_colors - 1)) + 1  # target index in the color palette
        return palette[ind]
    end
end

function value_color_map(val::Float64, showLimit::Float64)
    # Dont assign any color to arbitrary small values
    res = 0
    if val > showLimit || val < -showLimit
        res = val
    end
    return res
end
=#

function ev_value_to_size(combined_1d_arr)

    if combined_1d_arr == []
        return []
    else
        # change this number experimentally, if you modify the graph dimensions, so that it matches the square size for one timestep.
        max_marker_size = 5.2

        # take absolute value because for the size the sign doesn't affect
        abs_arr = abs.(combined_1d_arr)
        # normalize values with log base e. Mimium displayed value 0.01 will be almost the same (0.0099)
            # I add 1 so I don't get negative values. Initially the values shoudl range 0-22
        arr = log.(abs_arr .+ 1) 
        scale_factor = max_marker_size / maximum(arr)
        marker_size = arr * scale_factor

        return marker_size
    end
end

# Define the equivalent of the heatmap function
function ev_heatmap_daily(
    gdf::DataFrames.GroupedDataFrame,
    s::Dict{Symbol, Vector{Int64}}, #sets
    year_id::Int64, 
    day_id::Int64,
    time_step_size::Int64=5
)
    # exctract SubDataFrame from the GroupedDataFrame
    sdf = gdf[(year = year_id, day = day_id)]

    # extract flexibility values from DataFrame
    flexEV = [(sdf[:, Symbol("pev_kW_$i")] - sdf[:, Symbol("pinit_kW_$i")]) for i in s[:EV]]
    flexDN = deepcopy(flexEV)
    for vec in flexDN
        vec[vec .> 0] .= 0
    end
    flexUP = deepcopy(flexEV)
    for vec in flexUP
        vec[vec .< 0] .= 0
    end
    # interleave the two arrays above so that I have the following order of vectors: ev_1 flex_Dn in pos1, flex_Up in pos2, ev_2 flex_Dn in pos3 etc.
    flex_combine = collect(Iterators.flatten(zip(flexDN, flexUP))) 

    # combine all the vectors in a single vector
    flex_combine1D = vcat(flex_combine...)

    # create a vector with the corresponding timestep for each flexibility value 
    x_scatter = repeat(s[:T], 2*length(s[:EV]))
    # The y axis is also categorical (id's of evs). The flexibility value will only determine the size of the point
    y_up = s[:EV] .+ 0.22
    y_dn = s[:EV] .- 0.22
    y_combine = collect(Iterators.flatten(zip(y_dn, y_up)))
    y_scatter = repeat(y_combine, inner=length(s[:T]))

    # Combine into a new structure
    combined_data = [(x, y, z) for (x, y, z) in zip(x_scatter, y_scatter, flex_combine1D)]
    
    # Don't show values smaller than:
    show_limit = 0.01
    # Filter based on condition z greater than show_limit or less than -show_limit
    filtered_data = [(x, y, z) for (x, y, z) in combined_data if abs(z) > show_limit]
    x_filtered = [x[1] for x in filtered_data] 
    y_filtered = [y[2] for y in filtered_data] 
    z_filtered = [z[3] for z in filtered_data] 

    # define ticks for the scatteplot
    x_ticks = 1:time_step_size:s[:T][end]
    y_ticks = s[:EV]

    # set axis limits
    x_lims = (0.5, s[:T][end] + 0.5) 
    y_lims = (s[:EV][1] - 0.5, s[:EV][end] + 0.5) 

    background_color = :grey92
  
    # Add a color palette
    # Use 256 colors for the RdBu color palette
    #n_colors = 256
    #color_palette = reverse(ColorSchemes.diverging_bwr_40_95_c42_n256) # Create a diverging palette
    # Range of values that will be mapped to the palette
    #color_min = minimum(flex_combine1D)
    #color_max = maximum(flex_combine1D)  
    # Treating of very small values
    # don't show values smaller than:
    #show_limit = 0.01 
    #backgroung_color = ColorSchemes.RGB(0.94, 0.94, 0.94)
    # Assign the corresponding color to every point
    #colors = [value_to_color(flex_combine1D[i], color_min, color_max, color_palette, n_colors, show_limit, backgroung_color) for i in eachindex(flex_combine1D)]
    # Vector with final values that are used to color the markers 
    #pointZ = [value_color_map(flex_combine1D[i], show_limit) for i in eachindex(flex_combine1D)]

    scatter_plot = Plots.scatter(
        x_filtered,      
        y_filtered,
        marker_z = z_filtered,
        markercolor = :RdYlGn,
        #markersize = ev_value_to_size(z_filtered),  # Vector of square sizes, proportional to size parameter
        markershape = :square,  # Use square as scatterplot marker
        markerstrokewidth = 0,
        colorbar = :right,
        #colorbar_title = "Title",
        # adjust margins around the graph
        bottom_margin = (18, :px),
        left_margin = (25, :px),
        right_margin = (1, :px),
        top_margin = (5, :px),
        title = "EV Activated Flexibility for day $day_id (in kW)",
        xlabel = "Timestep",
        ylabel = "EV id",
        xticks = x_ticks,
        yticks = y_ticks,
        tick_dir = :out,
        xlims = x_lims,
        ylims = y_lims,
        legend = false,
        grid = false,
        size = (1280, 480),
        #background_color_outside = :grey79,
        background_color_inside = background_color,
        bordercolor = background_color        
    )
    # add grid lines
    h_lines = [i+0.5 for i in s[:EV][1:end-1]]
    v_lines = [i+0.5 for i in s[:T][1:end-1]]
    # between different evs
    Plots.hline!(scatter_plot, h_lines, linecolor = :white, linewidth = 2)
    # between timesteps
    Plots.vline!(scatter_plot, v_lines, linecolor = :white, linewidth = 0.5)
    # between pos and neg flexibility for the same ev
    Plots.hline!(scatter_plot, s[:EV], linecolor = :white, linewidth = 0.5)
end


#############################################
#################   PV   ####################
#############################################

function pv_value_to_size(combined_1d_arr)

    if combined_1d_arr == []
        return []
    else
        # change this number experimentally, if you modify the graph dimensions, so that it matches the square size for one timestep.
        max_marker_size = 5.2

        # take absolute value because for the size the sign doesn't affect
        abs_arr = abs.(combined_1d_arr)
        # normalize values with log base e. Mimium displayed value 0.01 will be almost the same (0.0099)
            # I add 1 so I don't get negative values. Initially the values shoudl range 0-22
        arr = log.(abs_arr .+ 1) 
        scale_factor = max_marker_size / maximum(arr)
        marker_size = arr * scale_factor

        return marker_size
    end
end

# Define the equivalent of the heatmap function
function pv_heatmap_daily(
    gdf::DataFrames.GroupedDataFrame,
    variable::Symbol, #  :pcurt or :qgen
    s::Dict{Symbol, Vector{Int64}}, #sets
    year_id::Int64, 
    day_id::Int64,
    time_step_size::Int64=5
)
    # exctract SubDataFrame from the GroupedDataFrame
    sdf = gdf[(year = year_id, day = day_id)]

    # extract flexibility values from DataFrame
    arr = [sdf[:, Symbol(string(variable, "_", i))] for i in s[:G]]
    # combine all the vectors in a single vector
    comb1d = vcat(arr...)

    # create a vector with the corresponding timestep for each flexibility value 
    x_scatter = repeat(s[:T], length(s[:G]))
    # The y axis is also categorical (id's of pvs). The flexibility value will only determine the size of the point
    y_scatter = repeat(s[:G], inner=length(s[:T]))

    # Combine into a new structure
    combined_data = [(x, y, z) for (x, y, z) in zip(x_scatter, y_scatter, comb1d)]
    
    # Don't show values smaller than:
    show_limit = 0.01
    # Filter based on condition z greater than show_limit or less than -show_limit
    filtered_data = [(x, y, z) for (x, y, z) in combined_data if abs(z) > show_limit]
    x_filtered = [x[1] for x in filtered_data] 
    y_filtered = [y[2] for y in filtered_data] 
    z_filtered = [z[3] for z in filtered_data] 

    # define ticks for the scatteplot
    x_ticks = 1:time_step_size:s[:T][end]
    y_ticks = 1:2:s[:G][end]

    # set axis limits
    x_lims = (0.5, s[:T][end] + 0.5) 
    y_lims = (s[:G][1] - 0.5, s[:G][end] + 0.5) 

    background_color = :grey92
    if variable == :pcurt
        graph_title = "PV Curtailed Power for day $day_id (in kW)"
        palette = :matter
    else
        graph_title = "PV Reactive Power for day $day_id (in kVAr)"
        palette = :diverging_rainbow_bgymr_45_85_c67_n256
    end
  
    scatter_plot = Plots.scatter(
        x_filtered,      
        y_filtered,
        marker_z = z_filtered,
        markercolor = palette,
        markersize = pv_value_to_size(z_filtered),  # Vector of square sizes, proportional to size parameter
        markershape = :square,  # Use square as scatterplot marker
        markerstrokewidth = 0,
        colorbar = :right,
        #colorbar_title = "Title",
        # adjust margins around the graph
        bottom_margin = (21, :px),
        left_margin = (25, :px),
        right_margin = (1, :px),
        top_margin = (5, :px),
        title = graph_title,
        xlabel = "Timestep",
        ylabel = "PV id",   
        xticks = x_ticks,
        yticks = y_ticks,
        tick_dir = :out,
        xlims = x_lims,
        ylims = y_lims,
        legend = false,
        grid = false,
        size = (1280, 400),
        #background_color_outside = :grey79,
        background_color_inside = background_color,
        bordercolor = background_color        
    )
    # add grid lines
    h_lines = [i+0.5 for i in s[:G][1:end-1]]
    v_lines = [i+0.5 for i in s[:T][1:end-1]]
    # between different pvs
    Plots.hline!(scatter_plot, h_lines, linecolor = :white, linewidth = 2)
    # between timesteps
    Plots.vline!(scatter_plot, v_lines, linecolor = :white, linewidth = 0.5)
end

##############################################
#################   BAT   ####################
##############################################

function battery_plot(
    bat_id::Int64,
    gdf::DataFrames.GroupedDataFrame,
    bat_info::Battery,
    s::Dict{Symbol, Vector{Int64}}, #sets
    year_id::Int64, 
    day_id::Int64
)
    # exctract SubDataFrame from the GroupedDataFrame
    sdf = gdf[(year = year_id, day = day_id)]

    x_plot = s[:T]
    E_plot = sdf[:, Symbol("energy_kWh_$bat_id")]
    Pch_plot = sdf[:, Symbol("pch_kW_$bat_id")]
    Pdis_plot = -sdf[:, Symbol("pdis_kW_$bat_id")]

    # set x-axis ticks
    x_ticks = 1:5:s[:T][end]
    # set axis limits
    x_lims = (s[:T][1], s[:T][end])
    y_lim_up = bat_info.capacity + 5
    y_lim_dn = floor(Int, minimum(Pdis_plot)) - 2 # I floor becasue it is negative

    # Power Plot -> left-side axis
    line_plot = Plots.plot(
        x_plot,      
        E_plot,
        title = "Battery $bat_id use for day $day_id",
        xlabel = "Timestep",
        ylabel = "Energy [kWh]", 
        xlim = x_lims,
        ylim = (y_lim_dn, y_lim_up),
        xticks = x_ticks,
        tick_dir = :out,
        # specify color
        linecolor = :dodgerblue1,
        # adjust margins around the graph
        #bottom_margin = (21, :px),
        #left_margin = (25, :px),
        #right_margin = (1, :px),
        #top_margin = (5, :px),
        legend = true,
        legend_position = :left,
        label= "energy",
        grid = true,
        size = (600, 400),
        #background_color_outside = :grey79,
        #background_color_inside = background_color,
        #bordercolor = background_color        
    )

    h_lines = bat_info.capacity * [1, bat_info.data.soc_min, bat_info.data.soc_max]
    Plots.hline!(line_plot, [h_lines[1]], linecolor = :red, linewidth = 1, label = "max capacity") 
    Plots.hline!(line_plot, h_lines[2:3], linecolor = :red4, linewidth = 0.5, label = "operating limits")

    Plots.plot!(
        Plots.twinx(), 
        [Pch_plot, Pdis_plot], 
        ylabel= "Power [kW]",
        xlim = x_lims,
        ylim = (y_lim_dn, y_lim_up),
        # specify colors
        cmap = [:springgreen3 :red2],
        # shade the area below the curve
        fill = true,
        fillalpha = 0.3,
        label= "power",  
    )    
end

##############################################
#################   TAPS   ###################
##############################################

function taps_plot(
    gdf_S::DataFrames.GroupedDataFrame,
    gdf_O::DataFrames.GroupedDataFrame,
    s::Dict{Symbol, Vector{Int64}}, #sets
    year_id::Int64, 
    day_id::Int64
)
    # exctract SubDataFrame from the GroupedDataFrame
    # sdf_S = gdf_S[(year = year_id, day = day_id)]
    sdf_O = gdf_O[(year = year_id, day = day_id)]

    x_plot = s[:T]
    tap_o1 = sdf_O[:, :tap_1]
    # tap_s1 = sdf_S[:, :tap_1]
    # tap_s2 = sdf_S[:, :tap_2]

    # Combine vectors just to calculate axis limits
    # combined = vcat(tap_s1, tap_s2, tap_o1)
    # max_y = maximum(combined) + 2
    # min_y = minimum(combined) - 2
    max_y = maximum(tap_o1) + 2
    min_y = minimum(tap_o1) - 2

    # set x-axis ticks
    x_ticks = 1:5:s[:T][end]
    # y_ticks = min_y:3:max_y
    # set axis limits
    x_lims = (s[:T][1], s[:T][end])
    y_lims = (min_y, max_y)

    # Power Plot -> left-side axis
    line_plot = Plots.plot(
        x_plot,      
        tap_o1, # [tap_o1, tap_s1, tap_s2],
        title = "Tap changes for day $day_id",
        xlabel = "Timestep",
        ylabel = "Tap position", 
        xlim = x_lims,
        ylim = y_lims,
        xticks = x_ticks,
        # yticks = y_ticks,
        tick_dir = :out,
        # Number of minor intervals between major ticks
        # yminorticks = 3,
        legend = true,
        label= "oltc tap position",
        grid = true,
        gridalpha = 0.2,
        size = (600, 400)    
    )
end


##############################################
################  Aggregate ##################
##############################################

function daily_overall_area_plot(
    gdfD::Dict{Symbol,DataFrames.GroupedDataFrame{DataFrames.DataFrame}},
    s::Dict{Symbol, Vector{Int64}}, #sets
    year_id::Int64, 
    day_id::Int64
)
    # exctract SubDataFrame from the GroupedDataFrame Dictionary
    sdf_ev = gdfD[:ev][(year = year_id, day = day_id)]
    sdf_pv = gdfD[:pv][(year = year_id, day = day_id)]
    sdf_bat = gdfD[:bat][(year = year_id, day = day_id)]

    x_plot = s[:T]
    # aggregate total "shifted" power from ev for each timestep. I plot only the curtailed power
        # maybe I should add also the other half cause it might be used by the optimization to help with grid contigencies
    flexDN = [(sdf_ev[:, Symbol("pinit_kW_$i")] - sdf_ev[:, Symbol("pev_kW_$i")]) for i in s[:EV]]
    flexDN_arr = hcat(flexDN...)
    flexDN_arr[flexDN_arr .< 0] .= 0
    flexDN_agg = [sum(flexDN_arr[i,:]) for i in s[:T]]

    # aggregate total curtailed power from pv for each timestep
    pcurt = [sdf_pv[:, Symbol("pcurt_$i")] for i in s[:G]]
    pcurt_arr = hcat(pcurt...)
    pcurt_agg = [sum(pcurt_arr[i,:]) for i in s[:T]]

    # aggregate total absolute reactive power from pv for each timestep
    qgen = [sdf_pv[:, Symbol("qgen_$i")] for i in s[:G]]
    qgen_arr = hcat(qgen...)
    qgen_arr_abs = abs.(qgen_arr)
    qgen_agg_abs = [sum(qgen_arr_abs[i,:]) for i in s[:T]]

    # aggregate total power from batteries for each timestep. 
        # There is no target final SoC. So whenever the battery is used is in order to help the network.
    # TODO: Fix that so the number of batteries is variable
    pbat_1 = sdf_bat[:, :pch_kW_1] + sdf_bat[:, :pdis_kW_1]
    # pbat_2 = sdf_bat[:, :pch_kW_2] + sdf_bat[:, :pdis_kW_2] 
    pbat = pbat_1 #+ pbat_2

    y4 = flexDN_agg + pcurt_agg + qgen_agg_abs + pbat

    # set x-axis ticks
    x_ticks = 1:5:s[:T][end]
    y_ticks = range(0, stop=maximum(y4), length=10)
    y_ticks_rounded = round.(Int, y_ticks)

    Plots.areaplot(
        x_plot, 
        [qgen_agg_abs pbat flexDN_agg pcurt_agg], 
        seriescolor = [1 2 3 4], 
        fillalpha = [0.35 0.35 0.35 0.35],
        title = "Aggregeated flexibility used for day $day_id",
        xlabel = "Timestep",
        ylabel = "Power [kW]",
        label = ["pv reactive power" "active power from ess" "shifted power from ev" "pv curtailed power"],
        legend = :topleft,
        xticks = x_ticks,
        yticks = y_ticks_rounded,
        tick_dir = :out,
        size = (600, 400)
    )
end

# more detailed graph
function daily_overall_area_plot2(
    gdfD::Dict{Symbol,DataFrames.GroupedDataFrame{DataFrames.DataFrame}},
    s::Dict{Symbol, Vector{Int64}}, #sets
    year_id::Int64, 
    day_id::Int64
)
    # exctract SubDataFrame from the GroupedDataFrame Dictionary
    sdf_ev = gdfD[:ev][(year = year_id, day = day_id)]
    sdf_pv = gdfD[:pv][(year = year_id, day = day_id)]
    sdf_bat = gdfD[:bat][(year = year_id, day = day_id)]

    x = s[:T] ./ 4 
    # aggregate total "shifted" power from ev for each timestep.
    # flex down 
    flexDN = [(-sdf_ev[:, Symbol("pinit_kW_$i")] + sdf_ev[:, Symbol("pev_kW_$i")]) for i in s[:EV]]
    flexDN_arr = hcat(flexDN...)
    flexDN_arr[flexDN_arr .> 0] .= 0
    flexDN_agg = [sum(flexDN_arr[i,:]) for i in s[:T]]

    # aggregate total "shifted" power from ev for each timestep. 
    # flex up
    flexUP = [(-sdf_ev[:, Symbol("pinit_kW_$i")] + sdf_ev[:, Symbol("pev_kW_$i")]) for i in s[:EV]]
    flexUP_arr = hcat(flexUP...)
    flexUP_arr[flexUP_arr .< 0] .= 0
    flexUP_agg = [sum(flexUP_arr[i,:]) for i in s[:T]]

    # aggregate total curtailed power from pv for each timestep
    pcurt = [sdf_pv[:, Symbol("pcurt_$i")] for i in s[:G]]
    pcurt_arr = hcat(pcurt...)
    pcurt_agg = -[sum(pcurt_arr[i,:]) for i in s[:T]]

    # aggregate total absolute reactive power from pv for each timestep
    qgen = [sdf_pv[:, Symbol("qgen_$i")] for i in s[:G]]
    qgen_arr = hcat(qgen...)
    qgen_arr_up = [x > 0 ? x : 0 for x in qgen_arr]
    qgen_agg_up = [sum(qgen_arr_up[i,:]) for i in s[:T]]

    # aggregate total absolute reactive power from pv for each timestep
    qgen = [sdf_pv[:, Symbol("qgen_$i")] for i in s[:G]]
    qgen_arr = hcat(qgen...)
    qgen_arr_dn = [x < 0 ? x : 0 for x in qgen_arr]
    qgen_agg_dn = [sum(qgen_arr_dn[i,:]) for i in s[:T]]

    # aggregate total power from batteries for each timestep. 
        # There is no target final SoC. So whenever the battery is used is in order to help the network.
    # TODO: Fix that so the number of batteries is variable
    # pbat_1 = sdf_bat[:, :pch_kW_1] + sdf_bat[:, :pdis_kW_1]
    # pbat_2 = sdf_bat[:, :pch_kW_2] + sdf_bat[:, :pdis_kW_2] 
    # pbat = pbat_1 #+ pbat_2
    pbat_ch = sdf_bat[:, :pch_kW_1]
    pbat_dis = -sdf_bat[:, :pdis_kW_1]

    max_y = flexUP_agg + qgen_agg_up + pbat_ch
    min_y = flexDN_agg + qgen_agg_dn + pbat_dis + pcurt_agg 

    # set x-axis ticks
    x_ticks = 0:24
    y_ticks = range(minimum(min_y), maximum(max_y), length=10)
    y_ticks_rounded = round.(Int, y_ticks)

    # Plots.areaplot(
    #     x_plot, 
    #     [qgen_agg_up qgen_agg_dn pbat_ch pbat_dis flexUP_agg flexDN_agg pcurt_agg], 
    #     seriescolor = [1 1 2 2 3 3 4], 
    #     fillalpha = [0.35 0.35 0.35 0.35 0.35 0.35 0.35],
    #     title = "Aggregeated flexibility used for day $day_id",
    #     xlabel = "Timestep",
    #     ylabel = "Power [kW]",
    #     label = ["pv reactive power" "" "active power from ess" "" "shifted power from ev" "" "pv curtailed power"],
    #     legend = :topleft,
    #     xticks = x_ticks,
    #     yticks = y_ticks_rounded,
    #     tick_dir = :out,
    #     size = (600, 400)
    # )

    y_pos_1 = abs.(pbat_dis) 
    y_neg_1 = -abs.(pbat_ch) 

    y_pos_2 = y_pos_1 .+ abs.(qgen_agg_up)
    y_neg_2 = y_neg_1 .- abs.(qgen_agg_dn)

    y_pos_3 = y_pos_2 .+ abs.(flexDN_agg)
    y_neg_3 = y_neg_2 .- abs.(flexUP_agg)

    y_neg_4 = y_neg_3 .- abs.(pcurt_agg)

    # Plot the stacked area plot
    Plots.plot(
        x, 
        y_pos_1, 
        fillrange=0, 
        seriestype=:path, 
        linealpha=0, 
        fillalpha=0.5, 
        label="BES", 
        seriescolor=3, 
        gridalpha=0.1,
        gridlinewidth=0.5,
        # gridcolor=:lightgrey, 
        gridstyle=:dash,
        framestyle = :box)
    Plots.plot!(x, y_pos_2, fillrange=y_pos_1, seriestype=:path, linealpha=0, fillalpha=0.5, label="RPC", seriescolor=1)
    Plots.plot!(x, y_pos_3, fillrange=y_pos_2, seriestype=:path, linealpha=0, fillalpha=0.5, label="EV", seriescolor=2)

    Plots.plot!(x, y_neg_1, fillrange=0, seriestype=:path, linealpha=0, fillalpha=0.5, label="", seriescolor=3)
    Plots.plot!(x, y_neg_2, fillrange=y_neg_1, seriestype=:path, linealpha=0, fillalpha=0.5, label="", seriescolor=1)
    Plots.plot!(x, y_neg_3, fillrange=y_neg_2, seriestype=:path, linealpha=0, fillalpha=0.5, label="", seriescolor=2)
    Plots.plot!(x, y_neg_4, fillrange=y_neg_3, seriestype=:path, linealpha=0, fillalpha=0.5, label="APC", seriescolor=4)

    Plots.hline!([0], linecolor = :black, linewidth = 0.5, label="")

    # Additional plot settings
    Plots.title!("Aggregate activated flexibility")
    Plots.xlabel!("Hour of Day")
    Plots.ylabel!("Power [kW]")
    Plots.xticks!(x_ticks)

    # Plots.grid!(true)  # Ensure the grid is enabled
    # Plots.gridcolor!(:gray)  # Set grid color to gray (or any color you prefer)
    # Plots.gridalpha!(0.7)  # Adjust the grid line transparency (0 is fully transparent, 1 is fully opaque)
    # Plots.gridlinewidth!(1)  # Increase the grid line width
end

function middle_element(vec::Vector)
    middle_index = div(length(vec), 2) + 1
    return middle_index
end

function all_days_bar_plot(
    ref::Dict{Symbol,<:Any},
    gdfD::Dict{Symbol,DataFrames.GroupedDataFrame{DataFrames.DataFrame}},
    s::Dict{Symbol, Vector{Int64}}, #sets
    params::Dict{String, Real}, #params
    year_id::Int64,
)
    days_ids = Int[]
    cluster_positions = Int[]
    current_cluster_end_position = 0
    vline_positions = Float64[]
    clusters_order = ["ov", "ov + ol", "ov + uv", "uv + ol", "uv", "ov + uv + ol"]
    for (j,cluster_name) in enumerate(clusters_order)
        cluster = ref[:clusters][year_id][cluster_name]
        append!(days_ids, cluster.repdays)
        middle_day_index = middle_element(cluster.repdays)
        if j == 1
            push!(cluster_positions, middle_day_index)
            current_cluster_end_position += length(cluster.repdays)
        else
            push!(cluster_positions, current_cluster_end_position + middle_day_index)
            current_cluster_end_position += length(cluster.repdays)
        end
        if j == 1
            push!(vline_positions, length(cluster.repdays) + 0.5)
        elseif j != length(clusters_order)
            push!(vline_positions, length(cluster.repdays) + vline_positions[j-1])
        end
    end

    n_days = length(days_ids)
    rpc = Vector{Float64}(undef, n_days)
    pbat = Vector{Float64}(undef, n_days)
    evflex = Vector{Float64}(undef, n_days)
    apc = Vector{Float64}(undef, n_days)

    for (i,day_id) in enumerate(days_ids)
        # exctract SubDataFrame from the GroupedDataFrame Dictionary
        sdf_ev = gdfD[:ev][(year = year_id, day = day_id)]
        sdf_pv = gdfD[:pv][(year = year_id, day = day_id)]
        sdf_bat = gdfD[:bat][(year = year_id, day = day_id)]

        # aggregate total "shifted" power from ev for each timestep. I plot only the curtailed power
            # maybe I should add also the other half cause it might be used by the optimization to help with grid contigencies
        flexUP = [(-sdf_ev[:, Symbol("pinit_kW_$i")] + sdf_ev[:, Symbol("pev_kW_$i")]) for i in s[:EV]]
        flexUP_arr = hcat(flexUP...)
        flexUP_arr[flexUP_arr .< 0] .= 0
        flexUP_arr = [sum(flexUP_arr[i,:]) for i in s[:T]]
        evflex[i] = sum(flexUP_arr) * params["Dt_hr_adjust"]

        # aggregate total curtailed power from pv for each timestep
        pcurt = [sdf_pv[:, Symbol("pcurt_$i")] for i in s[:G]]
        pcurt_arr = hcat(pcurt...)
        pcurt_agg = [sum(pcurt_arr[i,:]) for i in s[:T]]
        apc[i] = sum(pcurt_agg) * params["Dt_hr_adjust"]

        # aggregate total absolute reactive power from pv for each timestep
        qgen = [sdf_pv[:, Symbol("qgen_$i")] for i in s[:G]]
        qgen_arr = hcat(qgen...)
        qgen_arr_abs = abs.(qgen_arr)
        qgen_agg_abs = [sum(qgen_arr_abs[i,:]) for i in s[:T]]
        rpc[i] = sum(qgen_agg_abs) * params["Dt_hr_adjust"]

        # aggregate total power from batteries for each timestep. 
            # There is no target final SoC. So whenever the battery is used is in order to help the network.
        pbat_1 = sdf_bat[:, :pch_kW_1] + sdf_bat[:, :pdis_kW_1]
        pbat[i] = sum(pbat_1) * params["Dt_hr_adjust"]  
    end

    # Create a stacked bar plot
    StatsPlots.groupedbar(
        [apc evflex pbat rpc],
        label = ["APC" "EV" "BES" "RPC"],
        title = "Total flexibility used per day",
        # xlabel = "Clusters",
        ylabel = "Energy [kWh]",
        xticks = (1:n_days, days_ids),
        bar_width = 0.7,  
        bar_position = :stack,  # Stack the bars
        lw = 0,  # Remove bar borders
        framestyle = :box,
        legend = :topleft,
        #fillcolor = [4 2 3 1]  # Colors for each group
        colour = [4 2 3 1],  # Colors for each group
        xrotation = 45,
        bottom_margin = (7, :mm)  # Increase the bottom margin
    )
    Plots.vline!(vline_positions, lw=1, lc=:black, label="")
    labels = Tuple{String, Int, Symbol, Symbol}[]
    for cluster_name in clusters_order
        push!(labels, (replace(cluster_name, " "=>""), 9, :black, :top))
    end
    Plots.annotate!(cluster_positions, -80*ones(length(clusters_order)), labels)
    Plots.annotate!([(n_days//2, -180, ("Days/Clusters", 11, :bottom, :black))])
end

##############################################
################  Objective ##################
##############################################

function daily_costs_bar_plot(
    ref::Dict{Symbol,<:Any},
    df::DataFrames.DataFrame,
    year_id::Int64
)
    days_ids = Int[]
    cluster_positions = Int[]
    current_cluster_end_position = 0
    vline_positions = Float64[]
    clusters_order = ["ov", "ov + ol", "ov + uv", "uv + ol", "uv", "ov + uv + ol"]
    for (j,cluster_name) in enumerate(clusters_order)
        cluster = ref[:clusters][year_id][cluster_name]
        append!(days_ids, cluster.repdays)
        middle_day_index = middle_element(cluster.repdays)
        if j == 1
            push!(cluster_positions, middle_day_index)
            current_cluster_end_position += length(cluster.repdays)
        else
            push!(cluster_positions, current_cluster_end_position + middle_day_index)
            current_cluster_end_position += length(cluster.repdays)
        end
        if j == 1
            push!(vline_positions, length(cluster.repdays) + 0.5)
        elseif j != length(clusters_order)
            push!(vline_positions, length(cluster.repdays) + vline_positions[j-1])
        end
    end

    n_days = length(days_ids)

    order_map = Dict(day => i for (i, day) in enumerate(days_ids))
    df.id = [order_map[day] for day in sort(days_ids)]
    # Sort the DataFrame based on the temporary order column
    df_sorted = sort(df, :id)

    oltc = df_sorted[:, :oltc_1]
    pbat = df_sorted[:, :bat_1]
    apc = df_sorted[:, :apc_tot_cost]
    rpc = df_sorted[:, :rpc_tot_cost]
    evshift = df_sorted[:, :evshift_tot_cost]
  
    # set y-axis ticks
    tot_cost = oltc + pbat + apc + rpc + evshift
    # FInd max and round to the closest multiple of ten
    ymax = round(maximum(tot_cost) / 10) * 10 #
    y_ticks = range(0, stop=ymax, length=10)
    y_ticks_int = round.(y_ticks)

    # Create a stacked bar plot
    StatsPlots.groupedbar(
        [apc evshift pbat rpc],
        label = ["APC" "EV" "BES" "RPC"],
        title = "Flexibility costs per day",
        # xlabel = "Day ID",
        ylabel = "Cost [€]",
        xticks = (1:n_days, days_ids),
        yticks = y_ticks_int,
        # Number of minor intervals between major ticks
        yminorticks = 2,
        bar_width = 0.7,  
        bar_position = :stack,  # Stack the bars
        lw = 0,  # Remove bar borders
        framestyle = :box,
        legend = :topleft,
        #fillcolor = [4 2 3 1]  # Colors for each group
        colour = [4 2 3 1],  # Colors for each group
        xrotation = 45,
        bottom_margin = (8, :mm)  # Increase the bottom margin
    )
    Plots.vline!(vline_positions, lw=1, lc=:black, label="")  # Customize line width and color as needed
    labels = Tuple{String, Int, Symbol, Symbol}[]
    for cluster_name in clusters_order
        push!(labels, (replace(cluster_name, " "=>""), 9, :black, :top))
    end
    Plots.annotate!(cluster_positions, -3.4*ones(length(clusters_order)), labels)
    Plots.annotate!([(n_days//2, -8, ("Days/Clusters", 11, :bottom, :black))])
end