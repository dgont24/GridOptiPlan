## Optimization modelling of the tap changers (On-Load Tap Changers/ Voltage Regulators) in the network

# create variables for the model
function add_tap_changer_variables(
    model::JuMP.Model, 
    obj::Vector{TapChanger},
    sTapCh::Vector{Int},
    sT::Vector{Int},
    typename::String
)
    if typename == "OLTC"
        JuMP.@variables(model, begin
            vNOLTC[r in sTapCh, t in sT], Int
            vYOLTC[r in sTapCh, h in 1:obj[r].data.ns_binary_length, t in sT]
            vXOLTC[r in sTapCh, h in 1:obj[r].data.ns_binary_length, t in sT], Bin
        end)
    elseif typename == "SVR"
        JuMP.@variables(model, begin
            vNSVR[r in sTapCh, t in sT], Int
            vYSVR[r in sTapCh, h in 1:obj[r].data.ns_binary_length, t in sT]
            vXSVR[r in sTapCh, h in 1:obj[r].data.ns_binary_length, t in sT], Bin
        end)
    else
        throw(ArgumentError("The provided typename should be either OLTC or SVR"))
    end

    # We add two extra variables to avoid the absolute value in the tap changes summation of the objective function
    if typename == "OLTC"
        JuMP.@variable(model,
            vOTAPdiff[r in sTapCh, t in sT[2:end]]
        )
    elseif typename == "SVR"
        JuMP.@variable(model,
            vSTAPdiff[r in sTapCh, t in sT[2:end]]
        )
    else
        throw(ArgumentError("The provided typename should be either OLTC or SVR"))
    end

    return nothing
end

function add_tap_changer_constraints(
    model::JuMP.Model,
    obj::Vector{TapChanger},
    sTapCh::Vector{Int},
    sT::Vector{Int},
    vmin2::Float64,
    vmax2::Float64,
    typename::String
)
    # variables
    if typename == "OLTC"
        vNTAP = model[:vNOLTC]
        vYTAP = model[:vYOLTC]
        vXTAP = model[:vXOLTC]
    elseif typename == "SVR"
        vNTAP = model[:vNSVR]
        vYTAP = model[:vYSVR]
        vXTAP = model[:vXSVR]
    else
        error("Unknown typename: $typename")
    end

    # variables
    sending_nodes = [item.hv_node for item in obj]
    tap_nodes = [item.tap_node for item in obj]
    rows_to_select = union(sending_nodes, tap_nodes)
    vVOLT2 = model[:vVOLT2][rows_to_select, :]

    # add constraint for voltage at tap changer branches
    JuMP.@constraint(model, [r in sTapCh, t in sT],
        vVOLT2[obj[r].hv_node, t] == obj[r].data.tmin^2 * vVOLT2[obj[r].tap_node, t] + 2 * obj[r].data.dt * sum(2^(h-1) * vYTAP[r,h,t] for h in 1:obj[r].data.ns_binary_length)
    )

    # add constraint for tap changes to be integer values in the set [-NS/2, NS/2]
    JuMP.@constraint(model, [r in sTapCh, t in sT],
        vNTAP[r,t] == sum(2^(h-1) * vXTAP[r,h,t] for h in 1:obj[r].data.ns_binary_length) - div(obj[r].data.ns, 2)
    )

    # add McCormick linearization constraints for the variable y=x*v_m
    JuMP.@constraints(model, begin 
        [r in sTapCh, h in 1:obj[r].data.ns_binary_length, t in sT],
            vYTAP[r,h,t] >= obj[r].data.tmin^2 * vmin2 / sqrt(2) * vXTAP[r,h,t]
        [r in sTapCh, h in 1:obj[r].data.ns_binary_length, t in sT],
            vYTAP[r,h,t] <= obj[r].data.tmax^2 * vmax2 / sqrt(2) * vXTAP[r,h,t]
        [r in sTapCh, h in 1:obj[r].data.ns_binary_length, t in sT],
            vVOLT2[obj[r].tap_node, t] - vYTAP[r,h,t] >= obj[r].data.tmin^2 * vmin2 / sqrt(2) * (1 - vXTAP[r,h,t])
        [r in sTapCh, h in 1:obj[r].data.ns_binary_length, t in sT],
            vVOLT2[obj[r].tap_node, t] - vYTAP[r,h,t] <= obj[r].data.tmax^2 * vmax2 / sqrt(2) * (1 - vXTAP[r,h,t])
    end)

    # Constraint for SVR maximum tap. I don't add it here so that the optimization will decide how many taps are needed.
    #cSVR_maxTap[r in sSVR, t in sT],
    #    eSTAP[r,t] <= div(obj[r].data.ns, 2) 
    #cSVR_minTap[r in sSVR, t in sT],
    #    eSTAP[r,t] >= -div(obj[r].data.ns, 2)

    return nothing
end

function add_tap_changer_objective(
    model::JuMP.Model,
    obj::Vector{TapChanger},
    sTapCh::Vector{Int},
    sT::Vector{Int},
    typename::String
)
    # variables
    if typename == "OLTC"
        vNTAP = model[:vNOLTC]
        vTAPdiff = model[:vOTAPdiff]
    elseif typename == "SVR"
        vNTAP = model[:vNSVR]
        vTAPdiff = model[:vSTAPdiff]
    else
        error("Unknown typename: $typename")
    end

    JuMP.@constraints(model, begin
        [r in sTapCh, t in sT[2:end]],
            vTAPdiff[r,t] >= vNTAP[r,t] - vNTAP[r,t-1]
        [r in sTapCh, t in sT[2:end]],
            vTAPdiff[r,t] >= -vNTAP[r,t] + vNTAP[r,t-1]
    end)

    if typename == "OLTC"
        JuMP.@expression(model, costOLTC[r in sTapCh],
            (obj[r].data.cinv / obj[r].data.max_n) * sum(vTAPdiff[r,t] for t in sT[2:end])
        )
    else  
        JuMP.@expression(model, costSVR[r in sTapCh],
            (obj[r].data.cinv / obj[r].data.max_n) * sum(vTAPdiff[r,t] for t in sT[2:end])
        )
    end

    ## TODO: above in the cost function I should add the penalty for the first tap difference from the tap of the previous optimization stage
    # However, as the days are not subsequent, I don't think there is any need for that.

    return nothing
end

function add_tap_changer_formulation(
    model::JuMP.Model,
    obj::Vector{TapChanger},
    sTapCh::Vector{Int},
    sT::Vector{Int},
    vmin2::Float64,
    vmax2::Float64, 
    type::T
) where {T<:Union{IncludeOLTCConfig,IncludeSVRConfig}}

    if isa(type, IncludeOLTCConfig)
        typename = "OLTC"
    elseif isa(type, IncludeSVRConfig)
        typename = "SVR"
    else
        error("Unknown typename: $typename")
    end

    # add variables
    add_tap_changer_variables(model, obj, sTapCh, sT, typename)

    # add constraints
    add_tap_changer_constraints(model, obj, sTapCh, sT, vmin2, vmax2, typename)

    # add objective expression
    add_tap_changer_objective(model, obj, sTapCh, sT, typename)

    return nothing
end

function add_tap_changer_formulation(
    model::JuMP.Model,
    obj::Vector{TapChanger},
    sTapCh::Vector{Int},
    sT::Vector{Int},
    vmin2::Float64,
    vmax2::Float64,
    type::T
) where {T<:Union{NoOLTCConfig,NoSVRConfig}}
  
    return nothing
end


struct TapObject
    oltc_tap_begin::Int64
    oltc_tap_end::Int64
    svr_tap_begin::Vector{Int64}
    svr_tap_end::Vector{Int64}
end

function extract_daily_taps(
    model::JuMP.Model, 
    s::Dict #sets
)

    oltc_tapStart = 0
    oltc_tapEnd = 0
    if s[:OLTC] != [0]
        oltc_tapStart = value(model[:oltc_tap][1]) 
        oltc_tapEnd = value(model[:oltc_tap][end])
    end

    svr_tapStart = zeros(Int64, length(sets[:SVR]))
    svr_tapEND = zeros(Int64, length(sets[:SVR]))
    if sSVR != [0]
        for i in sSVR
            svr_tapStart[i] = value(model[:vNTAP][i,1])
            svr_tapEND[i] = value(model[:vNTAP][i,end])
        end
    end

    daily_taps = TapObject(oltc_tapStart, oltc_tapEnd, svr_tapStart, svr_tapEND)
    return daily_taps
end

function create_taps_vector(x::Int64, y::Int64, z::Int64, n::Int64)
    part1 = data[1:div(n, 3)]
    part2 = data[div(n, 3)+1:div(2n, 3)]
    part3 = data[div(2n, 3)+1:end]

    taps_vector = vcat(fill(x, length(part1)), fill(y, length(part2)), fill(z, length(part3)))
    return taps_vector
end

# Function to interpolate values for a given cluster
function interpolate_taps(
    tap_objects::Dict{String, TapObject}, 
    days::Vector{Int64},
    s::Dict #sets
)

    n_days = length(days)
    # Extract the values for min, median, max days
    min_obj = tap_objects["min"]
    median_obj = tap_objects["median"]
    max_obj = tap_objects["max"]
    
    # Interpolate values for each variable
    interpolated_oltc_tap_begin = create_taps_vector(min_obj.oltc_tap_begin, median_obj.oltc_tap_begin, max_obj.oltc_tap_begin, n_days)
    interpolated_oltc_tap_end = create_taps_vector(min_obj.oltc_tap_end, median_obj.oltc_tap_end, max_obj.oltc_tap_end, n_days)
    
    interpolated_svr_tap_begin = []
    interpolated_svr_tap_end = []
    if sSVR == 0
        push!(interpolated_svr_tap_begin, zeros(Int64, n_days))
        push!(interpolated_svr_tap_end, zeros(Int64, n_days))
    else
        for i in sSVR
            push!(interpolated_svr_tap_begin, create_taps_vector(min_obj.svr_tap_begin[i], median_obj.svr_tap_begin[i], max_obj.svr_tap_begin[i], n_days))
            push!(interpolated_svr_tap_end, create_taps_vector(min_obj.svr_tap_end[i], median_obj.svr_tap_end[i], max_obj.svr_tap_end[i], n_days))
        end
    end
    
    return interpolated_oltc_tap_begin, interpolated_oltc_tap_end, interpolated_svr_tap_begin, interpolated_svr_tap_end
end

function _yearly_tap_costs(
    s::Dict #sets
)

    columns = [:day_id, :oltc_start, :oltc_end]
    types = [Int64, Int64, Int64]
    
    # Dynamically add columns for svr id with specified suffixes
    for svr_id in sSVR
        push!(columns, Symbol("svr_$(svr_id)_start"))
        push!(types, Int64)  # Assuming all custom columns are Float64
        push!(columns, Symbol("svr_$(svr_id)_end"))
        push!(types, Int64)  # Assuming all custom columns are Float64
    end
    
    df = DataFrames.DataFrame([Vector{T}(undef, 0) for T in types], columns)
    return df
end
#=
function yearly_tap_costs(
    svr_df::DataFrame, 
    oltc_df::DataFrame,
    taps_df::DataFrame,
    s::Dict #sets
)

    sort!(taps_df, :day_id)
    svr_tap_cost = 0 
    if sSVR != 0
        svr_tap_cost = [svr_df[r,:cmain] / svr_df[r,:maxN] for r in sSVR]
    end
    oltc_tap_cost = 0 
    if s[:OLTC] != 0
        oltc_tap_cost = [oltc_df[1,:cmain] / oltc_df[1,:maxN] for r in sSVR]
    end
    
    for i in 1:nrow(taps_df)

    return
end
=#