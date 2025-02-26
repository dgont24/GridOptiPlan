## Optimization modelling of the batteries in the network

function add_battery_variables(
    model::JuMP.Model,
    sN_B::Vector{Int},
    sT::Vector{Int}
)
    
    JuMP.@variables(model, begin         
        EBAT[j in sN_B, t in sT]            # Battery available energy capacity 
        PBAT_ch[j in sN_B, t in sT]  >= 0   # BESS charging power
        PBAT_dis[j in sN_B, t in sT] >= 0   # BESS discharge power
        EBAT_THR[j in sN_B] >=0             # Battery daily energy throughput
    end)

    return nothing
end

function add_battery_constraints(
    model::JuMP.Model,
    obj::Dict{Int,Battery},
    sN_B::Vector{Int},
    sT::Vector{Int},
    dt::Float64
) 

    # variables
    EBAT = model[:EBAT]
    PBAT_ch = model[:PBAT_ch]
    PBAT_dis = model[:PBAT_dis]
    EBAT_THR = model[:EBAT_THR]

    JuMP.@constraints(model, begin
        cBatSoEStart[j in sN_B],
            EBAT[j,1] == obj[j].soc_start * obj[j].capacity_pu_degraded[]
        cBatSoE[j in sN_B, t in sT[2:end]],
            EBAT[j,t] == EBAT[j, t-1] + (obj[j].efficiency_degraded[] * PBAT_ch[j,t] - 
                            PBAT_dis[j,t] / obj[j].efficiency_degraded[]) * dt
        cBatSoE_LB[j in sN_B, t in sT],
            EBAT[j,t] >= obj[j].data.soc_min * obj[j].capacity_pu_degraded[]
        cBatSoE_UB[j in sN_B, t in sT],
            EBAT[j,t] <= obj[j].data.soc_max * obj[j].capacity_pu_degraded[]
        cBatPch_UB[j in sN_B, t in sT],
            PBAT_ch[j,t] <= obj[j].capacity_pu_degraded[] / 3
        cBatPdis_UB[j in sN_B, t in sT],
            PBAT_dis[j,t] <= obj[j].capacity_pu_degraded[] / 3    
        cBatThroughput[j in sN_B],
            EBAT_THR[j] == sum(PBAT_ch[j,t] * obj[j].efficiency_degraded[] + PBAT_dis[j,t] / obj[j].efficiency_degraded[] for t in sT) * dt
        cBatThroughput_UB[j in sN_B],
            # allow up to 1.5 cycles per day.
            EBAT_THR[j] <= 2 * (obj[j].data.soc_max - obj[j].data.soc_min) * obj[j].capacity_pu_degraded[]
        #TODO: Add correct maximum allowed charging and discharging power. It will be on the datasheet or be based on the C-rate of the battery.
    end)


    return nothing
end

function add_battery_objective(
    model::JuMP.Model,
    obj::Dict{Int,Battery},
    sN_B::Vector{Int}
)
    EBAT_THR = model[:EBAT_THR]
    
    # Add an extra slack variable to remove the non-linearities of the penalty function φ(.) 
    # JuMP.@variable(model, EBAT_slack[j in s[:sN_B], t in s[:sT]] >= 0)

    # this constraint is useful if you want to penalize having battery SoE close to each min and max limits.
    #JuMP.@constraints(model, begin
    #    cBatSoE_softLB[j in s[:sN_B], t in s[:sT]],
    #        EBAT_slack[j,t] >= bat_df[bat_df.battery_node .== j,:soc_min_strict][1] * bat_df[bat_df.battery_node .== j,:capacity_pu][1] - EBAT[j,t]
    #    cBatSoE_softUB[j in s[:sN_B], t in s[:sT]],
    #        EBAT_slack[j,t] >=  - bat_df[bat_df.battery_node .== j,:soc_max_strict][1] * bat_df[bat_df.battery_node .== j,:capacity_pu][1] + EBAT[j,t]
    #end)

    JuMP.@expression(model, bat_penalty[i in sN_B],
        obj[i].w * EBAT_THR[i]      # w1 * sum(EBAT_slack[bat_df[i,:battery_node],t] for t in s[:sT])
    )
    return nothing
end

function add_battery_formulation(
    model::JuMP.Model,
    obj::Dict{Int,Battery},
    sN_B::Vector{Int},
    sT::Vector{Int},
    dt::Float64,
    ::IncludeBatteryConfig
)
    
    add_battery_variables(model, sN_B, sT)
    add_battery_constraints(model, obj, sN_B, sT, dt)
    add_battery_objective(model, obj, sN_B)

    return nothing
end

function add_battery_formulation(
    model::JuMP.Model,
    obj::Dict{Int,Battery},
    sN_B::Vector{Int},
    sT::Vector{Int},
    dt::Float64,
    ::NoBatteryConfig
)

    return nothing
end

function get_battery_soc_from_rep_days(battery_daily_soc::Dict{Int,Vector{Float64}}, clusters::Dict{String,Cluster})
    days_to_curve = Dict{Int,Int}()
    for cluster in values(clusters)
        if length(cluster.days) <=3
            for i in keys(cluster.days)
                days_to_curve[i] = i
            end
        else
            for i in keys(cluster.days)
                if i in cluster.repdays
                    days_to_curve[i] = i
                else
                    if 3 <= length(cluster.repdays) <= 4    
                        days_to_curve[i] = cluster.repdays[2]
                    elseif length(cluster.repdays) == 5
                        days_to_curve[i] = cluster.repdays[3]
                    end
                end
            end
        end
    end

    yearly_soc = Float64[]
    for key in sort(collect(keys(days_to_curve)))
        append!(yearly_soc, battery_daily_soc[days_to_curve[key]])
    end

    return yearly_soc
end

# degradation
function yearly_battery_degredation(batteries_daily_soc::AbstractDict, clusters::Dict{String,Cluster}, batteries::Vector{Battery}, ref::Dict{Symbol,<:Any})
    nbat = length(batteries)
    representative_days = keys(batteries_daily_soc)
       
    # This is the capacity degradation. So the new capacity should be (1-deg)*init_capacity
    capacity_deg = zeros(nbat)
    for (i, battery) in enumerate(batteries)
        battery_soc_per_day = Dict(j => batteries_daily_soc[j][battery.id] for j in representative_days)
        # create yearly battery soc from representative day curves. Size: (96*365)
        battery_yearly_soc = get_battery_soc_from_rep_days(battery_soc_per_day, clusters)
        # calculate total yearly degredation percentage from soc data
        capacity_deg[i] = battery_capacity_fade(battery_yearly_soc, ref[:params]["Resolution"] , ref[:params]["Dt_hr_adjust"])
    end
   
    return capacity_deg
end

function get_battery_soc_from_year_data(batteries_daily_soc::AbstractDict, key::Int)
    days = sort(collect(keys(batteries_daily_soc)))
    yearly_soc = Float64[]
    for day in days
        append!(yearly_soc, batteries_daily_soc[day][key])
    end
    
    return yearly_soc
end
    
function yearly_battery_degredation(batteries_daily_soc::AbstractDict, batteries::Vector{Battery}, ref::Dict{Symbol,<:Any})       
    capacity_deg_percent = Dict(bat.id => 0.0 for bat in batteries)
    for battery in batteries
        battery_yearly_soc = get_battery_soc_from_year_data(batteries_daily_soc, battery.id)
        # calculate total yearly degredation percentage from soc data
        capacity_deg_percent[battery.id] = 100 * battery_capacity_fade(battery_yearly_soc, ref[:params]["Resolution"] , ref[:params]["Dt_hr_adjust"])
    end
   
    return capacity_deg_percent
end

function degradate_batteries!(batteries::Vector{Battery}, capacity_deg_percent::Vector{Float64})
    for (i,battery) in enumerate(batteries)
        # update battery's total capacity
        new_capacity = (1 - capacity_deg_percent[i]) * battery.capacity_pu_degraded[]
        battery.capacity_pu_degraded[] = new_capacity

        # update battery's efficiency
        new_efficiency = battery.efficiency_degraded[] - 0.2303 * capacity_deg_percent[i]
        battery.efficiency_degraded[] = new_efficiency
    end

    return nothing
end