function add_ev_formulation(
    model::JuMP.Model,
    ref::Dict{Symbol,<:Any},
    lref::Dict{Symbol,<:Any},
    ::LinearEVConfig,
)   
    sN_L::Vector{Int} = lref[:sets][:N_L]
    sN_EV::Vector{Int} = lref[:sets][:N_EV]
    sT::Vector{Int} = lref[:sets][:T]
    p_shift::Dict{Int,Float64} = ref[:pshift]
    c_ev::Float64 = ref[:costs]["c_ev"]

    JuMP.@variables(model, begin
        PEV[j in sN_L, t in sT]   >= 0            # EV Final Controlled Load
        flex_evUp[j in sN_EV, t in sT] >= 0   # variable for ev load shifting
        flex_evDn[j in sN_EV, t in sT] >= 0   # variable for ev load shifting
        #TODO: check if it is possible both of them to be greater that zero at the same time.
    end)

    # I set those values later using set_normalized_rhs
    p_ev = 0.0
    p_ev_sum = 0.0 

    # EV Load Control
    JuMP.@constraints(model, begin
        cEVpower[j in sN_EV, t in sT],
            PEV[j,t] == p_ev + flex_evUp[j,t] * p_shift[j] - flex_evDn[j,t] * p_shift[j]
        cEVpowerRest[j in setdiff(sN_L, sN_EV), t in sT],
            PEV[j,t] <= 0
        cEVpowerUB[j in sN_EV, t in sT],
            PEV[j,t] <= p_shift[j]
        cEVrestriction[j in sN_EV, t in sT],
            flex_evUp[j,t] + flex_evDn[j,t] <= 1
    end)

    # add an intra-day constraint of load shifting
    JuMP.@constraint(model, cEVbalance[j in sN_EV],
        sum(PEV[j,t] for t in sT) == p_ev_sum    #here I have to sum all the values of p_ev for the whole day for each vehicle.
    )

    JuMP.@expression(model, evshift_cost,
        c_ev * sum((flex_evDn[j,t] * p_shift[j]) for j in sN_EV, t in sT) * ref[:params]["Dt_hr_adjust"]   # load shifting total costs
    )

    return nothing
end

#=
function add_ev_formulation(
    model::JuMP.Model,
    p_shift::Dict{Int64, Float64},
    c_ev::Union{Int64, Float64},
    s::Dict, #sets
    ::IntegerEVConfig,
)
    JuMP.@variables(model, begin
        PEV[j in sN_L, t in sT]   >= 0    # EV Final Controlled Load
        x_ev[j in sN_EV, t in sT], Bin    # Binary variable for ev load shifting
        y_ev[j in sN_EV, t in sT], Bin    # Binary variable for ev load shifting
    end)

    # I set those values later using set_normalized_rhs
    p_ev = 0
    p_ev_sum = 0 

    # EV Load Control
    JuMP.@constraints(model, begin
        cEVpower[j in sN_EV, t in sT],
            PEV[j,t] == p_ev + x_ev[j,t]*p_shift[j] - y_ev[j,t]*p_shift[j] 
        cEVpowerRest[j in setdiff(sN_L, sN_EV), t in sT],
            PEV[j,t] <= 0
        cEVpowerUB[j in sN_EV, t in sT],
            PEV[j,t] <= p_shift[j]
        cEVrestriction[j in sN_EV, t in sT],
            x_ev[j,t] + y_ev[j,t] <= 1
    end)

    # add an intra-day constraint of load shifting
    JuMP.@constraint(model, cEVbalance[j in sN_EV],
        sum(PEV[j,t] for t in sT) == p_ev_sum    #here I have to sum all the values of p_ev for the whole day for each vehicle.
    )

    JuMP.@expression(model, evshift_cost,
        sum((c_ev * y_ev[j,t] * p_shift[j]) for j in sN_EV, t in sT)   # load shifting total costs
    )

    return nothing
end
=#