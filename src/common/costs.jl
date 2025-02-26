# COSTS

function load_cost_data(s_base::Float64)
    costs = Dict{String, Float64}()

    costs["c_p"] = 400.0 / 1e6 * s_base     # active power curtailment  # 123eur/MWh
    costs["c_q"] = 15.0 / 1e6 * s_base       # reactive power control    # 1.23eur/MWh
    costs["c_ev"] = 20.0 / 1e6 * s_base     # ev load shifting          # 1.04eur/MWh

    return costs
end


function calculate_capital_cost(
    ref::Dict{Symbol, <:Any},
    lref::Dict{Symbol, <:Any}, 
    line_states::Vector{Int}
)
    # OLTC, SVR
    cost_tap_changer = 0.0
    for element in union(lref[:oltc], lref[:svr])
        object_cost = tap_changer_total_cost(element, ref[:settings])
        cost_tap_changer +=  object_cost
    end

    # Lines
    cost_lines = 0.0
    sorted_ids = sort(collect(keys(ref[:upgrades_list][:lines])))
    for (i,id) in enumerate(sorted_ids)
        if line_states[i] != 0
            length = lref[:lines][id, :length]
            linecode = lref[:lines][id, :linecode]
            object_cost = line_total_cost(linecode, length, ref[:linecodes], ref[:settings])
            cost_lines += object_cost
        end
    end

    # Batteries
    cost_batteries = 0.0
    for element in values(lref[:batteries])
        object_cost = battery_total_cost(element, ref[:settings])
        cost_batteries += object_cost
    end

    total_capital_cost =  float(cost_tap_changer + cost_lines + cost_batteries)       

    return total_capital_cost
end


# Discount factor of present worth factor is used in calculating present values. It is equal to 1/(1 + i)^t where i the interest rate.
# i: interest rate, t: period of the incurred cash flow
# (P/F, i%, n)
present_worth_factor(i::Float64, t::Int) = 1 / (1 + i)^t
# Used for the levelization of all expenditures to the depreciation horizon. This levelization can also be called the annuity method.
# (A/P, i%, n)
capital_recovery_factor(i::Float64, t::Int) = i * (1 + i)^t / ((1 + i)^t - 1)
# (P/A, i%, n)
uniform_series_present_worth_factor(i::Float64, t::Int) = ((1 + i)^t - 1) / (i * (1 + i)^t)
# (F/A, i%, n),
series_compound_amount_factor(i::Float64, t::Int) = ((1 + i)^t - 1) / i

# Cost Values - 
# cinv: investment cost, 
# i: interest rate, 
# t: time of investment(year)
capex_npv(cinv::Float64, i::Float64, t::Int) = cinv * present_worth_factor(i, t)

# cop: operational costs, 
# i: interest rate, 
# dt: time of use(in years), 
# tend: final year that operational cost occurs
opex_npv(cop::Float64, i::Float64, dt::Int, tend::Int) = cop * series_compound_amount_factor(i, dt) * present_worth_factor(i, tend)

# cinv: investment cost, 
# i: interest rate, 
# ttot: total depreciation period of device(lifetime), 
# trem: remaing years of device lifetime after it is stopped being used, 
# tend: year that the device stops being used
residual_value_npv(cinv::Float64, i::Float64, ttot::Int, trem::Int, tend::Int) = 
    cinv * capital_recovery_factor(i, ttot) * uniform_series_present_worth_factor(i, trem) * present_worth_factor(i, tend)
    
function tap_changer_total_cost(element::TapChanger, info::OptiSettings)
    t_used = info.final_year - info.investment_year
    t_remaining = info.assets_life_expectancy - t_used

    capex = capex_npv(element.data.cinv, info.interest_rate, info.investment_year)
    opex = opex_npv(element.data.cmain, info.interest_rate, t_used, info.final_year)
    # TODO: calculate remaining cost based on the usage
    residual_value = residual_value_npv(element.data.cinv, info.interest_rate, info.assets_life_expectancy, t_remaining, info.final_year)
    totex_npv = capex - residual_value + opex

    return totex_npv
end

#TODO: Change the formula so that the residual value is calculated based of the State of Health of the battery and the end of the planning period
# as calculated by the operatinal level
function battery_total_cost(element::Battery, info::OptiSettings)
    t_used = info.final_year - info.investment_year
    t_remaining = element.data.lifetime_years - t_used

    capex = capex_npv(element.cap_cost, info.interest_rate, info.investment_year)
    opex = opex_npv(element.op_cost_yearly, info.interest_rate, t_used, info.final_year)
    residual_value = residual_value_npv(element.cap_cost, info.interest_rate, element.data.lifetime_years, t_remaining, info.final_year)
        
    totex_npv = capex - residual_value + opex

    return totex_npv
end

function line_total_cost(code::AbstractString, length::Float64, linecodes::DataFrames.DataFrame, info::OptiSettings)
    t_used = info.final_year - info.investment_year
    t_remaining = info.assets_life_expectancy - t_used

    row  = linecodes[linecodes.name .== code, :][1, :]
    cinv = length * (row.installation_cost_per_m + row.capital_cost_per_m)

    capex = capex_npv(cinv, info.interest_rate, info.investment_year)
    residual_value = residual_value_npv(cinv, info.interest_rate, info.assets_life_expectancy, t_remaining, info.final_year)
    totex_npv = capex - residual_value

    return totex_npv
end