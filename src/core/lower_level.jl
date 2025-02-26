function initialize_daily_optimization_problem!(
    model::JuMP.Model,
    ref::Dict{Symbol,<:Any},
    lref::Dict{Symbol,<:Any}
)
    s::Dict{Symbol, Vector{Int64}} = lref[:sets]
    # Set all the constants to zero for now, and they will be modified later for each day with set_normalized_rhs()
    p_load = 0.0
    p_gen = 0.0
    # penalty for violating voltage and and thermal constraints
    limit_penalty = 1000

    #DECISION VARIABLES
    JuMP.@variables(model, begin
        # votlage
        vVOLT2[union(s[:N_REF], s[:N], s[:N_SVR], s[:N_OLTC]), s[:T]]  # bus voltage magnitudes. BE AWARE those are the SQUARE MAGNITUDES !!        
        
        # power flows
        PFLOW[i in union(s[:N_REF], s[:N]), j in ref[:node_childs][i], t in s[:T]]    # active power flows between all connected nodes
            ## maybe in the future i change it so that PFLOW is reference by i,t only where i will be the line index.
            ## added a value p[(0,1)] for the power coming from the grid 
        QFLOW[i in union(s[:N_REF], s[:N]), j in ref[:node_childs][i], t in s[:T]]
        
        # nodal power injections
        PNODE[s[:N], s[:T]]           # Nodal total Active Power injection (pu)
        QNODE[s[:N], s[:T]]           # Nodal total Reactive Power injection (pu)
        
        # apc
        PGEN[s[:N_L], s[:T]]  >= 0    # Final Power generation for DERs (pu)
        PCURT[s[:N_G], s[:T]] >= 0    # Active Power curtailment (pu)
        # TODO: ADD a Variable DPCURT = PCURT[j, t] - PCURT[j, t-1], to limit the ramp rate of the ACP, so that is more similar to droop control
            # The constraint should be as follows: -DPCURT_max <= DPCURT <= DPCURT_max 
        
        # rpc
        QGEN[s[:N_L], s[:T]]          # Reactive Power Control of DERs
        QCTRL[s[:N_G], s[:T]]         # Absolute value of QGEN

        # slack variables 
        bus_slack[union(s[:N], s[:N_SVR], s[:N_OLTC]), s[:T]] >= 0
        line_slack[s[:LINES], s[:T]] >= 0
    end);

    # add ev formulation here because PEV is needed later
    add_ev_formulation(model, ref, lref, ref[:settings].ev_config)
    # add oltc, svr, batteries
    add_tap_changer_formulation(model, lref[:oltc], s[:OLTC], s[:T], ref[:params]["Vmin2"], ref[:params]["Vmax2"], lref[:oltc_config])
    add_tap_changer_formulation(model, lref[:svr], s[:SVR], s[:T], ref[:params]["Vmin2"], ref[:params]["Vmax2"], lref[:svr_config])
    add_battery_formulation(model, lref[:batteries], s[:N_B], s[:T], ref[:params]["Dt_hr_adjust"], lref[:bat_config])

    #CONSTRAINTS
    # fix reference bus voltage at value V_REF_PU
    for t in s[:T]
        JuMP.fix(vVOLT2[0, t], ref[:params]["V_REF_PU"]) 
    end
    
    # Voltage upper and lower limits for all buses except the reference bus
    JuMP.@constraints(model, begin
        # Maximum allowed bus voltage
        cVBusUB[l in union(s[:N], s[:N_SVR], s[:N_OLTC]), t in s[:T]],
            vVOLT2[l,t] <= (ref[:params]["Vmax2"] + bus_slack[l,t])
        # Minimum allowed bus voltage
        cVBusLB[l in union(s[:N], s[:N_SVR], s[:N_OLTC]), t in s[:T]],
            vVOLT2[l,t] >= (ref[:params]["Vmin2"] - bus_slack[l,t])
    end)

    # add voltage constraint for transformer edge
    JuMP.@constraint(model, cTrafoFlow[t in s[:T]],
        vVOLT2[lref[:sourcebus], t] - vVOLT2[1,t] == 
            2 * (ref[:trafo][1,:r_pu] * PFLOW[0, 1, t] + ref[:trafo][1,:x_pu] * QFLOW[0, 1, t])
    )

    # add voltage constraints for lines edges
    JuMP.@constraint(model, cLineFlows[l in s[:LINES], t in s[:T]],
        vVOLT2[lref[:lines][l,:newf_bus],t] - vVOLT2[lref[:lines][l,:to_bus],t] ==
            2 * (lref[:lines][l,:r_pu] * PFLOW[lref[:lines][l,:from_bus],lref[:lines][l,:to_bus], t]  + 
                lref[:lines][l,:x_pu] * QFLOW[lref[:lines][l,:from_bus],lref[:lines][l,:to_bus], t])
    );

    JuMP.@constraints(model, begin
        # Branches power balance
        cPBranchBalance[i in s[:N], t in s[:T]],
            sum(PFLOW[j, i, t] for j in ref[:node_parents][i]) + PNODE[i, t] == sum(PFLOW[i, j, t] for j in ref[:node_childs][i])
        cQBranchBalance[i in s[:N], t in s[:T]],
            sum(QFLOW[j, i, t] for j in ref[:node_parents][i]) + QNODE[i, t] == sum(QFLOW[i, j, t] for j in ref[:node_childs][i])

        # Nodal power balance
        cPLoadNode[j in s[:N_L], t in s[:T]],
            PNODE[j, t] - PGEN[j, t] + model[:PEV][j,t] == - p_load
        cQLoadNode[j in s[:N_L], t in s[:T]],
            QNODE[j, t] == QGEN[j, t]
        cPBatteryNode[j in s[:N_B], t in s[:T]],
            PNODE[j, t] + model[:PBAT_ch][j,t] - model[:PBAT_dis][j,t] == 0
        cPNodeNoLoad[j in setdiff(s[:N], s[:N_L], s[:N_B]), t in s[:T]],
            PNODE[j, t] == 0
        cQNodeNoLoad[j in setdiff(s[:N], s[:N_L]), t in s[:T]],
            QNODE[j, t] == 0

        # Active Power Curtailment
        cPowerGenActive[j in s[:N_G], t in s[:T]],
            PGEN[j,t] == p_gen - PCURT[j,t]
        cCurtailMax[j in s[:N_G], t in s[:T]],
            PCURT[j,t] <= p_gen
        cPowerGenNonActive[j in setdiff(s[:N_L], s[:N_G]), t in s[:T]],
            PGEN[j,t] <= 0
    
        # Reactive Power Control
        cReactiveGenLB[j in s[:N_G], t in s[:T]],
            QGEN[j,t] >= -ref[:p_g_nom][j] * tan(ref[:params]["max_phi"])
        cReactiveGenUB[j in s[:N_G], t in s[:T]],
            QGEN[j,t] <= ref[:p_g_nom][j] * tan(ref[:params]["max_phi"])    
        cReactiveCtrlPos[j in s[:N_G], t in s[:T]],
            QCTRL[j,t] >= QGEN[j,t]
        cReactiveCtrlNeg[j in s[:N_G], t in s[:T]],
            QCTRL[j,t] >= -QGEN[j,t]
        cReactiveGenRest[j in setdiff(s[:N_L], s[:N_G]), t in s[:T]],
            QGEN[j,t] == 0
    end);

    # lines thermal limits constraints
    JuMP.@constraints(model, begin
        cLineThermI[l in s[:LINES], t in s[:T]],
            PFLOW[lref[:lines][l,:from_bus],lref[:lines][l,:to_bus], t] <= (lref[:lines][l,:smax_pu] + line_slack[l,t]) - ref[:params]["a_phi2"] * QFLOW[lref[:lines][l,:from_bus],lref[:lines][l,:to_bus], t]
        cLineThermII[l in s[:LINES], t in s[:T]],
            PFLOW[lref[:lines][l,:from_bus],lref[:lines][l,:to_bus], t] <= (lref[:lines][l,:smax_pu] + line_slack[l,t]) + ref[:params]["a_phi2"] * QFLOW[lref[:lines][l,:from_bus],lref[:lines][l,:to_bus], t]
        cLineThermIII[l in s[:LINES], t in s[:T]],
            PFLOW[lref[:lines][l,:from_bus],lref[:lines][l,:to_bus], t] >= - (lref[:lines][l,:smax_pu] + line_slack[l,t]) - ref[:params]["a_phi2"] * QFLOW[lref[:lines][l,:from_bus],lref[:lines][l,:to_bus], t]
        cLineThermIV[l in s[:LINES], t in s[:T]],
            PFLOW[lref[:lines][l,:from_bus],lref[:lines][l,:to_bus], t] >= - (lref[:lines][l,:smax_pu] + line_slack[l,t]) + ref[:params]["a_phi2"] * QFLOW[lref[:lines][l,:from_bus],lref[:lines][l,:to_bus], t]
        cLineQThermI[l in s[:LINES], t in s[:T]],
            QFLOW[lref[:lines][l,:from_bus],lref[:lines][l,:to_bus], t] <= (lref[:lines][l,:smax_pu] + line_slack[l,t]) * ref[:params]["b_phi2"]
        cLineQThermII[l in s[:LINES], t in s[:T]],
            QFLOW[lref[:lines][l,:from_bus],lref[:lines][l,:to_bus], t] >= - (lref[:lines][l,:smax_pu] + line_slack[l,t]) * ref[:params]["b_phi2"]      
    end);

    # limit PFLOW[(0,1)] according to transformer smax_lines_pu
    JuMP.@constraints(model, begin
        cTrafoThermI[t in s[:T]],
            PFLOW[0,1,t] <= ref[:trafo][1,:s_pu] - ref[:params]["a_phi2"] * QFLOW[0,1,t]
        cTrafoThermII[t in s[:T]],
            PFLOW[0,1,t] <= ref[:trafo][1,:s_pu] + ref[:params]["a_phi2"] * QFLOW[0,1,t]
        cTrafoThermIII[t in s[:T]],
            PFLOW[0,1,t] >= - ref[:trafo][1,:s_pu] - ref[:params]["a_phi2"] * QFLOW[0,1,t]
        cTrafoThermIV[t in s[:T]],
            PFLOW[0,1,t] >= - ref[:trafo][1,:s_pu] + ref[:params]["a_phi2"] * QFLOW[0,1,t]
        cTrafoQThermI[t in s[:T]],
            QFLOW[0,1,t] <= ref[:trafo][1,:s_pu] * ref[:params]["b_phi2"] 
        cTrafoQThermII[t in s[:T]],
            QFLOW[0,1,t] >= - ref[:trafo][1,:s_pu] * ref[:params]["b_phi2"] 
    end);

    JuMP.@expressions(model, begin
        limit_relaxation_penalty, limit_penalty * (sum(bus_slack) + sum(line_slack))
        oltc_cost, sum(model[:costOLTC][r] for r in s[:OLTC]; init=0)
        svr_cost, sum(model[:costSVR][r] for r in s[:SVR]; init=0)
        bat_cost, sum(model[:bat_penalty][j] for j in s[:N_B]; init=0)
        apc_cost, sum(PCURT[j,t] for j in s[:N_G], t in s[:T]) * ref[:costs]["c_p"] * ref[:params]["Dt_hr_adjust"]
        qctrl_cost, sum(QCTRL[j,t] for j in s[:N_G], t in s[:T]) * ref[:costs]["c_q"] * ref[:params]["Dt_hr_adjust"]        
    end)

    # OBJECTIVE FUNCTION
    JuMP.@objective(model, Min,
    # minimize operational costs
        limit_relaxation_penalty +  # Minimize bus voltage limit and line thermal violations
        oltc_cost +                 # oltc total costs
        svr_cost +                  # svr total costs
        bat_cost +                  # battery penalties
        apc_cost +                  # APC Total Cost
        qctrl_cost +                # RPC Total Cost 
        model[:evshift_cost]         # load shifting total costs
    )

    return nothing
end

# TODO: Silence outputs (not possible?)
# Initialize gurobi environments
function init_gurobi_envs!(tot_envs::Int)
    println("Initialize Gurobi environment(s)")
    try
        for _ in 1:tot_envs
            push!(GRB_ENV_REF, Gurobi.Env())
        end
        println("End of Gurobi environment(s) initialization")
    catch
        println("Couldn't create specified number of gurobi environments")
        throw(error())
    end
end

create_optimizer(i::Int) = Gurobi.Optimizer(GRB_ENV_REF[i])

function init_model(env_id::Int, threads::Union{Int,Nothing}; verbose::Bool=false, direct::Bool=true, timelimit::Float64=100.0)
    model = if direct
        JuMP.direct_model(create_optimizer(env_id))
    else
        JuMP.Model(() -> create_optimizer(env_id))
    end
    verbose || JuMP.set_silent(model)
    !isnothing(threads) && JuMP.set_attribute(model, JuMP.MOI.NumberOfThreads(), threads)
    JuMP.set_time_limit_sec(model, timelimit)
    direct && JuMP.set_string_names_on_creation(model, false)
    
    return model  
end

function solve_multiyear_operational_problem!(
    model::JuMP.Model,
    env_id::Int,
    ref::Dict{Symbol,<:Any},
    lref::Dict{Symbol,<:Any};
    print_output::Bool
)
    # initialize global results
    cost = 0.0
    termination_status_code = 1
    operational_optimization_start_time = time()
    allowed_time_per_year = 4 * 60 # 4 min per year
    
    final_load_growth_year_is_solved = false
    final_load_growth_year_cost = 0.0

    years = (ref[:settings].first_flexibility_year):(ref[:settings].final_year)
    for (i,year) in enumerate(years)
        # skip year if it is not in the clusters => no violations occur
        if !haskey(ref[:clusters], year)
            continue
        end
        # solve the operational problem for the current year
        if year < ref[:settings].final_load_growth_year
            year_cost, year_termination_status_code = solve_yearly_operational_problem!(model, env_id, ref, lref, year, termination_status_code, print_output = print_output)
            if time() - operational_optimization_start_time > allowed_time_per_year * i
                return 1.0e9, 3
            end
        elseif year >= ref[:settings].final_load_growth_year && !final_load_growth_year_is_solved
            year_cost, year_termination_status_code = solve_yearly_operational_problem!(model, env_id, ref, lref, year, termination_status_code, print_output = print_output)
            final_load_growth_year_cost = year_cost
            final_load_growth_year_is_solved = true
        # if the year is after the final load growth year, the cost is the same as the final load growth year
        else
            year_cost = final_load_growth_year_cost
            year_termination_status_code = termination_status_code
        end 

        # update global results
        cost += year_cost * present_worth_factor(ref[:settings].interest_rate, year)
        termination_status_code = year_termination_status_code
    end
    
    return cost, termination_status_code
end

function solve_yearly_operational_problem!(
    model::JuMP.Model,
    env_id::Int,
    ref::Dict{Symbol,<:Any},
    lref::Dict{Symbol,<:Any},
    year::Int,
    termination_status_code::Int=1;
    print_output::Bool
)
    # Dict to store operational costs of each day
    days_to_run = keys(ref[:days][year])
    cost_of_day = Dict(day => 0.0 for day in days_to_run)
    isa(lref[:bat_config], IncludeBatteryConfig) && 
        (batteries_daily_soc = Dict(day => Dict{Int,Vector{Float64}}() for day in days_to_run))
    
    t_start_year = time()
    # @debug "Started_year Individual: $(env_id) Year: $(year) Time: $(now())"
    for (cluster_name, cluster) in ref[:clusters][year]
        for (day, type) in cluster.repdays_dict    
            # solve problem
            print_output && println("Solving operational problem for year $year, cluster: $cluster_name, and day $day ($type)")
            daily_cost, day_termination_status_code, violations = solve_daily_operational_problem!(year, day, model, env_id, ref, lref)

            # update results
            cost_of_day[day] = daily_cost
            termination_status_code = update_termination_status(day_termination_status_code, violations, termination_status_code)
            isa(lref[:bat_config], IncludeBatteryConfig) && 
                (batteries_daily_soc[day] = get_batteries_day_soc(model,collect(values(lref[:batteries])),ref[:params]["timesteps"]))
        end
    end
    t_run_year = round(time() - t_start_year, digits=2)
    # @debug "Finished_year Individual: $(env_id) Year: $(year) Time: $(now()) Solve_time: $(t_run_year)"

    # apply regression to get total cost
    total_cost = total_year_oper_cost(ref[:clusters][year], cost_of_day)

    # degradate batteries
    if isa(lref[:bat_config], IncludeBatteryConfig) && (termination_status_code == 1)
        bat_capacity_degradation_percent = yearly_battery_degredation(batteries_daily_soc, ref[:clusters][year], collect(values(lref[:batteries])), ref)
        degradate_batteries!(collect(values(lref[:batteries])), bat_capacity_degradation_percent)
    end

    return total_cost, termination_status_code
end

function get_batteries_day_soc(model::JuMP.Model, batteries::Vector{Battery}, timesteps::Int)
    if !JuMP.has_values(model)
        return Dict(0 => zeros(timesteps))
    end

    batteries_day_soc = Dict{Int,Vector{Float64}}()
    for battery in batteries
        ebat = JuMP.value.(model[:EBAT])[battery.node,:].data
        soc = ebat / battery.capacity_pu_degraded[]
        batteries_day_soc[battery.id] = soc
    end
    
    return batteries_day_soc
end

function save_daily_operational_problem_results(
    year::Int,
    day::Int,
    model::JuMP.Model,
    ref::Dict{Symbol,<:Any},
    lref::Dict{Symbol,<:Any},
    res_vars::Dict{Symbol,DataFrames.DataFrame},
    res_obj::DataFrames.DataFrame,
    evs_pu::Dict{Int, Vector{Float64}},
    gens_pu::Dict{Int, Vector{Float64}}
)
    # add new results to the dataframe
    results_variables!(year, day, model, ref, lref, res_vars, evs_pu, gens_pu)
    results_objective!(year, day, model, ref, lref, res_obj)

    return nothing
end

function solve_daily_operational_problem!(
    year::Int,
    day::Int,   
    model::JuMP.Model,
    env_id::Int,
    ref::Dict{Symbol,<:Any},
    lref::Dict{Symbol,<:Any},
    res_vars::Dict{Symbol,DataFrames.DataFrame}=Dict(:s=>DataFrames.DataFrame()), 
    res_obj::DataFrames.DataFrame=DataFrames.DataFrame();
    save_model::Bool=false
)
    # Active power injections at each node. Size of each array:(n_bus,96)
    loads_pu = node_daily_timeseries(ref[:load_profiles], lref[:sets][:N_L], day)
    evs_pu = node_daily_timeseries(ref[:ev_profiles], lref[:sets][:N_EV], day)
    gens_pu = node_daily_timeseries(ref[:gen_profiles], lref[:sets][:N_G], day)
    # Update the model constraints with the network consumption/generation values for the specific day
    update_const_terms(loads_pu, evs_pu, gens_pu, model, lref[:sets])

    # solve the optimization
    _solve_daily_optimization(model, env_id, day, year)

    # export results
    termination_status_code::Int = get_status_code(model)
    violated_constraints::Bool = model_is_relaxed(model)
    cost::Float64 = get_cost_value(model)
    if save_model && JuMP.is_solved_and_feasible(model)
        save_daily_operational_problem_results(year, day, model, ref, lref, res_vars, res_obj, evs_pu, gens_pu) 
    end

    return (cost, termination_status_code, violated_constraints)
end

function _solve_daily_optimization(model::JuMP.Model, env_id::Int, day::Int, year::Int)
    t_start = time()
    # @debug "Started_day. Individual: $(env_id). Year: $(year). Day: $(day)."
    JuMP.optimize!(model)
    t_run = round(time() - t_start, digits=2)
    # @debug "Finished_day. Individual: $(env_id). Year: $(year). Day: $(day). Solve time: $(t_run)"

    #if !is_solved_and_feasible(model; dual = true)
    #    error(
    #        """
    #        The model was not solved correctly:
    #        termination_status : $(termination_status(model))
    #        primal_status      : $(primal_status(model))
    #        dual_status        : $(dual_status(model))
    #        raw_status         : $(raw_status(model))
    #        """,
    #    )
    #end

    return nothing
end

#    "The problem is dual infeasible. If a primal feasible solution " *
#    "exists, the problem is unbounded. To check, set the objective " *
#    "to `@objective(model, Min, 0)` and re-solve. If the problem is " *
#    "feasible, the primal is unbounded. If the problem is " *
#    "infeasible, both the primal and dual are infeasible.",

function get_status_code(model::JuMP.Model)
    status = JuMP.termination_status(model)
    status_to_code = Dict{JuMP.MOI.TerminationStatusCode, Int}(
        JuMP.OPTIMAL => 1,
        JuMP.TIME_LIMIT => 3,
        JuMP.INFEASIBLE_OR_UNBOUNDED => 4,
        JuMP.INFEASIBLE => 4,
        JuMP.LOCALLY_INFEASIBLE => 4,
        JuMP.DUAL_INFEASIBLE => 4
    )

    return get(status_to_code, status, 5)
end

function model_is_relaxed(model::JuMP.Model)
    violations_cost = JuMP.value(model[:limit_relaxation_penalty])
    
    return violations_cost > 0
end

function update_termination_status(current_status_code::Int, constraint_violations::Bool, global_status_code::Int)    
    violations_code = constraint_violations ? 2 : 0
    new_status::Int = max(current_status_code, violations_code, global_status_code)

    return new_status
end

function get_cost_value(model)
    if JuMP.primal_status(model) == JuMP.FEASIBLE_POINT
        cost::Float64 = JuMP.value(model[:bat_cost]) + JuMP.value(model[:apc_cost]) + JuMP.value(model[:qctrl_cost]) + JuMP.value(model[:evshift_cost])
    else
        cost = rand((1.0e6:5.0e6))
    end
   
    return cost
end

function update_const_terms(
    p_load_pu::Dict{Int64, Vector{Float64}}, 
    p_ev_pu::Dict{Int64, Vector{Float64}}, 
    p_gen_pu::Dict{Int64, Vector{Float64}}, 
    model::JuMP.Model,
    s::Dict #sets
)
    for j in s[:N_L] 
        for t in s[:T]
            JuMP.set_normalized_rhs(model[:cPLoadNode][j,t], - p_load_pu[j][t])
        end
    end

    for j in s[:N_G] 
        for t in s[:T]
            JuMP.set_normalized_rhs(model[:cPowerGenActive][j,t], p_gen_pu[j][t])
            JuMP.set_normalized_rhs(model[:cCurtailMax][j,t], p_gen_pu[j][t])
        end
    end

    for j in s[:N_EV] 
        tot_ev_load = sum(p_ev_pu[j])
        JuMP.set_normalized_rhs(model[:cEVbalance][j], tot_ev_load)
    
        for t in s[:T]
            JuMP.set_normalized_rhs(model[:cEVpower][j,t], p_ev_pu[j][t])
        end
    end
    return nothing
end