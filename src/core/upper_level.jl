function find_bounds(upgrades::Dict{Symbol, Dict{Int64}}, nitems::Dict{Symbol,Int}, key::Symbol)
    lb = Int[]
    ub = Int[]
    
    # oltc, svr
    if key in [:oltc, :svr1, :svr2]
        push!(lb, 0)
        push!(ub, nitems[key])
    # lines
    elseif key == :lines
        # lb
        append!(lb, zeros(Int, nitems[key]))
        # ub
        dict_keys = sort(collect(keys(upgrades[key])))
        append!(ub, [Int(length(upgrades[key][i])) for i in dict_keys])
    # batteries
    elseif key == :batteries
        dict_keys = sort(collect(keys(upgrades[key])))
        # lb
        append!(lb, Int(upgrades[:batteries][i][1]) for i in dict_keys) 
        # ub
        append!(ub, Int(upgrades[:batteries][i][2]) for i in dict_keys)
    else
        throw(ArgumentError("The provided key type is not supported."))
    end

    return lb, ub
end

function create_search_space(
    upgrade_types::Vector{Symbol}, 
    upgrades::Dict{Symbol, Dict{Int64}}, 
    nitems::Dict{Symbol,Int}
)
    lb = Int[]
    ub = Int[]

    for type in upgrade_types
        temp_lb, temp_ub = find_bounds(upgrades, nitems, type)
        append!(lb, temp_lb)
        append!(ub, temp_ub)      
    end

    space = Metaheuristics.BoxConstrainedSpace(lb, ub)

    return space
end

# define settings for genetic algorithm
function initialize_search_algorithm(
    bounds::Metaheuristics.BoxConstrainedSpace; 
    type::String="GA", 
    population_size::Int=100,
    iter::Int=0,
    isparallel::Bool = false,
    time_limit_s::Float64 = 24*60*60
)
    if type != "GA"
        throw(ArgumentError("The specified algorithm is currently not supported"))
    end

    algorithm_options = Metaheuristics.Options(iterations = iter, parallel_evaluation=isparallel, time_limit=time_limit_s)

    ga = Metaheuristics.GA(; N=population_size,
            crossover=Metaheuristics.SBX(;bounds),
            mutation=Metaheuristics.PolynomialMutation(;bounds),
            environmental_selection=Metaheuristics.GenerationalReplacement(),
            options = algorithm_options
        )

    return ga
end

struct NetworkConfiguration
    x_vector::Vector{Int64}
    cap_cost::Float64
    oper_cost::Float64
    total_cost::Float64
    code::String
end

function NetworkConfiguration(x::Vector{Int64}, cap_cost::Float64, oper_cost::Float64; iter::Int, individual_id::Int)
    total_cost = cap_cost + oper_cost
    code = "$(iter)-$(individual_id)"
    return NetworkConfiguration(x, cap_cost, oper_cost, total_cost, code)
end

# Define a struct to store the best performing candidates
# maxcost: total cost of the most expensive candidate
# maxcost_code: code of the most expensive candidate
mutable struct BestUpgrades
    candidates_dict::Dict{String, NetworkConfiguration}
    maxcost::Float64
    maxcost_code::String

    function BestUpgrades(n::Int)
        candidates = [NetworkConfiguration([0], i * 1.0e9, i * 1.0e9, i * 1.0e9, "$i") for i in 1:n]
        list = Dict(candidate.code => candidate for candidate in candidates)
        maxcost = candidates[end].total_cost
        maxcost_code = candidates[end].code
        new(list, maxcost, maxcost_code)
    end
end

function get_five_best_costs(solutions::BestUpgrades)
    best_costs = [candidate.total_cost for candidate in values(solutions.candidates_dict)]
    best_costs_sorted = sort(best_costs)
    return best_costs_sorted[1:5]
end

function update_best_upgrades!(new_entry::NetworkConfiguration, solutions::BestUpgrades)
    # remove current most expensive solution and add new solution
    pop!(solutions.candidates_dict, solutions.maxcost_code)
    solutions.candidates_dict[new_entry.code] = new_entry

    # update most expensive solution
    new_maxcost = 0
    new_code = ""
    for (code,candidate) in solutions.candidates_dict
        if candidate.total_cost > new_maxcost
            new_maxcost = candidate.total_cost
            new_code = code
        end
    end
    solutions.maxcost = new_maxcost
    solutions.maxcost_code = new_code

    return nothing
end

function optimize(ref::Dict{Symbol,<:Any}, bounds::Metaheuristics.BoxConstrainedSpace, method::Metaheuristics.AbstractAlgorithm)
    @info "Start of Bi-Level Cost Optimization problem. Start_time: $(now())"
    # Counter for total iterations(current generation) of the GA
    iter_count = [1]
    start_time = time()

    # Store n best performing candidates
    best_upgrade_candidates = BestUpgrades(10)

    # Initialization of gurobi
    if method.options.parallel_evaluation
        init_gurobi_envs!(ref[:settings].population_size)
    else
        init_gurobi_envs!(1)
    end
    
    # Choose the appropriate cost function based on whether parallel evaluation is enabled
    closure_cost_f = if method.options.parallel_evaluation
        x -> f_parallel(x, ref, iter_count, start_time, best_upgrade_candidates)
    else
        x -> cost_f(x, ref, SerialMode())
    end

    # Perform optimization using the chosen cost function
    result = Metaheuristics.optimize(closure_cost_f, bounds, method)

    if result.termination_status_code == Metaheuristics.TIME_LIMIT
        @info """
        Time limit of: $(method.options.time_limit / 60) minutes reached. Exiting optimization.
        Last iteration was number was $(result.iteration).
        """
    end

    total_solve_time = canonicalize(Second(round(time() - start_time)))
    @info """
    Network Expansion Optimization Process Finished. 
        End Timestamp: $(now())
        Total solve time: $total_solve_time
        Best five solutions: $(round.(get_five_best_costs(best_upgrade_candidates), digits=3))    
    """

    return result, best_upgrade_candidates
end

function f_parallel(
    x_arr::Array{T,2},
    ref::Dict{Symbol,<:Any}, 
    iter_count::Vector{Int}, 
    time_init::Float64,
    best_solutions::BestUpgrades
) where {T<:Real}

    # print and update iterations counter
    @info "Start -> Solving Generation $(iter_count[1])... [start_time: $(Time(now()))]"
    time_start = time()
    
    # Solve MILP problems of current generation in parallel using multi-threading
    # x_arr -> Arr{Int}, size = (population_size, length(x_vector))
    fitness = zeros(size(x_arr,1))  
    Threads.@threads for i in axes(x_arr, 1)
        # run low-level MILP optimization
        opt_results = cost_f(x_arr[i,:], ref, i, ParallelMode())
        
        # return objective fitness value (capital+operational cost)
        total_cost = opt_results[1] + opt_results[2]
        fitness[i] =  total_cost

        # save selected result outputs
        ref[:output][i, iter_count[1], 1] = opt_results[1]  # capital_cost
        ref[:output][i, iter_count[1], 2] = opt_results[2]  # operational_cost
        ref[:output][i, iter_count[1], 3] = total_cost      # total cost
        ref[:output][i, iter_count[1], 4] = opt_results[3]  # termination status
        ref[:output][i, iter_count[1], 5] = opt_results[4]  # overall time
        ref[:output][i, iter_count[1], 6] = opt_results[5]  # network configuration
    end

    # update best candidates
    for (i, total_cost) in enumerate(ref[:output][:, iter_count[1], 3])
        if total_cost < best_solutions.maxcost
            new_network = NetworkConfiguration(
                ref[:output][i, iter_count[1], 6], 
                ref[:output][i, iter_count[1], 1], 
                ref[:output][i, iter_count[1], 2], 
                iter=iter_count[1], 
                individual_id=i
            )
            update_best_upgrades!(new_network, best_solutions)
        end
    end

    # Save results
    save_ga_results(ref[:settings].output_dir, ref[:output], iter_count[1])
    save_dict_to_json(Dict(iter_count[1] => best_solutions.candidates_dict), ref[:settings].output_dir, "best_solutions")

    solve_time = canonicalize(Second(round(time() - time_start)))
    runtime = canonicalize(Second(round(time() - time_init)))
    no_of_best_candiates_to_print = ref[:settings].population_size >= 5 ? 5 : ref[:settings].population_size
    current_iter_best = sort(ref[:output][:, iter_count[1], 3])[1:no_of_best_candiates_to_print]
    @info """
    Finished 
              solve time: $(solve_time)
              running time since start: $(runtime)
              best five solutions of current iteration: $(round.(current_iter_best, digits=3))
              best five overall solutions: $(round.(get_five_best_costs(best_solutions), digits=3)) 
    """

    iter_count[1] += 1
    return fitness  
end

function optimize_network_candidate(
    x::Vector{T}, 
    ref::Dict{Symbol,<:Any},
    env_id::Int;
    verbose::Bool
) where {T<:Real}

    # Decompose x vector
    states = decompose_x_vector(x, ref[:settings].upgrades, ref[:upgrades_nitems])
    network_configuration = reduce(vcat, [states[i] for i in ref[:settings].upgrades])

    # Create ref dictinoray with the unique data of the current individual network configuration
    local_ref = active_network_configuration_unique_data(states, ref)

    # Calculate Operational Costs -> run lower level daily optimization problem for the corresponding planning period
    # Initialize a new model
    model = init_model(env_id, ref[:settings].threads_milp)
    # Create model variables, constraints, objectives
    initialize_daily_optimization_problem!(model, ref, local_ref)
    
    t_start_multiyear = time()
    # @debug "Start_individual Ind: $(env_id) Time: $(now()) Net: $(x)"
    operational_cost, termination_status = solve_multiyear_operational_problem!(model, env_id, ref, local_ref, print_output=verbose)
    t_run_multiyear = round(time() - t_start_multiyear, digits=2)
    # @debug "Finished_individual Individual: $(env_id) Time: $(now()) Solve_time: $(t_run_multiyear)"

    # Calculate Capital Costs
    capital_cost = calculate_capital_cost(ref, local_ref, states[:lines])

    return (capital_cost, operational_cost, termination_status, network_configuration)
end

abstract type AbstaractOptimizationMode end
struct ParallelMode <: AbstaractOptimizationMode end
struct SerialMode <: AbstaractOptimizationMode end

# f(x) = c'*x
function cost_f(
    x::Vector{T}, 
    ref::Dict{Symbol,<:Any},
    ::SerialMode
) where {T<:Real}

    t_start = time()
    println("Solving given network configuration...")

    capital_cost, operational_cost, termination_status, network_configuration = optimize_network_candidate(x, ref, 1, verbose=true)
    cost = capital_cost + operational_cost

    # Total execution time 
    overall_time = time() - t_start
    
    println("Finished optimization for the network configuration: \n$(network_configuration)")
    println("Termination status: $(termination_status)")
    println("Total time: $(overall_time)")
    
    return cost
end

# f(x) = c'*x
function cost_f(
    x::Vector{T}, 
    ref::Dict{Symbol,<:Any},
    env_id::Int,
    ::ParallelMode,
) where {T<:Real}

    t_start = time()
    
    capital_cost, operational_cost, termination_status, network_configuration = optimize_network_candidate(x, ref, env_id, verbose=false)

    # Total execution time 
    overall_time = time() - t_start

    return (capital_cost, operational_cost, termination_status, overall_time, network_configuration)
end

function decompose_x_vector(x::Vector{T}, upgrade_types::Vector{Symbol}, nitems::Dict{Symbol,Int}) where {T<:Real}
    state = Dict{Symbol, Vector{Int}}()

    current_pos = 1
    # States dict should contain all keys in _valid_upgrades. If they technology is not used then the value is an empty vector.
    for type in _valid_upgrades
        if type in upgrade_types
            if type in [:oltc, :svr1, :svr2]
                state_vec = [x[current_pos]]
                state[type] = round.(Int, state_vec)
                current_pos += 1
            else
                # get total number of items for current upgrade type
                items_count = nitems[type]
                # obtain the state vector
                state_vec = x[current_pos:(current_pos + items_count - 1)]
                # pass the state vector to the dictionary 
                state[type] = round.(Int, state_vec)
                # increase index
                current_pos += items_count
            end
        else  
            state[type] = Int[]
        end
    end

    return state
end