#####
# Q: Where should the DiGraph be created? Before of after modifying line nodes by adding the SVR?
# A: for now the graph is only used in the power balance equations where the VR nodes are not considered. So before modifying any lines.

function create_graphDict(from_bus::Vector{Int}, to_bus::Vector{Int})
    g = create_di_graph(from_bus, to_bus)
    # child nodes
    node_childsDict = create_childs_Dict(g)
    # add the slack bus and connect it to the first node
    node_childsDict[0] = [1]
    # parent nodes
    node_parentsDict  = create_parents_Dict(g)
    # add the slack bus and connect it to the first node
    node_parentsDict[1] = [0]

    return node_parentsDict, node_childsDict
end

# perform a check on the give network upgrades vector and return a new one on the correct order
function validated_upgrades(upgs::Vector{String})
    upgrades::Vector{Symbol} = Symbol.(upgs)
    tot_upgrades::UnitRange{Int64} = 1:length(_valid_upgrades)
    sorted_upgrades_dict = Dict{Symbol,Int}(zip(_valid_upgrades, tot_upgrades))

    if isempty(upgrades)
        error("Provided updates list is empty.")
    end

    # Check if all given upgrade items are valid
    if !all(item ∈ _valid_upgrades for item in upgrades)
        error("List contains invalid upgrades. Possible items to add in the list are: :oltc, :svr, :lines, :batteries")
    end

    # Verify all items in upgrades are unique
    if !allunique(upgrades)
        error("List contains duplicate items")
    end

    # Sort the upgrades list
    sorted_list = sort(upgrades, by = x -> sorted_upgrades_dict[x])

    return sorted_list
end


"""
    items_per_category(upgrades)

Create a Dict{Symbol,Int} with the total number of possible upgrades for each category/type
"""
function items_per_category(upgrades::Dict{Symbol, Dict{Int64}})
    tot_items = Dict{Symbol,Int}()
    
    # Fill the dictionary with the lengths of each upgrade type
    for type in _valid_upgrades
        tot_items[type] = length(upgrades[type])
    end

    return tot_items 
end

function find_cluster_of_day(clusters::Dict{Int64, Dict{String, Cluster}}, year::Int, day::Int)
    clusters_dict = clusters[year]
    for cluster in values(clusters_dict)
        if day in keys(cluster.days)
            return cluster
        end
    end

    return nothing
end

function get_time_limit()
    time_limit_str = get(ENV,"TIME_LIMIT_HR", "25")
    time_limit = parse(Float64, time_limit_str)
    # Subtract 40 minutes to avoid timeout (assuming each iteration takes at max 40min)
    # Also add a minimum of 10 seconds
    time_limit_s = max(10, time_limit * 3600 - 40*60) 
    
    return time_limit_s
end

function print_search_space(
    space::Metaheuristics.SearchSpaces.BoxConstrainedSpace, 
    upgrades::Vector{Symbol}, 
    upgrades_nitems::Dict{Symbol, Int}
)
    lb = decompose_x_vector(space.lb, upgrades, upgrades_nitems)
    ub = decompose_x_vector(space.ub, upgrades, upgrades_nitems)

    bounds_info = IOBuffer()
    println(bounds_info, "Bounds:")
    for key in keys(lb)
        println(bounds_info, "$key:")
        println(bounds_info, "  lb: ", lb[key])
        println(bounds_info, "  ub: ", ub[key])
    end

    @info String(take!(bounds_info))

    return nothing
end

function print_algorithm_info(algo::Metaheuristics.Algorithm)
    @info """"
    GA Options: 
        iterations=$(algo.options.iterations), 
        population=$(algo.parameters.initializer.N),
        parallel_evaluation=$(algo.options.parallel_evaluation), 
        time_limit=$(Dates.canonicalize(Dates.Second(round(algo.options.time_limit)))), 
        f_calls_limit=$(algo.options.f_calls_limit)
    """
    return nothing
end

function fill_undefined_rows!(arr::Array{<:Any, 3}, i::Int, j::Int)
    for a in 1:size(arr, 1)
        for c in 1:size(arr, 3)
            arr[a, i:j, c] .= missing
        end
    end    
    return nothing
end
