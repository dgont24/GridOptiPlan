module DirectedGraph

import Graphs

export create_di_graph, create_childs_Dict, create_parents_Dict

function create_di_graph(fbus::Vector{Int}, tbus::Vector{Int})
    if length(fbus) != length(tbus)
        throw(ArgumentError("Start and end vectors should have the same size."))
    end
    
    # Create a directed graph
    g = Graphs.DiGraph(length(fbus) + 1)
    
    # Add edges from df
    for (i,j) in zip(fbus, tbus)
        Graphs.add_edge!(g, i, j)
    end
    # TODO: above add_edge just return true or false! So I can't know if all the network edges have been added correctly.
    return g
end

# Function to get parents of a node in a directed graph
function get_parents(graph, node::Int)
    parents = Int[]
    for incoming_edge in Graphs.edges(graph)
        if Graphs.dst(incoming_edge) == node
            push!(parents, Graphs.src(incoming_edge))
        end
    end
    return parents
end

# Function to get children of a node in a directed graph
function get_children(graph, node::Int)
    children = Int[]
    for outgoing_edge in Graphs.outneighbors(graph, node)
        push!(children, outgoing_edge)
    end
    return children
end

# Function to create children's dictionary for all nodes in a directed graph
function create_childs_Dict(graph)
    childsDict = Dict{Int64, Vector{Int64}}()
    for i in Graphs.vertices(graph)
        childsDict[i] = get_children(graph, i)
    end
    return childsDict
end

# Function to create parents' dictionary for all nodes in a directed graph
function create_parents_Dict(graph)
    parentsDict = Dict{Int64, Vector{Int64}}()
    for i in Graphs.vertices(graph)
        parentsDict[i] = get_parents(graph,i)
    end
    return parentsDict
end

end