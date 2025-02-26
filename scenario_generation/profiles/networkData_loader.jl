import CSV
using DataFrames


function renameColumns(df::DataFrame)
    # By convention make column names lowercase and replace spaces with "_"
    new_column_names = Dict(old => replace(lowercase(old), " " => "_") for old in names(df))
    rename!(df, new_column_names)
    return df
end

function leafNodes(branch_start, branch_end)
    leaf_nodes = Int64[]
    for x in branch_end
        if !(x in branch_start)
            push!(leaf_nodes, x) 
        end
    end
    return leaf_nodes
end

function junctionNodes(branches::DataFrame)
    # Identify junction nodes
    counts_df = combine(groupby(branches, :from_bus), nrow => :count)
    junction_nodes = counts_df[counts_df.count .> 1, :from_bus]
    return junction_nodes
end

function leafNodeHasLoad(leaf_nodes, load_buses)
    for i in leaf_nodes
        if i ∉ load_buses
            error("A leaf node is not connected to any load")
        end
    end
end

function noLoadLeafNodes(branch_start, branch_end, load_buses)
    dummy_nodes = Int64[]
    for x in branch_end
        if !(x in branch_start) && !(x in load_buses)
            push!(dummy_nodes, x) 
        end
    end
    return dummy_nodes
end

# Function to combine consecutive lines with the same linecode
function combineBranches(branches::DataFrame, load_buses)
    # Identify junction nodes
    junction_nodes = junctionNodes(branches)
    # Identify leaf nodes 
    leaf_nodes = leafNodes(branches[:, :from_bus], branches[:, :to_bus])
    # Check that all leaf_nodes are connected to loads
    leafNodeHasLoad(leaf_nodes, load_buses)
    
    # Initialize a list to store the combined branches
    combined_branches = DataFrame(
        from_bus = Int64[],
        to_bus = Int64[],
        length = Float64[],
        linecode = String[]
    )
    # Initialize variables to store information about the current branch
    current_from_bus = 1
    current_to_bus = 2
    current_length = 0.0
    current_linecode = ""

    # Initialize vector to store merged nodes
    merged_nodes = Int64[]
    # Add id column in branches DataFrame
    branches.id = 1:nrow(branches)
    for i in 1:nrow(branches)
        start_node = branches[i, :from_bus]
        middle_node = branches[i, :to_bus]

        if start_node in merged_nodes
            continue
        elseif (middle_node in junction_nodes) || (middle_node in leaf_nodes)
            current_from_bus = start_node
            current_to_bus = middle_node
            current_length = branches[i, :length]
            current_linecode = branches[i, :linecode]
        else 
            # index to next branch
            j = branches[branches.from_bus .== middle_node, :id][1]
            end_node = branches[j, :to_bus]

            current_from_bus = start_node
            current_to_bus = middle_node
            current_length = branches[i, :length]
            current_linecode = branches[i, :linecode] 

            canMerge = true
            while canMerge
                if canMergeBranch(branches[i, :linecode], branches[j, :linecode])
                    # Update the current branch
                    current_to_bus = end_node
                    current_length += branches[j, :length]
                    push!(merged_nodes, middle_node)

                    middle_node = end_node
                    if mergeableNode(middle_node, junction_nodes, merged_nodes, leaf_nodes)
                        j = branches[branches.from_bus .== middle_node, :id][1]
                        end_node = branches[j, :to_bus]
                    else
                        canMerge = false
                    end 
                else
                    canMerge = false
                end
            end        
        end    

        push!(combined_branches, (current_from_bus, current_to_bus, current_length, current_linecode))
    end
    return combined_branches
end

function mergeableNode(node, junction_nodes, merged_nodes, leaf_nodes)
    if (node in merged_nodes) || (node in leaf_nodes) || (node in junction_nodes)
        return false
    else
        return true
    end
end

function canMergeBranch(linecode1, linecode2)
    if linecode1 == linecode2
        return true
    else 
        return false
    end
end

function sanitizeBranches(branches::DataFrame, linecodes::DataFrame, load_buses)
    renameColumns(branches)
    df = branches
    #Remove Branches Not Connected To Loads
    #Define initial set of leaf nodes without any load
    dummy_nodes = noLoadLeafNodes(df[:, :from_bus], df[:, :to_bus], load_buses)
    #Define binary variable to check if there are any extra branches
    extraBranches = false
    if length(dummy_nodes) != 0
        extraBranches = true
    end
    #Iterate until all extra branches are removed from the dataset
    while extraBranches == true
        subset!(df, :to_bus => x -> (!in).(x,Ref(dummy_nodes)))
        dummy_nodes = noLoadLeafNodes(df[:, :from_bus], df[:, :to_bus], load_buses)
        if length(dummy_nodes) == 0
            extraBranches = false
        end
    end

    # Combine Sequential Branches With Same LineCode
    df = combineBranches(df, load_buses)

    # Add r,x columns
    df.r = linecodes[linecodes.name .== df.linecode, :r1]
    df.x = linecodes[linecodes.name .== df.linecode, :x1]
    df.cont_rating = linecode[linecodes.name .== df.linecode, :cont_rating]
    df.id = 1:nrow(df)

    return df
end

function sanitizeLinecodes(linecodes::DataFrame)
    renameColumns(linecodes)

    # adapt single phase lines current rating to three phase

    current_ratings = Dict( "2c_.007"       =>  1, 
                            "2c_.0225"      =>  2,
                            "35_SAC_XSC"    =>  2,
                            "4c_.06"        =>  2,
                            "4c_.1"         =>  2,
                            "4c_.35"        =>  2,
                            "4c_185"        =>  2,
                            "4c_70"         =>  2,
                            "4c_95_SAC_XC"  =>  2)
    return linecodes

end

function sanitizeLoads(loads::DataFrame)
    renameColumns(loads)
end

function sanitizeBuses(buses::DataFrame)
    renameColumns(buses)
end


# define root path
ROOT_DIR = pwd()
DATA_PATH = joinpath(ROOT_DIR,raw"data\European_LV_Test_Feeder_v2\European_LV_CSV")

# load csv files
original_branches = CSV.read(joinpath(DATA_PATH,"Lines.csv"), DataFrame)
original_buses = CSV.read(joinpath(DATA_PATH,"bus.csv"), DataFrame)
original_loads = CSV.read(joinpath(DATA_PATH,"Loads.csv"), DataFrame)
original_linecodes = CSV.read(joinpath(DATA_PATH,"LineCodes.csv"), DataFrame)

# pre-process dataframes
loads = sanitizeLoads(original_loads)
linecodes = sanitizeLinecodes(original_linecodes)
branches = sanitizeBranches(original_branches, linecodes, loads[:, :bus])
buses = sanitizeBuses(original_buses)


@enter
CSV.write("sanitizedBranchesFinal.csv", df)













