# SETS
# Creating a dictionary to hold all sets required for the MILP model
function create_optimization_sets(ref::Dict{Symbol,<:Any}, local_ref::Dict{Symbol,<:Any})

    battery_keys = sort(collect(keys(local_ref[:batteries])))::Vector{Int}

    sets = Dict(
            # time intervals
        :T => collect(1:ref[:params]["timesteps"]),
            # set of line buses/nodes and other special sets of buses
        :N => ref[:buses].id, 
        :N_REF => [0],
            # set of lines ids
        :LINES => local_ref[:lines].id,
            # set of buses with loads / evs / 
            # N_EV AND N_G SHOULD BE A SUBSET OF N_L.
        :N_L => ref[:loads].bus, 
        :N_EV => ref[:evs].bus, 
        :N_G => ref[:gens].bus,
            # set of all loads / evs / static generators
            # ONLY ONE GENERATOR/LOAD/EB CAN BE APPLIED AT A SINGLE BUS.
        :L => ref[:loads].id, 
        :EV => ref[:evs].id, 
        :G => ref[:gens].id,
            # oltc
        :N_OLTC => Int[obj.tap_node for obj in local_ref[:oltc]],
        :OLTC => Int[obj.id for obj in local_ref[:oltc]],
            # svr
        :N_SVR => Int[obj.tap_node for obj in local_ref[:svr]],
        :SVR => Int[obj.id for obj in local_ref[:svr]],  
            # batteries
        :N_B => Int[local_ref[:batteries][i].node for i in battery_keys],
        :B => Int[local_ref[:batteries][i].id for i in battery_keys]
    )
    
    return sets
end