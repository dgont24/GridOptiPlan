function save_ga_results(path::String, arr::AbstractArray{<:Any, 3}, iteration::Int)
    CSV.write(joinpath(path,"capital_cost.csv"),          DataFrames.DataFrame([arr[:,iteration,1]], :auto), append=true)
    CSV.write(joinpath(path,"operational_cost.csv"),      DataFrames.DataFrame([arr[:,iteration,2]], :auto), append=true)
    CSV.write(joinpath(path,"total_cost.csv"),            DataFrames.DataFrame([arr[:,iteration,3]], :auto), append=true)
    CSV.write(joinpath(path,"termination_status.csv"),    DataFrames.DataFrame([arr[:,iteration,4]], :auto), append=true)
    CSV.write(joinpath(path,"overall_time.csv"),          DataFrames.DataFrame([arr[:,iteration,5]], :auto), append=true)
    CSV.write(joinpath(path,"network_configuration.csv"), DataFrames.DataFrame([arr[:,iteration,6]], :auto), append=true)

    return nothing
end

# TODO: create function to check if the flex_ev_up and flex_ev_dn are not !=0 at the same time

# create output folder
function mk_output_dir(output_path::String)
    timestamp = Dates.format(Dates.now(), "YYYYmmdd-HHMMSS")
    dir_name = joinpath(output_path, "run_$timestamp")
    @assert !ispath(dir_name) "Somebody else already created the directory"
    mkpath(dir_name)
    return dir_name
end

function save_dict_to_json(dict::AbstractDict, path::String, filename::String)
    open(joinpath(path,filename*".json"), "w") do io
        JSON3.pretty(io, dict)
    end
    return nothing
end