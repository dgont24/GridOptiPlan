abstract type AbstractOptimizationConfig end

abstract type AbstractLimitViolationConfig <: AbstractOptimizationConfig end
struct SoftLimitsConfig <: AbstractLimitViolationConfig end # include slack variables and allow limit violations (thermal, voltage)
struct HardLimitsConfig <: AbstractLimitViolationConfig end # don't include slack variables. Limit violations not allowed

abstract type AbstractEVConfiguration <: AbstractOptimizationConfig end
struct LinearEVConfig <: AbstractEVConfiguration end
struct BinaryEVConfig <: AbstractEVConfiguration end

abstract type AbstractTypeConfiguration end

abstract type AbstractOLTCConfiguration <: AbstractTypeConfiguration end
struct IncludeOLTCConfig <: AbstractOLTCConfiguration end
struct NoOLTCConfig <: AbstractOLTCConfiguration end

abstract type AbstractSVRConfiguration <: AbstractTypeConfiguration end
struct IncludeSVRConfig <: AbstractSVRConfiguration end
struct NoSVRConfig <: AbstractSVRConfiguration end

abstract type AbstractBatteryConfiguration <: AbstractTypeConfiguration end
struct IncludeBatteryConfig <: AbstractBatteryConfiguration end
struct NoBatteryConfig <: AbstractBatteryConfiguration end


mutable struct OptiSettings
    threads_milp::Union{Nothing,Int}
    iterations::Int
    population_size::Int
    investment_year::Int
    first_flexibility_year::Int
    final_load_growth_year::Int
    final_year::Int
    interest_rate::Float64
    assets_life_expectancy::Int
    upgrades::Vector{Symbol}
    violations_config::AbstractLimitViolationConfig
    ev_config::AbstractEVConfiguration
    output_dir::String

    function OptiSettings(
        jsonpath::String,
        output_dir::String, 
        threads_milp::Union{Int,Nothing}=nothing,
        iter::Int=1,
        population::Int=1
    )
        # size of gurobi envs == population_size =>
        # I limit the population_size <= 150 bc I am not sure how many environments I am allowed to have.
        if population > 150
            throw(ArgumentError("Maximum allowed population size is 150"))
        end

        # get simulation settings from json file
        data = JSON3.read(jsonpath)

        # validate the given upgrades vector
        sorted_upgs::Vector{Symbol} = validated_upgrades(collect(data[:upgrades]))

        violations = string_to_config(get(data, :violations_config, "None"), AbstractLimitViolationConfig)
        ev = string_to_config(get(data, :ev_config, "None"), AbstractEVConfiguration)

        # Create output directory
        isdir(output_dir) || mkdir(output_dir)
       
        new( 
            threads_milp, iter, population,
            data[:investment_year], data[:first_flexibility_year], data[:final_load_growth_year], data[:final_year], 
            data[:interest_rate], data[:assets_life_expectancy],
            sorted_upgs, violations, ev, output_dir
        ) 
    end
end

function string_to_config(type_string::String, abstract_type::Type{T}) where {T<:AbstractOptimizationConfig}
    for config in InteractiveUtils.subtypes(abstract_type)
        if string(nameof(config)) == type_string
            return config()
        end
    end
    error("Invalid type of configuration")
end

abstract type AbstractType end

struct TapChangerType <:AbstractType
    code::Int
    name::String                   # type name
    reg::Float64                   # regulation range eg +-10%
    ns::Int64                      # total steps count
    tmin::Float64                  # Minimum tap ratio
    tmax::Float64                  # Maximum tap ratio
    dt::Float64                    # oltc turn ratio change per tap
    ns_binary_length::Int64        # the length of the binary representation of NS
    max_n::Int64                   # maximum operations without maintenance
    cinv::Float64                  # investement cost of the oltc in Euro
    cmain::Float64                 # maintenance costs. Usually add 20% of Cinv 
    
    function TapChangerType(row::T) where {T<:DataFrames.DataFrameRow}
        tmin = 1 - row.reg
        tmax = 1 + row.reg
        dt = 2 * row.reg / row.ns
        ns_binary_length = floor(Int, log2(row.ns)) # Should I remove +1 here? 8 is represented by 4 digits(1000). But for 0-7 you just need three
                                                    # I did so as to avoid the extra constraint 0<=tap<=tap_max
        new(row.code, row.name, row.reg, row.ns, tmin, tmax, dt, ns_binary_length, row.max_n, row.cinv, row.cmain)
    end
end

struct BatteryType <: AbstractType
    code::Int
    name::String
    c_rate::Float64
    soc_min::Float64
    soc_max::Float64
    efficiency::Float64
    costcapacity_eur_per_kwh::Float64
    costpower_eur_per_kw::Float64
    costmaint_eur_per_year_per_kWh::Float64
    lifetime_years::Int
    cycles::Int
    end_life_deg_percent::Float64

    function BatteryType(row::T) where {T<:DataFrames.DataFrameRow}
        new(
            row.code, row.name, row.c_rate, row.soc_min, row.soc_max, row.efficiency, 
            row.costcapacity_eur_per_kwh, row.costpower_eur_per_kw, row.costmaint_eur_per_year_per_kWh, 
            row.lifetime_years, row.cycles, row.end_life_deg_percent
        )
    end
end


abstract type AbstractEquipment end

struct TapChanger <:AbstractEquipment 
    id::Int
    tap_node::Int
    hv_node::Int
    lv_node::Int
    data::TapChangerType
end

struct Battery <:AbstractEquipment
    id::Int
    node::Int
    capacity::Float64
    capacity_pu_init::Float64
    capacity_pu_degraded::Base.RefValue{Float64}
    efficiency_degraded::Base.RefValue{Float64}
    soc_start::Float64
    w::Float64
    data::BatteryType
    cap_cost::Float64
    op_cost_yearly::Float64
end

struct DayData
    year::Int
    id::Int
    cluster::String
    fvalue::Float64
end

struct RepDayData
    year::Int
    id::Int
    cluster::String
    type::String
    fvalue::Float64
end

struct Cluster
    year::Int
    name::String
    days::Dict{Int,DayData}
    repdays_dict::Dict{Int, String}
    repdays::Vector{Int}
end