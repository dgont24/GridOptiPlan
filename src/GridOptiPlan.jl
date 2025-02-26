module GridOptiPlan

# Importing external modules
import DataFrames
import CSV
import JSON3
import InteractiveUtils
import Polynomials
import Plots
import StatsPlots
using Dates

import Metaheuristics
import JuMP
import Gurobi

const _valid_upgrades = [:oltc, :svr1, :svr2, :lines, :batteries]
const _normal_cluster_types = ["normal", "limited-ov/uv"]
const _violation_cluster_types = ["ov", "uv", "ov + uv", "ov + ol", "uv + ol", "ov + uv + ol"]
const GRB_ENV_REF = Gurobi.Env[]

# Including internal modules and core functionality
include("utils/DirectedGraph.jl")
using .DirectedGraph

include("common/types.jl")
include("common/parameters.jl")
include("common/costs.jl")
include("utils/utils.jl")
include("io/load_raw_data_consts.jl")
include("io/load_raw_data.jl")
include("io/results.jl")
include("io/plotting.jl")

include("core/data.jl")
include("core/tap_changer.jl")
include("core/evs.jl")
include("core/battery.jl")
include("core/regression.jl")
include("core/rainflow.jl")
include("core/upper_level.jl")
include("core/sets.jl")
include("core/lower_level.jl")

# Export
export OptiSettings, SoftLimitsConfig, HardLimitsConfig, LinearEVConfig, BinaryEVConfig
export load_data_from_files, create_search_space, initialize_search_algorithm, optimize, cost_f


end